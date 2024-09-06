
rm(list = ls())
gc()

#### Loading files and packages ####
# Load all required packages:
library(tidyr)
library(dplyr)
library(tibble)
library(lubridate)
library(stringr)
library(ggplot2)
library(data.table)
library(R.utils)
source("/home/ivm/config_te.R")

# File paths 
file1 <- endpoint_long_file_           # ENDPOINT LONGITUDINAL FILE PATH
file2 <- det_long_file_det_long_file_  # MEDS LONGITUDINAL FILE PATH
file3 <- min_endpoint_file_            # ENDPOINT FILE PATH for SEX and Birthyear
file4 <- vnr_dict_                     # VNRO dict      


long_ep <- fread(file1); gc()
raw <- fread(file2); gc()
meds <- raw %>% filter(SOURCE == "PURCH") %>% 
  # keep only relevant columns
  select(FINNGENID, EVENT_AGE, APPROX_EVENT_DAY, ATC = CODE1, VNRO = CODE3, PCKGNR = CODE4)   
sex <- fread(file3) %>% select(FINNGENID, SEX); gc()

# Read vnro dictionary
vnro_dict <- fread(file4)
vnro_dict <- vnro_dict %>% rename(VNRO = vnr, DOSE = vahvuus_num, PKOKO = pkoko_num) %>% select(VNRO, ATC, PKOKO, DOSE)
vnro_dict$VNRO <- as.character(vnro_dict$VNRO)
vnro_dict$VNRO <- str_pad(vnro_dict$VNRO, 6, pad = "0")

# Defining the opposite of %in% to exclude rows that meet exclusion criteria
'%!in%' <- function(x,y)!('%in%'(x,y))
gc()


#____________________________________________________________________
#### Creating SDF (Summary Dataframe) ####

# Summary Dataframe for Crushing numbers
df <- NULL

##### Add C0 to SDF #####
# Add column to summary df 0
df$c0 <- length(unique(meds$FINNGENID))
df <- as.data.frame(df)



#____________________________________________________________________
#### Drugs of Interest ####

# define the medication groups of interest
trial_meds <- c("Exposure" = "B01AF01",
                "Reference" = "B01AA03")

trial_meds <- as.data.frame(trial_meds) %>% rownames_to_column(var = "drugs") %>% rename(atcs = trial_meds) %>% mutate(drugs = if_else(drugs != "Exposure", "Reference", "Exposure"))

# create a new data frame
# extract meds of relevance and replace ATC with Drug/Group names
cohort <- meds %>%
  # filter for the endpoints of interest
  filter(str_detect(ATC, paste0("^", trial_meds$atcs, collapse = "|")))%>%
  # replace ATC Codes with drug type 
  mutate(
    ATC = case_when(
      str_detect(ATC, paste0("^", trial_meds$atcs, collapse = "|")) ~ trial_meds$drugs[match(str_extract(ATC, paste0("^", trial_meds$atcs, collapse = "|")), trial_meds$atcs)],
      T ~ ATC
    )
  )

##### Add C1 to SDF #####
# Patients who used exposure or a reference
df$c1 <- length(unique(cohort$FINNGENID))

# Exposure group 
exp <- cohort %>% 
  filter(ATC == "Exposure") %>%
  # choosing each IDs first purchase
  group_by(FINNGENID) %>% 
  arrange(EVENT_AGE) %>% 
  slice(1) %>% 
  ungroup()

# first purchase of exp (was the later drug to get MA)
min_year <- year(min(exp$APPROX_EVENT_DAY)) 

# Reference group 
ref <- cohort[cohort$FINNGENID %!in% exp$FINNGENID] %>% 
  # removing IDs who started treatment before Sitagliptin came to market
  filter(year(APPROX_EVENT_DAY) >= min_year) %>%
  # choosing each IDs first purchase
  group_by(FINNGENID) %>% 
  arrange(EVENT_AGE) %>% 
  slice(1) %>% 
  ungroup()


# remove from exposure group those who are taking reference meds up to 180 days to trial init
# reference and exposure purchases from exposure group
remove_from_exp <- cohort[cohort$FINNGENID %in% exp$FINNGENID] %>% 
  # add drug initiation age
  left_join(exp %>% select(FINNGENID, init_age = EVENT_AGE)) %>% 
  # add remove column, indicating these IDs don't qualify for the study
  mutate(
    remove = if_else(ATC == "Reference" & EVENT_AGE >= (init_age - 0.5) & EVENT_AGE <= init_age, 1, 0)
  ) %>% 
  filter(remove == 1) %>%
  distinct(FINNGENID)

# Exposure group withouth prior exposure 
exp <- exp[exp$FINNGENID %!in% remove_from_exp$FINNGENID, ]

# Cohort of exposure and reference
cohort <- rbind(exp, ref)

##### Add C2 to SDF #####
# Excluded due to prior use of referent/exposure/qualified in >1 exposure category
df$c2 <- n_distinct(cohort$FINNGENID)





#____________________________________________________________________
#### Age Function ####
# FUNCTION: Add ages before drug initiation to check for meeting different inclusion/exclusion criteria 

AGES <- function(group){
  group.ages <- group %>% mutate(
    past_365  = EVENT_AGE - 1,
    past_180  = EVENT_AGE - 0.5,
    past_90   = EVENT_AGE - 0.25,
    past_30   = EVENT_AGE - (1/12),
    past_14   = EVENT_AGE - (1/24)
  ) %>% select(FINNGENID, EVENT_AGE, starts_with("past")) %>% 
    rename(drug_init_age = EVENT_AGE)
}



#____________________________________________________________________
#### Exclude missing Age and Sex ####
# Excluding IDs with missing Age

if (any(is.na(cohort$EVENT_AGE))){
  cohort <- cohort %>% filter(!is.na(cohort))
}

# Excluding IDs with SEX NA
cohort <- cohort[cohort$FINNGENID %in% (sex %>% filter(!is.na(SEX)))$FINNGENID, ]

##### Add C3 to SDF #####
# Excluded based on missing Age / Gender
df$c3 <- n_distinct(cohort$FINNGENID)


# Adding relevant Ages for further inclusion/exclusion criteria:
cohort.ages <- AGES(cohort)


# Excluding IDs under Study Age
min_study_age <- 18
cohort.ages <- cohort.ages %>% filter(drug_init_age >= min_study_age)

##### Add C4 to SDF #####
# Excluded based on Inclusion Age >= min Study Age (18y)
df$c4 <- n_distinct(cohort.ages$FINNGENID)



#____________________________________________________________________
#### Inclusion Diagnosis: AF ####
# (measured min 2x up to 365 days prior to drug init)

incl_diagnosis <- "I9_AF"

incl <- (long_ep %>% filter(ENDPOINT == incl_diagnosis))[(long_ep %>% filter(ENDPOINT == incl_diagnosis))$FINNGENID %in% cohort.ages$FINNGENID] %>% left_join(cohort.ages) %>% filter(EVENT_AGE >= past_365 & EVENT_AGE <= drug_init_age)
incl_cohort <- cohort.ages[cohort.ages$FINNGENID %in% incl$FINNGENID, ]

##### Add C5 to SDF #####
# Excluded based on Inclusion - AF (measured min 2x in past 365 days)
df$c5 <- n_distinct(incl_cohort$FINNGENID)

# All endpoints of the AF cohort 
incl_longep <- long_ep[long_ep$FINNGENID %in% incl_cohort$FINNGENID]



#____________________________________________________________________
#### Inclusion: Minimum 2 outpatient visits in prior year ####
# Defined as: min 2 registry entries within a year prior to init

# updating meds (longitudinal drug purchase file) and long_ep (longitudinal endpoint file)
incl_meds <- meds[meds$FINNGENID %in% incl_cohort$FINNGENID]
incl_longep <- incl_longep[incl_longep$FINNGENID %in% incl_cohort$FINNGENID]


min_2_visits <- incl_meds %>% left_join(cohort.ages) %>% group_by(FINNGENID) %>%  mutate(prior_year_visits = ifelse(EVENT_AGE >= past_365 & EVENT_AGE < drug_init_age, 1, 0),
                                                                                         prior_year_visits = cumsum(prior_year_visits))

incl_cohort <- incl_cohort[incl_cohort$FINNGENID %in% unique((min_2_visits %>% filter(prior_year_visits >= 2))$FINNGENID), ]

##### Add C6 to SDF #####
# Excluded based in Inclusion min 2 registry entries in past year
df$c6 <- n_distinct(incl_cohort$FINNGENID)



#____________________________________________________________________
#### Inclusion: At least having one of the following endpoints/medication in different time periods before drug init ####

# updating meds (longitudinal drug purchase file) and long_ep (longitudinal endpoint file)
incl_meds <- meds[meds$FINNGENID %in% incl_cohort$FINNGENID]
incl_longep <- incl_longep[incl_longep$FINNGENID %in% incl_cohort$FINNGENID]


# either one of the following 
# all data prior to drug init
incl_diagnosis <- c(
  "C_STROKE", 
  "I9_TIA",
  "I9_PULMEMB"
)

incl <- incl_longep %>% left_join(incl_cohort) %>% filter(ENDPOINT %in% incl_diagnosis & EVENT_AGE <= drug_init_age) %>% select(FINNGENID) %>% distinct() 


# or two of the following 
# Age >= 75
incl2 <- incl_cohort %>% filter(drug_init_age >= 75) %>% select(FINNGENID) %>% distinct()

# all data prior to drug init
incl_diagnosis <- c(
  "I9_HEARTFAIL", 
  "Q17_LVOTO_NARROW",
  "E4_DIABETES"
)

incl3 <- incl_longep %>% left_join(incl_cohort) %>% filter(ENDPOINT %in% incl_diagnosis & EVENT_AGE <= drug_init_age) %>% select(FINNGENID) %>% distinct() 

# 0 - 180 prior to drug init
incl_diagnosis <- c(
  "I9_HYPERTENSION"
)

incl4 <- incl_longep %>% left_join(incl_cohort) %>% filter(ENDPOINT %in% incl_diagnosis & EVENT_AGE >= past_180 & EVENT_AGE <= drug_init_age) %>% select(FINNGENID) %>% distinct() 

at_least_two_incl <- rbind(incl2, incl3, incl4)
at_least_two_incl <- subset(at_least_two_incl, duplicated(at_least_two_incl$FINNGENID))



incl_cohort <- rbind(incl, at_least_two_incl) %>% distinct() %>% left_join(cohort.ages)


##### Add C7 to SDF #####
# Excluded based on Inclusion Stroke Risk Factor (either one of the following: Stroke, TIA, Embolus OR two of the following: Age >= 75, HF, LVEF, DM, HTN)
df$c7 <- n_distinct(incl_cohort$FINNGENID)





#____________________________________________________________________
#### EXCLUSION CRITERIA ####

# the Exclusions will be technically performed by time frames (e.g. all exclusion in 180 days window, all exclusions in 365 days window, etc.)

# updating meds (longitudinal drug purchase file) and long_ep (longitudinal endpoint file)
incl_meds <- meds[meds$FINNGENID %in% incl_cohort$FINNGENID]
incl_longep <- incl_longep[incl_longep$FINNGENID %in% incl_cohort$FINNGENID]


#### Exclusion: One of the following Endpoints in 180 days window #####

excl_diagnosis <- c(
  "I9_STENOSIS",
  "I9_VALVES",
  "I9_ENDOCARD",
  "K11_GIBLEEDING",
  "I9_PULMEMB",
  "I9_PVT",
  "I9_THROMBOTH",
  "I9_VTE",
  "I9_PHLETHROMBDVTLOW",
  "D3_ANAEMIAINCHRONICDISEASE",
  "D3_HAEMOLYTICANAEMIA",
  "XV_PREGNANCY_BIRTH",
  "Z21_CONTRACEPTIVE_MANAG",
  "AB1_HIV",
  "N14_CHRONKIDNEYDIS",
  "DIALYSIS",
  "K11_LIVER",
  "KRA_PSY_ALCOH",
  "KRA_PSY_SUBSTANCE"
)
excl <- incl_longep %>% left_join(incl_cohort) %>% filter(ENDPOINT %in% excl_diagnosis & EVENT_AGE >= past_180 & EVENT_AGE <= drug_init_age) %>% select(FINNGENID) %>% distinct() 
incl_cohort <- incl_cohort[incl_cohort$FINNGENID %!in% excl$FINNGENID, ]

##### Add C8 to SDF #####
# Excluded based on exclusion criteria for past 180 days: Mitral Stenosis, Prosthetic Heart valve, Thrombus, Endocarditis, Bleeding, PE/DVT, Anemia, Pregnancy/Use of Contraceptive, HIV, CKD, Dialysis, Liver disease, Alcohol/Substance abuse
df$c8 <- n_distinct(incl_cohort$FINNGENID)



#____________________________________________________________________
#### Exclusion: One of the following Endpoints in 90 days window #####

# updating meds (longitudinal drug purchase file) and long_ep (longitudinal endpoint file)
incl_meds <- meds[meds$FINNGENID %in% incl_cohort$FINNGENID]
incl_longep <- incl_longep[incl_longep$FINNGENID %in% incl_cohort$FINNGENID]

excl_diagnosis <- c(
  "C_STROKE"
)

excl <- incl_longep %>% left_join(incl_cohort) %>% filter(ENDPOINT %in% excl_diagnosis & EVENT_AGE >= past_90 & EVENT_AGE <= drug_init_age) %>% select(FINNGENID) %>% distinct() 
incl_cohort <- incl_cohort[incl_cohort$FINNGENID %!in% excl$FINNGENID, ]

##### Add C9 to SDF #####
# Excluded based on Exclusion - Stroke
df$c9 <- n_distinct(incl_cohort$FINNGENID)



#____________________________________________________________________
#### Exclusion: One of the following Medication in 30 days window  #####

# updating meds (longitudinal drug purchase file) and long_ep (longitudinal endpoint file)
incl_meds <- meds[meds$FINNGENID %in% incl_cohort$FINNGENID]
incl_longep <- incl_longep[incl_longep$FINNGENID %in% incl_cohort$FINNGENID]

excl_diagnosis <- c(
  "D3_COAGDEF_PURPUR_HAEMORRHAGIC",
  "C3_BRAIN",
  "I9_ANEURYSM",
  "D3_THROMBOCYTOPENIANAS"
)

excl <- incl_longep %>% left_join(incl_cohort) %>% filter(ENDPOINT %in% excl_diagnosis & EVENT_AGE >= past_30 & EVENT_AGE <= drug_init_age) %>% select(FINNGENID) %>% distinct() 
incl_cohort <- incl_cohort[incl_cohort$FINNGENID %!in% excl$FINNGENID, ]

##### Add C10 to SDF #####
# Excluded based on exclusion criteria for past 30 days: Haemorrhagic disorder, Intracranial neoplasms, Cerebral aneurysm, Thrombocytopenia
df$c10 <- n_distinct(incl_cohort$FINNGENID)



#____________________________________________________________________
#### Exclusion: One of the following Medication in 14 days window  #####

# updating meds (longitudinal drug purchase file) and long_ep (longitudinal endpoint file)
incl_meds <- meds[meds$FINNGENID %in% incl_cohort$FINNGENID]
incl_longep <- incl_longep[incl_longep$FINNGENID %in% incl_cohort$FINNGENID]

excl_atc <- c(
  "B01AD",
  "B01AC25",
  "J02AB02",
  "P01AB01",
  "J01XD01",
  "C01BD01",
  "A02BA01",
  "A02BC01",
  "N06AB03",
  "N06CA03",
  "J05AE02",
  "J05AE03",
  "N03AF01",
  "N03AB02",
  "N03AB52",
  "J04AB02",
  "N03AA02"
)

excl1 <- incl_meds %>% 
  left_join(incl_cohort) %>% 
  filter(EVENT_AGE >= past_14 & EVENT_AGE < drug_init_age) %>% 
  filter(str_detect(ATC, paste0("^", excl_atc, collapse = "|"))) %>% 
  distinct(FINNGENID, ATC) %>% 
  select(FINNGENID)

excl_diagnosis <- c(
  "I9_TIA"
)

excl2 <- incl_longep %>% left_join(incl_cohort) %>% filter(ENDPOINT %in% excl_diagnosis & EVENT_AGE >= past_14 & EVENT_AGE <= drug_init_age) %>% select(FINNGENID) %>% distinct() 

excl <- rbind(excl1, excl2)

incl_cohort <- incl_cohort[incl_cohort$FINNGENID %!in% excl$FINNGENID, ]

##### Add C11 to SDF #####
# Excluded based on Exclusion - CYP Inhibitors/Inducers, TIA
df$c11 <- n_distinct(incl_cohort$FINNGENID)


table((cohort[cohort$FINNGENID %in% incl_cohort$FINNGENID, ])$ATC)



#____________________________________________________________________
#### FINAL COHORT ####

cohort_final <- cohort[cohort$FINNGENID %in% incl_cohort$FINNGENID, ]
#table(final_cohort$ATC)


final_cohort <- data.frame(FINNGENID = cohort_final$FINNGENID,
                           medication = cohort_final$ATC)



#____________________________________________________________________
#### Calculating CCI ####

final_cohort_long <- final_cohort %>% 
  left_join(cohort.ages) %>% 
  left_join((raw %>% select("FINNGENID", "EVENT_AGE", "CODE1", "ICDVER"))) %>% distinct()

colnames(final_cohort_long)[which(names(final_cohort_long) == "FINNGENID")] <- "ID"
colnames(final_cohort_long)[which(names(final_cohort_long) == "CODE1")] <- "primary_ICD"
colnames(final_cohort_long)[which(names(final_cohort_long) == "ICDVER")] <- "ICD_version"
colnames(final_cohort_long)[which(names(final_cohort_long) == "EVENT_AGE")] <- "Event_age"
final_cohort_long <- final_cohort_long %>% filter(ICD_version %in% c(9, 10))

score_data_180prior <- ICCI::calc_cci(final_cohort_long, exp_start=final_cohort_long %>% pull(past_180), exp_end = final_cohort_long %>% pull(drug_init_age))
table(score_data_180prior$CCI)


#____________________________________________________________________
#### Prepare Data for Matching and Survival Analysis: ####

cohort.ppsmatching <- final_cohort %>% left_join(sex) %>% left_join(cohort.ages %>% select(FINNGENID, drug_init_age)) %>% left_join(score_data_180prior, by = c("FINNGENID" = "ID"))
colnames(cohort.ppsmatching)[which(names(cohort.ppsmatching) == "drug_init_age")] <- "age"
cohort.ppsmatching[cohort.ppsmatching == "Exposure"] <- 1  # Exposure
cohort.ppsmatching$medication[cohort.ppsmatching$medication != "1"] <- 0  # Reference
cohort.ppsmatching$medication <- as.numeric(cohort.ppsmatching$medication)
cohort.ppsmatching$SEX <- as.factor(cohort.ppsmatching$SEX)


#____________________________________________________________________
#### Exclude IDs with CCI >= 10 ####
cohort.ppsmatching <- cohort.ppsmatching %>% filter(CCI <= 9)

##### Add C12 to SDF #####
# Excluded based on CCI >= 10
df$c12 <- n_distinct(cohort.ppsmatching$FINNGENID)


df
table(cohort.ppsmatching$medication)

fwrite(cohort.ppsmatching, cohort_roc_elig_)
fwrite(df, cohort_roc_elig_flowchart_)

