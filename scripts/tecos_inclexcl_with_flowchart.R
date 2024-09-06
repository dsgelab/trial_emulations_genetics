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

# File paths (potentially needed)
file1 <- endpoint_long_file_  # ENDPOINT LONGITUDINAL FILE PATH
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
trial_meds <- c("Exposure" = "A10BH01",
                "Reference1" = "A10BB01",
                "Reference2" = "A10BB07",
                "Reference3" = "A10BB09",
                "Reference4" = "A10BB12")

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
    past_91   = EVENT_AGE - 0.25
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
# Excluded based on Inclusion Age >= min Study Age (50y)
df$c4 <- n_distinct(cohort.ages$FINNGENID)



#____________________________________________________________________
#### Inclusion Diagnosis: Type 2 Diabetes ####
# (measured up to 180 days prior to drug init)

incl_diagnosis <- "T2D"

incl <- (long_ep %>% filter(ENDPOINT == incl_diagnosis))[(long_ep %>% filter(ENDPOINT == incl_diagnosis))$FINNGENID %in% cohort.ages$FINNGENID] %>% left_join(cohort.ages) %>% filter(EVENT_AGE >= past_180 & EVENT_AGE <= drug_init_age)
incl_cohort <- cohort.ages[cohort.ages$FINNGENID %in% incl$FINNGENID, ]

##### Add C5 to SDF #####
# Excluded based on Inclusion T2D
df$c5 <- n_distinct(incl_cohort$FINNGENID)

# All endpoints of the T2D cohort 
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
#### Inclusion: Taking either Pioglitazone, Metformin, Insulins or SU 90 days prior to init ####

# updating meds (longitudinal drug purchase file) and long_ep (longitudinal endpoint file)
incl_meds <- meds[meds$FINNGENID %in% incl_cohort$FINNGENID]
incl_longep <- incl_longep[incl_longep$FINNGENID %in% incl_cohort$FINNGENID]

# 1 - 90 days prior to drug init
incl_atc <- c(
  "A10BA02",
  "A10BG03",
  "A10A",
  "A10BB"
)

incl <- incl_meds %>% 
  left_join(incl_cohort) %>% 
  filter(EVENT_AGE >= past_91 & EVENT_AGE < drug_init_age) %>% 
  filter(str_detect(ATC, paste0("^", incl_atc, collapse = "|"))) %>% 
  select(FINNGENID, ATC) %>% 
  distinct() 

incl_cohort <- incl_cohort[incl_cohort$FINNGENID %in% incl$FINNGENID, ]

##### Add C7 to SDF #####
# Excluded based on Inclusion - Rx for Metformin/Insulin/Piogliatzone/SU
df$c7 <- n_distinct(incl_cohort$FINNGENID)


#____________________________________________________________________
#### Inclusion: Taking either Pioglitazone, Metformin, Insulins or SU 90 days prior to init ####

# updating meds (longitudinal drug purchase file) and long_ep (longitudinal endpoint file)
incl_meds <- meds[meds$FINNGENID %in% incl_cohort$FINNGENID]
incl_longep <- incl_longep[incl_longep$FINNGENID %in% incl_cohort$FINNGENID]

# Age >= 50
incl1 <- incl_cohort %>% filter(drug_init_age >= 50) %>% select(FINNGENID)


# High CV Risk: 0 - 365 days prior to drug init
incl_diagnosis <- c(
  "I9_MI", 
  "I9_REVASC",
  "I9_ANGIO",
  "I9_CABG",
  "I9_STR_EXH",
  "I9_ATHSCLE",
  "I9_STENOSIS"
)

incl2 <- incl_longep %>% left_join(incl_cohort) %>% filter(ENDPOINT %in% incl_diagnosis & EVENT_AGE >= past_365 & EVENT_AGE <= drug_init_age) %>% select(FINNGENID) %>% distinct() 


incl_cohort <- rbind(incl1, incl2) %>% distinct() %>% left_join(cohort.ages)


##### Add C8 to SDF #####
# Exclusion based on Inclusion - Age 50 OR one major CV Risk factor in past 365 days (MI,  Revascularization, Coronary Angioplasty (PTCA), Bypass, Ischemic Stroke, Atherosclerosis, Amputation, Stenosis)
df$c8 <- n_distinct(incl_cohort$FINNGENID)




#____________________________________________________________________
#### EXCLUSION CRITERIA ####

# the Exclusions will be technically performed by time frames (e.g. all exclusion in 180 days window, all exclusions in 365 days window, etc.)

# updating meds (longitudinal drug purchase file) and long_ep (longitudinal endpoint file)
incl_meds <- meds[meds$FINNGENID %in% incl_cohort$FINNGENID]
incl_longep <- incl_longep[incl_longep$FINNGENID %in% incl_cohort$FINNGENID]


#### Exclusion: One of the following Endpoints in 365 days window #####

excl_diagnosis <- c(
  "T1D",
  "E4_HYPOGLYC"
)
excl <- incl_longep %>% left_join(incl_cohort) %>% filter(ENDPOINT %in% excl_diagnosis & EVENT_AGE >= past_365 & EVENT_AGE <= drug_init_age) %>% select(FINNGENID) %>% distinct() 
incl_cohort <- incl_cohort[incl_cohort$FINNGENID %!in% excl$FINNGENID, ]

##### Add C9 to SDF #####
# Exclusion based on Exclusion - T1D, Hypoglycaemia, Ketoacidosis
df$c9 <- n_distinct(incl_cohort$FINNGENID)



#____________________________________________________________________
#### Exclusion: One of the following Endpoints in 180 days window #####

# updating meds (longitudinal drug purchase file) and long_ep (longitudinal endpoint file)
incl_meds <- meds[meds$FINNGENID %in% incl_cohort$FINNGENID]
incl_longep <- incl_longep[incl_longep$FINNGENID %in% incl_cohort$FINNGENID]

excl_diagnosis <- c(
  "CHIRHEP_NAS",
  "XV_PREGNANCY_BIRTH",
  "Z21_CONTRACEPTIVE_MANAG"
)

excl <- incl_longep %>% left_join(incl_cohort) %>% filter(ENDPOINT %in% excl_diagnosis & EVENT_AGE >= past_180 & EVENT_AGE <= drug_init_age) %>% select(FINNGENID) %>% distinct() 
incl_cohort <- incl_cohort[incl_cohort$FINNGENID %!in% excl$FINNGENID, ]

##### Add C10 to SDF #####
# Excluded based on Exclusion - Pregnancy/Use of Contraceptive and Liver Cirrhosis
df$c10 <- n_distinct(incl_cohort$FINNGENID)



#____________________________________________________________________
#### Exclusion: One of the following Medication in 90 days window  #####

# updating meds (longitudinal drug purchase file) and long_ep (longitudinal endpoint file)
incl_meds <- meds[meds$FINNGENID %in% incl_cohort$FINNGENID]
incl_longep <- incl_longep[incl_longep$FINNGENID %in% incl_cohort$FINNGENID]

excl_atc <- c(
  "A10BH",
  "A10BJ",
  "A10BG01",
  "A10BG02",
  "A10BG04"
)

excl <- incl_meds %>% 
  left_join(incl_cohort) %>% 
  filter(EVENT_AGE >= past_91 & EVENT_AGE < drug_init_age) %>% 
  filter(str_detect(ATC, paste0("^", excl_atc, collapse = "|"))) %>% 
  distinct(FINNGENID, ATC) %>% 
  # since you exclude DPP4i, make sure to not exclude all Exposures (Sitagliptin [also a DPP4i])
  filter(ATC != "A10BH01") %>% 
  select(FINNGENID)
incl_cohort <- incl_cohort[incl_cohort$FINNGENID %!in% excl$FINNGENID, ]

##### Add C11 to SDF #####
# Excluded based on Exclusion - DPP4i, GLP-1, Thiazolidinediones (other than Pioglitazone)
df$c11 <- n_distinct(incl_cohort$FINNGENID)


table((cohort[cohort$FINNGENID %in% incl_cohort$FINNGENID, ])$ATC)



#____________________________________________________________________
#### FINAL COHORT ####

cohort_final <- cohort[cohort$FINNGENID %in% incl_cohort$FINNGENID, ]


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
# Excluded based on Exclusion - DPP4i, GLP-1, Thiazolidinediones (other than Pioglitazone)
df$c12 <- n_distinct(cohort.ppsmatching$FINNGENID)


table(cohort.ppsmatching$medication)

fwrite(cohort.ppsmatching, cohort_tec_elig_)
fwrite(df, cohort_tec_elig_flowchart_)

