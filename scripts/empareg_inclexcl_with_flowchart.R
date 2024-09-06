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
file1 <- endpoint_long_file_               # ENDPOINT LONGITUDINAL FILE PATH
file2 <- det_long_file_                    # MEDS LONGITUDINAL FILE PATH
file3 <- min_endpoint_file_                # ENDPOINT FILE PATH for SEX
file4 <- vnr_dict_                         # VNRO dictionary

long_ep <- fread(file1)
raw <- fread(file2)
meds <- raw %>% filter(SOURCE == "PURCH") %>% 
  # keep only relevant columns
  select(FINNGENID, EVENT_AGE, APPROX_EVENT_DAY, ATC = CODE1, VNRO = CODE3, PCKGNR = CODE4)   
sex <- fread(file3) %>% select(FINNGENID, SEX)

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
gc()


#____________________________________________________________________
#### Drugs of Interest ####

# define the medication groups of interest
trial_meds <- c("Exposure" = "A10BK03",
                "Reference" = "A10BH")

trial_meds <- as.data.frame(trial_meds) %>% rownames_to_column(var = "drugs") %>% rename(atcs = trial_meds)

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
df$c1 <- length(unique(cohort$FINNGENID)); gc()


# Exposure group 
exp <- cohort %>% 
  filter(ATC == "Exposure") %>%
  # choosing each IDs first purchase
  group_by(FINNGENID) %>% 
  arrange(EVENT_AGE) %>% 
  slice(1) %>% 
  ungroup()

# first purchase of exp (was the later drug to get MA)
min(exp$APPROX_EVENT_DAY) # 2015

# Reference group 
ref <- cohort[cohort$FINNGENID %!in% exp$FINNGENID] %>% 
  # removing IDs who started treatment before Empagliflozin came to market
  filter(APPROX_EVENT_DAY >= "2015-01-01") %>%
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
df$c2 <- n_distinct(cohort$FINNGENID); gc()





#____________________________________________________________________
#### Age Function ####
# FUNCTION: Add ages before drug initiation to check for meeting different inclusion/exclusion criteria 

AGES <- function(group){
  group.ages <- group %>% mutate(
    past_3650 = EVENT_AGE -10,
    past_1825 = EVENT_AGE - 5,
    past_730  = EVENT_AGE - 2,
    past_365  = EVENT_AGE - 1,
    past_180  = EVENT_AGE - 0.5,
    past_91   = EVENT_AGE - 0.25,
    past_61   = EVENT_AGE - 0.17,
    past_42   = EVENT_AGE - 0.115
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
df$c3 <- n_distinct(cohort$FINNGENID); gc()


# Adding relevant Ages for further inclusion/exclusion criteria:
cohort.ages <- AGES(cohort)


# Excluding IDs under Study Age
min_study_age <- 18
cohort.ages <- cohort.ages %>% filter(drug_init_age >= min_study_age)

##### Add C4 to SDF #####
# Excluded based on Inclusion Age >= min Study Age (18y)
df$c4 <- n_distinct(cohort.ages$FINNGENID); gc()



#____________________________________________________________________
#### Inclusion Diagnosis: Type 2 Diabetes ####
# (measured up to 180 days prior to drug init)

incl_diagnosis <- "T2D"

incl <- (long_ep %>% filter(ENDPOINT == incl_diagnosis))[(long_ep %>% filter(ENDPOINT == incl_diagnosis))$FINNGENID %in% cohort.ages$FINNGENID] %>% left_join(cohort.ages) %>% filter(EVENT_AGE >= past_180 & EVENT_AGE <= drug_init_age)
incl_cohort <- cohort.ages[cohort.ages$FINNGENID %in% incl$FINNGENID, ]

##### Add C5 to SDF #####
# Excluded based on Inclusion T2D
df$c5 <- n_distinct(incl_cohort$FINNGENID); gc()

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
df$c6 <- n_distinct(incl_cohort$FINNGENID); gc()



#____________________________________________________________________
#### Inclusion: At least having one of the following endpoints in different time periods before drug init ####

# updating meds (longitudinal drug purchase file) and long_ep (longitudinal endpoint file)
incl_meds <- meds[meds$FINNGENID %in% incl_cohort$FINNGENID]
incl_longep <- incl_longep[incl_longep$FINNGENID %in% incl_cohort$FINNGENID]

# 61 - 3650 days prior to drug init
incl_diagnosis <- c(
  "I9_ANGINA", 
  "I9_UAP",
  "C_STROKE", 
  "I9_MI"
)

incl <- incl_longep %>% left_join(incl_cohort) %>% filter(ENDPOINT %in% incl_diagnosis & EVENT_AGE >= past_3650 & EVENT_AGE <= past_61) %>% select(FINNGENID) %>% distinct() 

# 0 - 3650 days prior to drug init
incl_diagnosis <- c(
  "I9_STENOSIS", 
  "I9_REVASC"
)

incl2 <- incl_longep %>% left_join(incl_cohort) %>% filter(ENDPOINT %in% incl_diagnosis & EVENT_AGE >= past_3650 & EVENT_AGE <= drug_init_age) %>% select(FINNGENID) %>% distinct() 

# 0 - 180 days prior to drug init
incl_diagnosis <- c(
  "I9_PAD",
  "IX_CIRCULATORY",
  "I9_ANGIO",
  "I9_CABG",
  "I9_DOAAC"
)

incl3 <- incl_longep %>% left_join(incl_cohort) %>% filter(ENDPOINT %in% incl_diagnosis & EVENT_AGE >= past_180 & EVENT_AGE <= drug_init_age) %>% select(FINNGENID) %>% distinct() 

# 0 - 365 days prior to drug init
incl_atc <- c(
  "C09AA",
  "C09AC",
  "C08",
  "C01DA"
)

incl4 <- incl_meds %>% 
  left_join(incl_cohort) %>% 
  filter(EVENT_AGE >= past_180 & EVENT_AGE <= drug_init_age) %>% 
  filter(str_detect(ATC, paste0("^", incl_atc, collapse = "|"))) %>% 
  select(FINNGENID) %>% 
  distinct() 


incl_cohort <- rbind(incl, incl2, incl3, incl4) %>% distinct() %>% left_join(cohort.ages)


##### Add C7 to SDF #####
# Excluded based on high risk of cardio vascular event
df$c7 <- n_distinct(incl_cohort$FINNGENID); gc()





#____________________________________________________________________
#### EXCLUSION CRITERIA ####

# the Exclusions will be technically performed by time frames (e.g. all exclusion in 180 days window, all exclusions in 365 days window, etc.)

# updating meds (longitudinal drug purchase file) and long_ep (longitudinal endpoint file)
incl_meds <- meds[meds$FINNGENID %in% incl_cohort$FINNGENID]
incl_longep <- incl_longep[incl_longep$FINNGENID %in% incl_cohort$FINNGENID]


#### Exclusion: One of the following Endpoints in 180 days window #####

excl_diagnosis <- c(
  "E4_OBESITYXTRMALV",
  "K11_LIVER",
  "D3_BLOOD_IMMUN",
  "KRA_PSY_ALCOH",
  "KRA_PSY_SUBSTANCE",
  "XV_PREGNANCY_BIRTH",
  "Z21_CONTRACEPTIVE_MANAG"
)
excl <- incl_longep %>% left_join(incl_cohort) %>% filter(ENDPOINT %in% excl_diagnosis & EVENT_AGE >= past_180 & EVENT_AGE <= drug_init_age) %>% select(FINNGENID) %>% distinct() 
incl_cohort <- incl_cohort[incl_cohort$FINNGENID %!in% excl$FINNGENID, ]

##### Add C8 to SDF #####
# Excluded based on exclusion criteria for past 180 days: Obesity, Liver disease, Blood disease, Substance/Alcohol abuse, Pregnancy
df$c8 <- n_distinct(incl_cohort$FINNGENID); gc()



#____________________________________________________________________
#### Exclusion: One of the following Endpoints in 1825 days window #####

# updating meds (longitudinal drug purchase file) and long_ep (longitudinal endpoint file)
incl_meds <- meds[meds$FINNGENID %in% incl_cohort$FINNGENID]
incl_longep <- incl_longep[incl_longep$FINNGENID %in% incl_cohort$FINNGENID]

excl_diagnosis <- c(
  "C3_CANCER"
)

excl <- incl_longep %>% left_join(incl_cohort) %>% filter(ENDPOINT %in% excl_diagnosis & EVENT_AGE >= past_1825 & EVENT_AGE <= drug_init_age) %>% select(FINNGENID) %>% distinct() 
incl_cohort <- incl_cohort[incl_cohort$FINNGENID %!in% excl$FINNGENID, ]

##### Add C9 to SDF #####
# Excluded based on Exclusion History of Malignant Neoplasm past 5 years
df$c9 <- n_distinct(incl_cohort$FINNGENID); gc()



#____________________________________________________________________
#### Exclusion: One of the following Endpoints in 730 days window #####

# updating meds (longitudinal drug purchase file) and long_ep (longitudinal endpoint file)
incl_meds <- meds[meds$FINNGENID %in% incl_cohort$FINNGENID]
incl_longep <- incl_longep[incl_longep$FINNGENID %in% incl_cohort$FINNGENID]

excl_diagnosis <- c(
  "K11_MALABSORB"
)

excl <- incl_longep %>% left_join(incl_cohort) %>% filter(ENDPOINT %in% excl_diagnosis & EVENT_AGE >= past_730 & EVENT_AGE <= drug_init_age) %>% select(FINNGENID) %>% distinct() 
incl_cohort <- incl_cohort[incl_cohort$FINNGENID %!in% excl$FINNGENID, ]

##### Add C10 to SDF #####
# Excluded based on Exclusion Chronic Malabsorption
df$c10 <- n_distinct(incl_cohort$FINNGENID); gc()



#____________________________________________________________________
#### Exclusion: One of the following Medication in 90 days window  #####

# updating meds (longitudinal drug purchase file) and long_ep (longitudinal endpoint file)
incl_meds <- meds[meds$FINNGENID %in% incl_cohort$FINNGENID]
incl_longep <- incl_longep[incl_longep$FINNGENID %in% incl_cohort$FINNGENID]

excl_atc <- c(
  "A08"
)

excl <- incl_meds %>% 
  left_join(incl_cohort) %>% 
  filter(EVENT_AGE >= past_91 & EVENT_AGE <= drug_init_age) %>% 
  filter(str_detect(ATC, paste0("^", excl_atc, collapse = "|"))) %>% 
  select(FINNGENID) %>% 
  distinct() 
incl_cohort <- incl_cohort[incl_cohort$FINNGENID %!in% excl$FINNGENID, ]

##### Add C11 to SDF #####
# Excluded based on Anti-obseity drugs
df$c11 <- n_distinct(incl_cohort$FINNGENID); gc()



#____________________________________________________________________
#### Exclusion: One of the following Endpoints in 61 days window #####

# updating meds (longitudinal drug purchase file) and long_ep (longitudinal endpoint file)
incl_meds <- meds[meds$FINNGENID %in% incl_cohort$FINNGENID]
incl_longep <- incl_longep[incl_longep$FINNGENID %in% incl_cohort$FINNGENID]

excl_diagnosis <- c(
  "I9_TIA",
  "C_STROKE",
  "I9_REVASC"
)

excl <- incl_longep %>% left_join(incl_cohort) %>% filter(ENDPOINT %in% excl_diagnosis & EVENT_AGE >= past_61 & EVENT_AGE <= drug_init_age) %>% select(FINNGENID) %>% distinct() 
incl_cohort <- incl_cohort[incl_cohort$FINNGENID %!in% excl$FINNGENID, ]

##### Add C12 to SDF #####
# Excluded based on Acute coronary syndrome, stroke, or transient ischemic attack within 2 months
df$c12 <- n_distinct(incl_cohort$FINNGENID); gc()



#____________________________________________________________________
#### Exclusion: One of the following Medication in 42 days window  #####

# updating meds (longitudinal drug purchase file) and long_ep (longitudinal endpoint file)
incl_meds <- meds[meds$FINNGENID %in% incl_cohort$FINNGENID]
incl_longep <- incl_longep[incl_longep$FINNGENID %in% incl_cohort$FINNGENID]

excl_atc <- c(
  "H02"
)

excl <- incl_meds %>% 
  left_join(incl_cohort) %>% 
  filter(EVENT_AGE >= past_42 & EVENT_AGE <= drug_init_age) %>% 
  filter(str_detect(ATC, paste0("^", excl_atc, collapse = "|"))) %>% 
  select(FINNGENID) %>% 
  distinct() 
incl_cohort <- incl_cohort[incl_cohort$FINNGENID %!in% excl$FINNGENID, ]

##### Add C13 to SDF #####
# Excluded based on Exclusion - Use of corticosteroids
df$c13 <- n_distinct(incl_cohort$FINNGENID); gc()





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
table(score_data_180prior$CCI); gc()


#____________________________________________________________________
#### Prepare Data for Matching and Survival Analysis: ####

cohort.ppsmatching <- final_cohort %>% left_join(sex) %>% left_join(cohort.ages %>% select(FINNGENID, drug_init_age)) %>% left_join(score_data_180prior, by = c("FINNGENID" = "ID"))
colnames(cohort.ppsmatching)[which(names(cohort.ppsmatching) == "drug_init_age")] <- "age"
cohort.ppsmatching[cohort.ppsmatching == "Exposure"] <- 1  # Exposure
cohort.ppsmatching$medication[cohort.ppsmatching$medication != "1"] <- 0  # Reference
cohort.ppsmatching$medication <- as.numeric(cohort.ppsmatching$medication)
cohort.ppsmatching$SEX <- as.factor(cohort.ppsmatching$SEX)
table(cohort.ppsmatching$medication); gc()



#____________________________________________________________________
#### Exclude IDs with CCI >= 10 ####
cohort.ppsmatching <- cohort.ppsmatching %>% filter(CCI <= 9)


table(cohort.ppsmatching$medication)

##### Add C14 to SDF #####
# CCI >= 10
df$c14 <- n_distinct(cohort.ppsmatching$FINNGENID); gc()

fwrite(cohort.ppsmatching, cohort_emp_elig_)
fwrite(df, cohort_emp_elig_flowchart_)









