rm(list = ls())
gc()

#### Loading files and packages ####
# Load all required packages:
library(tidyr)
library(dplyr)
library(lubridate)
library(stringr)
library(forcats)
library(ggplot2)
library(data.table)
library(R.utils)
library(ICCI)
library(MatchIt)
library(survival)
library(ggfortify)
library(WeightIt)
library(cobalt)
source("/home/ivm/config_te.R")


# Defining the opposite of %in% to exclude rows that meet exclusion criteria
'%!in%' <- function(x,y)!('%in%'(x,y)) 

#### Read all required files ###
# Define file paths 
file01 <- cohort_ari_elig_
file1 <- endpoint_long_file_           # ENDPOINT LONGITUDINAL FILE PATH
file2 <- det_long_file_det_long_file_  # MEDS LONGITUDINAL FILE PATH
file3 <- min_endpoint_file_            # ENDPOINT FILE PATH for SEX and Birthyear
file4 <- vnr_dict_                     # VNRO dict                            


# Read file paths
trial_cohort <- fread(file01)  

long_ep <- fread(file1) %>% filter(FINNGENID %in% trial_cohort$FINNGENID); gc()
raw <- fread(file2) %>% filter(FINNGENID %in% trial_cohort$FINNGENID); gc()
mp <- fread(file3) %>% filter(FINNGENID %in% trial_cohort$FINNGENID); gc()
vnro_dict <- fread(file4)


##### Preprocess some files #####
# Meds
meds <- raw %>%
  # keep only meds (-> purchase registry) 
  filter(SOURCE == "PURCH") %>% 
  # keep only relevant columns
  select(FINNGENID, EVENT_AGE, APPROX_EVENT_DAY, ATC = CODE1, VNRO = CODE3, PCKGNR = CODE4) 

# VNRO Dictionary
vnro_dict <- vnro_dict %>% rename(VNRO = vnr, DOSE = vahvuus_num, PKOKO = pkoko_num) %>% select(VNRO, ATC, PKOKO, DOSE)
vnro_dict$VNRO <- as.character(vnro_dict$VNRO)
vnro_dict$VNRO <- str_pad(vnro_dict$VNRO, 6, pad = "0")

# Keep only birth year
mp <- mp %>% mutate(birthyear = year(APPROX_BIRTH_DATE)) %>% select(FINNGENID, birthyear)



#_____________________________________________
#### Add birthyear #####

trial_cohort <- trial_cohort %>% left_join(mp)



#_____________________________________________
#### Study follow-up #####

# Preprocess Trial cohort
trial_cohort <- trial_cohort %>% rename(group = medication)
trial_cohort$group <- as.character(trial_cohort$group)

# Only Cohort Endpoints after initiation of study
cohort_endpoints <- long_ep %>% filter(long_ep$FINNGENID %in% trial_cohort$FINNGENID) 
cohort_endpoints <- merge(x = cohort_endpoints,  y = trial_cohort[ , c("FINNGENID", "age")], all.x = T) %>% filter(!(EVENT_AGE < age))


##Stroke / Systemic Embolism
po <- "I9_ARTEMBTHR|I9_OTHINTRACRA|I9_STR_SAH|I9_PULMEMB"
cohort_po <- cohort_endpoints %>% filter(grepl(po, ENDPOINT)) %>% arrange(FINNGENID, EVENT_AGE)
cohort_po <- cohort_po[match(unique(cohort_po$FINNGENID), cohort_po$FINNGENID)]
cohort_po$primary_outcome <- 1
setnames(cohort_po, "EVENT_AGE", "primary_outcome_age")

# Join primary outcomes with matches
trial_cohort <- merge(x = trial_cohort, y = cohort_po[ , c("FINNGENID", "primary_outcome", "primary_outcome_age")], all.x = T)
trial_cohort$primary_outcome[is.na(trial_cohort$primary_outcome)] <- 0

# Time to event for primary outcome
trial_cohort <- trial_cohort %>% mutate(survival_years_po = primary_outcome_age - age) #%>% filter(!is.na(survival_years))
trial_cohort %>% group_by(group) %>% filter(!is.na(survival_years_po)) %>% summarise(mean_followup = mean(survival_years_po))


##### Death ####
# Death as end of study follow-up:
cohort_death <- cohort_endpoints %>% filter(ENDPOINT == "DEATH")
cohort_death$death <- 1
setnames(cohort_death, "EVENT_AGE", "death_age")

# Join death with matches
trial_cohort <- merge(x = trial_cohort, y = cohort_death[ , c("FINNGENID", "death", "death_age")], all.x = T)
trial_cohort$death[is.na(trial_cohort$death)] <- 0

# Time to event for death
trial_cohort <- trial_cohort %>% mutate(survival_years_death = death_age - age)


##### End of registry ####
# End of registry day
eor <- max(long_ep$APPROX_EVENT_DAY)  # 2023-05-10

# get age for end of the registry (2023-05-10):
age_endofregistry <- trial_cohort %>% select(FINNGENID, age) %>% 
  left_join(meds) %>% 
  filter(EVENT_AGE == age) %>% 
  select("FINNGENID", "age", "APPROX_EVENT_DAY") %>% 
  # calculate age at end of registry
  mutate(
    age_endofregistry = age + round(as.numeric(difftime(eor, APPROX_EVENT_DAY, unit = "weeks"))/52.25, 2)
  ) %>% 
  # select only FINNGENID and age at end of registry
  select(FINNGENID, age_endofregistry) %>% 
  # keep only distinct rows
  distinct(FINNGENID, .keep_all = T) 

# Join endofregistry with matches
trial_cohort <- trial_cohort %>% left_join(age_endofregistry)
trial_cohort <- trial_cohort %>% mutate(survival_years_endofregistry = age_endofregistry - age)


##### Discontinuation of Drug #####

## FIRST DISC PART (UNTIL vnros_to_remove) IS ONLY TO DETERMINE THE VNROS TO REMOVE - CAN BE DELETED, ONCE THERE IS A COMPLETE VNRO MAPPING
disc_tmp   <- trial_cohort %>% select(FINNGENID, group, age) %>% 
  left_join(meds) %>% 
  # only drugs after init
  filter(EVENT_AGE >= age) %>% 
  # filter for Apixaban and Warfarin
  filter(ATC == "B01AF02" | ATC == "B01AA03") %>%
  # keep only the medication of group assignment in case someone switches later
  filter((group == 1 & ATC == "B01AF02") | (group == 0 & ATC == "B01AA03")) %>% 
  # REMOVE VNRO CODES THAT DON'T HAVE ANY MAPPING CURRENTLY:
  left_join(vnro_dict)

# Remove VNROs that are unmapped
vnros_to_remove <- disc_tmp %>% filter(is.na(PKOKO)) %>% ungroup() %>% distinct(VNRO) %>% pull(VNRO)

disc <- trial_cohort %>% select(FINNGENID, group, age) %>% 
  left_join(meds) %>% 
  # only drugs after init
  filter(EVENT_AGE >= age) %>% 
  # filter for Apixaban and Warfarin
  filter(ATC == "B01AF02" | ATC == "B01AA03") %>%
  # keep only the medication of group assignment in case someone switches later
  filter((group == 1 & ATC == "B01AF02") | (group == 0 & ATC == "B01AA03"))  %>%  
  # REMOVE VNRO CODES THAT DON'T HAVE ANY MAPPING CURRENTLY:
  filter(!(VNRO %in% vnros_to_remove)) %>% 
  # select last purchase
  group_by(FINNGENID) %>% 
  arrange(EVENT_AGE) %>% 
  slice(n()) %>%
  # remove ATC column as it is redundant
  select(-ATC) %>% 
  rename(last_purch_age = EVENT_AGE)

disc <- disc %>% 
  # join the vnro dictionary
  left_join(vnro_dict) %>% 
  # calculate the disc_age (estimated time of last pill from package size and number of packages)
  # and add survival years until discontinuation
  mutate(disc_age = (PKOKO*PCKGNR)/365.25 + last_purch_age,
         survival_years_disc = disc_age - age) %>% 
  select(FINNGENID, disc_age, survival_years_disc)


disc %>% filter(if_any(everything(), is.na)) # check that survival years disc was set to 0

# Join last purchase with trial_cohort and discard IDs on unavailable VNROs from the start
trial_cohort <- trial_cohort %>% left_join(disc) %>% filter(!is.na(survival_years_disc))



##### Switch of Drug / or both meds at the same time #####

switch <- trial_cohort %>% select(FINNGENID, group, age) %>% 
  left_join(meds) %>% 
  # only drugs after init
  filter(EVENT_AGE >= age) %>% 
  # filter for Empagliflozine and DPP4i
  filter(str_detect(ATC, "^B01AF02|^B01AA03")) %>%
  # change the medication assignment -> SWITCHES/IDs with two meds    
  filter((group == 1 & ATC == "B01AA03") | (group == 0 & grepl("^B01AF02", ATC))) %>% 
  # select last purchase
  group_by(FINNGENID) %>% 
  arrange(EVENT_AGE) %>% 
  slice(1) %>%  # take first purchase  
  # remove ATC column as it is redundant
  select(-ATC) %>% 
  # EVENT_AGE is the switch age
  rename(switch_age = EVENT_AGE)

switch <- switch %>% 
  # join the vnro dictionary
  left_join(vnro_dict) %>% 
  # and add survival years until switch
  mutate(survival_years_switch = switch_age - age) %>% 
  select(FINNGENID, switch_age, survival_years_switch)

# Join switch with trial_cohort
trial_cohort <- trial_cohort %>% left_join(switch)


##### End of trial ####
# end of trial up to 5 years
trial_cohort <- trial_cohort %>% 
  mutate(
    age_endoftrial = age + 2,
    survival_years_endoftrial = 2
  )


trial_cohort %>% select(tidyselect::starts_with("survival"))

##### Survival time (as incident which has occured the first) #####
# Keep only the value with the lowest number in row
trial_cohort <- trial_cohort %>% 
  mutate(time = round(pmin(survival_years_death, survival_years_po, survival_years_endofregistry, survival_years_disc, survival_years_switch, survival_years_endoftrial, na.rm = T), 2))


##### Setting primary outcome (event) to 0 if it's censored

trial_cohort %>% group_by(group) %>% summarize(count_po = sum(primary_outcome == 1))
trial_cohort <- trial_cohort %>% mutate(primary_outcome = if_else(abs(time - round(survival_years_po, 2)) > 0.1 & !(time > 2)  | is.na(survival_years_po), 0, primary_outcome))
# primary outcomes per group AFTER survival analysis adjustment 
trial_cohort %>% group_by(group) %>% summarize(count_1 = sum(primary_outcome == 1))


# save trial_cohort
fwrite(trial_cohort, ari_surv_)

