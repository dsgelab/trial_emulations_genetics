# Clean Environment -------------------------------------------------------

rm(list = ls())
gc()


# Loading files and packages ----------------------------------------------

#### Read all required packages ###
library(tidyr)
library(dplyr)
library(tibble)
library(data.table)
library(stringr)
source("/home/ivm/config_te.R")


#### Read all required files ###
# Define file paths 
file01 <- cohort_ari_elig_
file1 <- endpoint_long_file_  # ENDPOINT LONGITUDINAL FILE PATH
file2 <- det_long_file_det_long_file_  # MEDS LONGITUDINAL FILE PATH
file3 <- min_endpoint_file_            # ENDPOINT FILE PATH for SEX and Birthyear
file4 <- vnr_dict_                     # VNRO dict 


# Read file paths
trial_cohort <- fread(file01)  

long_ep <- fread(file1) %>% filter(FINNGENID %in% trial_cohort$FINNGENID); gc()
raw <- fread(file2) %>% filter(FINNGENID %in% trial_cohort$FINNGENID); gc()
vnro_dict <- fread(file4)



#### Preprocess some files ####
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


# Defining the opposite of %in% 
'%!in%' <- function(x,y)!('%in%'(x,y))

# COVARIATES
covariates <- data.frame(covariates = c(
  "RX_ANTIHYP",
  "I9_CVD_HARD",
  "I9_CHD",
  'I9_HEARTFAIL_NS',
  "DM_NEPHROPATHY",
  "DM_NEUROPATHY",
  "DM_NEOVASCULAR_GLAUCOMA",
  "E4_HYPOGLYC",
  "N14_CHRONKIDNEYDIS",
  "J10_ASTHMA_INCLAVO",
  "J10_COPD",
  'J10_RESPIRATORY',
  'SMOKING_DEPEND',
  "RX_STATIN",
  'RX_OPIOIDS',
  'RX_NSAID',
  'RX_BZD',
  'RX_NO5C',
  'GLP1ANA',
  'DIAB_MED',
  'E4_DIABETES',
  "I9_E_DRONE",
  "I9_TIA"
))




# Covariates Analysis -----------------------------------------------------

### Time horizon for covariates assessment: 1 year prior to drug initiation

#### Endpoints ####

# Extract Endpoints from cohort
cohort_eps <- trial_cohort %>% 
  left_join(long_ep) %>%
  # keep only endpoints before trial init and no younger than 1 years before
  filter(EVENT_AGE <= age & EVENT_AGE >= age - 1) %>% 
  # remove some sources for consistency
  filter(EVENT_TYPE %!in% c("CANCER", "ERIK_OPER", "OPER")) %>%
  # take only the last endpoint occurrence per ID
  group_by(FINNGENID, ENDPOINT) %>% 
  arrange(EVENT_AGE) %>% 
  slice(n()) %>% 
  ungroup()


# Occurrence of covariate endpoints in cohorts
pop_desc <- cohort_eps %>%
  # filter for the endpoints of interest
  filter(ENDPOINT %in% covariates$covariates) %>%
  # group by patient ID and endpoint
  group_by(FINNGENID, ENDPOINT) %>%
  # summarize the age and time of endpoint for each endpoint of interest (last time if it occurred multiple times)
  summarize(endpoint = if_else(any(ENDPOINT %in% covariates$covariates), 1, 0)) %>% 
  # pivot the endpoint column to create new columns for each endpoint
  pivot_wider(names_from = ENDPOINT, values_from = endpoint) #, names_prefix = "endpoint_")

results <- data.frame(
  trial_cohort[,1:2] 
)

# Resulting dataframe with occurrence of each endpoint covariate for each ID
results <- results %>% left_join(pop_desc)
results[is.na(results)] <- 0


#### Medication ####

# Extract Medication from cohort
cohort_meds <- trial_cohort %>% left_join(meds) %>%
  # keep only endpoints before trial init
  filter(EVENT_AGE <= age & EVENT_AGE >= age - 1) %>%  
  # keep only relevant columns
  select(FINNGENID:age, EVENT_AGE, ATC) %>% 
  # keep only the 5-letter code of ATC
  mutate(ATC = substr(ATC, 1, 5)) %>%
  # take only the last medication occurrence per ID
  group_by(FINNGENID, ATC) %>% 
  arrange(EVENT_AGE) %>% 
  slice(n()) %>% 
  ungroup()


# Occurrence of covariate medication in cohorts
pop_desc_meds <- cohort_meds %>%
  # filter for the endpoints of interest
  filter(str_detect(ATC, paste0("^", covariates$covariates, collapse = "|"))) %>%
  # replace ATC code with the ATC group string that is defined as covariate (e.g. C03CA -> C03 and C10AA -> C10AA)
  mutate(
    ATC = case_when(
      str_detect(ATC, paste0("^", covariates$covariates, collapse = "|")) ~ covariates$covariates[match(str_extract(ATC, paste0("^", covariates$covariates, collapse = "|")), covariates$covariates)],
      T ~ ATC
    )
  ) %>% 
  distinct(FINNGENID, ATC, .keep_all = T) %>% 
  # add column with value 1 for the pivot_wider later -> 1 = ID had the medication
  mutate(column = 1) %>% 
  # remove EVENT_AGE for pivot_wider
  select(-EVENT_AGE, -SEX, -age) %>% 
  pivot_wider(names_from = ATC, values_from = column)


# Adding medication covariates to the resulting dataframe
results <- results %>% left_join(pop_desc_meds)
results[is.na(results)] <- 0

# Factorize
results <- results %>% mutate_at(vars(3:ncol(results)), as.factor)



# Save Outputs ------------------------------------------------------------

fwrite(results, ari_ps_cov_)
