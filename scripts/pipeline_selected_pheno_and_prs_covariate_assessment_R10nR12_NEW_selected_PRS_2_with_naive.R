# Packages ----------------------------------------------------------------

# Load all required packages:

library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(purrr)
library(tibble)
library(stringr)
library(forcats)
library(data.table)
library(MatchIt)
library(survival)
library(ggfortify)
library(cobalt)
library(gridExtra)
library(ggrepel)
library(cowplot)
library(patchwork)
library(grDevices)
source("/home/ivm/config_te.R")


# Functions ---------------------------------------------------------------

# PRS Preprocessing
prs_prep <- function(prs) {
  name <- deparse(substitute(prs))
  prs <- fread(prs)
  prs[[name]] <- scale(prs$SCORE1_AVG)
  prs <- prs %>% dplyr::rename("FINNGENID" = "IID") %>% select(FINNGENID, name)
  
  return(prs)
}


# Join PRSs to trials
add_prs <- function(trial) {
  
  prs_files <- Filter(function(x) is.data.frame(get(x, envir = parent.frame())), grep("^prs_", ls(envir = parent.frame()), value = TRUE))
  
  for (prs_name in prs_files) {
    prs <- get(prs_name)
    trial <- trial %>% left_join(prs)
  }
  return(trial)
}


# PS Matching
ps_matching <- function(data, output_file = "psmatches_rctd", col_start_cov = 10000, col_end_cov = 10000, caliper = 0.1) {
  # Get the names of all covariate columns (na.omit is in order to remove na as a result of no additional covariates specified)
  covariate_names <- na.omit(c("age", "SEX", "CCI", "birthyear", "age*birthyear", names(data)[col_start_cov:col_end_cov]))
  
  # Create a formula with covariates
  formula <- as.formula(paste("group ~", paste(covariate_names, collapse = " + ")))
  
  m.out <- matchit(formula,                                                    
                   data = data,
                   method = "nearest", 
                   distance = "logit", 
                   caliper = caliper)     
  
  # Summary of covariate balance before and after matching
  matching_summary <- summary(m.out)
  matching_summary_plot <- plot(summary(m.out))
  
  # PS of all the matches
  ps <- match.data(m.out) %>% select(FINNGENID, group, distance)
  # Plot the ps distribution of unmatched cohort
  ggplot(ps, aes(x = distance, fill = as.factor(group))) + 
    geom_density(alpha = 0.5)
  
  # get matches 
  ps_matches <- MatchIt::get_matches(m.out)
  ggplot(ps_matches, aes(x = distance, fill = group)) + 
    geom_density(alpha = 0.5)
  ps_matches <- ps_matches %>% select(intersect(names(data), names(ps_matches)), distance)
  assign(output_file, ps_matches, envir = .GlobalEnv)
  
  matching_summary
  matching_summary_plot
  
}


# PRS Analysis
gen_anal <- function(data, trial, type){
  # All PRS columns
  prs_columns <- grep("^prs_", names(data), value = TRUE)
  
  df <- data.frame(matrix(nrow=0, ncol=7))
  names(df) <- c("trial", "type", "prs_name", "diff", "ci_lower", "ci_upper", "Pval")
  
  for (prs in prs_columns){
    print(prs)
    average_diff <- glm(group ~ get(prs), family=binomial(), data = data)
    
    diff <- summary(average_diff)$coefficients[, "Estimate"][2]
    CI_lower <- confint(average_diff)[,"2.5 %"][2]
    CI_upper <- confint(average_diff)[,"97.5 %"][2]
    Pval <- summary(average_diff)$coefficients[, "Pr(>|z|)"][2]
    
    v <- data.frame(trial, type, prs, diff, CI_lower, CI_upper, Pval) 
    names(v) <- names(df)
    df <- rbind(df, v) %>% tibble() 
    
    
  }
  return(df)
  
}


# Assess Covariate Balance - aSMD
cov_balance <- function(before = emp_pur, name = "Empareg PRS Matching", cov_start = 110, cov_end = 187){
  
  # Balance Statistics of covariates of interest
  smd_before <- col_w_smd(before %>% select(cov_start:cov_end), treat = before$group)
  
  # Absolute SMDs of covariates of interest - unadjusted and PRS adjsuted 
  balance_stats <- data.frame(
    covariates = names(smd_before),
    unadjusted = abs(unname(smd_before))
  )
  
  # Reshape the dataframe for plotting
  balance_stats_longer <- pivot_longer(balance_stats, cols = c("unadjusted"), names_to = "metric")
  
  
  p <- ggplot(balance_stats_longer, aes(x = value, y = covariates, color = metric)) +
    geom_point(size = 1) +
    labs(x = "Absolute Standardized Mean Difference", y = "Covariates", title = paste0("Love Plot: ", name)) +
    scale_color_manual(values = c("unadjusted" = "blue")) +
    geom_vline(xintercept = 0.1, color = "black") +
    geom_vline(xintercept = 0.05, color = "black", linetype = "dashed") +
    theme_minimal()
  
  p
  
}


# One metric assessment - aSMD
one_metric_assess <- function(before = emp_pur, cov_start = 110, cov_end = 187, covariates = everything()){
  
  # Balance Statistics of covariates of interest
  smd_before <- col_w_smd(before %>% select(cov_start:cov_end), treat = before$group)
 
  # Absolute SMDs of covariates of interest - unadjusted and PRS adjsuted 
  balance_stats <- data.frame(
    unadjusted = smd_before
  ) %>% 
    mutate(unbalanced = ifelse(unadjusted > 0.1 | unadjusted < - 0.1, 1, 0)) %>%
    # arrange by largest SMD (regardless of direction)
    arrange(desc(abs(unadjusted))) %>% 
    filter(str_detect(paste0("^", covariates, collapse = "|"), rownames(.))) 
  
  return(balance_stats)
}


# One metric assessment for PRS - aSMD (removing criterium of 0.1)
one_metric_assess_prs <- function(before = emp_pur, cov_start = 110, cov_end = 187, covariates = everything()){
  
  # Balance Statistics of covariates of interest
  smd_before <- col_w_smd(before %>% select(cov_start:cov_end), treat = before$group)
  
  # Absolute SMDs of covariates of interest - unadjusted and PRS adjsuted 
  balance_stats <- data.frame(
    unadjusted = smd_before
  ) %>% 
    # arrange by largest SMD (regardless of direction)
    mutate(unbalanced = ifelse(unadjusted > 0.1 | unadjusted < - 0.1, 1, 0)) %>% 
    arrange(desc(abs(unadjusted))) %>% 
    filter(str_detect(paste0("^", covariates, collapse = "|"), rownames(.)))
  
  return(balance_stats)
}


# Covariate Analysis for specific covariates (stored in a string -> cov_columns) and p_value
cov_anal <- function(data, trial, type, cov_columns){
  print(cov_columns)
  df <- data.frame(matrix(nrow=0, ncol=7))
  names(df) <- c("trial", "type", "prs_name", "diff", "ci_lower", "ci_upper", "Pval")
  
  for (cov in cov_columns){
    print(cov)
    average_diff <- glm(group ~ get(cov), family=binomial(), data = data)
    
    diff <- summary(average_diff)$coefficients[, "Estimate"][2]
    CI_lower <- confint(average_diff)[,"2.5 %"][2]
    CI_upper <- confint(average_diff)[,"97.5 %"][2]
    Pval <- summary(average_diff)$coefficients[, "Pr(>|z|)"][2]
    
    
    v <- data.frame(trial, type, cov, diff, CI_lower, CI_upper, Pval) 
    names(v) <- names(df)
    df <- rbind(df, v) %>% tibble() 
    
    
  }
  return(df)
  
}


# Covariate Analysis for XXXXXXXX
cov_anal_2 <- function(data, trial, type, cov_columns, remove_cov = NULL){
  covariates <- data %>% select(cov_columns, -remove_cov) %>% names() %>% c()
  df <- data.frame(matrix(nrow=0, ncol=7))
  names(df) <- c("trial", "type", "prs_name", "diff", "ci_lower", "ci_upper", "Pval")
  
  for (cov in covariates){
    print(cov)
    average_diff <- glm(group ~ get(cov), family=binomial(), data = data)
    
    diff <- summary(average_diff)$coefficients[, "Estimate"][2]
    CI_lower <- confint(average_diff)[,"2.5 %"][2]
    CI_upper <- confint(average_diff)[,"97.5 %"][2]
    Pval <- summary(average_diff)$coefficients[, "Pr(>|z|)"][2]
    
    
    v <- data.frame(trial, type, cov, diff, CI_lower, CI_upper, Pval) 
    names(v) <- names(df)
    df <- rbind(df, v) %>% tibble() 
    
    
  }
  return(df)
  
}



# Paths -------------------------------------------------------------------

# File paths
# Trials

file1 <- emp_surv_            # Empareg
file2 <- ari_surv_            # Aristotle
file4 <- roc_surv_            # Rocket
file5 <- tec_surv_            # Tecos

# Covariates 
file1_cov <- emp_ps_cov_                         
file2_cov <- ari_ps_cov_
file4_cov <- roc_ps_cov_
file5_cov <- tec_ps_cov_


# Naive Trials
file1_naive <- emp_naive_
file5_naive <- tec_naive_
file4_naive <- roc_naive_
file2_naive <-  ari_naive_


# Reading trials ----------------------------------------------------------

# Read file paths
emp_cohort <- fread(file1)  
emp_cov <- fread(file1_cov)
emp_cohort_naive <- fread(file1_naive)

ari_cohort <- fread(file2)  
ari_cov <- fread(file2_cov)
ari_cohort_naive <- fread(file2_naive)

roc_cohort <- fread(file4)  
roc_cov <- fread(file4_cov)
roc_cohort_naive <- fread(file4_naive)

tec_cohort <- fread(file5)  
tec_cov <- fread(file5_cov)
tec_cohort_naive <- fread(file5_naive)



#_____________________________________________
#### Add covariates #####

##### Covariate #####

emp_cohort <- emp_cohort %>% 
  left_join(emp_cov, by = "FINNGENID") %>% 
  select(-medication) %>% 
  # make them a factor
  mutate_at(vars(22:ncol(.)), as.factor)
emp_cohort <- emp_cohort %>% select(-"GLP1ANA")


ari_cohort <- ari_cohort %>% 
  left_join(ari_cov, by = "FINNGENID") %>% 
  select(-medication) %>% 
  # make them a factor
  mutate_at(vars(22:ncol(.)), as.factor)


roc_cohort <- roc_cohort %>% 
  left_join(roc_cov, by = "FINNGENID") %>% 
  select(-medication) %>% 
  # make them a factor
  mutate_at(vars(22:ncol(.)), as.factor)


tec_cohort <- tec_cohort %>% 
  left_join(tec_cov, by = "FINNGENID") %>% 
  select(-medication) %>% 
  # make them a factor
  mutate_at(vars(22:ncol(.)), as.factor)
tec_cohort <- tec_cohort %>% select(-"GLP1ANA", -"SULFONYLUREAS")



# Read and process PRS ----------------------------------------------------

# PRS files 
# read, scale, select FID and scaled PRS column

# PRS Scores - Covariates
# PRS Scores - Covariates

prefix_prs <- prefix_prs2_
suffix_prs <- suffix_prs_

prs_t2d    <- prs_t2d_2_
prs_chd    <- prs_chd_2_
prs_af     <- prs_af_
prs_stroke <- prs_stroke_2_
prs_suba   <- prs_suba_
prs_bmi    <- prs_bmi_
prs_edu    <- prs_edu_
# HbA1c
prs_hba1c <- paste0(prefix_prs, prs_hba1c_, suffix_prs)
# BP
prs_sbp <- paste0(prefix_prs, prs_sbp_, suffix_prs)
# ALT
prs_alt <- paste0(prefix_prs, prs_alt_, suffix_prs)
# AST
prs_ast <- paste0(prefix_prs, prs_ast_, suffix_prs)
# CRP
prs_crp <- paste0(prefix_prs, prs_crp_, suffix_prs)
# Asthma
prs_asthma <- prs_asthma_
# HDL
prs_hdl <- paste0(prefix_prs, prs_hdl_, suffix_prs)
# LDL
prs_ldl <- paste0(prefix_prs, prs_ldl_, suffix_prs)
# Major Depressive Disorder
prs_mdd <- prs_mdd_
# Triglycerides
prs_tg <- paste0(prefix_prs, prs_tg_, suffix_prs)
# Liver Disease
prs_liver <- prs_liver_
# Heart Failure
prs_hf <- prs_hf_
# CKD
prs_ckd <- prs_ckd_



################################################################################



# Read and process PRS ----------------------------------------------------

# PRS files 
# read, scale, select FID and scaled PRS column
prs_chd    <- prs_prep(prs_chd)
prs_suba   <- prs_prep(prs_suba)
prs_t2d    <- prs_prep(prs_t2d)
prs_stroke <- prs_prep(prs_stroke)
prs_sbp    <- prs_prep(prs_sbp)
prs_bmi    <- prs_prep(prs_bmi)
prs_edu    <- prs_prep(prs_edu)
prs_liver  <- prs_prep(prs_liver)
prs_ast    <- prs_prep(prs_ast)
prs_alt    <- prs_prep(prs_alt)
prs_asthma <- prs_prep(prs_asthma)
prs_hba1c  <- prs_prep(prs_hba1c)
prs_ldl    <- prs_prep(prs_ldl)
prs_hdl    <- prs_prep(prs_hdl)
prs_mdd    <- prs_prep(prs_mdd)
prs_crp    <- prs_prep(prs_crp)
prs_hf     <- prs_prep(prs_hf)
prs_ckd    <- prs_prep(prs_ckd)
prs_af     <- prs_prep(prs_af)
prs_tg     <- prs_prep(prs_tg)



# Join Matching Covariates (PRS only) -------------------------------------

emp_cohort <- add_prs(emp_cohort)
emp_cohort_naive <- add_prs(emp_cohort_naive)

tec_cohort <- add_prs(tec_cohort)
tec_cohort_naive <- add_prs(tec_cohort_naive)

ari_cohort <- add_prs(ari_cohort)
ari_cohort_naive <- add_prs(ari_cohort_naive)

roc_cohort <- add_prs(roc_cohort)
roc_cohort_naive <- add_prs(roc_cohort_naive)


# Matching on most unbalanced covariates ----------------------------------
ps_matching(emp_cohort, "emp_psmatches", 22, 44, caliper = 0.1)
ps_matching(tec_cohort, "tec_psmatches", 22, 42, caliper = 0.1)
ps_matching(ari_cohort, "ari_psmatches", 22, 42, caliper = 0.01)
ps_matching(roc_cohort, "roc_psmatches", 22, 46, caliper = 0.1)



# Testing the significance of the difference in covariates and PRS --------

## EMPAREG

# PRS covariates naive
emp_cov_test_naive_prs <- gen_anal(emp_cohort_naive, "Empareg", "naive")
emp_cov_test_naive_prs %>% filter(Pval <= 0.05) 
emp_cov_test_naive_prs %>% filter(Pval > 0.05)
# PRS covariates crude
emp_cov_test_crude_prs <- gen_anal(emp_cohort, "Empareg", "unadjusted")
emp_cov_test_crude_prs %>% filter(Pval <= 0.05) 
emp_cov_test_crude_prs %>% filter(Pval > 0.05)
# PRS covariates after matching ON PHENOTYPE COVARIATES)
emp_cov_test_prs <- gen_anal(emp_psmatches, "Empareg", "PS matched")
emp_cov_test_prs %>% filter(Pval > 0.05) 


## Tecos

# PRS covariates naive
tec_cov_test_naive_prs <- gen_anal(tec_cohort_naive, "Tecos", "naive")
tec_cov_test_naive_prs %>% filter(Pval <= 0.05) 
tec_cov_test_naive_prs %>% filter(Pval > 0.05)
# PRS covariates crude
tec_cov_test_crude_prs <- gen_anal(tec_cohort, "Tecos", "unadjusted")
tec_cov_test_crude_prs %>% filter(Pval <= 0.05) 
tec_cov_test_crude_prs %>% filter(Pval > 0.05)
# PRS covariates after matching ON PHENOTYPE COVARIATES)
tec_cov_test_prs <- gen_anal(tec_psmatches, "Tecos", "PS matched")
tec_cov_test_prs %>% filter(Pval <= 0.05) 
tec_cov_test_prs %>% filter(Pval > 0.05) 


## Aristotle

# PRS covariates naive
ari_cov_test_naive_prs <- gen_anal(ari_cohort_naive, "Aristotle", "naive")
ari_cov_test_naive_prs %>% filter(Pval <= 0.05) 
ari_cov_test_naive_prs %>% filter(Pval > 0.05)
# PRS covariates crude
ari_cov_test_crude_prs <- gen_anal(ari_cohort, "Aristotle", "unadjusted")
ari_cov_test_crude_prs %>% filter(Pval <= 0.05) 
ari_cov_test_crude_prs %>% filter(Pval > 0.05)
# PRS covariates after matching ON PHENOTYPE COVARIATES)
ari_cov_test_prs <- gen_anal(ari_psmatches, "Aristotle", "PS matched")
ari_cov_test_prs %>% filter(Pval <= 0.05) 
ari_cov_test_prs %>% filter(Pval > 0.05) 


## Rocket

# PRS covariates naive
roc_cov_test_naive_prs <- gen_anal(roc_cohort_naive, "Rocket", "naive")
roc_cov_test_naive_prs %>% filter(Pval <= 0.05) 
roc_cov_test_naive_prs %>% filter(Pval > 0.05)
# PRS covariates crude
roc_cov_test_crude_prs <- gen_anal(roc_cohort, "Rocket", "unadjusted")
roc_cov_test_crude_prs %>% filter(Pval <= 0.05) 
roc_cov_test_crude_prs %>% filter(Pval > 0.05)
# PRS covariates after matching ON PHENOTYPE COVARIATES)
roc_cov_test_prs <- gen_anal(roc_psmatches, "Rocket", "PS matched")
roc_cov_test_prs %>% filter(Pval <= 0.05) 
roc_cov_test_prs %>% filter(Pval > 0.05) 


prs_analysis_all_trials <- rbind(emp_cov_test_naive_prs, emp_cov_test_crude_prs, emp_cov_test_prs, 
                                 tec_cov_test_naive_prs, tec_cov_test_crude_prs, tec_cov_test_prs,
                                 ari_cov_test_naive_prs, ari_cov_test_crude_prs, ari_cov_test_prs,
                                 roc_cov_test_naive_prs, roc_cov_test_crude_prs, roc_cov_test_prs)

# Replace values of PRS Names
prs_analysis_all_trials <- 
  prs_analysis_all_trials %>% 
  mutate(prs_name = case_when(
    prs_name == "prs_ckd" ~ "Chronic Kidney Disease",
    prs_name == "prs_ldl" ~ "LDL",
    prs_name == "prs_hdl" ~ "HDL",
    prs_name == "prs_sbp" ~ "Systolic Blood Pressure",
    prs_name == "prs_alt" ~ "ALT Enzyme",
    prs_name == "prs_ast" ~ "AST Enzyme",
    prs_name == "prs_hba1c" ~ "Hba1c",
    prs_name == "prs_crp" ~ "CRP",
    prs_name == "prs_bmi" ~ "Body Mass Index",
    prs_name == "prs_chd" ~ "Coronary Heart Disease",
    prs_name == "prs_asthma" ~ "Asthma",
    prs_name == "prs_liver" ~ "Liver Disease",
    prs_name == "prs_lipid" ~ "Hyperlipidemia",
    prs_name == "prs_mdd" ~ "Major Depression",
    prs_name == "prs_stroke" ~ "Stroke",
    prs_name == "prs_t2d" ~ "Type 2 Diabetes",
    prs_name == "prs_cad" ~ "CAD",
    prs_name == "prs_edu" ~ "Education",
    prs_name == "prs_tg" ~ "Triglycerides",
    prs_name == "prs_af" ~ "Atrial Fibrillation",
    prs_name == "prs_hf" ~ "Heart Failure",
    prs_name == "prs_vte" ~ "Venous Thromboembolism",
    prs_name == "prs_suba" ~ "Subarachnoid Hemorrhage"
  )) 

fwrite(prs_analysis_all_trials, prs_analysis_all_trials_)



prs_analysis_all_trials$type <- factor(prs_analysis_all_trials$type, levels = c("naive", "unadjusted", "PS matched"))

ari_prs_analysis <- prs_analysis_all_trials %>% filter(trial == "Aristotle") %>% 
  # sort by Pval of unadjusted (in the opposite order because it plots always inverse)
  arrange(desc(type),Pval) %>% group_by(type) %>%  mutate(rank = n() - row_number() +1)  %>% ungroup() %>% 
  group_by(prs_name) %>% 
  mutate(rank = if_else(type == "PS matched", rank[type == "unadjusted"][match(prs_name, prs_name[type == "unadjusted"])], rank),
         rank = if_else(type == "naive", rank[type == "unadjusted"][match(prs_name, prs_name[type == "unadjusted"])], rank))

emp_prs_analysis <- prs_analysis_all_trials %>% filter(trial == "Empareg") %>% 
  # sort by Pval of unadjusted (in the opposite order because it plots always inverse)
  arrange(desc(type),Pval) %>% group_by(type) %>%  mutate(rank = n() - row_number() +1)  %>% ungroup() %>% 
  group_by(prs_name) %>% 
  mutate(rank = if_else(type == "PS matched", rank[type == "unadjusted"][match(prs_name, prs_name[type == "unadjusted"])], rank),
         rank = if_else(type == "naive", rank[type == "unadjusted"][match(prs_name, prs_name[type == "unadjusted"])], rank))

roc_prs_analysis <- prs_analysis_all_trials %>% filter(trial == "Rocket") %>% 
  # sort by Pval of unadjusted (in the opposite order because it plots always inverse)
  arrange(desc(type),Pval) %>% group_by(type) %>%  mutate(rank = n() - row_number() +1)  %>% ungroup() %>% 
  group_by(prs_name) %>% 
  mutate(rank = if_else(type == "PS matched", rank[type == "unadjusted"][match(prs_name, prs_name[type == "unadjusted"])], rank),
         rank = if_else(type == "naive", rank[type == "naive"][match(prs_name, prs_name[type == "unadjusted"])], rank))

tec_prs_analysis <- prs_analysis_all_trials %>% filter(trial == "Tecos") %>% 
  # sort by Pval of unadjusted (in the opposite order because it plots always inverse)
  arrange(desc(type),Pval) %>% group_by(type) %>%  mutate(rank = n() - row_number() +1)  %>% ungroup() %>% 
  group_by(prs_name) %>% 
  mutate(rank = if_else(type == "PS matched", rank[type == "unadjusted"][match(prs_name, prs_name[type == "unadjusted"])], rank),
         rank = if_else(type == "naive", rank[type == "unadjusted"][match(prs_name, prs_name[type == "unadjusted"])], rank))

# improvement of the PRS imbalance 
# emp_prs_analysis %>% filter(prs_name == "CAD") %>% group_by(prs_name) %>% summarise(d = abs(sum(diff[type == "PS matched"])) - abs(sum(diff[type == "unadjusted"])))
emp_prs_analysis_diff <- emp_prs_analysis %>% select(trial, type, prs_name, diff, rank) %>% 
  pivot_wider(names_from = type, values_from = diff) %>% 
  mutate(improve = abs(`PS matched`) - abs(unadjusted),
         improve_p = improve/abs(unadjusted) *100) 

tec_prs_analysis_diff <- tec_prs_analysis %>% select(trial, type, prs_name, diff, rank) %>% 
  pivot_wider(names_from = type, values_from = diff) %>% 
  mutate(improve = abs(`PS matched`) - abs(unadjusted),
         improve_p = improve/abs(unadjusted) *100,
         improve_p_adj = ifelse(improve_p > 150, 150, improve_p))

ari_prs_analysis_diff <- ari_prs_analysis %>% select(trial, type, prs_name, diff, rank) %>% 
  pivot_wider(names_from = type, values_from = diff) %>% 
  mutate(improve = abs(`PS matched`) - abs(unadjusted),
         improve_p = improve/abs(unadjusted) *100,
         improve_p_adj = ifelse(improve_p > 150, 150, improve_p)) 

roc_prs_analysis_diff <- roc_prs_analysis %>% select(trial, type, prs_name, diff, rank) %>% 
  pivot_wider(names_from = type, values_from = diff) %>% 
  mutate(improve = abs(`PS matched`) - abs(unadjusted),
         improve_p = improve/abs(unadjusted) *100,
         improve_p_adj = ifelse(improve_p > 150, 150, improve_p)) 


ggplot(emp_prs_analysis_diff, aes(x = reorder(prs_name, rank), y = improve_p, fill = improve_p > 0))+
  geom_bar(stat = "identity", alpha = 0.1)+
  scale_fill_manual(values = c("darkred", "darkgreen"), guide = FALSE) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),  # Entfernen der y-Achsenbeschriftung
        axis.text.y = element_blank())+
  labs(x = "", y = "Balance change in % ", title = "")+
  coord_flip()  +
   theme(axis.text.x = element_text(angle = 0, size = 15),    
         axis.title.x = element_text(size = 15),
         axis.title.y = element_blank(),  # Entfernen der y-Achsenbeschriftung
         axis.text.y = element_blank(),
         panel.spacing = unit(2, "lines"),
         legend.position = "top",
         legend.text = element_text(size = 16),
         legend.title = element_blank(), #element_text(size = 11),
         plot.margin = margin(1, 1, 1, 1, "cm"), 
         plot.title = element_text(hjust = 0.5, size = 20)
   )


# Based on the followiing numbers you will adjust your scales in the plot 
max(prs_analysis_all_trials$ci_upper) 
min(prs_analysis_all_trials$ci_lower) 

 
pd <- position_dodge(width=1)

 ggplot(emp_prs_analysis, aes(x = reorder(prs_name, rank), y = diff, group = type)) +
  geom_point(aes( shape = type, color = type), position = pd, size = 2) +
  geom_linerange(aes(ymin = ci_lower, ymax = ci_upper, color = type), position = pd) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
   geom_point(data = subset(emp_prs_analysis, Pval <= (0.05/21)), pch =21, fill = NA, size = 5, colour = "black", stroke = 0.5) +
  theme_minimal() +
  coord_flip() + 
  ggtitle("Differences in PRS before and after covariate adjustment  -  Empareg") +   # change trial name for each trial!!
  guides(color = "none", shape = FALSE) +
  ylim(c(-0.5, 0.55)) +
  labs(y = "Normalized difference in PRS", x = "", shape = "Adjustment")+
  theme(axis.text.x = element_text(angle = 0, size = 15),    
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 14),
        panel.spacing = unit(2, "lines"),
        legend.position = "top",
        legend.text = element_text(size = 16),
        legend.title = element_blank(), 
        plot.margin = margin(1, 1, 1, 1, "cm"), 
        plot.title = element_text(hjust = 0.5, size = 20,),
  )+
   facet_grid(.~type, scales = "free_x", labeller = labeller(type = c("naive" = "Naive", "unadjusted" = "Unadjusted", "PS matched" = "Adjusted"))) +
   theme(strip.text = element_blank()) +  
   geom_text(aes(x = 21, y = min(diff), label = case_when(type == "unadjusted" ~ "Unadjusted",
                                                          type == "PS matched" ~ "Adjusted", 
                                                          type ==  "naive" ~ "Naive"),  hjust = 1.5, vjust = 0))

 
 

# Save Plots as PDF -------------------------------------------------------
 
a1 <- ggplot(emp_prs_analysis, aes(x = reorder(prs_name, rank), y = diff, group = type)) +
   geom_point(aes(color = type), position = pd, size = 2) +
   geom_linerange(aes(ymin = ci_lower, ymax = ci_upper, color = type), position = pd) +
   geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
   geom_point(data = subset(emp_prs_analysis, Pval <= (0.05/21)), pch =21, fill = NA, size = 4, colour = "black", stroke = 0.5) +
   theme_minimal() +
   coord_flip() + 
   guides(color = "none", shape = FALSE) +
   ylim(c(-0.3, 0.6)) +
   labs(y = "Standardized Mean Difference", x = "Polygenic Score", shape = "Adjustment")+
   theme(axis.text.x = element_text(angle = 0, size = 15),    
         axis.title.x = element_text(size = 15),
         axis.text.y = element_text(size = 15),
         axis.title.y = element_text(size = 15),
         panel.spacing = unit(2, "lines"),
         legend.position = "top",
         legend.text = element_text(size = 16),
         legend.title = element_blank(),
         plot.margin = margin(1, 1, 1, 1, "cm"), 
         plot.title = element_text(hjust = 0.5, size = 20,),
   ) +
  facet_grid(.~type, scales = "free_x", labeller = labeller(type = c("naive" = "Plain Observational", "unadjusted" = "After Eligibility Criteria", "PS matched" = "Propensity Score Adjusted"))) +
  theme(strip.text = element_text(size = 15)) 


pdf(emp_prs_shift_, width = 13, height = 7)
a1
dev.off()


a2 <- ggplot(tec_prs_analysis, aes(x = reorder(prs_name, rank), y = diff, group = type)) +
  geom_point(aes( shape = type, color = type), position = pd, size = 2) +
  geom_linerange(aes(ymin = ci_lower, ymax = ci_upper, color = type), position = pd) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(tec_prs_analysis, Pval <= (0.05/21)), pch =21, fill = NA, size = 5, colour = "black", stroke = 0.5) +
  theme_minimal() +
  coord_flip() + 
  guides(color = "none", shape = FALSE) +
  ylim(c(-0.3, 0.6)) +
  labs(y = "Standardized Mean Difference", x = "Polygenic Score", shape = "Adjustment")+
  theme(axis.text.x = element_text(angle = 0, size = 15),    
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.spacing = unit(2, "lines"),
        legend.position = "top",
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm"), 
        plot.title = element_text(hjust = 0.5, size = 20,),
  ) +
  facet_grid(.~type, scales = "free_x", labeller = labeller(type = c("naive" = "Plain Observational", "unadjusted" = "After Eligibility Criteria", "PS matched" = "Propensity Score Adjusted"))) +
  theme(strip.text = element_text(size = 15)) 



pdf(tec_prs_shift_, width = 13, height = 7)
a2
dev.off()




a3 <- ggplot(ari_prs_analysis, aes(x = reorder(prs_name, rank), y = diff, group = type)) +
  geom_point(aes( shape = type, color = type), position = pd, size = 2) +
  geom_linerange(aes(ymin = ci_lower, ymax = ci_upper, color = type), position = pd) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(ari_prs_analysis, Pval <= (0.05/21)), pch =21, fill = NA, size = 5, colour = "black", stroke = 0.5) +
  theme_minimal() +
  coord_flip() +
  guides(color = "none", shape = FALSE) +
  ylim(c(-0.3, 0.3)) +
  labs(y = "Standardized Mean Difference", x = "Polygenic Score", shape = "Adjustment")+
  theme(axis.text.x = element_text(angle = 0, size = 15),    
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.spacing = unit(2, "lines"),
        legend.position = "top",
        legend.text = element_text(size = 16),
        legend.title = element_blank(), 
        plot.margin = margin(1, 1, 1, 1, "cm"), 
        plot.title = element_text(hjust = 0.5, size = 20,),
  ) +
  facet_grid(.~type, scales = "free_x", labeller = labeller(type = c("naive" = "Plain Observational", "unadjusted" = "After Eligibility Criteria", "PS matched" = "Propensity Score Adjusted"))) +
  theme(strip.text = element_text(size = 15)) 
 


pdf(ari_prs_shift_, width = 13, height = 7)
a3
dev.off()


a4 <- ggplot(roc_prs_analysis, aes(x = reorder(prs_name, rank), y = diff, group = type)) +
  geom_point(aes( shape = type, color = type), position = pd, size = 2) +
  geom_linerange(aes(ymin = ci_lower, ymax = ci_upper, color = type), position = pd) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(roc_prs_analysis, Pval <= (0.05/21)), pch =21, fill = NA, size = 5, colour = "black", stroke = 0.5) +
  theme_minimal() +
  coord_flip() + 
  guides(color = "none", shape = FALSE) +
  ylim(c(-0.3, 0.3)) +
  labs(y = "Standardized Mean Difference", x = "Polygenic Score", shape = "Adjustment")+
  theme(axis.text.x = element_text(angle = 0, size = 15),    
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.spacing = unit(2, "lines"),
        legend.position = "top",
        legend.text = element_text(size = 16),
        legend.title = element_blank(), 
        plot.margin = margin(1, 1, 1, 1, "cm"), 
        plot.title = element_text(hjust = 0.5, size = 20,),
  ) +
  facet_grid(.~type, scales = "free_x", labeller = labeller(type = c("naive" = "Plain Observational", "unadjusted" = "After Eligibility Criteria", "PS matched" = "Propensity Score Adjusted"))) +
  theme(strip.text = element_text(size = 15)) 

pdf(roc_prs_shift_, width = 13, height = 7)
a4
dev.off()



# All trials
pdf(alltrials_prs_shift_, width = 26, height = 12)
plot_grid(a1, a2, a3, a4, ncol = 2)
dev.off()


# All trials (without Empareg) for supplementary
pdf(alltrials_woemp_prs_shift_, width = 13, height = 21)
plot_grid(a2, a3, a4, ncol = 1)
dev.off()
