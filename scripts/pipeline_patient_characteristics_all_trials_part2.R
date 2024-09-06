rm(list = ls())
rm(list = ls()[!ls() %in% c("long_ep", "meds", "raw", "mp", "vnro_dict")])
gc()

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
library(survminer)
library(ggfortify)
library(cobalt)
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
ps_matching <- function(data, output_file = "psmatches_rctd", col_start_cov = 10000, col_end_cov = 10000, caliper = NULL) {
  # Get the names of all covariate columns (na.omit is in order to remove na as a result of no additional covariates specified)
  covariate_names <- na.omit(c("age", "SEX", "CCI", "birthyear", "age*birthyear", names(data)[col_start_cov:col_end_cov]))
  
  # Create a formula with covariates
  formula <- as.formula(paste("group ~", paste(covariate_names, collapse = " + ")))
  
  m.out <- matchit(formula,                                                    
                   data = data,
                   method = "nearest", 
                   distance = "logit", 
                   caliper = caliper) 
  # cstat_pre <- Cstat(glm(formula = formula, data = data))
  
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
  
  return(list(
    matching_summary,
    Love_Plot = matching_summary_plot
    # Cstat1 = cstat_pre
  ))
}


# PRS Analysis
gen_anal <- function(data, trial, type){
  # All PRS columns
  prs_columns <- grep("^prs_", names(data), value = TRUE)
  
  df <- data.frame(matrix(nrow=0, ncol=6))
  names(df) <- c("trial", "type", "prs_name", "diff", "ci_lower", "ci_upper")
  
  for (prs in prs_columns){
    print(prs)
    average_diff <- glm(group ~ get(prs), family=binomial(), data = data)
    
    diff <- summary(average_diff)$coefficients[, "Estimate"][2]
    CI_lower <- confint(average_diff)[,"2.5 %"][2]
    CI_upper <- confint(average_diff)[,"97.5 %"][2]
    
    v <- data.frame(trial, type, prs, diff, CI_lower, CI_upper) 
    names(v) <- names(df)
    df <- rbind(df, v) %>% tibble() 
    
    
  }
  return(df)
  
}



# Assess Covariate Balance - aSMD
cov_balance <- function(before = emp_cohort, after = emp_psmatches, name = "Empareg PS Matching", cov_start = 22, cov_end = 98){
  
  # Balance Statistics of covariates of interest
  smd_before <- col_w_smd(before %>% select(cov_start:cov_end), treat = before$group)
  smd_after <- col_w_smd(after %>% select(cov_start:cov_end), treat = after$group)
  
  # Absolute SMDs of covariates of interest - unadjusted and PRS adjsuted 
  balance_stats <- data.frame(
    covariates = names(smd_before),
    unadjusted = abs(unname(smd_before)),
    PS_adjusted = abs(unname(smd_after))
  )
  
  # Reshape the dataframe for plotting
  balance_stats_longer <- pivot_longer(balance_stats, cols = c("unadjusted", "PS_adjusted"), names_to = "metric")
  
  
  p <- ggplot(balance_stats_longer, aes(x = value, y = covariates, color = metric)) +
    geom_point(size = 1) +
    labs(x = "Absolute Standardized Mean Difference", y = "Covariates", title = paste0("Love Plot: ", name)) +
    scale_color_manual(values = c("unadjusted" = "blue", "PS_adjusted" = "red")) +
    geom_vline(xintercept = 0.1, color = "black") +
    geom_vline(xintercept = 0.05, color = "black", linetype = "dashed") +
    theme_minimal()
  
  p
  
}



# Defining the opposite of %in% to exclude rows that meet exclusion criteria
'%!in%' <- function(x,y)!('%in%'(x,y)) 



# Paths -------------------------------------------------------------------

# File paths
# Trials

file1 <- emp_surv_            # Empareg
file2 <- ari_surv_      # Aristotle
file4 <- roc_surv_             # Rocket
file5 <- tec_surv_                # Tecos


# Covariates (RCTD)
file1_cov <- emp_ps_cov_                            # Cohort covariates
file2_cov <- ari_ps_cov_
file4_cov <- roc_ps_cov_
file5_cov <- tec_ps_cov_



# Read files --------------------------------------------------------------

# Read file paths
emp_cohort <- fread(file1)  
emp_cov <- fread(file1_cov)

ari_cohort <- fread(file2)  
ari_cov <- fread(file2_cov)

roc_cohort <- fread(file4)  
roc_cov <- fread(file4_cov)

tec_cohort <- fread(file5)  
tec_cov <- fread(file5_cov)



#_____________________________________________
#### Add covariates #####

##### Covariates (RCTD) #####

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
# remove covariate vectors that have only one value (all either 1 or 0)
# ari_cohort <- ari_cohort %>% select(-which(sapply(., function(col) length(unique(col)) == 1)))


roc_cohort <- roc_cohort %>% 
  left_join(roc_cov, by = "FINNGENID") %>% 
  select(-medication) %>% 
  # make them a factor
  mutate_at(vars(22:ncol(.)), as.factor)
# remove covariate vectors that have only one value (all either 1 or 0)
# roc_cohort <- roc_cohort %>% select(-which(sapply(., function(col) length(unique(col)) == 1)))


tec_cohort <- tec_cohort %>% 
  left_join(tec_cov, by = "FINNGENID") %>% 
  select(-medication) %>% 
  # make them a factor
  mutate_at(vars(22:ncol(.)), as.factor)
tec_cohort <- tec_cohort %>% select(-"GLP1ANA", -"SULFONYLUREAS")
#tec_cohort <- tec_cohort %>% select(-"A10BB")


# Matching ----------------------------------------------------------------

ps_matching(emp_cohort, "emp_psmatches", 22, 44, caliper = 0.1) # bis 98
ggplot(emp_psmatches, aes(x = distance, fill = as.factor(group))) + 
  geom_density(alpha = 0.5)

ps_matching(ari_cohort, "ari_psmatches", 22, 42, caliper =  0.01) # bis 89
ggplot(ari_psmatches, aes(x = distance, fill = as.factor(group))) + 
  geom_density(alpha = 0.5)


ps_matching(roc_cohort, "roc_psmatches", 22, 46, caliper = 0.1) # bis 90
ggplot(roc_psmatches, aes(x = distance, fill = as.factor(group))) + 
  geom_density(alpha = 0.5)

ps_matching(tec_cohort, "tec_psmatches", 22, 42, caliper = 0.1) # bis 97
ggplot(tec_psmatches, aes(x = distance, fill = as.factor(group))) + 
  geom_density(alpha = 0.5)



# Patient's characteristics -----------------------------------------------

## Function to calculate Cohen's d for continuous variables
cohens_d <- function(mean1, sd1, mean2, sd2) {
  pooled_sd <- sqrt((sd1^2 + sd2^2) / 2)
  d <- (mean1 - mean2) / pooled_sd
  return(d)
}

# Function to calculate Cohen's h for proportions
cohens_h <- function(p1, p2) {
  h <- 2 * (asin(sqrt(p1)) - asin(sqrt(p2)))
  return(h)
}

# Create a summary dataframe
create_summary <- function(data) {
  # Initialize empty list to store results
  results <- list()
  
  # Split data into two groups
  group1 <- data %>% filter(cohort_group == 1)
  group2 <- data %>% filter(cohort_group == 0)
  
  # Define the continuous variables
  continuous_vars <- c("age", "CCI", "birthyear")
  
  # Loop over each covariate
  for (covariate in colnames(data)[-c(1, 2)]) {
    if (covariate %in% continuous_vars) {
      # Calculate mean and SD for continuous variables in both groups
      mean1 <- mean(group1[[covariate]], na.rm = TRUE)
      sd1 <- sd(group1[[covariate]], na.rm = TRUE)
      mean2 <- mean(group2[[covariate]], na.rm = TRUE)
      sd2 <- sd(group2[[covariate]], na.rm = TRUE)
      
      # Calculate Cohen's d
      std_diff <- cohens_d(mean1, sd1, mean2, sd2)
      
      # Append to results
      results[[covariate]] <- data.frame(
        Covariate = covariate,
        Exposure_Absolute = NA,
        Exposure_Percentage = paste0(round(mean1, 2), " (", round(sd1, 2), ")"),
        Comparator_Absolute = NA,
        Comparator_Percentage = paste0(round(mean2, 2), " (", round(sd2, 2), ")"),
        Std_Diff = round(std_diff, 3)
      )
      
    } else {
      # Convert factors to numeric
      group1[[covariate]] <- as.numeric(as.character(group1[[covariate]]))
      group2[[covariate]] <- as.numeric(as.character(group2[[covariate]]))
      
      # Calculate absolute numbers and proportions for binary covariates in both groups
      abs1 <- sum(group1[[covariate]], na.rm = TRUE)
      prop1 <- mean(group1[[covariate]], na.rm = TRUE)
      abs2 <- sum(group2[[covariate]], na.rm = TRUE)
      prop2 <- mean(group2[[covariate]], na.rm = TRUE)
      
      # Calculate Cohen's h
      std_diff <- cohens_h(prop1, prop2)
      
      # Append to results
      results[[covariate]] <- data.frame(
        Covariate = covariate,
        Exposure_Absolute = abs1,
        Exposure_Percentage = paste0(round(prop1 * 100, 1), "%"),
        Comparator_Absolute = abs2,
        Comparator_Percentage = paste0(round(prop2 * 100, 1), "%"),
        Std_Diff = round(std_diff, 3)
      )
    }
  }
  
  # Combine all results into a single dataframe
  summary_df <- do.call(rbind, results)
  return(summary_df)
}


results <- fread(pat_char_)

# Empareg
emp <- emp_psmatches %>% 
  select(1:4) %>% 
  # female is coded as 1
  mutate(SEX = ifelse(SEX == "female", 1, 0)) %>% 
  left_join((results %>% filter(trial == "Empareg")), by = "FINNGENID") %>% 
  select(FINNGENID, cohort_group = medication, sex = SEX, age, everything(), -trial, -group)
create_summary(emp)


# Tecos
tec <- tec_psmatches %>% 
  select(1:4) %>% 
  # female is coded as 1
  mutate(SEX = ifelse(SEX == "female", 1, 0)) %>% 
  left_join((results %>% filter(trial == "Tecos")), by = "FINNGENID") %>% 
  select(FINNGENID, cohort_group = medication, sex = SEX, age, everything(), -trial, -group)
create_summary(tec)
  


# Aristotle
ari <- ari_psmatches %>% 
  select(1:4) %>% 
  # female is coded as 1
  mutate(SEX = ifelse(SEX == "female", 1, 0)) %>% 
  left_join((results %>% filter(trial == "Aristotle")), by = "FINNGENID") %>% 
  select(FINNGENID, cohort_group = medication, sex = SEX, age, everything(), -trial, -group)
create_summary(ari)


# Rocket
roc <- roc_psmatches %>% 
  select(1:4) %>% 
  # female is coded as 1
  mutate(SEX = ifelse(SEX == "female", 1, 0)) %>% 
  left_join((results %>% filter(trial == "Rocket")), by = "FINNGENID") %>% 
  select(FINNGENID, cohort_group = medication, sex = SEX, age, everything(), -trial, -group)
create_summary(roc)