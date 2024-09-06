rm(list = ls())
gc()

# Packages ----------------------------------------------------------------

# Load all required packages:

library(tidyverse)
library(data.table)
library(MatchIt)
library(survival)
library(survminer)
library(ggfortify)
library(cobalt)
library(DescTools)
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


# Cox Proportional Hazard Models and Kaplan Meyer Plots
# weights are a vector of 1 with the length of the dataset, unless otherwise specified
cox_kmp <- function(data = trial_cohort, prs = NULL, weights = rep(1, nrow(data)), trial = "Empareg", treated = "Empagliflozin", control = "DPP4i", title = "", save = F){
  # CoxPH 
  if (is.null(prs)){
    # CoxPH if PRS is NOT provided
    cox <- coxph(Surv(time, primary_outcome) ~ group, data = data, weights = weights)
    cox_summary <- summary(cox)
  } else {
    # CoxPH if PRS IS provided
    cox <- coxph(Surv(time, primary_outcome) ~ group*prs, data = data, weights = weights)
    cox_summary <- summary(cox)
  }
  
  # KMP
  kmfit <- survfit(Surv(time, primary_outcome) ~ group, data = data, weights = weights)
  survivalplot <- autoplot(kmfit,
                           xlab = "Time (years)",
                           ylab = "Survival probability",
                           surv.plot = "overall",
                           censor = FALSE
  ) +
    theme_minimal() +
    ggtitle(paste(trial, "Emulation", title)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.text = element_text(size = 15),
          legend.position = "top",
          axis.text.x = element_text(angle = 0, size = 15),    
          axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 15)) +
    scale_fill_discrete(name = "", labels = c(control, treated)) + 
    scale_color_discrete(name = "", labels = c(control, treated)) 
  
  if (save){
    dir_path <- "/home/ivm/Trial_Emulations_R12/Assets/"
    file_name <- paste("prs_survival_plot_", gsub("\\s+", "_", tolower(title)), ".", "png", sep = "")
    file_path <- file.path(dir_path, file_name)
    ggsave(filename = file_path, plot = survivalplot)
  }
  
  return(
    list(cox_summary = cox_summary, survivalplot = survivalplot)
  )
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
file2 <- ari_surv_            # Aristotle
file4 <- roc_surv_            # Rocket
file5 <- tec_surv_            # Tecos

# PRS Scores - Covariates
prs_t2d    <- paste0(prefix_prs_, prs_t2d_, suffix_prs_)
prs_chd    <- paste0(prefix_prs_, prs_chd_, suffix_prs_)
prs_stroke <- paste0(prefix_prs_, prs_stroke_, suffix_prs_)


# Covariates
file1_cov <- emp_ps_cov_                         
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


prs_t2d    <- prs_prep(prs_t2d)
prs_chd    <- prs_prep(prs_chd)
prs_stroke <- prs_prep(prs_stroke)




#_____________________________________________
#### Add covariates #####

##### Covariates #####

emp_cohort <- emp_cohort %>% 
  left_join(emp_cov, by = "FINNGENID") %>% 
  select(-medication) %>% 
  # make them a factor
  mutate_at(vars(22:ncol(.)), as.factor)
#emp_cohort <- emp_cohort %>% select(-"GLP1ANA")
emp_cohort <- add_prs(emp_cohort)


ari_cohort <- ari_cohort %>% 
  left_join(ari_cov, by = "FINNGENID") %>% 
  select(-medication) %>% 
  # make them a factor
  mutate_at(vars(22:ncol(.)), as.factor)
ari_cohort <- add_prs(ari_cohort)


roc_cohort <- roc_cohort %>% 
  left_join(roc_cov, by = "FINNGENID") %>% 
  select(-medication) %>% 
  # make them a factor
  mutate_at(vars(22:ncol(.)), as.factor)
roc_cohort <- add_prs(roc_cohort)


tec_cohort <- tec_cohort %>% 
  left_join(tec_cov, by = "FINNGENID") %>% 
  select(-medication) %>% 
  # make them a factor
  mutate_at(vars(22:ncol(.)), as.factor)
#tec_cohort <- tec_cohort %>% select(-"GLP1ANA", -"SULFONYLUREAS")
tec_cohort <- add_prs(tec_cohort)


# Matching ----------------------------------------------------------------

ps_matching(emp_cohort, "emp_psmatches", 22, 44, caliper = 0.1) 
ggplot(emp_psmatches, aes(x = distance, fill = as.factor(group))) + 
  geom_density(alpha = 0.5)

ps_matching(ari_cohort, "ari_psmatches", 22, 42, caliper =  0.01)
ggplot(ari_psmatches, aes(x = distance, fill = as.factor(group))) + 
  geom_density(alpha = 0.5)

ps_matching(roc_cohort, "roc_psmatches", 22, 46, caliper = 0.1) 
ggplot(roc_psmatches, aes(x = distance, fill = as.factor(group))) + 
  geom_density(alpha = 0.5)

ps_matching(tec_cohort, "tec_psmatches", 22, 42, caliper = 0.1) 
ggplot(tec_psmatches, aes(x = distance, fill = as.factor(group))) + 
  geom_density(alpha = 0.5)

# Covariate Balance before and after PRS Matching -------------------------

#### Graphical assessment ####

# Empareg
cov_balance(before = emp_cohort, after = emp_psmatches, name = "Empareg PS Matching", cov_start = 22, cov_end = 45)

# Empareg TMP
cov_balance(before = emp_tmp_cohort, after = emp_tmp_psmatches3, name = "Empareg PS Matching", cov_start = 22, cov_end = 98)

# Aristotle
cov_balance(before = ari_cohort, after = ari_psmatches, name = "Aristotle PS Matching", cov_start = 22, cov_end = 89)

# Rocket
cov_balance(before = roc_cohort, after = roc_psmatches, name = "Rocket PS Matching", cov_start = 22, cov_end = 90)

# Tecos
cov_balance(before = tec_cohort, after = tec_psmatches, name = "Tecos PS Matching", cov_start = 22, cov_end = 42)


# Survival Analysis (RCTD Matches) ----------------------------------------

# Empareg
cox_kmp(emp_psmatches, prs = emp_psmatches$prs_chd, title = "ps_matched_rctd", trial = "Empareg", treated = "Empagliflozin", control = "DPP4i", save = F)
cox_kmp(emp_cohort, title = "UNMATCHED", trial = "Empareg", treated = "Empagliflozin", control = "DPP4i", save = F)

# Aristotle
cox_kmp(ari_psmatches,  prs = ari_psmatches$prs_stroke, title = "ps_matched_rctd", trial = "Aristotle", treated = "Apixaban", control = "Warfarin", save = F)
cox_kmp(ari_cohort, title = "UNMATCHED", trial = "Aristotle", treated = "Apixaban", control = "Warfarin", save = F)

# Rocket
cox_kmp(roc_psmatches, prs = roc_psmatches$prs_stroke, title = "ps_matched_rctd",  trial = "Rocket", treated = "Rivaroxaban", control = "Warfarin", save = F)
cox_kmp(roc_cohort, title = "UNMATCHED", trial = "Rocket", treated = "Rivaroxaban", control = "Warfarin", save = F)

# Tecos
cox_kmp(tec_psmatches, prs = tec_psmatches$prs_chd, title = "ps_matched_rctd", trial = "Tecos", treated = "Sitagliptin", control = "2nd Gen SU", save = F)
cox_kmp(tec_cohort, title = "UNMATCHED", trial = "Tecos", treated = "Sitagliptin", control = "2nd Gen SU", save = F)




#### Save plots as PDF: ####

# Empareg
pdf(emp_surv_plot_, width = 13, height = 7)
cox_kmp(emp_psmatches2, prs = emp_psmatches2$prs_chd, title = "", trial = "Empareg", treated = "Empagliflozin", control = "DPP4 Inhibitors", save = F)[2]
dev.off()


# Tecos
pdf(tec_surv_plot_, width = 13, height = 7)
cox_kmp(tec_psmatches, title = "", trial = "Tecos", treated = "Sitagliptin", control = "2nd generation Sulfonylurea", save = F)[2]
dev.off()



# Aristotle
pdf(ari_surv_plot_, width = 13, height = 7)
cox_kmp(ari_psmatches, title = "", trial = "Aristotle", treated = "Apixaban", control = "Warfarin", save = F)[2]
dev.off()



# Rocket
pdf(roc_surv_plot_, width = 13, height = 7)
cox_kmp(roc_psmatches, title = "", trial = "Rocket", treated = "Rivaroxaban", control = "Warfarin", save = F)[2]
dev.off()


#### Save plots as PDF 2: ####
plot1 <- cox_kmp(emp_psmatches2, prs = emp_psmatches2$prs_chd, title = "", trial = "Empareg", treated = "Empagliflozin", control = "DPP4 Inhibitors", save = F)[2]
plot2 <- cox_kmp(tec_psmatches, title = "", trial = "Tecos", treated = "Sitagliptin", control = "2nd generation Sulfonylurea", save = F)[2]
plot3 <- cox_kmp(ari_psmatches, title = "", trial = "Aristotle", treated = "Apixaban", control = "Warfarin", save = F)[2]
plot4 <- cox_kmp(roc_psmatches, title = "", trial = "Rocket", treated = "Rivaroxaban", control = "Warfarin", save = F)[2]

pdf(alltrials_surv_plot_, width = 16, height = 10)
plot_grid(plot1[[1]], plot2[[1]], plot3[[1]], plot4[[1]], ncol = 2)
dev.off()


# Events and event rates --------------------------------------------------

emp_psmatches %>% group_by(group) %>% 
  summarise(events = sum(primary_outcome),
            total_time = sum(time),
            event_rate = events/total_time*100)

tec_psmatches %>% group_by(group) %>% 
  summarise(events = sum(primary_outcome),
            total_time = sum(time),
            event_rate = events/total_time * 100)

ari_psmatches %>% group_by(group) %>% 
  summarise(events = sum(primary_outcome),
            total_time = sum(time),
            event_rate = events/total_time * 100)

roc_psmatches %>% group_by(group) %>% 
  summarise(events = sum(primary_outcome),
            total_time = sum(time),
            event_rate = events/total_time * 100)
