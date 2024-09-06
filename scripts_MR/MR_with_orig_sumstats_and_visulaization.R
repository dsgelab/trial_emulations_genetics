library(data.table)
library(TwoSampleMR)
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(purrr)
library(tibble)
library(stringr)
library(forcats)
source("/home/ivm/config_te.R")

# This script serves the purpose of conducting Mendelian Randomization (MR) of "potential confounding" variables (exposure)
# on the trial arm and the trial outcome (outcome).


# Trial outcome GWAS ------------------------------------------------------

# CAD = trial_outcome
cad_gwas <- fread(lifted_cad_ss_)
cad_iv <- fread(cad_iv_)[,c(1,4,3)] # chrom, pos, rsID/SNP

trial_outcome_data <- data.frame(
  id.outcome = "CAD",
  outcome = "CAD (CARDIoGRAMplusC4D 2015)",
  SNP = cad_gwas$rsid,
  beta.outcome = cad_gwas$effect_weight,
  se.outcome = cad_gwas$se,
  effect_allele.outcome = cad_gwas$effect_allele,
  other_allele.outcome = cad_gwas$other_allele
) %>% 
  mutate(chr.outcome = cad_gwas$chromosome,
         pos.outcome = cad_gwas$chr_position,
         pval.outcome = cad_gwas$pval,
         samplesize.outcome = cad_gwas$samplesize,
         eaf.outcome = cad_gwas$af,
         mr_keep.outcome = TRUE,
         pval_origin.outcome = "reported",
         data_source.outcome = "") %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.outcome = as.integer(SNP %in% cad_iv$SNP))




# Trial arm GWAS ----------------------------------------------------------

# EMPA = trial_arm (initiators vs non-initiators)
empa_plain_gwas <- fread(empa_plain_gwas_)
empa_plain_iv <- fread(empa_plain_iv_)[,c(1,4,3)] # chrom, pos, rsID/SNP

trial_arm_plain_data <- data.frame(
  id.outcome = "Empagliflozin plain",
  outcome = "Empagliflozin Initiators (FG 2024)",
  SNP = empa_plain_gwas$rsid,
  beta.outcome = empa_plain_gwas$effect_weight,
  se.outcome = empa_plain_gwas$se,
  effect_allele.outcome = empa_plain_gwas$effect_allele,
  other_allele.outcome = empa_plain_gwas$other_allele
) %>%
  mutate(chr.outcome = empa_plain_gwas$chromosome,
         pos.outcome = empa_plain_gwas$chr_position,
         pval.outcome = empa_plain_gwas$pval,
         samplesize.outcome = empa_plain_gwas$samplesize,
         eaf.outcome = empa_plain_gwas$af,
         mr_keep.outcome = TRUE,
         pval_origin.outcome = "reported",
         data_source.outcome = "") %>%
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>%
  # add a column indicating whether it is an index SNP or not
  mutate(iv.outcome = as.integer(SNP %in% empa_plain_iv$SNP))


# EMPA = trial_arm (empareg after eligibility criteria)
empa_gwas <- fread(empa_gwas_)
empa_iv <- NULL

trial_arm_data <- data.frame(
  id.outcome = "Empagliflozin",
  outcome = "Empagliflozin Treatment after eligibility criteria (FG 2024)",
  SNP = empa_gwas$rsid,
  beta.outcome = empa_gwas$effect_weight,
  se.outcome = empa_gwas$se,
  effect_allele.outcome = empa_gwas$effect_allele,
  other_allele.outcome = empa_gwas$other_allele
) %>% 
  mutate(chr.outcome = empa_gwas$chromosome,
         pos.outcome = empa_gwas$chr_position,
         pval.outcome = empa_gwas$pval,
         samplesize.outcome = empa_gwas$samplesize,
         eaf.outcome = empa_gwas$af,
         mr_keep.outcome = TRUE,
         pval_origin.outcome = "reported",
         data_source.outcome = "") %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.outcome = as.integer(SNP %in% empa_iv$SNP))




# MR with the exposure variables ------------------------------------------

#### MR BMI - CAD ####
bmi_gwas <- fread(lifted_bmi_ss_)
bmi_iv <- fread(bmi_iv_)[,c(1,4,3)]

bmi_data <- data.frame(
  id.exposure = "BMI",
  exposure = "Body Mass Index (GIANT 2018)",
  SNP = bmi_gwas$rsid,
  beta.exposure = bmi_gwas$effect_weight,
  se.exposure = bmi_gwas$se,
  effect_allele.exposure = bmi_gwas$effect_allele,
  other_allele.exposure = bmi_gwas$other_allele
) %>% 
  mutate(chr.exposure = bmi_gwas$chromosome,
         pos.exposure = bmi_gwas$chr_position,
         pval.exposure = bmi_gwas$pval,
         samplesize.exposure = bmi_gwas$samplesize,
         eaf.exposure = bmi_gwas$af,
         mr_keep.exposure = TRUE,
         pval_origin.exposure = "reported",
         data_source.exposure = "") %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.exposure = as.integer(SNP %in% bmi_iv$SNP)) %>%
  # keep only IVs
  filter(iv.exposure == 1) 

# remove IVs if they are also IVs in the second GWAS
bmi_data_cad <- bmi_data %>%
  filter(!(SNP %in% cad_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
cad_data_bmi_ivs <- trial_outcome_data %>% filter(SNP %in% bmi_data_cad$SNP)


bmi_on_cad <- harmonise_data(bmi_data_cad, cad_data_bmi_ivs)
bmi_on_cad_result <- mr(bmi_on_cad)
bmi_on_cad_plot <- mr_scatter_plot(bmi_on_cad_result, bmi_on_cad)
bmi_on_cad_plot


#### MR BMI - Empa plain ####

# remove IVs if they are also IVs in the second GWAS
bmi_data_empa <- bmi_data %>%
  filter(!(SNP %in% empa_plain_iv$SNP))

empa_plain_data_bmi_ivs <- trial_arm_plain_data %>% filter(SNP %in% bmi_data_empa$SNP)

bmi_on_empa <- harmonise_data(bmi_data_empa, empa_plain_data_bmi_ivs)
bmi_on_empa_plain_result <- mr(bmi_on_empa)
bmi_on_empa_plain_plot <- mr_scatter_plot(bmi_on_empa_plain_result, bmi_on_empa)
bmi_on_empa_plain_plot


#### MR BMI - Empa ####

# remove IVs if they are also IVs in the second GWAS
bmi_data_empa <- bmi_data %>%
  filter(!(SNP %in% empa_iv$SNP))

empa_data_bmi_ivs <- trial_arm_data %>% filter(SNP %in% bmi_data_empa$SNP)

bmi_on_empa <- harmonise_data(bmi_data_empa, empa_data_bmi_ivs)
bmi_on_empa_result <- mr(bmi_on_empa)
bmi_on_empa_plot <- mr_scatter_plot(bmi_on_empa_result, bmi_on_empa)
bmi_on_empa_plot



#### MR LDL - CAD ####
ldl_gwas <- fread(lifted_ldl_ss_)
ldl_iv <- fread(ldl_iv_)[,c(1,4,3)]

ldl_data <- data.frame(
  id.exposure = "LDL",
  exposure = "LDL (PAN-UKBB)",
  SNP = ldl_gwas$rsid,
  beta.exposure = ldl_gwas$effect_weight,
  se.exposure = ldl_gwas$se,
  effect_allele.exposure = ldl_gwas$effect_allele,
  other_allele.exposure = ldl_gwas$other_allele
) %>% 
  mutate(chr.exposure = ldl_gwas$chromosome,
         pos.exposure = ldl_gwas$chr_position,
         pval.exposure = ldl_gwas$pval,
         samplesize.exposure = ldl_gwas$samplesize,
         eaf.exposure = ldl_gwas$af,
         mr_keep.exposure = TRUE,
         pval_origin.exposure = "reported",
         data_source.exposure = "") %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.exposure = as.integer(SNP %in% ldl_iv$SNP)) %>%
  # keep only IVs
  filter(iv.exposure == 1) 

# remove IVs if they are also IVs in the second GWAS
ldl_data_cad <- ldl_data %>%
  filter(!(SNP %in% cad_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
cad_data_ldl_ivs <- trial_outcome_data %>% filter(SNP %in% ldl_data_cad$SNP)


ldl_on_cad <- harmonise_data(ldl_data_cad, cad_data_ldl_ivs)
ldl_on_cad_result <- mr(ldl_on_cad)
ldl_on_cad_plot <- mr_scatter_plot(ldl_on_cad_result, ldl_on_cad)
ldl_on_cad_plot


#### MR LDL - Empa plain ####

# remove IVs if they are also IVs in the second GWAS
ldl_data_empa <- ldl_data %>%
  filter(!(SNP %in% empa_plain_iv$SNP))

empa_plain_data_ldl_ivs <- trial_arm_plain_data %>% filter(SNP %in% ldl_data_empa$SNP)

ldl_on_empa <- harmonise_data(ldl_data_empa, empa_plain_data_ldl_ivs)
ldl_on_empa_plain_result <- mr(ldl_on_empa)
ldl_on_empa_plain_plot <- mr_scatter_plot(ldl_on_empa_plain_result, ldl_on_empa)
ldl_on_empa_plain_plot


#### MR LDL - Empa ####

# remove IVs if they are also IVs in the second GWAS
ldl_data_empa <- ldl_data %>%
  filter(!(SNP %in% empa_iv$SNP))

empa_data_ldl_ivs <- trial_arm_data %>% filter(SNP %in% ldl_data_empa$SNP)

ldl_on_empa <- harmonise_data(ldl_data_empa, empa_data_ldl_ivs)
ldl_on_empa_result <- mr(ldl_on_empa)
ldl_on_empa_plot <- mr_scatter_plot(ldl_on_empa_result, ldl_on_empa)
ldl_on_empa_plot



#### MR HDL - CAD ####
hdl_gwas <- fread(lifted_hdl_ss_)
hdl_iv <- fread(hdl_iv_)[,c(1,4,3)]

hdl_data <- data.frame(
  id.exposure = "HDL",
  exposure = "HDL (PAN-UKBB)",
  SNP = hdl_gwas$rsid,
  beta.exposure = hdl_gwas$effect_weight,
  se.exposure = hdl_gwas$se,
  effect_allele.exposure = hdl_gwas$effect_allele,
  other_allele.exposure = hdl_gwas$other_allele
) %>% 
  mutate(chr.exposure = hdl_gwas$chromosome,
         pos.exposure = hdl_gwas$chr_position,
         pval.exposure = hdl_gwas$pval,
         samplesize.exposure = hdl_gwas$samplesize,
         eaf.exposure = hdl_gwas$af,
         mr_keep.exposure = TRUE,
         pval_origin.exposure = "reported",
         data_source.exposure = "") %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.exposure = as.integer(SNP %in% hdl_iv$SNP)) %>%
  # keep only IVs
  filter(iv.exposure == 1) 

# remove IVs if they are also IVs in the second GWAS
hdl_data_cad <- hdl_data %>%
  filter(!(SNP %in% cad_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
cad_data_hdl_ivs <- trial_outcome_data %>% filter(SNP %in% hdl_data_cad$SNP)


hdl_on_cad <- harmonise_data(hdl_data_cad, cad_data_hdl_ivs)
hdl_on_cad_result <- mr(hdl_on_cad)
hdl_on_cad_plot <- mr_scatter_plot(hdl_on_cad_result, hdl_on_cad)
hdl_on_cad_plot


#### MR HDL - Empa plain ####

# remove IVs if they are also IVs in the second GWAS
hdl_data_empa <- hdl_data %>%
  filter(!(SNP %in% empa_plain_iv$SNP))

empa_plain_data_hdl_ivs <- trial_arm_plain_data %>% filter(SNP %in% hdl_data_empa$SNP)

hdl_on_empa <- harmonise_data(hdl_data_empa, empa_plain_data_hdl_ivs)
hdl_on_empa_plain_result <- mr(hdl_on_empa)
hdl_on_empa_plain_plot <- mr_scatter_plot(hdl_on_empa_plain_result, hdl_on_empa)
hdl_on_empa_plain_plot


#### MR HDL - Empa ####

# remove IVs if they are also IVs in the second GWAS
hdl_data_empa <- hdl_data %>%
  filter(!(SNP %in% empa_iv$SNP))

empa_data_hdl_ivs <- trial_arm_data %>% filter(SNP %in% hdl_data_empa$SNP)

hdl_on_empa <- harmonise_data(hdl_data_empa, empa_data_hdl_ivs)
hdl_on_empa_result <- mr(hdl_on_empa)
hdl_on_empa_plot <- mr_scatter_plot(hdl_on_empa_result, hdl_on_empa)
hdl_on_empa_plot



#### MR CKD - CAD ####
ckd_gwas <- fread(lifted_ckd_ss_)
ckd_iv <- fread(ckd_iv_)[,c(1,4,3)]

ckd_data <- data.frame(
  id.exposure = "CKD",
  exposure = "CKD (2018)",
  SNP = ckd_gwas$rsid,
  beta.exposure = ckd_gwas$effect_weight,
  se.exposure = ckd_gwas$se,
  effect_allele.exposure = ckd_gwas$effect_allele,
  other_allele.exposure = ckd_gwas$other_allele
) %>% 
  mutate(chr.exposure = ckd_gwas$chromosome,
         pos.exposure = ckd_gwas$chr_position,
         pval.exposure = ckd_gwas$pval,
         samplesize.exposure = ckd_gwas$samplesize,
         eaf.exposure = ckd_gwas$af,
         mr_keep.exposure = TRUE,
         pval_origin.exposure = "reported",
         data_source.exposure = "") %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.exposure = as.integer(SNP %in% ckd_iv$SNP)) %>%
  # keep only IVs
  filter(iv.exposure == 1) 

# remove IVs if they are also IVs in the second GWAS
ckd_data_cad <- ckd_data %>%
  filter(!(SNP %in% cad_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
cad_data_ckd_ivs <- trial_outcome_data %>% filter(SNP %in% ckd_data_cad$SNP)


ckd_on_cad <- harmonise_data(ckd_data_cad, cad_data_ckd_ivs)
ckd_on_cad_result <- mr(ckd_on_cad)
ckd_on_cad_plot <- mr_scatter_plot(ckd_on_cad_result, ckd_on_cad)
ckd_on_cad_plot


#### MR CKD - Empa plain ####

# remove IVs if they are also IVs in the second GWAS
ckd_data_empa <- ckd_data %>%
  filter(!(SNP %in% empa_plain_iv$SNP))

empa_plain_data_ckd_ivs <- trial_arm_plain_data %>% filter(SNP %in% ckd_data_empa$SNP)

ckd_on_empa <- harmonise_data(ckd_data_empa, empa_plain_data_ckd_ivs)
ckd_on_empa_plain_result <- mr(ckd_on_empa)
ckd_on_empa_plain_plot <- mr_scatter_plot(ckd_on_empa_plain_result, ckd_on_empa)
ckd_on_empa_plain_plot


#### MR CKD - Empa ####

# remove IVs if they are also IVs in the second GWAS
ckd_data_empa <- ckd_data %>%
  filter(!(SNP %in% empa_iv$SNP))

empa_data_ckd_ivs <- trial_arm_data %>% filter(SNP %in% ckd_data_empa$SNP)

ckd_on_empa <- harmonise_data(ckd_data_empa, empa_data_ckd_ivs)
ckd_on_empa_result <- mr(ckd_on_empa)
ckd_on_empa_plot <- mr_scatter_plot(ckd_on_empa_result, ckd_on_empa)
ckd_on_empa_plot



#### MR EDU - CAD ####
edu_gwas <- fread(lifted_edu_ss_)
edu_iv <- fread(edu_iv_)[,c(1,4,3)]

edu_data <- data.frame(
  id.exposure = "Educational Attainment",
  exposure = "Educational Attainment (Lee 2018)",
  SNP = edu_gwas$rsid,
  beta.exposure = edu_gwas$effect_weight,
  se.exposure = edu_gwas$se,
  effect_allele.exposure = edu_gwas$effect_allele,
  other_allele.exposure = edu_gwas$other_allele
) %>% 
  mutate(chr.exposure = edu_gwas$chromosome,
         pos.exposure = edu_gwas$chr_position,
         pval.exposure = edu_gwas$pval,
         samplesize.exposure = edu_gwas$samplesize,
         eaf.exposure = NA,
         mr_keep.exposure = TRUE,
         pval_origin.exposure = "reported",
         data_source.exposure = "") %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.exposure = as.integer(SNP %in% edu_iv$SNP)) %>%
  # keep only IVs
  filter(iv.exposure == 1) 

# remove IVs if they are also IVs in the second GWAS
edu_data_cad <- edu_data %>%
  filter(!(SNP %in% cad_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
cad_data_edu_ivs <- trial_outcome_data %>% filter(SNP %in% edu_data_cad$SNP)


edu_on_cad <- harmonise_data(edu_data_cad, cad_data_edu_ivs)
edu_on_cad_result <- mr(edu_on_cad)
edu_on_cad_plot <- mr_scatter_plot(edu_on_cad_result, edu_on_cad)
edu_on_cad_plot


#### MR EDU - Empa plain ####

# remove IVs if they are also IVs in the second GWAS
edu_data_empa <- edu_data %>%
  filter(!(SNP %in% empa_plain_iv$SNP))

empa_plain_data_edu_ivs <- trial_arm_plain_data %>% filter(SNP %in% edu_data_empa$SNP)

edu_on_empa <- harmonise_data(edu_data_empa, empa_plain_data_edu_ivs)
edu_on_empa_plain_result <- mr(edu_on_empa)
edu_on_empa_plain_plot <- mr_scatter_plot(edu_on_empa_plain_result, edu_on_empa)
edu_on_empa_plain_plot


#### MR EDU - Empa ####

# remove IVs if they are also IVs in the second GWAS
edu_data_empa <- edu_data %>%
  filter(!(SNP %in% empa_iv$SNP))

empa_data_edu_ivs <- trial_arm_data %>% filter(SNP %in% edu_data_empa$SNP)

edu_on_empa <- harmonise_data(edu_data_empa, empa_data_edu_ivs)
edu_on_empa_result <- mr(edu_on_empa)
edu_on_empa_plot <- mr_scatter_plot(edu_on_empa_result, edu_on_empa)
edu_on_empa_plot



#### MR MDD - CAD ####
mdd_gwas <- fread(lifted_mdd_ss_)
mdd_iv <- fread(mdd_iv_)[,c(1,4,3)]

mdd_data <- data.frame(
  id.exposure = "Major Depressive Disorder",
  exposure = "Major Depressive Disorder (Wray 2018)",
  SNP = mdd_gwas$rsid,
  beta.exposure = mdd_gwas$effect_weight,
  se.exposure = mdd_gwas$se,
  effect_allele.exposure = mdd_gwas$effect_allele,
  other_allele.exposure = mdd_gwas$other_allele
) %>% 
  mutate(chr.exposure = mdd_gwas$chromosome,
         pos.exposure = mdd_gwas$chr_position,
         pval.exposure = mdd_gwas$pval,
         samplesize.exposure = mdd_gwas$samplesize,
         eaf.exposure = NA,
         mr_keep.exposure = TRUE,
         pval_origin.exposure = "reported",
         data_source.exposure = "") %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.exposure = as.integer(SNP %in% mdd_iv$SNP)) %>%
  # keep only IVs
  filter(iv.exposure == 1) 

# remove IVs if they are also IVs in the second GWAS
mdd_data_cad <- mdd_data %>%
  filter(!(SNP %in% cad_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
cad_data_mdd_ivs <- trial_outcome_data %>% filter(SNP %in% mdd_data_cad$SNP)


mdd_on_cad <- harmonise_data(mdd_data_cad, cad_data_mdd_ivs)
mdd_on_cad_result <- mr(mdd_on_cad)
mdd_on_cad_plot <- mr_scatter_plot(mdd_on_cad_result, mdd_on_cad)
mdd_on_cad_plot


#### MR MDD - Empa plain ####

# remove IVs if they are also IVs in the second GWAS
mdd_data_empa <- mdd_data %>%
  filter(!(SNP %in% empa_plain_iv$SNP))

empa_plain_data_mdd_ivs <- trial_arm_plain_data %>% filter(SNP %in% mdd_data_empa$SNP)

mdd_on_empa <- harmonise_data(mdd_data_empa, empa_plain_data_mdd_ivs)
mdd_on_empa_plain_result <- mr(mdd_on_empa)
mdd_on_empa_plain_plot <- mr_scatter_plot(mdd_on_empa_plain_result, mdd_on_empa)
mdd_on_empa_plain_plot


#### MR MDD - Empa ####

# remove IVs if they are also IVs in the second GWAS
mdd_data_empa <- mdd_data %>%
  filter(!(SNP %in% empa_iv$SNP))

empa_data_mdd_ivs <- trial_arm_data %>% filter(SNP %in% mdd_data_empa$SNP)

mdd_on_empa <- harmonise_data(mdd_data_empa, empa_data_mdd_ivs)
mdd_on_empa_result <- mr(mdd_on_empa)
mdd_on_empa_plot <- mr_scatter_plot(mdd_on_empa_result, mdd_on_empa)
mdd_on_empa_plot


#### MR Asthma - CAD ####
asthma_gwas <- fread(lifted_asthma_ss_)
asthma_iv <- fread(asthma_iv_)[,c(1,4,3)]

asthma_data <- data.frame(
  id.exposure = "Asthma",
  exposure = "Asthma (Han 2020)",
  SNP = asthma_gwas$rsid,
  beta.exposure = asthma_gwas$effect_weight,
  se.exposure = asthma_gwas$se,
  effect_allele.exposure = asthma_gwas$effect_allele,
  other_allele.exposure = asthma_gwas$other_allele
) %>% 
  mutate(chr.exposure = asthma_gwas$chromosome,
         pos.exposure = asthma_gwas$chr_position,
         pval.exposure = asthma_gwas$pval,
         samplesize.exposure = asthma_gwas$samplesize,
         eaf.exposure = asthma_gwas$af,
         mr_keep.exposure = TRUE,
         pval_origin.exposure = "reported",
         data_source.exposure = "") %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.exposure = as.integer(SNP %in% asthma_iv$SNP)) %>%
  # keep only IVs
  filter(iv.exposure == 1) 

# remove IVs if they are also IVs in the second GWAS
asthma_data_cad <- asthma_data %>%
  filter(!(SNP %in% cad_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
cad_data_asthma_ivs <- trial_outcome_data %>% filter(SNP %in% asthma_data_cad$SNP)


asthma_on_cad <- harmonise_data(asthma_data_cad, cad_data_asthma_ivs)
asthma_on_cad_result <- mr(asthma_on_cad)
asthma_on_cad_plot <- mr_scatter_plot(asthma_on_cad_result, asthma_on_cad)
asthma_on_cad_plot


#### MR Asthma - Empa plain ####

# remove IVs if they are also IVs in the second GWAS
asthma_data_empa <- asthma_data %>%
  filter(!(SNP %in% empa_plain_iv$SNP))

empa_plain_data_asthma_ivs <- trial_arm_plain_data %>% filter(SNP %in% asthma_data_empa$SNP)

asthma_on_empa <- harmonise_data(asthma_data_empa, empa_plain_data_asthma_ivs)
asthma_on_empa_plain_result <- mr(asthma_on_empa)
asthma_on_empa_plain_plot <- mr_scatter_plot(asthma_on_empa_plain_result, asthma_on_empa)
asthma_on_empa_plain_plot


#### MR Asthma - Empa ####

# remove IVs if they are also IVs in the second GWAS
asthma_data_empa <- asthma_data %>%
  filter(!(SNP %in% empa_iv$SNP))

empa_data_asthma_ivs <- trial_arm_data %>% filter(SNP %in% asthma_data_empa$SNP)

asthma_on_empa <- harmonise_data(asthma_data_empa, empa_data_asthma_ivs)
asthma_on_empa_result <- mr(asthma_on_empa)
asthma_on_empa_plot <- mr_scatter_plot(asthma_on_empa_result, asthma_on_empa)
asthma_on_empa_plot



#### MR T2D - CAD ####
t2d_gwas <- fread(write_t2d_ss_)
t2d_gwas$chromosome <- as.numeric(sub("^chr", "", t2d_gwas$chr_name)) 
t2d_iv <- fread(t2d_iv_)[,c(1,4,3)]

t2d_data <- data.frame(
  id.exposure = "T2D",
  exposure = "T2D (Suzuki 2024)",
  SNP = t2d_gwas$rsID,
  beta.exposure = t2d_gwas$effect_weight,
  se.exposure = t2d_gwas$se,
  effect_allele.exposure = t2d_gwas$effect_allele,
  other_allele.exposure = t2d_gwas$other_allele
) %>% 
  mutate(chr.exposure = t2d_gwas$chromosome,
         pos.exposure = t2d_gwas$chr_position,
         pval.exposure = t2d_gwas$pval,
         samplesize.exposure = t2d_gwas$samplesize,
         eaf.exposure = 0.5, # random number, is anyways not important
         mr_keep.exposure = TRUE,
         pval_origin.exposure = "reported",
         data_source.exposure = "") %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.exposure = as.integer(SNP %in% t2d_iv$SNP)) %>%
  # keep only IVs
  filter(iv.exposure == 1) 

# remove IVs if they are also IVs in the second GWAS
t2d_data_cad <- t2d_data %>%
  filter(!(SNP %in% cad_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
cad_data_t2d_ivs <- trial_outcome_data %>% filter(SNP %in% t2d_data_cad$SNP)


t2d_on_cad <- harmonise_data(t2d_data_cad, cad_data_t2d_ivs)
t2d_on_cad_result <- mr(t2d_on_cad)
t2d_on_cad_plot <- mr_scatter_plot(t2d_on_cad_result, t2d_on_cad)
t2d_on_cad_plot


#### MR T2D - Empa plain ####

# remove IVs if they are also IVs in the second GWAS
t2d_data_empa <- t2d_data %>%
  filter(!(SNP %in% empa_plain_iv$SNP))

empa_plain_data_t2d_ivs <- trial_arm_plain_data %>% filter(SNP %in% t2d_data_empa$SNP)

t2d_on_empa <- harmonise_data(t2d_data_empa, empa_plain_data_t2d_ivs)
t2d_on_empa_plain_result <- mr(t2d_on_empa)
t2d_on_empa_plain_plot <- mr_scatter_plot(t2d_on_empa_plain_result, t2d_on_empa)
t2d_on_empa_plain_plot


#### MR T2D - Empa ####

# remove IVs if they are also IVs in the second GWAS
t2d_data_empa <- t2d_data %>%
  filter(!(SNP %in% empa_iv$SNP))

empa_data_t2d_ivs <- trial_arm_data %>% filter(SNP %in% t2d_data_empa$SNP)

t2d_on_empa <- harmonise_data(t2d_data_empa, empa_data_t2d_ivs)
t2d_on_empa_result <- mr(t2d_on_empa)
t2d_on_empa_plot <- mr_scatter_plot(t2d_on_empa_result, t2d_on_empa)
t2d_on_empa_plot



#### MR Hba1c - CAD ####
hba1c_gwas <- fread(lifted_hba1c_ss_)
hba1c_iv <- fread(hba1c_iv_)[,c(1,4,3)]

hba1c_data <- data.frame(
  id.exposure = "Hba1c",
  exposure = "Hba1c (PAN-UKBB)",
  SNP = hba1c_gwas$rsid,
  beta.exposure = hba1c_gwas$effect_weight,
  se.exposure = hba1c_gwas$se,
  effect_allele.exposure = hba1c_gwas$effect_allele,
  other_allele.exposure = hba1c_gwas$other_allele
) %>% 
  mutate(chr.exposure = hba1c_gwas$chromosome,
         pos.exposure = hba1c_gwas$chr_position,
         pval.exposure = hba1c_gwas$pval,
         samplesize.exposure = hba1c_gwas$samplesize,
         eaf.exposure = NA,
         mr_keep.exposure = TRUE,
         pval_origin.exposure = "reported",
         data_source.exposure = "") %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.exposure = as.integer(SNP %in% hba1c_iv$SNP)) %>%
  # keep only IVs
  filter(iv.exposure == 1) 

# remove IVs if they are also IVs in the second GWAS
hba1c_data_cad <- hba1c_data %>%
  filter(!(SNP %in% cad_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
cad_data_hba1c_ivs <- trial_outcome_data %>% filter(SNP %in% hba1c_data_cad$SNP)


hba1c_on_cad <- harmonise_data(hba1c_data_cad, cad_data_hba1c_ivs)
hba1c_on_cad_result <- mr(hba1c_on_cad)
hba1c_on_cad_plot <- mr_scatter_plot(hba1c_on_cad_result, hba1c_on_cad)
hba1c_on_cad_plot


#### MR Hba1c - Empa plain ####

# remove IVs if they are also IVs in the second GWAS
hba1c_data_empa <- hba1c_data %>%
  filter(!(SNP %in% empa_plain_iv$SNP))

empa_plain_data_hba1c_ivs <- trial_arm_plain_data %>% filter(SNP %in% hba1c_data_empa$SNP)

hba1c_on_empa <- harmonise_data(hba1c_data_empa, empa_plain_data_hba1c_ivs)
hba1c_on_empa_plain_result <- mr(hba1c_on_empa)
hba1c_on_empa_plain_plot <- mr_scatter_plot(hba1c_on_empa_plain_result, hba1c_on_empa)
hba1c_on_empa_plain_plot


#### MR Hba1c - Empa ####

# remove IVs if they are also IVs in the second GWAS
hba1c_data_empa <- hba1c_data %>%
  filter(!(SNP %in% empa_iv$SNP))

empa_data_hba1c_ivs <- trial_arm_data %>% filter(SNP %in% hba1c_data_empa$SNP)

hba1c_on_empa <- harmonise_data(hba1c_data_empa, empa_data_hba1c_ivs)
hba1c_on_empa_result <- mr(hba1c_on_empa)
hba1c_on_empa_plot <- mr_scatter_plot(hba1c_on_empa_result, hba1c_on_empa)
hba1c_on_empa_plot



#### MR ALT - CAD ####
alt_gwas <- fread(lifted_alt_ss_)
alt_iv <- fread(alt_iv_)[,c(1,4,3)]

alt_data <- data.frame(
  id.exposure = "ALT",
  exposure = "ALT (PAN-UKBB)",
  SNP = alt_gwas$rsid,
  beta.exposure = alt_gwas$effect_weight,
  se.exposure = alt_gwas$se,
  effect_allele.exposure = alt_gwas$effect_allele,
  other_allele.exposure = alt_gwas$other_allele
) %>% 
  mutate(chr.exposure = alt_gwas$chromosome,
         pos.exposure = alt_gwas$chr_position,
         pval.exposure = alt_gwas$pval,
         samplesize.exposure = alt_gwas$samplesize,
         eaf.exposure = alt_gwas$af,
         mr_keep.exposure = TRUE,
         pval_origin.exposure = "reported",
         data_source.exposure = "") %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.exposure = as.integer(SNP %in% alt_iv$SNP)) %>%
  # keep only IVs
  filter(iv.exposure == 1) 

# remove IVs if they are also IVs in the second GWAS
alt_data_cad <- alt_data %>%
  filter(!(SNP %in% cad_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
cad_data_alt_ivs <- trial_outcome_data %>% filter(SNP %in% alt_data_cad$SNP)


alt_on_cad <- harmonise_data(alt_data_cad, cad_data_alt_ivs)
alt_on_cad_result <- mr(alt_on_cad)
alt_on_cad_plot <- mr_scatter_plot(alt_on_cad_result, alt_on_cad)
alt_on_cad_plot


#### MR ALT - Empa plain ####

# remove IVs if they are also IVs in the second GWAS
alt_data_empa <- alt_data %>%
  filter(!(SNP %in% empa_plain_iv$SNP))

empa_plain_data_alt_ivs <- trial_arm_plain_data %>% filter(SNP %in% alt_data_empa$SNP)

alt_on_empa <- harmonise_data(alt_data_empa, empa_plain_data_alt_ivs)
alt_on_empa_plain_result <- mr(alt_on_empa)
alt_on_empa_plain_plot <- mr_scatter_plot(alt_on_empa_plain_result, alt_on_empa)
alt_on_empa_plain_plot


#### MR ALT - Empa ####

# remove IVs if they are also IVs in the second GWAS
alt_data_empa <- alt_data %>%
  filter(!(SNP %in% empa_iv$SNP))

empa_data_alt_ivs <- trial_arm_data %>% filter(SNP %in% alt_data_empa$SNP)

alt_on_empa <- harmonise_data(alt_data_empa, empa_data_alt_ivs)
alt_on_empa_result <- mr(alt_on_empa)
alt_on_empa_plot <- mr_scatter_plot(alt_on_empa_result, alt_on_empa)
alt_on_empa_plot




#### MR HF - CAD ####
hf_gwas <- fread(lifted_hf_ss_)
hf_iv <- fread(hf_iv_)[,c(1,4,3)]

hf_data <- data.frame(
  id.exposure = "hf",
  exposure = "HF (HERMES 2019)",
  SNP = hf_gwas$rsid,
  beta.exposure = hf_gwas$effect_weight,
  se.exposure = hf_gwas$se,
  effect_allele.exposure = hf_gwas$effect_allele,
  other_allele.exposure = hf_gwas$other_allele
) %>% 
  mutate(chr.exposure = hf_gwas$chromosome,
         pos.exposure = hf_gwas$chr_position,
         pval.exposure = hf_gwas$pval,
         samplesize.exposure = hf_gwas$samplesize,
         eaf.exposure = hf_gwas$af,
         mr_keep.exposure = TRUE,
         pval_origin.exposure = "reported",
         data_source.exposure = "") %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.exposure = as.integer(SNP %in% hf_iv$SNP)) %>%
  # keep only IVs
  filter(iv.exposure == 1) 

# remove IVs if they are also IVs in the second GWAS
hf_data_cad <- hf_data %>%
  filter(!(SNP %in% cad_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
cad_data_hf_ivs <- trial_outcome_data %>% filter(SNP %in% hf_data_cad$SNP)


hf_on_cad <- harmonise_data(hf_data_cad, cad_data_hf_ivs)
hf_on_cad_result <- mr(hf_on_cad)
hf_on_cad_plot <- mr_scatter_plot(hf_on_cad_result, hf_on_cad)
hf_on_cad_plot


#### MR HF - Empa plain ####

# remove IVs if they are also IVs in the second GWAS
hf_data_empa <- hf_data %>%
  filter(!(SNP %in% empa_plain_iv$SNP))

empa_plain_data_hf_ivs <- trial_arm_plain_data %>% filter(SNP %in% hf_data_empa$SNP)

hf_on_empa <- harmonise_data(hf_data_empa, empa_plain_data_hf_ivs)
hf_on_empa_plain_result <- mr(hf_on_empa)
hf_on_empa_plain_plot <- mr_scatter_plot(hf_on_empa_plain_result, hf_on_empa)
hf_on_empa_plain_plot


#### MR HF - Empa ####

# remove IVs if they are also IVs in the second GWAS
hf_data_empa <- hf_data %>%
  filter(!(SNP %in% empa_iv$SNP))

empa_data_hf_ivs <- trial_arm_data %>% filter(SNP %in% hf_data_empa$SNP)

hf_on_empa <- harmonise_data(hf_data_empa, empa_data_hf_ivs)
hf_on_empa_result <- mr(hf_on_empa)
hf_on_empa_plot <- mr_scatter_plot(hf_on_empa_result, hf_on_empa)
hf_on_empa_plot




#### MR TG2 - CAD ####
tg2_gwas <- fread(tg2_gwas_)
tg2_iv <- fread(tg2_iv_)[,c(1,4,3)]

tg2_data <- data.frame(
  id.exposure = "TG",
  exposure = "TG (PAN-UKBB)",
  SNP = tg2_gwas$rsid,
  beta.exposure = tg2_gwas$effect_weight,
  se.exposure = tg2_gwas$se,
  effect_allele.exposure = tg2_gwas$effect_allele,
  other_allele.exposure = tg2_gwas$other_allele
) %>% 
  mutate(chr.exposure = tg2_gwas$chromosome,
         pos.exposure = tg2_gwas$chr_position,
         pval.exposure = tg2_gwas$pval,
         samplesize.exposure = tg2_gwas$samplesize,
         eaf.exposure = tg2_gwas$af,
         mr_keep.exposure = TRUE,
         pval_origin.exposure = "reported",
         data_source.exposure = "") %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.exposure = as.integer(SNP %in% tg2_iv$SNP)) %>%
  # keep only IVs
  filter(iv.exposure == 1) 

# remove IVs if they are also IVs in the second GWAS
tg2_data_cad <- tg2_data %>%
  filter(!(SNP %in% cad_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
cad_data_tg2_ivs <- trial_outcome_data %>% filter(SNP %in% tg2_data_cad$SNP)


tg2_on_cad <- harmonise_data(tg2_data_cad, cad_data_tg2_ivs)
tg2_on_cad_result <- mr(tg2_on_cad)
tg2_on_cad_plot <- mr_scatter_plot(tg2_on_cad_result, tg2_on_cad)
tg2_on_cad_plot


#### MR TG2 - Empa plain ####

# remove IVs if they are also IVs in the second GWAS
tg2_data_empa <- tg2_data %>%
  filter(!(SNP %in% empa_plain_iv$SNP))

empa_plain_data_tg2_ivs <- trial_arm_plain_data %>% filter(SNP %in% tg2_data_empa$SNP)

tg2_on_empa <- harmonise_data(tg2_data_empa, empa_plain_data_tg2_ivs)
tg2_on_empa_plain_result <- mr(tg2_on_empa)
tg2_on_empa_plain_plot <- mr_scatter_plot(tg2_on_empa_plain_result, tg2_on_empa)
tg2_on_empa_plain_plot


#### MR TG2 - Empa ####

# remove IVs if they are also IVs in the second GWAS
tg2_data_empa <- tg2_data %>%
  filter(!(SNP %in% empa_iv$SNP))

empa_data_tg2_ivs <- trial_arm_data %>% filter(SNP %in% tg2_data_empa$SNP)

tg2_on_empa <- harmonise_data(tg2_data_empa, empa_data_tg2_ivs)
tg2_on_empa_result <- mr(tg2_on_empa)
tg2_on_empa_plot <- mr_scatter_plot(tg2_on_empa_result, tg2_on_empa)
tg2_on_empa_plot




#### MR SAH - CAD ####
sah_gwas <- fread(lifted_sah_ss_)
sah_iv <- fread(sah_iv_)[,c(1,4,3)]

sah_data <- data.frame(
  id.exposure = "Subarach. Hemo.",
  exposure = "SAH (Bakker 2020)",
  SNP = sah_gwas$rsid,
  beta.exposure = sah_gwas$effect_weight,
  se.exposure = sah_gwas$se,
  effect_allele.exposure = sah_gwas$effect_allele,
  other_allele.exposure = sah_gwas$other_allele
) %>% 
  mutate(chr.exposure = sah_gwas$chromosome,
         pos.exposure = sah_gwas$chr_position,
         pval.exposure = sah_gwas$pval,
         samplesize.exposure = sah_gwas$samplesize,
         eaf.exposure = sah_gwas$af,
         mr_keep.exposure = TRUE,
         pval_origin.exposure = "reported",
         data_source.exposure = "") %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.exposure = as.integer(SNP %in% sah_iv$SNP)) %>%
  # keep only IVs
  filter(iv.exposure == 1) 

# remove IVs if they are also IVs in the second GWAS
sah_data_cad <- sah_data %>%
  filter(!(SNP %in% cad_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
cad_data_sah_ivs <- trial_outcome_data %>% filter(SNP %in% sah_data_cad$SNP)


sah_on_cad <- harmonise_data(sah_data_cad, cad_data_sah_ivs)
sah_on_cad_result <- mr(sah_on_cad)
sah_on_cad_plot <- mr_scatter_plot(sah_on_cad_result, sah_on_cad)
sah_on_cad_plot


#### MR SAH - Empa plain ####

# remove IVs if they are also IVs in the second GWAS
sah_data_empa <- sah_data %>%
  filter(!(SNP %in% empa_plain_iv$SNP))

empa_plain_data_sah_ivs <- trial_arm_plain_data %>% filter(SNP %in% sah_data_empa$SNP)

sah_on_empa <- harmonise_data(sah_data_empa, empa_plain_data_sah_ivs)
sah_on_empa_plain_result <- mr(sah_on_empa)
sah_on_empa_plain_plot <- mr_scatter_plot(sah_on_empa_plain_result, sah_on_empa)
sah_on_empa_plain_plot


#### MR SAH - Empa ####

# remove IVs if they are also IVs in the second GWAS
sah_data_empa <- sah_data %>%
  filter(!(SNP %in% empa_iv$SNP))

empa_data_sah_ivs <- trial_arm_data %>% filter(SNP %in% sah_data_empa$SNP)

sah_on_empa <- harmonise_data(sah_data_empa, empa_data_sah_ivs)
sah_on_empa_result <- mr(sah_on_empa)
sah_on_empa_plot <- mr_scatter_plot(sah_on_empa_result, sah_on_empa)
sah_on_empa_plot




#### MR SBP - CAD ####
sbp_gwas <- fread(lifted_sbp_ss_)
sbp_iv <- fread(sbp_iv_)[,c(1,4,3)]

sbp_data <- data.frame(
  id.exposure = "Systolic Blood Pressure",
  exposure = "SBP (UKBB ICBP)",
  SNP = sbp_gwas$rsid,
  beta.exposure = sbp_gwas$effect_weight,
  se.exposure = sbp_gwas$se,
  effect_allele.exposure = sbp_gwas$effect_allele,
  other_allele.exposure = sbp_gwas$other_allele
) %>% 
  mutate(chr.exposure = sbp_gwas$chromosome,
         pos.exposure = sbp_gwas$chr_position,
         pval.exposure = sbp_gwas$pval,
         samplesize.exposure = sbp_gwas$samplesize,
         eaf.exposure = sbp_gwas$af,
         mr_keep.exposure = TRUE,
         pval_origin.exposure = "reported",
         data_source.exposure = "") %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.exposure = as.integer(SNP %in% sbp_iv$SNP)) %>%
  # keep only IVs
  filter(iv.exposure == 1) 

# remove IVs if they are also IVs in the second GWAS
sbp_data_cad <- sbp_data %>%
  filter(!(SNP %in% cad_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
cad_data_sbp_ivs <- trial_outcome_data %>% filter(SNP %in% sbp_data_cad$SNP)


sbp_on_cad <- harmonise_data(sbp_data_cad, cad_data_sbp_ivs)
sbp_on_cad_result <- mr(sbp_on_cad)
sbp_on_cad_plot <- mr_scatter_plot(sbp_on_cad_result, sbp_on_cad)
sbp_on_cad_plot


#### MR SBP - Empa plain ####

# remove IVs if they are also IVs in the second GWAS
sbp_data_empa <- sbp_data %>%
  filter(!(SNP %in% empa_plain_iv$SNP))

empa_plain_data_sbp_ivs <- trial_arm_plain_data %>% filter(SNP %in% sbp_data_empa$SNP)

sbp_on_empa <- harmonise_data(sbp_data_empa, empa_plain_data_sbp_ivs)
sbp_on_empa_plain_result <- mr(sbp_on_empa)
sbp_on_empa_plain_plot <- mr_scatter_plot(sbp_on_empa_plain_result, sbp_on_empa)
sbp_on_empa_plain_plot


#### MR SBP - Empa ####

# remove IVs if they are also IVs in the second GWAS
sbp_data_empa <- sbp_data %>%
  filter(!(SNP %in% empa_iv$SNP))

empa_data_sbp_ivs <- trial_arm_data %>% filter(SNP %in% sbp_data_empa$SNP)

sbp_on_empa <- harmonise_data(sbp_data_empa, empa_data_sbp_ivs)
sbp_on_empa_result <- mr(sbp_on_empa)
sbp_on_empa_plot <- mr_scatter_plot(sbp_on_empa_result, sbp_on_empa)
sbp_on_empa_plot




#### MR LIV - CAD ####
liv_gwas <- fread(write_liv_ss_)
liv_iv <- fread(liv_iv_)[,c(1,4,3)]

liv_data <- data.frame(
  id.exposure = "Liver Disease",
  exposure = "Liver Disease (Fairfield 2022)",
  SNP = liv_gwas$rsid,
  beta.exposure = liv_gwas$effect_weight,
  se.exposure = liv_gwas$se,
  effect_allele.exposure = liv_gwas$effect_allele,
  other_allele.exposure = liv_gwas$other_allele
) %>% 
  mutate(chr.exposure = liv_gwas$chromosome,
         pos.exposure = liv_gwas$chr_position,
         pval.exposure = liv_gwas$pval,
         samplesize.exposure = liv_gwas$samplesize,
         eaf.exposure = liv_gwas$af,
         mr_keep.exposure = TRUE,
         pval_origin.exposure = "reported",
         data_source.exposure = "") %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.exposure = as.integer(SNP %in% liv_iv$SNP)) %>%
  # keep only IVs
  filter(iv.exposure == 1) 

# remove IVs if they are also IVs in the second GWAS
liv_data_cad <- liv_data %>%
  filter(!(SNP %in% cad_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
cad_data_liv_ivs <- trial_outcome_data %>% filter(SNP %in% liv_data_cad$SNP)


liv_on_cad <- harmonise_data(liv_data_cad, cad_data_liv_ivs)
liv_on_cad_result <- mr(liv_on_cad)
liv_on_cad_plot <- mr_scatter_plot(liv_on_cad_result, liv_on_cad)
liv_on_cad_plot


#### MR LIV - Empa plain ####

# remove IVs if they are also IVs in the second GWAS
liv_data_empa <- liv_data %>%
  filter(!(SNP %in% empa_plain_iv$SNP))

empa_plain_data_liv_ivs <- trial_arm_plain_data %>% filter(SNP %in% liv_data_empa$SNP)

liv_on_empa <- harmonise_data(liv_data_empa, empa_plain_data_liv_ivs)
liv_on_empa_plain_result <- mr(liv_on_empa)
liv_on_empa_plain_plot <- mr_scatter_plot(liv_on_empa_plain_result, liv_on_empa)
liv_on_empa_plain_plot


#### MR LIV - Empa ####

# remove IVs if they are also IVs in the second GWAS
liv_data_empa <- liv_data %>%
  filter(!(SNP %in% empa_iv$SNP))

empa_data_liv_ivs <- trial_arm_data %>% filter(SNP %in% liv_data_empa$SNP)

liv_on_empa <- harmonise_data(liv_data_empa, empa_data_liv_ivs)
liv_on_empa_result <- mr(liv_on_empa)
liv_on_empa_plot <- mr_scatter_plot(liv_on_empa_result, liv_on_empa)
liv_on_empa_plot




#### MR Stroke - CAD ####
stroke_gwas <- fread(lifted_stroke_ss_)
stroke_iv <- fread(stroke_iv_)[,c(1,4,3)]

stroke_data <- data.frame(
  id.exposure = "Stroke",
  exposure = "Stroke (MEGASTROKE 2018)",
  SNP = stroke_gwas$rsid,
  beta.exposure = stroke_gwas$effect_weight,
  se.exposure = stroke_gwas$se,
  effect_allele.exposure = stroke_gwas$effect_allele,
  other_allele.exposure = stroke_gwas$other_allele
) %>% 
  mutate(chr.exposure = stroke_gwas$chromosome,
         pos.exposure = stroke_gwas$chr_position,
         pval.exposure = stroke_gwas$pval,
         samplesize.exposure = stroke_gwas$samplesize,
         eaf.exposure = stroke_gwas$af,
         mr_keep.exposure = TRUE,
         pval_origin.exposure = "reported",
         data_source.exposure = "") %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.exposure = as.integer(SNP %in% stroke_iv$SNP)) %>%
  # keep only IVs
  filter(iv.exposure == 1) 

# remove IVs if they are also IVs in the second GWAS
stroke_data_cad <- stroke_data %>%
  filter(!(SNP %in% cad_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
cad_data_stroke_ivs <- trial_outcome_data %>% filter(SNP %in% stroke_data_cad$SNP)


stroke_on_cad <- harmonise_data(stroke_data_cad, cad_data_stroke_ivs)
stroke_on_cad_result <- mr(stroke_on_cad)
stroke_on_cad_plot <- mr_scatter_plot(stroke_on_cad_result, stroke_on_cad)
stroke_on_cad_plot


#### MR Stroke - Empa plain ####

# remove IVs if they are also IVs in the second GWAS
stroke_data_empa <- stroke_data %>%
  filter(!(SNP %in% empa_plain_iv$SNP))

empa_plain_data_stroke_ivs <- trial_arm_plain_data %>% filter(SNP %in% stroke_data_empa$SNP)

stroke_on_empa <- harmonise_data(stroke_data_empa, empa_plain_data_stroke_ivs)
stroke_on_empa_plain_result <- mr(stroke_on_empa)
stroke_on_empa_plain_plot <- mr_scatter_plot(stroke_on_empa_plain_result, stroke_on_empa)
stroke_on_empa_plain_plot


#### MR Stroke - Empa ####

# remove IVs if they are also IVs in the second GWAS
stroke_data_empa <- stroke_data %>%
  filter(!(SNP %in% empa_iv$SNP))

empa_data_stroke_ivs <- trial_arm_data %>% filter(SNP %in% stroke_data_empa$SNP)

stroke_on_empa <- harmonise_data(stroke_data_empa, empa_data_stroke_ivs)
stroke_on_empa_result <- mr(stroke_on_empa)
stroke_on_empa_plot <- mr_scatter_plot(stroke_on_empa_result, stroke_on_empa)
stroke_on_empa_plot




#### MR CRP - CAD ####
crp_gwas <- fread(lifted_crp_ss_)
crp_iv <- fread(crp_iv_)[,c(1,4,3)]

crp_data <- data.frame(
  id.exposure = "CRP",
  exposure = "CRP (PAN-UKBB)",
  SNP = crp_gwas$rsid,
  beta.exposure = crp_gwas$effect_weight,
  se.exposure = crp_gwas$se,
  effect_allele.exposure = crp_gwas$effect_allele,
  other_allele.exposure = crp_gwas$other_allele
) %>% 
  mutate(chr.exposure = crp_gwas$chromosome,
         pos.exposure = crp_gwas$chr_position,
         pval.exposure = crp_gwas$pval,
         samplesize.exposure = crp_gwas$samplesize,
         eaf.exposure = crp_gwas$af,
         mr_keep.exposure = TRUE,
         pval_origin.exposure = "reported",
         data_source.exposure = "") %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.exposure = as.integer(SNP %in% crp_iv$SNP)) %>%
  # keep only IVs
  filter(iv.exposure == 1) 

# remove IVs if they are also IVs in the second GWAS
crp_data_cad <- crp_data %>%
  filter(!(SNP %in% cad_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
cad_data_crp_ivs <- trial_outcome_data %>% filter(SNP %in% crp_data_cad$SNP)


crp_on_cad <- harmonise_data(crp_data_cad, cad_data_crp_ivs)
crp_on_cad_result <- mr(crp_on_cad)
crp_on_cad_plot <- mr_scatter_plot(crp_on_cad_result, crp_on_cad)
crp_on_cad_plot


#### MR CRP - Empa plain ####

# remove IVs if they are also IVs in the second GWAS
crp_data_empa <- crp_data %>%
  filter(!(SNP %in% empa_plain_iv$SNP))

empa_plain_data_crp_ivs <- trial_arm_plain_data %>% filter(SNP %in% crp_data_empa$SNP)

crp_on_empa <- harmonise_data(crp_data_empa, empa_plain_data_crp_ivs)
crp_on_empa_plain_result <- mr(crp_on_empa)
crp_on_empa_plain_plot <- mr_scatter_plot(crp_on_empa_plain_result, crp_on_empa)
crp_on_empa_plain_plot


#### MR CRP - Empa ####

# remove IVs if they are also IVs in the second GWAS
crp_data_empa <- crp_data %>%
  filter(!(SNP %in% empa_iv$SNP))

empa_data_crp_ivs <- trial_arm_data %>% filter(SNP %in% crp_data_empa$SNP)

crp_on_empa <- harmonise_data(crp_data_empa, empa_data_crp_ivs)
crp_on_empa_result <- mr(crp_on_empa)
crp_on_empa_plot <- mr_scatter_plot(crp_on_empa_result, crp_on_empa)
crp_on_empa_plot




#### MR AF - CAD ####
af_gwas <- fread(lifted_af_ss_)
af_iv <- fread(af_iv_)[,c(1,4,3)]

af_data <- data.frame(
  id.exposure = "Atrial Fibrillation",
  exposure = "Atrial Fibrillation (Roselli 2018)",
  SNP = af_gwas$rsid,
  beta.exposure = af_gwas$effect_weight,
  se.exposure = af_gwas$se,
  effect_allele.exposure = af_gwas$effect_allele,
  other_allele.exposure = af_gwas$other_allele
) %>% 
  mutate(chr.exposure = af_gwas$chromosome,
         pos.exposure = af_gwas$chr_position,
         pval.exposure = af_gwas$pval,
         samplesize.exposure = af_gwas$samplesize,
         eaf.exposure = af_gwas$af,
         mr_keep.exposure = TRUE,
         pval_origin.exposure = "reported",
         data_source.exposure = "") %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.exposure = as.integer(SNP %in% af_iv$SNP)) %>%
  # keep only IVs
  filter(iv.exposure == 1) 

# remove IVs if they are also IVs in the second GWAS
af_data_cad <- af_data %>%
  filter(!(SNP %in% cad_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
cad_data_af_ivs <- trial_outcome_data %>% filter(SNP %in% af_data_cad$SNP)


af_on_cad <- harmonise_data(af_data_cad, cad_data_af_ivs)
af_on_cad_result <- mr(af_on_cad)
af_on_cad_plot <- mr_scatter_plot(af_on_cad_result, af_on_cad)
af_on_cad_plot


#### MR AF - Empa plain ####

# remove IVs if they are also IVs in the second GWAS
af_data_empa <- af_data %>%
  filter(!(SNP %in% empa_plain_iv$SNP))

empa_plain_data_af_ivs <- trial_arm_plain_data %>% filter(SNP %in% af_data_empa$SNP)

af_on_empa <- harmonise_data(af_data_empa, empa_plain_data_af_ivs)
af_on_empa_plain_result <- mr(af_on_empa)
af_on_empa_plain_plot <- mr_scatter_plot(af_on_empa_plain_result, af_on_empa)
af_on_empa_plain_plot


#### MR AF - Empa ####

# remove IVs if they are also IVs in the second GWAS
af_data_empa <- af_data %>%
  filter(!(SNP %in% empa_iv$SNP))

empa_data_af_ivs <- trial_arm_data %>% filter(SNP %in% af_data_empa$SNP)

af_on_empa <- harmonise_data(af_data_empa, empa_data_af_ivs)
af_on_empa_result <- mr(af_on_empa)
af_on_empa_plot <- mr_scatter_plot(af_on_empa_result, af_on_empa)
af_on_empa_plot




#### MR AST - CAD ####
ast_gwas <- fread(lifted_ast_ss_)
ast_iv <- fread(ast_iv_)[,c(1,4,3)]

ast_data <- data.frame(
  id.exposure = "AST",
  exposure = "AST (PAN-UKBB)",
  SNP = ast_gwas$rsid,
  beta.exposure = ast_gwas$effect_weight,
  se.exposure = ast_gwas$se,
  effect_allele.exposure = ast_gwas$effect_allele,
  other_allele.exposure = ast_gwas$other_allele
) %>% 
  mutate(chr.exposure = ast_gwas$chromosome,
         pos.exposure = ast_gwas$chr_position,
         pval.exposure = ast_gwas$pval,
         samplesize.exposure = ast_gwas$samplesize,
         eaf.exposure = ast_gwas$af,
         mr_keep.exposure = TRUE,
         pval_origin.exposure = "reported",
         data_source.exposure = "") %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.exposure = as.integer(SNP %in% ast_iv$SNP)) %>%
  # keep only IVs
  filter(iv.exposure == 1) 

# remove IVs if they are also IVs in the second GWAS
ast_data_cad <- ast_data %>%
  filter(!(SNP %in% cad_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
cad_data_ast_ivs <- trial_outcome_data %>% filter(SNP %in% ast_data_cad$SNP)


ast_on_cad <- harmonise_data(ast_data_cad, cad_data_ast_ivs)
ast_on_cad_result <- mr(ast_on_cad)
ast_on_cad_plot <- mr_scatter_plot(ast_on_cad_result, ast_on_cad)
ast_on_cad_plot


#### MR AST - Empa plain ####

# remove IVs if they are also IVs in the second GWAS
ast_data_empa <- ast_data %>%
  filter(!(SNP %in% empa_plain_iv$SNP))

empa_plain_data_ast_ivs <- trial_arm_plain_data %>% filter(SNP %in% ast_data_empa$SNP)

ast_on_empa <- harmonise_data(ast_data_empa, empa_plain_data_ast_ivs)
ast_on_empa_plain_result <- mr(ast_on_empa)
ast_on_empa_plain_plot <- mr_scatter_plot(ast_on_empa_plain_result, ast_on_empa)
ast_on_empa_plain_plot


#### MR AST - Empa ####

# remove IVs if they are also IVs in the second GWAS
ast_data_empa <- ast_data %>%
  filter(!(SNP %in% empa_iv$SNP))

empa_data_ast_ivs <- trial_arm_data %>% filter(SNP %in% ast_data_empa$SNP)

ast_on_empa <- harmonise_data(ast_data_empa, empa_data_ast_ivs)
ast_on_empa_result <- mr(ast_on_empa)
ast_on_empa_plot <- mr_scatter_plot(ast_on_empa_result, ast_on_empa)
ast_on_empa_plot





# Save results ------------------------------------------------------------

# Number of MR results
ls()[grep("_result$", ls())]
# bind them together
results_all_mr <- do.call(rbind, lapply(ls()[grep("_result$", ls())], get))
fwrite(results_all_mr, results_all_mr_)

results_all_mr <- results_all_mr %>% 
  select("Exposure" = "id.exposure",
         "Outcome" = "id.outcome",
         "Method" = "method",
         nsnp, b, se, pval)

fwrite(results_all_mr, results_all_mr_redCol_)




# Save MR plots -----------------------------------------------------------

#### Each separately ####

main_path <- main_path_mr_plots_

# Loop over each plot
for (plot_name in ls(pattern = "_plot$")) {
  # Remove "_plot" from the object name to get the file name
  file_name <- sub("_plot$", "", plot_name)
  
  # Get the dataframe object
  df <- get(plot_name)
  
  # Save as pdf
  pdf_file <- paste0(main_path, file_name, ".pdf")
  pdf(pdf_file,  width = 16, height = 10)
  print(df)
  dev.off()
  
  # Save as png
  png_file <- paste0(main_path, file_name, ".png")
  png(png_file,  width = 16*300, height = 10*300, units = "px", res = 300)
  print(df)
  dev.off()
  
  # Save as svg
  svg_file <- paste0(main_path, file_name, ".svg")
  svg(svg_file,  width = 16, height = 10)
  print(df)
  dev.off()
  
  # Print confirmation
  cat("Plots saved for", file_name, "\n")
}


#### All together ####
# Get a list of all objects in global environment ending with "_plot"
plot_objects <- ls(pattern = "_plot$")

# Extract prefixes and suffixes
file_info <- strsplit(plot_objects, "_")
prefixes <- sapply(file_info, function(x) x[1])
suffixes <- sapply(file_info, function(x) paste(x[-c(1:2)], collapse = "_"))

# Unique prefixes and suffixes
unique_prefixes <- unique(prefixes)
unique_suffixes <- unique(suffixes)

pdf_file <- paste0(main_path, "all_mr_plots.pdf")
png_file <- paste0(main_path, "all_mr_plots.png")
svg_file <- paste0(main_path, "all_mr_plots.svg")

all_plots_grid <- NULL

for (prefix in unique_prefixes) {
  for (suffix in unique_suffixes) {
    # Find corresponding plot
    plot_name <- paste0(prefix, "_on_", suffix)
    if (plot_name %in% plot_objects) {
      # Get the existing plot object
      plot_obj <- get(plot_name)
      
      # Add the plot to the grid
      all_plots_grid <- rbind(all_plots_grid, plot_obj)
      
      cat("Plot for", plot_name, "\n")
    }
  }
}

pdf(pdf_file, width = 16*3, height = 10*19)
grid.arrange(grobs = all_plots_grid, ncol = length(unique_suffixes))
dev.off()

# no png -> file would be way to big

svg(svg_file, width = 16*3, height = 10*19)
grid.arrange(grobs = all_plots_grid, ncol = length(unique_suffixes))
dev.off()

#### without Liver ####
# There is no plot for Liver Disease because of too few SNPs
# Remove "liv" from unique_prefixes
unique_prefixes <- unique_prefixes[unique_prefixes != "liv"]

pdf_file <- paste0(main_path, "all_mr_plots_noLiv.pdf")
svg_file <- paste0(main_path, "all_mr_plots_noLiv.svg")

all_plots_grid <- NULL

for (prefix in unique_prefixes) {
  for (suffix in unique_suffixes) {
    # Find corresponding plot
    plot_name <- paste0(prefix, "_on_", suffix)
    if (plot_name %in% plot_objects) {
      # Get the existing plot object
      plot_obj <- get(plot_name)
      
      # Add the plot to the grid
      all_plots_grid <- rbind(all_plots_grid, plot_obj)
      
      cat("Plot for", plot_name, "\n")
    }
  }
}

pdf(pdf_file, width = 16*3, height = 10*18)
grid.arrange(grobs = all_plots_grid, ncol = length(unique_suffixes))
dev.off()

# no png -> file would be way to big

svg(svg_file, width = 16*3, height = 10*18)
grid.arrange(grobs = all_plots_grid, ncol = length(unique_suffixes))
dev.off()




###########################################################################
###########################################################################
###########################################################################
###########################################################################
# When starting from saved results ----------------------------------------

results_all_mr <- fread(results_all_mr_)

# split dataframe in individual results dataframes based on combinations of id.exposure and id.outcome
group_var <- tolower(gsub(" ", "_", paste(results_all_mr$id.exposure, results_all_mr$id.outcome, sep = "_")))

# Iterate over the list of dataframes and save them individually
split_results_all_mr <- split(results_all_mr, group_var)
for (i in seq_along(split_results_all_mr)) {
  # Generate a unique name for the dataframe
  df_name <- paste("df_", names(split_results_all_mr)[i], sep = "")
  # Assign the dataframe to the generated name
  assign(df_name, split_results_all_mr[[i]])
}

results_all_mr <- results_all_mr %>% 
  select("Exposure" = "id.exposure",
         "Outcome" = "id.outcome",
         "Method" = "method",
         nsnp, b, se, pval)


# Operate an all MR results -----------------------------------------------

results_all_mr <- results_all_mr %>% 
  filter(Exposure != "Liver Disease") %>%  # remove Liver Disease because only one SNP after clumping
  mutate(ci_lower = b - 1.96*se,
         ci_upper = b + 1.96*se,
         # create rank 
         rank = case_when(Exposure == "BMI" ~ 18,
                          Exposure == "T2D" ~ 17,
                          Exposure == "Major Depressive Disorder" ~ 16,
                          Exposure == "hf" ~ 15,
                          Exposure == "CKD" ~ 14,
                          Exposure == "ALT" ~ 13,
                          Exposure == "TG"  ~ 12,
                          Exposure == "Subarach. Hemo." ~ 11,
                          Exposure == "Systolic Blood Pressure" ~ 10,
                          Exposure == "Asthma" ~ 9,
                          Exposure == "AST" ~ 8,
                          Exposure == "Hba1c" ~ 7,
                          Exposure == "Stroke" ~ 6,
                          Exposure == "CRP" ~ 5,
                          Exposure == "LDL" ~ 4,
                          Exposure == "Atrial Fibrillation" ~ 3,
                          Exposure == "Educational Attainment" ~ 2,
                          Exposure == "HDL" ~ 1),
         # change Exposure var names
         Exposure = case_when(Exposure == "BMI" ~ "Body Mass Index",
                              Exposure == "T2D" ~ "Type 2 Diabetes",
                              Exposure == "Major Depressive Disorder" ~ "Major Depression",
                              Exposure == "hf" ~ "Heart Failure",
                              Exposure == "CKD" ~ "Chronic Kidney Disease",
                              Exposure == "ALT" ~ "ALT Enzyme",
                              Exposure == "TG"  ~ "Triglycerides",
                              Exposure == "Subarach. Hemo." ~ "Subarachnoid Hemorrhage",
                              Exposure == "Systolic Blood Pressure" ~ "Systolic Blood Pressure",
                              Exposure == "Asthma" ~ "Asthma",
                              Exposure == "AST" ~ "AST Enzyme",
                              Exposure == "Hba1c" ~ "Hba1c",
                              Exposure == "Stroke" ~ "Stroke",
                              Exposure == "CRP" ~ "CRP",
                              Exposure == "LDL" ~ "LDL",
                              Exposure == "Atrial Fibrillation" ~ "Atrial Fibrillation",
                              Exposure == "Educational Attainment" ~ "Education",
                              Exposure == "HDL" ~ "HDL")) %>% 
  arrange(rank) %>% 
  mutate(stripe = rep(c(rep(0, 15), rep(1, 15)), 9))

results_all_mr$Outcome <- factor(results_all_mr$Outcome, levels = c("CAD", "Empagliflozin plain", "Empagliflozin"))

results_all_mr_ivw <- results_all_mr %>% filter(Method == "Inverse variance weighted")
fwrite(results_all_mr_ivw, results_all_mr_ivw_redCol_)

results_all_mr_ivw <- results_all_mr_ivw %>% mutate(OR = exp(b),
                                                    ci_lower_OR = exp(ci_lower),
                                                    ci_upper_OR = exp(ci_upper)) 

#### Visualize all MR results ####


pd <- position_dodge(width=1)
min(results_all_mr_ivw$ci_lower)
max(results_all_mr_ivw$ci_upper)

# For IVW
ivw <- ggplot(results_all_mr_ivw, aes(x = OR, y = reorder(Exposure, rank), group = Method)) +
  geom_point(aes(color = Exposure), position = pd, size = 2) +
  geom_linerange(aes(xmin = ci_lower_OR, xmax = ci_upper_OR, color = Exposure), position = pd) +
  geom_vline(xintercept = 1, color = "black", linetype = "dashed") +
  geom_point(data = subset(results_all_mr_ivw, pval <= 0.05), pch =21, fill = NA, size = 4, colour = "black", stroke = 0.5) +
  geom_rect(aes(ymin = rank-0.5, ymax = rank+0.5, xmin = -Inf, xmax = Inf, fill = as.factor(stripe)), color= NA, alpha = 0.1) +
  theme_minimal() +
  scale_x_continuous(limits = c(0.4, 3.1), breaks = seq(0.5, 3.0, by = 0.5)) +
  facet_grid(.~Outcome, scales = "free_x", labeller = labeller(Outcome = c("CAD" = "B1: Coronary heart disease\n      (C --> Y)", "Empagliflozin plain" = "B2: Emplaglifozin\n      (C --> X)", "Empagliflozin" = "B3: Empagliflozin after elig. crit.\n      (C --> X)"))) +
  theme(strip.text = element_text(size = 15, hjust = 0)) +
  scale_fill_manual(values = c("white", "grey50"), name="", guide="none") +
  labs(y = "Genetically-predicted confounders (C)", 
       x = "Odds Ratio")+
  theme(axis.text.x = element_text(angle = 0, size = 15),    
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.spacing = unit(c(3, 1), "lines"),
        legend.position = "",
        legend.text = element_blank(),
        legend.title = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm"), 
        plot.title = element_text(hjust = 0.5, size = 20,),
        plot.subtitle = element_text(hjust = 0.5, vjust = -1, size = 15,),
  )



#### Save plot ####
pdf(ivw_pdf_, width = 16, height = 10)
ivw
dev.off()

png(ivw_png_, width = 16*300, height = 10*300, units = "px", res = 300)
ivw
dev.off()

svg(ivw_svg_, width = 16, height = 10)
ivw
dev.off()

# For all
all <- ggplot(results_all_mr, aes(x = b, y = reorder(Exposure, rank), group = Method)) +
  geom_point(aes(color = Exposure), position = pd, size = 2) +
  geom_linerange(aes(xmin = ci_lower, xmax = ci_upper, color = Exposure), position = pd) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(results_all_mr, pval <= 0.05),aes(x = b), pch =21, fill = NA, size = 4, colour = "black", stroke = 0.5, position = pd) +
  geom_rect(aes(ymin = rank-0.5, ymax = rank+0.5, xmin = -Inf, xmax = Inf, fill = as.factor(stripe)), color= NA, alpha = 0.1) +
  theme_minimal() +
  #xlim(c(-0.8, 1.15)) +
  facet_grid(Method~Outcome, scales = "free_x", labeller = labeller(Outcome = c("CAD" = "CAD", "Empagliflozin plain" = "Empagliflozin Initiation", "Empagliflozin" = "Empareg after eligibility criteria"))) +
  theme(strip.text = element_text(size = 15)) +
  scale_fill_manual(values = c("white", "grey50"), name="", guide="none") +
  ggtitle("Mendelian Randomization Results (all methods)") +
  labs(y = "Tested covariates as Exposure", x = "Beta")+
  theme(axis.text.x = element_text(angle = 0, size = 15),    
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.spacing = unit(2, "lines"),
        legend.position = "",
        legend.text = element_blank(),
        legend.title = element_blank(), 
        plot.margin = margin(1, 1, 1, 1, "cm"), 
        plot.title = element_text(hjust = 0.5, size = 20,),
  )



# MR Empa - CAD -----------------------------------------------------------
# only possible with Emmpa GWAS of init vs non init

# replace 'outcome' to exposure in colnames
empa_data <- trial_arm_data
colnames(empa_data) <- gsub('outcome', 'exposure', colnames(empa_data))

empa_data <- empa_data %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.exposure = as.integer(SNP %in% empa_iv$SNP)) %>%
  # keep only IVs
  filter(iv.exposure == 1) 

# remove IVs if they are also IVs in the second GWAS
empa_data_cad <- empa_data %>%
  filter(!(SNP %in% cad_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
cad_data_empa_ivs <- trial_outcome_data %>% filter(SNP %in% empa_data_cad$SNP)


empa_on_cad <- harmonise_data(empa_data_cad, cad_data_empa_ivs)
empa_on_cad_result <- mr(empa_on_cad)
empa_on_cad_plot <- mr_scatter_plot(empa_on_cad_result, empa_on_cad)
empa_on_cad_plot




# MR CAD - Empa -----------------------------------------------------------
# possible with both: Emmpa GWAS of init vs non init and Empareg GWAS after elig. crit.

# replace 'outcome' to exposure in colnames
cad_data <- trial_outcome_data
colnames(cad_data) <- gsub('outcome', 'exposure', colnames(cad_data))

cad_data <- cad_data %>% 
  # remove duplicate SNPs (from mapping of chrom_pos to rsID)
  filter(SNP != "" & !(base::duplicated(SNP))) %>% 
  # add a column indicating whether it is an index SNP or not
  mutate(iv.exposure = as.integer(SNP %in% cad_iv$SNP)) %>%
  # keep only IVs
  filter(iv.exposure == 1) 


# WITH EMPA AFTER ELIGIBILITY CRITERIA
# remove IVs if they are also IVs in the second GWAS
cad_data_empa <- cad_data %>%
  filter(!(SNP %in% empa_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
empa_data_cad_ivs <- trial_arm_data %>% filter(SNP %in% cad_data_empa$SNP)


cad_on_empa <- harmonise_data(cad_data_empa, empa_data_cad_ivs)
cad_on_empa_result <- mr(cad_on_empa)
cad_on_empa_plot <- mr_scatter_plot(cad_on_empa_result, cad_on_empa)
cad_on_empa_plot


# WITH EMPA INITIATORS
# remove IVs if they are also IVs in the second GWAS
cad_data_empa_plain <- cad_data %>%
  filter(!(SNP %in% empa_plain_iv$SNP))

# FILTER the IVs from Trial Outcome GWAS
empa_plain_data_cad_ivs <- trial_arm_plain_data %>% filter(SNP %in% cad_data_empa_plain$SNP)


cad_on_empa_plain <- harmonise_data(cad_data_empa_plain, empa_plain_data_cad_ivs)
cad_on_empa_plain_result <- mr(cad_on_empa_plain)
cad_on_empa_plain_plot <- mr_scatter_plot(cad_on_empa_plain_result, cad_on_empa_plain)
cad_on_empa_plain_plot


# Number of MR results
ls()[grep("_result$", ls())]
# bind them together
results_cad_empa_mr <- do.call(rbind, lapply(ls()[grep("^cad.*_result$", ls())], get))
fwrite(results_cad_empa_mr, results_cad_empa_mr_)

results_cad_empa_mr <- results_cad_empa_mr %>% 
  select("Exposure" = "id.exposure",
         "Outcome" = "id.outcome",
         "Method" = "method",
         nsnp, b, se, pval)
fwrite(results_cad_empa_mr, results_cad_empa_mr_redCol_)

results_cad_empa_mr <- results_cad_empa_mr %>% 
  mutate(ci_lower = b - 1.96*se,
         ci_upper = b + 1.96*se)
results_cad_empa_mr_ivw <- results_cad_empa_mr %>% filter(Method == "Inverse variance weighted")


pd <- position_dodge(width=1)
ivw <- ggplot(results_cad_empa_mr_ivw, aes(x = b, y = Exposure, group = Method)) +
  geom_point(aes(color = Exposure), position = pd, size = 2) +
  geom_linerange(aes(xmin = ci_lower, xmax = ci_upper, color = Exposure), position = pd) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(results_cad_empa_mr_ivw, pval <= 0.05), pch =21, fill = NA, size = 4, colour = "black", stroke = 0.5) +
  theme_minimal() +
  xlim(c(-0.8, 1.15)) +
  facet_grid(.~Outcome, scales = "free_x", labeller = labeller(Outcome = c("Empagliflozin plain" = "Empagliflozin Initiation (X)", "Empagliflozin" = "Empareg after elig. crit. (X)"))) +
  theme(strip.text = element_text(size = 15)) +
  scale_fill_manual(values = c("white", "grey50"), name="", guide="none") +
  labs(y = "Exposures (C)", 
    x = "Beta") +
  theme(axis.text.x = element_text(angle = 0, size = 15),    
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.spacing = unit(2, "lines"),
        legend.position = "",
        legend.text = element_blank(),
        legend.title = element_blank(), # element_text(size = 11),
        plot.margin = margin(1, 1, 1, 1, "cm"), 
        plot.title = element_text(hjust = 0.5, size = 20,),
        plot.subtitle = element_text(hjust = 0.5, vjust = -1, size = 15,),
  )



all <- ggplot(results_cad_empa_mr, aes(x = b, y = Exposure, group = Method)) +
  geom_point(aes(color = Exposure), position = pd, size = 2) +
  geom_linerange(aes(xmin = ci_lower, xmax = ci_upper, color = Exposure), position = pd) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  geom_point(data = subset(results_cad_empa_mr, pval <= 0.05),aes(x = b), pch =21, fill = NA, size = 4, colour = "black", stroke = 0.5, position = pd) +
  geom_text(aes(label = round(b, 2)), vjust = -1, hjust = 1)+
  theme_minimal() +
  #xlim(c(-0.8, 1.15)) +
  facet_grid(Method~Outcome, scales = "free_x", labeller = labeller(Outcome = c("Empagliflozin plain" = "Empagliflozin Initiation", "Empagliflozin" = "Empareg after eligibility criteria"))) +
  theme(strip.text = element_text(size = 15)) +
  scale_fill_manual(values = c("white", "grey50"), name="", guide="none") +
  ggtitle("Mendelian Randomization Results (all methods)") +
  labs(y = "Tested covariates as Exposure", x = "Beta")+
  theme(axis.text.x = element_text(angle = 0, size = 15),    
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.spacing = unit(2, "lines"),
        legend.position = "",
        legend.text = element_blank(),
        legend.title = element_blank(), # element_text(size = 11),
        plot.margin = margin(1, 1, 1, 1, "cm"), 
        plot.title = element_text(hjust = 0.5, size = 20,),
  )
