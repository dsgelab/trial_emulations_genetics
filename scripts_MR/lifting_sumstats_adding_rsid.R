library(tidyverse)
library(data.table)
library(rtracklayer)
source("/home/ivm/config_te.R")

#######################################################
### for liftover ###
chain_file <- chain_file_liftover
chainObject <- import.chain(chain_file)
#######################################################

### rsID dictionary ###
rsid_dict <- fread(rsid_dict_)

# ALT ---------------------------------------------------------------------

# Source: Pan-UKBB (https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1450719288)

# Read Sumstats
alt_ss <- fread(alt_ss_)
glimpse(alt_ss)
# Keep only EUR ancestry
alt_ss <- alt_ss %>% dplyr::select(chr, pos, ref, alt, af_EUR, beta_EUR, se_EUR, pval_EUR) %>% filter(!is.na(beta_EUR))

alt_ss <-  alt_ss %>% dplyr::rename("chr_position" = "pos",
                             "effect_allele" = "alt",
                             "other_allele" = "ref",
                             "effect_weight" = "beta_EUR", 
                             "se" = "se_EUR",
                             "pval" = "pval_EUR",
                             "af" = "af_EUR")

alt_ss <-  alt_ss %>% mutate(chr_name = paste0("chr",chr))  
head(alt_ss)

# create track of old bp hg19 build
alt_ss$old_hg19_bp = alt_ss$chr_position   # chr_position is BP

alt_ss$startField = alt_ss$chr_position    # the start position for the SNP
alt_ss$endField = alt_ss$chr_position      # the end position is the start because it's a SNP 

# Convert to genomic ranges 
# GRanges
grGWAS_SNPs <- makeGRangesFromDataFrame(
  alt_ss, # df that I want to convert to GRanges
  seqnames.field = "chr_name",  # name of the chromosome information
  start.field = "startField",   # start of the SNP
  end.field = "endField",       # end of the SNP
  keep.extra.columns = T        # keep all other columns after conversion 
)

lifted_alt_ss <- as.tibble(liftOver(grGWAS_SNPs, chainObject))

lifted_alt_ss <- lifted_alt_ss %>% dplyr::select("chr_name" = seqnames, "chr_position" = start, effect_allele, other_allele, 
                                                 effect_weight, se, pval, af, old_hg19_bp) %>% 
  mutate(chromosome = sub("^chr", "", chr_name))

head(lifted_alt_ss)
lifted_alt_ss <- lifted_alt_ss %>% mutate(variant = paste(gsub('chr', '', chr_name), chr_position, other_allele, effect_allele, sep = "_"),
                                          samplesize = as.numeric(400822)) # samplesize from column X from https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1450719288

# Add rsIDs for clumping
lifted_alt_ss <- lifted_alt_ss %>% mutate(chrom_pos = paste0(chromosome, "_", chr_position)) %>% 
  left_join(rsid_dict, by = c("chrom_pos" = "V2")) %>% dplyr::rename("rsid" = "V1") %>% select(-V3)

fwrite(lifted_alt_ss, lifted_alt_ss_)


# create List for clumping 
lifted_alt_snplist <- lifted_alt_ss %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
lifted_alt_snplist <- lifted_alt_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(lifted_alt_snplist, lifted_alt_snplist_, sep = "\t")




# LDL ---------------------------------------------------------------------

# Source: Pan-UKBB (https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1450719288)

# Read Sumstats
ldl_ss <- fread(ldl_ss_)
glimpse(ldl_ss)

ldl_ss <- ldl_ss %>% dplyr::select(chr, pos, ref, alt, af_EUR, beta_EUR, se_EUR, pval_EUR) %>% filter(!is.na(beta_EUR))

ldl_ss <-  ldl_ss %>% dplyr::rename("chr_position" = "pos",
                                    "effect_allele" = "alt",
                                    "other_allele" = "ref",
                                    "effect_weight" = "beta_EUR", 
                                    "se" = "se_EUR",
                                    "pval" = "pval_EUR",
                                    "af" = "af_EUR")

ldl_ss <-  ldl_ss %>% mutate(chr_name = paste0("chr",chr))  
head(ldl_ss)

# create track of old bp hg19 build
ldl_ss$old_hg19_bp = ldl_ss$chr_position   # chr_position is BP

ldl_ss$startField = ldl_ss$chr_position    # the start position for the SNP
ldl_ss$endField = ldl_ss$chr_position      # the end position is the start because it's a SNP 

# Convert to genomic ranges 
# GRanges
grGWAS_SNPs <- makeGRangesFromDataFrame(
  ldl_ss, # df that I want to convert to GRanges
  seqnames.field = "chr_name",  # name of the chromosome information
  start.field = "startField",   # start of the SNP
  end.field = "endField",       # end of the SNP
  keep.extra.columns = T        # keep all other columns after conversion 
)

lifted_ldl_ss <- as.tibble(liftOver(grGWAS_SNPs, chainObject))

lifted_ldl_ss <- lifted_ldl_ss %>% dplyr::select("chr_name" = seqnames, "chr_position" = start, effect_allele, other_allele, 
                                                 effect_weight, se, pval, af, old_hg19_bp) %>% 
  mutate(chromosome = sub("^chr", "", chr_name))

head(lifted_ldl_ss)
lifted_ldl_ss <- lifted_ldl_ss %>% mutate(variant = paste(gsub('chr', '', chr_name), chr_position, other_allele, effect_allele, sep = "_"),
                                          samplesize = as.numeric(400223)) # samplesize from column X from https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1450719288

# Add rsIDs for clumping
lifted_ldl_ss <- lifted_ldl_ss %>% mutate(chrom_pos = paste0(chromosome, "_", chr_position)) %>% 
  left_join(rsid_dict, by = c("chrom_pos" = "V2")) %>% dplyr::rename("rsid" = "V1") %>% select(-V3)
fwrite(lifted_ldl_ss, lifted_ldl_ss_)


# create List for clumping 
lifted_ldl_snplist <- lifted_ldl_ss %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
lifted_ldl_snplist <- lifted_ldl_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(lifted_ldl_snplist, lifted_ldl_snplist_, sep = "\t")




# HDL ---------------------------------------------------------------------

# Source: Pan-UKBB (https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1450719288)

# Read Sumstats
hdl_ss <- fread(hdl_ss_)
glimpse(hdl_ss)

hdl_ss <- hdl_ss %>% dplyr::select(chr, pos, ref, alt, af_EUR, beta_EUR, se_EUR, pval_EUR) %>% filter(!is.na(beta_EUR))

hdl_ss <-  hdl_ss %>% dplyr::rename("chr_position" = "pos",
                                    "effect_allele" = "alt",
                                    "other_allele" = "ref",
                                    "effect_weight" = "beta_EUR", 
                                    "se" = "se_EUR",
                                    "pval" = "pval_EUR",
                                    "af" = "af_EUR")

hdl_ss <-  hdl_ss %>% mutate(chr_name = paste0("chr",chr))  
head(hdl_ss)

# create track of old bp hg19 build
hdl_ss$old_hg19_bp = hdl_ss$chr_position   # chr_position is BP

hdl_ss$startField = hdl_ss$chr_position    # the start position for the SNP
hdl_ss$endField = hdl_ss$chr_position      # the end position is the start because it's a SNP 

# Convert to genomic ranges 
# GRanges
grGWAS_SNPs <- makeGRangesFromDataFrame(
  hdl_ss, # df that I want to convert to GRanges
  seqnames.field = "chr_name",  # name of the chromosome information
  start.field = "startField",   # start of the SNP
  end.field = "endField",       # end of the SNP
  keep.extra.columns = T        # keep all other columns after conversion 
)

lifted_hdl_ss <- as.tibble(liftOver(grGWAS_SNPs, chainObject))

lifted_hdl_ss <- lifted_hdl_ss %>% dplyr::select("chr_name" = seqnames, "chr_position" = start, effect_allele, other_allele, 
                                                 effect_weight, se, pval, af, old_hg19_bp) %>% 
  mutate(chromosome = sub("^chr", "", chr_name))

head(lifted_hdl_ss)
lifted_hdl_ss <- lifted_hdl_ss %>% mutate(variant = paste(gsub('chr', '', chr_name), chr_position, other_allele, effect_allele, sep = "_"),
                                          samplesize = as.numeric(367021)) # samplesize from column X from https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1450719288

# Add rsIDs for clumping
lifted_hdl_ss <- lifted_hdl_ss %>% mutate(chrom_pos = paste0(chromosome, "_", chr_position)) %>% 
  left_join(rsid_dict, by = c("chrom_pos" = "V2")) %>% dplyr::rename("rsid" = "V1") %>% select(-V3)
fwrite(lifted_hdl_ss, lifted_hdl_ss_)


# create List for clumping 
lifted_hdl_snplist <- lifted_hdl_ss %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
lifted_hdl_snplist <- lifted_hdl_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(lifted_hdl_snplist, lifted_hdl_snplist_, sep = "\t")




# Hba1c -------------------------------------------------------------------

# Source: Pan-UKBB (https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1450719288)

# Read Sumstats
hba1c_ss <- fread(hba1c_ss_)
glimpse(hba1c_ss)

hba1c_ss <- hba1c_ss %>% dplyr::select(chr, pos, ref, alt, af_EUR, beta_EUR, se_EUR, pval_EUR) %>% filter(!is.na(beta_EUR))

hba1c_ss <-  hba1c_ss %>% dplyr::rename("chr_position" = "pos",
                                    "effect_allele" = "alt",
                                    "other_allele" = "ref",
                                    "effect_weight" = "beta_EUR", 
                                    "se" = "se_EUR",
                                    "pval" = "pval_EUR",
                                    "af" = "af_EUR")


hba1c_ss <-  hba1c_ss %>% mutate(chr_name = paste0("chr",chr))  
head(hba1c_ss)

# create track of old bp hg19 build
hba1c_ss$old_hg19_bp = hba1c_ss$chr_position   # chr_position is BP

hba1c_ss$startField = hba1c_ss$chr_position    # the start position for the SNP
hba1c_ss$endField = hba1c_ss$chr_position      # the end position is the start because it's a SNP 

# Convert to genomic ranges 
# GRanges
grGWAS_SNPs <- makeGRangesFromDataFrame(
  hba1c_ss, # df that I want to convert to GRanges
  seqnames.field = "chr_name",  # name of the chromosome information
  start.field = "startField",   # start of the SNP
  end.field = "endField",       # end of the SNP
  keep.extra.columns = T        # keep all other columns after conversion 
)

lifted_hba1c_ss <- as.tibble(liftOver(grGWAS_SNPs, chainObject))

lifted_hba1c_ss <- lifted_hba1c_ss %>% dplyr::select("chr_name" = seqnames, "chr_position" = start, effect_allele, other_allele, 
                                                 effect_weight, se, pval, af, old_hg19_bp) %>% 
  mutate(chromosome = sub("^chr", "", chr_name))

head(lifted_hba1c_ss)
lifted_hba1c_ss <- lifted_hba1c_ss %>% mutate(variant = paste(gsub('chr', '', chr_name), chr_position, other_allele, effect_allele, sep = "_"),
                                          samplesize = as.numeric(400825)) # samplesize from column X from https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1450719288


# Add rsIDs for clumping
lifted_hba1c_ss <- lifted_hba1c_ss %>% mutate(chrom_pos = paste0(chromosome, "_", chr_position)) %>% 
  left_join(rsid_dict, by = c("chrom_pos" = "V2")) %>% dplyr::rename("rsid" = "V1") %>% select(-V3)

fwrite(lifted_hba1c_ss, lifted_hba1c_ss_)


# create List for clumping 
lifted_hba1c_snplist <- lifted_hba1c_ss %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
lifted_hba1c_snplist <- lifted_hba1c_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(lifted_hba1c_snplist, , sep = "\t")



# Asthma ------------------------------------------------------------------

# Source: Han et.al. 2020 (https://www.nature.com/articles/s41467-020-15649-3) (https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010042/)

# Read Sumstats
asthma_ss <- fread(asthma_ss_)
glimpse(asthma_ss)
asthma_ss <-  asthma_ss %>% mutate(effect_weight = log(OR),
                                   SE = (log(OR_95U) - log(OR)) / 1.96) 

asthma_ss <-  asthma_ss %>% dplyr::rename("chr" = "CHR",
                                          "chr_position" = "BP",
                                        "effect_allele" = "EA",
                                        "other_allele" = "NEA",
                                        "effect_weight" = "effect_weight", 
                                        "se" = "SE",
                                        "pval" = "P",
                                        "af" = "EAF")  

asthma_ss <-  asthma_ss %>% mutate(chr_name = paste0("chr",chr))  
head(asthma_ss)

# create track of old bp hg19 build
asthma_ss$old_hg19_bp = asthma_ss$chr_position   # chr_position is BP

asthma_ss$startField = asthma_ss$chr_position    # the start position for the SNP
asthma_ss$endField = asthma_ss$chr_position      # the end position is the start because it's a SNP 

# Convert to genomic ranges 
# GRanges
grGWAS_SNPs <- makeGRangesFromDataFrame(
  asthma_ss, # df that I want to convert to GRanges
  seqnames.field = "chr_name",  # name of the chromosome information
  start.field = "startField",   # start of the SNP
  end.field = "endField",       # end of the SNP
  keep.extra.columns = T        # keep all other columns after conversion 
)

lifted_asthma_ss <- as.tibble(liftOver(grGWAS_SNPs, chainObject))

lifted_asthma_ss <- lifted_asthma_ss %>% dplyr::select("chr_name" = seqnames, "chr_position" = start, effect_allele, other_allele, 
                                                     effect_weight, se, pval, af, OR, OR_95L, OR_95U, INFO, old_hg19_bp, "samplesize" = N) %>% 
  mutate(chromosome = sub("^chr", "", chr_name))

head(lifted_asthma_ss)
lifted_asthma_ss <- lifted_asthma_ss %>% mutate(variant = paste(gsub('chr', '', chr_name), chr_position, other_allele, effect_allele, sep = "_"))
              
                                
# Add rsIDs for clumping
lifted_asthma_ss <- lifted_asthma_ss %>% mutate(chrom_pos = paste0(chromosome, "_", chr_position)) %>% 
  left_join(rsid_dict, by = c("chrom_pos" = "V2")) %>% dplyr::rename("rsid" = "V1") %>% select(-V3)

fwrite(lifted_asthma_ss, lifted_asthma_ss_)


# create List for clumping 
lifted_asthma_snplist <- lifted_asthma_ss %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
lifted_asthma_snplist <- lifted_asthma_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(lifted_asthma_snplist, lifted_asthma_snplist_, sep = "\t")




# MDD ---------------------------------------------------------------------

# Source: Han et.al. 2020 (https://www.nature.com/articles/s41588-018-0090-3#MOESM3) (https://figshare.com/articles/dataset/mdd2018/14672085?file=28169502)

# Read Sumstats
mdd_ss <- fread(mdd_ss_)
glimpse(mdd_ss)

mdd_ss <-  mdd_ss %>% mutate(effect_weight = log(OR),
                             #new_SE = SE / OR,
                             #new_SE = (log(OR+SE*1.96) - log(OR)) / 1.96,
                             N = Nca + Nco,
                             EAF = (FRQ_A_59851*Nca + FRQ_U_113154*Nco)/N ) 

mdd_ss <-  mdd_ss %>% dplyr::rename("chr" = "CHR",
                                          "chr_position" = "BP",
                                          "effect_allele" = "A1",
                                          "other_allele" = "A2",
                                          "effect_weight" = "effect_weight", 
                                          "se" = "SE",
                                          "pval" = "P",
                                          "af" = "EAF")

# samplesize: 480359 (cas: 135458 con: 344901) https://figshare.com/articles/dataset/mdd2018/14672085?file=28169502

mdd_ss <-  mdd_ss %>% mutate(chr_name = paste0("chr",chr))  
head(mdd_ss)

# create track of old bp hg19 build
mdd_ss$old_hg19_bp = mdd_ss$chr_position   # chr_position is BP

mdd_ss$startField = mdd_ss$chr_position    # the start position for the SNP
mdd_ss$endField = mdd_ss$chr_position      # the end position is the start because it's a SNP 

# Convert to genomic ranges 
# GRanges
grGWAS_SNPs <- makeGRangesFromDataFrame(
  mdd_ss, # df that I want to convert to GRanges
  seqnames.field = "chr_name",  # name of the chromosome information
  start.field = "startField",   # start of the SNP
  end.field = "endField",       # end of the SNP
  keep.extra.columns = T        # keep all other columns after conversion 
)

lifted_mdd_ss <- as.tibble(liftOver(grGWAS_SNPs, chainObject))

lifted_mdd_ss <- lifted_mdd_ss %>% dplyr::select("chr_name" = seqnames, "chr_position" = start, effect_allele, other_allele, 
                                                       effect_weight, se, pval, af, OR, INFO, old_hg19_bp, "samplesize" = N, FRQ_A_59851, FRQ_U_113154 ) %>% 
  mutate(chromosome = sub("^chr", "", chr_name))

head(lifted_mdd_ss)
lifted_mdd_ss <- lifted_mdd_ss %>% mutate(variant = paste(gsub('chr', '', chr_name), chr_position, other_allele, effect_allele, sep = "_"))


# Add rsIDs for clumping
lifted_mdd_ss <- lifted_mdd_ss %>% mutate(chrom_pos = paste0(chromosome, "_", chr_position)) %>% 
  left_join(rsid_dict, by = c("chrom_pos" = "V2")) %>% dplyr::rename("rsid" = "V1") %>% select(-V3)

fwrite(lifted_mdd_ss, lifted_mdd_ss_)


# create List for clumping 
lifted_mdd_snplist <- lifted_mdd_ss %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
lifted_mdd_snplist <- lifted_mdd_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(lifted_mdd_snplist, lifted_mdd_snplist_, sep = "\t")




# EDU ---------------------------------------------------------------------

# Source: Lee et al 2018 (https://www.nature.com/articles/s41588-018-0147-3) (https://www.thessgac.com/papers/3)

# Read Sumstats
edu_ss <- fread(edu_ss_)
glimpse(edu_ss)

edu_ss <-  edu_ss %>% dplyr::rename("chr" = "CHR",
                                    "chr_position" = "POS",
                                    "effect_allele" = "A1",
                                    "other_allele" = "A2",
                                    "effect_weight" = "Beta", 
                                    "se" = "SE",
                                    "pval" = "Pval",
                                    "af" = "EAF",
                                    "snp" = "MarkerName")

edu_ss <-  edu_ss %>% mutate(chr_name = paste0("chr",chr))  
head(edu_ss)

# create track of old bp hg19 build
edu_ss$old_hg19_bp = edu_ss$chr_position   # chr_position is BP

edu_ss$startField = edu_ss$chr_position    # the start position for the SNP
edu_ss$endField = edu_ss$chr_position      # the end position is the start because it's a SNP 

# Convert to genomic ranges 
# GRanges
grGWAS_SNPs <- makeGRangesFromDataFrame(
  edu_ss, # df that I want to convert to GRanges
  seqnames.field = "chr_name",  # name of the chromosome information
  start.field = "startField",   # start of the SNP
  end.field = "endField",       # end of the SNP
  keep.extra.columns = T        # keep all other columns after conversion 
)

lifted_edu_ss <- as.tibble(liftOver(grGWAS_SNPs, chainObject))

lifted_edu_ss <- lifted_edu_ss %>% dplyr::select("chr_name" = seqnames, "chr_position" = start, effect_allele, other_allele, 
                                                       effect_weight, se, pval, af, old_hg19_bp, snp) %>% 
  mutate(chromosome = sub("^chr", "", chr_name))

head(lifted_edu_ss)
lifted_edu_ss <- lifted_edu_ss %>% mutate(variant = paste(gsub('chr', '', chr_name), chr_position, other_allele, effect_allele, sep = "_"),
                                          samplesize = as.numeric(1131881)) # samplesize from https://www.nature.com/articles/s41588-018-0147-3

# Add rsIDs for clumping
lifted_edu_ss <- lifted_edu_ss %>% mutate(chrom_pos = paste0(chromosome, "_", chr_position)) %>% 
  left_join(rsid_dict, by = c("chrom_pos" = "V2")) %>% dplyr::rename("rsid" = "V1") %>% select(-V3)

fwrite(lifted_edu_ss, lifted_edu_ss_)


# create List for clumping 
lifted_edu_snplist <- lifted_edu_ss %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
lifted_edu_snplist <- lifted_edu_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(lifted_edu_snplist, lifted_edu_snplist_, sep = "\t")






# BMI ---------------------------------------------------------------------

# Source: GIANT 2018 (https://academic.oup.com/hmg/article/27/20/3641/5067845#134313169) 

# Read Sumstats
bmi_ss <- fread(bmi_ss_)

glimpse(bmi_ss)
bmi_ss <-  bmi_ss %>% filter(!is.na(CHR))
bmi_ss <-  bmi_ss %>% dplyr::rename("chr" = "CHR",
                                    "chr_position" = "POS",
                                    "effect_allele" = "Tested_Allele",
                                    "other_allele" = "Other_Allele",
                                    "effect_weight" = "BETA", 
                                    "se" = "SE",
                                    "pval" = "P",
                                    "af" = "Freq_Tested_Allele",
                                    "snp" = "SNP")

bmi_ss <-  bmi_ss %>% mutate(chr_name = paste0("chr",chr))  
head(bmi_ss)

# create track of old bp hg19 build
bmi_ss$old_hg19_bp = bmi_ss$chr_position   # chr_position is BP

bmi_ss$startField = bmi_ss$chr_position    # the start position for the SNP
bmi_ss$endField = bmi_ss$chr_position      # the end position is the start because it's a SNP 

# Convert to genomic ranges 
# GRanges
grGWAS_SNPs <- makeGRangesFromDataFrame(
  bmi_ss, # df that I want to convert to GRanges
  seqnames.field = "chr_name",  # name of the chromosome information
  start.field = "startField",   # start of the SNP
  end.field = "endField",       # end of the SNP
  keep.extra.columns = T        # keep all other columns after conversion 
)

lifted_bmi_ss <- as.tibble(liftOver(grGWAS_SNPs, chainObject))

lifted_bmi_ss <- lifted_bmi_ss %>% dplyr::select("chr_name" = seqnames, "chr_position" = start, effect_allele, other_allele, 
                                                 effect_weight, se, pval, af, old_hg19_bp, snp, samplesize = N, INFO) %>% 
  mutate(chromosome = sub("^chr", "", chr_name))

head(lifted_bmi_ss)
lifted_bmi_ss <- lifted_bmi_ss %>% mutate(variant = paste(gsub('chr', '', chr_name), chr_position, other_allele, effect_allele, sep = "_"))


# Add rsIDs for clumping
lifted_bmi_ss <- lifted_bmi_ss %>% mutate(chrom_pos = paste0(chromosome, "_", chr_position)) %>% 
  left_join(rsid_dict, by = c("chrom_pos" = "V2")) %>% dplyr::rename("rsid" = "V1") %>% select(-V3)

fwrite(lifted_bmi_ss, lifted_bmi_ss_)


# create List for clumping 
lifted_bmi_snplist <- lifted_bmi_ss %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
lifted_bmi_snplist <- lifted_bmi_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(lifted_bmi_snplist, lifted_bmi_snplist_, sep = "\t")



# CAD ---------------------------------------------------------------------

# Source: Cardiogram 2015 (https://www.nature.com/articles/ng.3396#change-history) (https://www.thessgac.com/papers/3)

# Read Sumstats
cad_ss <- fread(cad_ss_)
glimpse(cad_ss)

cad_ss <-  cad_ss %>% dplyr::rename("chr" = "chr",
                                    "chr_position" = "bp_hg19",
                                    "effect_allele" = "effect_allele",
                                    "other_allele" = "noneffect_allele",
                                    "effect_weight" = "beta", 
                                    "se" = "se_dgc",
                                    "pval" = "p_dgc",
                                    "af" = "effect_allele_freq",
                                    "snp" = "markername")

cad_ss <-  cad_ss %>% mutate(chr_name = paste0("chr",chr))  
head(cad_ss)

# create track of old bp hg19 build
cad_ss$old_hg19_bp = cad_ss$chr_position   # chr_position is BP

cad_ss$startField = cad_ss$chr_position    # the start position for the SNP
cad_ss$endField = cad_ss$chr_position      # the end position is the start because it's a SNP 

# Convert to genomic ranges 
# GRanges
grGWAS_SNPs <- makeGRangesFromDataFrame(
  cad_ss, # df that I want to convert to GRanges
  seqnames.field = "chr_name",  # name of the chromosome information
  start.field = "startField",   # start of the SNP
  end.field = "endField",       # end of the SNP
  keep.extra.columns = T        # keep all other columns after conversion 
)

lifted_cad_ss <- as.tibble(liftOver(grGWAS_SNPs, chainObject))

lifted_cad_ss <- lifted_cad_ss %>% dplyr::select("chr_name" = seqnames, "chr_position" = start, effect_allele, other_allele, 
                                                 effect_weight, se, pval, af, old_hg19_bp, snp, het_pvalue, n_studies, median_info) %>% 
  mutate(chromosome = sub("^chr", "", chr_name))

head(lifted_cad_ss)
lifted_cad_ss <- lifted_cad_ss %>% mutate(variant = paste(gsub('chr', '', chr_name), chr_position, other_allele, effect_allele, sep = "_"),
                                          samplesize = as.numeric(184305)) # samplesize from https://www.nature.com/articles/ng.3396#change-history


# Add rsIDs for clumping
lifted_cad_ss <- lifted_cad_ss %>% mutate(chrom_pos = paste0(chromosome, "_", chr_position)) %>% 
  left_join(rsid_dict, by = c("chrom_pos" = "V2")) %>% dplyr::rename("rsid" = "V1") %>% select(-V3)

fwrite(lifted_cad_ss, lifted_cad_ss_)


# create List for clumping 
lifted_cad_snplist <- lifted_cad_ss %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
lifted_cad_snplist <- lifted_cad_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(lifted_cad_snplist, lifted_cad_snplist_, sep = "\t")




# CKD ---------------------------------------------------------------------

# Source: CKDGen (https://ckdgen.imbi.uni-freiburg.de/)

# Read Sumstats
ckd_ss <- fread(ckd_ss_)
glimpse(ckd_ss)

ckd_ss <-  ckd_ss %>% dplyr::rename("chr" = "Chr",
                                    "chr_position" = "Pos_b37",
                                    "effect_allele" = "Allele1",
                                    "other_allele" = "Allele2",
                                    "effect_weight" = "Effect", 
                                    "se" = "StdErr",
                                    "pval" = `P-value`,
                                    "af" = "Freq1",
                                    "snp" = "RSID")

ckd_ss <-  ckd_ss %>% mutate(chr_name = paste0("chr",chr))  
head(ckd_ss)

# create track of old bp hg19 build
ckd_ss$old_hg19_bp = ckd_ss$chr_position   # chr_position is BP

ckd_ss$startField = ckd_ss$chr_position    # the start position for the SNP
ckd_ss$endField = ckd_ss$chr_position      # the end position is the start because it's a SNP 

# Convert to genomic ranges 
# GRanges
grGWAS_SNPs <- makeGRangesFromDataFrame(
  ckd_ss, # df that I want to convert to GRanges
  seqnames.field = "chr_name",  # name of the chromosome information
  start.field = "startField",   # start of the SNP
  end.field = "endField",       # end of the SNP
  keep.extra.columns = T        # keep all other columns after conversion 
)

lifted_ckd_ss <- as.tibble(liftOver(grGWAS_SNPs, chainObject))

lifted_ckd_ss <- lifted_ckd_ss %>% dplyr::select("chr_name" = seqnames, "chr_position" = start, effect_allele, other_allele, 
                                                 effect_weight, se, pval, af, old_hg19_bp, snp, samplesize = n_total_sum) %>% 
  mutate(chromosome = sub("^chr", "", chr_name))

head(lifted_ckd_ss)
lifted_ckd_ss <- lifted_ckd_ss %>% mutate(variant = paste(gsub('chr', '', chr_name), chr_position, other_allele, effect_allele, sep = "_"))


# Add rsIDs for clumping
lifted_ckd_ss <- lifted_ckd_ss %>% mutate(chrom_pos = paste0(chromosome, "_", chr_position)) %>% 
  left_join(rsid_dict, by = c("chrom_pos" = "V2")) %>% dplyr::rename("rsid" = "V1") %>% select(-V3)

fwrite(lifted_ckd_ss, lifted_ckd_ss_)


# create List for clumping 
lifted_ckd_snplist <- lifted_ckd_ss %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
lifted_ckd_snplist <- lifted_ckd_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(lifted_ckd_snplist, lifted_ckd_snplist_, sep = "\t")





# T2D ---------------------------------------------------------------------

# these sumstats are already in hg38 and just need to be clumped
t2d_ss <- fread(t2d_ss_)
head(t2d_ss)
t2d_ss <- t2d_ss %>% dplyr::select("chr_name" = seqnames, "chr_position" = pos, effect_allele, other_allele, 
                                                 "effect_weight" = BETA, se, "pval" = Pval, old_hg19_bp, snp, "chrom_pos" = V2, "rsid" = rsID) %>% 
  mutate(samplesize = 2535601 )
fwrite(t2d_ss, write_t2d_ss_)


# create List for clumping 
t2d_ss_snplist <- t2d_ss %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
t2d_ss_snplist <- t2d_ss_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(t2d_ss_snplist, t2d_ss_snplist_, sep = "\t")


# HF  ---------------------------------------------------------------------

# Source: HERMES HF GWAS (https://cvd.hugeamp.org/dinspector.html?dataset=GWAS_HERMES_eu)

# Read Sumstats
hf_ss <- fread(hf_ss_)
glimpse(hf_ss)

hf_ss <-  hf_ss %>% dplyr::rename("chr" = "CHR",
                                    "chr_position" = "BP",
                                    "effect_allele" = "A1",
                                    "other_allele" = "A2",
                                    "effect_weight" = "b", 
                                    "se" = "se",
                                    "pval" = "p",
                                    "af" = "freq",
                                    "snp" = "SNP") %>% 
  mutate(effect_allele = toupper(effect_allele),
         other_allele = toupper(other_allele))

hf_ss <-  hf_ss %>% mutate(chr_name = paste0("chr",chr))  
head(hf_ss)

# create track of old bp hg19 build
hf_ss$old_hg19_bp = hf_ss$chr_position   # chr_position is BP

hf_ss$startField = hf_ss$chr_position    # the start position for the SNP
hf_ss$endField = hf_ss$chr_position      # the end position is the start because it's a SNP 

# Convert to genomic ranges 
# GRanges
grGWAS_SNPs <- makeGRangesFromDataFrame(
  hf_ss, # df that I want to convert to GRanges
  seqnames.field = "chr_name",  # name of the chromosome information
  start.field = "startField",   # start of the SNP
  end.field = "endField",       # end of the SNP
  keep.extra.columns = T        # keep all other columns after conversion 
)

lifted_hf_ss <- as.tibble(liftOver(grGWAS_SNPs, chainObject))

lifted_hf_ss <- lifted_hf_ss %>% dplyr::select("chr_name" = seqnames, "chr_position" = start, effect_allele, other_allele, 
                                                 effect_weight, se, pval, af, old_hg19_bp, snp, samplesize = N) %>% 
  mutate(chromosome = sub("^chr", "", chr_name))

head(lifted_hf_ss)
lifted_hf_ss <- lifted_hf_ss %>% mutate(variant = paste(gsub('chr', '', chr_name), chr_position, other_allele, effect_allele, sep = "_"))


# rsIDs for clumping is already present, just rename column for unity
all(!is.na(lifted_hf_ss$snp))
lifted_hf_ss <- lifted_hf_ss %>% dplyr::rename("rsid" = "snp") 
fwrite(lifted_hf_ss, lifted_hf_ss_)


# create List for clumping 
lifted_hf_snplist <- lifted_hf_ss %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
lifted_hf_snplist <- lifted_hf_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(lifted_hf_snplist, lifted_hf_snplist_, sep = "\t")




# TG  ---------------------------------------------------------------------

# Source: Surakka et al (https://www.nature.com/articles/ng.3300)

# Read Sumstats
tg_ss <- fread(tg_ss_)
glimpse(tg_ss)
tg_ss <- tg_ss %>% mutate(chr_name = sub(":.*", "", rs_number),
                          chr_position = as.numeric(sub(".*:", "", rs_number)),
                          chr = as.numeric(sub("chr", "", chr_name)),
                          # calculate the se from the pvalue and beta
                          se = abs(beta / (qnorm(1-`p-value`/2))))

tg_ss <-  tg_ss %>% dplyr::rename("effect_allele" = "reference_allele",  # CAVE: not sure if this is the right way around!
                                    "other_allele" = "other_allele",
                                    "effect_weight" = "beta", 
                                    "se" = "se",
                                    "pval" = `p-value`,
                                    "af" = "eaf")

# create track of old bp hg19 build
tg_ss$old_hg19_bp = tg_ss$chr_position   # chr_position is BP

tg_ss$startField = tg_ss$chr_position    # the start position for the SNP
tg_ss$endField = tg_ss$chr_position      # the end position is the start because it's a SNP 

# Convert to genomic ranges 
# GRanges
grGWAS_SNPs <- makeGRangesFromDataFrame(
  tg_ss, # df that I want to convert to GRanges
  seqnames.field = "chr_name",  # name of the chromosome information
  start.field = "startField",   # start of the SNP
  end.field = "endField",       # end of the SNP
  keep.extra.columns = T        # keep all other columns after conversion 
)

lifted_tg_ss <- as.tibble(liftOver(grGWAS_SNPs, chainObject))

lifted_tg_ss <- lifted_tg_ss %>% dplyr::select("chr_name" = seqnames, "chr_position" = start, effect_allele, other_allele, 
                                                 effect_weight, se, pval, af, old_hg19_bp) %>% 
  mutate(chromosome = sub("^chr", "", chr_name))

head(lifted_tg_ss)
lifted_tg_ss <- lifted_tg_ss %>% mutate(variant = paste(gsub('chr', '', chr_name), chr_position, other_allele, effect_allele, sep = "_"),
                                        samplesize = as.numeric(361194)) # samplesize from https://www.nature.com/articles/s41467-019-13690-5


# Add rsIDs for clumping
lifted_tg_ss <- lifted_tg_ss %>% mutate(chrom_pos = paste0(chromosome, "_", chr_position)) %>% 
  left_join(rsid_dict, by = c("chrom_pos" = "V2")) %>% dplyr::rename("rsid" = "V1") %>% select(-V3)

fwrite(lifted_tg_ss, lifted_tg_ss_)


# create List for clumping 
lifted_tg_snplist <- lifted_tg_ss %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
lifted_tg_snplist <- lifted_tg_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(lifted_tg_snplist, lifted_tg_snplist_, sep = "\t")




# SAH ---------------------------------------------------------------------

# Source: sahGen (https://sahgen.imbi.uni-freiburg.de/)

# Read Sumstats
sah_ss <- fread(sah_ss_)
glimpse(sah_ss)

sah_ss <-  sah_ss %>% dplyr::rename("chr" = "CHR",
                                    "chr_position" = "BP",
                                    "effect_allele" = "A_EFF",
                                    "other_allele" = "A_NONEFF",
                                    "effect_weight" = "BETA", 
                                    "se" = "SE",
                                    "pval" = "P",
                                    "af" = "Freq_EFF",
                                    "snp" = "SNP")

sah_ss <-  sah_ss %>% mutate(chr_name = paste0("chr",chr))  
head(sah_ss)

# create track of old bp hg19 build
sah_ss$old_hg19_bp = sah_ss$chr_position   # chr_position is BP

sah_ss$startField = sah_ss$chr_position    # the start position for the SNP
sah_ss$endField = sah_ss$chr_position      # the end position is the start because it's a SNP 

# Convert to genomic ranges 
# GRanges
grGWAS_SNPs <- makeGRangesFromDataFrame(
  sah_ss, # df that I want to convert to GRanges
  seqnames.field = "chr_name",  # name of the chromosome information
  start.field = "startField",   # start of the SNP
  end.field = "endField",       # end of the SNP
  keep.extra.columns = T        # keep all other columns after conversion 
)

lifted_sah_ss <- as.tibble(liftOver(grGWAS_SNPs, chainObject))

lifted_sah_ss <- lifted_sah_ss %>% dplyr::select("chr_name" = seqnames, "chr_position" = start, effect_allele, other_allele, 
                                                 effect_weight, se, pval, af, old_hg19_bp, snp, samplesize = Neff) %>% 
  mutate(chromosome = sub("^chr", "", chr_name))

head(lifted_sah_ss)
lifted_sah_ss <- lifted_sah_ss %>% mutate(variant = paste(gsub('chr', '', chr_name), chr_position, other_allele, effect_allele, sep = "_"))


# rename snp into rsid
lifted_sah_ss <- lifted_sah_ss %>% dplyr::rename("rsid" = "snp") 

fwrite(lifted_sah_ss, lifted_sah_ss_)


# create List for clumping 
lifted_sah_snplist <- lifted_sah_ss %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
lifted_sah_snplist <- lifted_sah_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(lifted_sah_snplist, lifted_sah_snplist_, sep = "\t")




# SBP ---------------------------------------------------------------------

# Source: UKBB ICBP (https://www.kp4cd.org/node/679)

# Read Sumstats
sbp_ss <- fread(sbp_ss_)
glimpse(sbp_ss)
sbp_ss <- sbp_ss %>% mutate(chr = sapply(str_split(MarkerName, pattern = ":"), `[`, 1),
                            chr_position = sapply(str_split(MarkerName, pattern = ":"), `[`, 2))

sbp_ss <-  sbp_ss %>% dplyr::rename("chr" = "chr",
                                    "chr_position" = "chr_position",
                                    "effect_allele" = "Allele1",
                                    "other_allele" = "Allele2",
                                    "effect_weight" = "Effect", 
                                    "se" = "StdErr",
                                    "pval" = "P",
                                    "af" = "Freq1") %>% 
  mutate(effect_allele = toupper(effect_allele),
         other_allele = toupper(other_allele)) 

sbp_ss <-  sbp_ss %>% mutate(chr_name = paste0("chr",chr)) 
head(sbp_ss)

# create track of old bp hg19 build
sbp_ss$old_hg19_bp = sbp_ss$chr_position   # chr_position is BP

sbp_ss$startField = sbp_ss$chr_position    # the start position for the SNP
sbp_ss$endField = sbp_ss$chr_position      # the end position is the start because it's a SNP 

# Convert to genomic ranges 
# GRanges
grGWAS_SNPs <- makeGRangesFromDataFrame(
  sbp_ss, # df that I want to convert to GRanges
  seqnames.field = "chr_name",  # name of the chromosome information
  start.field = "startField",   # start of the SNP
  end.field = "endField",       # end of the SNP
  keep.extra.columns = T        # keep all other columns after conversion 
)

lifted_sbp_ss <- as.tibble(liftOver(grGWAS_SNPs, chainObject))

lifted_sbp_ss <- lifted_sbp_ss %>% dplyr::select("chr_name" = seqnames, "chr_position" = start, effect_allele, other_allele, 
                                                 effect_weight, se, pval, af, old_hg19_bp, samplesize = N_effective, TotalSampleSize) %>% 
  mutate(chromosome = sub("^chr", "", chr_name))

head(lifted_sbp_ss)
lifted_sbp_ss <- lifted_sbp_ss %>% mutate(variant = paste(gsub('chr', '', chr_name), chr_position, other_allele, effect_allele, sep = "_"))


# Add rsIDs for clumping
lifted_sbp_ss <- lifted_sbp_ss %>% mutate(chrom_pos = paste0(chromosome, "_", chr_position)) %>% 
  left_join(rsid_dict, by = c("chrom_pos" = "V2")) %>% dplyr::rename("rsid" = "V1") %>% select(-V3)
# cave!: multiple many-to-many reationships -> remove duplicate rsIDs later 
# sum(base::duplicated((lifted_sbp_ss %>% filter(rsid != ""))$rsid))
fwrite(lifted_sbp_ss, lifted_sbp_ss_)


# create List for clumping 
lifted_sbp_snplist <- lifted_sbp_ss %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
lifted_sbp_snplist <- lifted_sbp_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(lifted_sbp_snplist, lifted_sbp_snplist_, sep = "\t")




# LIV ---------------------------------------------------------------------

# Source: https://journals.lww.com/hepcomm/fulltext/2022/02000/genome_wide_association_study_of_nafld_using.6.aspx

# Read Sumstats
liv_ss <- fread(liv_ss_)
glimpse(liv_ss)
liv_ss <-  liv_ss %>% mutate(chr_name = chr,
                             chr = sub("^chr", "", chr_name)) %>% 
  dplyr::rename("chromosome" = "chr",
                                    "chr_position" = "SNP",
                                    "effect_allele" = "effect_allele",
                                    "other_allele" = "other_allele",
                                    "effect_weight" = "beta", 
                                    "se" = "standard_error",
                                    "pval" = "p_value",
                                    "af" = "effect_allele_frequency",
                                    "rsid" = "variant_id") %>% 
  mutate(variant = paste(gsub('chr', '', chr_name), chr_position, other_allele, effect_allele, sep = "_"),
         samplesize = 300000) 

# rsIDs for clumping is already present
all(!is.na(liv_ss$rsid))
fwrite(liv_ss, write_liv_ss_)

# create List for clumping 
liv_ss_snplist <- liv_ss %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
liv_ss_snplist <- liv_ss_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(liv_ss_snplist, liv_ss_snplist_, sep = "\t")





# Stroke ------------------------------------------------------------------

# Source: Megastroke (https://www.nature.com/articles/s41588-018-0058-3)

# Read Sumstats
stroke_ss <- fread(stroke_ss_)
glimpse(stroke_ss)

stroke_ss <-  stroke_ss %>% dplyr::rename("effect_allele" = "Allele1",
                                    "other_allele" = "Allele2",
                                    "effect_weight" = "Effect", 
                                    "se" = "StdErr",
                                    "pval" = `P-value`,
                                    "af" = "Freq1",
                                    "rsid" = "MarkerName") %>% 
  mutate(effect_allele = toupper(effect_allele),
         other_allele = toupper(other_allele),
         samplesize = 446696)  # https://www.megastroke.org/index.html


stroke_ss <-  stroke_ss %>% mutate(chr_name = paste0("chr",chr))  
head(stroke_ss)

# Add chrom_pos, because only rsID is here 
lifted_stroke_ss <- stroke_ss %>% 
  left_join(rsid_dict, by = c("rsid" = "V1")) %>% dplyr::rename("chrom_pos" = "V2") %>% select(-V3)
lifted_stroke_ss <- lifted_stroke_ss %>% mutate(chromosome = sapply(str_split(chrom_pos, pattern = "_"), `[`, 1),
                                                chr_position = sapply(str_split(chrom_pos, pattern = "_"), `[`, 2))

fwrite(lifted_stroke_ss, lifted_stroke_ss_)


# create List for clumping 
lifted_stroke_snplist <- lifted_stroke_ss %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
lifted_stroke_snplist <- lifted_stroke_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(lifted_stroke_snplist, lifted_stroke_snplist_, sep = "\t")





# CRP ---------------------------------------------------------------------

# Source: Pan-UKBB (https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1450719288)

# Read Sumstats
crp_ss <- fread(crp_ss_)
glimpse(crp_ss)
crp_ss <- crp_ss %>% dplyr::select(chr, pos, ref, alt, af_EUR, beta_EUR, se_EUR, pval_EUR) %>% filter(!is.na(beta_EUR))

crp_ss <-  crp_ss %>% dplyr::rename("chr_position" = "pos",
                                    "effect_allele" = "alt",
                                    "other_allele" = "ref",
                                    "effect_weight" = "beta_EUR", 
                                    "se" = "se_EUR",
                                    "pval" = "pval_EUR",
                                    "af" = "af_EUR")

crp_ss <-  crp_ss %>% mutate(chr_name = paste0("chr",chr))  
head(crp_ss)

# create track of old bp hg19 build
crp_ss$old_hg19_bp = crp_ss$chr_position   # chr_position is BP

crp_ss$startField = crp_ss$chr_position    # the start position for the SNP
crp_ss$endField = crp_ss$chr_position      # the end position is the start because it's a SNP 

# Convert to genomic ranges 
# GRanges
grGWAS_SNPs <- makeGRangesFromDataFrame(
  crp_ss, # df that I want to convert to GRanges
  seqnames.field = "chr_name",  # name of the chromosome information
  start.field = "startField",   # start of the SNP
  end.field = "endField",       # end of the SNP
  keep.extra.columns = T        # keep all other columns after conversion 
)

lifted_crp_ss <- as.tibble(liftOver(grGWAS_SNPs, chainObject))

lifted_crp_ss <- lifted_crp_ss %>% dplyr::select("chr_name" = seqnames, "chr_position" = start, effect_allele, other_allele, 
                                                 effect_weight, se, pval, af, old_hg19_bp) %>% 
  mutate(chromosome = sub("^chr", "", chr_name))

head(lifted_crp_ss)
lifted_crp_ss <- lifted_crp_ss %>% mutate(variant = paste(gsub('chr', '', chr_name), chr_position, other_allele, effect_allele, sep = "_"),
                                          samplesize = as.numeric(419693)) # samplesize from column X from https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1450719288

# Add rsIDs for clumping
lifted_crp_ss <- lifted_crp_ss %>% mutate(chrom_pos = paste0(chromosome, "_", chr_position)) %>% 
  left_join(rsid_dict, by = c("chrom_pos" = "V2")) %>% dplyr::rename("rsid" = "V1") %>% select(-V3)

fwrite(lifted_crp_ss, lifted_crp_ss_)


# create List for clumping 
lifted_crp_snplist <- lifted_crp_ss %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
lifted_crp_snplist <- lifted_crp_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(lifted_crp_snplist, lifted_crp_snplist_, sep = "\t")





# AST ---------------------------------------------------------------------

# Source: Pan-UKBB (https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1450719288)

# Read Sumstats
ast_ss <- fread(ast_ss_)
glimpse(ast_ss)
ast_ss <- ast_ss %>% dplyr::select(chr, pos, ref, alt, af_EUR, beta_EUR, se_EUR, pval_EUR) %>% filter(!is.na(beta_EUR))

ast_ss <-  ast_ss %>% dplyr::rename("chr_position" = "pos",
                                    "effect_allele" = "alt",
                                    "other_allele" = "ref",
                                    "effect_weight" = "beta_EUR", 
                                    "se" = "se_EUR",
                                    "pval" = "pval_EUR",
                                    "af" = "af_EUR")

ast_ss <-  ast_ss %>% mutate(chr_name = paste0("chr",chr))  
head(ast_ss)

# create track of old bp hg19 build
ast_ss$old_hg19_bp = ast_ss$chr_position   # chr_position is BP

ast_ss$startField = ast_ss$chr_position    # the start position for the SNP
ast_ss$endField = ast_ss$chr_position      # the end position is the start because it's a SNP 

# Convert to genomic ranges 
# GRanges
grGWAS_SNPs <- makeGRangesFromDataFrame(
  ast_ss, # df that I want to convert to GRanges
  seqnames.field = "chr_name",  # name of the chromosome information
  start.field = "startField",   # start of the SNP
  end.field = "endField",       # end of the SNP
  keep.extra.columns = T        # keep all other columns after conversion 
)

lifted_ast_ss <- as.tibble(liftOver(grGWAS_SNPs, chainObject))

lifted_ast_ss <- lifted_ast_ss %>% dplyr::select("chr_name" = seqnames, "chr_position" = start, effect_allele, other_allele, 
                                                 effect_weight, se, pval, af, old_hg19_bp) %>% 
  mutate(chromosome = sub("^chr", "", chr_name))

head(lifted_ast_ss)
lifted_ast_ss <- lifted_ast_ss %>% mutate(variant = paste(gsub('chr', '', chr_name), chr_position, other_allele, effect_allele, sep = "_"),
                                          samplesize = as.numeric(399482)) # samplesize from column X from https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1450719288

# Add rsIDs for clumping
lifted_ast_ss <- lifted_ast_ss %>% mutate(chrom_pos = paste0(chromosome, "_", chr_position)) %>% 
  left_join(rsid_dict, by = c("chrom_pos" = "V2")) %>% dplyr::rename("rsid" = "V1") %>% select(-V3)

fwrite(lifted_ast_ss, lifted_ast_ss_)


# create List for clumping 
lifted_ast_snplist <- lifted_ast_ss %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
lifted_ast_snplist <- lifted_ast_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(lifted_ast_snplist, lifted_ast_snplist_, sep = "\t")





# AF ----------------------------------------------------------------------

# Source: AF (https://www.nature.com/articles/s41588-018-0133-9)

# Read Sumstats
af_ss <- fread(af_ss_)
glimpse(af_ss)

af_ss <-  af_ss %>% dplyr::rename("chr" = "chr",
                                    "chr_position" = "pos",
                                    "effect_allele" = "Allele1",
                                    "other_allele" = "Allele2",
                                    "effect_weight" = "Effect", 
                                    "se" = "StdErr",
                                    "pval" = `P-value`,
                                    "snp" = "MarkerName")


af_ss <-  af_ss %>% mutate(chr_name = paste0("chr",chr),
                           aef = 0.672) 
head(af_ss)

# create track of old bp hg19 build
af_ss$old_hg19_bp = af_ss$chr_position   # chr_position is BP

af_ss$startField = af_ss$chr_position    # the start position for the SNP
af_ss$endField = af_ss$chr_position      # the end position is the start because it's a SNP 

# Convert to genomic ranges 
# GRanges
grGWAS_SNPs <- makeGRangesFromDataFrame(
  af_ss, # df that I want to convert to GRanges
  seqnames.field = "chr_name",  # name of the chromosome information
  start.field = "startField",   # start of the SNP
  end.field = "endField",       # end of the SNP
  keep.extra.columns = T        # keep all other columns after conversion 
)

lifted_af_ss <- as.tibble(liftOver(grGWAS_SNPs, chainObject))

lifted_af_ss <- lifted_af_ss %>% dplyr::select("chr_name" = seqnames, "chr_position" = start, effect_allele, other_allele, 
                                                 effect_weight, se, pval, "af" = aef, old_hg19_bp, snp) %>% 
  mutate(chromosome = sub("^chr", "", chr_name))

head(lifted_af_ss)
lifted_af_ss <- lifted_af_ss %>% mutate(variant = paste(gsub('chr', '', chr_name), chr_position, other_allele, effect_allele, sep = "_"),
                                        samplesize = as.numeric(351017)) # samplesize from the paper above (fig1)

# Add rsIDs for clumping
lifted_af_ss <- lifted_af_ss %>% dplyr::rename("rsid" = "snp") 

fwrite(lifted_af_ss, lifted_af_ss_)


# create List for clumping 
lifted_af_snplist <- lifted_af_ss %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
lifted_af_snplist <- lifted_af_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(lifted_af_snplist, lifted_af_snplist_, sep = "\t")

