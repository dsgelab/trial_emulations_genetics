# This script has the purpose to format trial GWASs and create files for clumping

### rsID dictionary ###
rsid_dict <- fread(rsid_dict_)



# EMPA GWAS (initiators vs non-initiators) --------------------------------

empa_gwas <- fread(empa_plain_gwas_raw_)
empa_gwas <- empa_gwas %>% dplyr::select("chromosome" = `#chrom`, "chr_position" = pos, effect_allele = alt, other_allele = ref, 
                                         "effect_weight" = beta, "se" = sebeta, pval, "af" = af_alt)

empa_gwas <- empa_gwas %>%  mutate(variant = paste(chromosome, chr_position, other_allele, effect_allele, sep = "_"),
                                   samplesize = as.numeric(426776)) 

# Add rsIDs for clumping
empa_gwas <- empa_gwas %>% mutate(chrom_pos = paste0(chromosome, "_", chr_position)) %>% 
  left_join(rsid_dict, by = c("chrom_pos" = "V2")) %>% dplyr::rename("rsid" = "V1") %>% select(-V3)

fwrite(empa_gwas, empa_plain_gwas_)

# create List for clumping 
empa_snplist <- empa_gwas %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
empa_snplist <- empa_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(empa_snplist, empa_plain_snplist_, sep = "\t")



# EMPAREG GWAS (after eligibility criteria) -------------------------------
# there is no significant SNP (P > 5e-8), so clumping should produce 0 clumps / 0 IVs

empa_gwas <- fread(empa_gwas_raw_)
empa_gwas <- empa_gwas %>% dplyr::select("chromosome" = `#chrom`, "chr_position" = pos, effect_allele = alt, other_allele = ref, 
                                         "effect_weight" = beta, "se" = sebeta, pval, "af" = af_alt)

empa_gwas <- empa_gwas %>%  mutate(variant = paste(chromosome, chr_position, other_allele, effect_allele, sep = "_"),
                                   samplesize = as.numeric(11364)) 

# Add rsIDs for clumping
empa_gwas <- empa_gwas %>% mutate(chrom_pos = paste0(chromosome, "_", chr_position)) %>% 
  left_join(rsid_dict, by = c("chrom_pos" = "V2")) %>% dplyr::rename("rsid" = "V1") %>% select(-V3)

fwrite(empa_gwas, empa_gwas_)

# create List for clumping 
empa_snplist <- empa_gwas %>% filter(!(base::duplicated(rsid) & !is.na(rsid)))
empa_snplist <- empa_snplist %>% dplyr::select(SNP = rsid, P = pval) %>% filter(!(SNP == ""))
# save in txt with only SNP and P
fwrite(empa_snplist, empa_snplist_, sep = "\t")
