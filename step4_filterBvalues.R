# -----------------------------
# Load necessary packages
# -----------------------------
library(tidyverse)
library(vroom, verbose = FALSE)

# -----------------------------
# Load necessary functions
# -----------------------------

`%!in%` = Negate(`%in%`)

# -----------------------------
# Read input files
# -----------------------------
setwd("~/sc1/cancerproj")

#all variants in the catalog
gwas_catalog = read_delim("~/sc1/cancerproj/processGWAScatalog/processed_gwas-associations.tsv", delim = "\t")

#rsIDs of variants in gwas catalog
rsids = read_delim("~/sc1/cancerproj/processGWAScatalog/rsid_positions_hg38.txt", delim = "\t") %>%
  select(chr, start, end, rsID, alleles) %>%
  mutate(
    cancer_type = gwas_catalog$cancer_type[match(rsID, gwas_catalog$SNP_ID_CURRENT)],
    associated_allele = gwas_catalog$allele[match(rsID, gwas_catalog$SNP_ID_CURRENT)]
  ) 

#all sites in the vcf with mcvicker b-values on hg38 coordinates
sites_in_vcf = read_delim("~/sc1/cancerproj/snps_in_vcf_bvals.bed", delim = "\t", col_names = c("chr", "start", "end", "bval")) %>%
  filter(end%!in%rsids$end) %>% #remove the sites that are associated with cancer
  mutate(bval = bval*1000) %>% #rescale b-value
  filter(bval > 800 & bval < 980) #keep neutral b-values that are between 800 and 980 (remove the top 2% bc bvals aren't well calibrated https://elifesciences.org/articles/76065#data)

write_tsv(sites_in_vcf, file = "snps_in_vcf_bvals_neutralGr800Lt980_noCancerVars.bed") #write out
