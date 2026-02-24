###ASK FOR 20G FOR 4 HOURS TO RUN (was done in about 3hrs)

# -----------------------------
# Load necessary packages
# -----------------------------
library(tidyverse)
library(vroom, verbose = FALSE)

# -----------------------------
# Read input files
# -----------------------------
setwd("~/sc1/cancerproj")

gwas_catalog = read_delim("~/sc1/cancerproj/processGWAScatalog/processed_gwas-associations.tsv", delim = "\t")
rsids = read_delim("~/sc1/cancerproj/processGWAScatalog/rsid_positions_hg38.txt", delim = "\t") %>%
  mutate(
    cancer_type = gwas_catalog$cancer_type[match(rsID, gwas_catalog$SNP_ID_CURRENT)],
    allele = gwas_catalog$allele[match(rsID, gwas_catalog$SNP_ID_CURRENT)]
  )

# -----------------------------
# Directories & parameters
# -----------------------------
base_dir = "~/sc1/cancerproj/alleleFrequencies" 
chroms = paste0("chr", 1:22)
pop_pattern = "\\.frq$"
out_dir = "~/sc1/cancerproj/processedAlleleFrequencies" 
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(23)  # For reproducible null SNP sampling

# -----------------------------
# Loop over chromosomes
# -----------------------------
for(chr in chroms){
  
  chr_path = file.path(base_dir, chr)
  chr_out = file.path(out_dir, chr)
  dir.create(chr_out, showWarnings = FALSE)
  
  frq_files = list.files(chr_path, pattern = pop_pattern, full.names = TRUE)
  
  # -----------------------------
  # Generate 100k null SNPs once per chromosome
  # -----------------------------
  df_first = vroom(frq_files[1], delim = "\t", col_types = cols(
    CHROM = col_character(),
    POS = col_integer(),
    N_ALLELES = col_integer(),
    N_CHR = col_integer(),
    `{ALLELE:FREQ}` = col_character()
  ),
  show_col_types = FALSE
  ) %>%
    separate_wider_delim(`{ALLELE:FREQ}`, delim="\t", names=c("af1","af2")) %>%
    separate_wider_delim(af1, delim=":", names=c("REF","p")) %>%
    separate_wider_delim(af2, delim=":", names=c("ALT","q")) %>%
    mutate(p = as.numeric(p), q = as.numeric(q))
  
  non_gwas_hits = df_first %>%
    anti_join(rsids, by = join_by(CHROM == chr, POS == end))
  
  null_snps_sample = non_gwas_hits %>%
    slice_sample(n = min(100000, nrow(non_gwas_hits))) %>%
    pull(POS)
  
  # -----------------------------
  # Loop over population files
  # -----------------------------
  for(f in frq_files){
    
    pop_name = basename(f) %>% str_remove("_chr[0-9]+\\.frq$")
    message("Processing ", chr, " - ", pop_name)
    
    df = vroom(f, delim = "\t", col_types = cols(
      CHROM = col_character(),
      POS = col_integer(),
      N_ALLELES = col_integer(),
      N_CHR = col_integer(),
      `{ALLELE:FREQ}` = col_character()
    ),
    show_col_types = FALSE
    ) %>%
      separate_wider_delim(`{ALLELE:FREQ}`, delim="\t", names=c("af1","af2")) %>%
      separate_wider_delim(af1, delim=":", names=c("REF","p")) %>%
      separate_wider_delim(af2, delim=":", names=c("ALT","q")) %>%
      mutate(p = as.numeric(p), q = as.numeric(q), population = pop_name)
    
    # -----------------------------
    # GWAS hits
    # -----------------------------
    af_gwasCatalog = df %>%
      inner_join(rsids, by = join_by(CHROM == chr, POS == end)) %>%
      select(-c(N_ALLELES, N_CHR, start))
    
    # -----------------------------
    # Null SNPs (100k per chromosome)
    # -----------------------------
    af_nullSNPs = df %>%
      filter(POS %in% null_snps_sample) %>%
      select(-c(N_ALLELES, N_CHR))
    
    # -----------------------------
    # Save outputs
    # -----------------------------
    write_tsv(af_gwasCatalog, file.path(chr_out, paste0(pop_name, "_gwas.txt")))
    write_tsv(af_nullSNPs, file.path(chr_out, paste0(pop_name, "_null.txt")))
  }
  
  message("Done with ", chr)
}
