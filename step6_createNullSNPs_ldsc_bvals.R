library(tidyverse)
library(vroom)
library(data.table)

# -----------------------------
# Load static data ONCE
# -----------------------------
gwas_catalog <- vroom(
  "~/sc1/cancerproj/processGWAScatalog/processed_gwas-associations.tsv",
  delim = "\t",
  col_select = c(SNP_ID_CURRENT, cancer_type, allele)
)

rsids <- vroom(
  "~/sc1/cancerproj/processGWAScatalog/rsid_positions_hg38.txt",
  delim = "\t",
  col_select = c(chr, start, end, rsID, alleles)
) %>%
  mutate(
    cancer_type = gwas_catalog$cancer_type[
      match(rsID, gwas_catalog$SNP_ID_CURRENT)
    ],
    associated_allele = gwas_catalog$allele[
      match(rsID, gwas_catalog$SNP_ID_CURRENT)
    ]
  )

bvals <- vroom(
  "~/sc1/cancerproj/snps_in_vcf_bvals_neutralGr800Lt980_noCancerVars.bed"
)

setDT(rsids)
setDT(bvals)

# -----------------------------
# Iterate over chromosomes
# -----------------------------
chromosomes <- 1:22

for (chr in chromosomes) {
  chr_dir <- paste0(
    "~/sc1/cancerproj/plinkFiles_ldScore/ldsc_chr", chr
  )
  message("Processing chromosome ", chr)
  
  ld_files <- list.files(
    path = chr_dir,
    pattern = paste0("_chr", chr, "\\.l2\\.ldscore\\.gz$"),
    full.names = TRUE
  )
  
  for (file in ld_files) {
    pop <- basename(file) %>%
      str_remove(paste0("_Unrels_strictMask_chr", chr, "\\.l2\\.ldscore\\.gz"))
    message("  Population: ", pop)
    
    # -----------------------------
    # Load LDSC data (only needed columns)
    # -----------------------------
    ldsc <- vroom(file, col_select = c(CHR, BP, L2)) %>%
      mutate(CHR = paste0("chr", CHR))
    setDT(ldsc)
    
    # GWAS SNPs
    ldsc_GWAS <- ldsc[rsids, on = .(CHR == chr, BP = end), nomatch=0]
    
    # Null pool
    ldsc_null <- ldsc[!rsids, on = .(CHR == chr, BP = end)][bvals, on = .(CHR == chr, BP = end), nomatch=0]
    
    rm(ldsc); gc()
    
    if (nrow(ldsc_GWAS) == 0 | nrow(ldsc_null) == 0) {
      message("    Skipping (no GWAS or null SNPs)")
      rm(ldsc_GWAS, ldsc_null); gc()
      next
    }
    
    # -----------------------------
    # Rolling join matching
    # -----------------------------
    ldsc_GWAS[, gwas_id := .I]
    ldsc_null[, null_id := .I]
    
    setkey(ldsc_GWAS, CHR, L2)
    setkey(ldsc_null, CHR, L2)
    
    candidate_pairs <- ldsc_null[ldsc_GWAS, 
                                 .(CHR, BP_gwas = i.BP, L2_gwas = i.L2,
                                   BP_null = BP, L2_null = L2,
                                   gwas_id = i.gwas_id, null_id = null_id,
                                   cancer_type = i.cancer_type),
                                 roll = 5]
    
    if (nrow(candidate_pairs) == 0) {
      message("    No matches within threshold")
      rm(ldsc_GWAS, ldsc_null, candidate_pairs); gc()
      next
    }
    
    candidate_pairs[, diff := abs(L2_gwas - L2_null)]
    
    # Pick top 10 per GWAS SNP
    setorder(candidate_pairs, gwas_id, diff)
    top_matches <- candidate_pairs[, head(.SD, 10), by = gwas_id]
    
    # Ensure each null SNP used once
    setorder(top_matches, diff)
    matched <- top_matches[, head(.SD, 1), by = null_id]
    
    # -----------------------------
    # Split by cancer_type, make folders, and save
    # -----------------------------
    cancer_types <- unique(matched$cancer_type)
    
    for (cancer in cancer_types) {
      # sanitize cancer_type for folder and filenames
      cancer_clean <- cancer %>%
        str_replace_all("[^A-Za-z0-9]", "_")
      
      # create folder if it doesn't exist
      if (!dir.exists(cancer_clean)) dir.create(cancer_clean)
      
      matched_subset <- matched[cancer_type == cancer]
      matched_bed <- matched_subset[, .(CHR, start = BP_null, end = BP_null)]
      
      # Filenames include pop, cancer_type, and chromosome
      matched_file <- file.path(cancer_clean, paste0(pop, "_", cancer_clean, "_chr", chr, "_null_matched_bvals_ldsc.txt"))
      bed_file <- file.path(cancer_clean, paste0(pop, "_", cancer_clean, "_chr", chr, "_null_matched_bvals_ldsc.bed"))
      
      write_tsv(matched_subset, matched_file)
      write_tsv(matched_bed, bed_file)
    }
    
    message("    Files written for population ", pop, " chromosome ", chr)
    
    # Cleanup memory
    rm(ldsc_GWAS, ldsc_null, candidate_pairs, top_matches, matched); gc()
  }
}
