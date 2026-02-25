#!/bin/sh
#SBATCH --job-name=frequencies
#SBATCH --output=/scratch1/jazlynmo/cancerproj/logfiles/allelefreq_chr%a.out
#SBATCH --error=/scratch1/jazlynmo/cancerproj/logfiles/allelefreq_chr%a.err
#SBATCH --time=30:00:00
#SBATCH -p qcbr
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-22
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazlynmo@usc.edu

module purge
module load gcc/13.3.0 vcftools/0.1.16

CHROM="${SLURM_ARRAY_TASK_ID}"
#CHROM=22

VCF_DIR="/project2/jazlynmo_738/DataRepository/Human/1000GenomeNYGC_hg38/vcfs_strictMask"
SAMPLE_ID_DIR="/project2/jazlynmo_738/DataRepository/Human/1000GenomeNYGC_hg38/FileInformation/pops"
REF_DIR="/project2/jazlynmo_738/DataRepository/Reference_Files/CancerGenomics"
BASE_OUTPUT_DIR="/scratch1/jazlynmo/cancerproj/alleleFrequencies"

VCF_FILE="${VCF_DIR}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHROM}.filtered.shapeit2-duohmm-phased.nodupmarkers.snps.strict.vcf.gz"
POPS_FILE="${SAMPLE_ID_DIR}/pops.txt"

# Chromosome directory, make if needed
CHR_OUT_DIR="${BASE_OUTPUT_DIR}/chr${CHROM}"

if [ ! -d "${CHR_OUT_DIR}" ]; then
    echo "Creating directory ${CHR_OUT_DIR}"
    mkdir -p "${CHR_OUT_DIR}"
fi

echo "Starting allele frequency calculation for chr${CHROM}"

# Cancer types
CANCER_TYPES=("breastCancer" "colorectalCancer" "kidneyCancer" "lungCancer" "ovarianCancer" "pancreaticCancer" "prostateCancer")

# Loop over cancer types
for CANCER in "${CANCER_TYPES[@]}"
do
    echo "Processing cancer type: ${CANCER}"

    # Cancer-specific directory, make if needed
    CANCER_OUT_DIR="${CHR_OUT_DIR}/${CANCER}"

    if [ ! -d "${CANCER_OUT_DIR}" ]; then
        echo "Creating directory ${CANCER_OUT_DIR}"
        mkdir -p "${CANCER_OUT_DIR}"
    fi

    # Loop over populations
    while read -r POP
    do
        echo "Processing population: ${POP}"

        vcftools \
            --gzvcf "${VCF_FILE}" \
            --freq \
            --bed "${REF_DIR}/gwasSNPs/${CANCER}/${CANCER}_chr${CHROM}.bed" \
            --keep "${SAMPLE_ID_DIR}/${POP}_unrels_no_3rd_N50.txt" \
            --out "${CANCER_OUT_DIR}/${POP}_chr${CHROM}_${CANCER}_gwas"

        vcftools \
            --gzvcf "${VCF_FILE}" \
            --freq \
            --bed "${REF_DIR}/nullSNPs_bval_ldsc/${CANCER}/${POP}_${CANCER}_chr${CHROM}_null_matched_bvals_ldsc.bed" \
            --keep "${SAMPLE_ID_DIR}/${POP}_unrels_no_3rd_N50.txt" \
            --out "${CANCER_OUT_DIR}/${POP}_chr${CHROM}_${CANCER}_null"

        echo "Finished ${POP} for ${CANCER}"

    done < "${POPS_FILE}"

done

echo "All cancer types and populations finished for chr${CHROM}"
