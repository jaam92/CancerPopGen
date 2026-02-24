#!/bin/sh
#SBATCH --job-name=frequencies
#SBATCH --output=/scratch1/jazlynmo/cancerproj/logfiles/allelefreq_chr%a.out  #UPDATE TO YOUR SCRATCH
#SBATCH --error=/scratch1/jazlynmo/cancerproj/logfiles/allelefreq_chr%a.err  #UPDATE TO YOUR SCRATCH
#SBATCH --time=20:00:00 #twenty hour run time
#SBATCH -p qcbr #update the partition based on cluster you are using
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 #one cpu per task
#SBATCH --array=1-22 #number of jobs is number of chromosomes 
#SBATCH --mem-per-cpu=1G #this is equivalent to 1G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jazlynmo@usc.edu #UPDATE TO YOUR EMAIL

#load everything up
module purge
module load gcc/13.3.0 vcftools/0.1.16 #load vcftools module with dependencies first

# Variables
CHROM="${SLURM_ARRAY_TASK_ID}"

VCF_DIR="/project2/jazlynmo_738/DataRepository/Human/1000GenomeNYGC_hg38/vcfs_strictMask"
SAMPLE_ID_DIR="/project2/jazlynmo_738/DataRepository/Human/1000GenomeNYGC_hg38/FileInformation/pops"
BASE_OUTPUT_DIR="/scratch1/jazlynmo/cancerproj/alleleFrequencies" #UPDATE TO YOUR SCRATCH

VCF_FILE="${VCF_DIR}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHROM}.filtered.shapeit2-duohmm-phased.nodupmarkers.snps.strict.vcf.gz"
POPS_FILE="${SAMPLE_ID_DIR}/pops.txt"

# Chromosome-specific output directory
CHR_OUT_DIR="${BASE_OUTPUT_DIR}/chr${CHROM}"
###mkdir -p "${CHR_OUT_DIR}" #COMMENT THIS OUT IF YOU HAVE MADE THESE ALREADY AND REMOVE COMMENTS TO MAKE

echo "Starting allele frequency calculation for chr${CHROM}"
echo "Output directory: ${CHR_OUT_DIR}"


# Loop over populations
while read -r POP 
do
    echo "Processing population: ${POP}"

    vcftools \
        --gzvcf "${VCF_FILE}" \
        --freq \
        --keep "${SAMPLE_ID_DIR}/${POP}_unrels_no_3rd_N50.txt" \
        --out "${CHR_OUT_DIR}/${POP}_chr${CHROM}"

    echo "Finished ${POP}" #output the population you are on for a given chromosome

done < "${POPS_FILE}"

echo "All populations finished for chr${CHROM}"

sleep 180
