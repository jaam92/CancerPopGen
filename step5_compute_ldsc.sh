#!/bin/bash
#SBATCH --job-name=ldsc
#SBATCH --output=/scratch1/jazlynmo/cancerproj/logfiles/ldsc_chr%a.out
#SBATCH --error=/scratch1/jazlynmo/cancerproj/logfiles/ldsc_chr%a.err
#SBATCH --time=01:30:00
#SBATCH -p qcb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-22
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazlynmo@usc.edu

set -euo pipefail

############################
# Dry run control
############################
DRY_RUN=${DRY_RUN:-0}

############################
# Load environment
############################
module purge
eval "$(conda shell.bash hook)"
conda activate ldsc39

echo "Python:"
which python
python --version

############################
# Variables
############################
CHROM="${SLURM_ARRAY_TASK_ID}"

PLINK_DIR="/scratch1/jazlynmo/cancerproj/plinkFiles_ldScore/splitChroms"
LDSC_DIR="/project2/jazlynmo_738/DataRepository/Software/ldsc"
BASE_OUTPUT_DIR="/scratch1/jazlynmo/cancerproj/plinkFiles_ldScore"
POPS_FILE="/scratch1/jazlynmo/cancerproj/plinkFiles_ldScore/pops.txt"

LDSC_SCRIPT="${LDSC_DIR}/ldsc.py"

CHR_OUT_DIR="${BASE_OUTPUT_DIR}/ldsc_chr${CHROM}"

# Check if directory exists if it does then don't make it again
if [[ ! -d "${CHR_OUT_DIR}" ]]; then
    echo "Creating output directory: ${CHR_OUT_DIR}"
    mkdir -p "${CHR_OUT_DIR}"
else
    echo "Output directory already exists: ${CHR_OUT_DIR} — skipping mkdir"
fi


echo "======================================"
echo "Chromosome ${CHROM}"
echo "Output dir: ${CHR_OUT_DIR}"
echo "======================================"

############################
# Loop populations
############################
while read -r POP || [[ -n "$POP" ]]
do
    [[ -z "$POP" || "$POP" =~ ^# ]] && continue

    BFILE="${PLINK_DIR}/${POP}_Unrels_strictMask_chr${CHROM}"

    echo "---- Population: ${POP}"

    if [[ ! -f "${BFILE}.bed" ]]; then
        echo "Missing ${BFILE}.bed — skipping"
        continue
    fi

    CMD=(python "${LDSC_SCRIPT}"
        --bfile "${BFILE}"
        --l2
        --ld-wind-cm 1
        --out "${CHR_OUT_DIR}/${POP}_Unrels_strictMask_chr${CHROM}"
    )

    echo "${CMD[@]}"

    if [ "$DRY_RUN" -eq 0 ]; then
        "${CMD[@]}"
    else
        echo "[DRY RUN]"
    fi

done < "${POPS_FILE}"

echo "Finished chromosome ${CHROM}"

