#!/bin/bash
# ============================================================
# Step 1: Local Preprocessing and QC with trim_galore
# Converted from SLURM array job to local sequential processing
# ============================================================

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

DATA_DIR="${PROJECT_DIR}/0_data"
OUT_DIR="${PROJECT_DIR}/01_preprocessing_qc"
THREADS=8

echo "=========================================="
echo "Starting local preprocessing and QC"
echo "=========================================="
echo "Data directory: ${DATA_DIR}"
echo "Output directory: ${OUT_DIR}"
echo "Threads: ${THREADS}"
echo ""

# Check if conda environment is activated
if ! command -v trim_galore &> /dev/null; then
    echo "ERROR: trim_galore not found!"
    echo "Please activate conda environment: conda activate ratlas_env"
    exit 1
fi

# Create output directories
mkdir -p "${OUT_DIR}/trimmed_fastq"
mkdir -p "${OUT_DIR}/fastqc"
mkdir -p "${OUT_DIR}/logs"

# Check if data directory exists
if [ ! -d "${DATA_DIR}" ]; then
    echo "ERROR: Data directory not found: ${DATA_DIR}"
    echo "Please place your FASTQ files in ${DATA_DIR}/"
    exit 1
fi

# ============================================================
# Process all paired-end samples
# ============================================================
echo "[$(date)] Processing samples..."
echo ""

SAMPLE_COUNT=0
SUCCESS_COUNT=0

# Loop through all R1 files
for INPUT_R1 in ${DATA_DIR}/*_R1.fastq.gz; do
    if [ ! -f "${INPUT_R1}" ]; then
        echo "Warning: No R1 FASTQ files found"
        continue
    fi

    # Extract sample name
    SAMPLE=$(basename ${INPUT_R1} _R1.fastq.gz)
    INPUT_R2="${DATA_DIR}/${SAMPLE}_R2.fastq.gz"

    if [ ! -f "${INPUT_R2}" ]; then
        echo "Warning: R2 file not found for ${SAMPLE}, skipping..."
        continue
    fi

    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))

    echo "=========================================="
    echo "Processing sample ${SAMPLE_COUNT}: ${SAMPLE}"
    echo "=========================================="
    echo "R1: ${INPUT_R1}"
    echo "R2: ${INPUT_R2}"

    # Run trim_galore
    echo "[$(date)] Running trim_galore..."

    trim_galore \
        --paired \
        --quality 20 \
        --length 20 \
        --fastqc \
        --fastqc_args "--outdir ${OUT_DIR}/fastqc" \
        --cores 4 \
        -o "${OUT_DIR}/trimmed_fastq" \
        "${INPUT_R1}" "${INPUT_R2}" \
        2>&1 | tee "${OUT_DIR}/logs/${SAMPLE}_trim.log"

    if [ $? -eq 0 ]; then
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        echo "[$(date)] Sample ${SAMPLE} completed successfully"
    else
        echo "[$(date)] ERROR: Sample ${SAMPLE} failed"
    fi

    echo ""
done

# ============================================================
# Generate MultiQC report
# ============================================================
if command -v multiqc &> /dev/null; then
    echo "[$(date)] Generating MultiQC report..."
    multiqc "${OUT_DIR}/" -o "${OUT_DIR}/" -n "multiqc_report" --force 2>&1 | tee "${OUT_DIR}/logs/multiqc.log"
    echo "MultiQC report: ${OUT_DIR}/multiqc_report.html"
else
    echo "Warning: multiqc not found, skipping report generation"
fi

# ============================================================
# Summary
# ============================================================
echo ""
echo "=========================================="
echo "Preprocessing Complete!"
echo "=========================================="
echo "Total samples processed: ${SAMPLE_COUNT}"
echo "Successfully completed: ${SUCCESS_COUNT}"
echo "Failed: $((SAMPLE_COUNT - SUCCESS_COUNT))"
echo ""
echo "Results:"
echo "  - Trimmed reads: ${OUT_DIR}/trimmed_fastq/"
echo "  - FastQC reports: ${OUT_DIR}/fastqc/"
echo "  - Logs: ${OUT_DIR}/logs/"
if [ -f "${OUT_DIR}/multiqc_report.html" ]; then
    echo "  - MultiQC report: ${OUT_DIR}/multiqc_report.html"
fi
echo ""
echo "Next step: Run alignment (02_alignment.sh)"
echo "=========================================="
