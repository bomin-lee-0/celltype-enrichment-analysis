#!/bin/bash
# ============================================================
# Step 1: Pre-processing and Quality Control
# - Raw FASTQ QC
# - Adapter trimming & quality filtering
# - Post-trim QC
# - MultiQC report generation
# ============================================================

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Load configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
CONFIG_FILE="${PROJECT_DIR}/config.yaml"

# Parse YAML config (simple parsing - for production use yq or python)
DATA_DIR="${PROJECT_DIR}/0_data"
OUT_DIR="${PROJECT_DIR}/01_preprocessing_qc"
THREADS=8

echo "=========================================="
echo "Starting preprocessing and QC pipeline"
echo "=========================================="
echo "Data directory: ${DATA_DIR}"
echo "Output directory: ${OUT_DIR}"
echo ""

# Check if data directory exists
if [ ! -d "${DATA_DIR}" ]; then
    echo "Error: Data directory not found: ${DATA_DIR}"
    echo "Please place your FASTQ files in ${DATA_DIR}/"
    exit 1
fi

# Create output directories
mkdir -p "${OUT_DIR}/fastqc_raw"
mkdir -p "${OUT_DIR}/fastqc_trimmed"
mkdir -p "${OUT_DIR}/trimming_logs"
mkdir -p "${OUT_DIR}/trimmed_fastq"

# ============================================================
# Step 1.1: Raw FASTQ Quality Control
# ============================================================
echo "[$(date)] Step 1.1: Running FastQC on raw reads..."

FASTQ_FILES=(${DATA_DIR}/*.fastq.gz)
if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    echo "Warning: No FASTQ files found in ${DATA_DIR}"
    echo "Please add your FASTQ files and run again"
    exit 1
fi

fastqc -t ${THREADS} ${DATA_DIR}/*.fastq.gz -o "${OUT_DIR}/fastqc_raw/" 2>&1 | tee "${OUT_DIR}/fastqc_raw.log"

echo "[$(date)] Raw FastQC complete"

# ============================================================
# Step 1.2: Adapter Trimming and Quality Filtering
# ============================================================
echo "[$(date)] Step 1.2: Running Trim Galore..."

# Process paired-end samples
for R1 in ${DATA_DIR}/*_R1.fastq.gz; do
    if [ -f "${R1}" ]; then
        BASENAME=$(basename ${R1} _R1.fastq.gz)
        R2="${DATA_DIR}/${BASENAME}_R2.fastq.gz"

        if [ -f "${R2}" ]; then
            echo "Processing paired-end sample: ${BASENAME}"

            trim_galore \
                --paired \
                --quality 20 \
                --length 20 \
                --fastqc \
                --cores 4 \
                --output_dir "${OUT_DIR}/trimmed_fastq/" \
                "${R1}" "${R2}" \
                2>&1 | tee "${OUT_DIR}/trimming_logs/${BASENAME}_trim.log"
        else
            echo "Warning: R2 file not found for ${R1}"
        fi
    fi
done

echo "[$(date)] Trimming complete"

# ============================================================
# Step 1.3: Post-trim Quality Control
# ============================================================
echo "[$(date)] Step 1.3: Running FastQC on trimmed reads..."

if [ -d "${OUT_DIR}/trimmed_fastq" ] && [ "$(ls -A ${OUT_DIR}/trimmed_fastq/*.fq.gz 2>/dev/null)" ]; then
    fastqc -t ${THREADS} ${OUT_DIR}/trimmed_fastq/*_val_*.fq.gz -o "${OUT_DIR}/fastqc_trimmed/" 2>&1 | tee "${OUT_DIR}/fastqc_trimmed.log"
    echo "[$(date)] Post-trim FastQC complete"
else
    echo "Warning: No trimmed FASTQ files found"
fi

# ============================================================
# Step 1.4: MultiQC Report
# ============================================================
echo "[$(date)] Step 1.4: Generating MultiQC report..."

multiqc "${OUT_DIR}/" -o "${OUT_DIR}/" -n "multiqc_report" --force 2>&1 | tee "${OUT_DIR}/multiqc.log"

echo "[$(date)] MultiQC report generated: ${OUT_DIR}/multiqc_report.html"

# ============================================================
# Summary
# ============================================================
echo ""
echo "=========================================="
echo "Preprocessing and QC Complete!"
echo "=========================================="
echo "Results:"
echo "  - Raw FastQC reports: ${OUT_DIR}/fastqc_raw/"
echo "  - Trimmed reads: ${OUT_DIR}/trimmed_fastq/"
echo "  - Trimmed FastQC reports: ${OUT_DIR}/fastqc_trimmed/"
echo "  - MultiQC summary: ${OUT_DIR}/multiqc_report.html"
echo ""
echo "Next step: Run alignment (02_alignment.sh)"
echo "=========================================="
