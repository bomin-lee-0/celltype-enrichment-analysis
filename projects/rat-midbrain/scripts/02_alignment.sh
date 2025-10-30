#!/bin/bash
# ============================================================
# Step 2: Read Alignment
# - Align trimmed reads to rat genome (rn7) using Bowtie2
# - Sort and index BAM files
# - Mark and remove duplicates
# - Filter for high-quality alignments
# - Generate alignment statistics
# ============================================================

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

INPUT_DIR="${PROJECT_DIR}/01_preprocessing_qc/trimmed_fastq"
OUT_DIR="${PROJECT_DIR}/02_alignment"
THREADS=8
MAX_FRAGMENT=2000

echo "=========================================="
echo "Starting alignment pipeline"
echo "=========================================="
echo "Input directory: ${INPUT_DIR}"
echo "Output directory: ${OUT_DIR}"
echo ""

# Create output directories
mkdir -p "${OUT_DIR}/aligned_bam"
mkdir -p "${OUT_DIR}/dedup_bam"
mkdir -p "${OUT_DIR}/filtered_bam"
mkdir -p "${OUT_DIR}/logs"

# ============================================================
# Check for Bowtie2 index
# ============================================================
echo "[$(date)] Checking for Bowtie2 reference index..."

# You need to set this to your actual Bowtie2 index path
# Example: BOWTIE2_INDEX="/path/to/rn7_index/rn7"
BOWTIE2_INDEX="${BOWTIE2_INDEX:-/path/to/rn7_index/rn7}"

if [ ! -f "${BOWTIE2_INDEX}.1.bt2" ]; then
    echo "ERROR: Bowtie2 index not found!"
    echo "Please set BOWTIE2_INDEX environment variable or update this script"
    echo "Example: export BOWTIE2_INDEX=/path/to/rn7_index/rn7"
    echo ""
    echo "To build the index, run:"
    echo "  bowtie2-build rn7.fa rn7"
    exit 1
fi

echo "Using Bowtie2 index: ${BOWTIE2_INDEX}"

# ============================================================
# Process each sample
# ============================================================

# Find all R1 trimmed files
for R1 in ${INPUT_DIR}/*_R1_val_1.fq.gz; do
    if [ -f "${R1}" ]; then
        # Extract sample name
        BASENAME=$(basename ${R1} _R1_val_1.fq.gz)
        R2="${INPUT_DIR}/${BASENAME}_R2_val_2.fq.gz"

        if [ ! -f "${R2}" ]; then
            echo "Warning: R2 file not found for ${BASENAME}, skipping..."
            continue
        fi

        echo ""
        echo "=========================================="
        echo "Processing sample: ${BASENAME}"
        echo "=========================================="

        # --------------------------------------------------------
        # Step 2.1: Alignment with Bowtie2
        # --------------------------------------------------------
        echo "[$(date)] Step 2.1: Aligning reads with Bowtie2..."

        bowtie2 \
            -p ${THREADS} \
            --very-sensitive \
            --maxins ${MAX_FRAGMENT} \
            -x ${BOWTIE2_INDEX} \
            -1 ${R1} \
            -2 ${R2} \
            2> "${OUT_DIR}/logs/${BASENAME}_bowtie2.log" \
            | samtools view -@ ${THREADS} -bS - \
            | samtools sort -@ ${THREADS} -o "${OUT_DIR}/aligned_bam/${BASENAME}.sorted.bam" -

        samtools index "${OUT_DIR}/aligned_bam/${BASENAME}.sorted.bam"

        echo "[$(date)] Alignment complete: ${BASENAME}.sorted.bam"

        # --------------------------------------------------------
        # Step 2.2: Mark and remove duplicates
        # --------------------------------------------------------
        echo "[$(date)] Step 2.2: Removing duplicates..."

        samtools markdup \
            -r \
            -@ ${THREADS} \
            "${OUT_DIR}/aligned_bam/${BASENAME}.sorted.bam" \
            "${OUT_DIR}/dedup_bam/${BASENAME}.dedup.bam" \
            2> "${OUT_DIR}/logs/${BASENAME}_markdup.log"

        samtools index "${OUT_DIR}/dedup_bam/${BASENAME}.dedup.bam"

        echo "[$(date)] Duplicates removed: ${BASENAME}.dedup.bam"

        # --------------------------------------------------------
        # Step 2.3: Filter for high-quality alignments
        # --------------------------------------------------------
        echo "[$(date)] Step 2.3: Filtering for high-quality alignments..."

        # Filter:
        # -F 1804: exclude unmapped, secondary, qc-fail, duplicate, supplementary
        # -f 2: include properly paired
        # -q 20: minimum MAPQ 20
        samtools view \
            -@ ${THREADS} \
            -b \
            -F 1804 \
            -f 2 \
            -q 20 \
            "${OUT_DIR}/dedup_bam/${BASENAME}.dedup.bam" \
            > "${OUT_DIR}/filtered_bam/${BASENAME}.filtered.bam"

        samtools index "${OUT_DIR}/filtered_bam/${BASENAME}.filtered.bam"

        echo "[$(date)] Filtering complete: ${BASENAME}.filtered.bam"

        # --------------------------------------------------------
        # Step 2.4: Generate alignment statistics
        # --------------------------------------------------------
        echo "[$(date)] Step 2.4: Generating alignment statistics..."

        samtools flagstat "${OUT_DIR}/aligned_bam/${BASENAME}.sorted.bam" \
            > "${OUT_DIR}/logs/${BASENAME}_aligned_flagstat.txt"

        samtools flagstat "${OUT_DIR}/dedup_bam/${BASENAME}.dedup.bam" \
            > "${OUT_DIR}/logs/${BASENAME}_dedup_flagstat.txt"

        samtools flagstat "${OUT_DIR}/filtered_bam/${BASENAME}.filtered.bam" \
            > "${OUT_DIR}/logs/${BASENAME}_filtered_flagstat.txt"

        echo "[$(date)] Statistics saved"
    fi
done

# ============================================================
# Generate summary report
# ============================================================
echo ""
echo "[$(date)] Generating alignment summary..."

SUMMARY_FILE="${OUT_DIR}/alignment_summary.txt"
echo "========================================" > ${SUMMARY_FILE}
echo "Alignment Summary" >> ${SUMMARY_FILE}
echo "Generated: $(date)" >> ${SUMMARY_FILE}
echo "========================================" >> ${SUMMARY_FILE}
echo "" >> ${SUMMARY_FILE}

for FLAGSTAT in ${OUT_DIR}/logs/*_filtered_flagstat.txt; do
    if [ -f "${FLAGSTAT}" ]; then
        SAMPLE=$(basename ${FLAGSTAT} _filtered_flagstat.txt)
        echo "Sample: ${SAMPLE}" >> ${SUMMARY_FILE}
        echo "---" >> ${SUMMARY_FILE}
        cat ${FLAGSTAT} >> ${SUMMARY_FILE}
        echo "" >> ${SUMMARY_FILE}
    fi
done

# Save alignment commands for reference
cat > "${OUT_DIR}/bowtie2_alignment_commands.txt" << EOF
# Bowtie2 alignment commands used
# Generated: $(date)

# Alignment parameters:
# - Index: ${BOWTIE2_INDEX}
# - Mode: --very-sensitive
# - Max fragment length: ${MAX_FRAGMENT}
# - Threads: ${THREADS}

# Filtering parameters:
# - Remove unmapped, secondary, QC-fail, duplicates, supplementary
# - Require properly paired
# - Minimum MAPQ: 20

# Example command:
bowtie2 -p ${THREADS} --very-sensitive --maxins ${MAX_FRAGMENT} \\
    -x ${BOWTIE2_INDEX} \\
    -1 sample_R1_val_1.fq.gz \\
    -2 sample_R2_val_2.fq.gz \\
    | samtools view -bS - \\
    | samtools sort -o sample.sorted.bam -

samtools markdup -r sample.sorted.bam sample.dedup.bam
samtools view -b -F 1804 -f 2 -q 20 sample.dedup.bam > sample.filtered.bam
EOF

# ============================================================
# Summary
# ============================================================
echo ""
echo "=========================================="
echo "Alignment Complete!"
echo "=========================================="
echo "Results:"
echo "  - Aligned BAMs: ${OUT_DIR}/aligned_bam/"
echo "  - Deduplicated BAMs: ${OUT_DIR}/dedup_bam/"
echo "  - Filtered BAMs: ${OUT_DIR}/filtered_bam/"
echo "  - Alignment summary: ${OUT_DIR}/alignment_summary.txt"
echo "  - Logs: ${OUT_DIR}/logs/"
echo ""
echo "Next step: Run peak calling (03_peak_calling.sh)"
echo "=========================================="
