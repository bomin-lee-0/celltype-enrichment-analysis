#!/bin/bash
# ============================================================
# Step 3: Peak Calling with MACS2
# - Call H3K27ac peaks using MACS2
# - Generate narrowPeak files
# - Create summary statistics
# ============================================================

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

INPUT_DIR="${PROJECT_DIR}/02_alignment/filtered_bam"
OUT_DIR="${PROJECT_DIR}/03_peak_calling"
GENOME_SIZE="2.5e8"  # Rat genome size
Q_VALUE=0.05

echo "=========================================="
echo "Starting peak calling pipeline"
echo "=========================================="
echo "Input directory: ${INPUT_DIR}"
echo "Output directory: ${OUT_DIR}"
echo ""

# Create output directories
mkdir -p "${OUT_DIR}/macs2_output"
mkdir -p "${OUT_DIR}/logs"

# ============================================================
# Identify input/control samples
# ============================================================
echo "[$(date)] Identifying samples..."

# Find input/control BAM files (usually named with "Input" or "IgG")
INPUT_BAMS=(${INPUT_DIR}/Input*.filtered.bam ${INPUT_DIR}/IgG*.filtered.bam 2>/dev/null || true)

if [ ${#INPUT_BAMS[@]} -eq 0 ] || [ ! -f "${INPUT_BAMS[0]}" ]; then
    echo "WARNING: No input/control BAM files found"
    echo "MACS2 will run without control"
    USE_CONTROL=false
    CONTROL_BAM=""
else
    echo "Found input/control files: ${INPUT_BAMS[@]}"
    USE_CONTROL=true
    # Use first input file as control (or combine them if multiple)
    CONTROL_BAM="${INPUT_BAMS[0]}"
fi

# ============================================================
# Call peaks for each ChIP sample
# ============================================================

for BAM in ${INPUT_DIR}/*.filtered.bam; do
    if [ -f "${BAM}" ]; then
        BASENAME=$(basename ${BAM} .filtered.bam)

        # Skip if this is an input/control file
        if [[ ${BASENAME} == Input* ]] || [[ ${BASENAME} == IgG* ]]; then
            echo "Skipping control sample: ${BASENAME}"
            continue
        fi

        echo ""
        echo "=========================================="
        echo "Calling peaks for: ${BASENAME}"
        echo "=========================================="

        # Construct MACS2 command
        MACS2_CMD="macs2 callpeak \
            -t ${BAM} \
            -f BAM \
            -g ${GENOME_SIZE} \
            -n ${BASENAME} \
            -q ${Q_VALUE} \
            --keep-dup auto \
            --outdir ${OUT_DIR}/macs2_output"

        # Add control if available
        if [ "${USE_CONTROL}" = true ]; then
            MACS2_CMD="${MACS2_CMD} -c ${CONTROL_BAM}"
        fi

        # Run MACS2
        echo "[$(date)] Running MACS2..."
        echo "Command: ${MACS2_CMD}"

        eval ${MACS2_CMD} 2>&1 | tee "${OUT_DIR}/logs/${BASENAME}_macs2.log"

        echo "[$(date)] Peak calling complete for ${BASENAME}"

        # Count peaks
        if [ -f "${OUT_DIR}/macs2_output/${BASENAME}_peaks.narrowPeak" ]; then
            PEAK_COUNT=$(wc -l < "${OUT_DIR}/macs2_output/${BASENAME}_peaks.narrowPeak")
            echo "  Peaks identified: ${PEAK_COUNT}"
        fi
    fi
done

# ============================================================
# Generate peak summary
# ============================================================
echo ""
echo "[$(date)] Generating peak calling summary..."

SUMMARY_FILE="${OUT_DIR}/peak_summary.txt"
echo "========================================" > ${SUMMARY_FILE}
echo "Peak Calling Summary" >> ${SUMMARY_FILE}
echo "Generated: $(date)" >> ${SUMMARY_FILE}
echo "========================================" >> ${SUMMARY_FILE}
echo "" >> ${SUMMARY_FILE}
echo "MACS2 parameters:" >> ${SUMMARY_FILE}
echo "  Genome size: ${GENOME_SIZE}" >> ${SUMMARY_FILE}
echo "  Q-value cutoff: ${Q_VALUE}" >> ${SUMMARY_FILE}
echo "  Control used: ${USE_CONTROL}" >> ${SUMMARY_FILE}
if [ "${USE_CONTROL}" = true ]; then
    echo "  Control file: ${CONTROL_BAM}" >> ${SUMMARY_FILE}
fi
echo "" >> ${SUMMARY_FILE}
echo "Peak counts per sample:" >> ${SUMMARY_FILE}
echo "---" >> ${SUMMARY_FILE}

for PEAK_FILE in ${OUT_DIR}/macs2_output/*_peaks.narrowPeak; do
    if [ -f "${PEAK_FILE}" ]; then
        SAMPLE=$(basename ${PEAK_FILE} _peaks.narrowPeak)
        PEAK_COUNT=$(wc -l < "${PEAK_FILE}")
        printf "%-30s %10d peaks\n" "${SAMPLE}:" "${PEAK_COUNT}" >> ${SUMMARY_FILE}
    fi
done

echo "" >> ${SUMMARY_FILE}
echo "Output files for each sample:" >> ${SUMMARY_FILE}
echo "  *_peaks.narrowPeak    : Peak locations in BED6+4 format" >> ${SUMMARY_FILE}
echo "  *_summits.bed         : Peak summit positions" >> ${SUMMARY_FILE}
echo "  *_peaks.xls           : Detailed peak information" >> ${SUMMARY_FILE}
echo "  *_model.r             : R script for peak model visualization" >> ${SUMMARY_FILE}

# Save MACS2 commands for reference
cat > "${OUT_DIR}/macs2_commands.txt" << EOF
# MACS2 peak calling commands used
# Generated: $(date)

# Parameters:
# - Format: BAM
# - Genome size: ${GENOME_SIZE}
# - Q-value: ${Q_VALUE}
# - Keep duplicates: auto
$([ "${USE_CONTROL}" = true ] && echo "# - Control: ${CONTROL_BAM}")

# Example command:
macs2 callpeak \\
  -t sample.filtered.bam \\
  $([ "${USE_CONTROL}" = true ] && echo "-c ${CONTROL_BAM} \\")
  -f BAM \\
  -g ${GENOME_SIZE} \\
  -n sample_name \\
  -q ${Q_VALUE} \\
  --keep-dup auto \\
  --outdir 03_peak_calling/macs2_output
EOF

# ============================================================
# Summary
# ============================================================
echo ""
echo "=========================================="
echo "Peak Calling Complete!"
echo "=========================================="
echo "Results:"
echo "  - Peak files: ${OUT_DIR}/macs2_output/"
echo "  - Summary: ${OUT_DIR}/peak_summary.txt"
echo "  - Logs: ${OUT_DIR}/logs/"
echo ""
cat ${SUMMARY_FILE}
echo ""
echo "Next step: Calculate QC metrics (04_qc_metrics.sh)"
echo "=========================================="
