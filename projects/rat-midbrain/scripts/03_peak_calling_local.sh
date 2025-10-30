#!/bin/bash
# ============================================================
# Step 3: Local Peak Calling with MACS2
# Calls peaks by combining replicates for each cell type
# Matches original server script exactly
# ============================================================

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

BAM_DIR="${PROJECT_DIR}/02_alignment/filtered_bam"
OUT_DIR="${PROJECT_DIR}/03_peak_calling"
GENOME="2.5e8"

echo "=========================================="
echo "Starting local peak calling"
echo "=========================================="
echo "BAM directory: ${BAM_DIR}"
echo "Output directory: ${OUT_DIR}"
echo "Genome size: ${GENOME}"
echo ""

# Check if macs2 is available
if ! command -v macs2 &> /dev/null; then
    echo "ERROR: macs2 not found!"
    echo "Please activate conda environment: conda activate ratlas_env"
    exit 1
fi

# Create output directories
mkdir -p "${OUT_DIR}/macs2_output/NeuN"
mkdir -p "${OUT_DIR}/macs2_output/Nurr"
mkdir -p "${OUT_DIR}/macs2_output/Olig"
mkdir -p "${OUT_DIR}/macs2_output/Neg"
mkdir -p "${OUT_DIR}/logs"

# Check if BAM files exist
if [ ! -d "${BAM_DIR}" ]; then
    echo "ERROR: BAM directory not found: ${BAM_DIR}"
    echo "Please run alignment step first"
    exit 1
fi

# ============================================================
# Peak Calling: NeuN
# ============================================================
echo ""
echo "=========================================="
echo "Calling peaks for: NeuN"
echo "=========================================="

if [ -f "${BAM_DIR}/IGF131357.filtered.bam" ] && \
   [ -f "${BAM_DIR}/IGF131358.filtered.bam" ] && \
   [ -f "${BAM_DIR}/IGF131359.filtered.bam" ] && \
   [ -f "${BAM_DIR}/IGF131373.filtered.bam" ]; then

    macs2 callpeak \
        -t "${BAM_DIR}/IGF131357.filtered.bam" \
           "${BAM_DIR}/IGF131358.filtered.bam" \
           "${BAM_DIR}/IGF131359.filtered.bam" \
        -c "${BAM_DIR}/IGF131373.filtered.bam" \
        -f BAM -g ${GENOME} -n NeuN_all \
        --outdir "${OUT_DIR}/macs2_output/NeuN" \
        --nomodel --shift 0 --extsize 200 -B --SPMR \
        2>&1 | tee "${OUT_DIR}/logs/NeuN_macs2.log"

    echo "[$(date)] NeuN peak calling complete"
    PEAK_COUNT=$(wc -l < "${OUT_DIR}/macs2_output/NeuN/NeuN_all_peaks.narrowPeak")
    echo "Peaks identified: ${PEAK_COUNT}"
else
    echo "ERROR: NeuN BAM files not found"
fi

# ============================================================
# Peak Calling: Nurr
# ============================================================
echo ""
echo "=========================================="
echo "Calling peaks for: Nurr"
echo "=========================================="

if [ -f "${BAM_DIR}/IGF131366.filtered.bam" ] && \
   [ -f "${BAM_DIR}/IGF131367.filtered.bam" ] && \
   [ -f "${BAM_DIR}/IGF131376.filtered.bam" ]; then

    macs2 callpeak \
        -t "${BAM_DIR}/IGF131366.filtered.bam" \
           "${BAM_DIR}/IGF131367.filtered.bam" \
        -c "${BAM_DIR}/IGF131376.filtered.bam" \
        -f BAM -g ${GENOME} -n Nurr_all \
        --outdir "${OUT_DIR}/macs2_output/Nurr" \
        --nomodel --shift 0 --extsize 200 -B --SPMR \
        2>&1 | tee "${OUT_DIR}/logs/Nurr_macs2.log"

    echo "[$(date)] Nurr peak calling complete"
    PEAK_COUNT=$(wc -l < "${OUT_DIR}/macs2_output/Nurr/Nurr_all_peaks.narrowPeak")
    echo "Peaks identified: ${PEAK_COUNT}"
else
    echo "ERROR: Nurr BAM files not found"
fi

# ============================================================
# Peak Calling: Olig
# ============================================================
echo ""
echo "=========================================="
echo "Calling peaks for: Olig"
echo "=========================================="

if [ -f "${BAM_DIR}/IGF131360.filtered.bam" ] && \
   [ -f "${BAM_DIR}/IGF131361.filtered.bam" ] && \
   [ -f "${BAM_DIR}/IGF131362.filtered.bam" ] && \
   [ -f "${BAM_DIR}/IGF131374.filtered.bam" ]; then

    macs2 callpeak \
        -t "${BAM_DIR}/IGF131360.filtered.bam" \
           "${BAM_DIR}/IGF131361.filtered.bam" \
           "${BAM_DIR}/IGF131362.filtered.bam" \
        -c "${BAM_DIR}/IGF131374.filtered.bam" \
        -f BAM -g ${GENOME} -n Olig_all \
        --outdir "${OUT_DIR}/macs2_output/Olig" \
        --nomodel --shift 0 --extsize 200 -B --SPMR \
        2>&1 | tee "${OUT_DIR}/logs/Olig_macs2.log"

    echo "[$(date)] Olig peak calling complete"
    PEAK_COUNT=$(wc -l < "${OUT_DIR}/macs2_output/Olig/Olig_all_peaks.narrowPeak")
    echo "Peaks identified: ${PEAK_COUNT}"
else
    echo "ERROR: Olig BAM files not found"
fi

# ============================================================
# Peak Calling: Neg
# ============================================================
echo ""
echo "=========================================="
echo "Calling peaks for: Neg"
echo "=========================================="

if [ -f "${BAM_DIR}/IGF131369.filtered.bam" ] && \
   [ -f "${BAM_DIR}/IGF131370.filtered.bam" ] && \
   [ -f "${BAM_DIR}/IGF131371.filtered.bam" ] && \
   [ -f "${BAM_DIR}/IGF131377.filtered.bam" ]; then

    macs2 callpeak \
        -t "${BAM_DIR}/IGF131369.filtered.bam" \
           "${BAM_DIR}/IGF131370.filtered.bam" \
           "${BAM_DIR}/IGF131371.filtered.bam" \
        -c "${BAM_DIR}/IGF131377.filtered.bam" \
        -f BAM -g ${GENOME} -n Neg_all \
        --outdir "${OUT_DIR}/macs2_output/Neg" \
        --nomodel --shift 0 --extsize 200 -B --SPMR \
        2>&1 | tee "${OUT_DIR}/logs/Neg_macs2.log"

    echo "[$(date)] Neg peak calling complete"
    PEAK_COUNT=$(wc -l < "${OUT_DIR}/macs2_output/Neg/Neg_all_peaks.narrowPeak")
    echo "Peaks identified: ${PEAK_COUNT}"
else
    echo "ERROR: Neg BAM files not found"
fi

# ============================================================
# Generate summary
# ============================================================
echo ""
echo "[$(date)] Generating peak calling summary..."

SUMMARY_FILE="${OUT_DIR}/peak_summary.txt"
cat > ${SUMMARY_FILE} << EOF
========================================
Peak Calling Summary
Generated: $(date)
========================================

MACS2 parameters:
  Genome size: ${GENOME}
  --nomodel --shift 0 --extsize 200
  -B --SPMR

Peak counts by cell type:
---
EOF

for CELLTYPE in NeuN Nurr Olig Neg; do
    PEAK_FILE="${OUT_DIR}/macs2_output/${CELLTYPE}/${CELLTYPE}_all_peaks.narrowPeak"
    if [ -f "${PEAK_FILE}" ]; then
        PEAK_COUNT=$(wc -l < "${PEAK_FILE}")
        printf "%-20s: %10d peaks\n" "${CELLTYPE}" "${PEAK_COUNT}" >> ${SUMMARY_FILE}
    else
        printf "%-20s: Not found\n" "${CELLTYPE}" >> ${SUMMARY_FILE}
    fi
done

cat >> ${SUMMARY_FILE} << EOF

Output files per cell type:
  *_all_peaks.narrowPeak    : Peak locations
  *_all_summits.bed         : Peak summits
  *_all_peaks.xls           : Detailed information
  *_all_treat_pileup.bdg    : Treatment bedGraph
  *_all_control_lambda.bdg  : Control bedGraph

========================================
EOF

# ============================================================
# Summary
# ============================================================
echo ""
echo "=========================================="
echo "Peak Calling Complete!"
echo "=========================================="
cat ${SUMMARY_FILE}
echo ""
echo "Next step: Calculate QC metrics (04_qc_metrics_local.sh)"
echo "=========================================="
