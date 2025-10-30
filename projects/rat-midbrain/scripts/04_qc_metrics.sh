#!/bin/bash
# ============================================================
# Step 4: Calculate QC Metrics
# - Calculate FRiP (Fraction of Reads in Peaks) scores
# - Generate QC summary report
# ============================================================

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

BAM_DIR="${PROJECT_DIR}/02_alignment/filtered_bam"
PEAK_DIR="${PROJECT_DIR}/03_peak_calling/macs2_output"
OUT_DIR="${PROJECT_DIR}/04_qc_metrics"
MIN_FRIP=0.2

echo "=========================================="
echo "Starting QC metrics calculation"
echo "=========================================="
echo "BAM directory: ${BAM_DIR}"
echo "Peak directory: ${PEAK_DIR}"
echo "Output directory: ${OUT_DIR}"
echo ""

# Create output directory
mkdir -p "${OUT_DIR}"

# Initialize output file
FRIP_FILE="${OUT_DIR}/frip_scores.csv"
echo "Sample,Total_Reads,Reads_in_Peaks,FRiP_Score,QC_Status" > ${FRIP_FILE}

# ============================================================
# Calculate FRiP for each sample
# ============================================================

for BAM in ${BAM_DIR}/*.filtered.bam; do
    if [ -f "${BAM}" ]; then
        BASENAME=$(basename ${BAM} .filtered.bam)

        # Skip control samples
        if [[ ${BASENAME} == Input* ]] || [[ ${BASENAME} == IgG* ]]; then
            echo "Skipping control sample: ${BASENAME}"
            continue
        fi

        # Check if peak file exists
        PEAK_FILE="${PEAK_DIR}/${BASENAME}_peaks.narrowPeak"
        if [ ! -f "${PEAK_FILE}" ]; then
            echo "Warning: Peak file not found for ${BASENAME}, skipping..."
            continue
        fi

        echo ""
        echo "=========================================="
        echo "Calculating FRiP for: ${BASENAME}"
        echo "=========================================="

        # --------------------------------------------------------
        # Count total reads
        # --------------------------------------------------------
        echo "[$(date)] Counting total reads..."
        TOTAL_READS=$(samtools view -c ${BAM})
        echo "Total reads: ${TOTAL_READS}"

        # --------------------------------------------------------
        # Count reads in peaks
        # --------------------------------------------------------
        echo "[$(date)] Counting reads in peaks..."
        READS_IN_PEAKS=$(bedtools intersect -a ${BAM} -b ${PEAK_FILE} -bed | wc -l)
        echo "Reads in peaks: ${READS_IN_PEAKS}"

        # --------------------------------------------------------
        # Calculate FRiP score
        # --------------------------------------------------------
        if [ ${TOTAL_READS} -gt 0 ]; then
            FRIP=$(echo "scale=4; ${READS_IN_PEAKS} / ${TOTAL_READS}" | bc)
            echo "FRiP score: ${FRIP}"

            # Determine QC status
            PASS=$(echo "${FRIP} >= ${MIN_FRIP}" | bc)
            if [ ${PASS} -eq 1 ]; then
                QC_STATUS="PASS"
                echo "QC Status: PASS (FRiP >= ${MIN_FRIP})"
            else
                QC_STATUS="WARNING"
                echo "QC Status: WARNING (FRiP < ${MIN_FRIP})"
            fi
        else
            FRIP="0.0000"
            QC_STATUS="FAIL"
            echo "QC Status: FAIL (no reads)"
        fi

        # Save to CSV
        echo "${BASENAME},${TOTAL_READS},${READS_IN_PEAKS},${FRIP},${QC_STATUS}" >> ${FRIP_FILE}
    fi
done

# ============================================================
# Generate summary report
# ============================================================
echo ""
echo "[$(date)] Generating QC summary..."

SUMMARY_FILE="${OUT_DIR}/qc_summary.txt"
cat > ${SUMMARY_FILE} << EOF
========================================
ChIP-seq QC Summary Report
Generated: $(date)
========================================

FRiP Score Threshold: >= ${MIN_FRIP}

Detailed Results:
---
EOF

# Read and format CSV data
tail -n +2 ${FRIP_FILE} | while IFS=',' read -r sample total reads_in_peaks frip status; do
    cat >> ${SUMMARY_FILE} << EOF

Sample: ${sample}
  Total reads:       ${total}
  Reads in peaks:    ${reads_in_peaks}
  FRiP score:        ${frip}
  QC status:         ${status}
EOF
done

cat >> ${SUMMARY_FILE} << EOF

---
Summary Statistics:
EOF

# Calculate summary statistics
TOTAL_SAMPLES=$(tail -n +2 ${FRIP_FILE} | wc -l)
PASS_SAMPLES=$(tail -n +2 ${FRIP_FILE} | grep -c "PASS" || echo "0")
WARNING_SAMPLES=$(tail -n +2 ${FRIP_FILE} | grep -c "WARNING" || echo "0")
FAIL_SAMPLES=$(tail -n +2 ${FRIP_FILE} | grep -c "FAIL" || echo "0")

cat >> ${SUMMARY_FILE} << EOF

Total samples analyzed: ${TOTAL_SAMPLES}
  PASS (FRiP >= ${MIN_FRIP}):    ${PASS_SAMPLES}
  WARNING (FRiP < ${MIN_FRIP}):  ${WARNING_SAMPLES}
  FAIL (no reads):               ${FAIL_SAMPLES}

========================================

Notes:
- FRiP (Fraction of Reads in Peaks) measures enrichment quality
- High-quality ChIP-seq typically has FRiP >= 0.2 (20%)
- Samples with FRiP < 0.1 may indicate failed experiments
- Consider removing low-quality samples from downstream analysis

========================================
EOF

# ============================================================
# Summary
# ============================================================
echo ""
echo "=========================================="
echo "QC Metrics Complete!"
echo "=========================================="
echo "Results:"
echo "  - FRiP scores: ${FRIP_FILE}"
echo "  - QC summary: ${SUMMARY_FILE}"
echo ""
cat ${SUMMARY_FILE}
echo ""
echo "Next step: Process peaks (05_peak_processing.sh)"
echo "=========================================="
