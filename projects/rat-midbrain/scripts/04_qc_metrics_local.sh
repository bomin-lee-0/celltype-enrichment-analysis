#!/bin/bash
# ============================================================
# Step 4: Local FRiP Score Calculation
# Fragment-level FRiP calculation matching original server script
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

echo "=========================================="
echo "Starting local FRiP score calculation"
echo "=========================================="
echo "BAM directory: ${BAM_DIR}"
echo "Peak directory: ${PEAK_DIR}"
echo "Output directory: ${OUT_DIR}"
echo ""

# Check required tools
if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools not found!"
    echo "Please activate conda environment: conda activate ratlas_env"
    exit 1
fi

if ! command -v bedtools &> /dev/null; then
    echo "ERROR: bedtools not found!"
    echo "Please activate conda environment: conda activate ratlas_env"
    exit 1
fi

# Create output directory
mkdir -p "${OUT_DIR}"

# ============================================================
# Define peak files for each cell type
# ============================================================
NEUN="${PEAK_DIR}/NeuN/NeuN_all_peaks.narrowPeak"
NURR="${PEAK_DIR}/Nurr/Nurr_all_peaks.narrowPeak"
OLIG="${PEAK_DIR}/Olig/Olig_all_peaks.narrowPeak"
NEG="${PEAK_DIR}/Neg/Neg_all_peaks.narrowPeak"

# Check if peak files exist
if [ ! -f "${NEUN}" ]; then echo "Warning: NeuN peaks not found"; fi
if [ ! -f "${NURR}" ]; then echo "Warning: Nurr peaks not found"; fi
if [ ! -f "${OLIG}" ]; then echo "Warning: Olig peaks not found"; fi
if [ ! -f "${NEG}" ]; then echo "Warning: Neg peaks not found"; fi

# ============================================================
# Calculate FRiP for each sample
# ============================================================
echo ""
echo "[$(date)] Calculating FRiP scores (fragment-level)..."
echo ""

# Initialize output file
FRIP_FILE="${OUT_DIR}/frip_scores.csv"
echo "sample,total_fragments,fragments_in_peaks,frip_score,celltype" > ${FRIP_FILE}

SAMPLE_COUNT=0
SUCCESS_COUNT=0

# Process all BAM files
for bam in ${BAM_DIR}/*.filtered.bam; do
    if [ ! -f "${bam}" ]; then
        continue
    fi

    sample=$(basename "$bam" .filtered.bam)
    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))

    echo "=========================================="
    echo "Sample ${SAMPLE_COUNT}: ${sample}"
    echo "=========================================="

    # Auto-match cell type based on sample ID
    peaks=""
    celltype=""

    if [[ "$sample" == IGF13135* ]]; then
        peaks="$NEUN"
        celltype="NeuN"
    elif [[ "$sample" == IGF13136[67]* ]]; then
        peaks="$NURR"
        celltype="Nurr"
    elif [[ "$sample" == IGF13136[012]* ]]; then
        peaks="$OLIG"
        celltype="Olig"
    elif [[ "$sample" == IGF13136[9]* ]] || [[ "$sample" == IGF13137[01]* ]]; then
        peaks="$NEG"
        celltype="Neg"
    elif [[ "$sample" == IGF13137[3467]* ]]; then
        echo "  Skipping input/control sample"
        continue
    else
        echo "  Warning: Could not determine cell type, skipping"
        continue
    fi

    echo "  Cell type: ${celltype}"

    # Check if peak file exists
    if [ ! -f "${peaks}" ]; then
        echo "  ERROR: Peak file not found for ${celltype}"
        continue
    fi

    # Count total fragments
    echo "  Counting total fragments..."
    total_fragments=$(samtools view -c "$bam")
    echo "  Total fragments: ${total_fragments}"

    # Count fragments in peaks
    echo "  Counting fragments in peaks..."
    fragments_in_peaks=$(bedtools intersect -a "$bam" -b "$peaks" -bed | wc -l)
    echo "  Fragments in peaks: ${fragments_in_peaks}"

    # Calculate FRiP
    if [ ${total_fragments} -gt 0 ]; then
        frip=$(awk -v a="$fragments_in_peaks" -v b="$total_fragments" 'BEGIN {printf "%.4f", a/b}')
        echo "  FRiP score: ${frip}"

        # Save to CSV
        echo "${sample},${total_fragments},${fragments_in_peaks},${frip},${celltype}" >> ${FRIP_FILE}
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        echo "  ERROR: No fragments found"
        echo "${sample},0,0,0.0000,${celltype}" >> ${FRIP_FILE}
    fi

    echo ""
done

# Sort the CSV by sample name
sort -t',' -k1,1 ${FRIP_FILE} -o ${FRIP_FILE}

# ============================================================
# Generate summary report
# ============================================================
echo "[$(date)] Generating FRiP summary..."

SUMMARY_FILE="${OUT_DIR}/qc_summary.txt"
cat > ${SUMMARY_FILE} << EOF
========================================
FRiP Score QC Summary
Generated: $(date)
========================================

FRiP (Fraction of Reads in Peaks):
  - Measures ChIP enrichment quality
  - Good quality: FRiP >= 0.2 (20%)
  - Acceptable: FRiP >= 0.1 (10%)
  - Poor: FRiP < 0.1 (<10%)

Results by Cell Type:
---

EOF

# Generate cell type summaries
for CELLTYPE in NeuN Nurr Olig Neg; do
    echo "${CELLTYPE}:" >> ${SUMMARY_FILE}

    # Extract samples for this cell type
    grep ",${CELLTYPE}$" ${FRIP_FILE} 2>/dev/null | while IFS=',' read -r sample total in_peaks frip ct; do
        if [ "$sample" != "sample" ]; then
            status="PASS"
            frip_float=$(echo "$frip" | awk '{print $1+0}')
            if (( $(echo "$frip_float < 0.2" | bc -l) )); then
                if (( $(echo "$frip_float < 0.1" | bc -l) )); then
                    status="FAIL"
                else
                    status="WARNING"
                fi
            fi
            printf "  %-15s FRiP: %6s  [%s]\n" "$sample" "$frip" "$status" >> ${SUMMARY_FILE}
        fi
    done

    # Calculate mean FRiP for cell type
    MEAN_FRIP=$(grep ",${CELLTYPE}$" ${FRIP_FILE} 2>/dev/null | tail -n +2 | awk -F',' '{sum+=$4; count++} END {if(count>0) printf "%.4f", sum/count; else print "N/A"}')
    echo "  Mean FRiP: ${MEAN_FRIP}" >> ${SUMMARY_FILE}
    echo "" >> ${SUMMARY_FILE}
done

cat >> ${SUMMARY_FILE} << EOF
========================================
Summary Statistics:
EOF

TOTAL=$(tail -n +2 ${FRIP_FILE} | wc -l)
PASS=$(tail -n +2 ${FRIP_FILE} | awk -F',' '$4 >= 0.2 {count++} END {print count+0}')
WARNING=$(tail -n +2 ${FRIP_FILE} | awk -F',' '$4 >= 0.1 && $4 < 0.2 {count++} END {print count+0}')
FAIL=$(tail -n +2 ${FRIP_FILE} | awk -F',' '$4 < 0.1 {count++} END {print count+0}')

cat >> ${SUMMARY_FILE} << EOF

Total samples analyzed: ${TOTAL}
  PASS (FRiP >= 0.2):    ${PASS}
  WARNING (0.1-0.2):     ${WARNING}
  FAIL (FRiP < 0.1):     ${FAIL}

========================================
EOF

# ============================================================
# Summary
# ============================================================
echo ""
echo "=========================================="
echo "FRiP Calculation Complete!"
echo "=========================================="
echo "Samples processed: ${SAMPLE_COUNT}"
echo "Successfully calculated: ${SUCCESS_COUNT}"
echo ""
echo "Results:"
echo "  - FRiP scores: ${FRIP_FILE}"
echo "  - Summary: ${SUMMARY_FILE}"
echo ""
cat ${SUMMARY_FILE}
echo ""
echo "Next step: Process peaks (05_peak_processing.sh)"
echo "=========================================="
