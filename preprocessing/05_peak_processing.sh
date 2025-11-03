#!/bin/bash
# ============================================================
# Step 5: Peak Processing
# - Merge replicate peaks per cell type
# - Identify cell-type-specific unique peaks
# - Generate Venn diagrams and statistics
# ============================================================

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

PEAK_DIR="${PROJECT_DIR}/03_peak_calling/macs2_output"
OUT_DIR="${PROJECT_DIR}/05_peak_processing"

echo "=========================================="
echo "Starting peak processing pipeline"
echo "=========================================="
echo "Peak directory: ${PEAK_DIR}"
echo "Output directory: ${OUT_DIR}"
echo ""

# Create output directories
mkdir -p "${OUT_DIR}/merged_peaks"
mkdir -p "${OUT_DIR}/unique_peaks"
mkdir -p "${OUT_DIR}/logs"

# ============================================================
# Define cell type groups
# Modify these based on your actual sample names
# ============================================================
declare -A CELLTYPES
CELLTYPES[Dopaminergic]="Nurr1"
CELLTYPES[Oligodendrocyte]="Olig2"
CELLTYPES[Astrocyte]="GFAP"

echo "Cell types defined:"
for CELLTYPE in "${!CELLTYPES[@]}"; do
    echo "  ${CELLTYPE}: ${CELLTYPES[$CELLTYPE]}"
done
echo ""

# ============================================================
# Step 5.1: Merge replicate peaks per cell type
# ============================================================
echo "[$(date)] Step 5.1: Merging replicate peaks per cell type..."

for CELLTYPE in "${!CELLTYPES[@]}"; do
    MARKER=${CELLTYPES[$CELLTYPE]}

    echo ""
    echo "Processing cell type: ${CELLTYPE} (marker: ${MARKER})"

    # Find all peak files for this cell type
    PEAK_FILES=(${PEAK_DIR}/${MARKER}*_peaks.narrowPeak)

    if [ ${#PEAK_FILES[@]} -eq 0 ] || [ ! -f "${PEAK_FILES[0]}" ]; then
        echo "Warning: No peak files found for ${MARKER}, skipping..."
        continue
    fi

    echo "Found ${#PEAK_FILES[@]} replicate(s) for ${MARKER}"

    # Concatenate all replicates
    cat "${PEAK_FILES[@]}" | sort -k1,1 -k2,2n > "${OUT_DIR}/merged_peaks/${MARKER}_all_peaks.bed"

    # Merge overlapping peaks
    bedtools merge \
        -i "${OUT_DIR}/merged_peaks/${MARKER}_all_peaks.bed" \
        > "${OUT_DIR}/merged_peaks/${MARKER}_merged.bed"

    PEAK_COUNT=$(wc -l < "${OUT_DIR}/merged_peaks/${MARKER}_merged.bed")
    echo "Merged peaks for ${MARKER}: ${PEAK_COUNT}"
done

echo "[$(date)] Peak merging complete"

# ============================================================
# Step 5.2: Identify cell-type-specific unique peaks
# ============================================================
echo ""
echo "[$(date)] Step 5.2: Identifying cell-type-specific unique peaks..."

# Get list of merged peak files
MERGED_FILES=(${OUT_DIR}/merged_peaks/*_merged.bed)

if [ ${#MERGED_FILES[@]} -lt 2 ]; then
    echo "Warning: Need at least 2 cell types to identify unique peaks"
    echo "Skipping unique peak identification"
else
    for CELLTYPE in "${!CELLTYPES[@]}"; do
        MARKER=${CELLTYPES[$CELLTYPE]}
        QUERY_FILE="${OUT_DIR}/merged_peaks/${MARKER}_merged.bed"

        if [ ! -f "${QUERY_FILE}" ]; then
            continue
        fi

        echo ""
        echo "Finding unique peaks for: ${MARKER}"

        # Get all other cell type files
        OTHER_FILES=""
        for OTHER_FILE in ${OUT_DIR}/merged_peaks/*_merged.bed; do
            OTHER_MARKER=$(basename ${OTHER_FILE} _merged.bed)
            if [ "${OTHER_MARKER}" != "${MARKER}" ]; then
                OTHER_FILES="${OTHER_FILES} ${OTHER_FILE}"
            fi
        done

        if [ -n "${OTHER_FILES}" ]; then
            # Concatenate other cell types
            cat ${OTHER_FILES} | sort -k1,1 -k2,2n | bedtools merge -i - \
                > "${OUT_DIR}/merged_peaks/others_for_${MARKER}.bed"

            # Find peaks unique to this cell type (no overlap with others)
            bedtools intersect \
                -a ${QUERY_FILE} \
                -b "${OUT_DIR}/merged_peaks/others_for_${MARKER}.bed" \
                -v \
                > "${OUT_DIR}/unique_peaks/${MARKER}_unique.bed"

            UNIQUE_COUNT=$(wc -l < "${OUT_DIR}/unique_peaks/${MARKER}_unique.bed")
            TOTAL_COUNT=$(wc -l < ${QUERY_FILE})
            PERCENT=$(echo "scale=2; ${UNIQUE_COUNT} * 100 / ${TOTAL_COUNT}" | bc)

            echo "  Total peaks: ${TOTAL_COUNT}"
            echo "  Unique peaks: ${UNIQUE_COUNT} (${PERCENT}%)"
        fi
    done
fi

echo "[$(date)] Unique peak identification complete"

# ============================================================
# Step 5.3: Generate peak statistics
# ============================================================
echo ""
echo "[$(date)] Step 5.3: Generating peak statistics..."

STATS_FILE="${OUT_DIR}/peak_statistics.txt"
cat > ${STATS_FILE} << EOF
========================================
Peak Processing Statistics
Generated: $(date)
========================================

Merged Peaks (per cell type):
---
EOF

for CELLTYPE in "${!CELLTYPES[@]}"; do
    MARKER=${CELLTYPES[$CELLTYPE]}
    MERGED_FILE="${OUT_DIR}/merged_peaks/${MARKER}_merged.bed"

    if [ -f "${MERGED_FILE}" ]; then
        COUNT=$(wc -l < ${MERGED_FILE})
        printf "%-20s %10d peaks\n" "${MARKER}:" "${COUNT}" >> ${STATS_FILE}
    fi
done

cat >> ${STATS_FILE} << EOF

Unique Peaks (cell-type-specific):
---
EOF

for CELLTYPE in "${!CELLTYPES[@]}"; do
    MARKER=${CELLTYPES[$CELLTYPE]}
    UNIQUE_FILE="${OUT_DIR}/unique_peaks/${MARKER}_unique.bed"

    if [ -f "${UNIQUE_FILE}" ]; then
        COUNT=$(wc -l < ${UNIQUE_FILE})
        printf "%-20s %10d peaks\n" "${MARKER}:" "${COUNT}" >> ${STATS_FILE}
    fi
done

cat >> ${STATS_FILE} << EOF

========================================
Peak Overlap Analysis:
---
EOF

# Calculate pairwise overlaps
for CELLTYPE1 in "${!CELLTYPES[@]}"; do
    MARKER1=${CELLTYPES[$CELLTYPE1]}
    FILE1="${OUT_DIR}/merged_peaks/${MARKER1}_merged.bed"

    if [ ! -f "${FILE1}" ]; then
        continue
    fi

    for CELLTYPE2 in "${!CELLTYPES[@]}"; do
        MARKER2=${CELLTYPES[$CELLTYPE2]}
        FILE2="${OUT_DIR}/merged_peaks/${MARKER2}_merged.bed"

        if [ ! -f "${FILE2}" ] || [ "${MARKER1}" = "${MARKER2}" ]; then
            continue
        fi

        OVERLAP=$(bedtools intersect -a ${FILE1} -b ${FILE2} | wc -l)
        echo "${MARKER1} ∩ ${MARKER2}: ${OVERLAP} peaks" >> ${STATS_FILE}
    done
done

cat >> ${STATS_FILE} << EOF

========================================
EOF

# Save bedtools commands for reference
cat > "${OUT_DIR}/bedtools_commands.txt" << EOF
# Bedtools commands used for peak processing
# Generated: $(date)

# Merge replicate peaks:
cat sample_rep1_peaks.narrowPeak sample_rep2_peaks.narrowPeak | sort -k1,1 -k2,2n > all_peaks.bed
bedtools merge -i all_peaks.bed > merged_peaks.bed

# Find cell-type-specific peaks:
# (peaks in cell type A that don't overlap with cell type B)
bedtools intersect -a cellA_merged.bed -b cellB_merged.bed -v > cellA_unique.bed

# Calculate peak overlaps:
bedtools intersect -a peaks1.bed -b peaks2.bed | wc -l
EOF

# ============================================================
# Summary
# ============================================================
echo ""
echo "=========================================="
echo "Peak Processing Complete!"
echo "=========================================="
echo "Results:"
echo "  - Merged peaks: ${OUT_DIR}/merged_peaks/"
echo "  - Unique peaks: ${OUT_DIR}/unique_peaks/"
echo "  - Statistics: ${OUT_DIR}/peak_statistics.txt"
echo ""
cat ${STATS_FILE}
echo ""
echo "Next step: Annotate peaks (06_annotation.R)"
echo "=========================================="
