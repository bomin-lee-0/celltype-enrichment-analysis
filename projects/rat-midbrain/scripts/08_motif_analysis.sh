#!/bin/bash
# ============================================================
# Step 8: Motif Analysis with HOMER
# - De novo motif discovery in unique peaks
# - Known motif enrichment analysis
# - Generate motif reports
# ============================================================

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

PEAK_DIR="${PROJECT_DIR}/05_peak_processing/unique_peaks"
OUT_DIR="${PROJECT_DIR}/08_motif_analysis"
GENOME="rn7"
SIZE="given"  # Use actual peak sizes

echo "=========================================="
echo "Starting motif analysis pipeline"
echo "=========================================="
echo "Peak directory: ${PEAK_DIR}"
echo "Output directory: ${OUT_DIR}"
echo "Genome: ${GENOME}"
echo ""

# Create output directory
mkdir -p "${OUT_DIR}/motif_results"
mkdir -p "${OUT_DIR}/logs"

# ============================================================
# Check for genome files
# ============================================================
echo "[$(date)] Checking for HOMER genome..."

# Check if HOMER genome is installed
if ! homer2 listGenomes 2>/dev/null | grep -q "${GENOME}"; then
    echo "WARNING: HOMER genome ${GENOME} not found"
    echo "To install, run: homer2 installGenome ${GENOME}"
    echo ""
    echo "For now, continuing with available genome..."
    echo "You may need to adjust the genome parameter"
fi

# ============================================================
# Run motif discovery for each cell type
# ============================================================

# Find all unique peak files
PEAK_FILES=(${PEAK_DIR}/*_unique.bed)

if [ ${#PEAK_FILES[@]} -eq 0 ] || [ ! -f "${PEAK_FILES[0]}" ]; then
    echo "Error: No unique peak files found in ${PEAK_DIR}"
    echo "Please run 05_peak_processing.sh first"
    exit 1
fi

echo "Found ${#PEAK_FILES[@]} peak file(s) to analyze"
echo ""

for PEAK_FILE in "${PEAK_FILES[@]}"; do
    if [ -f "${PEAK_FILE}" ]; then
        BASENAME=$(basename ${PEAK_FILE} _unique.bed)

        echo ""
        echo "=========================================="
        echo "Analyzing motifs for: ${BASENAME}"
        echo "=========================================="

        PEAK_COUNT=$(wc -l < ${PEAK_FILE})
        echo "Number of peaks: ${PEAK_COUNT}"

        if [ ${PEAK_COUNT} -lt 10 ]; then
            echo "Warning: Too few peaks (${PEAK_COUNT}) for meaningful motif analysis"
            echo "Skipping ${BASENAME}"
            continue
        fi

        OUTPUT_SUBDIR="${OUT_DIR}/motif_results/${BASENAME}"
        mkdir -p ${OUTPUT_SUBDIR}

        # --------------------------------------------------------
        # Run HOMER motif discovery
        # --------------------------------------------------------
        echo "[$(date)] Running HOMER findMotifsGenome.pl..."

        findMotifsGenome.pl \
            ${PEAK_FILE} \
            ${GENOME} \
            ${OUTPUT_SUBDIR} \
            -size ${SIZE} \
            -mask \
            -p 4 \
            2>&1 | tee "${OUT_DIR}/logs/${BASENAME}_homer.log"

        echo "[$(date)] Motif discovery complete for ${BASENAME}"

        # Check results
        if [ -f "${OUTPUT_SUBDIR}/homerResults.html" ]; then
            echo "  Results saved: ${OUTPUT_SUBDIR}/homerResults.html"

            # Count discovered motifs
            if [ -f "${OUTPUT_SUBDIR}/homerMotifs.all.motifs" ]; then
                MOTIF_COUNT=$(grep -c "^>" "${OUTPUT_SUBDIR}/homerMotifs.all.motifs" || echo "0")
                echo "  De novo motifs discovered: ${MOTIF_COUNT}"
            fi

            # Count known motif matches
            if [ -f "${OUTPUT_SUBDIR}/knownResults.txt" ]; then
                KNOWN_COUNT=$(tail -n +2 "${OUTPUT_SUBDIR}/knownResults.txt" | wc -l)
                echo "  Known motifs analyzed: ${KNOWN_COUNT}"
            fi
        else
            echo "  Warning: HOMER results not found"
        fi
    fi
done

# ============================================================
# Generate motif summary
# ============================================================
echo ""
echo "[$(date)] Generating motif analysis summary..."

SUMMARY_FILE="${OUT_DIR}/motif_summary.txt"
cat > ${SUMMARY_FILE} << EOF
========================================
Motif Analysis Summary
Generated: $(date)
========================================

HOMER Parameters:
  Genome: ${GENOME}
  Size: ${SIZE}
  Masking: enabled

Results by Cell Type:
---
EOF

for PEAK_FILE in "${PEAK_FILES[@]}"; do
    if [ -f "${PEAK_FILE}" ]; then
        BASENAME=$(basename ${PEAK_FILE} _unique.bed)
        OUTPUT_SUBDIR="${OUT_DIR}/motif_results/${BASENAME}"

        echo "" >> ${SUMMARY_FILE}
        echo "Sample: ${BASENAME}" >> ${SUMMARY_FILE}

        PEAK_COUNT=$(wc -l < ${PEAK_FILE})
        echo "  Peaks analyzed: ${PEAK_COUNT}" >> ${SUMMARY_FILE}

        if [ -f "${OUTPUT_SUBDIR}/homerMotifs.all.motifs" ]; then
            MOTIF_COUNT=$(grep -c "^>" "${OUTPUT_SUBDIR}/homerMotifs.all.motifs" || echo "0")
            echo "  De novo motifs: ${MOTIF_COUNT}" >> ${SUMMARY_FILE}
        else
            echo "  De novo motifs: N/A" >> ${SUMMARY_FILE}
        fi

        if [ -f "${OUTPUT_SUBDIR}/knownResults.txt" ]; then
            echo "  Top known motifs:" >> ${SUMMARY_FILE}
            head -n 6 "${OUTPUT_SUBDIR}/knownResults.txt" | tail -n 5 | awk '{print "    - " $1}' >> ${SUMMARY_FILE}
        fi

        if [ -f "${OUTPUT_SUBDIR}/homerResults.html" ]; then
            echo "  Report: ${OUTPUT_SUBDIR}/homerResults.html" >> ${SUMMARY_FILE}
        fi
    fi
done

cat >> ${SUMMARY_FILE} << EOF

========================================
Key Motifs to Look For:

Based on the rat midbrain context, expect to see:
  - Nurr1 (NR4A2): Dopaminergic neuron marker
  - FOXA2: Midbrain development
  - PU.1 (SPI1): Microglial marker
  - OLIG2: Oligodendrocyte marker
  - GFAP-related TFs: Astrocyte markers

========================================
EOF

# Save HOMER commands for reference
cat > "${OUT_DIR}/homer_commands.txt" << EOF
# HOMER motif analysis commands used
# Generated: $(date)

# Parameters:
# - Genome: ${GENOME}
# - Size: ${SIZE} (use actual peak sizes)
# - Masking: enabled (mask repetitive sequences)
# - Threads: 4

# Example command:
findMotifsGenome.pl \\
  05_peak_processing/unique_peaks/Nurr1_unique.bed \\
  ${GENOME} \\
  08_motif_analysis/motif_results/Nurr1/ \\
  -size ${SIZE} \\
  -mask \\
  -p 4

# Output files:
# - homerResults.html: Main results page
# - homerMotifs.all.motifs: All de novo motifs
# - knownResults.txt: Known motif enrichment
# - motif1.logo.png, motif2.logo.png, ...: Motif logos

# To install HOMER genome:
# perl /path/to/homer/configureHomer.pl -install ${GENOME}
EOF

# ============================================================
# Summary
# ============================================================
echo ""
echo "=========================================="
echo "Motif Analysis Complete!"
echo "=========================================="
echo "Results:"
echo "  - Motif results: ${OUT_DIR}/motif_results/"
echo "  - Summary: ${OUT_DIR}/motif_summary.txt"
echo "  - Logs: ${OUT_DIR}/logs/"
echo ""
echo "To view results, open the HTML files:"
for PEAK_FILE in "${PEAK_FILES[@]}"; do
    if [ -f "${PEAK_FILE}" ]; then
        BASENAME=$(basename ${PEAK_FILE} _unique.bed)
        HTML_FILE="${OUT_DIR}/motif_results/${BASENAME}/homerResults.html"
        if [ -f "${HTML_FILE}" ]; then
            echo "  - ${HTML_FILE}"
        fi
    fi
done
echo ""
cat ${SUMMARY_FILE}
echo ""
echo "=========================================="
echo "Pipeline Complete!"
echo "=========================================="
