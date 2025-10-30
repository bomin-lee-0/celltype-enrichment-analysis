#!/bin/bash
# ============================================================
# Non-intersect Merged Peaks
# Merges all peaks across cell types without intersection
# ============================================================

cd /scratch/prj/bcn_marzi_lab/ratlas/Bomin

echo "=========================================="
echo "Creating non-intersect merged peaks"
echo "=========================================="

# Cell type별 narrowPeak 파일
NEUN="macs2_output/NeuN/NeuN_all_peaks.narrowPeak"
NURR="macs2_output/Nurr/Nurr_all_peaks.narrowPeak"
OLIG="macs2_output/Olig/Olig_all_peaks.narrowPeak"
NEG="macs2_output/Neg/Neg_all_peaks.narrowPeak"

# 출력 폴더
mkdir -p merged_unique

echo ""
echo "Input files:"
TOTAL_PEAKS=0
for peak_file in "$NEUN" "$NURR" "$OLIG" "$NEG"; do
    if [ -f "$peak_file" ]; then
        COUNT=$(wc -l < "$peak_file")
        TOTAL_PEAKS=$((TOTAL_PEAKS + COUNT))
        echo "  ✓ $(basename $peak_file): ${COUNT} peaks"
    else
        echo "  ✗ $(basename $peak_file): NOT FOUND"
    fi
done
echo "  Total peaks (before merge): ${TOTAL_PEAKS}"

# ============================================================
# Non-intersect (all-cell merged) 버전
# ============================================================
echo ""
echo "[$(date)] Merging all peaks (non-intersect)..."

cat "$NEUN" "$NURR" "$OLIG" "$NEG" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - > merged_unique/all_merged_peaks.bed

MERGED_COUNT=$(wc -l < merged_unique/all_merged_peaks.bed)

echo "[$(date)] Merging complete"

# ============================================================
# Summary
# ============================================================
echo ""
echo "=========================================="
echo "Non-intersect Merged Peaks Complete"
echo "=========================================="
echo ""
echo "Summary:"
echo "  Total peaks before merge: ${TOTAL_PEAKS}"
echo "  Merged peaks: ${MERGED_COUNT}"
echo "  Reduction: $((TOTAL_PEAKS - MERGED_COUNT)) peaks"
echo "  Merge rate: $(awk -v m="$MERGED_COUNT" -v t="$TOTAL_PEAKS" 'BEGIN {printf "%.1f%%", m*100/t}')"
echo ""
echo "Output file:"
echo "  merged_unique/all_merged_peaks.bed"
echo "=========================================="
