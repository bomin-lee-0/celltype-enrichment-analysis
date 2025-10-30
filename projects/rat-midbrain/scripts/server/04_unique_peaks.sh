#!/bin/bash
# ============================================================
# Cell Type-Specific Unique Peaks Identification
# Finds peaks unique to each cell type (no overlap with others)
# ============================================================

cd /scratch/prj/bcn_marzi_lab/ratlas/Bomin

echo "=========================================="
echo "Identifying cell type-specific unique peaks"
echo "=========================================="

# Cell type별 narrowPeak 파일
NEUN="macs2_output/NeuN/NeuN_all_peaks.narrowPeak"
NURR="macs2_output/Nurr/Nurr_all_peaks.narrowPeak"
OLIG="macs2_output/Olig/Olig_all_peaks.narrowPeak"
NEG="macs2_output/Neg/Neg_all_peaks.narrowPeak"

# 출력 폴더 생성
mkdir -p merged_unique

echo ""
echo "Input files:"
for peak_file in "$NEUN" "$NURR" "$OLIG" "$NEG"; do
    if [ -f "$peak_file" ]; then
        COUNT=$(wc -l < "$peak_file")
        echo "  ✓ $(basename $peak_file): ${COUNT} peaks"
    else
        echo "  ✗ $(basename $peak_file): NOT FOUND"
    fi
done

# ============================================================
# 전체 병합 (All cell types merged)
# ============================================================
echo ""
echo "[$(date)] Merging all peaks across cell types..."

cat "$NEUN" "$NURR" "$OLIG" "$NEG" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - > merged_unique/merged_all.bed

ALL_MERGED=$(wc -l < merged_unique/merged_all.bed)
echo "  All merged peaks: ${ALL_MERGED}"

# ============================================================
# NeuN 고유 피크
# ============================================================
echo ""
echo "[$(date)] Finding unique peaks for NeuN..."
echo "  Excluding: Nurr, Olig, Neg"

cat "$NURR" "$OLIG" "$NEG" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - > merged_unique/tmp.bed

bedtools intersect -a "$NEUN" -b merged_unique/tmp.bed -v \
  > merged_unique/unique_NeuN.bed

NEUN_UNIQUE=$(wc -l < merged_unique/unique_NeuN.bed)
NEUN_TOTAL=$(wc -l < "$NEUN")
NEUN_PCT=$(awk -v u="$NEUN_UNIQUE" -v t="$NEUN_TOTAL" 'BEGIN {printf "%.1f", u*100/t}')
echo "  NeuN unique: ${NEUN_UNIQUE} / ${NEUN_TOTAL} (${NEUN_PCT}%)"

# ============================================================
# Nurr 고유 피크
# ============================================================
echo ""
echo "[$(date)] Finding unique peaks for Nurr..."
echo "  Excluding: NeuN, Olig, Neg"

cat "$NEUN" "$OLIG" "$NEG" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - > merged_unique/tmp.bed

bedtools intersect -a "$NURR" -b merged_unique/tmp.bed -v \
  > merged_unique/unique_Nurr.bed

NURR_UNIQUE=$(wc -l < merged_unique/unique_Nurr.bed)
NURR_TOTAL=$(wc -l < "$NURR")
NURR_PCT=$(awk -v u="$NURR_UNIQUE" -v t="$NURR_TOTAL" 'BEGIN {printf "%.1f", u*100/t}')
echo "  Nurr unique: ${NURR_UNIQUE} / ${NURR_TOTAL} (${NURR_PCT}%)"

# ============================================================
# Olig 고유 피크
# ============================================================
echo ""
echo "[$(date)] Finding unique peaks for Olig..."
echo "  Excluding: NeuN, Nurr, Neg"

cat "$NEUN" "$NURR" "$NEG" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - > merged_unique/tmp.bed

bedtools intersect -a "$OLIG" -b merged_unique/tmp.bed -v \
  > merged_unique/unique_Olig.bed

OLIG_UNIQUE=$(wc -l < merged_unique/unique_Olig.bed)
OLIG_TOTAL=$(wc -l < "$OLIG")
OLIG_PCT=$(awk -v u="$OLIG_UNIQUE" -v t="$OLIG_TOTAL" 'BEGIN {printf "%.1f", u*100/t}')
echo "  Olig unique: ${OLIG_UNIQUE} / ${OLIG_TOTAL} (${OLIG_PCT}%)"

# ============================================================
# Neg 고유 피크
# ============================================================
echo ""
echo "[$(date)] Finding unique peaks for Neg..."
echo "  Excluding: NeuN, Nurr, Olig"

cat "$NEUN" "$NURR" "$OLIG" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - > merged_unique/tmp.bed

bedtools intersect -a "$NEG" -b merged_unique/tmp.bed -v \
  > merged_unique/unique_Neg.bed

NEG_UNIQUE=$(wc -l < merged_unique/unique_Neg.bed)
NEG_TOTAL=$(wc -l < "$NEG")
NEG_PCT=$(awk -v u="$NEG_UNIQUE" -v t="$NEG_TOTAL" 'BEGIN {printf "%.1f", u*100/t}')
echo "  Neg unique: ${NEG_UNIQUE} / ${NEG_TOTAL} (${NEG_PCT}%)"

# Clean up temporary file
rm merged_unique/tmp.bed

# ============================================================
# Summary
# ============================================================
echo ""
echo "=========================================="
echo "Unique Peak Identification Complete"
echo "=========================================="
echo ""
echo "Summary:"
echo "---"
printf "%-10s %10s %10s %10s\n" "Cell Type" "Total" "Unique" "% Unique"
printf "%-10s %10s %10s %10s\n" "----------" "------" "------" "--------"
printf "%-10s %10d %10d %9.1f%%\n" "NeuN" "$NEUN_TOTAL" "$NEUN_UNIQUE" "$NEUN_PCT"
printf "%-10s %10d %10d %9.1f%%\n" "Nurr" "$NURR_TOTAL" "$NURR_UNIQUE" "$NURR_PCT"
printf "%-10s %10d %10d %9.1f%%\n" "Olig" "$OLIG_TOTAL" "$OLIG_UNIQUE" "$OLIG_PCT"
printf "%-10s %10d %10d %9.1f%%\n" "Neg" "$NEG_TOTAL" "$NEG_UNIQUE" "$NEG_PCT"
echo ""
echo "All merged: ${ALL_MERGED} peaks"
echo ""
echo "Output files:"
echo "  - merged_unique/unique_NeuN.bed"
echo "  - merged_unique/unique_Nurr.bed"
echo "  - merged_unique/unique_Olig.bed"
echo "  - merged_unique/unique_Neg.bed"
echo "  - merged_unique/merged_all.bed"
echo "=========================================="
