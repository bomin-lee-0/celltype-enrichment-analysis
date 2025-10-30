#!/bin/bash
# ============================================================
# FRiP Score Calculation (Fragment-level)
# Final successful version from server
# ============================================================

cd /scratch/prj/bcn_marzi_lab/ratlas/Bomin

echo "=========================================="
echo "Starting FRiP score calculation"
echo "Fragment-level analysis"
echo "=========================================="

# Cell type별 MACS2 peak 파일 지정
NEUN="macs2_output/NeuN/NeuN_all_peaks.narrowPeak"
NURR="macs2_output/Nurr/Nurr_all_peaks.narrowPeak"
OLIG="macs2_output/Olig/Olig_all_peaks.narrowPeak"
NEG="macs2_output/Neg/Neg_all_peaks.narrowPeak"

# Check if peak files exist
echo "Checking peak files..."
for peak_file in "$NEUN" "$NURR" "$OLIG" "$NEG"; do
    if [ -f "$peak_file" ]; then
        echo "  ✓ $(basename $peak_file)"
    else
        echo "  ✗ $(basename $peak_file) NOT FOUND"
    fi
done
echo ""

# 출력 파일 초기화
echo "sample,frip_score" > frip_scores.csv

# alignment_copy 내 filtered BAM 기준 (fragment 단위)
echo "Processing BAM files..."
echo ""

SAMPLE_COUNT=0
for bam in alignment_copy/*.bam; do
    sample=$(basename "$bam" .bam)
    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))

    echo "=========================================="
    echo "Sample ${SAMPLE_COUNT}: ${sample}"
    echo "=========================================="

    # Cell type 자동 매칭
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
    elif [[ "$sample" == *"input"* ]] || [[ "$sample" == *"Input"* ]]; then
        echo "  Skipping input/control sample"
        echo ""
        continue
    else
        echo "  Could not determine cell type, skipping"
        echo ""
        continue
    fi

    echo "  Cell type: ${celltype}"

    # BAM fragment 수 계산
    echo "  Counting total fragments..."
    total_fragments=$(samtools view -c "$bam")
    echo "  Total fragments: ${total_fragments}"

    # Peaks 안에 포함된 fragment 수 계산 (BEDTools)
    echo "  Counting fragments in peaks..."
    fragments_in_peaks=$(bedtools intersect -a "$bam" -b "$peaks" -bed | wc -l)
    echo "  Fragments in peaks: ${fragments_in_peaks}"

    # FRiP 계산
    if [ ${total_fragments} -gt 0 ]; then
        frip=$(awk -v a="$fragments_in_peaks" -v b="$total_fragments" 'BEGIN {if (b>0) print a/b; else print 0}')
        echo "  FRiP score: ${frip}"

        # Determine QC status
        frip_percent=$(awk -v f="$frip" 'BEGIN {printf "%.1f%%", f*100}')
        if (( $(echo "$frip >= 0.2" | bc -l) )); then
            status="PASS"
        elif (( $(echo "$frip >= 0.1" | bc -l) )); then
            status="WARNING"
        else
            status="FAIL"
        fi
        echo "  QC Status: ${status} (${frip_percent})"

        # Save to CSV
        echo "${sample},${frip}" >> frip_scores.csv
    else
        echo "  ERROR: No fragments found"
        echo "${sample},0" >> frip_scores.csv
    fi

    echo ""
done

# 정렬 후 덮어쓰기 (최종)
echo "Sorting results..."
sort -t',' -k1,1 frip_scores.csv -o frip_scores.csv

# ============================================================
# Summary
# ============================================================
echo ""
echo "=========================================="
echo "FRiP Calculation Complete"
echo "=========================================="
echo "Results saved to: frip_scores.csv"
echo ""
echo "Summary by cell type:"
echo "---"

for CELLTYPE in NeuN Nurr Olig Neg; do
    echo ""
    echo "${CELLTYPE}:"
    grep -E "IGF.*" frip_scores.csv | while IFS=',' read sample frip; do
        case "$sample" in
            IGF13135*) ct="NeuN" ;;
            IGF13136[67]*) ct="Nurr" ;;
            IGF13136[012]*) ct="Olig" ;;
            IGF13136[9]*|IGF13137[01]*) ct="Neg" ;;
            *) ct="" ;;
        esac

        if [ "$ct" = "$CELLTYPE" ]; then
            frip_pct=$(awk -v f="$frip" 'BEGIN {printf "%.2f%%", f*100}')
            printf "  %-15s %s\n" "$sample" "$frip_pct"
        fi
    done
done

echo ""
echo "=========================================="
