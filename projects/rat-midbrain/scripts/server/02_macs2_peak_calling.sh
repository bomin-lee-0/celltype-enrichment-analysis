#!/bin/bash
# ============================================================
# MACS2 Peak Calling for Server
# Calls peaks by combining replicates for each cell type
# Uses exact parameters from original analysis
# ============================================================

cd /scratch/prj/bcn_marzi_lab/ratlas/Bomin

GENOME="2.5e8"

echo "=========================================="
echo "Starting MACS2 peak calling"
echo "Genome size: ${GENOME}"
echo "=========================================="

# Create output directories
mkdir -p macs2_output/NeuN
mkdir -p macs2_output/Nurr
mkdir -p macs2_output/Olig
mkdir -p macs2_output/Neg

# ============================================================
# NeuN - Neuronal marker
# ============================================================
echo ""
echo "[$(date)] Calling peaks for: NeuN"
echo "Replicates: IGF131357, IGF131358, IGF131359"
echo "Control: IGF131373"

macs2 callpeak \
  -t alignment_copy/IGF131357.bam \
     alignment_copy/IGF131358.bam \
     alignment_copy/IGF131359.bam \
  -c alignment_copy/IGF131373_input.bam \
  -f BAM -g ${GENOME} -n NeuN_all \
  --outdir macs2_output/NeuN \
  --nomodel --shift 0 --extsize 200 -B --SPMR

echo "[$(date)] NeuN peak calling complete"
wc -l macs2_output/NeuN/NeuN_all_peaks.narrowPeak

# ============================================================
# Nurr - Dopaminergic neurons
# ============================================================
echo ""
echo "[$(date)] Calling peaks for: Nurr"
echo "Replicates: IGF131366, IGF131367"
echo "Control: IGF131376"

macs2 callpeak \
  -t alignment_copy/IGF131366.bam \
     alignment_copy/IGF131367.bam \
  -c alignment_copy/IGF131376_input.bam \
  -f BAM -g ${GENOME} -n Nurr_all \
  --outdir macs2_output/Nurr \
  --nomodel --shift 0 --extsize 200 -B --SPMR

echo "[$(date)] Nurr peak calling complete"
wc -l macs2_output/Nurr/Nurr_all_peaks.narrowPeak

# ============================================================
# Olig - Oligodendrocytes
# ============================================================
echo ""
echo "[$(date)] Calling peaks for: Olig"
echo "Replicates: IGF131360, IGF131361, IGF131362"
echo "Control: IGF131374"

macs2 callpeak \
  -t alignment_copy/IGF131360.bam \
     alignment_copy/IGF131361.bam \
     alignment_copy/IGF131362.bam \
  -c alignment_copy/IGF131374_input.bam \
  -f BAM -g ${GENOME} -n Olig_all \
  --outdir macs2_output/Olig \
  --nomodel --shift 0 --extsize 200 -B --SPMR

echo "[$(date)] Olig peak calling complete"
wc -l macs2_output/Olig/Olig_all_peaks.narrowPeak

# ============================================================
# Neg - Negative control
# ============================================================
echo ""
echo "[$(date)] Calling peaks for: Neg"
echo "Replicates: IGF131369, IGF131370, IGF131371"
echo "Control: IGF131377"

macs2 callpeak \
  -t alignment_copy/IGF131369.bam \
     alignment_copy/IGF131370.bam \
     alignment_copy/IGF131371.bam \
  -c alignment_copy/IGF131377_input.bam \
  -f BAM -g ${GENOME} -n Neg_all \
  --outdir macs2_output/Neg \
  --nomodel --shift 0 --extsize 200 -B --SPMR

echo "[$(date)] Neg peak calling complete"
wc -l macs2_output/Neg/Neg_all_peaks.narrowPeak

# ============================================================
# Summary
# ============================================================
echo ""
echo "=========================================="
echo "Peak Calling Complete"
echo "=========================================="
echo "Peak counts:"
echo "  NeuN: $(wc -l < macs2_output/NeuN/NeuN_all_peaks.narrowPeak)"
echo "  Nurr: $(wc -l < macs2_output/Nurr/Nurr_all_peaks.narrowPeak)"
echo "  Olig: $(wc -l < macs2_output/Olig/Olig_all_peaks.narrowPeak)"
echo "  Neg:  $(wc -l < macs2_output/Neg/Neg_all_peaks.narrowPeak)"
echo "=========================================="
