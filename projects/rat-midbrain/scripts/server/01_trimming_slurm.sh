#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 8
#SBATCH --mem=10G
#SBATCH --output=/scratch/prj/bcn_marzi_lab/ratlas/Bomin/trimmed_output/trim.%a.out
#SBATCH --array=1-15
#SBATCH --time=00:30:00

# ============================================================
# SLURM Array Job for Trimming FASTQ files
# Processes multiple samples in parallel on HPC cluster
# ============================================================

module load trim_galore

# Change to working directory
cd /scratch/prj/bcn_marzi_lab/ratlas/Bomin

# Read sample name from sample.txt based on array task ID
sample=$(head -n $SLURM_ARRAY_TASK_ID sample.txt | tail -n 1)

# Define input file paths
INPUT_R1="/scratch/prj/bcn_marzi_lab/ratlas/raw_reads/${sample}.fastq.gz"
INPUT_R2="/scratch/prj/bcn_marzi_lab/ratlas/raw_reads/${sample/R1/R2}.fastq.gz"

echo "=========================================="
echo "Processing sample: ${sample}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "R1: ${INPUT_R1}"
echo "R2: ${INPUT_R2}"
echo "=========================================="

# Create output directories if they don't exist
mkdir -p trimmed_output/fastqc

# Run trim_galore with paired-end mode
trim_galore --paired --fastqc \
    --fastqc_args "--outdir trimmed_output/fastqc" \
    -o trimmed_output \
    "$INPUT_R1" "$INPUT_R2"

echo "=========================================="
echo "Sample ${sample} completed"
echo "=========================================="
