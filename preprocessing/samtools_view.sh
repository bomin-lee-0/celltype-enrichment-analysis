#!/bin/bash

#SBATCH -p cpu
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --output=/scratch/prj/bcn_marzi_lab/ratlas/Mingshen/alignment/samtools_view.out
#SBATCH --array=1-17
#SBATCH --time=24:00:00
#SBATCH --job-name=samtools_view

module load samtools

DATA_DIR=/scratch/prj/bcn_marzi_lab/ratlas/Mingshen
BAM_DIR=/scratch/prj/bcn_marzi_lab/ratlas/Mingshen/alignment/BAM


cd $DATA_DIR

file=$(head -n $SLURM_ARRAY_TASK_ID sample_17.txt | tail -1)

samtools view -b -F 4 -q 30 ${BAM_DIR}/sorted/$file.dedup.bam > ${BAM_DIR}/final/$file.dedup.q30.bam
