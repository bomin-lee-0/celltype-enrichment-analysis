# Cell-type-specific H3K27ac enhancer landscape in the rat midbrain and its genetic enrichment for Parkinson’s disease

This repository contains the full analysis pipeline used to identify cell-type-specific enhancer landscapes in the rat midbrain and evaluate their genetic enrichment for Parkinson’s disease (PD) GWAS signals.  
All analyses were conducted using **H3K27ac ChIP-seq data** derived from purified nuclei of midbrain cell types, and performed on the King’s College London **CREATE HPC cluster**.

---

## 🧠 Overview

Parkinson’s disease is a neurodegenerative disorder caused by complex interactions between genetic and environmental factors.  
Most PD-associated genetic variants lie in **non-coding regions**, suggesting that they may affect gene regulation via enhancers active in specific cell types.  
This study aims to determine **which midbrain cell types** show the strongest genetic enrichment for PD risk variants by integrating epigenomic and statistical genetics analyses.

---

## ⚙️ Workflow Summary

| Step | Description | Tools |
|------|--------------|-------|
| **1. Pre-QC & Trimming** | FASTQ QC → adapter/quality trimming → post-trim QC | `FastQC`, `Trim Galore`/`cutadapt`, `MultiQC` |
| **2. Alignment** | Map trimmed reads to rat genome (rn7), sort, mark/remove duplicates, filter | `Bowtie2` |
| **3. Peak Calling** | Identification of H3K27ac peaks per sample | `MACS2` |
| **4. QC Metrics** | FRiP (Fraction of Reads in Peaks) on fragment-level | `bedtools` |
| **5. Peak Processing** | Merge peaks and extract cell-type–specific unique peaks | `bedtools merge`, `bedtools intersect` |
| **6. Annotation** | Genomic/functional annotation of peaks | `ChIPseeker`, `TxDb` |
| **7. Functional Enrichment** | GO and KEGG enrichment analyses | `clusterProfiler` |
| **8. Motif Analysis** | De novo motif discovery from unique peaks | `HOMER` |
---

## 📁 Repository Structure

```
rat-midbrain-celltype-enrichment/
│
├── 01_preprocessing_qc/
│   ├── fastqc_raw/
│   ├── fastqc_trimmed/
│   ├── multiqc_report.html
│   └── trimming_logs/
│
├── 02_alignment/
│   ├── bowtie2_alignment_commands.txt
│   ├── aligned_bam/
│   ├── dedup_bam/
│   ├── filtered_bam/
│   └── alignment_summary.txt
│
├── 03_peak_calling/
│   ├── macs2_output/
│   ├── macs2_commands.txt
│
├── 04_qc_metrics/
│   └── frip_scores.csv
│
├── 05_peak_processing/
│   ├── merged_peaks/
│   ├── unique_peaks/
│   └── bedtools_commands.txt
│
├── 06_annotation/
│   ├── annotate_peak.R
│   ├── annotation_summary.csv
│   ├── plots/
│
├── 07_enrichment/
│   ├── clusterProfiler_GO_KEGG.R
│   ├── results/
│   ├── figures/
│
└── 08_motif_analysis/
    ├── homer_commands.txt
    ├── motif_results/ 
```

---

## 🧩 Example Commands

### 🔹 1) Pre-QC, Trimming, Post-QC
```bash
# Raw QC
fastqc -t 8 raw/sample_R1.fastq.gz raw/sample_R2.fastq.gz -o 01_preprocessing_qc/fastqc_raw/

# Trimming (choose one: Trim Galore or cutadapt)


# Post-trim QC
fastqc -t 8 01_preprocessing_qc/sample_R1_val_1.fq.gz 01_preprocessing_qc/sample_R2_val_2.fq.gz -o 01_preprocessing_qc/fastqc_trimmed/

# MultiQC summary
multiqc 01_preprocessing_qc/ -o 01_preprocessing_qc/
```

---

### 🔹 2) Alignment 
```bash
```

---

### 🔹 3) Peak Calling (MACS2)
```bash
macs2 callpeak \
  -t 02_alignment/filtered_bam/sample.filtered.bam \
  -c 02_alignment/filtered_bam/input.filtered.bam \
  -f BAM -g 2.5e8 -n sample_name \
  --outdir 03_peak_calling/macs2_output
```

---

### 🔹 4) FRiP Score (Fragment-level)
```bash
# Example: fragments within peaks (pair-to-fragment implied by filtering above)
bedtools intersect -a 02_alignment/filtered_bam/sample.filtered.bam -b 03_peak_calling/macs2_output/sample_name_peaks.narrowPeak -bed \
| wc -l
```

---

### 🔹 5) Peak Annotation (R / ChIPseeker)
```r
library(ChIPseeker)
library(TxDb.Rnorvegicus.UCSC.rn7.refGene)

txdb <- TxDb.Rnorvegicus.UCSC.rn7.refGene
peakAnno <- annotatePeak("05_peak_processing/unique_peaks/Nurr_unique.bed",
                         TxDb = txdb, tssRegion = c(-3000, 3000))
plotAnnoPie(peakAnno)
```

---

### 🔹 6) GO / KEGG Enrichment (clusterProfiler)
```r
library(clusterProfiler)
ego <- enrichGO(gene = geneList,
                OrgDb = org.Rn.eg.db,
                keyType = "ENTREZID",
                ont = "BP")
dotplot(ego, showCategory = 20)
```

---

### 🔹 7) Motif Discovery (HOMER)
```bash
findMotifsGenome.pl 05_peak_processing/unique_peaks/Nurr_unique.bed rn7 08_motif_analysis/motif_results/ -size given
```

---



## 🧬 Key Results

- 🧠 **Dopaminergic (Nurr1⁺) neurons** exhibited the highest PD genetic enrichment among midbrain cell types.  
- **Oligodendrocytes** showed moderate enrichment, consistent with glial involvement in PD.  
- FRiP and QC metrics confirmed high data quality (FRiP > 0.2 for all samples).  
- GO/KEGG analyses revealed pathways related to **synaptic signaling, axonogenesis, and immune response**.  
- Motif analysis identified enrichment of **Nurr1**, **FOXA2**, and **PU.1** motifs.  

---

## ⚙️ Computational Environment

```bash
# Conda environment
conda create -n ratlas_env python=3.9 -y
conda activate ratlas_env
conda install -y -c bioconda bowtie2 samtools macs2 bedtools homer
conda install -y -c conda-forge r-base=4.3 r-ggplot2
conda install -y -c bioconda bioconductor-chipseeker bioconductor-clusterprofiler
```

**Software versions**
- Python 3.9  
- R 4.3.1  
- Bowtie2 2.5.x  
- Samtools 1.18+  
- MACS2 2.2.9.1  
- bedtools 2.31.0  
- HOMER v4.11  
---

## 📚 References

1. Nott, A. *et al.* (2019). **Cell type–specific enhancer–promoter interactome maps of the human brain.** *Science*, 366(6469), 1134–1139.  
2. Agarwal, D. *et al.* (2020). **A single-cell atlas of the human midbrain identifies cell-specific genetic risk for Parkinson’s disease.** *Nature Neuroscience*, 23, 939–951.  
3. Bulik-Sullivan, B. *et al.* (2015). **An atlas of genetic correlations across human diseases and traits.** *Nature Genetics*, 47(11), 1236–1241.  
4. Alexei, A. *et al.* (2018). **Epigenomic annotation of human cortical cell types.** *Nature Neuroscience*, 21(1), 125–132.  

---

## ✨ Author

**Bomin Lee**  
MSc Psychology & Neuroscience of the Mind–Body Interface  
King’s College London  

For questions or collaboration: **bomin.lee@kcl.ac.uk**

---

## 🧩 License

This repository is distributed for academic and educational use only.  
All data and analyses are part of the **Marzi Lab** at King’s College London and may not be used for commercial purposes.
```
