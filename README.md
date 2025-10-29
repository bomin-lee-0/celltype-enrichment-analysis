# Cell-Type Enrichment Analysis

[![Research Focus](https://img.shields.io/badge/Research-Parkinson's%20Disease-blue)](https://github.com)
[![Analysis Type](https://img.shields.io/badge/Analysis-GWAS%20Enrichment-green)](https://github.com)
[![Cell Types](https://img.shields.io/badge/Cell%20Types-Multi--Type-orange)](https://github.com)

## Overview

This repository contains comprehensive cell-type-specific enrichment analyses for Parkinson's disease (PD). The analyses integrate genomic, epigenomic, and statistical genetics approaches to identify which brain cell types show the strongest genetic contribution to PD pathogenesis.

## Projects

This repository is organized into independent but related projects, each focusing on different aspects of cell-type-specific analysis:

### 1. Parkinson's Disease GWAS Enrichment Analysis (`projects/parkinsons/`)

A complete analysis pipeline for quantifying cell-type-specific genetic enrichment in Parkinson's disease using human GWAS data and enhancer regions.

**Key Features:**
- Analysis of 17.4M SNPs from large-scale PD GWAS (37,688 cases, 1.4M controls)
- Four brain cell types analyzed: Oligodendrocytes, Dopaminergic Neurons, General Neurons, Microglia
- LDSC (Linkage Disequilibrium Score Regression) partitioned heritability analysis
- Automated batch analysis system with visualization
- hg19/GRCh37 genome build

**Status:** Production-ready with full LDSC implementation

[View detailed documentation →](projects/parkinsons/README.md)

### 2. Rat Midbrain H3K27ac Enhancer Landscape (`projects/rat-midbrain/`)

A comprehensive ChIP-seq analysis pipeline for identifying cell-type-specific enhancer landscapes in the rat midbrain and evaluating their genetic enrichment for Parkinson's disease.

**Key Features:**
- H3K27ac ChIP-seq data from purified midbrain nuclei
- Complete bioinformatics pipeline: QC, alignment, peak calling, annotation
- Cell-type-specific unique peak identification
- Functional enrichment and motif analysis
- rn7 genome build

**Status:** Active development

[View detailed documentation →](projects/rat-midbrain/README.md)

## Repository Structure

```
celltype-enrichment-analysis/
├── projects/
│   ├── parkinsons/          # Human GWAS enrichment analysis
│   │   ├── 0.Data/          # GWAS and enhancer data
│   │   ├── 1.Scripts/       # Analysis scripts (LDSC, visualization)
│   │   ├── 2.Results/       # Analysis outputs
│   │   ├── 3.Documentation/ # Workflow documentation
│   │   └── main.py          # Main execution script
│   │
│   └── rat-midbrain/        # Rat midbrain ChIP-seq analysis
│       ├── 01_preprocessing_qc/
│       ├── 02_alignment/
│       ├── 03_peak_calling/
│       ├── 04_qc_metrics/
│       ├── 05_peak_processing/
│       ├── 06_annotation/
│       ├── 07_enrichment/
│       └── 08_motif_analysis/
│
├── README.md                # This file
└── LICENSE                  # License information
```

## Scientific Rationale

### Cell-Type-Specific Genetic Architecture

Parkinson's disease is a complex neurodegenerative disorder with strong genetic and environmental components. Most PD-associated genetic variants are located in non-coding regions, suggesting they affect gene regulation through enhancers active in specific cell types.

### Multi-Cell-Type Comparison

This research addresses the fundamental question: **Which brain cell types contribute most to PD genetic risk?**

| Cell Type | Function | PD Relevance |
|-----------|----------|--------------|
| **Dopaminergic Neurons** | Dopamine production, motor control | Primary lesion site, neuronal loss |
| **Oligodendrocytes** | Myelin formation, white-matter maintenance | White-matter damage, connectivity deficits |
| **Microglia** | Immune response, neuroprotection | Neuroinflammation, immune activation |
| **General Neurons** | Neural signal transmission | General neural network disruption |

## Quick Start

### Parkinson's GWAS Analysis

```bash
cd projects/parkinsons/

# Run complete LDSC analysis
python main.py --all

# Or run by stages
python main.py --step coordinate  # Coordinate transformation
python main.py --step ldsc        # LDSC analysis
python main.py --step visualize   # Visualization
```

### Rat Midbrain Analysis

```bash
cd projects/rat-midbrain/

# Follow the step-by-step workflow in the project README
# See projects/rat-midbrain/README.md for detailed commands
```

## Key Methodologies

### Statistical Genetics
- **LDSC Partitioned Heritability:** Quantifies genetic contribution per cell type
- **Mann-Whitney U Test:** Compares p-value distributions
- **Fisher's Exact Test:** Tests enrichment of significant SNPs
- **Enrichment Ratio:** Measures fold-enrichment in enhancer regions

### Epigenomics
- **H3K27ac ChIP-seq:** Identifies active enhancers
- **Peak Calling:** MACS2-based peak detection
- **Cell-Type-Specific Peaks:** Unique and shared enhancer identification
- **Functional Annotation:** GO/KEGG enrichment analysis

### Genomics
- **Coordinate Liftover:** Cross-species and build conversion (rn7 ↔ hg38 ↔ hg19)
- **Large-Scale Data Processing:** Efficient handling of 17M+ SNPs
- **Quality Control:** Rigorous QC at every step

## Key Results

### Human GWAS Findings
- Dopaminergic neurons show **highest genetic enrichment** (enrichment > 3.0)
- Oligodendrocytes show **significant enrichment** (enrichment > 2.0)
- Microglia show **moderate enrichment** (enrichment > 1.5)
- General neurons show **no significant enrichment**

### Rat Midbrain Findings
- Dopaminergic (Nurr1+) neurons exhibit **highest PD genetic enrichment**
- FRiP scores > 0.2 confirm high data quality
- Enrichment of Nurr1, FOXA2, and PU.1 transcription factor motifs
- GO pathways: synaptic signaling, axonogenesis, immune response

## Clinical and Therapeutic Implications

### Tiered Therapeutic Strategy
1. **Priority 1:** Dopaminergic neuron protection and regeneration
2. **Priority 2:** Oligodendrocyte support and myelin repair
3. **Priority 3:** Neuroinflammation modulation (microglia-targeted)

### Personalized Medicine
- Cell-type-specific risk profiling
- Therapy response prediction based on enrichment patterns
- Biomarker development for targeted cell types

## Technical Requirements

### Python Environment
```bash
conda create -n celltype_analysis python=3.9 -y
conda activate celltype_analysis
pip install numpy pandas scipy matplotlib seaborn
```

### R Environment (for ChIP-seq)
```bash
conda install -y -c conda-forge r-base=4.3
conda install -y -c bioconda bioconductor-chipseeker bioconductor-clusterprofiler
```

### External Tools
- LDSC (for partitioned heritability)
- UCSC liftOver (for coordinate conversion)
- Bowtie2, MACS2, bedtools, HOMER (for ChIP-seq)

## Data Sources

- **GWAS Data:** GWAS Catalog (GCST009325) - Nalls et al. (2019)
- **Reference Genomes:** hg19, hg38, rn7 (UCSC Genome Browser)
- **LDSC References:** 1000 Genomes, BaselineLD v2.2, HapMap3
- **ChIP-seq Data:** H3K27ac from purified midbrain nuclei

## References

### Primary Publications
1. Nalls, M.A., et al. (2019). Identification of novel risk loci for Parkinson's disease. *The Lancet Neurology*, 18(12), 1091–1102.
2. Finucane, H.K., et al. (2015). Partitioning heritability by functional annotation. *Nature Genetics*, 47(11), 1228–1235.
3. Nott, A. et al. (2019). Cell type-specific enhancer-promoter interactome maps. *Science*, 366(6469), 1134–1139.
4. Agarwal, D. et al. (2020). A single-cell atlas of the human midbrain. *Nature Neuroscience*, 23, 939–951.

### Resources
- GWAS Catalog: https://www.ebi.ac.uk/gwas/
- LDSC Software: https://github.com/bulik/ldsc
- 1000 Genomes: https://www.internationalgenome.org/
- UCSC Genome Browser: https://genome.ucsc.edu/

## Author

**Bomin Lee**
MSc Psychology & Neuroscience of the Mind-Body Interface
King's College London

For questions or collaboration: bomin.lee@kcl.ac.uk

## License

This repository is distributed for academic and educational use. All analyses are part of research conducted at King's College London and may not be used for commercial purposes without permission.

## Acknowledgments

- Marzi Lab, King's College London
- CREATE HPC Cluster
- GWAS Catalog and 1000 Genomes Project
- All contributors to the open-source bioinformatics tools used in this research

---

**Last Updated:** 2025-10-29
