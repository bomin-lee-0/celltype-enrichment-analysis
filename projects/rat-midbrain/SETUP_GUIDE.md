# Rat Midbrain ChIP-seq Analysis - Local Setup Guide

로컬 환경에서 ChIP-seq 분석 파이프라인을 실행하기 위한 가이드입니다.

## 📋 Table of Contents
1. [환경 설정](#환경-설정)
2. [데이터 준비](#데이터-준비)
3. [파이프라인 실행](#파이프라인-실행)
4. [주요 스크립트 설명](#주요-스크립트-설명)
5. [문제 해결](#문제-해결)

---

## 🔧 환경 설정

### 1. Conda 환경 생성

```bash
# 프로젝트 디렉토리로 이동
cd projects/rat-midbrain

# Conda 환경 생성 (처음 한 번만 실행)
conda env create -f environment.yml

# 환경 활성화
conda activate ratlas_env
```

### 2. 참조 게놈 다운로드 (필요시)

#### Rat genome (rn7) 다운로드:
```bash
# 게놈 FASTA 다운로드
wget https://hgdownload.soe.ucsc.edu/goldenPath/rn7/bigZips/rn7.fa.gz
gunzip rn7.fa.gz

# Bowtie2 인덱스 생성 (alignment 필요시)
bowtie2-build rn7.fa rn7_index/rn7
```

#### config.yaml 수정:
```yaml
genome:
  fasta: "/path/to/rn7.fa"
  bowtie2_index: "/path/to/rn7_index/rn7"
```

---

## 📁 데이터 준비

### 1. FASTQ 파일 배치

FASTQ 파일을 `0_data/` 폴더에 배치합니다:

```bash
projects/rat-midbrain/0_data/
├── IGF131357_R1.fastq.gz
├── IGF131357_R2.fastq.gz
├── IGF131358_R1.fastq.gz
├── IGF131358_R2.fastq.gz
├── ...
└── IGF131377_R2.fastq.gz
```

### 2. 샘플 정보

**Cell type별 샘플 구성:**
- **NeuN** (뉴런): IGF131357, IGF131358, IGF131359
  - Input: IGF131373
- **Nurr** (도파민성 뉴런): IGF131366, IGF131367
  - Input: IGF131376
- **Olig** (올리고덴드로사이트): IGF131360, IGF131361, IGF131362
  - Input: IGF131374
- **Neg** (음성 대조군): IGF131369, IGF131370, IGF131371
  - Input: IGF131377

---

## 🚀 파이프라인 실행

### 방법 1: 개별 스크립트 실행 (추천)

각 단계를 순차적으로 실행:

```bash
# 환경 활성화
conda activate ratlas_env

# Step 1: Preprocessing (FASTQ 파일 필요)
bash scripts/01_preprocessing_qc_local.sh

# Step 2: Alignment (참조 게놈 필요)
bash scripts/02_alignment.sh

# Step 3: Peak Calling (BAM 파일 필요)
bash scripts/03_peak_calling_local.sh

# Step 4: QC Metrics (FRiP 계산)
bash scripts/04_qc_metrics_local.sh

# Step 5: Peak Processing
bash scripts/05_peak_processing.sh

# Step 6: Annotation (R 필요)
Rscript scripts/06_annotation.R

# Step 7: Enrichment Analysis (R 필요)
Rscript scripts/07_enrichment.R

# Step 8: Motif Analysis
bash scripts/08_motif_analysis.sh
```

### 방법 2: Python 메인 파이프라인 실행

```bash
# 전체 파이프라인 실행
python main.py

# 특정 단계부터 실행
python main.py --start-from peak_calling

# 특정 단계까지 실행
python main.py --stop-at qc_metrics

# 특정 단계 건너뛰기
python main.py --skip preprocessing alignment

# 사용 가능한 단계 확인
python main.py --list-steps
```

---

## 📝 주요 스크립트 설명

### 로컬 버전 스크립트 (SLURM 불필요)

#### `01_preprocessing_qc_local.sh`
- **기능**: FASTQ QC 및 어댑터 트리밍
- **입력**: `0_data/*_R1.fastq.gz`, `*_R2.fastq.gz`
- **출력**: `01_preprocessing_qc/trimmed_fastq/`
- **원본과의 차이**:
  - SLURM array job → for loop
  - `module load` 제거 (conda 사용)
  - 로컬 경로 사용

#### `03_peak_calling_local.sh`
- **기능**: Cell type별 peak calling (replicates 통합)
- **입력**: `02_alignment/filtered_bam/*.bam`
- **출력**: `03_peak_calling/macs2_output/`
- **원본과의 차이**:
  - 경로만 로컬로 수정
  - MACS2 파라미터 동일 유지

#### `04_qc_metrics_local.sh`
- **기능**: FRiP score 계산 (fragment-level)
- **입력**: BAM 파일 + Peak 파일
- **출력**: `04_qc_metrics/frip_scores.csv`
- **원본과의 차이**:
  - 경로만 로컬로 수정
  - 로직 완전 동일

---

## ⚠️ 문제 해결

### 1. "command not found" 에러

```bash
# Conda 환경이 활성화되었는지 확인
conda activate ratlas_env

# 특정 도구 설치 확인
which trim_galore
which macs2
which bedtools
```

### 2. 메모리 부족

```bash
# MACS2에서 메모리 에러 발생시, 샘플을 나눠서 실행
# 또는 --buffer-size 옵션 조정
```

### 3. Bowtie2 인덱스 에러

```bash
# BOWTIE2_INDEX 환경 변수 설정
export BOWTIE2_INDEX=/path/to/rn7_index/rn7

# 또는 스크립트 내에서 직접 수정
# scripts/02_alignment.sh 파일의 BOWTIE2_INDEX 변수 수정
```

### 4. R 패키지 설치 문제

```R
# R 콘솔에서 수동 설치
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "ChIPseeker",
    "clusterProfiler",
    "TxDb.Rnorvegicus.UCSC.rn7.refGene",
    "org.Rn.eg.db"
))
```

### 5. HOMER genome 설치

```bash
# HOMER genome 설치 (motif analysis 필요시)
perl $(which configureHomer.pl) -install rn7

# 또는 최신 버전 확인
homer2 listGenomes
```

---

## 📊 결과 확인

### 주요 출력 파일:

```
01_preprocessing_qc/
├── multiqc_report.html          ← QC 전체 요약
└── trimmed_fastq/               ← 트리밍된 FASTQ

03_peak_calling/
└── macs2_output/
    ├── NeuN/NeuN_all_peaks.narrowPeak
    ├── Nurr/Nurr_all_peaks.narrowPeak
    ├── Olig/Olig_all_peaks.narrowPeak
    └── Neg/Neg_all_peaks.narrowPeak

04_qc_metrics/
├── frip_scores.csv              ← FRiP score 결과
└── qc_summary.txt               ← QC 요약

06_annotation/
├── annotation_summary.csv       ← Peak annotation
└── plots/                       ← 시각화

07_enrichment/
├── results/                     ← GO/KEGG 결과
└── figures/                     ← Enrichment plots

08_motif_analysis/
└── motif_results/               ← HOMER motif 결과
    └── */homerResults.html      ← 결과 HTML
```

---

## 🎯 빠른 시작 (BAM 파일이 이미 있는 경우)

이미 alignment까지 완료된 BAM 파일이 있다면:

```bash
# 1. BAM 파일을 올바른 위치에 배치
cp /path/to/bams/*.bam projects/rat-midbrain/02_alignment/filtered_bam/

# 2. Peak calling부터 시작
conda activate ratlas_env
bash scripts/03_peak_calling_local.sh
bash scripts/04_qc_metrics_local.sh
bash scripts/05_peak_processing.sh
Rscript scripts/06_annotation.R
Rscript scripts/07_enrichment.R
bash scripts/08_motif_analysis.sh
```

---

## 📧 문의

문제가 발생하거나 질문이 있으면:
- GitHub Issues 생성
- Email: bomin.lee@kcl.ac.uk

---

## 📚 참고 자료

- [MACS2 Documentation](https://github.com/macs3-project/MACS)
- [ChIPseeker Documentation](https://bioconductor.org/packages/ChIPseeker)
- [HOMER Documentation](http://homer.ucsd.edu/homer/)
- [Trim Galore Documentation](https://github.com/FelixKrueger/TrimGalore)
