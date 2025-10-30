# Server Scripts for HPC Cluster

이 폴더에는 King's College London CREATE HPC 클러스터에서 실행할 수 있는 서버용 스크립트가 포함되어 있습니다.

## 📁 파일 구성

```
server/
├── sample.txt                          # 샘플 리스트 (SLURM array job용)
├── 01_trimming_slurm.sh               # FASTQ trimming (SLURM array job)
├── 02_macs2_peak_calling.sh           # MACS2 peak calling
├── 03_frip_calculation.sh             # FRiP score calculation
├── 04_unique_peaks.sh                 # Cell type-specific unique peaks
├── 05_merged_peaks_non_intersect.sh   # All-cell merged peaks
└── 06_peak_annotation.R               # Peak annotation with ChIPseeker
```

---

## 🚀 실행 순서

### 0. 준비 사항

```bash
# 서버에 로그인
ssh username@create.hpc.kcl.ac.uk

# 작업 디렉토리로 이동
cd /scratch/prj/bcn_marzi_lab/ratlas/Bomin

# 필요한 디렉토리 생성
mkdir -p raw_reads
mkdir -p trimmed_output/fastqc
mkdir -p alignment_copy
mkdir -p macs2_output
mkdir -p merged_unique
mkdir -p annotation
```

### 1. Trimming (SLURM Array Job)

```bash
# sample.txt 파일이 작업 디렉토리에 있는지 확인
cat sample.txt

# SLURM array job 제출 (15개 샘플 병렬 처리)
sbatch 01_trimming_slurm.sh

# 작업 상태 확인
squeue -u $USER

# 로그 확인
tail -f trimmed_output/trim.*.out
```

**주의사항:**
- `sample.txt`에는 R1 파일명만 나열 (R2는 자동으로 매칭)
- `--array=1-15`: 샘플 개수에 맞게 조정
- FASTQ 파일은 `raw_reads/` 폴더에 있어야 함

---

### 2. Alignment (이미 완료된 경우 skip)

Alignment는 별도의 Bowtie2 스크립트로 수행.
결과 BAM 파일들을 `alignment_copy/` 폴더에 배치.

---

### 3. Peak Calling with MACS2

```bash
# 스크립트 실행
bash 02_macs2_peak_calling.sh

# 또는 SLURM job으로 제출
sbatch -p cpu -c 4 --mem=16G --time=02:00:00 02_macs2_peak_calling.sh

# 결과 확인
ls -lh macs2_output/*/
```

**출력:**
- `NeuN_all_peaks.narrowPeak` (3 replicates)
- `Nurr_all_peaks.narrowPeak` (2 replicates)
- `Olig_all_peaks.narrowPeak` (3 replicates)
- `Neg_all_peaks.narrowPeak` (3 replicates)

---

### 4. FRiP Score Calculation

```bash
# 스크립트 실행
bash 03_frip_calculation.sh

# 결과 확인
cat frip_scores.csv
```

**출력:**
- `frip_scores.csv`: Sample별 FRiP score

**FRiP Score 해석:**
- **PASS**: FRiP ≥ 0.2 (20%)
- **WARNING**: 0.1 ≤ FRiP < 0.2
- **FAIL**: FRiP < 0.1

---

### 5. Cell Type-Specific Unique Peaks

```bash
# 각 cell type의 고유한 peaks 찾기
bash 04_unique_peaks.sh

# 결과 확인
ls -lh merged_unique/unique_*.bed
```

**출력:**
- `unique_NeuN.bed`
- `unique_Nurr.bed`
- `unique_Olig.bed`
- `unique_Neg.bed`

---

### 6. Merged Peaks (Non-intersect)

```bash
# 모든 cell type의 peaks를 합침
bash 05_merged_peaks_non_intersect.sh

# 결과 확인
wc -l merged_unique/all_merged_peaks.bed
```

**출력:**
- `all_merged_peaks.bed`: 모든 peaks 통합

---

### 7. Peak Annotation

```bash
# R 모듈 로드 (필요시)
module load R/4.3

# 스크립트 실행
Rscript 06_peak_annotation.R

# 또는 SLURM job으로 제출
sbatch -p cpu -c 4 --mem=32G --time=01:00:00 --wrap="Rscript 06_peak_annotation.R"

# 결과 확인
ls -lh annotation/
```

**출력:**
- `peakAnnoList.rds`: 전체 annotation 객체
- `annotation_summary.csv`: Annotation 통계
- `*_annotated.csv`: Cell type별 상세 annotation
- `*_genes.txt`: Cell type별 gene list
- `*.pdf`: 시각화 plots

---

## 📊 예상 결과

### Peak Counts (예시)
```
Cell Type    Total Peaks    Unique Peaks    % Unique
---------    -----------    ------------    --------
NeuN         ~50,000        ~20,000         40%
Nurr         ~40,000        ~15,000         37%
Olig         ~45,000        ~18,000         40%
Neg          ~30,000        ~10,000         33%
```

### FRiP Scores (예시)
```
Sample         Cell Type    FRiP Score    Status
------         ---------    ----------    ------
IGF131357      NeuN         0.25          PASS
IGF131358      NeuN         0.23          PASS
IGF131359      NeuN         0.21          PASS
IGF131366      Nurr         0.28          PASS
IGF131367      Nurr         0.26          PASS
...
```

---

## 🔧 문제 해결

### SLURM Job이 실패할 때

```bash
# 로그 파일 확인
cat trimmed_output/trim.1.out

# 에러 메시지 검색
grep -i error trimmed_output/*.out

# 특정 작업 취소
scancel JOB_ID

# 모든 작업 취소
scancel -u $USER
```

### 모듈 로드 에러

```bash
# 사용 가능한 모듈 확인
module avail

# R 모듈 로드
module load R/4.3

# trim_galore 모듈
module load trim_galore
```

### 메모리 부족

```bash
# 더 많은 메모리로 재제출
sbatch -p cpu -c 8 --mem=64G 02_macs2_peak_calling.sh
```

---

## 📝 경로 수정 가이드

스크립트를 다른 프로젝트에 사용하려면:

1. **작업 디렉토리 변경:**
   ```bash
   # 모든 스크립트에서 다음 부분 수정
   cd /scratch/prj/YOUR_PROJECT/YOUR_NAME
   ```

2. **샘플 ID 변경:**
   - `sample.txt` 파일 수정
   - FRiP 스크립트의 샘플 매칭 로직 수정

3. **SLURM 출력 경로:**
   ```bash
   #SBATCH --output=/your/path/trim.%a.out
   ```

---

## 📧 문의

문제가 발생하면:
- Email: bomin.lee@kcl.ac.uk
- HPC support: hpc@kcl.ac.uk

---

## 🔗 참고

- [CREATE HPC Documentation](https://docs.er.kcl.ac.uk/)
- [SLURM Documentation](https://slurm.schedmd.com/)
- [MACS2 Documentation](https://github.com/macs3-project/MACS)
