#!/usr/bin/env python3
"""
세포타입별 enhancer intersect SNP들의 맨하탄 플롯
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import gzip
from pathlib import Path
import os

# 프로젝트 루트 디렉토리 설정
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
DATA_DIR = PROJECT_ROOT / "0.Data"
RESULTS_DIR = PROJECT_ROOT / "2.Results"

def load_celltype_annotations():
    """세포타입별 annotation 데이터 로드"""
    
    print("🧬 세포타입별 annotation 데이터 로딩 중...")
    
    celltypes = {
        'Microglia': 'Neg_cleaned',
        'Neuron': 'NeuN_cleaned', 
        'Oligodendrocyte': 'Olig_cleaned',
        'Dopaminergic': 'Nurr_cleaned'
    }
    
    all_annotations = {}
    
    for celltype, file_prefix in celltypes.items():
        print(f"  {celltype} 로딩 중...")
        
        # 염색체 1번만 먼저 테스트
        # LDSC annotation 파일이 없으면 BED 파일로부터 영역 정보 로드
        annot_file = DATA_DIR / "Results" / "annotations" / f"{file_prefix}.1.annot.gz"
        bed_file = DATA_DIR / "processed" / "hg19_coordinates" / f"cleaned_data_{file_prefix}_hg19.bed"

        try:
            if annot_file.exists():
                # LDSC annotation 파일 사용
                df = pd.read_csv(annot_file, sep='\t', compression='gzip')

                # 마지막 열이 enhancer annotation
                enhancer_col = df.columns[-1]  # 마지막 열

                # Intersect된 SNP들만 선택
                intersect_snps = df[df[enhancer_col] == 1]['SNP'].tolist()

                print(f"    염색체 1: {len(intersect_snps)}개 SNP가 {celltype} enhancer와 intersect")

                all_annotations[celltype] = {
                    'intersect_snps': set(intersect_snps),
                    'total_snps': len(df),
                    'intersect_count': len(intersect_snps)
                }
            elif bed_file.exists():
                # BED 파일에서 enhancer 영역 로드
                bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end', 'name'])

                print(f"    BED 파일에서 {len(bed_df)}개 enhancer 영역 로드")
                print(f"    주의: SNP-enhancer overlap은 GWAS 데이터 로딩 후 계산됩니다")

                all_annotations[celltype] = {
                    'bed_regions': bed_df,
                    'intersect_snps': set(),  # 나중에 계산
                    'total_snps': 0,
                    'intersect_count': 0
                }
            else:
                print(f"    경고: annotation 파일과 BED 파일을 찾을 수 없습니다")
                all_annotations[celltype] = {
                    'intersect_snps': set(),
                    'total_snps': 0,
                    'intersect_count': 0
                }

        except Exception as e:
            print(f"    오류: {e}")
            all_annotations[celltype] = {
                'intersect_snps': set(),
                'total_snps': 0,
                'intersect_count': 0
            }
    
    return all_annotations

def compute_snp_overlaps(gwas_df, annotations):
    """BED 영역과 GWAS SNP의 overlap 계산"""

    print("\n🔍 SNP-enhancer overlap 계산 중...")

    for celltype, data in annotations.items():
        if 'bed_regions' in data and data['bed_regions'] is not None:
            print(f"  {celltype} overlap 계산 중...")

            bed_df = data['bed_regions']
            intersect_snps = []

            # 각 SNP에 대해 enhancer 영역과 overlap 확인
            for _, snp in gwas_df.iterrows():
                snp_chr = snp['CHR']
                snp_pos = snp['BP']

                # 해당 SNP이 어떤 enhancer 영역에 포함되는지 확인
                overlaps = bed_df[
                    (bed_df['chr'].astype(str) == str(snp_chr)) &
                    (bed_df['start'] <= snp_pos) &
                    (bed_df['end'] >= snp_pos)
                ]

                if len(overlaps) > 0:
                    intersect_snps.append(snp['SNP'])

            data['intersect_snps'] = set(intersect_snps)
            data['intersect_count'] = len(intersect_snps)

            print(f"    {len(intersect_snps)}개 SNP가 {celltype} enhancer와 overlap")

    return annotations

def load_gwas_with_positions():
    """GWAS 데이터와 위치 정보 함께 로드"""
    
    print("📊 GWAS 데이터 로딩 중...")

    # GWAS 데이터
    gwas_file = DATA_DIR / "GWAS" / "GCST009325.h.tsv.gz"
    
    # 샘플링해서 로드 (메모리 절약)
    sample_size = 100000  # 테스트용: 10만개로 축소
    
    print(f"  GWAS 데이터 샘플링 로딩 중 (n={sample_size:,})...")

    df = pd.read_csv(gwas_file, sep='\t', compression='gzip', nrows=sample_size)

    # 컬럼명 확인 및 표준화
    print(f"  컬럼: {', '.join(df.columns[:10])}...")

    # P-value 컬럼 찾기
    p_col = None
    for col in ['p_value', 'P', 'pvalue', 'PVAL']:
        if col in df.columns:
            p_col = col
            break

    if p_col:
        df['P'] = df[p_col]
    elif 'Z' in df.columns:
        # Z-score가 있으면 P-value 계산
        df['P'] = 2 * (1 - norm.cdf(np.abs(df['Z'])))
    else:
        print("  경고: P-value 컬럼을 찾을 수 없습니다!")
        df['P'] = 0.5

    df['-log10P'] = -np.log10(np.maximum(df['P'], 1e-50))

    # 위치 정보 추가
    chr_col = next((c for c in ['chromosome', 'CHR', 'chr'] if c in df.columns), None)
    bp_col = next((c for c in ['base_pair_location', 'BP', 'bp', 'POS'] if c in df.columns), None)
    snp_col = next((c for c in ['variant_id', 'rsid', 'SNP', 'ID'] if c in df.columns), None)

    if chr_col:
        df['CHR'] = df[chr_col]
    else:
        df['CHR'] = 1

    if bp_col:
        df['BP'] = df[bp_col]
    else:
        df['BP'] = range(len(df))

    if snp_col:
        df['SNP'] = df[snp_col]
    elif 'rsid' in df.columns:
        df['SNP'] = df['rsid']
    else:
        df['SNP'] = [f"SNP_{i}" for i in range(len(df))]
    
    print(f"  로딩 완료: {len(df):,} SNPs")
    print(f"  P-value 범위: {df['P'].min():.2e} - {df['P'].max():.2e}")
    
    return df

def create_celltype_manhattan_plots(gwas_df, annotations):
    """세포타입별 맨하탄 플롯 생성"""
    
    print("🗽 세포타입별 맨하탄 플롯 생성 중...")
    
    # 4개 세포타입 서브플롯
    fig, axes = plt.subplots(2, 2, figsize=(20, 16))
    axes = axes.flatten()
    
    celltypes = ['Microglia', 'Neuron', 'Oligodendrocyte', 'Dopaminergic']
    colors_intersect = ['red', 'darkgreen', 'purple', 'orange']
    colors_background = ['lightcoral', 'lightgreen', 'plum', 'moccasin']
    
    # 유의성 기준선
    genome_wide = -np.log10(5e-8)
    suggestive = -np.log10(1e-5)
    
    for i, celltype in enumerate(celltypes):
        ax = axes[i]
        
        print(f"  {celltype} 플롯 생성 중...")
        
        # 해당 세포타입의 intersect SNP들
        intersect_snps = annotations[celltype]['intersect_snps']
        intersect_count = len(intersect_snps)
        
        # GWAS 데이터에서 intersect 여부 표시
        gwas_df['is_intersect'] = gwas_df['SNP'].isin(intersect_snps)
        
        # Non-intersect SNP들 (배경)
        non_intersect = gwas_df[~gwas_df['is_intersect']]
        intersect_data = gwas_df[gwas_df['is_intersect']]
        
        print(f"    Non-intersect SNPs: {len(non_intersect):,}")
        print(f"    Intersect SNPs: {len(intersect_data):,}")
        
        # 배경 SNP들 (작게, 투명하게)
        if len(non_intersect) > 0:
            # 너무 많으면 샘플링
            if len(non_intersect) > 50000:
                non_intersect_sample = non_intersect.sample(n=50000, random_state=42)
            else:
                non_intersect_sample = non_intersect
                
            ax.scatter(non_intersect_sample['BP'], non_intersect_sample['-log10P'], 
                      c=colors_background[i], alpha=0.3, s=1, label='Background SNPs')
        
        # Intersect SNP들 (크게, 진하게)
        if len(intersect_data) > 0:
            ax.scatter(intersect_data['BP'], intersect_data['-log10P'], 
                      c=colors_intersect[i], alpha=0.8, s=20, 
                      label=f'{celltype} Enhancer SNPs (n={len(intersect_data)})')
        
        # 유의성 기준선
        ax.axhline(y=genome_wide, color='red', linestyle='--', alpha=0.7, 
                  linewidth=1.5, label='p=5×10⁻⁸')
        ax.axhline(y=suggestive, color='blue', linestyle='--', alpha=0.5, 
                  linewidth=1, label='p=1×10⁻⁵')
        
        # 상위 intersect SNP들 라벨링
        if len(intersect_data) > 0:
            top_intersect = intersect_data.nlargest(3, '-log10P')
            for _, snp in top_intersect.iterrows():
                if snp['-log10P'] > suggestive:
                    ax.annotate(f"{snp['SNP']}", 
                               xy=(snp['BP'], snp['-log10P']),
                               xytext=(5, 5), textcoords='offset points',
                               fontsize=8, alpha=0.8,
                               bbox=dict(boxstyle='round,pad=0.2', 
                                       facecolor=colors_intersect[i], alpha=0.7))
        
        # 축 설정
        ax.set_xlabel('SNP Position (Chr 1)', fontsize=12, fontweight='bold')
        ax.set_ylabel('-log₁₀(P-value)', fontsize=12, fontweight='bold')
        ax.set_title(f'{celltype} Enhancer Manhattan Plot\\n'
                    f'({intersect_count} intersect SNPs)', 
                    fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right', fontsize=10)
        
        # Y축 범위 설정
        max_y = max(gwas_df['-log10P'].max(), genome_wide + 2)
        ax.set_ylim(0, max_y)
    
    plt.tight_layout()

    # 저장
    plots_dir = RESULTS_DIR / "Plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    plt.savefig(plots_dir / 'celltype_manhattan_plots.png',
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(plots_dir / 'celltype_manhattan_plots.pdf',
                bbox_inches='tight', facecolor='white')

    print(f"  저장 완료: {plots_dir / 'celltype_manhattan_plots.png'}")
    
    plt.show()
    
    # 통계 요약
    print(f"\\n📈 세포타입별 맨하탄 플롯 통계:")
    print("="*60)
    
    for celltype in celltypes:
        intersect_snps = annotations[celltype]['intersect_snps']
        intersect_data = gwas_df[gwas_df['SNP'].isin(intersect_snps)]
        
        if len(intersect_data) > 0:
            significant_intersect = intersect_data[intersect_data['-log10P'] > suggestive]
            genome_wide_intersect = intersect_data[intersect_data['-log10P'] > genome_wide]
            
            print(f"\\n🧬 {celltype}:")
            print(f"  Total intersect SNPs: {len(intersect_data)}")
            print(f"  Suggestive (p<1e-5): {len(significant_intersect)}")
            print(f"  Genome-wide (p<5e-8): {len(genome_wide_intersect)}")
            
            if len(intersect_data) > 0:
                print(f"  Max -log10(P): {intersect_data['-log10P'].max():.2f}")
                
                # 상위 3개 SNP
                top_snps = intersect_data.nlargest(3, '-log10P')
                print(f"  Top SNPs:")
                for _, snp in top_snps.iterrows():
                    print(f"    {snp['SNP']}: p={snp['P']:.2e}")
        else:
            print(f"\\n🧬 {celltype}: No intersect SNPs found")

def create_comparison_manhattan(gwas_df, annotations):
    """4개 세포타입 비교 맨하탄 플롯"""
    
    print("\\n🔄 세포타입 비교 맨하탄 플롯 생성 중...")
    
    # 단일 플롯에서 4개 세포타입 비교
    fig, ax = plt.subplots(figsize=(16, 10))
    
    celltypes = ['Microglia', 'Neuron', 'Oligodendrocyte', 'Dopaminergic']
    colors = ['red', 'green', 'blue', 'orange']
    
    # 각 세포타입별로 점 표시
    for i, celltype in enumerate(celltypes):
        intersect_snps = annotations[celltype]['intersect_snps']
        intersect_data = gwas_df[gwas_df['SNP'].isin(intersect_snps)]
        
        if len(intersect_data) > 0:
            ax.scatter(intersect_data['BP'], intersect_data['-log10P'], 
                      c=colors[i], alpha=0.7, s=15, label=f'{celltype} (n={len(intersect_data)})')
    
    # 유의성 기준선
    genome_wide = -np.log10(5e-8)
    suggestive = -np.log10(1e-5)
    
    ax.axhline(y=genome_wide, color='red', linestyle='--', alpha=0.8, 
              linewidth=2, label='Genome-wide (p=5×10⁻⁸)')
    ax.axhline(y=suggestive, color='blue', linestyle='--', alpha=0.6, 
              linewidth=1.5, label='Suggestive (p=1×10⁻⁵)')
    
    # 축 설정
    ax.set_xlabel('SNP Position (Chr 1)', fontsize=14, fontweight='bold')
    ax.set_ylabel('-log₁₀(P-value)', fontsize=14, fontweight='bold')
    ax.set_title('Cell Type-Specific Enhancer Manhattan Plot Comparison', 
                 fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()

    # 저장
    plots_dir = RESULTS_DIR / "Plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    plt.savefig(plots_dir / 'celltype_comparison_manhattan.png',
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(plots_dir / 'celltype_comparison_manhattan.pdf',
                bbox_inches='tight', facecolor='white')

    print(f"  저장 완료: {plots_dir / 'celltype_comparison_manhattan.png'}")

    plt.show()

if __name__ == "__main__":
    print("🗽 세포타입별 enhancer 맨하탄 플롯 생성!")
    print("="*60)

    # 1. 세포타입별 annotation 로드
    annotations = load_celltype_annotations()

    # 2. GWAS 데이터 로드
    gwas_data = load_gwas_with_positions()

    # 3. BED 파일 사용 시 SNP-enhancer overlap 계산
    annotations = compute_snp_overlaps(gwas_data, annotations)

    # 4. 세포타입별 맨하탄 플롯 생성
    create_celltype_manhattan_plots(gwas_data, annotations)

    # 5. 비교 맨하탄 플롯 생성
    create_comparison_manhattan(gwas_data, annotations)

    print(f"\n✅ 세포타입별 맨하탄 플롯 생성 완료!")
    print(f"   저장된 파일:")
    plots_dir = RESULTS_DIR / "Plots"
    print(f"   - {plots_dir / 'celltype_manhattan_plots.png'} (4-panel)")
    print(f"   - {plots_dir / 'celltype_comparison_manhattan.png'} (overlay)")