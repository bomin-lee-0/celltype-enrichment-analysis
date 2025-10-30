#!/usr/bin/env python3
"""
파킨슨병 GWAS - 세포타입별 Enhancer Enrichment 분석
Main execution script
"""

import os
import sys
import argparse
from pathlib import Path

# Add scripts directories to path
sys.path.append(str(Path(__file__).parent / "1.Scripts" / "LDSC"))
sys.path.append(str(Path(__file__).parent / "1.Scripts" / "Visualization"))
sys.path.append(str(Path(__file__).parent / "1.Scripts" / "Utils"))

def main():
    parser = argparse.ArgumentParser(
        description="파킨슨병 GWAS 세포타입별 Enhancer Enrichment 분석",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
사용 예시:
  # 1. 좌표 변환
  python main.py --step coordinate
  
  # 2. LDSC 분석
  python main.py --step ldsc
  
  # 3. 시각화
  python main.py --step visualize
  
  # 전체 파이프라인 실행
  python main.py --all
        """
    )
    
    parser.add_argument(
        '--step', 
        choices=['coordinate', 'ldsc', 'visualize'],
        help='실행할 단계 선택'
    )
    parser.add_argument(
        '--all', 
        action='store_true',
        help='전체 파이프라인 실행'
    )
    
    args = parser.parse_args()
    
    if args.all:
        print("=" * 60)
        print("파킨슨병 GWAS 세포타입별 Enhancer Enrichment 분석")
        print("전체 파이프라인 실행")
        print("=" * 60)

        # Step 1: Coordinate conversion
        print("\n[1/3] 좌표계 변환...")
        try:
            from setup_liftover import convert_all_enhancer_files
            convert_all_enhancer_files()
            print("✅ 좌표계 변환 완료")
        except Exception as e:
            print(f"⚠️  좌표계 변환 건너뜀 (이미 완료되었거나 오류 발생): {e}")

        # Step 2: LDSC analysis (optional - may fail if LDSC not configured)
        print("\n[2/3] LDSC 분석...")
        try:
            from ldsc_analysis_system import main as ldsc_main
            ldsc_main()
            print("✅ LDSC 분석 완료")
        except Exception as e:
            print(f"⚠️  LDSC 분석 건너뜀 (LDSC 설정 필요 또는 오류): {e}")
            print("   (시각화는 계속 진행됩니다)")

        # Step 3: Visualization
        print("\n[3/3] 시각화...")
        try:
            import sys
            sys.path.insert(0, str(Path(__file__).parent / "1.Scripts" / "Visualization"))
            import celltype_manhattan_plot
            # Run main code from the module
            celltype_manhattan_plot.annotations = celltype_manhattan_plot.load_celltype_annotations()
            celltype_manhattan_plot.gwas_data = celltype_manhattan_plot.load_gwas_with_positions()
            celltype_manhattan_plot.annotations = celltype_manhattan_plot.compute_snp_overlaps(
                celltype_manhattan_plot.gwas_data, celltype_manhattan_plot.annotations
            )
            celltype_manhattan_plot.create_celltype_manhattan_plots(
                celltype_manhattan_plot.gwas_data, celltype_manhattan_plot.annotations
            )
            celltype_manhattan_plot.create_comparison_manhattan(
                celltype_manhattan_plot.gwas_data, celltype_manhattan_plot.annotations
            )
            print("✅ 시각화 완료")
        except Exception as e:
            print(f"❌ 시각화 실패: {e}")
            import traceback
            traceback.print_exc()

        print("\n✅ 전체 파이프라인 완료!")
        
    elif args.step == 'coordinate':
        print("좌표계 변환 실행...")
        try:
            from setup_liftover import convert_all_enhancer_files
            convert_all_enhancer_files()
            print("✅ 좌표계 변환 완료")
        except Exception as e:
            print(f"❌ 좌표계 변환 실패: {e}")
            import traceback
            traceback.print_exc()

    elif args.step == 'ldsc':
        print("LDSC 분석 실행...")
        try:
            from ldsc_analysis_system import main as ldsc_main
            ldsc_main()
            print("✅ LDSC 분석 완료")
        except Exception as e:
            print(f"❌ LDSC 분석 실패: {e}")
            print("LDSC 환경이 설정되어 있는지 확인하세요.")
            import traceback
            traceback.print_exc()

    elif args.step == 'visualize':
        print("시각화 실행...")
        try:
            import sys
            sys.path.insert(0, str(Path(__file__).parent / "1.Scripts" / "Visualization"))
            import celltype_manhattan_plot
            # Run main code
            annotations = celltype_manhattan_plot.load_celltype_annotations()
            gwas_data = celltype_manhattan_plot.load_gwas_with_positions()
            annotations = celltype_manhattan_plot.compute_snp_overlaps(gwas_data, annotations)
            celltype_manhattan_plot.create_celltype_manhattan_plots(gwas_data, annotations)
            celltype_manhattan_plot.create_comparison_manhattan(gwas_data, annotations)
            print("✅ 시각화 완료")
        except Exception as e:
            print(f"❌ 시각화 실패: {e}")
            import traceback
            traceback.print_exc()
        
    else:
        parser.print_help()

if __name__ == "__main__":
    main()