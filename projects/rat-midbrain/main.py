#!/usr/bin/env python3
"""
Main pipeline for rat midbrain cell-type-specific enhancer analysis

This script orchestrates the entire ChIP-seq analysis pipeline:
1. Pre-processing and QC
2. Alignment
3. Peak calling
4. QC metrics
5. Peak processing
6. Annotation
7. Enrichment analysis
8. Motif analysis
"""

import os
import sys
import argparse
import subprocess
import logging
from pathlib import Path
from datetime import datetime

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('pipeline.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


class RatMidbrainPipeline:
    """Main pipeline orchestrator"""

    def __init__(self, project_dir, skip_steps=None):
        self.project_dir = Path(project_dir)
        self.scripts_dir = self.project_dir / "scripts"
        self.skip_steps = skip_steps or []

        # Define pipeline steps
        self.steps = [
            {
                'name': 'preprocessing',
                'script': '01_preprocessing_qc.sh',
                'description': 'Pre-processing and Quality Control'
            },
            {
                'name': 'alignment',
                'script': '02_alignment.sh',
                'description': 'Read Alignment'
            },
            {
                'name': 'peak_calling',
                'script': '03_peak_calling.sh',
                'description': 'Peak Calling with MACS2'
            },
            {
                'name': 'qc_metrics',
                'script': '04_qc_metrics.sh',
                'description': 'QC Metrics Calculation'
            },
            {
                'name': 'peak_processing',
                'script': '05_peak_processing.sh',
                'description': 'Peak Processing'
            },
            {
                'name': 'annotation',
                'script': '06_annotation.R',
                'description': 'Peak Annotation'
            },
            {
                'name': 'enrichment',
                'script': '07_enrichment.R',
                'description': 'Functional Enrichment Analysis'
            },
            {
                'name': 'motif',
                'script': '08_motif_analysis.sh',
                'description': 'Motif Analysis'
            }
        ]

    def run_step(self, step):
        """Run a single pipeline step"""
        logger.info("=" * 60)
        logger.info(f"Step: {step['description']}")
        logger.info("=" * 60)

        script_path = self.scripts_dir / step['script']

        if not script_path.exists():
            logger.error(f"Script not found: {script_path}")
            return False

        # Make script executable
        os.chmod(script_path, 0o755)

        # Determine how to run the script
        if script_path.suffix == '.sh':
            cmd = ['bash', str(script_path)]
        elif script_path.suffix == '.R':
            cmd = ['Rscript', str(script_path)]
        else:
            logger.error(f"Unknown script type: {script_path.suffix}")
            return False

        # Run the script
        try:
            logger.info(f"Running: {' '.join(cmd)}")
            result = subprocess.run(
                cmd,
                cwd=self.project_dir,
                check=True,
                capture_output=False
            )
            logger.info(f"Step '{step['name']}' completed successfully")
            return True

        except subprocess.CalledProcessError as e:
            logger.error(f"Step '{step['name']}' failed with exit code {e.returncode}")
            return False

        except Exception as e:
            logger.error(f"Error running step '{step['name']}': {str(e)}")
            return False

    def run_pipeline(self, start_from=None, stop_at=None):
        """Run the full pipeline or a subset of steps"""
        logger.info("=" * 60)
        logger.info("Rat Midbrain ChIP-seq Analysis Pipeline")
        logger.info("=" * 60)
        logger.info(f"Project directory: {self.project_dir}")
        logger.info(f"Start time: {datetime.now()}")
        logger.info("")

        # Filter steps based on start_from and stop_at
        steps_to_run = []
        started = start_from is None

        for step in self.steps:
            if start_from and step['name'] == start_from:
                started = True

            if started and step['name'] not in self.skip_steps:
                steps_to_run.append(step)

            if stop_at and step['name'] == stop_at:
                break

        if not steps_to_run:
            logger.error("No steps to run")
            return False

        logger.info("Pipeline steps to run:")
        for i, step in enumerate(steps_to_run, 1):
            logger.info(f"  {i}. {step['description']}")
        logger.info("")

        # Run each step
        success = True
        for i, step in enumerate(steps_to_run, 1):
            logger.info("")
            logger.info(f"[{i}/{len(steps_to_run)}] Starting: {step['description']}")

            if not self.run_step(step):
                logger.error(f"Pipeline stopped at step: {step['description']}")
                success = False
                break

        # Summary
        logger.info("")
        logger.info("=" * 60)
        if success:
            logger.info("Pipeline completed successfully!")
        else:
            logger.info("Pipeline failed")
        logger.info(f"End time: {datetime.now()}")
        logger.info("=" * 60)

        return success


def main():
    parser = argparse.ArgumentParser(
        description='Run rat midbrain ChIP-seq analysis pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run full pipeline
  python main.py

  # Run from specific step
  python main.py --start-from alignment

  # Run until specific step
  python main.py --stop-at peak_calling

  # Run specific range
  python main.py --start-from alignment --stop-at qc_metrics

  # Skip specific steps
  python main.py --skip preprocessing motif

Available steps:
  preprocessing, alignment, peak_calling, qc_metrics,
  peak_processing, annotation, enrichment, motif
        """
    )

    parser.add_argument(
        '--project-dir',
        default='projects/rat-midbrain',
        help='Project directory (default: projects/rat-midbrain)'
    )

    parser.add_argument(
        '--start-from',
        choices=['preprocessing', 'alignment', 'peak_calling', 'qc_metrics',
                 'peak_processing', 'annotation', 'enrichment', 'motif'],
        help='Start pipeline from this step'
    )

    parser.add_argument(
        '--stop-at',
        choices=['preprocessing', 'alignment', 'peak_calling', 'qc_metrics',
                 'peak_processing', 'annotation', 'enrichment', 'motif'],
        help='Stop pipeline at this step'
    )

    parser.add_argument(
        '--skip',
        nargs='+',
        choices=['preprocessing', 'alignment', 'peak_calling', 'qc_metrics',
                 'peak_processing', 'annotation', 'enrichment', 'motif'],
        help='Skip these steps'
    )

    parser.add_argument(
        '--list-steps',
        action='store_true',
        help='List all pipeline steps and exit'
    )

    args = parser.parse_args()

    # Create pipeline instance
    pipeline = RatMidbrainPipeline(args.project_dir, skip_steps=args.skip)

    # List steps if requested
    if args.list_steps:
        print("\nPipeline steps:")
        print("-" * 60)
        for i, step in enumerate(pipeline.steps, 1):
            print(f"{i}. {step['name']:<20} - {step['description']}")
        print("-" * 60)
        return 0

    # Run pipeline
    success = pipeline.run_pipeline(
        start_from=args.start_from,
        stop_at=args.stop_at
    )

    return 0 if success else 1


if __name__ == '__main__':
    sys.exit(main())
