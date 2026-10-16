"""Reproducible indexed matcher and optional R/HTML reporting entry point."""
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
import yaml
from indexed_match import run
from metorigin_index import build_index


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    parser.add_argument('--plots', action='store_true')
    parser.add_argument('--html', action='store_true')
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--rscript', default=os.environ.get('RSCRIPT') or shutil.which('Rscript') or 'Rscript',
                        help='Rscript executable; defaults to RSCRIPT or PATH.')
    args = parser.parse_args()
    config_path = Path(args.config).resolve()
    project = config_path.parent.parent
    cfg = yaml.safe_load(config_path.read_text(encoding='utf-8'))
    index = project / cfg['metorigin_index']
    if not index.exists():
        build_index(project / cfg['metorigin_file'], index)
    run_dir = run(config_path, validate=not args.offline)
    env = dict(os.environ, METORIGIN_RUN_DIR=str(run_dir), METORIGIN_INCLUDE_DOMAINMASST='false')
    if args.plots:
        subprocess.run([args.rscript,'scripts/Graph_All_compact_workflow.R'],cwd=project,env=env,check=True)
    if args.html:
        subprocess.run([sys.executable,'scripts/generate_interactive_origin_report.py','--run',str(run_dir)],cwd=project,check=True)


if __name__ == '__main__':main()
