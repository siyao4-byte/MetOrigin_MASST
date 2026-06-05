#!/usr/bin/env python
"""
Project front-door workflow for MAPS + MetOrigin + domainMASST.

Modes
-----
1. HGMD mode:
   Find an HGMD_xxxx folder under the HGMD LCMS base path, copy:
     - final-annotation-df.csv -> data/maps/final-annotation-df.csv
     - mzmine/data_iimn_gnps.mgf -> data/mgf/data_iimn_gnps.mgf

2. User-file/default mode:
   Use the files already present in data/maps/ and data/mgf/.

After input setup, optionally run:
   - MetOrigin matching
   - domainMASST preparation or full domainMASST batch
   - compact R graph workflow
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_HGMD_BASE = Path(r"Y:\MA_BPA_Microbiome\LCMS data")
DEFAULT_MAPS_INPUT = PROJECT_ROOT / "data" / "maps" / "final-annotation-df.csv"
DEFAULT_MGF_INPUT = PROJECT_ROOT / "data" / "mgf" / "data_iimn_gnps.mgf"
CONFIG_FILE = PROJECT_ROOT / "config" / "config.yaml"
PYTHON_EXE = PROJECT_ROOT / ".venv310" / "Scripts" / "python.exe"
RSCRIPT_EXE = Path(r"C:\Program Files\R\R-4.5.1\bin\Rscript.exe")
GRAPH_DIRS_TO_SAVE = [
    "Graph1_MetOrigin_R_outputs",
    "Graph2_domainMASST_R_outputs",
    "Graph3_integrated_origin_R_outputs",
]


def prompt_yes_no(question: str, default: bool = False) -> bool:
    suffix = " [Y/n] " if default else " [y/N] "
    ans = input(question + suffix).strip().lower()
    if not ans:
        return default
    return ans in {"y", "yes", "true", "1"}


def normalize_hgmd_id(value: str) -> str:
    value = value.strip()
    if not value:
        return ""
    if value.upper().startswith("HGMD_"):
        return value.upper()
    return f"HGMD_{value}"


def sanitize_run_label(value: str) -> str:
    allowed = []
    for ch in str(value).strip():
        if ch.isalnum() or ch in {"_", "-", "."}:
            allowed.append(ch)
        else:
            allowed.append("_")
    label = "".join(allowed).strip("._-")
    return label or datetime.now().strftime("USER_FILES_%Y%m%d_%H%M%S")


def find_hgmd_folder(hgmd_id: str, base_path: Path) -> Path:
    if not base_path.exists():
        raise FileNotFoundError(f"HGMD base path does not exist: {base_path}")

    matches = sorted(
        p for p in base_path.iterdir()
        if p.is_dir() and p.name.upper().startswith(hgmd_id.upper())
    )

    if not matches:
        raise FileNotFoundError(f"No folder starting with {hgmd_id} found under {base_path}")

    if len(matches) == 1:
        return matches[0]

    print("\nMultiple matching HGMD folders found:")
    for i, path in enumerate(matches, start=1):
        print(f"  {i}. {path.name}")

    while True:
        choice = input("Choose folder number: ").strip()
        if choice.isdigit() and 1 <= int(choice) <= len(matches):
            return matches[int(choice) - 1]
        print("Please enter a valid number.")


def find_annotation_file(hgmd_folder: Path) -> Path:
    candidates = list(hgmd_folder.rglob("final-annotation-df.csv"))
    if not candidates:
        raise FileNotFoundError(f"No final-annotation-df.csv found under {hgmd_folder}")
    if len(candidates) == 1:
        return candidates[0]

    candidates = sorted(candidates, key=lambda p: (len(p.parts), str(p)))
    print("\nMultiple final-annotation-df.csv files found; using the shortest path:")
    for path in candidates[:5]:
        print(f"  {path}")
    return candidates[0]


def find_mgf_file(hgmd_folder: Path) -> Path:
    mzmine_dir = hgmd_folder / "mzmine"
    if not mzmine_dir.exists():
        raise FileNotFoundError(f"Expected mzmine folder was not found: {mzmine_dir}")

    candidates = sorted(mzmine_dir.rglob("data_iimn_gnps.mgf"))
    if not candidates:
        raise FileNotFoundError(f"No data_iimn_gnps.mgf found under {mzmine_dir}")
    if len(candidates) > 1:
        print("\nMultiple data_iimn_gnps.mgf files found under mzmine; using the shortest path:")
        for path in candidates[:5]:
            print(f"  {path.relative_to(hgmd_folder)} ({path.stat().st_size:,} bytes)")
        candidates = sorted(candidates, key=lambda p: (len(p.parts), str(p)))
    return candidates[0]


def backup_existing_inputs(paths: Iterable[Path]) -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_dir = PROJECT_ROOT / "backups" / f"input_files_before_setup_{stamp}"
    copied_any = False

    for path in paths:
        if path.exists():
            backup_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, backup_dir / path.name)
            copied_any = True

    return backup_dir if copied_any else backup_dir


def copy_inputs(annotation_file: Path, mgf_file: Path) -> None:
    DEFAULT_MAPS_INPUT.parent.mkdir(parents=True, exist_ok=True)
    DEFAULT_MGF_INPUT.parent.mkdir(parents=True, exist_ok=True)

    backup_dir = backup_existing_inputs([DEFAULT_MAPS_INPUT, DEFAULT_MGF_INPUT])
    if backup_dir.exists():
        print(f"Existing input files backed up to: {backup_dir}")

    shutil.copy2(annotation_file, DEFAULT_MAPS_INPUT)
    shutil.copy2(mgf_file, DEFAULT_MGF_INPUT)

    print("\nInput files installed:")
    print(f"  MAPS annotation: {DEFAULT_MAPS_INPUT}")
    print(f"  MGF:             {DEFAULT_MGF_INPUT}")


def copy_file_if_exists(src: Path, dst: Path) -> None:
    if not src.exists():
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def copy_dir_if_exists(src: Path, dst: Path) -> None:
    if not src.exists():
        return
    shutil.copytree(src, dst, dirs_exist_ok=True)


def save_run_bundle(run_label: str, include_outputs: bool = True, include_graphs: bool = True) -> Path:
    run_label = sanitize_run_label(run_label)
    runs_root = PROJECT_ROOT / "runs"
    run_dir = runs_root / run_label
    runs_root.mkdir(parents=True, exist_ok=True)

    if run_dir.exists():
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        previous_dir = runs_root / f"{run_label}_previous_{stamp}"
        shutil.move(str(run_dir), str(previous_dir))
        print(f"Existing run folder moved to: {previous_dir}")

    run_dir.mkdir(parents=True, exist_ok=True)

    copy_file_if_exists(DEFAULT_MAPS_INPUT, run_dir / "data" / "maps" / DEFAULT_MAPS_INPUT.name)
    copy_file_if_exists(DEFAULT_MGF_INPUT, run_dir / "data" / "mgf" / DEFAULT_MGF_INPUT.name)

    saved_dirs = []

    if include_outputs:
        copy_dir_if_exists(PROJECT_ROOT / "outputs", run_dir / "outputs")
        if (PROJECT_ROOT / "outputs").exists():
            saved_dirs.append("outputs")

    if include_graphs:
        for rel_dir in GRAPH_DIRS_TO_SAVE:
            copy_dir_if_exists(PROJECT_ROOT / rel_dir, run_dir / rel_dir)
            if (PROJECT_ROOT / rel_dir).exists():
                saved_dirs.append(rel_dir)

    manifest = run_dir / "RUN_MANIFEST.txt"
    manifest.write_text(
        "\n".join([
            f"run_label: {run_label}",
            f"created_at: {datetime.now().isoformat(timespec='seconds')}",
            f"project_root: {PROJECT_ROOT}",
            "",
            "Saved folders:",
            "  data/maps/final-annotation-df.csv",
            "  data/mgf/data_iimn_gnps.mgf",
            *[f"  {x}/" for x in saved_dirs],
            "",
        ]),
        encoding="utf-8",
    )

    print(f"\nSaved run bundle to: {run_dir}")
    return run_dir


def run_command(cmd: list[str], label: str) -> None:
    print("\n" + "=" * 90)
    print(label)
    print(" ".join(f'"{x}"' if " " in x else x for x in cmd))
    print("=" * 90)
    subprocess.run(cmd, cwd=PROJECT_ROOT, check=True)


def run_metorigin() -> None:
    if not PYTHON_EXE.exists():
        raise FileNotFoundError(f"Python environment not found: {PYTHON_EXE}")
    run_command(
        [str(PYTHON_EXE), "scripts/metorigin_match.py", "--config", str(CONFIG_FILE)],
        "Running MetOrigin matching",
    )


def run_domainmasst(run_batch: bool) -> None:
    cmd = [
        "powershell",
        "-ExecutionPolicy", "Bypass",
        "-File", "scripts/run_domainmasst_unmatched.ps1",
    ]
    if run_batch:
        cmd.append("-RunMasst")
    run_command(cmd, "Running domainMASST batch" if run_batch else "Preparing domainMASST unmatched MS2 input")


def run_graphs() -> None:
    if not RSCRIPT_EXE.exists():
        raise FileNotFoundError(f"Rscript not found: {RSCRIPT_EXE}")
    run_command(
        [str(RSCRIPT_EXE), "scripts/Graph_All_compact_workflow.R"],
        "Running compact R graph workflow",
    )


def print_current_inputs() -> None:
    print("\nCurrent default input files:")
    for label, path in [("MAPS annotation", DEFAULT_MAPS_INPUT), ("MGF", DEFAULT_MGF_INPUT)]:
        status = f"{path.stat().st_size:,} bytes" if path.exists() else "MISSING"
        print(f"  {label}: {path} ({status})")


def configure_hgmd_mode(args: argparse.Namespace) -> None:
    hgmd_id = normalize_hgmd_id(args.hgmd or input("Enter HGMD number or ID, e.g. 1047 or HGMD_1047: "))
    if not hgmd_id:
        raise ValueError("HGMD ID was empty.")

    base_path = Path(args.hgmd_base)
    hgmd_folder = find_hgmd_folder(hgmd_id, base_path)
    annotation_file = find_annotation_file(hgmd_folder)
    mgf_file = find_mgf_file(hgmd_folder)

    print("\nHGMD input candidate:")
    print(f"  Requested ID:     {hgmd_id}")
    print(f"  Folder:           {hgmd_folder}")
    print(f"  Annotation CSV:   {annotation_file}")
    print(f"  MGF:              {mgf_file}")

    if not args.yes and not prompt_yes_no("Is this the correct HGMD input?", default=False):
        raise SystemExit("Stopped before copying files.")

    copy_inputs(annotation_file, mgf_file)


def configure_default_mode(args: argparse.Namespace) -> None:
    print("\nNo HGMD ID was specified. Using default input files already in the data folder.")
    print_current_inputs()

    missing = [p for p in [DEFAULT_MAPS_INPUT, DEFAULT_MGF_INPUT] if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing default input file(s): " + ", ".join(str(p) for p in missing))

    if not args.yes and not prompt_yes_no("Continue with these default input files?", default=True):
        raise SystemExit("Stopped before analysis.")


def main() -> None:
    parser = argparse.ArgumentParser(description="Set up inputs and optionally run MAPS/MetOrigin/domainMASST/graph workflow.")
    parser.add_argument("--hgmd", help="HGMD number or ID, e.g. 1047 or HGMD_1047. Omit to use data/ defaults.")
    parser.add_argument("--hgmd_base", default=str(DEFAULT_HGMD_BASE), help="Base path containing HGMD_xxxx folders.")
    parser.add_argument("--yes", action="store_true", help="Accept prompts automatically.")
    parser.add_argument("--run_metorigin", action="store_true", help="Run MetOrigin matching after input setup.")
    parser.add_argument("--prepare_domainmasst", action="store_true", help="Prepare MetOrigin-unmatched MS2 MGF for domainMASST.")
    parser.add_argument("--run_domainmasst", action="store_true", help="Prepare and run the domainMASST batch search.")
    parser.add_argument("--run_graphs", action="store_true", help="Run the compact R graph workflow.")
    parser.add_argument("--run_label", help="Folder name under runs/. Defaults to the HGMD ID or USER_FILES_timestamp.")
    parser.add_argument("--save_run_outputs", action="store_true", help="Save current data/outputs/graph folders under runs/<run_label>/ even for default user-file mode.")
    parser.add_argument("--no_save_run", action="store_true", help="Do not save a run bundle at the end.")
    parser.add_argument("--interactive_run_prompts", action="store_true", help="Ask which analysis steps to run if no run flags are supplied.")
    args = parser.parse_args()

    run_label = sanitize_run_label(args.run_label or args.hgmd or datetime.now().strftime("USER_FILES_%Y%m%d_%H%M%S"))

    if args.hgmd:
        configure_hgmd_mode(args)
    else:
        configure_default_mode(args)

    no_run_flags = not any([args.run_metorigin, args.prepare_domainmasst, args.run_domainmasst, args.run_graphs])
    if no_run_flags and args.interactive_run_prompts:
        args.run_metorigin = prompt_yes_no("Run MetOrigin matching now?", default=False)
        args.run_domainmasst = prompt_yes_no("Run the full domainMASST batch search now?", default=False)
        if not args.run_domainmasst:
            args.prepare_domainmasst = prompt_yes_no("Prepare domainMASST unmatched MS2 input only?", default=False)
        args.run_graphs = prompt_yes_no("Run the compact R graph workflow now?", default=False)

    if args.run_metorigin:
        run_metorigin()
    if args.run_domainmasst:
        run_domainmasst(run_batch=True)
    elif args.prepare_domainmasst:
        run_domainmasst(run_batch=False)
    if args.run_graphs:
        run_graphs()

    analysis_requested = any([args.run_metorigin, args.prepare_domainmasst, args.run_domainmasst, args.run_graphs])
    should_save_run = not args.no_save_run and (args.save_run_outputs or analysis_requested)

    if should_save_run:
        save_run_bundle(
            run_label,
            include_outputs=args.save_run_outputs or args.run_metorigin or args.prepare_domainmasst or args.run_domainmasst or args.run_graphs,
            include_graphs=args.save_run_outputs or args.run_graphs,
        )

    print("\nWorkflow setup complete.")
    print_current_inputs()


if __name__ == "__main__":
    main()
