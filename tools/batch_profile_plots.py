#!/usr/bin/env python3
"""
batch_profile_plots.py
Recursively finds every sample_diameter_profile.json under a results directory
and generates a sample_profile_I_Apm.html in the same folder.

Usage:
    python batch_profile_plots.py                        # uses ./results/
    python batch_profile_plots.py /path/to/results       # explicit path
    python batch_profile_plots.py --force /path/to/results  # overwrite existing
"""
import argparse
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PLOTTER    = SCRIPT_DIR / "plot_runlog_compact.py"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results_dir", nargs="?", default=None,
                    help="Root results directory to scan (default: ./results/)")
    ap.add_argument("--force", action="store_true",
                    help="Overwrite existing plots")
    args = ap.parse_args()

    results_dir = Path(args.results_dir) if args.results_dir else SCRIPT_DIR / "results"
    if not results_dir.is_dir():
        sys.exit(f"ERROR: results directory not found: {results_dir}")

    json_files = sorted(results_dir.rglob("sample_diameter_profile.json"))
    if not json_files:
        sys.exit(f"No sample_diameter_profile.json files found under {results_dir}")

    print(f"Scanning: {results_dir}")
    print(f"Found {len(json_files)} profile(s)\n")

    found = len(json_files)
    done = skipped = failed = 0

    for json_file in json_files:
        out = json_file.parent / "sample_profile_I_Apm.html"
        label = json_file.parent.name

        if out.exists() and not args.force:
            print(f"[skip]  {label}")
            skipped += 1
            continue

        print(f"[plot]  {label}")
        result = subprocess.run(
            [sys.executable, str(PLOTTER),
             "sample-profile",
             "--json",    str(json_file),
             "--x-units", "mm",
             "--y",       "I_Apm",
             "--out",     str(out)],
            capture_output=True, text=True
        )
        if result.returncode == 0:
            done += 1
        else:
            msg = (result.stderr or result.stdout).strip()
            print(f"        ERROR: {msg}")
            failed += 1

    print(f"\nDone.  found={found}  plotted={done}  skipped={skipped}  failed={failed}")
    if skipped:
        print("(use --force to overwrite existing plots)")


if __name__ == "__main__":
    main()
