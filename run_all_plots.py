#!/usr/bin/env python
"""Generate all ISMIP7 AIS plots for a selected lab and model."""

import argparse
import os
import subprocess
import sys


HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lab", choices=("NORCE", "NCAR"), default="NORCE",
                        help="which lab's data to plot (default: NORCE)")
    parser.add_argument("--model", choices=("CISM", "CISM8"),
                        default="CISM",
                        help="which model resolution to plot (default: CISM)")
    parser.add_argument("--abs-only", action="store_true",
                        help="write only absolute scalar time series")
    parser.add_argument("--anom-only", action="store_true",
                        help="write only anomaly scalar time series")
    parser.add_argument("--exp", action="append", dest="experiments",
                        metavar="ID",
                        help="limit experiment maps; repeat for multiple IDs")
    args = parser.parse_args()

    if args.abs_only and args.anom_only:
        parser.error("--abs-only and --anom-only cannot be used together")
    if args.model == "CISM8" and args.lab != "NORCE":
        parser.error("--model CISM8 is currently available only for --lab NORCE")

    dataset_args = ["--lab", args.lab, "--model", args.model]
    summary_args = dataset_args.copy()
    if args.abs_only:
        summary_args.append("--abs-only")
    if args.anom_only:
        summary_args.append("--anom-only")

    experiment_args = dataset_args.copy()
    for experiment in args.experiments or []:
        experiment_args.extend(("--exp", experiment))

    commands = (
        ("plot_scalar_summary.py", summary_args),
        ("plot_scalar_exps.py", experiment_args),
        ("plot_initial.py", dataset_args),
    )
    for script, script_args in commands:
        print(f"\n=== {script} ===", flush=True)
        subprocess.run(
            [sys.executable, os.path.join(HERE, script), *script_args],
            cwd=HERE,
            check=True,
        )


if __name__ == "__main__":
    main()
