#!/usr/bin/env python
"""
plot_scalar_summary.py -- ISMIP7 NORCE AIS scalar summary plots (approach A).

Plots scalar ISMIP7 variables for all 11 experiments (C001-C011) of the
NORCE CISM ensemble in one figure per variable. Output: one PNG per
variable in Plots/Summary/.

Time series handling:
  - historical (C001/C002, 1970-2014) is concatenated with each scenario
    (ssp126/ssp370/ssp585/ctrl, 2015-...) per ESM member, giving continuous
    1970-2300 lines.
  - ocx (C011, 1990-2025) is plotted as a standalone line.

Two plot variants are produced per variable:
  - absolute values:            {var}_CISM_NORCE_{end}.png
  - anomaly relative to the last time step of the historical experiment
    (end of 2014) of the same ESM member: {var}-anom_CISM_NORCE_{end}.png
    ocx (C011) has no historical counterpart, but its anomaly is also
    referenced to the end of 2014 (interpolated from its own series).

Each variant is written for several end years (2030, 2100, 2200 and the
full range ending 2300), e.g. limnsw-anom_CISM_NORCE_2100.png.

Additionally, a sea-level contribution figure is derived from limnsw
anomalies: -361.8 Gt of mass above flotation = +1 mm sea-level
contribution (limnsw-sea_level_CISM_NORCE_{end}.png).

Run with the 'plotting' environment:
  /nird/datapeak/NS11016K/miniforge3_26/envs/plotting/bin/python plot_scalar_summary.py
Optional flags: --abs-only, --anom-only, --lab {NORCE,NCAR,CISM8}
"""

import os
import glob
import argparse

import numpy as np
import netCDF4 as nc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from config import SCALAR_YLIM_ANOM, SLC_YLIM

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Workspace root (this script lives in norce-ismip7-ais-analysis/)
HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "Plots", "Summary")

# Per-lab configuration: data root, model/lab tag, experiment table
#   id: (exp_name, ESM, member, start_year, end_year)
LABS = {
    "NORCE": {
        "data_root": os.path.join(HERE, "..", "AIS", "NORCE", "CISM", "CORE"),
        "model_lab": "CISM_NORCE",
        "experiments": {
            "C001": ("historical", "CESM2-WACCM", "m001", 1970, 2014),
            "C002": ("historical", "MRI-ESM2-0", "m002", 1970, 2014),
            "C003": ("ssp370",    "CESM2-WACCM", "m001", 2015, 2100),
            "C004": ("ssp370",    "MRI-ESM2-0", "m002", 2015, 2100),
            "C005": ("ssp126",    "CESM2-WACCM", "m001", 2015, 2300),
            "C006": ("ssp126",    "MRI-ESM2-0", "m002", 2015, 2300),
            "C007": ("ssp585",    "CESM2-WACCM", "m001", 2015, 2300),
            "C008": ("ssp585",    "MRI-ESM2-0", "m002", 2015, 2300),
            "C009": ("ctrl",      "CESM2-WACCM", "m001", 2015, 2300),
            "C010": ("ctrl",      "MRI-ESM2-0", "m002", 2015, 2300),
            "C011": ("ocx",       "ERA",        "m001", 1990, 2025),
        },
    },
    "NCAR": {
        "data_root": os.path.join(HERE, "..", "AIS", "NCAR", "CISM", "CORE"),
        "model_lab": "CISM_NCAR",
        "experiments": {
            "C001": ("historical", "CESM2-WACCM", "m001", 2000, 2014),
            "C007": ("ssp585",     "CESM2-WACCM", "m001", 2015, 2300),
        },
    },
    "CISM8": {
        "data_root": os.path.join(HERE, "..", "AIS", "NORCE", "CISM8", "CORE"),
        "model_lab": "CISM8_NORCE",
        "experiments": {
            "C001": ("historical", "CESM2-WACCM", "m001", 1970, 2014),
            "C002": ("historical", "MRI-ESM2-0", "m002", 1970, 2014),
            "C003": ("ssp370",    "CESM2-WACCM", "m001", 2015, 2100),
            "C004": ("ssp370",    "MRI-ESM2-0", "m002", 2015, 2100),
            "C005": ("ssp126",    "CESM2-WACCM", "m001", 2015, 2300),
            "C006": ("ssp126",    "MRI-ESM2-0", "m002", 2015, 2300),
            "C007": ("ssp585",    "CESM2-WACCM", "m001", 2015, 2300),
            "C008": ("ssp585",    "MRI-ESM2-0", "m002", 2015, 2300),
            "C009": ("ctrl",      "CESM2-WACCM", "m001", 2015, 2300),
            "C010": ("ctrl",      "MRI-ESM2-0", "m002", 2015, 2300),
            "C011": ("ocx",       "ERA",        "m001", 1990, 2025),
        },
    },
}

# Scenario colors (historical/scenario families); ocx gets its own color
SCENARIO_COLORS = {
    "historical": "black",
    "ssp126":     "tab:blue",
    "ssp370":     "tab:orange",
    "ssp585":     "tab:red",
    "ctrl":       "tab:green",
    "ocx":        "tab:purple",
}

# Line style per ESM member
MEMBER_STYLES = {
    "m001": "-",
    "m002": "--",
}

# Variables to plot: var -> (long_name, plot_unit, scale, ftype)
#   scale converts raw file units to plot units:
#     kg      -> Gt       : * 1e-12
#     kg s-1  -> Gt/yr    : * 1e-12 * 86400 * 365.25
#     m^2     -> 10^6 km^2: * 1e-12
#   ftype: 'ST' = snapshot (Jan 1), 'FL' = annual mean flux (Jul 1)
VARIABLES = {
    "lim":            ("total ice mass",                    "Gt",        1e-12, "ST"),
    "limnsw":         ("mass above flotation",              "Gt",        1e-12, "ST"),
    "iareagr":        ("grounded ice area",                 "10$^6$ km$^2$", 1e-12, "ST"),
    "iareafl":        ("floating ice area",                 "10$^6$ km$^2$", 1e-12, "ST"),
    "tendacabf":      ("total SMB flux",                    "Gt/yr",     1e-12 * 86400 * 365.25, "FL"),
    "tendlibmassbfgr": ("total grounded basal mass balance flux", "Gt/yr", 1e-12 * 86400 * 365.25, "FL"),
    "tendlibmassbffl": ("total floating basal mass balance flux", "Gt/yr", 1e-12 * 86400 * 365.25, "FL"),
    "tendlicalvf":    ("total calving flux",                "Gt/yr",     1e-12 * 86400 * 365.25, "FL"),
}

# All-zero variables, skipped on purpose
SKIP_VARS = {"tendlifmassbf", "tendligroundf"}

# Sea-level contribution conversion: -361.8 Gt limnsw = +1 mm SLC
GT_PER_MM_SL = 361.8

# Model/lab tag and end years for output file names
MODEL_LAB = "CISM_NORCE"  # default; overridden by --lab
END_YEARS = (2030, 2100, 2200, 2300)

# Active lab configuration (set in main())
DATA_ROOT = None
EXPERIMENTS = {}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def find_file(exp_id, var):
    """Return path of the scalar nc file for var in experiment exp_id, or None."""
    pattern = os.path.join(DATA_ROOT, exp_id, f"{var}_AIS_*_{exp_id}_*.nc")
    matches = sorted(glob.glob(pattern))
    return matches[0] if matches else None


def read_scalar(path, var):
    """Read (years, values) from a scalar nc file.

    Time is stored as days since 1850-01-01 (standard calendar); converted
    to decimal years here. Returns (years ndarray, values ndarray) or
    (None, None) if the file/variable is unusable.
    """
    try:
        with nc.Dataset(path) as ds:
            if var not in ds.variables or "time" not in ds.variables:
                print(f"  WARNING: {var} or time missing in {os.path.basename(path)}")
                return None, None
            t = ds.variables["time"][:]
            d = ds.variables[var][:]
            units = ds.variables["time"].units  # e.g. 'days since 1850-01-01'
    except OSError as e:
        print(f"  WARNING: cannot read {path}: {e}")
        return None, None

    # days since <origin> -> decimal years (standard calendar, 365.25 d/yr avg)
    origin_year = float(units.split("since")[1].split("-")[0])
    years = origin_year + np.asarray(t, dtype=float) / 365.25
    vals = np.ma.filled(np.asarray(d, dtype=float), np.nan)
    return years, vals


def concatenate_hist_scenario(hist, scen):
    """Concatenate historical and scenario series, dropping overlapping years.

    hist, scen: (years, values) tuples. Returns merged (years, values).
    """
    hy, hv = hist
    sy, sv = scen
    # keep historical points strictly before the first scenario point
    keep = hy < sy[0]
    years = np.concatenate([hy[keep], sy])
    vals = np.concatenate([hv[keep], sv])
    return years, vals


def build_series(var):
    """Build the list of line series for one variable.

    Each series dict has keys: label, color, style, years, vals, ref.
    'ref' is the anomaly reference value: the last time step of the
    historical experiment of the same ESM member (end of 2014). For ocx
    (no historical counterpart) the value at end of 2014 is interpolated
    from its own series; for scenarios without a historical run, the
    first time step of the series itself is used.
    """
    series = []
    hist = {}  # member -> (years, vals)

    for exp_id, (exp, esm, member, _, _) in EXPERIMENTS.items():
        path = find_file(exp_id, var)
        if path is None:
            print(f"  WARNING: no file for {var} in {exp_id}, skipped")
            continue
        years, vals = read_scalar(path, var)
        if years is None:
            continue

        if exp == "historical":
            hist[member] = (years, vals)
        elif exp == "ocx":
            # standalone line, anomaly referenced to end of 2014
            # (interpolated from its own series)
            series.append({
                "label": f"{exp_id} {exp}",
                "color": SCENARIO_COLORS["ocx"],
                "style": MEMBER_STYLES.get(member, "-"),
                "years": years, "vals": vals,
                "ref": np.interp(2015.0, years, vals),
            })
        else:
            if member in hist:
                merged = concatenate_hist_scenario(hist[member], (years, vals))
                # anomaly reference: last historical time step (end of 2014)
                ref = hist[member][1][-1]
            else:
                print(f"  WARNING: no historical for member {member}, "
                      f"plotting {exp_id} standalone")
                merged = (years, vals)
                ref = vals[0]
            series.append({
                "label": f"{exp_id} {exp} {esm}",
                "color": SCENARIO_COLORS.get(exp, None),
                "style": MEMBER_STYLES.get(member, "-"),
                "years": merged[0], "vals": merged[1],
                "ref": ref,
            })

    return series


def plot_variable(var, long_name, unit, scale, ftype, series, anom=False):
    """Create and save summary figures for a variable.

    series: list of dicts with keys: label, color, style, years, vals, ref.
    If anom is True, plot (vals - ref) instead of absolute values.
    One figure is written per end year in END_YEARS (data beyond the end
    year is masked out).
    """
    for end_year in END_YEARS:
        fig, ax = plt.subplots(figsize=(10, 6))

        for s in series:
            vals = (s["vals"] - s["ref"]) if anom else s["vals"]
            keep = s["years"] <= end_year
            ax.plot(s["years"][keep], vals[keep] * scale, color=s["color"],
                    linestyle=s["style"], linewidth=1.5, label=s["label"])

        ylims = SCALAR_YLIM_ANOM.get(var, {}) if anom else {}
        ylim = ylims.get(end_year)
        if ylim is not None:
            ax.set_ylim(*ylim)

        if anom:
            ax.axhline(0.0, color="gray", linewidth=0.8)
            ax.set_title(f"{var}: {long_name} anomaly ({ftype}-type)",
                         fontsize=12)
            ax.set_ylabel(f"$\\Delta$ {long_name} [{unit}]")
            fname = f"{var}-anom_{MODEL_LAB}_{end_year}.png"
        else:
            ax.set_title(f"{var}: {long_name}  ({ftype}-type)", fontsize=13)
            ax.set_ylabel(f"{long_name} [{unit}]")
            fname = f"{var}_{MODEL_LAB}_{end_year}.png"

        ax.set_xlabel("year")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=8, ncol=2)
        fig.tight_layout()

        out = os.path.join(OUT_DIR, fname)
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(f"  wrote {out}")


# ---------------------------------------------------------------------------
# Sea-level contribution
# ---------------------------------------------------------------------------


def plot_sea_level_contribution(series, end_year):
    """Plot sea-level contribution from limnsw anomalies.

    Mass above flotation lost by the ice sheet raises sea level:
    -361.8 Gt limnsw = +1 mm SLC. Anomalies are relative to the last
    historical time step (end of 2014) of the same ESM member; ocx is
    referenced to end of 2014 (interpolated from its own series).
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    for s in series:
        # kg -> Gt (1e-12), then Gt -> mm SLC via -361.8 Gt = +1 mm
        slc = -(s["vals"] - s["ref"]) * 1e-12 / GT_PER_MM_SL
        keep = s["years"] <= end_year
        ax.plot(s["years"][keep], slc[keep], color=s["color"],
                linestyle=s["style"], linewidth=1.5, label=s["label"])

    ax.axhline(0.0, color="gray", linewidth=0.8)
    ax.set_title("Sea-level contribution", fontsize=13)
    ax.set_xlabel("year")
    ax.set_ylabel("sea-level contribution [mm]")
    if end_year in SLC_YLIM:
        ax.set_ylim(*SLC_YLIM[end_year])
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8, ncol=2)
    fig.tight_layout()

    out = os.path.join(OUT_DIR, f"limnsw-sea_level_{MODEL_LAB}_{end_year}.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    global MODEL_LAB, DATA_ROOT, EXPERIMENTS

    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--abs-only", action="store_true",
                       help="plot only absolute values")
    group.add_argument("--anom-only", action="store_true",
                       help="plot only anomalies rel. to end of historical")
    parser.add_argument("--lab", choices=sorted(LABS), default="NORCE",
                        help="which lab's data to plot (default: NORCE)")
    args = parser.parse_args()

    lab_cfg = LABS[args.lab]
    MODEL_LAB = lab_cfg["model_lab"]
    DATA_ROOT = os.path.normpath(lab_cfg["data_root"])
    EXPERIMENTS = lab_cfg["experiments"]
    print(f"Lab: {args.lab} ({MODEL_LAB}), data root: {DATA_ROOT}")

    do_abs = not args.anom_only
    do_anom = not args.abs_only

    os.makedirs(OUT_DIR, exist_ok=True)

    for var, (long_name, unit, scale, ftype) in VARIABLES.items():
        print(f"Processing {var} ...")

        series = build_series(var)
        if not series:
            print(f"  WARNING: no data at all for {var}, no figure written")
            continue

        if do_abs:
            plot_variable(var, long_name, unit, scale, ftype, series, anom=False)
        if do_anom:
            plot_variable(var, long_name, unit, scale, ftype, series, anom=True)

        if var == "limnsw":
            for end_year in END_YEARS:
                plot_sea_level_contribution(series, end_year)

    print("Done.")


if __name__ == "__main__":
    main()
