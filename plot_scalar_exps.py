#!/usr/bin/env python
"""
plot_scalar_exps.py -- ISMIP7 AIS per-experiment plots (approach B).

One plot per experiment: 2D maps of lithk (ice thickness) change relative
to the end of the historical experiment (end of 2014) of the same ESM
member. Output: one PNG per experiment and end year in Plots/Exps/, e.g.

  lithk-anom_C008_CISM_NORCE_2300.png

The historical reference field is taken from the historical run of the
same ESM member (C001/C002); the anomaly is (lithk_scenario -
lithk_hist_end). Experiments without a historical counterpart (ocx) are
skipped.

Run with the 'plotting' environment:
  /nird/datapeak/NS11016K/miniforge3_26/envs/plotting/bin/python plot_scalar_exps.py
Optional flags: --lab {NORCE,NCAR}, --model {CISM,CISM8},
                --exp C008 (repeatable)
"""

import os
import glob
import argparse

import numpy as np
import netCDF4 as nc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, ListedColormap, BoundaryNorm

from config import LITHK_VMAX

# Mask variables: var name -> (short tag, long name)
MASK_VARS = {
    "sftgrf": ("grounded-mask", "grounded ice mask"),
    "sftflf": ("floating-mask", "floating ice mask"),
    "sftgif": ("ice-mask", "ice mask"),
}
MASK_THRESH = 0.5  # fractional fields (sftgrf/sftflf) -> binary at 0.5
OROG_SNAPSHOT_YEARS = (2015, 2030, 2100, 2200, 2300)
OROG_LEVELS = np.arange(0, 4500, 500)

# ---------------------------------------------------------------------------
# Configuration (kept consistent with plot_scalar_summary.py)
# ---------------------------------------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "Plots", "Exps")

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

# End years for output file names (2300 = full range)
END_YEARS = (2100, 2200, 2300)

# Historical member per ESM: member id -> historical experiment id
# (derived from the experiment table at runtime)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def find_file(data_root, exp_id, var):
    """Return path of the nc file for var in experiment exp_id, or None."""
    pattern = os.path.join(data_root, exp_id, f"{var}_AIS_*_{exp_id}_*.nc")
    matches = sorted(glob.glob(pattern))
    return matches[0] if matches else None


def read_field(path, var):
    """Read (years, x, y, field) from a gridded nc file (time, y, x)."""
    with nc.Dataset(path) as ds:
        t = ds.variables["time"][:]
        x = ds.variables["x"][:]
        y = ds.variables["y"][:]
        d = ds.variables[var][:]
        units = ds.variables["time"].units
    origin_year = float(units.split("since")[1].split("-")[0])
    years = origin_year + np.asarray(t, dtype=float) / 365.25
    return years, np.asarray(x), np.asarray(y), np.ma.filled(d, np.nan)


def read_lithk(path):
    """Read (years, x, y, lithk) from a lithk nc file."""
    return read_field(path, "lithk")


def hist_member_map(experiments):
    """Map ESM member -> historical experiment id."""
    out = {}
    for exp_id, (exp, _, member, _, _) in experiments.items():
        if exp == "historical":
            out[member] = exp_id
    return out


def plot_lithk_anom(exp_id, exp, esm, member, end_year, ref_field,
                    years, x, y, lithk, model_lab, out_dir):
    """Plot lithk anomaly (relative to end of historical) for one experiment.

    Only cells with ice in either field (ref or current) are shown; the
    rest is masked. A diverging colormap centered at 0 is used.
    """
    keep_time = years <= end_year
    if not keep_time.any():
        print(f"  WARNING: no data until {end_year} for {exp_id}, skipped")
        return
    idx = np.nonzero(keep_time)[0][-1]
    field = lithk[idx] - ref_field

    # mask cells without ice in both fields (thickness 0 / NaN)
    ice = (ref_field > 0) | (lithk[idx] > 0)
    field = np.where(ice, field, np.nan)

    fig, ax = plt.subplots(figsize=(8, 8))
    vmax = LITHK_VMAX[end_year]
    field_max = np.nanmax(np.abs(field))
    if np.isfinite(field_max) and field_max > vmax:
        print(f"  WARNING: {exp_id} max |lithk anomaly| {field_max:.1f} m "
              f"exceeds vmax {vmax:.0f} m for {end_year} -- "
              "update LITHK_VMAX in config.py")
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    im = ax.pcolormesh(x / 1000, y / 1000, field, cmap="RdBu", norm=norm,
                       shading="auto", rasterized=True)
    ax.set_aspect("equal")
    ax.set_title(f"{exp_id} {exp} {esm}: lithk change [m] "
                 f"({int(years[idx])})", fontsize=11)
    ax.set_xlabel("x [km]")
    ax.set_ylabel("y [km]")
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("$\\Delta$ lithk [m]")

    fig.tight_layout()
    fname = f"lithk-anom_{exp_id}_{model_lab}_{end_year}.png"
    out = os.path.join(out_dir, fname)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def plot_mask_change(exp_id, exp, esm, member, end_year, var, tag, long_name,
                     ref_mask, years, x, y, mask, model_lab, out_dir):
    """Plot change of a binary ice mask relative to end of historical.

    Categories: 0 = no change, 1 = mask gained (new), -1 = mask lost.
    """
    keep_time = years <= end_year
    if not keep_time.any():
        print(f"  WARNING: no data until {end_year} for {exp_id}, skipped")
        return
    idx = np.nonzero(keep_time)[0][-1]
    cur = mask[idx] > MASK_THRESH
    change = np.where(cur & ~ref_mask, 1, np.where(~cur & ref_mask, -1, 0))

    fig, ax = plt.subplots(figsize=(8, 8))
    cmap = ListedColormap(["#d73027", "#f7f7f7", "#1a9850"])  # lost, same, gained
    norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], cmap.N)
    im = ax.pcolormesh(x / 1000, y / 1000, change, cmap=cmap, norm=norm,
                       shading="auto", rasterized=True)
    ax.set_aspect("equal")
    ax.set_title(f"{exp_id} {exp} {esm}: {long_name} change "
                 f"({int(years[idx])})", fontsize=11)
    ax.set_xlabel("x [km]")
    ax.set_ylabel("y [km]")
    cbar = fig.colorbar(im, ax=ax, shrink=0.8, ticks=[-1, 0, 1])
    cbar.ax.set_yticklabels(["lost", "unchanged", "gained"])

    fig.tight_layout()
    fname = f"{tag}-anom_{exp_id}_{model_lab}_{end_year}.png"
    out = os.path.join(out_dir, fname)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def plot_orog_snapshots(data_root, exp_id, exp, esm, model_lab, out_dir):
    """Plot masked surface elevation snapshots for an experiment."""
    paths = {var: find_file(data_root, exp_id, var)
             for var in ("orog", "sftflf", "sftgif")}
    missing = [var for var, path in paths.items() if path is None]
    if missing:
        print(f"WARNING: missing {', '.join(missing)} for {exp_id}, "
              "orog snapshots skipped")
        return

    with nc.Dataset(paths["orog"]) as orog_ds, \
            nc.Dataset(paths["sftflf"]) as floating_ds, \
            nc.Dataset(paths["sftgif"]) as ice_ds:
        time_var = orog_ds.variables["time"]
        times = time_var[:]
        dates = nc.num2date(
            times, time_var.units,
            calendar=getattr(time_var, "calendar", "standard"))
        # ST fields are dated Jan 1 of the year after the represented model year.
        sample_years = np.asarray([date.year - 1 for date in dates])
        if (len(floating_ds.variables["time"]) != len(times)
                or len(ice_ds.variables["time"]) != len(times)):
            print(f"WARNING: inconsistent time dimensions for {exp_id}, "
                  "orog snapshots skipped")
            return

        x = np.asarray(orog_ds.variables["x"][:])
        y = np.asarray(orog_ds.variables["y"][:])
        for snapshot_year in OROG_SNAPSHOT_YEARS:
            valid = np.flatnonzero(sample_years <= snapshot_year)
            if not valid.size:
                print(f"  WARNING: no orog data through {snapshot_year} for "
                      f"{exp_id}, skipped")
                continue
            index = valid[-1]
            field = np.ma.filled(orog_ds.variables["orog"][index], np.nan)
            floating = np.ma.filled(
                floating_ds.variables["sftflf"][index], np.nan) > MASK_THRESH
            ice = np.ma.filled(ice_ds.variables["sftgif"][index], np.nan) > MASK_THRESH
            field = np.where(floating, OROG_LEVELS[0], field)
            field = np.where(ice, field, np.nan)
            actual_year = int(sample_years[index])

            fig, ax = plt.subplots(figsize=(8, 8))
            cmap = plt.get_cmap("terrain", len(OROG_LEVELS) - 1)
            cmap.set_bad("lightgrey")
            norm = BoundaryNorm(OROG_LEVELS, cmap.N, clip=True)
            im = ax.pcolormesh(x / 1000, y / 1000, field, cmap=cmap,
                               norm=norm, shading="auto", rasterized=True)
            ax.set_aspect("equal")
            title = (f"{exp_id} {exp} {esm}: surface elevation [m] "
                     f"({actual_year})")
            if actual_year != snapshot_year:
                title += f"; latest available for {snapshot_year}"
            ax.set_title(title, fontsize=11)
            ax.set_xlabel("x [km]")
            ax.set_ylabel("y [km]")
            cbar = fig.colorbar(im, ax=ax, shrink=0.8)
            cbar.set_label("surface elevation [m]")
            fig.tight_layout()

            fname = f"orog_{exp_id}_{model_lab}_{snapshot_year}.png"
            out = os.path.join(out_dir, fname)
            fig.savefig(out, dpi=150)
            plt.close(fig)
            print(f"  wrote {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lab", choices=("NORCE", "NCAR"), default="NORCE",
                        help="which lab's data to plot (default: NORCE)")
    parser.add_argument("--model", choices=("CISM", "CISM8"),
                        default="CISM",
                        help="which model resolution to plot (default: CISM)")
    parser.add_argument("--exp", action="append", default=None,
                        help="experiment id to plot (repeatable, e.g. "
                             "--exp C008); default: all non-historical exps")
    args = parser.parse_args()

    if args.model == "CISM8" and args.lab != "NORCE":
        parser.error("--model CISM8 is currently available only for --lab NORCE")
    dataset_key = "CISM8" if args.model == "CISM8" else args.lab
    lab_cfg = LABS[dataset_key]
    data_root = os.path.normpath(lab_cfg["data_root"])
    model_lab = lab_cfg["model_lab"]
    experiments = lab_cfg["experiments"]
    out_dir = os.path.join(OUT_DIR, model_lab)
    print(f"Lab: {args.lab}, model: {args.model} ({model_lab}), "
          f"data root: {data_root}")

    os.makedirs(out_dir, exist_ok=True)

    # historical reference fields per member (end of historical run)
    hist_map = hist_member_map(experiments)
    ref_fields = {}
    ref_masks = {}
    for member, hist_id in hist_map.items():
        path = find_file(data_root, hist_id, "lithk")
        if path is None:
            print(f"WARNING: no lithk file for historical {hist_id}")
            continue
        years, x, y, lithk = read_lithk(path)
        ref_fields[member] = lithk[-1]
        print(f"Reference field: {hist_id} at year {int(years[-1])}")
        # mask references (grounded / floating / ice)
        masks = {}
        for var, (tag, long_name) in MASK_VARS.items():
            mpath = find_file(data_root, hist_id, var)
            if mpath is None:
                print(f"WARNING: no {var} file for historical {hist_id}")
                continue
            _, _, _, m = read_field(mpath, var)
            masks[var] = m[-1] > MASK_THRESH
        ref_masks[member] = masks

    # experiments to plot: all non-historical by default
    if args.exp:
        targets = {e: experiments[e] for e in args.exp if e in experiments}
    else:
        targets = {e: info for e, info in experiments.items()
                   if info[0] != "historical"}

    for exp_id, (exp, esm, member, _, _) in targets.items():
        if member not in ref_fields:
            print(f"WARNING: no historical reference for {exp_id} "
                  f"(member {member}), skipped")
            continue
        path = find_file(data_root, exp_id, "lithk")
        if path is None:
            print(f"WARNING: no lithk file for {exp_id}, skipped")
            continue
        print(f"Processing {exp_id} {exp} ...")
        plot_orog_snapshots(data_root, exp_id, exp, esm, model_lab, out_dir)
        years, x, y, lithk = read_lithk(path)
        ref_field = ref_fields[member]
        for end_year in END_YEARS:
            plot_lithk_anom(exp_id, exp, esm, member, end_year, ref_field,
                            years, x, y, lithk, model_lab, out_dir)
        # mask change plots
        for var, (tag, long_name) in MASK_VARS.items():
            mpath = find_file(data_root, exp_id, var)
            if mpath is None:
                print(f"WARNING: no {var} file for {exp_id}, skipped")
                continue
            if var not in ref_masks.get(member, {}):
                print(f"WARNING: no {var} historical reference for {exp_id}, "
                      "skipped")
                continue
            _, _, _, m = read_field(mpath, var)
            ref_mask = ref_masks[member][var]
            for end_year in END_YEARS:
                plot_mask_change(exp_id, exp, esm, member, end_year, var,
                                 tag, long_name, ref_mask, years, x, y, m,
                                 model_lab, out_dir)

    print("Done.")


if __name__ == "__main__":
    main()
