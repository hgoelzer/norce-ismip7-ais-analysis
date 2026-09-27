#!/usr/bin/env python
"""
plot_init_velocity.py -- ISMIP7 NORCE AIS initial state maps.

Plots the initial state (first time step of the historical experiments
C001/C002):

  - ice velocity magnitude (from xvelmean/yvelmean, m s-1 -> m/yr) on a
    log-like scale with levels 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000,
    3000 m/yr
  - surface elevation (orog, m) on discrete levels; ice shelves in the
    lowest (dark blue) color, ice-free ocean in light grey
  - initial surface elevation (orog) difference to the observational
    forcing (BMA3_CISM3_sm6_v3, usurf), red/blue diverging colormap;
    only for C001

One PNG per historical experiment, variable and lab:

  init-velocity_C001_CISM_NORCE.png
  init-orog_C001_CISM_NORCE.png
  init-orog-diff_C001_CISM_NORCE.png

Run with the 'plotting' environment:
  /nird/datapeak/NS11016K/miniforge3_26/envs/plotting/bin/python plot_init_velocity.py
Optional flags: --lab {NORCE,NCAR,CISM8}
"""

import os
import glob
import argparse

import numpy as np
import netCDF4 as nc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm

from plot_scalar_summary import LABS

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "Plots", "Init")

# observational bedrock topography (per grid spacing)
OBS_DIR = os.path.normpath(os.path.join(HERE, "..", "Obs"))
OBS_FILES = {16000.0: "BMA3_CISM3_sm6_v3_16000m.nc",
             8000.0: "BMA3_CISM3_sm6_v3_08000m.nc"}

# log-like levels [m/yr]
LEVELS = np.array([0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000])

# surface elevation levels [m]
OROG_LEVELS = np.arange(0, 4500, 500)

# mask threshold for fractional mask variables (sftflf)
MASK_THRESH = 0.5

# colors for the orog plot
OCEAN_COLOR = "lightgrey"   # ice-free ocean
# ice shelves are drawn in the lowest level color (dark blue for 'terrain')

# m s-1 -> m/yr
S_TO_YR = 86400 * 365.25


def find_file(data_root, exp_id, var):
    pattern = os.path.join(data_root, exp_id, f"{var}_AIS_*_{exp_id}_*.nc")
    matches = sorted(glob.glob(pattern))
    return matches[0] if matches else None


def read_first_field(data_root, exp_id, var):
    """Return (x, y, field) of the first time step of a gridded variable."""
    path = find_file(data_root, exp_id, var)
    if path is None:
        return None
    with nc.Dataset(path) as ds:
        x = np.asarray(ds.variables["x"][:])
        y = np.asarray(ds.variables["y"][:])
        field = np.ma.filled(ds.variables[var][0], np.nan)
    return x, y, field


def read_first_speed(data_root, exp_id):
    """Return (x, y, speed) of the first time step, speed in m/yr."""
    px = find_file(data_root, exp_id, "xvelmean")
    py = find_file(data_root, exp_id, "yvelmean")
    if px is None or py is None:
        return None
    with nc.Dataset(px) as ds:
        x = np.asarray(ds.variables["x"][:])
        y = np.asarray(ds.variables["y"][:])
        vx = np.ma.filled(ds.variables["xvelmean"][0], np.nan)
    with nc.Dataset(py) as ds:
        vy = np.ma.filled(ds.variables["yvelmean"][0], np.nan)
    speed = np.sqrt(vx**2 + vy**2) * S_TO_YR
    return x, y, speed


def plot_init_field(exp_id, esm, member, x, y, field, levels, cmap, label,
                    tag, title, model_lab, out_dir, bad_color=None):
    """Plot a 2D field on discrete levels (NaN cells in bad_color)."""
    fig, ax = plt.subplots(figsize=(8, 8))
    cmap = plt.get_cmap(cmap, len(levels) - 1)
    if bad_color is not None:
        cmap.set_bad(bad_color)
    norm = BoundaryNorm(levels, cmap.N, clip=True)
    im = ax.pcolormesh(x / 1000, y / 1000, field, cmap=cmap, norm=norm,
                       shading="auto", rasterized=True)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("x [km]")
    ax.set_ylabel("y [km]")
    cbar = fig.colorbar(im, ax=ax, shrink=0.8,
                        ticks=levels[::2] if len(levels) > 10 else levels)
    cbar.set_label(label)

    fig.tight_layout()
    fname = f"init-{tag}_{exp_id}_{model_lab}.png"
    out = os.path.join(out_dir, fname)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def plot_init_velocity(exp_id, esm, member, x, y, speed, model_lab, out_dir):
    """Plot speed on log-like discrete scale."""
    plot_init_field(
        exp_id, esm, member, x, y, speed, LEVELS, "viridis",
        "speed [m/yr]", "velocity",
        f"{exp_id} historical {esm}: initial velocity magnitude [m/yr]",
        model_lab, out_dir)


def read_obs_usurf(dx):
    """Return (x, y, usurf) of the observational surface elevation for
    grid spacing dx."""
    fname = OBS_FILES.get(float(dx))
    if fname is None:
        return None
    path = os.path.join(OBS_DIR, fname)
    with nc.Dataset(path) as ds:
        x = np.asarray(ds.variables["x1"][:])
        y = np.asarray(ds.variables["y1"][:])
        usurf = np.ma.filled(ds.variables["usurf"][0], np.nan)
    return x, y, usurf


def plot_orog_diff(exp_id, esm, member, x, y, orog, model_lab, out_dir):
    """Plot model minus observed surface elevation (red/blue)."""
    dx = float(x[1] - x[0])
    obs = read_obs_usurf(dx)
    if obs is None:
        print(f"WARNING: no obs file for grid spacing {dx:.0f} m, "
              f"{exp_id} skipped")
        return
    xo, yo, usurf_obs = obs
    if not (np.allclose(xo, x) and np.allclose(yo, y)):
        print(f"WARNING: obs grid does not match model grid for {exp_id}, "
              f"skipped")
        return
    diff = orog - usurf_obs
    vmax = 500.0
    print(f"  orog diff range {np.nanmin(diff):.1f} .. "
          f"{np.nanmax(diff):.1f} m (rms {np.sqrt(np.nanmean(diff**2)):.1f} m),"
          f" vmax {vmax:.0f} m")

    fig, ax = plt.subplots(figsize=(8, 8))
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("lightgrey")
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    im = ax.pcolormesh(x / 1000, y / 1000, diff, cmap=cmap, norm=norm,
                       shading="auto", rasterized=True)
    ax.set_aspect("equal")
    ax.set_title(f"{exp_id} historical {esm}: initial surface elevation "
                 "difference to observations [m]", fontsize=11)
    ax.set_xlabel("x [km]")
    ax.set_ylabel("y [km]")
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("model - obs [m]")

    fig.tight_layout()
    fname = f"init-orog-diff_{exp_id}_{model_lab}.png"
    out = os.path.join(out_dir, fname)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lab", choices=sorted(LABS), default="NORCE",
                        help="which lab's data to plot (default: NORCE)")
    args = parser.parse_args()

    lab_cfg = LABS[args.lab]
    data_root = os.path.normpath(lab_cfg["data_root"])
    model_lab = lab_cfg["model_lab"]
    experiments = lab_cfg["experiments"]
    print(f"Lab: {args.lab} ({model_lab}), data root: {data_root}")

    os.makedirs(OUT_DIR, exist_ok=True)

    for exp_id, (exp, esm, member, _, _) in experiments.items():
        if exp != "historical":
            continue
        # velocity magnitude
        res = read_first_speed(data_root, exp_id)
        if res is None:
            print(f"WARNING: no velocity files for {exp_id}, skipped")
        else:
            x, y, speed = res
            print(f"Processing {exp_id} ({esm}) ... "
                  f"max speed {np.nanmax(speed):.1f} m/yr")
            plot_init_velocity(exp_id, esm, member, x, y, speed, model_lab,
                               OUT_DIR)
        # surface elevation, with ice shelves in the lowest (dark blue)
        # color and ice-free ocean in light grey
        res = read_first_field(data_root, exp_id, "orog")
        res_fl = read_first_field(data_root, exp_id, "sftflf")
        res_ice = read_first_field(data_root, exp_id, "sftgif")
        if res is None:
            print(f"WARNING: no orog file for {exp_id}, skipped")
        else:
            x, y, orog = res
            if res_fl is not None and res_ice is not None:
                floating = res_fl[2] > MASK_THRESH
                ice = res_ice[2] > MASK_THRESH
                # ice shelves -> 0 m (lowest level bin, dark blue);
                # ice-free ocean -> NaN (light grey)
                orog = np.where(floating, OROG_LEVELS[0], orog)
                orog = np.where(ice, orog, np.nan)
            else:
                print(f"WARNING: no sftflf/sftgif for {exp_id}, "
                      "plotting plain orog")
            print(f"Processing {exp_id} ({esm}) ... "
                  f"max surface elevation {np.nanmax(orog):.1f} m")
            plot_init_field(
                exp_id, esm, member, x, y, orog, OROG_LEVELS, "terrain",
                "surface elevation [m]", "orog",
                f"{exp_id} historical {esm}: initial surface elevation [m]",
                model_lab, OUT_DIR, bad_color=OCEAN_COLOR)
        # surface elevation difference to observations (C001 only)
        if exp_id == "C001":
            res = read_first_field(data_root, exp_id, "orog")
            if res is None:
                print(f"WARNING: no orog file for {exp_id}, skipped")
            else:
                x, y, orog = res
                print(f"Processing {exp_id} ({esm}) orog difference ...")
                plot_orog_diff(exp_id, esm, member, x, y, orog, model_lab,
                               OUT_DIR)

    print("Done.")


if __name__ == "__main__":
    main()
