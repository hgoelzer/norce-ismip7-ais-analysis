# Copilot Instructions — NORCE ISMIP7 AIS Analysis

Plotting and analysis of the ISMIP7-compliant NetCDF output produced by
`../norce-ismip7-ais-processing/` (see its `copilot-instructions.md` for the
processing side). Three plotting scripts, one all-plots wrapper, and one
shared scaling config.

## Layout

- `plot_scalar_summary.py` — approach A: all experiments of a lab/model in one
  figure per scalar variable (absolute + anomaly + sea-level contribution).
  Holds the `LABS` dict (per-lab data root, model_lab tag, experiment table).
- `plot_scalar_exps.py` — approach B: per-experiment 2D maps of lithk and
  grounded/floating/ice-mask changes vs end of historical, plus orography
  snapshots for model years 2015, 2030, 2100, 2200, and 2300.
  Has its own copy of `LABS` (keep in sync with `plot_scalar_summary.py`).
- `plot_initial.py` — initial-state maps for C001 (first historical time step;
  C002 omitted because its initial fields are nearly identical):
  velocity magnitude, surface elevation (with ice shelves dark blue, ice-free
  ocean light grey), and surface-elevation difference to the observational
  forcing (`../Obs/BMA3_CISM3_sm6_v3_{16000,08000}m.nc`, C001 only).
  `--observations` plots all six BMA3 fields on both grids once, independently
  of the selected lab/model, including `init-orog_OBS_{16,8}km.png` in the
  same discrete terrain style as model orography plots.
- `run_all_plots.py` — runs the three plotting scripts sequentially for one
  selected lab/model combination; accepts repeatable `--exp` filters for maps.
- `config.py` — shared scaling tables ONLY (`LITHK_VMAX`, `SCALAR_YLIM_ANOM`,
  `SLC_YLIM`). **No `LABS` here** — import `LABS` from `plot_scalar_summary`.
- Output dirs: `Plots/Summary/`, `Plots/Exps/`, `Plots/Init/`.

## Environments & commands

- Plotting python: `/nird/datapeak/NS11016K/miniforge3_26/envs/plotting/bin/python`
  (matplotlib, netCDF4, numpy, cftime). The `nc` env used for processing has
  **no matplotlib**.
- Select datasets with `--lab NORCE|NCAR` and `--model CISM|CISM8` (both
  default to NORCE/CISM); supported pairs are NORCE/CISM, NCAR/CISM, and
  NORCE/CISM8.
- Run: `<plotting python> plot_scalar_summary.py [--lab ...] [--model ...]
  [--abs-only|--anom-only]`, `... plot_scalar_exps.py [--lab ...] [--model ...]
  [--exp C008]`, `... plot_initial.py [--lab ...] [--model ...]`.
- Generate all plot families in one call: `<plotting python> run_all_plots.py
  [--lab ...] [--model ...] [--abs-only|--anom-only] [--exp C008]`.

## Data conventions

- Scalar files: one float32 data var on dim `time`; time = days since
  1850-01-01, calendar standard. ST-type (lim, limnsw, iareagr, iareafl)
  sampled Jan 1; FL-type (tend* fluxes) annual means at Jul 1 with `time_bnds`.
- Gridded files: dims (time, y, x); 16 km grid 381×381, x/y −3040000..3040000 m
  (CISM8: 761×761 at 8 km); crs epsg:3031. Years = origin + days/365.25.
- Velocity: `xvelmean`/`yvelmean` in m s⁻¹ → m/yr via ×86400×365.25.
- Masks: `sftgrf`/`sftflf` are fractional 0–1 (threshold at 0.5); `sftgif` is
  binary 0/1.
- Unit conversions: kg→Gt ×1e-12; kg s⁻¹→Gt/yr ×1e-12×86400×365.25;
  m²→10⁶ km² ×1e-12; sea level: −361.8 Gt limnsw = +1 mm SLC.
- Experiments: C001/C002 historical (CESM2-WACCM / MRI-ESM2-0, 1970–2014),
  C003/C004 ssp370 (→2100), C005/C006 ssp126 (→2300), C007/C008 ssp585
  (→2300), C009/C010 ctrl (→2300), C011 ocx (ERA, 1990–2025). NCAR lab has
  only C001 + C007 (2000–2014 / 2015–2300). CISM8 has all 11 (C007 ends
  2224, the rest reach 2300).

## Conventions & pitfalls

- Anomaly reference = last historical time step (end 2014) of the same ESM
  member; ocx has no historical predecessor → reference interpolated at 2015
  from its own series.
- Global y-axis scaling (from `config.py`) applies ONLY to anomaly plots;
  absolute plots auto-scale.
- `LITHK_VMAX` and `SLC_YLIM` in `config.py` are guided by the largest
  ranges across all labs (LITHK: NCAR C007; SLC: NORCE C007 high end,
  CISM8 C007 low end — it gains mass first, SLC −43 mm at 2100, then
  +93 mm by 2224). The scripts print a WARNING when a field exceeds a
  bound — update `config.py` then and regenerate the affected plots for
  ALL labs so the common scale holds.
- File naming: `{var}[-anom]_{model_lab}_{endyear}.png`; titles carry
  `{exp_id} {exp} {esm}: ...` and no "rel. to end of historical" text.
- An experiment may end before 2300 (e.g. CISM8 C007 ends 2224) — the
  scripts plot the latest available year ≤ end_year.
- `plt.get_cmap(name, N).colors` fails for LinearSegmentedColormap — use
  `plt.get_cmap(name, N)` directly with `BoundaryNorm`.
- `cp` is aliased interactive on the login node — use `\cp -f`.

## Workflow conventions

- User drives step-by-step; confirm before regenerating large plot batches.
- This folder is a git repo (remote
  `git@github.com:hgoelzer/norce-ismip7-ais-analysis.git`, branch `main`).
  PNG outputs are git-ignored; commit scripts + `config.py` changes.
- Commit messages: concise imperative summary + bullet body of the
  substantive changes.
