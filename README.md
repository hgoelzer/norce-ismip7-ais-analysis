# ISMIP7 AIS Analysis Tools

This repository creates diagnostic plots from the ISMIP7-compliant Antarctic
ice-sheet NetCDF files produced by
[`../norce-ismip7-ais-processing/`](../norce-ismip7-ais-processing/). It does
not run the ice-sheet model or alter the input data. The scripts read the
processed files and write PNG figures under `Plots/`.

## Setup

Run the scripts from this directory with the `plotting` Python environment,
which provides `matplotlib`, `netCDF4`, `numpy`, and `cftime`:

```bash
PYTHON=/nird/datapeak/NS11016K/miniforge3_26/envs/plotting/bin/python
"$PYTHON" plot_scalar_summary.py --help
```

By default, input data are found relative to this repository:

- NORCE 16 km CISM: `../AIS/NORCE/CISM/CORE/`
- NCAR CISM: `../AIS/NCAR/CISM/CORE/`
- NORCE 8 km CISM8: `../AIS/NORCE/CISM8/CORE/`
- Observational fields for the initial-state maps and model comparison:
  `../Obs/BMA3_CISM3_sm6_v3_{16000,08000}m.nc`

Each data root should contain experiment directories such as `C001/`, with
processed variable files named like `{var}_AIS_*_{experiment}_*.nc`. The
analysis scripts expect ISMIP7 variables including scalar `lim`, `limnsw`,
`iareagr`, `iareafl`, `tendacabf`, `tendlibmassbfgr`, `tendlibmassbffl`, and
`tendlicalvf`; gridded `lithk`, `sftgrf`, `sftflf`, `sftgif`, `orog`,
`xvelmean`, and `yvelmean` are used where available. Missing files or fields
are reported as warnings and the affected plots are skipped.

## Run the tools

All commands below are run from this directory. Use `--lab {NORCE,NCAR}` to
select the institution and `--model {CISM,CISM8}` to select the model. Both
default to `NORCE` and `CISM`; supported combinations are NORCE/CISM,
NCAR/CISM, and NORCE/CISM8. CISM8 is currently available only for NORCE.

```bash
# Generate all three plot families for one dataset
"$PYTHON" run_all_plots.py --lab NORCE --model CISM8

# Run all plot families for NCAR/CISM; restrict only the experiment maps
"$PYTHON" run_all_plots.py --lab NCAR --exp C007

# Scalar time series, anomalies, and sea-level contribution
"$PYTHON" plot_scalar_summary.py

# Only absolute-value or anomaly time series (sea-level plots are also written)
"$PYTHON" plot_scalar_summary.py --abs-only
"$PYTHON" plot_scalar_summary.py --anom-only --lab NORCE --model CISM8

# Per-experiment spatial changes; repeat --exp to select multiple experiments
"$PYTHON" plot_scalar_exps.py --lab NORCE
"$PYTHON" plot_scalar_exps.py --lab NCAR --exp C007
"$PYTHON" plot_scalar_exps.py --lab NORCE --model CISM8 --exp C007 --exp C008

# Initial historical-state maps
"$PYTHON" plot_initial.py --lab NORCE

# Plot the observation fields once at both available resolutions
"$PYTHON" plot_initial.py --observations
```

`run_all_plots.py` runs the scalar summary, experiment maps, and initial-state
maps sequentially with the selected lab/model. It stops if any script fails.
Its `--abs-only` and `--anom-only` options affect the scalar summary only;
repeat `--exp ID` to limit the experiment maps without skipping the other
plot families.

The scripts create output folders as needed:

- `Plots/Summary/{model_lab}/` for scalar time-series figures.
- `Plots/Exps/{model_lab}/` for per-experiment map figures.
- `Plots/Init/{model_lab}/` for model-specific initial-state figures.
- Model-independent observation figures remain in `Plots/Init/`.

## Expected figures

### Scalar summaries

`plot_scalar_summary.py` joins each scenario to its historical run for the same
ESM member. The available experiment sets are:

| Lab / model | Configured experiments |
| --- | --- |
| NORCE / CISM | C001-C011: historical, ssp370, ssp126, ssp585, ctrl, and ocx |
| NCAR / CISM | C001 and C007 only: historical and ssp585 |
| NORCE / CISM8 | C001-C011: historical, ssp370, ssp126, ssp585, ctrl, and ocx |

For each available variable, the default run writes absolute and anomaly plots
for end years 2030, 2100, 2200, and 2300. Anomalies use the last historical
value (end of 2014) for the same member; ocx has no historical counterpart and
uses its own value interpolated at 2015. Scenario lines are joined to their
historical series where one is configured. Data are truncated at each
requested end year, so a shorter run only shows its available period.

Filenames follow `{variable}[-anom]_{model_lab}_{end_year}.png`, for example
`limnsw-anom_CISM_NORCE_2100.png`. A further
`sea-level_{model_lab}_{end_year}.png` figure shows sea-level
contribution derived from the limnsw anomaly: a loss of 361.8 Gt corresponds
to a 1 mm rise. There are eight scalar variables; with complete input data,
the default NORCE or CISM8 run produces up to 68 summary PNGs, while NCAR's
two configured experiments produce fewer series in the same figure set.
The `--abs-only` and `--anom-only` flags select the corresponding scalar
series; sea-level contribution figures are still generated.

Scalar mass is converted from kg to Gt, flux from kg/s to Gt/year, and area
from m2 to million km2. Snapshot fields are shown as ST-type values and
annual-mean flux fields as FL-type values. Sea-level contribution is reported
in mm.

### Per-experiment maps

`plot_scalar_exps.py` compares gridded ice thickness and masks against the
last historical field from the same member. For each configured non-historical
experiment with the required files and reference, it writes maps for 2100,
2200, and 2300:

- Ice-thickness anomaly (`lithk`), in metres relative to the end of
  historical.
- Grounded, floating, and total ice-mask changes (`sftgrf`, `sftflf`,
  `sftgif`). Fractional masks use a 0.5 threshold; red indicates lost mask
  coverage and green indicates gained coverage.
- Surface elevation (`orog`) snapshots for model years 2015, 2030, 2100, 2200,
  and 2300, using the same terrain levels and shelf/ocean masking as the
  initial-state maps.

Names follow `lithk-anom_{experiment}_{model_lab}_{end_year}.png` and
`{grounded-mask|floating-mask|ice-mask}-anom_{experiment}_{model_lab}_{end_year}.png`,
for example `lithk-anom_C007_CISM_NORCE_2300.png`. Orography snapshots are
named `orog_{experiment}_{model_lab}_{year}.png`. If an experiment ends before
a requested snapshot year, its latest available field is used and the actual
sample year is shown in the title. A complete NORCE or CISM8 dataset can
produce up to 153 map PNGs; NCAR's single configured scenario produces up to
17. The plotted field is the latest available time step at or before the
named year.

### Initial-state and observation maps

`plot_initial.py` uses the first time step of C001 only; C002 is omitted
because its initial fields are nearly identical. It produces velocity
magnitude in m/year and surface elevation in metres. Velocity is converted
from m/s. Ice shelves are shown in the lowest elevation color and ice-free
ocean is grey. It also compares model surface elevation with observational
`usurf` and plots model-minus-observation in metres. The difference is
meaningful only when the observation file for the model grid is available and
the grids match.

Output names look like `init-velocity_C001_CISM_NORCE.png`,
`init-orog_C001_CISM_NORCE.png`, and
`init-orog-diff_C001_CISM_NORCE.png`. Each configured lab/model produces up to
three model-specific initial-state figures for C001.

Run `plot_initial.py --observations` once to plot the six observed fields
(`topg`, `usurf`, `thk`, `lsurf`, `grounded_mask`, and `floating_mask`) from
both the 16 km and 8 km BMA3 files. This mode does not read a model dataset, so
it does not need to be repeated when adding another lab or model. It writes
`obs-{variable}_{16|8}km.png` under `Plots/Init/`; continuous fields use shared
color scales across the two grids and masks use categorical 0/1 colors. It also
writes `init-orog_OBS_16km.png` and `init-orog_OBS_8km.png`, using the same
discrete terrain levels and floating-shelf/ice-free-ocean treatment as the
model `init-orog` plots.

## Reading results and scaling

- Scalar anomaly curves are relative to the end of historical, not to the
  first year of each scenario. Negative limnsw anomalies correspond to
  positive sea-level contribution.
- In the spatial change maps, ice-thickness change is signed (loss or gain);
  mask maps show where a thresholded category is lost or gained.
- Shared map and anomaly-axis limits are defined in `config.py` so labs can be
  compared on consistent scales. Sea-level plot limits cover the measured
  extrema across configured datasets. Absolute scalar plots auto-scale. The map
  script warns if an ice-thickness anomaly exceeds the configured color range;
  review `LITHK_VMAX` and regenerate maps if needed. Anomaly axis limits are
  fixed for comparison, so values outside them may be clipped.
- Output PNGs are generated artifacts and are excluded from version control.

For the processing workflow that creates the input NetCDF files, see the
[`norce-ismip7-ais-processing` README](../norce-ismip7-ais-processing/README.md).