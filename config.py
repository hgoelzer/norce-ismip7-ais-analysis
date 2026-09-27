"""Shared scaling config for ISMIP7 NORCE AIS analysis plots.

All bounds are taken from the largest ranges found across all
models/labs (NORCE CISM 16km, NCAR CISM, NORCE CISM8 8km), rounded to
round numbers. If another model/lab exceeds a bound, update the value
here so all plots keep a common scale.
"""

# ---------------------------------------------------------------------------
# Per-experiment lithk anomaly maps (plot_scalar_exps.py)
# ---------------------------------------------------------------------------
# Guide: NCAR C007 (ssp585, CESM2-WACCM m001), max |lithk anomaly|:
#   2100:  757.9 m, 2200: 2365.6 m, 2300: 3264.2 m

# end year -> vmax [m] for lithk anomaly colorbar (symmetric: -vmax..+vmax)
LITHK_VMAX = {
    2100: 800.0,
    2200: 2400.0,
    2300: 3300.0,
}

# ---------------------------------------------------------------------------
# Scalar summary plots (plot_scalar_summary.py)
# ---------------------------------------------------------------------------
# y-axis limits per variable and end year, in PLOT units (Gt, Gt/yr,
# 10^6 km^2, mm). Each entry is (ymin, ymax); None = auto.
# Values are the extremes over all labs, rounded outward.

# absolute values: var -> {end_year: (ymin, ymax)}
SCALAR_YLIM_ABS = {
    "lim": {
        2030: (23.8e6, 24.1e6),
        2100: (23.8e6, 24.1e6),
        2200: (23.0e6, 24.1e6),
        2300: (22.3e6, 24.1e6),
    },
    "limnsw": {
        2030: (20.3e6, 20.6e6),
        2100: (20.3e6, 20.6e6),
        2200: (20.3e6, 20.6e6),
        2300: (20.0e6, 20.6e6),
    },
    "iareagr": {
        2030: (12.2, 12.3),
        2100: (12.1, 12.3),
        2200: (11.4, 12.3),
        2300: (10.7, 12.3),
    },
    "iareafl": {
        2030: (1.2, 1.5),
        2100: (1.2, 1.5),
        2200: (0.1, 1.6),
        2300: (0.0, 1.6),
    },
    "tendacabf": {
        2030: (-3000, 4200),
        2100: (-3000, 4200),
        2200: (-3000, 4200),
        2300: (-3000, 4200),
    },
    "tendlibmassbfgr": {
        2030: (-4000, 0),
        2100: (-17000, 0),
        2200: (-30000, 0),
        2300: (-30000, 0),
    },
    "tendlibmassbffl": {
        2030: (0, 2200),
        2100: (0, 9500),
        2200: (0, 19000),
        2300: (0, 23000),
    },
    "tendlicalvf": {
        2030: (-2400, 0),
        2100: (-2400, 0),
        2200: (-2800, 0),
        2300: (-3100, 0),
    },
}

# anomalies rel. to end of historical: var -> {end_year: (ymin, ymax)}
SCALAR_YLIM_ANOM = {
    "lim": {
        2030: (-11000, 5000),
        2100: (-170000, 5000),
        2200: (-1.1e6, 5000),
        2300: (-1.8e6, 10000),
    },
    "limnsw": {
        2030: (-4500, 5000),
        2100: (-5000, 25000),
        2200: (-110000, 30000),
        2300: (-480000, 30000),
    },
    "iareagr": {
        2030: (-0.02, 0.02),
        2100: (-0.1, 0.02),
        2200: (-0.8, 0.02),
        2300: (-1.5, 0.02),
    },
    "iareafl": {
        2030: (-0.02, 0.02),
        2100: (-0.03, 0.02),
        2200: (-1.3, 0.2),
        2300: (-1.4, 0.2),
    },
    "tendacabf": {
        2030: (-1500, 1100),
        2100: (-1500, 1100),
        2200: (-6200, 1100),
        2300: (-6200, 1100),
    },
    "tendlibmassbfgr": {
        2030: (-1600, 1200),
        2100: (-15000, 1200),
        2200: (-29000, 1200),
        2300: (-29000, 1200),
    },
    "tendlibmassbffl": {
        2030: (-700, 900),
        2100: (-700, 8500),
        2200: (-700, 18000),
        2300: (-700, 22000),
    },
    "tendlicalvf": {
        2030: (-150, 900),
        2100: (-150, 900),
        2200: (-1200, 1100),
        2300: (-1300, 1400),
    },
}

# sea-level contribution [mm]: end_year -> (ymin, ymax)
# Guide: NORCE C007 (ssp585) upper bound; CISM8 C007 (ssp585, CESM2,
# 2015-2224) dips to -43 mm at 2100 and reaches +93 mm at 2224.
SLC_YLIM = {
    2030: (-15, 15),
    2100: (-50, 70),
    2200: (-300, 100),
    2300: (-1350, 100),
}
