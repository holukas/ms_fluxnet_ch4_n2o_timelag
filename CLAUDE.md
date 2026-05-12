# CLAUDE.md

This file documents the codebase for time lag detection in eddy covariance flux tower data.

## What This Project Does

**ms_fluxnet_ch4_n2o_timelag** detects time lags in eddy covariance measurements of CH4 (methane) and N2O (nitrous oxide). The challenge: closed-path gas analyzers introduce a physical delay as air samples travel through tubing. Incorrect lags lead to systematically biased flux estimates. The tool uses RFlux, which implements a recent method from Vitale et al. (2024) called pre-whitening with block-bootstrap (PWB) — a significant improvement over traditional cross-correlation methods, especially when fluxes are weak.

### Why This Matters

Standard cross-correlation works fine when correlations are strong. But CH4 and N2O at CH-Cha site correlate near −0.07 to −0.13 — weak enough that multiple spurious peaks appear, all of similar height, confusing standard lag-detection. The PWB approach filters out autocorrelation from the time series first (pre-whitening), then uses 99 bootstrap resamples to quantify uncertainty. The result: a detected lag plus confidence bounds indicating whether to trust it.

### Sources

- **Vitale et al. (2024)**: "A pre-whitening with block-bootstrap cross-correlation procedure for temporal alignment of data sampled by eddy covariance systems." *Environmental and Ecological Statistics* 31:219–244. The paper benchmarks PWB on CH-Cha data (Chamau, Switzerland).
- **RFlux manual** (v3.2.0, Vitale et al. 2022): See `info/RFlux-manual.pdf`
- **RFlux on GitHub**: https://github.com/icos-etc/RFlux (maintained by ICOS-ETC)

---

## Layout & Flow

The project flows in two languages:

**R** (the analysis): Raw flux data → RFlux's PWB algorithm → CSV output with detected lags and uncertainty bounds.

**Python** (the visualization): Those CSVs → comparison plots against EddyPro (the industry standard software) → PDF/PNG.

Here's the directory tree:

```
.
├── ms_fluxnet_ch4_n2o_timelag.R    # Main script. Run this first.
├── comparison_python_plot/
│   ├── plot.py                      # Reads R output, makes plots
│   ├── tlag_results_partN.csv      # R generates these (7 parts, Sept 2021 data)
│   └── *.csv                        # EddyPro reference files for comparison
├── input/01-raw_data_ascii_done/   # Your raw data
│   └── CH-CHA_*.csv.gz             # ~4,000 compressed files, Sept 2021
├── output/                          # Where results go (git-ignored)
├── tlag_plots/                     # Generated plots
├── info/                           # Reference docs
│   ├── RFlux-manual.pdf
│   └── Vitale et al. - 2024 - ...pdf
├── renv.lock & .venv               # Dependency snapshots
└── info_tlag_results.md            # Column definitions for output
```

---

## Getting the Environments Set Up

### R

Uses **renv** for isolation. Key packages: RFlux, forecast, xts, and (archived) egcm. RFlux requires R ≥ 4.0. On first run, `renv::init()` installs everything pinned to `renv.lock`. The script handles this automatically.

### Python

Uses **uv** for speed. Python ≥ 3.9, with pandas, matplotlib, scipy, numpy for data and plotting. Config is in `pyproject.toml`.

---

## Running the Analyses

### The R Script

```bash
# In R or RStudio:
source("ms_fluxnet_ch4_n2o_timelag.R")
```

What happens:
1. Installs RFlux and dependencies (via renv) if not already present
2. Reads every `CH-CHA_*.csv.gz` file from `input/01-raw_data_ascii_done/`
3. Runs PWB time lag detection on CH4 and N2O separately
4. Writes output CSV files to `output/` (or the configured output directory)

**If it crashes mid-run**: Check the last completed file timestamp. The script's designed to resume from there, not restart from file 1. RFlux logs this internally.

**Output columns** are documented in `info_tlag_results.md`. Key ones: `ch4_tlag_sec` (raw detected lag), `ch4_hdi_range_sec` (uncertainty), `ch4_pwbopt_sec` (final recommended lag after S1/S2/S3 logic).

### The Python Plots

```bash
# First time:
uv sync

# Then run:
uv run python comparison_python_plot/plot.py
```

Generates `timelags.pdf` and `timelags.png` — side-by-side comparisons of PWB lags vs. EddyPro's output.

### Linting & Formatting

```bash
ruff check comparison_python_plot/
black --line-length 100 comparison_python_plot/
```

---

## Understanding the Output

### What Each Column Means

**Raw detection:**
- `ch4_tlag_sec` / `n2o_tlag_sec` — Lag detected for that 30-min file (in seconds)
- `ch4_hdi_lci_sec` / `ch4_hdi_uci_sec` — Lower and upper bounds of the 95% confidence interval from bootstrap
- `ch4_hdi_range_sec` — Width of that interval (0 = all bootstrap samples agreed; >6 s = disagree wildly)
- `ch4_cor` / `n2o_cor` — Cross-correlation magnitude at the peak (usually −0.07 to −0.13 for low fluxes; that's *expected*)

**Final recommended lag (use this for flux calculations):**
- `ch4_pwbopt_sec` / `n2o_pwbopt_sec` — The lag after S1/S2/S3 logic (explained below)
- `ch4_flag` / `n2o_flag` — How that lag was chosen (S1, S2, or S3)

### The S1/S2/S3 Logic

PWB assigns a reliability flag to each 30-min period:

**S1_optimal**: The detected lag has an HDI range < 0.5 s. Inherently reliable; use it.

**S2_optimal**: HDI range ≥ 0.5 s, but the lag is within 0.5 s of the previous S1/S2 optimal. Accept it for continuity — short-term stability (30–60 min windows) justifies this. (If CH4 was 1.8 s at 10:00, and 1.9 s at 10:30 with a wide HDI, take it; don't jump to −2 s.)

**S3_unreliable**: The detection is too uncertain (wide HDI) or too far from the prior. Fall back to the last trusted lag. This avoids the mirroring effect: spurious negative lags (−4.85 s when expected ~1.8 s) that standard methods get fooled by on weak signals.

Bottom line: **Always use `*pwbopt_sec`, not the raw `*tlag_sec`.**

### Interpreting the Comparison Plots

`timelags.pdf` and `timelags.png` show detected lags (PWB) vs. EddyPro's lags (the baseline), plotted across time.

Expected patterns:
- Both tracks (PWB and EddyPro) cluster around 1.7–2.0 s for most of the day.
- They often track each other closely.
- When PWB shows a S3_unreliable flag (wide HDI), it typically *doesn't* jump dramatically; instead, it stays near the prior good value. EddyPro might jump wildly on weak signal days; PWB doesn't, which is the point.
- If PWB and EddyPro diverge significantly, check the corresponding HDI range. A wide HDI (>1 s) indicates PWB uncertainty; narrow HDI indicates PWB confidence and suggests EddyPro might be chasing noise.

Red flags:
- If PWB lags cluster at implausible values (e.g., always −3 s), something's wrong with the input data or tube setup description.
- If every single day has S3 flags, the signal is too weak; consider whether the analyzer settings or tube was changed during the campaign.

---

## How the Algorithm Works

PWB has three steps:

1. **Pre-whiten**: Fit an autoregressive (AR) model to each time series, then filter it out. This removes autocorrelation that would otherwise inflate noise in the cross-correlation.
2. **Bootstrap**: Resample the pre-whitened data 99 times in ~20-timestep blocks (keeping short-timescale structure intact). For each resample, compute the cross-correlation and find the peak lag.
3. **Synthesize**: Collect the 99 lag estimates. The most common one is the detected lag. The range of the 95% highest-density interval tells you the uncertainty.

The S1/S2/S3 decision tree then handles edge cases: narrow intervals are accepted, wide but stable intervals are accepted for continuity, and unstable intervals fall back to prior values.

### RFlux Functions You'll See

- `tlag_detection()` — Does the PWB calculation. Takes scalar (CH4/N2O), sonic temperature, vertical wind, sampling frequency, bootstrap count, plot flag.
- `despiking()` — Removes outliers (spikes) from time series.
- `cleanFlux()` — A full quality-control pipeline (useful reference when modifying preprocessing).

---

## Practical Notes for Edits

1. **Data format**: Files are gzip-compressed CSV. Columns: U, V, W (wind), T_SONIC (sonic temperature), CO2, H2O (for reference), plus instrument diagnostics. Decompress with `gunzip CH-CHA_202109010100.csv.gz` to inspect the structure.

2. **Bootstrap replicates**: The script uses `Rboot = 3` (PWB minimum). The paper used 99 for the benchmark. During debugging, 3 runs fast; increase it for tighter confidence intervals.

3. **N2O vs. CH4 behavior**: N2O has stronger, more stable signal. Most N2O periods are S1_optimal with HDI < 0.15 s. CH4 is noisier. Seeing several S3_unreliable flags for CH4 is normal and correct.

4. **Afternoon dropouts**: Many 30-min files have NA in the output. That's the `valid_enough()` guard: if the series is constant or mostly missing (instrument drift, calibration, outage), skip it. Don't force processing through those gaps.

5. **Physical reality check**: At CH-Cha, expected lags are 1.7–2.0 s (tube length, flow rate, wall absorption). If detected lags stray far (e.g., consistently −2 s) and have tight HDI, suspect a problem with the analyzer setup or tube documentation, not the algorithm.

---

## What to Read Next

- `info_tlag_results.md` — Deep dive on output columns
- `info/RFlux-manual.pdf` — Function signatures and examples
- `info/Vitale et al. - 2024 - ...pdf` — The peer-reviewed paper (equations, validation, limitations)
