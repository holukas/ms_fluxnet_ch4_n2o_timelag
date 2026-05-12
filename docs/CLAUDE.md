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
├── scripts/
│   ├── 01_tlag_detection_pwb.R          # Time lag detection (PWB + S1/S2/S3)
│   ├── 02_tlag_compare_pwbopt_strategies.R  # Compare standard vs pre-filtered approach
│   └── visualization/
│       ├── plot.py                      # EddyPro vs PWB comparison
│       └── plot_comparison_strategies.py # EddyPro vs PWB pre-filtered
├── input/
│   ├── 01-raw_data_ascii/               # Original ASCII data files
│   ├── 02-run_eddypro_for_rotation/     # EddyPro config and processing notes
│   ├── 03-rotated_data_from_eddypro_level5_part[1-7]/  # Rotation-corrected data
│   │   └── *.txt                        # High-frequency wind + scalar (20 Hz)
│   └── README.md                        # Input data documentation
├── output/
│   ├── tlag_results_part[1-7].csv       # Script 1 output
│   ├── tlag_results_prefiltered_all.csv # Script 2 output (combined)
│   ├── tlag_results_prefiltered_part[1-7].csv  # Script 2 output (per-part)
│   └── tlag_plots/                      # Diagnostic plots (PDF)
├── references/                          # Reference materials
│   ├── RFlux-manual.pdf
│   ├── Vitale_et_al_2024.pdf
│   └── Rflux_source_code/               # RFlux source (v3.1.0)
├── docs/
│   ├── CLAUDE.md                        # This file
│   ├── README.md                        # Project overview
│   └── WORKFLOW.md                      # Workflow diagrams
├── renv.lock                            # R dependency snapshot
└── pyproject.toml                       # Python dependencies (uv)
```

---

## Getting the Environments Set Up

### R

Uses **renv** for isolation. Key packages: RFlux, forecast, xts, and (archived) egcm. RFlux requires R ≥ 4.0. On first run, `renv::init()` installs everything pinned to `renv.lock`. The script handles this automatically.

### Python

Uses **uv** for speed. Python ≥ 3.9, with pandas, matplotlib, scipy, numpy for data and plotting. Config is in `pyproject.toml`.

---

## Running the Analyses

### Step 1: Time Lag Detection (01_tlag_detection_pwb.R)

```bash
# In R or RStudio (from project root):
source("scripts/01_tlag_detection_pwb.R")
```

This script:
1. Loads rotation-corrected wind and scalar data from `03-rotated_data_from_eddypro_level5_partN/*.txt`
2. For each 30-min period:
   - Runs RFlux `tlag_detection()` with 99 bootstrap replicates (PWB method)
   - Extracts: raw lag, HDI bounds (lci, uci), correlation magnitude
3. Applies **S1/S2/S3 selection logic** (Vitale et al. Sect. 2.3):
   - **S1 (optimal)**: HDI < 0.5 s — low uncertainty, accept
   - **S2 (optimal)**: HDI ≥ 0.5 s BUT within 0.5 s of prior — accept for continuity
   - **S3 (unreliable)**: High uncertainty or far from prior — use last good lag
4. Generates 6-panel diagnostic plots (PDF) for each gas/period
5. Outputs CSV with lags, HDI, correlations, and PWBOPT flags

**Configuration**: Change `part` variable (line 20) to process different dataset parts (1–7).

**Resume capability**: If the script crashes, run it again. It resumes from the last completed file using a checkpoint system.

**Output columns** documented in `info_tlag_results.md`. Key ones:
- `ch4_tlag_sec` — raw detected lag
- `ch4_hdi_range_sec` — uncertainty (95% HDI width)
- `ch4_pwbopt_sec` — final recommended lag (after S1/S2/S3)
- `ch4_flag` — reliability flag (S1_optimal, S2_optimal, or S3_unreliable)

### Step 2 (Optional): Test Pre-filtering Strategy (02_tlag_compare_pwbopt_strategies.R)

```bash
# In R or RStudio (from project root, after running script 1):
source("scripts/02_tlag_compare_pwbopt_strategies.R")
```

This script tests whether pre-filtering high-uncertainty lags before applying S1/S2/S3 improves reliability. The approach:

1. Reads all parts from `output/tlag_results_part[1-7].csv` (from script 1)
2. Combines parts to maintain temporal continuity
3. Applies S1/S2/S3 twice:
   - **Standard**: All detected lags → S1/S2/S3
   - **Pre-filtered**: Remove lags with HDI > threshold (default 1.0 s) → S1/S2/S3
4. Outputs combined results with both approaches for comparison

**Output columns** (example for CH4):
- `ch4_pwbopt_std`, `ch4_flag_std` — Standard approach (all lags)
- `ch4_pwbopt_prefilter`, `ch4_flag_prefilter` — Pre-filtered approach
- Same for N2O

**Configuration** (edit at top of script):
- `hdi_prefilter_sec` (line 13): Pre-filter threshold in seconds (default 1.0)
- `parts_to_process` (line 12): Which parts to process (default 1:7)

**Outputs**:
- `output/tlag_results_prefiltered_all.csv` — Combined results (all parts, both approaches)
- `output/tlag_results_prefiltered_part[1-7].csv` — Per-part results for reference
- Console summary: % reliable for each approach, # filtered out per gas

### Step 3: Visualization

#### Option A: Compare PWB vs. EddyPro (Original Plot)

```bash
cd scripts/visualization
uv sync
uv run python plot.py
```

Generates `timelags.png` and `timelags.pdf`:
- 2×2 grid: CH4 & N2O × CM (EddyPro) & PWB
- Scatter plots with jitter + KDE distributions
- Mode lines (black dashed) for each method
- Event overlays (tillage, fertilization, precipitation)

**Configuration** (edit at top of `plot.py`):
- `FILE1["path"]`: EddyPro CSV path
- `FILE2["paths"]`: List of PWB result CSVs
- `SAVE_PATH`: Output filename
- `JITTER`, `YLIM`, `DPI`: Plot styling
- Event dates: `TILLAGE_DATE`, `FERTILIZATION_DATE`, `PRECIPITATION_DATE`

#### Option B: Compare EddyPro vs. PWB Pre-filtered

```bash
cd scripts/visualization
uv run python plot_comparison_strategies.py
```

Generates `timelags_strategies_comparison.png`:
- 2×2 grid: CH4 & N2O × EddyPro & PWB pre-filtered
- Same layout as plot.py (scatter + KDE)
- Color-coded: blue (EddyPro), orange (PWB pre-filtered)
- Event overlays + mode lines

**Configuration** (edit at top of `plot_comparison_strategies.py`):
- `FILE_EDDYPRO["path"]`: EddyPro CSV path
- `FILE_PREFILTERED["paths"]`: Combined prefiltered result CSV
- `SAVE_PATH`: Output filename
- `JITTER`, `YLIM`, `DPI`: Plot styling
- Event dates (same as above)

### Code Quality

```bash
cd scripts/visualization
ruff check .
black --line-length 100 .
```

---

## Workflow & File Dependencies

```
                 ┌─────────────────────────────────────┐
                 │ 03-rotated_data_from_eddypro_level5 │
                 │ (EddyPro L5 rotation-corrected)      │
                 └──────────────────┬──────────────────┘
                                    │
                    ┌───────────────▼────────────────┐
                    │ 01_tlag_detection_pwb.R        │
                    │ (RFlux PWB + S1/S2/S3)         │
                    └───────────────┬────────────────┘
                                    │
                    ┌───────────────▼────────────────┐
                    │ output/tlag_results_part*.csv  │
                    │ (Raw lags + PWBOPT flags)      │
                    └───────────────┬────────────────┘
                                    │
                ┌───────────────────┬───────────────────┐
                │                   │                   │
        ┌───────▼──────────┐  ┌──────▼──────────┐  ┌────▼────────────┐
        │ plot.py          │  │ 02_tlag_compare │  │ plot_comparison │
        │ (PWB vs Eddypro) │  │ _pwbopt_*.R     │  │ _strategies.py  │
        │                  │  │ (Test pre-filter)   │ (EddyPro vs PWB)│
        └───────┬──────────┘  └────┬─────────────┘  └────┬────────────┘
                │                  │                     │
        ┌───────▼──────────┐  ┌─────▼──────────────┐  ┌──▼─────────────┐
        │ timelags.png     │  │ output/tlag_results│  │ timelags_      │
        │ timelags.pdf     │  │ _prefiltered_all   │  │ strategies_    │
        │ (PWB vs EddyPro) │  │ .csv (combined)    │  │ comparison.png │
        │                  │  │ + per-part CSVs    │  │ (EddyPro vs    │
        │                  │  │                    │  │  PWB pre-filt) │
        └──────────────────┘  └────────────────────┘  └────────────────┘
```

**Read-only inputs** (git-tracked):
- Rotation-corrected wind/scalar data in `03-rotated_data_from_eddypro_level5_partN/`

**Generated files** (git-ignored):
- `output/tlag_results_part[1-7].csv` — Script 1 output (7 parts)
- `output/tlag_results_prefiltered_all.csv` — Script 2 output (combined, both approaches)
- `output/tlag_results_prefiltered_part[1-7].csv` — Script 2 output (per-part, reference)
- `tlag_plots/*.pdf` — Diagnostic plots from script 1
- `timelags.png`, `timelags.pdf` — Comparison plot (Option A: PWB vs EddyPro)
- `timelags_strategies_comparison.png` — Comparison plot (Option B: EddyPro vs PWB pre-filtered)

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

### Prefiltered Output Columns

Script 2 outputs `tlag_results_prefiltered_all.csv` (and per-part variants) with columns for both approaches:

**Standard approach** (all detected lags, S1/S2/S3 applied directly):
- `ch4_pwbopt_std` / `n2o_pwbopt_std` — Recommended lag (baseline)
- `ch4_flag_std` / `n2o_flag_std` — Flag (S1_optimal, S2_optimal, S3_unreliable)

**Pre-filtered approach** (lags with HDI > 1.0 s removed, then S1/S2/S3):
- `ch4_pwbopt_prefilter` / `n2o_pwbopt_prefilter` — Recommended lag (filtered variant)
- `ch4_flag_prefilter` / `n2o_flag_prefilter` — Flag (S1_optimal, S2_optimal, S3_unreliable)

**Original detection results** (shared by both approaches):
- `ch4_tlag_sec`, `ch4_hdi_lci_sec`, `ch4_hdi_uci_sec`, `ch4_hdi_range_sec`, `ch4_cor`
- `n2o_tlag_sec`, `n2o_hdi_lci_sec`, `n2o_hdi_uci_sec`, `n2o_hdi_range_sec`, `n2o_cor`
- `file_name` — Input file name
- `timestamp` (implicit in file_name, extractable with regex: `\d{8}-\d{4}`)

### The S1/S2/S3 Logic

PWB assigns a reliability flag to each 30-min period:

**S1_optimal**: The detected lag has an HDI range < 0.5 s. Inherently reliable; use it.

**S2_optimal**: HDI range ≥ 0.5 s, but the lag is within 0.5 s of the previous S1/S2 optimal. Accept it for continuity — short-term stability (30–60 min windows) justifies this. (If CH4 was 1.8 s at 10:00, and 1.9 s at 10:30 with a wide HDI, take it; don't jump to −2 s.)

**S3_unreliable**: The detection is too uncertain (wide HDI) or too far from the prior. Fall back to the last trusted lag. This avoids the mirroring effect: spurious negative lags (−4.85 s when expected ~1.8 s) that standard methods get fooled by on weak signals.

Bottom line: **Always use `*pwbopt_sec`, not the raw `*tlag_sec`.**

### Interpreting the Comparison Plots

#### EddyPro vs. PWB (plot.py)

`timelags.pdf` and `timelags.png` show detected lags (PWB) vs. EddyPro's lags (the baseline), plotted across time.

**Layout**: 2×2 grid (CH4/N2O × CM/PWB), each cell has scatter + KDE distribution
- **Scatter plot** (left of each pair): Time series with jittered points (to avoid overplotting)
- **KDE plot** (right of each pair): Kernel density estimate of lag distribution
- **Black dashed line**: Mode (most frequent lag value)
- **Event lines**: Tillage (−), Fertilization (:), Precipitation (−−)

**What to look for**:
- Close agreement between CM and PWB clusters → methods agree
- PWB mode sharper than CM → more consistent lag detection
- CM jumps erratically on weak-signal days → PWB smooths out spurious peaks
- N2O KDE peaks narrower than CH4 → N2O has stronger, more stable signal

#### EddyPro vs. PWB Pre-filtered (plot_comparison_strategies.py)

`timelags_strategies_comparison.png` compares EddyPro's detection with the PWB pre-filtered approach:
- **Left (a, c)**: EddyPro (CM) in blue
- **Right (b, d)**: PWB pre-filtered (lags with HDI > 1.0 s removed) in orange

**What to look for**:
- Close agreement between blue and orange clusters → both methods agree
- Orange sparse in certain periods → pre-filtering removed those lags (high uncertainty)
- Orange KDE peak sharper than blue → pre-filtering produces more consistent lags
- Orange avoids spurious jumps that blue might make on weak-signal days

**Interpreting patterns**:
- If plots look very similar: pre-filtering threshold (1.0 s) is effective at removing only truly uncertain lags
- If orange frequently disappears: threshold may be too aggressive; consider increasing it in script 2
- If orange mode differs from blue: pre-filtering changed the consensus lag (evaluate which is more stable)
- If orange peaks are narrower: pre-filtering increased internal consistency

Expected patterns:
- Both methods cluster around 1.7–2.0 s for most of the day on strong-signal periods
- EddyPro may jump erratically on weak-signal afternoons; PWB pre-filtered should stay stable
- N2O (bottom panels) typically has narrower distributions than CH4 (stronger, more stable signal)

Red flags:
- If PWB pre-filtered lags cluster at implausible values (e.g., always −3 s), check input data or analyzer setup
- If pre-filtering removes >50% of lags, consider lowering the threshold (`hdi_prefilter_sec`)

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

## Common Tasks & Troubleshooting

### Processing Multiple Parts

Scripts 1 and 2 process one part at a time. To process all 7 parts:

```r
# Option A: Run script manually 7 times, changing 'part' each time
part <- 1; source("scripts/01_tlag_detection_pwb.R")
part <- 2; source("scripts/01_tlag_detection_pwb.R")
# ... etc

# Option B: Wrap in a loop
for (part in 1:7) {
  source("scripts/01_tlag_detection_pwb.R")
  cat(sprintf("Part %d complete.\n", part))
}
```

Script 2 processes all parts in one run (configured at lines 12).

### Testing Different Pre-filter Thresholds

To compare multiple thresholds in script 2:

```r
for (threshold in c(0.5, 0.75, 1.0, 1.5, 2.0)) {
  hdi_prefilter_sec <- threshold
  source("scripts/02_tlag_compare_pwbopt_strategies.R")
  cat(sprintf("Threshold %.2f s complete.\n", threshold))
}
```

Then examine the output CSVs to see which threshold best balances reliability and data retention.

### Customizing Plots

**Change event dates** (plot.py or plot_comparison_strategies.py):
```python
TILLAGE_DATE = "2021-08-20"        # Line 37
FERTILIZATION_DATE = "2021-07-29"  # Line 38
PRECIPITATION_DATE = "2021-08-16"  # Line 39
```

**Change colors**:
```python
COLOR_STD = "#0072B2"      # Blue (hex color code)
COLOR_PREFILTER = "#E05C2A" # Orange
COLOR_F1 = "#0072B2"       # CM (plot.py)
COLOR_F2 = "#E05C2A"       # PWB (plot.py)
```

**Change output size/quality**:
```python
DPI = 300                  # Resolution (line 41)
FIG_W = 28 / 2.54          # Width in inches (line 129)
FIG_H = 14 / 2.54          # Height in inches (line 130)
```

### Troubleshooting

**Script 1 crashes mid-run**:
- Check `output/tlag_results_checkpoint_partN.rds` — it tracks progress
- Run script again; it resumes from the last completed file
- If checkpoint is corrupted, delete it and restart the part

**Plot.py says "file not found"**:
- Ensure `output/tlag_results_part*.csv` files exist (run script 1 first)
- Verify file paths in `FILE1` and `FILE2` (lines 19–42)
- Check working directory: run script from `scripts/visualization/` directory

**Prefiltered plot looks identical to standard**:
- Pre-filter threshold may be too high
- Try lowering `hdi_prefilter_sec` in script 2 (line 13)
- Or adjust threshold in plot script labels if only visualization changed

**HDI ranges are all very wide (>5 s)**:
- Signal may be too weak; check input data quality
- N2O typically has narrower HDI than CH4 (expected)
- If expected 1–3 s lags have >5 s HDI, investigate analyzer setup or atmospheric conditions

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
