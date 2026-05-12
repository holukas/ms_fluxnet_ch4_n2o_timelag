# Time Lag Detection for Closed-Path Eddy Covariance Measurements

Detects time lags in closed-path gas analyzer measurements using the pre-whitening with block-bootstrap (PWB) method — an improvement over standard cross-correlation for weak-signal cases.

## The Problem

Closed-path analyzers measure gases through tubing, introducing a physical delay between sampling and detection. Standard cross-correlation works well for strong signals but fails on weak ones. When flux correlations are near zero (typical for trace gases), multiple spurious peaks appear with similar heights, leading to incorrect lags and biased flux estimates.

## How PWB Works

1. Pre-whiten the time series (filter out autocorrelation via AR model)
2. Block-bootstrap resample 99 times to quantify lag uncertainty
3. Use S1/S2/S3 logic to accept reliable lags and reject artifacts

Output: detected lag, 95% confidence interval (HDI), and a reliability flag.

## Running It

R time lag detection:

```bash
source("scripts/01_tlag_detection_pwb.R")
```

Python visualization:

```bash
cd scripts/visualization
uv sync
uv run python plot.py
```

For detailed instructions and workflow, see [CLAUDE.md](CLAUDE.md).

## Input & Output

**Input**: 30-minute windows of high-frequency eddy covariance data at ≥20 Hz sampling. Wind components must be rotation-corrected (planar fit or double rotation). Includes sonic temperature and gas concentrations (CH₄, N₂O, or both). PWB calculates the cross-correlation between turbulent fluctuations of the vertical wind (w') and scalar (e.g., N₂O'), obtained via Reynolds decomposition.

**Output**: CSV with detected lag, HDI bounds, and S1/S2/S3 flags for each window. See [CLAUDE.md](CLAUDE.md#understanding-the-output) for detailed column descriptions.

## References

- Vitale et al. (2024): A pre-whitening with block-bootstrap cross-correlation procedure for temporal alignment of data sampled by eddy covariance systems. *Environmental and Ecological Statistics* 31:219–244. [10.1007/s10651-024-00615-9](https://doi.org/10.1007/s10651-024-00615-9)
- RFlux package: https://github.com/icos-etc/RFlux (Vitale et al. 2022)
- See `references/` folder for full papers and RFlux manual.

## Requirements

- R ≥ 4.0 (dependencies managed via renv)
- Python ≥ 3.9 (dependencies managed via uv)
