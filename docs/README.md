# Time Lag Detection for Eddy Covariance Measurements

Detects time lags between vertical wind and scalar concentrations in eddy covariance flux measurements using the pre-whitening with block-bootstrap (PWB) method from Vitale et al. (2024). Designed for weak-signal cases where standard cross-correlation fails.

## The Problem

Gas analyzers (especially closed-path systems) introduce physical delays via tubing or optical path differences. Scalar time series also lag due to instrumental response time. Standard cross-correlation works well for strong signals but struggles when flux correlations are near zero. With weak signals, multiple spurious correlation peaks appear with similar heights, making it impossible to identify the true lag by eye, leading to incorrect lag estimates and systematically biased flux calculations.

## The Solution: PWB (Vitale et al., 2024)

The pre-whitening with block-bootstrap method improves lag detection on weak signals:

1. **Pre-whiten** — Filter out autocorrelation from each time series via AR model
2. **Block-bootstrap** — Resample 99 times in ~20-step blocks to preserve short-timescale structure
3. **Quantify uncertainty** — Extract lag from the mode of bootstrap distribution; 95% HDI width indicates confidence
4. **Apply S1/S2/S3 logic** — Accept S1 (high confidence), accept S2 (within continuity window), reject S3 (spurious)

Output: detected lag, 95% highest-density interval (HDI), and reliability flag per 30-min window.

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

**Input**: 30-minute windows of high-frequency eddy covariance data (≥20 Hz sampling). Includes rotation-corrected wind components (U, V, W via planar-fit or double-rotation), sonic temperature, and scalar concentrations (e.g., CH₄, N₂O). The algorithm computes cross-correlation between vertical wind turbulent fluctuations (w') and scalar fluctuations (e.g., c') via Reynolds decomposition.

**Output**: CSV with detected lag, 95% HDI bounds (uncertainty), correlation magnitude, and S1/S2/S3 reliability flags. One row per 30-min averaging window. See [CLAUDE.md](CLAUDE.md#understanding-the-output) for detailed column descriptions.

## References

- Vitale et al. (2024): A pre-whitening with block-bootstrap cross-correlation procedure for temporal alignment of data sampled by eddy covariance systems. *Environmental and Ecological Statistics* 31:219–244. [10.1007/s10651-024-00615-9](https://doi.org/10.1007/s10651-024-00615-9)
- RFlux package: https://github.com/icos-etc/RFlux (Vitale et al. 2022)
- See `references/` folder for full papers and RFlux manual.

## Requirements

- R ≥ 4.0 (dependencies managed via renv)
- Python ≥ 3.9 (dependencies managed via uv)
