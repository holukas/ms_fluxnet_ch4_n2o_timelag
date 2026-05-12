# ==============================================================================
# Script: RFlux Time Lag Detection Initialization & Testing
# Purpose: Sets up an isolated R environment, installs the RFlux package from 
#          source, and tests the pre-whitening time lag detection algorithm 
#          on sample eddy covariance data.
# References: 
#   - RFlux Manual (Vitale et al., 2022)
#   - Vitale et al. (2024) - A pre-whitening with block-bootstrap cross-correlation...
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Environment & Dependency Management
# ------------------------------------------------------------------------------

# Initialize an isolated virtual environment to prevent package version conflicts.
install.packages("renv")
renv::init()

# The 'egcm' package is required by RFlux but is archived/unavailable for R 4.5.2 on CRAN.
# We bypass CRAN and install it directly from the GitHub mirror.
remotes::install_github("cran/egcm")

# ------------------------------------------------------------------------------
# 2. RFlux Package Installation
# ------------------------------------------------------------------------------

# Install devtools to compile RFlux from its source code repository.
install.packages("devtools")

# Install RFlux directly from the ICOS-ETC GitHub repository.
devtools::install_github("icos-etc/RFlux")

# Display the package documentation and citation information for reference.
help(package=RFlux)
citation("RFlux")

# Load required libraries into the workspace.
library(RFlux)
library(forecast) # Used internally by RFlux for ARIMA/AR modeling
library(xts)      # Required for handling high-frequency time series objects

# ------------------------------------------------------------------------------
# 3. Data Preparation
# ------------------------------------------------------------------------------

# Load the built-in raw data from the package. 
# 'closed_path_rawdata' is an example of raw, high-frequency eddy covariance data 
# for a closed-path system.
data("closed_path_rawdata")
rawdata <- closed_path_rawdata

# Count the number of observations (typically 36,000 for a 30-min file at 20Hz).
N <- nrow(rawdata)

# Create a dummy timestamp sequence for the xts object.
# This aligns the raw data points into a strict chronological sequence for analysis.
timestamp_orig <- strptime("0000.00", format="%H%M.%OS", tz="GMT") + seq(1, N, 1)
data.xts <- xts(rawdata, order.by = timestamp_orig)

# ------------------------------------------------------------------------------
# 4. Time Lag Detection using Pre-Whitening
# ------------------------------------------------------------------------------

# Test the tlag_detection function between water vapor (H2O) and vertical wind speed (W).
# We supply the scalar variable, the sonic temperature, the vertical wind, 
# the number of bootstrap iterations (Rboot), and tell it to plot (plot.it).

tlag_h2o_out <- tlag_detection(
  scalar_var = rawdata$H2O,      # The scalar atmospheric variable (gas concentration)
  tsonic_var = rawdata$T_SONIC,  # The sonic temperature
  w_var      = rawdata$W,        # The vertical wind velocity component
  mfreq      = 20,               # The acquisition frequency in Hz
  Rboot      = 3,              # Number of bootstrap replicates (required by the new PWB method)
  plot.it    = TRUE              # Plots the resulting clean CCF (renamed from show.plot)
)

# Print the detected optimal time lag (using the PWB method output)
print(paste("Optimal Time Lag (PWB method):", tlag_h2o_out$pwb, "timesteps"))

# Print the estimated correlation at the optimal time lag
print(paste("Correlation Estimate:", tlag_h2o_out$cor_pwb))
