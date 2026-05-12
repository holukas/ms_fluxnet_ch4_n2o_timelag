# ==============================================================================
# 02_tlag_compare_pwbopt_strategies.R
#
# Purpose:
#   Compares two PWBOPT selection strategies to test whether pre-filtering
#   high-uncertainty lags improves lag detection reliability and impacts
#   flux estimation.
#
# Strategies compared:
#   (1) Standard (baseline):
#       Apply S1/S2/S3 logic to all detected lags (as per Vitale et al. 2024)
#
#   (2) Pre-filtered (alternative):
#       Remove lags with HDI > threshold, then apply S1/S2/S3 logic
#       (more conservative: rejects spurious detections upfront)
#
# Workflow:
#   1. Read tlag_results_partN.csv (output from 01_tlag_detection_pwb.R)
#   2. For each part:
#      a. Apply standard S1/S2/S3 to all raw lags
#      b. Apply HDI pre-filter (default: 1.0 second) to remove uncertain lags
#      c. Apply S1/S2/S3 to pre-filtered lags
#      d. Output both sets of results for comparison
#   3. Generate summary statistics (% reliable, # filtered out)
#
# Inputs:
#   tlag_results_part[1-7].csv  (from 01_tlag_detection_pwb.R)
#
# Outputs:
#   tlag_results_prefiltered_part[1-7].csv  -- results with both approaches
#   (columns: ch4/n2o_pwbopt_std, ch4/n2o_flag_std,
#            ch4/n2o_pwbopt_prefilter, ch4/n2o_flag_prefilter)
#
# Configuration:
#   hdi_prefilter_sec (line 13): Adjustable threshold for pre-filtering
#                                (default 1.0 s; test different values)
#   parts_to_process  (line 12): Which parts to process
# ==============================================================================

setwd("F:/Sync/luhk_work/20 - CODING/29 - WORKBENCH/ms_fluxnet_ch4_n2o_timelag")
library(data.table)

# ==============================================================================
# SETTINGS: Adjust these to test different pre-filter thresholds
# ==============================================================================
parts_to_process  <- 1:7        # Which parts to process
hdi_prefilter_sec <- 1.0        # Pre-filter threshold (seconds)
                                # Lags with HDI > this are set to NA before S1/S2/S3

# PWBOPT thresholds (same as main script)
HDI_THRESH_S      <- 0.5
DEV_THRESH_S      <- 0.5

# ==============================================================================
# PWBOPT LOGIC (from main script)
# ==============================================================================
apply_pwbopt <- function(tlag_sec, hdi_range_sec,
                         hdi_thresh = HDI_THRESH_S,
                         dev_thresh = DEV_THRESH_S) {
  n            <- length(tlag_sec)
  flag         <- rep("S3_unreliable", n)
  optimal      <- rep(NA_real_, n)
  last_optimal <- NA_real_

  for (i in seq_len(n)) {
    tl  <- tlag_sec[i]
    hdi <- hdi_range_sec[i]

    if (is.na(tl) || is.na(hdi)) {
      optimal[i] <- last_optimal
      next
    }

    if (hdi < hdi_thresh) {                                         # S1
      flag[i]      <- "S1_optimal"
      optimal[i]   <- tl
      last_optimal <- tl

    } else if (!is.na(last_optimal) &&
               abs(tl - last_optimal) <= dev_thresh) {             # S2
      flag[i]      <- "S2_optimal"
      optimal[i]   <- tl
      last_optimal <- tl

    } else {                                                        # S3
      optimal[i] <- last_optimal
    }
  }

  data.frame(pwbopt_sec = optimal, flag = flag, stringsAsFactors = FALSE)
}

# ==============================================================================
# PRE-FILTER FUNCTION
# ==============================================================================
apply_hdi_prefilter <- function(tlag_sec, hdi_range_sec, threshold) {
  # Set lags with HDI > threshold to NA
  tlag_filtered <- tlag_sec
  tlag_filtered[hdi_range_sec > threshold & !is.na(hdi_range_sec)] <- NA_real_
  tlag_filtered
}

# ==============================================================================
# MAIN LOOP: Load all parts and combine (preserving continuity)
# ==============================================================================
all_frames <- list()

for (part in parts_to_process) {
  input_csv <- sprintf("output/tlag_results_part%d.csv", part)

  if (!file.exists(input_csv)) {
    cat(sprintf("Part %d: %s not found, skipping.\n", part, input_csv))
    next
  }

  cat(sprintf("Loading Part %d...\n", part))
  all_frames[[part]] <- fread(input_csv)
}

# Combine all parts into one dataset (maintains temporal continuity for S1/S2/S3)
results <- rbindlist(all_frames, fill = TRUE)

cat(sprintf("\nApplying PWBOPT logic to %d rows (all parts combined)...\n", nrow(results)))

# --- Standard approach (no pre-filter) ---
ch4_std <- apply_pwbopt(results$ch4_tlag_sec, results$ch4_hdi_range_sec)
n2o_std <- apply_pwbopt(results$n2o_tlag_sec, results$n2o_hdi_range_sec)

# --- Pre-filter approach ---
ch4_filt <- apply_hdi_prefilter(results$ch4_tlag_sec,
                                 results$ch4_hdi_range_sec,
                                 hdi_prefilter_sec)
n2o_filt <- apply_hdi_prefilter(results$n2o_tlag_sec,
                                 results$n2o_hdi_range_sec,
                                 hdi_prefilter_sec)

ch4_pf <- apply_pwbopt(ch4_filt, results$ch4_hdi_range_sec)
n2o_pf <- apply_pwbopt(n2o_filt, results$n2o_hdi_range_sec)

# --- Combine results ---
results[, ch4_pwbopt_std       := ch4_std$pwbopt_sec]
results[, ch4_flag_std         := ch4_std$flag]
results[, ch4_pwbopt_prefilter := ch4_pf$pwbopt_sec]
results[, ch4_flag_prefilter   := ch4_pf$flag]

results[, n2o_pwbopt_std       := n2o_std$pwbopt_sec]
results[, n2o_flag_std         := n2o_std$flag]
results[, n2o_pwbopt_prefilter := n2o_pf$pwbopt_sec]
results[, n2o_flag_prefilter   := n2o_pf$flag]

# Save combined result
output_csv <- "output/tlag_results_prefiltered_all.csv"
fwrite(results, output_csv)
cat(sprintf("Saved combined: %s\n", output_csv))

# Also save per-part for reference
for (part in parts_to_process) {
  if (part %in% names(all_frames)) {
    part_idx <- which(all_frames[[part]]$file_name %in% results$file_name)
    part_results <- results[part_idx]
    output_csv <- sprintf("output/tlag_results_prefiltered_part%d.csv", part)
    fwrite(part_results, output_csv)
    cat(sprintf("  -> Part %d: %d rows\n", part, nrow(part_results)))
  }
}

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================
cat("\n", strrep("=", 80), "\n", sep = "")
cat("COMPARISON: Standard vs. Pre-filtered (HDI > ", hdi_prefilter_sec, " s)\n", sep = "")
cat(strrep("=", 80), "\n\n", sep = "")

pct_reliable <- function(flags) {
  round(100 * mean(flags %in% c("S1_optimal", "S2_optimal"), na.rm = TRUE), 1)
}

for (part in parts_to_process) {
  if (!part %in% names(all_frames)) next

  # Extract this part's results from the combined dataset
  part_idx <- which(all_frames[[part]]$file_name %in% results$file_name)
  part_results <- results[part_idx]

  cat(sprintf("PART %d:\n", part))
  cat("  CH4 Standard:     "); cat(sprintf("%.1f%% reliable\n", pct_reliable(part_results$ch4_flag_std)))
  cat("  CH4 Pre-filtered: "); cat(sprintf("%.1f%% reliable\n", pct_reliable(part_results$ch4_flag_prefilter)))
  cat("  N2O Standard:     "); cat(sprintf("%.1f%% reliable\n", pct_reliable(part_results$n2o_flag_std)))
  cat("  N2O Pre-filtered: "); cat(sprintf("%.1f%% reliable\n", pct_reliable(part_results$n2o_flag_prefilter)))

  # Count how many were pre-filtered
  ch4_prefiltered <- sum(part_results$ch4_hdi_range_sec > hdi_prefilter_sec &
                         !is.na(part_results$ch4_hdi_range_sec), na.rm = TRUE)
  n2o_prefiltered <- sum(part_results$n2o_hdi_range_sec > hdi_prefilter_sec &
                         !is.na(part_results$n2o_hdi_range_sec), na.rm = TRUE)

  cat(sprintf("  CH4 pre-filtered out: %d / %d\n", ch4_prefiltered, nrow(part_results)))
  cat(sprintf("  N2O pre-filtered out: %d / %d\n", n2o_prefiltered, nrow(part_results)))
  cat("\n")
}

cat(strrep("=", 80), "\n", sep = "")
cat("Output files: tlag_results_prefiltered_part[1-7].csv\n")
cat("Each file contains:\n")
cat("  - Standard approach:  ch4_pwbopt_std, ch4_flag_std (etc.)\n")
cat("  - Pre-filtered:       ch4_pwbopt_prefilter, ch4_flag_prefilter (etc.)\n")
