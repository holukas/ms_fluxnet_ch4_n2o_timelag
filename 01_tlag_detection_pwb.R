# ==============================================================================
# 01_tlag_detection_pwb.R
#
# Purpose:
#   Detects time lags in CH4 and N2O eddy covariance measurements using the
#   pre-whitening with block-bootstrap (PWB) method from Vitale et al. (2024).
#   Processes 30-min periods from rotation-corrected wind/scalar data files.
#
# Workflow:
#   1. Load high-frequency data from EddyPro-rotated files (one part at a time)
#   2. Run RFlux tlag_detection() with 99 bootstrap replicates on CH4 and N2O
#   3. Extract raw lags, HDI bounds, and correlations
#   4. Apply S1/S2/S3 selection logic (Vitale et al. Sect. 2.3):
#      - S1: HDI < 0.5s → reliable, accept
#      - S2: HDI ≥ 0.5s but within 0.5s of prior → reliable, accept for continuity
#      - S3: High uncertainty or far from prior → unreliable, use last good lag
#   5. Generate 6-panel diagnostic plots (PDF) for each gas/period
#   6. Output CSV with detected lags, HDI bounds, correlations, PWBOPT flags
#
# Inputs:
#   03-rotated_data_from_eddypro_level5_partN/*.txt  (rotation-corrected, EddyPro L5)
#
# Outputs:
#   tlag_results_partN.csv                 -- raw lags + PWBOPT selection
#   tlag_plots/<stem>_ch4.pdf              -- diagnostic plots for CH4
#   tlag_plots/<stem>_n2o.pdf              -- diagnostic plots for N2O
#   tlag_results_checkpoint_partN.rds      -- auto-resume checkpoint
#
# Notes:
#   - Change 'part' variable (line 20) to process different dataset parts (1-7)
#   - Script resumes from last completed file if interrupted (checkpoint system)
#   - Requires RFlux (github.com/icos-etc/RFlux)
# ==============================================================================

setwd("F:/Sync/luhk_work/20 - CODING/29 - WORKBENCH/ms_fluxnet_ch4_n2o_timelag")
renv::restore()
library(data.table)
library(RFlux)

# ------------------------------------------------------------------------------
# 1. SETTINGS
# ------------------------------------------------------------------------------
part <- 5   # <-- change this number to switch between dataset parts

folder_path      <- sprintf("03-rotated_data_from_eddypro_level5_part%d", part)
output_csv       <- sprintf("tlag_results_part%d.csv", part)
plot_dir         <- "tlag_plots"
checkpoint_file  <- sprintf("tlag_results_checkpoint_part%d.rds", part)
file_pattern     <- "\\.txt$"
mfreq            <- 20    # Hz

n_boot           <- 99    # Vitale et al. 2024: N_B = 99 bootstrap samples
wdt              <- floor(mfreq / 2) + 1  # = 11 steps; paper: hz/2 + 1 (Sect. 2.2)

# PWBOPT thresholds (Section 2.3, Vitale et al. 2024)
HDI_THRESH_S     <- 0.5   # seconds: HDI range below this = low uncertainty (S1)
DEV_THRESH_S     <- 0.5   # seconds: max deviation from preceding optimal lag (S2)

# Minimum fraction of non-NA values required before attempting detection.
# Lowered to 0.3 to avoid discarding marginal but usable periods (e.g. during
# sensor maintenance or rain events) while still guarding against the NA p-value
# crash in bvr.test() for near-empty series.
MIN_VALID_FRAC   <- 0.3

# ------------------------------------------------------------------------------
# 2. CREATE OUTPUT DIRECTORIES
# ------------------------------------------------------------------------------
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 3. PWBOPT SELECTION LOGIC (S1 -> S2 -> S3)
# ------------------------------------------------------------------------------
apply_pwbopt <- function(tlag_sec, hdi_range_sec,
                         hdi_thresh = HDI_THRESH_S,
                         dev_thresh  = DEV_THRESH_S) {
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

# ------------------------------------------------------------------------------
# 4. HELPERS
# ------------------------------------------------------------------------------

# Convert a sample-unit field from tlag_detection result to seconds
extr <- function(res, field) {
  if (is.null(res)) return(NA_real_)
  val <- res[[field]]
  if (is.null(val)) return(NA_real_)
  as.numeric(val) / mfreq
}

# HDI width in seconds
hdi_range <- function(res) {
  if (is.null(res)) return(NA_real_)
  (res$pwb_uci - res$pwb_lci) / mfreq
}

# Guard: returns TRUE only when a series has enough valid, non-constant data.
valid_enough <- function(x) {
  if (mean(!is.na(x)) < MIN_VALID_FRAC) return(FALSE)
  if (sd(x, na.rm = TRUE) < .Machine$double.eps) return(FALSE)
  TRUE
}

# Empty result row — defined once outside the loop
make_empty <- function(file_name) data.frame(
  file_name         = file_name,
  ch4_tlag_sec      = NA_real_, ch4_hdi_lci_sec   = NA_real_,
  ch4_hdi_uci_sec   = NA_real_, ch4_hdi_range_sec = NA_real_,
  ch4_cor           = NA_real_,
  n2o_tlag_sec      = NA_real_, n2o_hdi_lci_sec   = NA_real_,
  n2o_hdi_uci_sec   = NA_real_, n2o_hdi_range_sec = NA_real_,
  n2o_cor           = NA_real_,
  stringsAsFactors  = FALSE
)

# Detect time lag for one gas, writing the plot to a PDF.
# The pdf() device is opened and closed inside this function so that an error
# in tlag_detection never leaves an orphaned open device.
detect_gas <- function(scalar, tsonic, w, gas_label, pdf_path) {
  pdf(pdf_path, width = 14, height = 10)
  res <- tryCatch({
    tlag_detection(scalar_var = scalar, tsonic_var = tsonic, w_var = w,
                   mfreq = mfreq, Rboot = n_boot, wdt = wdt, plot.it = TRUE)
  }, error = function(e) {
    message(sprintf("  -> %s error: %s", gas_label, e$message))
    NULL
  })
  # Always close the device, whether detection succeeded or failed
  if (dev.cur() > 1) dev.off()
  res
}

# ------------------------------------------------------------------------------
# 5. RESUME FROM CHECKPOINT IF AVAILABLE
# ------------------------------------------------------------------------------
file_list <- sort(list.files(folder_path, pattern = file_pattern, full.names = TRUE))
if (length(file_list) == 0) stop("No files found in the specified folder.")
cat(sprintf("Found %d files to process.\n", length(file_list)))

if (file.exists(checkpoint_file)) {
  results_list <- readRDS(checkpoint_file)
  # Last completed index = last non-NULL entry
  completed    <- which(!sapply(results_list, is.null))
  start_i      <- if (length(completed) > 0) max(completed) + 1 else 1
  cat(sprintf("Checkpoint found — resuming from file %d / %d\n",
              start_i, length(file_list)))
} else {
  results_list <- vector("list", length(file_list))
  start_i      <- 1
}

# ------------------------------------------------------------------------------
# 6. MAIN LOOP
# ------------------------------------------------------------------------------
for (i in seq(start_i, length(file_list))) {
  current_file <- file_list[i]
  file_name    <- basename(current_file)
  file_stem    <- tools::file_path_sans_ext(file_name)
  cat(sprintf("[%d/%d] %s\n", i, length(file_list), file_name))
  
  # --- Load data ---
  rawdata <- tryCatch({
    dat <- fread(current_file,
                 skip         = 9,
                 header       = TRUE,
                 sep          = "auto",
                 check.names  = FALSE,
                 na.strings   = c("-9999.0", "-9999.0000000000000", "-9999"),
                 select       = c(1, 2, 3, 4, 7, 8),
                 showProgress = FALSE)
    setnames(dat, new = c("u", "v", "w", "ts", "ch4", "n2o"))
    dat
  }, error = function(e) {
    message(sprintf("  -> Read error: %s", e$message))
    NULL
  })
  
  if (is.null(rawdata) || nrow(rawdata) < 25) {
    cat("  -> Skipped (insufficient rows).\n")
    results_list[[i]] <- make_empty(file_name)
    saveRDS(results_list, checkpoint_file)
    next
  }
  
  W       <- as.numeric(rawdata$w)
  T_SONIC <- as.numeric(rawdata$ts)
  CH4     <- as.numeric(rawdata$ch4)
  N2O     <- as.numeric(rawdata$n2o)
  
  shared_ok <- valid_enough(W) && valid_enough(T_SONIC)
  
  # --- CH4 ---
  if (shared_ok && valid_enough(CH4)) {
    tlag_ch4 <- detect_gas(CH4, T_SONIC, W,
                           gas_label = "CH4",
                           pdf_path  = file.path(plot_dir,
                                                 paste0(file_stem, "_ch4.pdf")))
  } else {
    cat("  -> CH4 skipped (too many NAs or constant series).\n")
    tlag_ch4 <- NULL
  }
  
  # --- N2O ---
  if (shared_ok && valid_enough(N2O)) {
    tlag_n2o <- detect_gas(N2O, T_SONIC, W,
                           gas_label = "N2O",
                           pdf_path  = file.path(plot_dir,
                                                 paste0(file_stem, "_n2o.pdf")))
  } else {
    cat("  -> N2O skipped (too many NAs or constant series).\n")
    tlag_n2o <- NULL
  }
  
  # --- Collect results ---
  results_list[[i]] <- data.frame(
    file_name         = file_name,
    ch4_tlag_sec      = extr(tlag_ch4, "pwb"),
    ch4_hdi_lci_sec   = extr(tlag_ch4, "pwb_lci"),
    ch4_hdi_uci_sec   = extr(tlag_ch4, "pwb_uci"),
    ch4_hdi_range_sec = hdi_range(tlag_ch4),
    ch4_cor           = if (!is.null(tlag_ch4)) tlag_ch4$cor_pwb else NA_real_,
    n2o_tlag_sec      = extr(tlag_n2o, "pwb"),
    n2o_hdi_lci_sec   = extr(tlag_n2o, "pwb_lci"),
    n2o_hdi_uci_sec   = extr(tlag_n2o, "pwb_uci"),
    n2o_hdi_range_sec = hdi_range(tlag_n2o),
    n2o_cor           = if (!is.null(tlag_n2o)) tlag_n2o$cor_pwb else NA_real_,
    stringsAsFactors  = FALSE
  )
  
  cat(sprintf("  -> Done. CH4: %s s | N2O: %s s\n",
              ifelse(is.na(results_list[[i]]$ch4_tlag_sec), "NA",
                     sprintf("%.3f", results_list[[i]]$ch4_tlag_sec)),
              ifelse(is.na(results_list[[i]]$n2o_tlag_sec), "NA",
                     sprintf("%.3f", results_list[[i]]$n2o_tlag_sec))))
  
  # Save checkpoint after every file so a crash can be resumed
  saveRDS(results_list, checkpoint_file)
}

# ------------------------------------------------------------------------------
# 7. COMBINE AND APPLY PWBOPT
# ------------------------------------------------------------------------------
results <- rbindlist(results_list)

cat("\nApplying PWBOPT selection logic...\n")

ch4_opt <- apply_pwbopt(results$ch4_tlag_sec, results$ch4_hdi_range_sec)
n2o_opt <- apply_pwbopt(results$n2o_tlag_sec, results$n2o_hdi_range_sec)

results[, ch4_pwbopt_sec := ch4_opt$pwbopt_sec]
results[, ch4_flag       := ch4_opt$flag]
results[, n2o_pwbopt_sec := n2o_opt$pwbopt_sec]
results[, n2o_flag       := n2o_opt$flag]

# ------------------------------------------------------------------------------
# 8. SAVE CSV + PRINT SUMMARY
# ------------------------------------------------------------------------------
pct <- function(flags) {
  round(100 * mean(flags %in% c("S1_optimal", "S2_optimal"), na.rm = TRUE), 1)
}

cat("\n--- PWBOPT Flag Summary ---\n")
cat("CH4:\n"); print(table(results$ch4_flag, useNA = "ifany"))
cat("N2O:\n"); print(table(results$n2o_flag, useNA = "ifany"))
cat(sprintf("\nCH4 reliable (S1+S2): %.1f%%\n", pct(results$ch4_flag)))
cat(sprintf("N2O reliable (S1+S2): %.1f%%\n",   pct(results$n2o_flag)))

write.csv(results, file = output_csv, row.names = FALSE)
cat(sprintf("\nResults    : %s\n", output_csv))
cat(sprintf("Plots      : %s/\n", plot_dir))
cat(sprintf("Checkpoint : %s\n",  checkpoint_file))

print(head(results[, .(file_name,
                       ch4_tlag_sec, ch4_hdi_range_sec, ch4_pwbopt_sec, ch4_flag,
                       n2o_tlag_sec, n2o_hdi_range_sec, n2o_pwbopt_sec, n2o_flag)]))