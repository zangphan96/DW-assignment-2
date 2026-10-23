# Reproduction of the CATEGORICAL experiments (Section 4.2, Figures 2 and 4,
# Table 1) of
#   Stekhoven & Buehlmann (2012), "MissForest - non-parametric missing value
#   imputation for mixed-type data"
#
# Only missForest is run.
#
# Simulation plan (Section 3 of the paper):
#   * values removed completely at random (MCAR) at 10%, 20%, 30%
#   * 50 independent simulations per dataset and missingness rate
#   * true error: PFC (proportion of falsely classified entries)
#   * OOB error estimate of missForest recorded in every run
#   * missForest default parameters: ntree = 100, mtry = floor(sqrt(p)), maxiter = 10
#
# Datasets (all three categorical datasets of the paper):
#   SPECT        n = 267, p = 22  (22 binary features, diagnosis removed)
#   Promoter     n = 106, p = 57  (57 nucleotide positions a/c/g/t, class removed)
#   Lymphography n = 148, p = 19  (18 attributes + class, all categorical)
#
# Outputs:
#   results/raw_results.csv   one row per simulation run
#   Summary table             printed in the console
#   Figure 2 and Figure 4     drawn in the plot window


## ---------------------------------------------------------------------------
## 0. Settings
## ---------------------------------------------------------------------------
QUICK      <- FALSE           # TRUE = 3 simulations per setting, for testing
N_SIM      <- if (QUICK) 3 else 50
MISS_RATES <- c(0.1, 0.2, 0.3)       # fractions of entries deleted (MCAR)
NTREE      <- 100                    # trees per random forest (paper default)
MAXITER    <- 10                     # max. number of missForest iterations
# The paper used the 'randomForest' package. missForest 1.6+ defaults to
# 'ranger' instead, so request randomForest explicitly to match the paper.
BACKEND    <- "randomForest"
N_CORES    <- max(1, parallel::detectCores() - 1)  # simulations run in parallel
SEED       <- 2012                   # base seed; each run derives its own seed
OUT_DIR    <- "results"
DATA_DIR         <- "data"
# Expected locations of the raw UCI data files (relative to the working dir)
SPECT_TRAIN_FILE <- file.path(DATA_DIR, "SPECT.train")
SPECT_TEST_FILE  <- file.path(DATA_DIR, "SPECT.test")
PROMOTER_FILE    <- file.path(DATA_DIR, "promoters.data")
LYMPHO_FILE      <- file.path(DATA_DIR, "lymphography.data")

# Average missForest runtimes reported in Table 1 of the paper [s]
PAPER_RUNTIME <- c(SPECT = 5.5, Promoter = 38, Lymphography = 7.0)

# Install any missing packages, then load them
for (pkg in c("missForest", "randomForest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(missForest)
library(parallel)
# Print the version: behaviour (e.g. default backend) differs between versions
cat("missForest version:", as.character(packageVersion("missForest")), "\n")

# Create the output folder (no warning if it already exists)
dir.create(OUT_DIR, showWarnings = FALSE)

# Stop early with a clear message if a data file is missing
files <- c(SPECT_TRAIN_FILE, SPECT_TEST_FILE, PROMOTER_FILE, LYMPHO_FILE)
if (any(!file.exists(files)))
  stop("Data file(s) not found:\n  ",
       paste(files[!file.exists(files)], collapse = "\n  "),
       "\nDownload them (links in the settings above) or fix the paths.")

## ---------------------------------------------------------------------------
## 1. Data loaders (each returns a complete data.frame of factors)
## ---------------------------------------------------------------------------

# SPECT heart data (Kurgan et al. 2001): train + test = 267 patients
load_spect <- function() {
  # Stack the training (80) and test (187) files into one dataset
  d <- rbind(read.csv(SPECT_TRAIN_FILE, header = FALSE),
             read.csv(SPECT_TEST_FILE,  header = FALSE))
  d <- d[, -1]                        # column 1 = diagnosis; keep 22 features
  names(d) <- paste0("F", 1:22)
  # Values are 0/1 integers; convert to factors so missForest treats them as
  # categorical (classification forests, PFC error)
  d[] <- lapply(d, factor)
  d
}

# E. coli promoter gene sequences (Harley & Reynolds 1987; Towell et al. 1990)
# Each line looks like:  +,S10,<tabs>tactagcaatacgcttgcg...   (57 letters)
load_promoter <- function() {
  lines <- readLines(PROMOTER_FILE)
  lines <- lines[nzchar(trimws(lines))]           # drop empty lines
  # Take the 3rd comma-separated field (the sequence), strip tabs/spaces and
  # lower-case it so all nucleotides are a/c/g/t
  seqs  <- sapply(strsplit(lines, ","), function(x) tolower(gsub("[[:space:]]", "", x[3])))
  stopifnot(all(nchar(seqs) == 57))               # every sequence must be complete
  # Split each sequence into single letters: one row per sequence,
  # one column per nucleotide position
  d <- as.data.frame(do.call(rbind, strsplit(seqs, "")), stringsAsFactors = FALSE)
  names(d) <- paste0("p", c(-50:-1, 1:7))         # positions -50 ... +7
  # Fix the level set to all 4 nucleotides, even if one is absent at a position
  d[] <- lapply(d, factor, levels = c("a", "c", "g", "t"))
  d                                               # class (+/-) not used: p = 57
}

# Lymphography domain data (Zwitter & Soklic); 148 patients, 19 variables.
# All values are coded as integers; every variable is treated as categorical.
load_lymphography <- function() {
  d <- read.csv(LYMPHO_FILE, header = FALSE)
  # Variable names from the UCI description (lymphography.names)
  names(d) <- c("class", "lymphatics", "block_affere", "bl_lymph_c",
                "bl_lymph_s", "by_pass", "extravasates", "regeneration",
                "early_uptake", "lym_nodes_dimin", "lym_nodes_enlar",
                "changes_lym", "defect_node", "changes_node", "changes_stru",
                "special_forms", "dislocation", "exclusion_no", "no_nodes")
  d[] <- lapply(d, factor)                        # class kept: p = 19
  d
}

# List of datasets to process: a display name plus the function that loads it.
# The order here is also the order used in the table and the figures.
DATASETS <- list(
  list(name = "SPECT",        loader = load_spect),
  list(name = "Promoter",     loader = load_promoter),
  list(name = "Lymphography", loader = load_lymphography)
)

## ---------------------------------------------------------------------------
## 2. One simulation run: delete values, impute, measure the error
## ---------------------------------------------------------------------------
# Arguments:
#   sim      index of the simulation (1 ... N_SIM)
#   xtrue    the complete dataset
#   rate     fraction of entries to delete
#   seed     base seed for this dataset/rate; sim is added to it so every run
#            is reproducible and different from the others
#   ntree, maxiter, backend   missForest settings
# Returns a one-row data.frame with the true PFC, the OOB PFC and the runtime.
one_run <- function(sim, xtrue, rate, seed, ntree, maxiter, backend) {
  set.seed(seed + sim)

  # missForest 1.6.1 bug: its initial mode imputation uses summary(), which
  # includes an "NA's" row. When a column's NA count exactly ties its most
  # frequent level, "NA's" can be drawn as the fill value and the run crashes
  # ("NA not permitted in predictors"). Rare, but likely at 30% missing.
  # Workaround: redraw the MCAR mask (a fresh, equally valid MCAR sample).
  for (attempt in 1:20) {
    xmis <- missForest::prodNA(xtrue, noNA = rate)      # MCAR removal

    # Build the argument list for missForest (mtry left at its default,
    # floor(sqrt(p)), as in the paper)
    args <- list(xmis, ntree = ntree, maxiter = maxiter, verbose = FALSE)
    # Older missForest versions have no 'backend' argument (always randomForest)
    if ("backend" %in% names(formals(missForest::missForest)))
      args$backend <- backend

    # Time the imputation (wall-clock seconds); errors return NULL so the
    # loop can retry with a new mask instead of stopping the whole script
    t0  <- proc.time()[["elapsed"]]
    imp <- tryCatch(do.call(missForest::missForest, args), error = function(e) NULL)
    rt  <- proc.time()[["elapsed"]] - t0
    if (!is.null(imp)) break                            # success: leave the loop
  }
  if (is.null(imp)) stop("missForest failed 20 times in a row for sim ", sim)

  # True PFC over the artificially removed entries only
  err <- missForest::mixError(imp$ximp, xmis, xtrue)

  # imp$OOBerror is missForest's own estimate of the imputation error,
  # computed from the out-of-bag predictions (no access to the true values)
  data.frame(sim       = sim,
             PFC_true  = unname(err[["PFC"]]),
             PFC_OOB   = unname(imp$OOBerror[["PFC"]]),
             runtime_s = rt,
             retries   = attempt - 1)                   # 0 = worked first time
}

## ---------------------------------------------------------------------------
## 3. Run all experiments (simulations in parallel)
## ---------------------------------------------------------------------------
# Start the worker processes and load missForest on each of them
cl <- makeCluster(N_CORES)
clusterEvalQ(cl, library(missForest))

all_res <- list()                       # collects one data.frame per setting
for (ds in DATASETS) {
  xtrue <- ds$loader()
  # Sanity check: the dataset must be complete and fully categorical
  stopifnot(!anyNA(xtrue), all(sapply(xtrue, is.factor)))
  cat(sprintf("\n=== %s: n = %d, p = %d categorical variables\n",
              ds$name, nrow(xtrue), ncol(xtrue)))

  for (rate in MISS_RATES) {
    t0 <- Sys.time()
    # Run the N_SIM simulations on the workers. The seed is shifted by the
    # missingness rate so the three rates use different random numbers.
    runs <- parLapply(cl, seq_len(N_SIM), one_run,
                      xtrue = xtrue, rate = rate,
                      seed = SEED + round(rate * 1000),
                      ntree = NTREE, maxiter = MAXITER, backend = BACKEND)
    runs <- do.call(rbind, runs)        # list of rows -> one data.frame
    # Label the rows with the setting they belong to
    runs$dataset <- ds$name
    runs$miss    <- rate
    all_res[[length(all_res) + 1]] <- runs
    # Progress report for this setting
    cat(sprintf("  %2.0f%% missing: %d runs in %.1f min (mean PFC = %.3f)\n",
                100 * rate, N_SIM,
                as.numeric(difftime(Sys.time(), t0, units = "mins")),
                mean(runs$PFC_true)))
  }
  # Save after every dataset so nothing is lost if a later one fails
  write.csv(do.call(rbind, all_res), file.path(OUT_DIR, "raw_results.csv"),
            row.names = FALSE)
}
stopCluster(cl)                         # shut down the worker processes

# Combine all settings into one long data.frame (one row per run)
res <- do.call(rbind, all_res)
# Factor with fixed level order so tables and plots follow the DATASETS order
res$dataset <- factor(res$dataset, levels = sapply(DATASETS, `[[`, "name"))
res$diff    <- res$PFC_true - res$PFC_OOB          # err_true - err_OOB

## ---------------------------------------------------------------------------
## 4. Summary table
## ---------------------------------------------------------------------------
# One group per dataset x missingness rate combination
groups <- split(res, list(res$dataset, res$miss), drop = TRUE)

# Aggregate each group: means over the simulations, plus the standard error
# of the mean true PFC (used for the error bars in Figure 2)
summary_tab <- do.call(rbind, lapply(groups, function(g) data.frame(
  dataset          = g$dataset[1],
  miss             = g$miss[1],
  n_sim            = nrow(g),
  PFC_true         = mean(g$PFC_true),
  PFC_true_se      = sd(g$PFC_true) / sqrt(nrow(g)),
  PFC_OOB          = mean(g$PFC_OOB),
  mean_diff        = mean(g$diff),
  runtime_s        = mean(g$runtime_s),
  paper_runtime_s  = PAPER_RUNTIME[[as.character(g$dataset[1])]]
)))

# Sort by dataset, then by rate (split() returns groups in a different order)
summary_tab <- summary_tab[order(summary_tab$dataset, summary_tab$miss), ]
rownames(summary_tab) <- NULL

# Readable version of the table: percentages for the rates, 3 decimals elsewhere
show_tab <- data.frame(
  "Dataset"            = summary_tab$dataset,
  "Missing"            = paste0(100 * summary_tab$miss, "%"),
  "Runs"               = summary_tab$n_sim,
  "PFC (true)"         = round(summary_tab$PFC_true, 3),
  "SE"                 = round(summary_tab$PFC_true_se, 4),
  "PFC (OOB)"          = round(summary_tab$PFC_OOB, 3),
  "True - OOB"         = round(summary_tab$mean_diff, 3),
  "Time [s]"           = round(summary_tab$runtime_s, 1),
  "Paper time [s]"     = summary_tab$paper_runtime_s,
  check.names = FALSE)                 # keep spaces/brackets in column names

cat("\n===== Summary (mean over simulations) =====\n")
print(show_tab, row.names = FALSE)
# Save the summary in summary.csv file
write.csv(show_tab, file.path(OUT_DIR, "summary.csv"), row.names = FALSE)

if (interactive()) View(show_tab, "missForest summary")   # table tab in RStudio

## ---------------------------------------------------------------------------
## 5a. Figure 2: average PFC per dataset and missingness rate
## ---------------------------------------------------------------------------
old_par <- par(no.readonly = TRUE)      # remember graphics settings to restore later
# One panel per dataset, side by side
par(mfrow = c(1, length(DATASETS)), mar = c(3, 4, 2.5, 0.5))
for (ds in DATASETS) {
  s  <- summary_tab[summary_tab$dataset == ds$name, ]   # rows for this dataset
  # Upper y limit: leave room above the error bars for the value labels
  up <- max(0.1, 1.3 * max(s$PFC_true + 2 * s$PFC_true_se))
  # Black bars as in the paper; barplot() returns the bar centres (x positions)
  bp <- barplot(s$PFC_true, names.arg = paste0(100 * s$miss, "%"),
                col = "black", ylim = c(0, up), width = 0.6, space = 1,
                ylab = "PFC", main = paste(ds$name, "data"), las = 1)
  # +/- 2 standard errors
  arrows(bp, s$PFC_true - 2 * s$PFC_true_se, bp, s$PFC_true + 2 * s$PFC_true_se,
         angle = 90, code = 3, length = 0.04, col = "grey50")
  # Write the mean PFC above each error bar
  text(bp, s$PFC_true + 2 * s$PFC_true_se, sprintf("%.3f", s$PFC_true),
       pos = 3, cex = 0.8)
  abline(h = 0)                         # baseline
}
par(old_par)

## ---------------------------------------------------------------------------
## 5b. Figure 4 (right panel): err_true - err_OOB, with means of both
## ---------------------------------------------------------------------------
k   <- length(MISS_RATES)
# x positions: 3 boxes per dataset, a gap of 1 between datasets
pos <- unlist(lapply(seq_along(DATASETS) - 1, function(i) i * (k + 1) + seq_len(k)))
lab <- rep(MISS_RATES, length(DATASETS))   # tick labels under each box: 0.1/0.2/0.3

# y range must cover the differences, both sets of means and the zero line;
# extra space on top for the legend
ylim <- range(c(res$diff, summary_tab$PFC_true, summary_tab$PFC_OOB, 0))
ylim <- ylim + c(-0.05, 0.12) * diff(ylim)

par(mar = c(5, 5, 2, 1))   # RStudio: use the arrows in the Plots pane to go back to Figure 2
# First pass: draws axes and the boxes. The formula diff ~ miss + dataset
# orders the boxes by rate within dataset, matching the 'pos' positions.
boxplot(diff ~ miss + dataset, data = res, at = pos, names = lab,
        ylim = ylim, xlim = c(0.5, max(pos) + 0.5), col = "white",
        cex.axis = 0.8, las = 1, xlab = "",
        ylab = expression(err[true] - widehat(err)[OOB]))
abline(h = pretty(ylim), col = "grey85")          # grid lines ...
boxplot(diff ~ miss + dataset, data = res, at = pos, add = TRUE,  # ... behind boxes
        col = "white", axes = FALSE)
abline(h = 0, lty = 2)                  # zero line: OOB estimate equals true error
# Dataset names centred under each group of boxes
mtext(sapply(DATASETS, `[[`, "name"), side = 1, line = 2.5,
      at = sapply(seq_along(DATASETS) - 1, function(i) i * (k + 1) + (k + 1) / 2))
mtext("Missing values", side = 1, line = 3.8, cex = 0.8)

# Means of the true error (circle) and of the OOB estimate (plus), as in the paper
# (summary_tab is sorted by dataset, then rate, so it lines up with 'pos')
points(pos, summary_tab$PFC_true, pch = 16)
points(pos, summary_tab$PFC_OOB,  pch = 3, cex = 1.2)

legend("topleft", bg = "white", cex = 0.85,
       legend = c(expression(mean(widehat(err)[OOB])), expression(mean(err[true]))),
       pch = c(3, 16))
par(old_par)                            # restore the original graphics settings

cat("\nDone. Raw results saved to", file.path(OUT_DIR, "raw_results.csv"), "\n")
