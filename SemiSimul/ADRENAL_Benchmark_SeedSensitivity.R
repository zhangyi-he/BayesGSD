#' @title Multi-seed type I error estimate of the K=9 matched-budget proposed
#'   framework benchmark (Section 3.3.1)
#' @description Re-runs the K=9 row of Table 1 (single-criterion efficacy
#'   p = 0.99, binding futility q = 0.90) at the matched simulation budget
#'   R = 5,000 under H0, across 1000 independent RNG seeds (configurable via the
#'   N_SEEDS environment variable), to put empirical
#'   bars around the Proposed framework's single-seed point estimate
#'   (3.70% at seed 21). The R = 5,000 estimate sits 0.45 percentage points
#'   below the high-precision Proposed estimate (R = 10^6) of 4.15% and
#'   below the BATSS (4.28%) and adaptr (4.26%) matched-budget estimates,
#'   which the manuscript currently dismisses as "within combined Monte
#'   Carlo SE" without an empirical SE from a multi-seed re-run. This
#'   script provides that empirical SE.
#'
#'   The K = 9 cell is the largest-K row of Table 1, where the binomial
#'   Monte Carlo SE at R = 5,000 is also largest (0.28 percentage points
#'   at p ≈ 0.04). It is therefore the most sensitive cell to seed.
#'
#'   Output: prints per-seed type I error estimates and a summary; reports
#'   the empirical sd across seeds; computes a z-statistic for the
#'   seed-21 estimate against the high-precision Proposed reference and
#'   the BATSS / adaptr matched-budget estimates; writes the structured
#'   result to `ADRENAL_Benchmark_SeedSensitivity.rda` and a
#'   sessionInfo dump.
#' @author Zhangyi He, Feng Yu, Suzie Cro, Laurent Billot

# Defer project-root resolution and setup sourcing to the shared bootstrap.
local({
  root <- Sys.getenv("BAYESGSD_ROOT", unset = "")
  cur  <- if (nzchar(root) && dir.exists(root)) root else getwd()
  while (!file.exists(file.path(cur, "Code", "Code v1.0", "adabay_bern_bootstrap.R"))) {
    if (cur == dirname(cur))
      stop("Could not locate adabay_bern_bootstrap.R; set BAYESGSD_ROOT.")
    cur <- dirname(cur)
  }
  source(file.path(cur, "Code", "Code v1.0", "adabay_bern_bootstrap.R"))
})

suppressPackageStartupMessages({
  library(RBesT)
  library(parallel)
})
source("./Code/Code v1.0/adabay_bern.R")

# ----------------------------------------------------------------------------
# Inputs
# ----------------------------------------------------------------------------

RDA_PATH   <- file.path(OUTPUT_DIR, "ADRENAL_Benchmark_SeedSensitivity.rda")
CACHE_PATH <- file.path(OUTPUT_DIR, "ADRENAL_Benchmark_UniformPrior.rda")

# Reference values from the existing matched-budget benchmark cache
# (verified by Rscript inspection of bench_log$K9_H0 prior to this run).
# Reference values read from the matched-budget benchmark cache produced by
# `ADRENAL_Benchmark_UniformPrior.R` at the current R_BENCH. The values are
# loaded lazily at run time so this script stays in sync with the upstream
# cache even if its R_BENCH changes (e.g., 5,000 -> 10,000).
REF <- local({
  cache <- file.path(OUTPUT_DIR, "ADRENAL_Benchmark_UniformPrior.rda")
  if (!file.exists(cache)) {
    stop("Upstream benchmark cache missing: ", cache,
         ". Run ADRENAL_Benchmark_UniformPrior.R first.")
  }
  e <- new.env(parent = emptyenv()); load(cache, envir = e)
  cell <- e$bench_log$K9_H0
  list(
    proposed_seed21_Rbench = cell$proposed_small$efficacy_prob,
    proposed_seed21_R1M    = cell$proposed_big$efficacy_prob,
    batss_seed21_Rbench    = cell$batss$efficacy_prob,
    adaptr_seed21_Rbench   = cell$adaptr_small$efficacy_prob,
    R_bench_upstream       = cell$proposed_small$R
  )
})

# K = 9 design with equally spaced looks (matches Table 1 row K = 9).
K            <- 9L
N_MAX        <- 3800L
make_looks   <- function(K) as.integer(ceiling(seq(0, N_MAX, length.out = K + 1L))[-1])
LOOKS        <- make_looks(K)
P_EFF        <- 0.99
Q_FUT        <- 0.90
P_FUT_LOW    <- 1 - Q_FUT          # 0.10
P_CONTROL    <- 0.33

# Sweep settings.
SEED_GRID <- 21:(20L + as.integer(Sys.getenv("N_SEEDS", "1000")))   # 1000 seeds; seed 21 reproduces the cached benchmark value.
R_TRIALS  <- 5000L    # Matched-budget benchmark size (R_BENCH) of Section 3.3.1.

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

simulate_H0 <- function(seed) {
  simulationSettings <- list(
    numberOfTrials   = R_TRIALS,
    numberOfSubjects = LOOKS,
    allocationRatio  = 1,
    rates            = c(P_CONTROL, P_CONTROL),
    hypothesis       = "null",
    seed             = seed
  )
  posteriorSettings <- list(
    effectThresholds = list(
      efficacy    = list(size = 0, type = "absolute.risk"),
      equivalence = NULL,
      futility    = list(size = 0, type = "absolute.risk")
    ),
    alternative = "less"
  )
  priorSettings <- initialisePriorSettings()
  runTrialMonitoring(simulationSettings, posteriorSettings, priorSettings = priorSettings)
}

eval_K9 <- function(trialSimulation) {
  fut_prob_mat <- matrix(P_FUT_LOW, nrow = K, ncol = 1)
  fut_prob_mat[K, 1L] <- 0     # no futility at the final look (Section 2.1)
  trialDesign <- list(
    numberOfSubjects = LOOKS,
    effectThresholds = list(
      efficacy    = list(size = 0, type = "absolute.risk"),
      equivalence = NULL,
      futility    = list(size = 0, type = "absolute.risk")
    ),
    probabilityThresholds = list(
      efficacy    = matrix(P_EFF, nrow = K, ncol = 1),
      equivalence = NULL,
      futility    = fut_prob_mat
    )
  )
  oc <- getOperatingCharacteristics(trialSimulation, trialDesign)
  list(type1 = oc$type1ErrorRate, EN = oc$expSampleSize)
}

# ----------------------------------------------------------------------------
# Sweep
# ----------------------------------------------------------------------------

cat(sprintf("=== K = %d matched-budget seed sensitivity (R = %d, single-criterion efficacy p = %.2f, binding futility q = %.2f) ===\n",
            K, R_TRIALS, P_EFF, Q_FUT))
cat(sprintf("Seeds: %s\n\n", paste(SEED_GRID, collapse = ", ")))

result <- data.frame(
  seed  = SEED_GRID,
  type1 = NA_real_,
  EN    = NA_real_,
  stringsAsFactors = FALSE
)

t_start <- proc.time()
# The per-seed run is single-threaded (R = 5,000 is below the chunking
# threshold), so we parallelise across seeds instead, keeping each run
# single-threaded to avoid oversubscription.
options(BayesGSD.cores = 1L)
res_list <- mclapply(SEED_GRID, function(seed) {
  ts <- simulate_H0(seed)
  oc <- eval_K9(ts)
  c(type1 = oc$type1, EN = oc$EN)
}, mc.cores = max(1L, N_CORES))
ok <- vapply(res_list, function(x) is.numeric(x) && length(x) == 2L && all(is.finite(x)), logical(1))
if (!all(ok))
  stop(sprintf("mclapply failed for %d seed(s): %s", sum(!ok), paste(SEED_GRID[!ok], collapse = ", ")))
result$type1 <- vapply(res_list, function(x) unname(x[["type1"]]), numeric(1))
result$EN    <- vapply(res_list, function(x) unname(x[["EN"]]),    numeric(1))
t_total <- (proc.time() - t_start)[["elapsed"]]
cat(sprintf("Completed %d seeds in %.1f s (%.2f min) on %d core(s)\n",
            length(SEED_GRID), t_total, t_total / 60, max(1L, N_CORES)))

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------

emp_mean <- mean(result$type1)
emp_sd   <- sd(result$type1)
emp_min  <- min(result$type1)
emp_max  <- max(result$type1)

# Binomial SE under the population type I (use the high-precision reference 4.15%).
binom_se <- sqrt(REF$proposed_seed21_R1M * (1 - REF$proposed_seed21_R1M) / R_TRIALS)

R_label <- sprintf("R = %s", formatC(R_TRIALS, big.mark = ",", format = "d"))
if (!identical(as.integer(REF$R_bench_upstream), as.integer(R_TRIALS))) {
  warning(sprintf("R_TRIALS (%d) differs from upstream cache R (%d); cross-comparison values are at the upstream R.",
                  R_TRIALS, REF$R_bench_upstream))
}

cat(sprintf("\n=== Empirical seed-to-seed distribution at %s ===\n", R_label))
cat(sprintf("  mean      = %.4f%% (%.4f pp)\n", emp_mean * 100, emp_mean))
cat(sprintf("  sd        = %.4f pp  (binomial SE at p = %.4f: %.4f pp)\n",
            emp_sd * 100, REF$proposed_seed21_R1M, binom_se * 100))
cat(sprintf("  range     = [%.4f%%, %.4f%%]\n", emp_min * 100, emp_max * 100))

cat(sprintf("\n=== Cross-method comparison against the %s distribution ===\n", R_label))
compare <- function(label, value) {
  z <- (value - emp_mean) / emp_sd
  p_lo <- mean(result$type1 <= value)
  p_hi <- mean(result$type1 >= value)
  cat(sprintf("  %-40s %.4f%%  (z = %+.2f from empirical mean, pr[obs <= ref] = %.3f, pr[obs >= ref] = %.3f)\n",
              label, value * 100, z, p_lo, p_hi))
}
compare(sprintf("Proposed seed-21, %s:", R_label),  REF$proposed_seed21_Rbench)
compare("Proposed seed-21, R = 1,000,000:",         REF$proposed_seed21_R1M)
compare(sprintf("BATSS seed-21, %s:",    R_label),  REF$batss_seed21_Rbench)
compare(sprintf("adaptr seed-21, %s:",   R_label),  REF$adaptr_seed21_Rbench)

cat("\n=== Interpretation ===\n")
gap_Rb_vs_1M    <- REF$proposed_seed21_Rbench - REF$proposed_seed21_R1M
gap_batss_vs_1M <- REF$batss_seed21_Rbench    - REF$proposed_seed21_R1M
cat(sprintf("  Proposed seed-21 (%s) minus Proposed seed-21 (R=10^6) = %+.4f pp (%.2f empirical sd)\n",
            R_label, gap_Rb_vs_1M * 100, gap_Rb_vs_1M / emp_sd))
cat(sprintf("  BATSS seed-21 (%s) minus Proposed seed-21 (R=10^6)    = %+.4f pp (%.2f empirical sd)\n",
            R_label, gap_batss_vs_1M * 100, gap_batss_vs_1M / emp_sd))

# ----------------------------------------------------------------------------
# Persist
# ----------------------------------------------------------------------------

K_DESIGN <- K
save(result, K_DESIGN, LOOKS, R_TRIALS, SEED_GRID,
     P_EFF, Q_FUT, P_FUT_LOW, REF,
     emp_mean, emp_sd, emp_min, emp_max, binom_se,
     file = RDA_PATH)
cat(sprintf("\nWrote %s\n", RDA_PATH))

bayesgsd_save_session("ADRENAL_Benchmark_SeedSensitivity")
