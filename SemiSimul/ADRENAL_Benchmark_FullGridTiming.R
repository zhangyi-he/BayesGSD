#' @title Re-measure the full 438-design-grid timing for the proposed framework
#'   on the current 8-core machine
#' @description Times two stages of the calibration workflow under each
#'   hypothesis (H0 and H1):
#'     (1) the single simulation pass at the union of all candidate look times
#'         (`runTrialMonitoring` with 1,000,000 virtual trials, Beta(1,1)
#'         priors, equally spaced K = 2..10 schedules);
#'     (2) the 438-design grid sweep against that cache
#'         (`getOperatingCharacteristics` over K in 1..10, p in {0.90, 0.95,
#'         0.99}, q in {0.85, 0.90, 0.95, 1.00}, single- and dual-criterion
#'         efficacy, binding and non-binding futility -- the exact same axes
#'         as `ADRENAL_Evaluation_OperatingCharacteristics.R`).
#'   The sweep time pins the "approximately seven minutes" claim referenced
#'   in the abstract, Section 3.3, and Section 4 of the manuscript. Re-run
#'   this script after the OC simulation budget changes so the timing
#'   figures referenced in the prose reflect the current configuration.
#' @author Zhangyi He, Feng Yu, Suzie Cro, Laurent Billot

# Project root, parallelism, INLA threads, sessionInfo helper.
# See Code/Code v1.0/adabay_bern_setup.R for behaviour and override env vars.
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

OUTPUT_DIR <- "./Output/Output v1.0"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Same union-of-look-times as Section 3.4 (1..9 equally spaced interim
# analyses; max sample size 3,800).
all_look_sizes <- sort(unique(c(
  ceiling(seq(0, 3800, length.out =  3))[-1],
  ceiling(seq(0, 3800, length.out =  4))[-1],
  ceiling(seq(0, 3800, length.out =  5))[-1],
  ceiling(seq(0, 3800, length.out =  6))[-1],
  ceiling(seq(0, 3800, length.out =  7))[-1],
  ceiling(seq(0, 3800, length.out =  8))[-1],
  ceiling(seq(0, 3800, length.out =  9))[-1],
  ceiling(seq(0, 3800, length.out = 10))[-1],
  ceiling(seq(0, 3800, length.out = 11))[-1]
)))
cat("Full-grid look times (n =", length(all_look_sizes), "):", paste(all_look_sizes, collapse=","), "\n")

P_CONTROL <- 0.33
DELTA_ALT <- -0.05

run_full_grid_one <- function(hypothesis) {
  stopifnot(hypothesis %in% c("H0", "H1"))
  rates <- if (hypothesis == "H0") c(P_CONTROL, P_CONTROL) else c(P_CONTROL, P_CONTROL + DELTA_ALT)
  simulationSettings <- list(
    numberOfTrials   = 1000000L,
    numberOfSubjects = all_look_sizes,
    allocationRatio  = 1,
    rates            = rates,
    hypothesis       = if (hypothesis == "H0") "null" else "alternative",
    seed             = 21L
  )
  posteriorSettings <- list(
    effectThresholds = list(
      efficacy    = list(size = c(-0.05, 0.00), type = "absolute.risk"),
      equivalence = NULL,
      futility    = list(size = c(0.00),        type = "absolute.risk")
    ),
    alternative = "less"
  )
  priorSettings <- initialisePriorSettings()

  t0 <- proc.time()
  trialSimulation <- runTrialMonitoring(
    simulationSettings = simulationSettings,
    posteriorSettings  = posteriorSettings,
    priorSettings      = priorSettings
  )
  t_elapsed <- (proc.time() - t0)[["elapsed"]]
  list(
    elapsed = t_elapsed,
    n_looks = length(all_look_sizes),
    R       = simulationSettings$numberOfTrials,
    hypothesis = hypothesis
  )
}

cat("\n=== Stage 1a: simulation pass under H0 ===\n")
res_H0 <- run_full_grid_one("H0")
cat(sprintf("  elapsed = %.1f s\n", res_H0$elapsed))

cat("\n=== Stage 1b: simulation pass under H1 ===\n")
res_H1 <- run_full_grid_one("H1")
cat(sprintf("  elapsed = %.1f s\n", res_H1$elapsed))

# Derive the headline design-count from the grid product rather than
# hard-coding it: K = 1 fixed-sample contributes (efficacy thresholds in
# {0.90, 0.95, 0.99}) x (single + dual criterion) = 6 designs; the GSDs
# K in {2..10} contribute 9 x 3 (p) x 4 (q) x 2 (binding/non-binding) x
# 2 (single/dual) = 432.
n_p          <- 3L
n_q          <- 4L
n_K          <- 9L
n_crit       <- 2L
n_binding    <- 2L
n_designs    <- 1L * n_p * n_crit +              # K=1: q and binding inactive
                n_K * n_p * n_q * n_binding * n_crit
stopifnot(n_designs == 438L)   # sanity check vs the manuscript number

# ----------------------------------------------------------------------------
# Stage 2: 438-design grid sweep against the simulation cache (per hypothesis)
# ----------------------------------------------------------------------------
# Mirrors the sweep in ADRENAL_Evaluation_OperatingCharacteristics.R but is
# timing-focused: each design's getOperatingCharacteristics call is included
# in the wall-clock, with no figure rendering or downstream caching.

P_EFF_GRID  <- c(0.90, 0.95, 0.99)
P_FUT_LOW   <- c(0.15, 0.10, 0.05, 0.00)   # 1 - q for q in {0.85, 0.90, 0.95, 1.00}
K_GRID      <- 1:10
CRIT_LEVELS <- c("single", "dual")

make_looks  <- function(K) as.integer(ceiling(seq(0, 3800, length.out = K + 1L))[-1])

build_design <- function(K, p_eff, p_fut_low, crit, binding) {
  looks <- make_looks(K)
  if (binding) {
    futility_thr <- list(size = 0, type = "absolute.risk")
    fut_prob_mat <- matrix(p_fut_low, nrow = K, ncol = 1)
    fut_prob_mat[K, 1L] <- 0  # Section 2.1: no futility at the final look
  } else {
    futility_thr <- NULL
    fut_prob_mat <- NULL
  }
  if (crit == "single") {
    eff_thr      <- list(size = 0, type = "absolute.risk")
    eff_prob_mat <- matrix(p_eff, nrow = K, ncol = 1)
  } else {
    eff_thr      <- list(size = c(-0.05, 0.00), type = "absolute.risk")
    eff_prob_mat <- matrix(rep(c(0.50, p_eff), each = K), nrow = K, ncol = 2)
  }
  list(
    numberOfSubjects = looks,
    effectThresholds = list(efficacy = eff_thr, equivalence = NULL,
                            futility = futility_thr),
    probabilityThresholds = list(efficacy = eff_prob_mat, equivalence = NULL,
                                 futility = fut_prob_mat)
  )
}

# Helper: run the 438-design sweep for one (hypothesis, trialSimulation) cache
# and return the wall-clock plus a count of designs actually evaluated.
sweep_438 <- function(trialSimulation) {
  n_ok <- 0L
  t0 <- proc.time()
  for (K in K_GRID) {
    for (p_eff in P_EFF_GRID) {
      for (crit in CRIT_LEVELS) {
        if (K == 1L) {
          # K=1 fixed-sample: futility and binding/non-binding inactive (6 designs).
          d <- build_design(K, p_eff, 0, crit, binding = TRUE)
          oc <- getOperatingCharacteristics(trialSimulation, d)
          n_ok <- n_ok + 1L
        } else {
          for (p_fut in P_FUT_LOW) {
            for (binding in c(TRUE, FALSE)) {
              d <- build_design(K, p_eff, p_fut, crit, binding = binding)
              oc <- getOperatingCharacteristics(trialSimulation, d)
              n_ok <- n_ok + 1L
            }
          }
        }
      }
    }
  }
  t_elapsed <- (proc.time() - t0)[["elapsed"]]
  stopifnot(n_ok == n_designs)
  list(elapsed = t_elapsed, n_designs = n_ok)
}

# The simulation results from Stage 1 are not held in scope across the two
# hypotheses (run_full_grid_one returns only the timing). Re-run the
# simulation once more per hypothesis to get a fresh `trialSimulation` for
# the sweep step, and time that step independently. The Stage 1 timings
# above remain the canonical simulation-pass numbers; the second-run wall-
# clock is approximately equal and is reported here only as a sanity check.

run_sim_and_sweep <- function(hypothesis) {
  rates <- if (hypothesis == "H0") c(P_CONTROL, P_CONTROL) else c(P_CONTROL, P_CONTROL + DELTA_ALT)
  simulationSettings <- list(
    numberOfTrials   = 1000000L,
    numberOfSubjects = all_look_sizes,
    allocationRatio  = 1,
    rates            = rates,
    hypothesis       = if (hypothesis == "H0") "null" else "alternative",
    seed             = 21L
  )
  posteriorSettings <- list(
    effectThresholds = list(
      efficacy    = list(size = c(-0.05, 0.00), type = "absolute.risk"),
      equivalence = NULL,
      futility    = list(size = c(0.00),        type = "absolute.risk")
    ),
    alternative = "less"
  )
  priorSettings <- initialisePriorSettings()
  t_sim_start <- proc.time()
  trialSimulation <- runTrialMonitoring(
    simulationSettings = simulationSettings,
    posteriorSettings  = posteriorSettings,
    priorSettings      = priorSettings
  )
  t_sim <- (proc.time() - t_sim_start)[["elapsed"]]
  sweep <- sweep_438(trialSimulation)
  list(t_sim = t_sim, t_sweep = sweep$elapsed, n_designs = sweep$n_designs,
       t_total_this_hyp = t_sim + sweep$elapsed)
}

cat("\n=== Stage 2: 438-design grid sweep under H0 ===\n")
res2_H0 <- run_sim_and_sweep("H0")
cat(sprintf("  simulation pass    = %.1f s (sanity-check vs stage 1a = %.1f s)\n",
            res2_H0$t_sim, res_H0$elapsed))
cat(sprintf("  438-design sweep   = %.1f s\n", res2_H0$t_sweep))
cat(sprintf("  H0 total (sim+sweep) = %.1f s\n", res2_H0$t_total_this_hyp))

cat("\n=== Stage 2: 438-design grid sweep under H1 ===\n")
res2_H1 <- run_sim_and_sweep("H1")
cat(sprintf("  simulation pass    = %.1f s (sanity-check vs stage 1b = %.1f s)\n",
            res2_H1$t_sim, res_H1$elapsed))
cat(sprintf("  438-design sweep   = %.1f s\n", res2_H1$t_sweep))
cat(sprintf("  H1 total (sim+sweep) = %.1f s\n", res2_H1$t_total_this_hyp))

t_sim_total   <- res2_H0$t_sim    + res2_H1$t_sim
t_sweep_total <- res2_H0$t_sweep  + res2_H1$t_sweep
t_grand_total <- t_sim_total + t_sweep_total

cat("\n=== Combined 438-design grid wall-clock ===\n")
cat(sprintf("  simulation pass total (H0 + H1) = %.1f s (%.2f min)\n",
            t_sim_total, t_sim_total / 60))
cat(sprintf("  438-design sweep total (H0 + H1) = %.1f s (%.2f min)\n",
            t_sweep_total, t_sweep_total / 60))
cat(sprintf("  grand total                     = %.1f s (%.2f min)\n",
            t_grand_total, t_grand_total / 60))

timing <- list(
  H0 = res_H0,
  H1 = res_H1,
  sweep = list(
    H0 = res2_H0, H1 = res2_H1,
    sim_total_s   = t_sim_total,
    sweep_total_s = t_sweep_total,
    grand_total_s = t_grand_total,
    grand_total_min = t_grand_total / 60
  ),
  n_cores   = N_CORES,
  n_designs = n_designs,
  R         = res_H0$R,
  machine   = Sys.info()[["nodename"]]
)
save(timing, file = file.path(OUTPUT_DIR, "ADRENAL_Benchmark_FullGridTiming.rda"))

cat("\n=== Full-grid summary ===\n")
cat(sprintf("Simulation pass: %.1f s, 438-design sweep: %.1f s, total: %.1f s (%.2f min) on %d cores\n",
            t_sim_total, t_sweep_total, t_grand_total, t_grand_total / 60, N_CORES))
cat("Done.\n")


# Persist sessionInfo and package versions alongside the cache.
bayesgsd_save_session("ADRENAL_Benchmark_FullGridTiming")
