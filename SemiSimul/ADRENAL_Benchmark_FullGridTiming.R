#' @title Re-measure the full 438-design-grid timing for the proposed framework
#'   on the current 8-core machine
#' @description Runs `runTrialMonitoring` once for each hypothesis (H0 and H1)
#'   over the full union of sample sizes used in Section 3.4 (designs with
#'   between 1 and 9 equally spaced interim analyses, total 1,000,000 virtual
#'   trials, Beta(1,1) priors).  Reports the wall-clock cost.  Re-run this
#'   script after the OC simulation budget changes so the timing figures
#'   referenced in Sections 3.3 and 4 reflect the current configuration.
#' @author Zhangyi He, Feng Yu, Suzie Cro, Laurent Billot

# Project root, parallelism, INLA threads, sessionInfo helper.
# See Code/Code v1.0/bayseqSim_bern_setup.R for behaviour and override env vars.
local({
  bootstrap_root <- Sys.getenv("BAYESGSD_ROOT", unset = "")
  if (!nzchar(bootstrap_root) || !dir.exists(bootstrap_root)) {
    cur <- getwd()
    while (cur != dirname(cur)) {
      if (dir.exists(file.path(cur, "Code")) && dir.exists(file.path(cur, "Article"))) {
        bootstrap_root <- cur; break
      }
      cur <- dirname(cur)
    }
  }
  if (!nzchar(bootstrap_root) || !dir.exists(bootstrap_root)) {
    stop("Could not resolve BAYESGSD_ROOT. Set it before running.")
  }
  source(file.path(bootstrap_root, "Code", "Code v1.0", "bayseqSim_bern_setup.R"))
})



suppressPackageStartupMessages({
  library(RBesT)
  library(parallel)
})

source("./Code/Code v1.0/bayseqSim_bern.R")

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

cat("\n=== Full-grid runTrialMonitoring under H0 ===\n")
res_H0 <- run_full_grid_one("H0")
cat(sprintf("  elapsed = %.1f s\n", res_H0$elapsed))

cat("\n=== Full-grid runTrialMonitoring under H1 ===\n")
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
timing <- list(H0 = res_H0, H1 = res_H1, n_cores = N_CORES,
               n_designs = n_designs, R = res_H0$R,
               machine = Sys.info()[["nodename"]])
stopifnot(n_designs == 438L)   # sanity check vs the manuscript number
save(timing, file = file.path(OUTPUT_DIR, "ADRENAL_Benchmark_FullGridTiming.rda"))

cat("\n=== Full-grid summary ===\n")
cat(sprintf("H0: %.1f s, H1: %.1f s, on %d cores\n",
            res_H0$elapsed, res_H1$elapsed, N_CORES))
cat("Done.\n")


# Persist sessionInfo and package versions alongside the cache.
bayesgsd_save_session("ADRENAL_Benchmark_FullGridTiming")
