#' @title Uniform-prior benchmark of the proposed framework against BATSS and
#'   adaptr (Section 3.3.1)
#' @description Head-to-head benchmark of the proposed semi-simulation
#'   framework against BATSS and adaptr on the ADRENAL re-design under an
#'   independent Beta(1,1) prior on the per-arm response rates. All three
#'   methods share the same trial parameters: efficacy threshold p = 0.99,
#'   futility threshold q = 0.90, equal allocation, maximum sample size
#'   N = 3,800, and five designs with K in {1, 3, 5, 7, 9} equally spaced
#'   looks (K = 1 is the fixed-sample baseline).
#'
#'   The three reference implementations use the Beta(1,1) prior as follows:
#'
#'     - Proposed framework: independent Beta(1,1) on the per-arm rates via
#'       the framework default `initialisePriorSettings()`.
#'
#'     - BATSS: INLA's default weakly informative Gaussian priors on the
#'       regression coefficients of the model y ~ group, because BATSS tests
#'       the treatment-vs-control contrast on the coefficient of group rather
#'       than on the per-arm rates directly.
#'
#'     - adaptr: built-in Beta(1,1) posterior sampler via the default
#'       `setup_trial_binom()` (no custom `fun_draws` override).
#'
#'   INLA is pinned to one thread per worker (num.threads = "1:1") so that
#'   BATSS and the proposed framework each use exactly 8 cores via a single
#'   layer of mclapply parallelism.
#'
#'   Results are written into a single cache at
#'   ./Output/Output v1.0/ADRENAL_Benchmark_UniformPrior.rda, with one list
#'   entry per (K, hypothesis) tag and slots
#'     $proposed_small, $proposed_big, $batss, $adaptr_small
#'   for the matched-budget proposed run, the high-precision proposed run,
#'   the matched-budget BATSS run, and the matched-budget adaptr run
#'   respectively. Each cell is filled lazily so that interrupted runs
#'   resume from the last completed cell.
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
  library(BATSS)
  library(INLA)
  library(RBesT)
  library(parallel)
  library(adaptr)
})

INLA::inla.setOption(num.threads = "1:1")
cat("INLA num.threads pinned to:", INLA::inla.getOption("num.threads"), "\n")

source("./Code/Code v1.0/adabay_bern.R")

OUTPUT_DIR <- "./Output/Output v1.0"
RDA_PATH   <- file.path(OUTPUT_DIR, "ADRENAL_Benchmark_UniformPrior.rda")

# ----------------------------------------------------------------------------
# Trial parameters (match the Section 3.3 headline experiment)
# ----------------------------------------------------------------------------

P_CONTROL    <- 0.33
DELTA_ALT    <- -0.05
P_TREAT_H1   <- P_CONTROL + DELTA_ALT
N_MAX        <- 3800L

P_EFF        <- 0.99
Q_FUT_MS     <- 0.90
P_FUT_BATSS  <- 1 - Q_FUT_MS   # BATSS futility on P(treatment best | data) < b
P_FUT_LOW    <- 1 - Q_FUT_MS   # proposed futility on lower-tail posterior
P_INF_LOW    <- 1 - Q_FUT_MS   # adaptr inferiority on P(arm best | data) < b

R_BENCH    <- 5000L      # matched budget across all three methods
R_PROP_BIG <- 1000000L   # high-precision budget for the proposed framework
R_ADAPTR_BIG <- 1000000L # high-precision budget for adaptr (it is fast enough
                         # to reach 10^6; BATSS is not -- see the manuscript note)

# Single simulation seed used across all runs.
SEED <- 21L

make_looks <- function(K) {
  as.integer(ceiling(seq(from = 0, to = N_MAX, length.out = K + 1L))[-1])
}

DESIGNS <- list(
  list(name = "K1", K = 1L, looks = make_looks(1L)),
  list(name = "K3", K = 3L, looks = make_looks(3L)),
  list(name = "K5", K = 5L, looks = make_looks(5L)),
  list(name = "K7", K = 7L, looks = make_looks(7L)),
  list(name = "K9", K = 9L, looks = make_looks(9L))
)
cat("\nDesigns:\n")
for (d in DESIGNS) cat(sprintf("  %s K=%d looks=%s\n",
                                d$name, d$K, paste(d$looks, collapse=",")))

# ----------------------------------------------------------------------------
# Runners
# ----------------------------------------------------------------------------

logit <- function(p) log(p / (1 - p))

run_proposed <- function(design, hypothesis, R, seed) {
  rates <- if (hypothesis == "H0") c(P_CONTROL, P_CONTROL) else c(P_CONTROL, P_TREAT_H1)
  simulationSettings <- list(
    numberOfTrials = R, numberOfSubjects = design$looks,
    allocationRatio = 1, rates = rates,
    hypothesis = if (hypothesis == "H0") "null" else "alternative",
    seed = seed
  )
  posteriorSettings <- list(
    effectThresholds = list(
      efficacy = list(size = 0, type = "absolute.risk"),
      equivalence = NULL,
      futility = list(size = 0, type = "absolute.risk")
    ),
    alternative = "less"
  )
  priorSettings <- initialisePriorSettings()        # default Beta(1,1) on each arm

  # Match the benchmark's parallel budget: chunk_size <= R / N_CORES ensures
  # that mclapply is actually entered at all R values. Without this, the
  # default chunk_size = 50,000 leaves the small-R proposed run single-
  # threaded while BATSS and adaptr use all N_CORES, understating the
  # proposed framework's parallel speedup. Restore the option on exit.
  old_chunk <- getOption("BayesGSD.chunk_size")
  on.exit(options(BayesGSD.chunk_size = old_chunk), add = TRUE)
  options(BayesGSD.chunk_size = max(1L, as.integer(ceiling(R / N_CORES))))

  t0 <- proc.time()
  trialSimulation <- runTrialMonitoring(
    simulationSettings = simulationSettings,
    posteriorSettings  = posteriorSettings,
    priorSettings      = priorSettings
  )
  t_sim <- (proc.time() - t0)[["elapsed"]]

  # Manuscript decision criteria (Section 2.1): futility is evaluated only at interim looks
  # (k = 1, ..., K - 1), so we deactivate the final-stage futility row by
  # setting its threshold to 0 (pp < 0 is never true). At K = 1 there are
  # no interim looks, so the matrix has a single row of zeros and futility
  # is fully off, consistent with §3.2's "fixed-sample design has no
  # futility rule" remark.
  fut_thr <- matrix(P_FUT_LOW, nrow = design$K, ncol = 1)
  fut_thr[design$K, 1L] <- 0
  trialDesign <- list(
    numberOfSubjects = design$looks,
    effectThresholds = list(
      efficacy = list(size = 0, type = "absolute.risk"),
      equivalence = NULL,
      futility = list(size = 0, type = "absolute.risk")
    ),
    probabilityThresholds = list(
      efficacy = matrix(P_EFF, nrow = design$K, ncol = 1),
      equivalence = NULL,
      futility = fut_thr
    )
  )
  t0 <- proc.time()
  oc <- getOperatingCharacteristics(
    trialSimulation = trialSimulation,
    trialDesign     = trialDesign
  )
  t_oc <- (proc.time() - t0)[["elapsed"]]

  efficacy_prob <- if (hypothesis == "H0") oc$type1ErrorRate else (1 - oc$type2ErrorRate)
  list(
    efficacy_prob = efficacy_prob,
    expected_n    = oc$expSampleSize,
    elapsed_sim   = t_sim,
    elapsed_oc    = t_oc,
    elapsed       = t_sim + t_oc,
    R             = R,
    prior         = "Beta(1,1)"
  )
}

run_batss <- function(design, hypothesis, R) {
  beta_intercept <- logit(P_CONTROL)
  beta_group     <- if (hypothesis == "H0") 0 else logit(P_TREAT_H1) - logit(P_CONTROL)
  if (design$K == 1L) {
    interim <- list(recruited = numeric(0))
  } else {
    interim <- list(recruited = head(design$looks, design$K - 1L))
  }
  t0 <- proc.time()
  fit <- batss.glm(
    model = y ~ group, var = list(y = rbinom, group = alloc.balanced),
    var.control = list(y = list(size = 1L)),
    family = "binomial", link = "logit",
    beta = c(beta_intercept, beta_group),
    which = 2L, alternative = "less",
    R = R, N = N_MAX, interim = interim,
    prob0 = c(C = 0.5, T = 0.5),
    delta.eff = 0, delta.fut = 0,
    eff.arm = eff.arm.simple, eff.arm.control = list(b = P_EFF),
    fut.arm = fut.arm.simple, fut.arm.control = list(b = P_FUT_BATSS),
    H0 = (hypothesis == "H0"),
    computation = "parallel", mc.cores = N_CORES
  )
  t_elapsed <- (proc.time() - t0)[["elapsed"]]
  block <- if (hypothesis == "H0") fit$H0 else fit$H1
  list(
    efficacy_prob = block$efficacy$par$overall[1],
    futility_prob = block$futility$par$overall[1],
    expected_n    = mean(block$sample$C + block$sample$T),
    elapsed       = t_elapsed,
    R             = R,
    inla_threads  = INLA::inla.getOption("num.threads"),
    prior         = "INLA default weakly informative Gaussian on (Intercept, group)"
  )
}

run_adaptr <- function(design, hypothesis, R, seed) {
  rates <- if (hypothesis == "H0") c(P_CONTROL, P_CONTROL) else c(P_CONTROL, P_TREAT_H1)
  spec <- setup_trial_binom(
    arms              = c("C", "T"),
    true_ys           = rates,
    control           = "C",
    fixed_probs       = c(0.5, 0.5),
    data_looks        = design$looks,
    highest_is_best   = FALSE,
    superiority       = P_EFF,
    inferiority       = P_INF_LOW,
    n_draws           = 5000L
  )
  t0 <- proc.time()
  res <- run_trials(spec, n_rep = R, cores = N_CORES, base_seed = seed,
                    progress = NULL)
  t_elapsed <- (proc.time() - t0)[["elapsed"]]
  perf <- check_performance(res, raw_ests = FALSE)
  # `prob_superior` is two-sided (counts a trial whenever either arm crosses
  # the superiority threshold), so under H0 with a symmetric null its rate is
  # roughly twice the one-sided type I error rate we want. We use
  # `prob_select_arm_T`, the proportion of trials that end with the treatment
  # declared the winner, which matches the proposed framework's efficacy
  # probability.
  prob_eff  <- perf$est[perf$metric == "prob_select_arm_T"]
  size_mean <- perf$est[perf$metric == "size_mean"]
  list(
    efficacy_prob = prob_eff,
    expected_n    = size_mean,
    elapsed       = t_elapsed,
    R             = R,
    seed          = seed,
    prior         = "Beta(1,1) non-informative"
  )
}

# ----------------------------------------------------------------------------
# Main loop with incremental save
# ----------------------------------------------------------------------------

bench_log <- if (file.exists(RDA_PATH)) {
  e <- new.env(); load(RDA_PATH, envir = e)
  cat("Resuming from existing", RDA_PATH, "with",
      length(e$bench_log), "entries.\n")
  e$bench_log
} else list()

save_cache <- function() {
  save(bench_log, file = RDA_PATH)
}

# A slot is considered "cached" only if it exists AND was computed at the
# target R. If R differs (e.g., user bumped R_BENCH from 5,000 to 10,000),
# the slot is invalidated and re-run.
.cached_at_R <- function(slot, R_target) {
  !is.null(slot) && identical(as.integer(slot$R), as.integer(R_target))
}

# ----------------------------------------------------------------------------
# Warm-up (timing hygiene). The proposed framework's one-time initialisation
# -- Rcpp JIT compilation of the tail-probability kernel and the first
# mclapply worker fork -- is triggered lazily on the first call to
# runTrialMonitoring(). Left unmanaged it is charged in full to the first
# timed cell (K = 1, H0), inflating that single wall-clock figure to ~20 s and
# making the R = 5,000 run appear slower than the R = 10^6 run. We absorb it
# here with a small throwaway run so that every timed cell reflects steady-
# state cost. The result is discarded; the deterministic per-cell seeds are
# unaffected.
invisible(run_proposed(DESIGNS[[1]], "H0", R = 200L, seed = SEED))

for (design in DESIGNS) {
  for (hypothesis in c("H0", "H1")) {
    tag <- sprintf("%s_%s", design$name, hypothesis)
    cat(sprintf("\n=== %s ===\n", tag))

    cur <- bench_log[[tag]]
    if (is.null(cur)) cur <- list(design = design, hypothesis = hypothesis)

    if (!.cached_at_R(cur$proposed_small, R_BENCH)) {
      cat(sprintf("  proposed @ R=%d  ... ", R_BENCH))
      cur$proposed_small <- run_proposed(design, hypothesis, R = R_BENCH, seed = SEED)
      cat(sprintf("%.1f s; eff=%.4f; E[N]=%.0f\n",
                  cur$proposed_small$elapsed,
                  cur$proposed_small$efficacy_prob,
                  cur$proposed_small$expected_n))
      bench_log[[tag]] <- cur
      save_cache()
    } else {
      cat(sprintf("  proposed @ R=%d  ... cached\n", R_BENCH))
    }

    if (!.cached_at_R(cur$proposed_big, R_PROP_BIG)) {
      cat(sprintf("  proposed @ R=%d ... ", R_PROP_BIG))
      cur$proposed_big <- run_proposed(design, hypothesis, R = R_PROP_BIG, seed = SEED)
      cat(sprintf("%.1f s; eff=%.4f; E[N]=%.0f\n",
                  cur$proposed_big$elapsed,
                  cur$proposed_big$efficacy_prob,
                  cur$proposed_big$expected_n))
      bench_log[[tag]] <- cur
      save_cache()
    } else {
      cat(sprintf("  proposed @ R=%d ... cached\n", R_PROP_BIG))
    }

    if (!.cached_at_R(cur$batss, R_BENCH)) {
      cat(sprintf("  BATSS    @ R=%d (INLA pinned 1:1) ... ", R_BENCH))
      cur$batss <- tryCatch(run_batss(design, hypothesis, R = R_BENCH),
                            error = function(e) {
                              cat("ERROR:", conditionMessage(e), "\n")
                              list(efficacy_prob = NA_real_, futility_prob = NA_real_,
                                   expected_n = NA_real_, elapsed = NA_real_,
                                   R = R_BENCH,
                                   inla_threads = INLA::inla.getOption("num.threads"),
                                   prior = "INLA default weakly informative Gaussian on (Intercept, group)")
                            })
      cat(sprintf("%.1f s; eff=%.4f; E[N]=%.0f\n",
                  cur$batss$elapsed,
                  cur$batss$efficacy_prob,
                  cur$batss$expected_n))
      bench_log[[tag]] <- cur
      save_cache()
    } else {
      cat(sprintf("  BATSS    @ R=%d ... cached\n", R_BENCH))
    }

    if (!.cached_at_R(cur$adaptr_small, R_BENCH)) {
      cat(sprintf("  adaptr   @ R=%d (seed=%d) ... ", R_BENCH, SEED))
      cur$adaptr_small <- run_adaptr(design, hypothesis, R = R_BENCH, seed = SEED)
      cat(sprintf("%.1f s; eff=%.4f; E[N]=%.0f\n",
                  cur$adaptr_small$elapsed,
                  cur$adaptr_small$efficacy_prob,
                  cur$adaptr_small$expected_n))
      bench_log[[tag]] <- cur
      save_cache()
    } else {
      cat(sprintf("  adaptr   @ R=%d ... cached\n", R_BENCH))
    }

    # High-precision adaptr at R = 10^6 (the same budget as proposed_big).
    # adaptr's conjugate beta-binomial update is light enough to reach this
    # budget; BATSS is not (see manuscript note). Cached incrementally so a
    # long run resumes from the last completed cell.
    if (!.cached_at_R(cur$adaptr_big, R_ADAPTR_BIG)) {
      cat(sprintf("  adaptr   @ R=%d (seed=%d) ... ", R_ADAPTR_BIG, SEED))
      cur$adaptr_big <- run_adaptr(design, hypothesis, R = R_ADAPTR_BIG, seed = SEED)
      cat(sprintf("%.1f s; eff=%.4f; E[N]=%.0f\n",
                  cur$adaptr_big$elapsed,
                  cur$adaptr_big$efficacy_prob,
                  cur$adaptr_big$expected_n))
      bench_log[[tag]] <- cur
      save_cache()
    } else {
      cat(sprintf("  adaptr   @ R=%d ... cached\n", R_ADAPTR_BIG))
    }
  }
}

cat("\n=== Final summary ===\n")
for (tag in names(bench_log)) {
  r <- bench_log[[tag]]
  cat(sprintf(
    "%s: BATSS=%.4f  Prop_2k=%.4f  Prop_100k=%.4f  adaptr_2k=%.4f  E[N]: prop=%.0f BATSS=%.0f  t(s): prop=%.1f BATSS=%.1f adaptr=%.1f\n",
    tag,
    r$batss$efficacy_prob,
    r$proposed_small$efficacy_prob,
    r$proposed_big$efficacy_prob,
    r$adaptr_small$efficacy_prob,
    r$proposed_small$expected_n, r$batss$expected_n,
    r$proposed_small$elapsed, r$batss$elapsed, r$adaptr_small$elapsed
  ))
}
cat("\nDone. Uniform-prior benchmark complete.\n")

# Persist sessionInfo and package versions alongside the cache.
bayesgsd_save_session("ADRENAL_Benchmark_UniformPrior")
