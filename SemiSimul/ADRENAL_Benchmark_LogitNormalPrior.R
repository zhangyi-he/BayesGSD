#' @title Logit-normal-prior benchmark of the proposed framework against
#'   BATSS and adaptr (Section 3.3.2)
#' @description Replicates the head-to-head benchmark of Section 3.3 under an
#'   informative independent logit-normal prior on each arm
#'     logit(theta_c), logit(theta_t) ~ N(mu, sigma^2)
#'   with mu = logit(0.33), sigma = 0.5. The three reference implementations
#'   handle the non-conjugate prior in three different ways:
#'
#'     - Proposed framework: the logit-normal prior is approximated by a
#'       finite beta mixture via the prior-approximation step of
#'       Section 2.3 (RBesT::automixfit, forward Kullback-Leibler).
#'
#'     - BATSS: the independent logit-normal on the arms is translated into
#'       the closest INLA-representable independent normal prior on the
#'       regression coefficients (beta0, beta1) of the model y ~ group,
#'           beta0 = logit(theta_c)                  ~ N(mu, sigma^2)
#'           beta1 = logit(theta_t) - logit(theta_c) ~ N(0, 2 sigma^2),
#'       dropping the cross-covariance Cov(beta0, beta1) = -sigma^2, which
#'       INLA's `control.fixed` cannot express. This is the closest faithful
#'       prior approximation available within BATSS.
#'
#'     - adaptr: the logit-normal prior is supplied directly through a
#'       custom `fun_draws` that draws from each arm's posterior by
#'       inverse-CDF sampling on a fine grid of logit(theta), with no
#'       beta-mixture intermediate.
#'
#'   INLA is pinned to one thread per worker (num.threads = "1:1") so that
#'   BATSS and the proposed framework each use exactly 8 cores via a single
#'   layer of mclapply parallelism, matching the configuration of the
#'   headline benchmark in Section 3.3.
#'
#'   The script writes results into a single cache at
#'   ./Output/Output v1.0/ADRENAL_Benchmark_LogitNormalPrior.rda, with one
#'   list entry per (K, hypothesis) tag and slots
#'     $proposed_small, $proposed_big, $batss, $adaptr_small
#'   for the matched-budget proposed run, the high-precision proposed run,
#'   the matched-budget BATSS run, and the matched-budget adaptr run
#'   respectively. Each cell is filled lazily so that interrupted runs
#'   resume from the last completed cell.
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
  library(BATSS)
  library(INLA)
  library(RBesT)
  library(parallel)
  library(adaptr)
})

INLA::inla.setOption(num.threads = "1:1")
cat("INLA num.threads pinned to:", INLA::inla.getOption("num.threads"), "\n")

source("./Code/Code v1.0/bayseqSim_bern.R")

OUTPUT_DIR <- "./Output/Output v1.0"
RDA_PATH   <- file.path(OUTPUT_DIR, "ADRENAL_Benchmark_LogitNormalPrior.rda")

# ----------------------------------------------------------------------------
# Logit-normal prior on each arm (target prior for all three methods)
# ----------------------------------------------------------------------------

PRIOR_MU        <- log(0.33 / (1 - 0.33))
PRIOR_SIGMA     <- 0.5
N_PRIOR_SAMPLES <- 50000L

# Single simulation seed used across the prior-sample draws and all runs.
SEED <- 21L

logit_normal_samples <- function(n, mu, sigma) {
  z <- rnorm(n, mu, sigma)
  1 / (1 + exp(-z))
}

set.seed(SEED)
samples_c <- logit_normal_samples(N_PRIOR_SAMPLES, PRIOR_MU, PRIOR_SIGMA)
samples_t <- logit_normal_samples(N_PRIOR_SAMPLES, PRIOR_MU, PRIOR_SIGMA)

cat("Logit-normal prior on each arm: mu =", PRIOR_MU,
    " (logit), sigma =", PRIOR_SIGMA, "\n")
cat("Sample summary on response-rate scale (control arm):\n")
print(summary(samples_c))

# ----- Proposed framework: beta-mixture approximation -----------------------
priorSettings <- initialisePriorSettings(
  priorSamples       = list(samples_c, samples_t),
  numberOfComponents = 4
)
cat("\nBeta-mixture fit, control arm:\n");   print(priorSettings[[1]])
cat("\nBeta-mixture fit, treatment arm:\n"); print(priorSettings[[2]])

# ----- BATSS: marginal-matched independent normal priors on (beta0, beta1) --
BATSS_PRIOR_FIXED <- list(
  mean = list("(Intercept)" = PRIOR_MU,           group = 0),
  prec = list("(Intercept)" = 1/(PRIOR_SIGMA^2),  group = 1/(2*PRIOR_SIGMA^2))
)
cat("\nBATSS control.fixed translated from independent logit-normal:\n")
str(BATSS_PRIOR_FIXED)

# ----- adaptr: direct logit-normal posterior via inverse-CDF sampling -------
# Per-arm posterior on theta given binomial data (y events out of n) under
# logit(theta) ~ N(mu, sigma^2) is not in a named family, so we sample from
# it numerically by inverse-CDF on a fine grid of z = logit(theta) centred
# at the MLE with width set by the Fisher information. The grid is wide
# enough to cover the prior tails when n is small and tightens automatically
# as n grows.
custom_fun_draws <- function(arms, allocs, ys, control, n_draws) {
  draws <- list()
  for (a in arms) {
    ii <- which(allocs == a)
    n_events <- sum(ys[ii])
    n_total  <- length(ii)
    p_hat  <- (n_events + 0.5) / (n_total + 1)
    z_hat  <- log(p_hat / (1 - p_hat))
    se_hat <- 1 / sqrt(n_total * p_hat * (1 - p_hat))
    half_width <- max(8 * se_hat, 4 * PRIOR_SIGMA)
    z_grid <- seq(z_hat - half_width, z_hat + half_width, length.out = 1024L)
    p_grid <- 1 / (1 + exp(-z_grid))
    log_prior <- dnorm(z_grid, PRIOR_MU, PRIOR_SIGMA, log = TRUE)
    log_lik   <- n_events * log(p_grid) + (n_total - n_events) * log(1 - p_grid)
    log_post  <- log_prior + log_lik
    log_post  <- log_post - max(log_post)
    post      <- exp(log_post)
    cdf       <- cumsum(post) / sum(post)
    u   <- runif(n_draws)
    idx <- pmin(findInterval(u, cdf) + 1L, length(z_grid))
    draws[[a]] <- 1 / (1 + exp(-z_grid[idx]))
  }
  do.call(cbind, draws)
}

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

R_BENCH    <- 5000L
R_PROP_BIG <- 1000000L

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

  # Match the benchmark's parallel budget at every R; see comment in the
  # Beta(1,1) benchmark script.
  old_chunk <- getOption("BayesGSD.chunk_size")
  on.exit(options(BayesGSD.chunk_size = old_chunk), add = TRUE)
  options(BayesGSD.chunk_size = max(1L, as.integer(ceiling(R / N_CORES))))

  t0 <- proc.time()
  trialSimulation <- runTrialMonitoring(
    simulationSettings = simulationSettings,
    posteriorSettings  = posteriorSettings,
    priorSettings      = priorSettings        # logit-normal beta-mixture
  )
  t_sim <- (proc.time() - t0)[["elapsed"]]

  # Manuscript Eq. 3.2.3: futility evaluated only at interim looks
  # (k = 1, ..., K - 1); deactivate the final-stage row.
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
    prior         = "logit-normal beta-mixture"
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
    computation = "parallel", mc.cores = N_CORES,
    control.fixed = BATSS_PRIOR_FIXED        # logit-normal-matched prior
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
    prior         = "logit-normal-matched (independent N on beta0, beta1)"
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
  spec$fun_draws <- custom_fun_draws

  t0 <- proc.time()
  res <- run_trials(spec, n_rep = R, cores = N_CORES, base_seed = seed,
                    progress = NULL, export = c("PRIOR_MU", "PRIOR_SIGMA"))
  t_elapsed <- (proc.time() - t0)[["elapsed"]]
  perf <- check_performance(res, raw_ests = FALSE)
  prob_eff  <- perf$est[perf$metric == "prob_select_arm_T"]
  size_mean <- perf$est[perf$metric == "size_mean"]
  list(
    efficacy_prob = prob_eff,
    expected_n    = size_mean,
    elapsed       = t_elapsed,
    R             = R,
    seed          = seed,
    prior         = "logit-normal (direct, inverse-CDF on a fine grid)"
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
  save(bench_log, priorSettings, BATSS_PRIOR_FIXED, file = RDA_PATH)
}

for (design in DESIGNS) {
  for (hypothesis in c("H0", "H1")) {
    tag <- sprintf("%s_%s", design$name, hypothesis)
    cat(sprintf("\n=== %s ===\n", tag))

    cur <- bench_log[[tag]]
    if (is.null(cur)) cur <- list(design = design, hypothesis = hypothesis)

    if (is.null(cur$proposed_small)) {
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

    if (is.null(cur$proposed_big)) {
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

    if (is.null(cur$batss)) {
      cat(sprintf("  BATSS    @ R=%d (INLA pinned 1:1) ... ", R_BENCH))
      cur$batss <- tryCatch(run_batss(design, hypothesis, R = R_BENCH),
                            error = function(e) {
                              cat("ERROR:", conditionMessage(e), "\n")
                              list(efficacy_prob = NA_real_, futility_prob = NA_real_,
                                   expected_n = NA_real_, elapsed = NA_real_,
                                   R = R_BENCH,
                                   inla_threads = INLA::inla.getOption("num.threads"),
                                   prior = "logit-normal-matched (independent N on beta0, beta1)")
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

    if (is.null(cur$adaptr_small)) {
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
cat("\nDone. Logit-normal-prior benchmark complete.\n")

# Persist sessionInfo and package versions alongside the cache.
bayesgsd_save_session("ADRENAL_Benchmark_LogitNormalPrior")
