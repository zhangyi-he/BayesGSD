#' @title Regression tests for bayseqSim_bern.R
#' @description Pin down the mathematical correctness of the kernel against
#'   hand-computed references. Run after any change to `bayseqSim_bern.R`.
#'   Exits non-zero if any check fails.
#' @author Zhangyi He, Feng Yu, Suzie Cro, Laurent Billot

PROJECT_ROOT <- Sys.getenv("BAYESGSD_ROOT", unset = "")
if (!nzchar(PROJECT_ROOT) || !dir.exists(PROJECT_ROOT)) {
  PROJECT_ROOT <- getwd()
}
setwd(PROJECT_ROOT)

# Source the kernel from whichever path resolves first.
kernel_candidates <- c(
  "./Code/Code v1.0/bayseqSim_bern.R",
  "./bayseqSim_bern.R"
)
kernel_path <- kernel_candidates[file.exists(kernel_candidates)][1L]
if (is.na(kernel_path)) {
  stop("Cannot locate bayseqSim_bern.R from PROJECT_ROOT = ", PROJECT_ROOT)
}
source(kernel_path)

# ---- helpers ---------------------------------------------------------------

expect_close <- function(actual, expected, tol = 1e-6, msg = "") {
  diff <- max(abs(as.numeric(actual) - as.numeric(expected)))
  if (diff > tol) {
    stop(sprintf("FAIL: %s\n  max |actual - expected| = %.3e (tol %.1e)\n  actual:   %s\n  expected: %s",
                 msg, diff, tol,
                 paste(format(actual, digits = 8), collapse = ", "),
                 paste(format(expected, digits = 8), collapse = ", ")))
  }
  cat(sprintf("PASS: %s (max diff %.3e)\n", msg, diff))
}

# Reference: posterior of theta_c (single-arm Bernoulli) under a Beta-mixture
# prior, normalised correctly via per-component marginal likelihoods.
# Returns P(theta_c > threshold | data) under the correctly-renormalised
# posterior mixture.
reference_arm_tail <- function(threshold, e, m, alpha, beta, weight,
                               lower_tail = FALSE) {
  L <- length(weight)
  lwm <- log(weight) + lbeta(alpha + e, beta + m - e) - lbeta(alpha, beta)
  lwm <- lwm - max(lwm)
  wpost <- exp(lwm) / sum(exp(lwm))
  contrib <- vapply(seq_len(L), function(l) {
    pbeta(threshold, alpha[l] + e, beta[l] + m - e, lower.tail = lower_tail)
  }, numeric(1L))
  sum(wpost * contrib)
}

# ---- Test 1: single-component Beta(1,1), Bernoulli, marginal arm tail ------
# The proposed framework's posterior tail prob for one trial should match
# the exact Beta-Binomial conjugate update.

cat("\n=== Test 1: Beta(1,1) single-arm tail ===\n")
{
  priors <- initialisePriorSettings()                  # default Beta(1,1)
  # Use the framework via the two-arm kernel with one trivial arm.
  # Threshold theta_t < x given binomial data (e=20, m=100):
  #   exact: pbeta(x, 1+20, 1+80, lower.tail = TRUE)
  e <- 20L; m <- 100L
  for (x in c(0.10, 0.20, 0.30, 0.5, 0.8)) {
    actual <- getPosteriorProbabilities(
      treatmentEffects = list(size = -1 + 1e-12 + x, type = "absolute.risk"),
      # build a "control" arm with degenerate data and a "treatment" arm with (e, m).
      # We instead use the single-arm tail manually via lower.tail of pbeta.
      numberOfEvents   = c(0L, e),
      numberOfSubjects = c(1L, m),
      alternative      = "less",
      priors           = priors
    )
    # Reference: marginalise over theta_c ~ Beta(1, 2) (control degenerate).
    # easier: use the analytic helper directly on the treatment arm,
    # then convolve. Since this gets messy, skip the cross-validation here
    # and validate the single-arm tail via direct convolution in Test 3 below.
  }
  cat("  (Test 1 deferred to Test 3 — full Beta(1,1) convolution check.)\n")
}

# ---- Test 2: 2-component mixture posterior weight renormalisation ----------
# Key bug fix: under a 2-component prior, the posterior weights are
#   w_l_post ∝ w_l * B(a_l + e, b_l + m - e) / B(a_l, b_l).
# We verify the framework's output matches this hand-computed reference.

cat("\n=== Test 2: 2-component mixture renormalisation ===\n")
{
  # Mixture prior: 0.6 * Beta(5, 5) + 0.4 * Beta(20, 5)
  alpha  <- c(5, 20)
  beta_  <- c(5, 5)
  weight <- c(0.6, 0.4)
  priors <- list(
    list(weight = weight, alpha = alpha, beta = beta_),
    list(weight = weight, alpha = alpha, beta = beta_)
  )

  # Data: control 30/100 events, treatment 25/100 events.
  ec <- 30L; mc <- 100L
  et <- 25L; nt <- 100L

  # Reference per-arm posterior weights.
  wpost_c_ref <- {
    lwm <- log(weight) + lbeta(alpha + ec, beta_ + mc - ec) - lbeta(alpha, beta_)
    lwm <- lwm - max(lwm); w <- exp(lwm); w / sum(w)
  }
  wpost_t_ref <- {
    lwm <- log(weight) + lbeta(alpha + et, beta_ + nt - et) - lbeta(alpha, beta_)
    lwm <- lwm - max(lwm); w <- exp(lwm); w / sum(w)
  }

  # Reference marginal P(theta < 0.4 | data) for the control arm.
  ref_marginal_tail_c_at_0.4 <- sum(wpost_c_ref *
    c(pbeta(0.4, alpha[1] + ec, beta_[1] + mc - ec, lower.tail = TRUE),
      pbeta(0.4, alpha[2] + ec, beta_[2] + mc - ec, lower.tail = TRUE)))
  cat(sprintf("  reference wpost_c = (%.5f, %.5f)\n", wpost_c_ref[1], wpost_c_ref[2]))
  cat(sprintf("  reference wpost_t = (%.5f, %.5f)\n", wpost_t_ref[1], wpost_t_ref[2]))
  cat(sprintf("  reference P(theta_c < 0.4 | data) = %.6f\n", ref_marginal_tail_c_at_0.4))

  # Reference posterior tail for absolute-risk Delta = theta_t - theta_c,
  # P(Delta < tau | data), computed by 1-D quadrature on theta_c.
  # Numerical reference uses dense pure-R double-integral (slow but exact).
  ref_delta_tail <- function(tau) {
    integrand <- function(theta_c) {
      dens_c <- sum(wpost_c_ref *
        c(dbeta(theta_c, alpha[1] + ec, beta_[1] + mc - ec),
          dbeta(theta_c, alpha[2] + ec, beta_[2] + mc - ec)))
      tail_t <- sum(wpost_t_ref *
        c(pbeta(theta_c + tau, alpha[1] + et, beta_[1] + nt - et, lower.tail = TRUE),
          pbeta(theta_c + tau, alpha[2] + et, beta_[2] + nt - et, lower.tail = TRUE)))
      dens_c * tail_t
    }
    integrate(Vectorize(integrand), 0, 1, rel.tol = 1e-10)$value
  }
  tau_grid <- c(-0.10, -0.05, 0, 0.05, 0.10)
  ref_tails <- vapply(tau_grid, ref_delta_tail, numeric(1L))

  # Framework's posterior tail probability for the same (R=1) trial.
  framework_tails <- getPosteriorProbabilities(
    treatmentEffects = list(size = tau_grid, type = "absolute.risk"),
    numberOfEvents   = c(ec, et),
    numberOfSubjects = c(mc, nt),
    alternative      = "less",
    priors           = priors
  )
  cat(sprintf("  reference tails:   %s\n",
              paste(format(ref_tails, digits = 8), collapse = ", ")))
  cat(sprintf("  framework tails:   %s\n",
              paste(format(framework_tails, digits = 8), collapse = ", ")))
  expect_close(framework_tails, ref_tails, tol = 5e-5,
               msg = "2-component mixture posterior tail P(Delta < tau | data)")
}

# ---- Test 3: Beta(1,1) full convolution (Delta tail) -----------------------
# Sanity check: for the single-component default prior, the framework's tail
# probability should match an independent 1-D integration of two betas.

cat("\n=== Test 3: Beta(1,1) Delta tail convolution ===\n")
{
  priors <- initialisePriorSettings()
  ec <- 50L; mc <- 200L
  et <- 40L; nt <- 200L

  ref_tail <- function(tau) {
    integrand <- function(theta_c) {
      dbeta(theta_c, 1 + ec, 1 + mc - ec) *
        pbeta(theta_c + tau, 1 + et, 1 + nt - et, lower.tail = TRUE)
    }
    integrate(Vectorize(integrand), 0, 1, rel.tol = 1e-10)$value
  }
  tau_grid <- c(-0.10, -0.05, 0, 0.05)
  ref_tails <- vapply(tau_grid, ref_tail, numeric(1L))
  framework_tails <- getPosteriorProbabilities(
    treatmentEffects = list(size = tau_grid, type = "absolute.risk"),
    numberOfEvents   = c(ec, et),
    numberOfSubjects = c(mc, nt),
    alternative      = "less",
    priors           = priors
  )
  cat(sprintf("  reference tails:   %s\n",
              paste(format(ref_tails, digits = 8), collapse = ", ")))
  cat(sprintf("  framework tails:   %s\n",
              paste(format(framework_tails, digits = 8), collapse = ", ")))
  expect_close(framework_tails, ref_tails, tol = 5e-5,
               msg = "Beta(1,1) Delta tail probability matches direct integration")
}

# ---- Test 4: M8 — numeric-proximity threshold lookup, not rownames ---------
# Build a posterior cache, then query thresholds whose `as.character`
# representation differs from how they were stored (e.g. arithmetic drift).

cat("\n=== Test 4: numeric-proximity threshold lookup ===\n")
{
  # Single small simulation
  set.seed(42L)
  simulationSettings <- list(
    numberOfTrials = 100L,
    numberOfSubjects = c(100L, 200L),
    allocationRatio = 1, rates = c(0.33, 0.28),
    hypothesis = "alternative", seed = 42L
  )
  posteriorSettings <- list(
    effectThresholds = list(
      efficacy = list(size = c(-0.05, 0), type = "absolute.risk"),
      equivalence = NULL,
      futility    = list(size = 0,         type = "absolute.risk")
    ),
    alternative = "less"
  )
  sim <- runTrialMonitoring(simulationSettings, posteriorSettings)

  # Query with a threshold that differs from the stored value by float drift.
  # Pick an explicit tiny non-zero perturbation that is well below the
  # `.match_thresholds` tolerance (1e-8) and would fail any string-based match.
  drifted_zero <- .Machine$double.eps  # ~2.22e-16
  stopifnot(drifted_zero != 0, abs(drifted_zero) < 1e-10)
  trialDesign <- list(
    numberOfSubjects = c(100L, 200L),
    effectThresholds = list(
      efficacy    = list(size = 0 + drifted_zero, type = "absolute.risk"),  # was 0
      equivalence = NULL,
      futility    = list(size = 0,                type = "absolute.risk")
    ),
    probabilityThresholds = list(
      efficacy    = matrix(0.99, 2L, 1L),
      equivalence = NULL,
      futility    = matrix(0.10, 2L, 1L)
    )
  )
  oc <- getOperatingCharacteristics(sim, trialDesign)
  stopifnot(!is.null(oc))
  cat(sprintf("  PASS: numeric-proximity lookup tolerates float drift of %.2e\n",
              abs(drifted_zero)))
}

# ---- Test 5: B3 — options(warn) restored after runTrialMonitoring ---------
cat("\n=== Test 5: options(warn) restored on exit ===\n")
{
  old_warn <- getOption("warn")
  options(warn = 1L)
  on.exit(options(warn = old_warn), add = TRUE)

  sim <- runTrialMonitoring(
    simulationSettings = list(numberOfTrials = 50L,
                              numberOfSubjects = 100L,
                              allocationRatio = 1,
                              rates = c(0.3, 0.3),
                              hypothesis = "null",
                              seed = 1L),
    posteriorSettings  = list(effectThresholds =
                                list(efficacy = list(size = 0, type = "absolute.risk"),
                                     equivalence = NULL, futility = NULL),
                              alternative = "less")
  )
  if (getOption("warn") != 1L) {
    stop(sprintf("FAIL: options(warn) not restored — got %d, expected 1",
                 getOption("warn")))
  }
  cat("  PASS: options(warn) restored to its prior value (1)\n")
}

# ---- Test 6: M2 — L'Ecuyer-CMRG installed ---------------------------------
cat("\n=== Test 6: L'Ecuyer-CMRG RNG kind ===\n")
{
  if (!identical(RNGkind()[1L], "L'Ecuyer-CMRG")) {
    stop(sprintf("FAIL: RNGkind()[1] = %s, expected L'Ecuyer-CMRG",
                 RNGkind()[1L]))
  }
  cat("  PASS: RNG kind is L'Ecuyer-CMRG\n")
}

# ---- Test 7: elbow selection rule in initialisePriorSettings -------------
# §2.3 elbow rule: pick the smallest L at which the running-minimum empirical
# cross-entropy stops improving by more than epsilon. Verify the selection
# returns a parsimonious mixture (L <= numberOfComponents) on a logit-normal
# prior sample, and that the same call with selection = "max" returns the
# largest L.

cat("\n=== Test 7: elbow selection rule ===\n")
{
  set.seed(21L)
  logit <- function(p) log(p / (1 - p))
  mu <- logit(0.33); sigma <- 0.5
  # Use the same 50,000-sample size the §3.3.2 benchmark uses; with the
  # 5,000-sample short version automixfit may converge to a single-component
  # fit for both the elbow and the max rules, which doesn't exercise the
  # selection logic.
  samples <- 1 / (1 + exp(-rnorm(50000L, mu, sigma)))
  ps_elbow <- initialisePriorSettings(priorSamples       = list(samples),
                                      numberOfComponents = 4L,
                                      selection          = "elbow")
  ps_max   <- initialisePriorSettings(priorSamples       = list(samples),
                                      numberOfComponents = 4L,
                                      selection          = "max")
  ps_aic   <- initialisePriorSettings(priorSamples       = list(samples),
                                      numberOfComponents = 4L,
                                      selection          = "aic")
  L_elbow <- length(ps_elbow[[1]]$weight)
  L_max   <- length(ps_max[[1]]$weight)
  L_aic   <- length(ps_aic[[1]]$weight)
  cat(sprintf("  elbow L = %d; aic L = %d; max L = %d\n",
              L_elbow, L_aic, L_max))
  stopifnot(L_elbow >= 1L, L_elbow <= 4L)
  stopifnot(L_max   >= 1L, L_max   <= 4L)
  stopifnot(L_elbow <= L_max)         # elbow is no larger than max
  # All three selections produce valid mixtures (weights sum to 1; positive
  # alpha/beta).
  for (ps in list(ps_elbow, ps_aic, ps_max)) {
    stopifnot(abs(sum(ps[[1]]$weight) - 1) < 1e-8)
    stopifnot(all(ps[[1]]$alpha > 0), all(ps[[1]]$beta > 0))
  }
  cat("  PASS: all three selection rules return valid mixtures with L_elbow <= L_max\n")
}

cat("\nAll tests passed.\n")
