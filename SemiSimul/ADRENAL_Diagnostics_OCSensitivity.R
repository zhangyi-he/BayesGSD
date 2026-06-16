#' @title Operating-characteristic sensitivity to the mixture order L
#'   (Appendix A.1, diagnostic iv)
#' @description Appendix A.1 reports that the beta-mixture approximation of the
#'   logit-normal prior matches the target prior's tail probabilities at the
#'   decision thresholds (Table A.1), but notes that prior-tail agreement is a
#'   necessary-but-not-sufficient check: what ultimately matters is whether the
#'   mixture order L moves the *operating characteristics* that drive design
#'   decisions. This script supplies that stronger diagnostic.
#'
#'   It takes the §3.5 calibrated three-look design (maximum power subject to
#'   the type I error rate <= 2.5%), replaces its
#'   Beta(1,1) prior with the informative logit-normal prior of §3.3.2
#'   approximated by a beta mixture of order L in {1, 2, 3}, and re-evaluates
#'   the type I error rate, power, and expected sample sizes under each L at
#'   R = 10^6 virtual trials. Because the simulated event counts share the
#'   single seed across all L, the only thing that varies between rows is the
#'   prior approximation, isolating the effect of L on the operating
#'   characteristics.
#'
#'   Output: prints the OC-by-L table; writes the structured result to
#'   `ADRENAL_Diagnostics_OCSensitivity.rda` plus a sessionInfo dump.
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

RDA_PATH <- file.path(OUTPUT_DIR, "ADRENAL_Diagnostics_OCSensitivity.rda")

# ----------------------------------------------------------------------------
# Logit-normal target prior (matches §3.3.2 / Appendix A.1)
# ----------------------------------------------------------------------------

PRIOR_MU        <- log(0.33 / (1 - 0.33))
PRIOR_SIGMA     <- 0.5
N_PRIOR_SAMPLES <- 50000L
SEED            <- SEED_DEFAULT   # project canonical seed (adabay_bern_setup.R)

logit_normal_samples <- function(n, mu, sigma) {
  z <- rnorm(n, mu, sigma)
  1 / (1 + exp(-z))
}

set.seed(SEED)
samples_c <- logit_normal_samples(N_PRIOR_SAMPLES, PRIOR_MU, PRIOR_SIGMA)
samples_t <- logit_normal_samples(N_PRIOR_SAMPLES, PRIOR_MU, PRIOR_SIGMA)

# ----------------------------------------------------------------------------
# §3.5 calibrated three-look design (maximum power s.t. type I <= 2.5%)
# ----------------------------------------------------------------------------

# N_MAX, make_looks() and the calibrated maximum-power design thresholds
# (bayesgsd_p_stage(), bayesgsd_q()) are the single source of truth defined in
# adabay_bern_setup.R (sourced via the bootstrap above).
K          <- 3L
LOOKS      <- make_looks(K)            # 1267, 2534, 3800
# Single-criterion efficacy thresholds of the headline three-look design (source
# of truth: BAYESGSD_DESIGNS in adabay_bern_setup.R): interim 0.999, final 0.977.
P_STAGE    <- bayesgsd_p_stage(K)
Q_FUT      <- BAYESGSD_DESIGNS$q_futility
P_FUT_LOW  <- 1 - Q_FUT                # 0.10; futility at interim looks only
P_CONTROL  <- 0.33
DELTA_ALT  <- -0.05
R_TRIALS   <- 1000000L

L_GRID <- 1:3

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

simulate <- function(hypothesis, priorSettings) {
  rates <- if (hypothesis == "H0") c(P_CONTROL, P_CONTROL) else c(P_CONTROL, P_CONTROL + DELTA_ALT)
  simulationSettings <- list(
    numberOfTrials   = R_TRIALS,
    numberOfSubjects = LOOKS,
    allocationRatio  = 1,
    rates            = rates,
    hypothesis       = if (hypothesis == "H0") "null" else "alternative",
    seed             = SEED
  )
  posteriorSettings <- list(
    effectThresholds = list(
      efficacy    = list(size = 0, type = "absolute.risk"),
      equivalence = NULL,
      futility    = list(size = 0, type = "absolute.risk")
    ),
    alternative = "less"
  )
  runTrialMonitoring(simulationSettings, posteriorSettings, priorSettings = priorSettings)
}

eval_design <- function(trialSimulation) {
  fut_prob_mat <- matrix(P_FUT_LOW, nrow = K, ncol = 1)
  fut_prob_mat[K, 1L] <- 0   # Section 2.1 decision rule: futility only at interim looks
  trialDesign <- list(
    numberOfSubjects = LOOKS,
    effectThresholds = list(
      efficacy    = list(size = 0, type = "absolute.risk"),
      equivalence = NULL,
      futility    = list(size = 0, type = "absolute.risk")
    ),
    probabilityThresholds = list(
      efficacy    = matrix(P_STAGE, nrow = K, ncol = 1),
      equivalence = NULL,
      futility    = fut_prob_mat
    )
  )
  oc <- getOperatingCharacteristics(trialSimulation, trialDesign)
  err <- if (is.finite(oc$type1ErrorRate)) oc$type1ErrorRate else oc$type2ErrorRate
  list(error_rate = err, EN = oc$expSampleSize)
}

# ----------------------------------------------------------------------------
# Sweep L = 1, 2, 3
# ----------------------------------------------------------------------------

cat(sprintf("=== OC sensitivity to mixture order L (calibrated three-look design) ===\n"))
cat(sprintf("Design: K=%d, looks=%s, p_stage=(%s), q=%.2f\n",
            K, paste(LOOKS, collapse = ","),
            paste(formatC(P_STAGE, digits = 3, format = "f"), collapse = ", "), Q_FUT))
cat(sprintf("Prior: logit-normal mu=logit(0.33), sigma=%.2f, approximated by beta mixture of order L.\n", PRIOR_SIGMA))
cat(sprintf("R = %s virtual trials per (L, hypothesis); seed = %d.\n\n", format(R_TRIALS, big.mark = ","), SEED))

result <- data.frame(L = L_GRID, ess = NA_real_, type1 = NA_real_,
                     power = NA_real_, EN_H0 = NA_real_, EN_H1 = NA_real_)

for (i in seq_along(L_GRID)) {
  L <- L_GRID[i]
  priorSettings <- initialisePriorSettings(
    priorSamples       = list(samples_c, samples_t),
    numberOfComponents = L,
    selection          = "max"     # force exactly L components
  )
  result$ess[i] <- priorSettings[[1]]$ESS

  simH0 <- simulate("H0", priorSettings)
  simH1 <- simulate("H1", priorSettings)
  ocH0  <- eval_design(simH0)
  ocH1  <- eval_design(simH1)

  result$type1[i] <- ocH0$error_rate
  result$power[i] <- 1 - ocH1$error_rate
  result$EN_H0[i] <- ocH0$EN
  result$EN_H1[i] <- ocH1$EN

  cat(sprintf("L=%d: type I = %.4f%%, power = %.4f%%, E[N|H0] = %.0f, E[N|H1] = %.0f\n",
              L, ocH0$error_rate * 100, (1 - ocH1$error_rate) * 100,
              ocH0$EN, ocH1$EN))
}

# ----------------------------------------------------------------------------
# Summary: spread across L (the diagnostic of interest)
# ----------------------------------------------------------------------------

cat("\n=== Spread across L (max - min) ===\n")
cat(sprintf("  type I : %.4f pp\n", (max(result$type1) - min(result$type1)) * 100))
cat(sprintf("  power  : %.4f pp\n", (max(result$power) - min(result$power)) * 100))
cat(sprintf("  E[N|H0]: %.1f patients\n", max(result$EN_H0) - min(result$EN_H0)))
cat(sprintf("  E[N|H1]: %.1f patients\n", max(result$EN_H1) - min(result$EN_H1)))
cat(sprintf("\nBinomial MC SE at R=%s, p~0.025: %.4f pp (reference)\n",
            format(R_TRIALS, big.mark = ","), sqrt(0.025 * 0.975 / R_TRIALS) * 100))

LHAT <- 2L  # within-tolerance selection from Appendix A.1
save(result, K, LOOKS, P_STAGE, Q_FUT, P_FUT_LOW, R_TRIALS, SEED,
     PRIOR_MU, PRIOR_SIGMA, L_GRID, LHAT, file = RDA_PATH)
cat(sprintf("\nWrote %s\n", RDA_PATH))

bayesgsd_save_session("ADRENAL_Diagnostics_OCSensitivity")
