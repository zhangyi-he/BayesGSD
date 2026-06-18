#' @title Seed sensitivity of the §3.5 calibrated designs
#' @description The calibrated three-look and five-look designs reported in
#'   Section 3.5 sit close to the regulatory boundary $\alpha \leq 2.5\%$. To
#'   pin the question "would these still be the chosen designs under a different
#'   RNG seed?", this script re-evaluates their binding and non-binding type I
#'   error rates and power at a small sweep of independent seeds, holding the
#'   simulation budget at R = 1,000,000 virtual trials per (seed, hypothesis)
#'   cell.
#'
#'   The two designs are on the 0.001-mesh grid of
#'   ADRENAL_Calibration_GridSearch.R / _Frontier.R and are stated here
#'   explicitly, each the maximum-power design subject to the type I error rate
#'   <= 2.5%, with common futility threshold q = 0.90.
#'
#'   Each seed is an independent pair of R = 10^6 simulations (H0 and H1) and is
#'   checkpointed as soon as it completes, so a session interruption loses at
#'   most one seed and a re-launch resumes from the cache. Env overrides:
#'   R_TRIALS (smaller R for a smoke test) and KEEP_CKPT (retain the checkpoint).
#'
#'   Output: prints per-seed binding type I, non-binding type I (DSMC override),
#'   and power for each design; reports mean / SD / range and the fraction of
#'   seeds at which the binding type I exceeds the $\alpha = 2.5\%$ target;
#'   writes `ADRENAL_Calibration_SeedSensitivity.rda` plus a sessionInfo dump.
#' @author Zhangyi He, Feng Yu, Suzie Cro, Laurent Billot

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

suppressPackageStartupMessages({ library(RBesT); library(parallel) })
source("./Code/Code v1.0/adabay_bern.R")

RDA_PATH <- file.path(OUTPUT_DIR, "ADRENAL_Calibration_SeedSensitivity.rda")
CKPT     <- file.path(OUTPUT_DIR, "ADRENAL_SeedSensitivity_ckpt.rda")

# N_MAX, make_looks() and the calibrated maximum-power designs
# (bayesgsd_p_stage(), bayesgsd_q()) are the single source of truth defined in
# adabay_bern_setup.R (sourced via the bootstrap above).

# The two headline Section 3.5 calibrated designs (source of truth:
# BAYESGSD_DESIGNS in adabay_bern_setup.R), each the maximum-power design subject
# to the type I error rate <= 2.5%: interim 0.999, final 0.977, q = 0.90.
designs <- list(
  list(K = 3L, p_stage = bayesgsd_p_stage(3L), q = BAYESGSD_DESIGNS$q_futility),
  list(K = 5L, p_stage = bayesgsd_p_stage(5L), q = BAYESGSD_DESIGNS$q_futility)
)
design_K <- vapply(designs, function(d) d$K, integer(1))

ALPHA_TARGET <- 0.025
union_looks  <- sort(unique(c(make_looks(3L), make_looks(5L))))

SEED_GRID <- 21:30
R_TRIALS  <- as.integer(Sys.getenv("R_TRIALS", unset = "1000000"))
P_CONTROL <- 0.33
DELTA_ALT <- -0.05

simulate_at_seed <- function(hypothesis, seed) {
  rates <- if (hypothesis == "H0") c(P_CONTROL, P_CONTROL) else c(P_CONTROL, P_CONTROL + DELTA_ALT)
  simulationSettings <- list(
    numberOfTrials   = R_TRIALS,
    numberOfSubjects = union_looks,
    allocationRatio  = 1,
    rates            = rates,
    hypothesis       = if (hypothesis == "H0") "null" else "alternative",
    seed             = seed
  )
  posteriorSettings <- list(
    effectThresholds = list(efficacy = list(size = 0, type = "absolute.risk"),
                            equivalence = NULL,
                            futility = list(size = 0, type = "absolute.risk")),
    alternative = "less"
  )
  runTrialMonitoring(simulationSettings, posteriorSettings,
                     priorSettings = initialisePriorSettings())
}

eval_design <- function(trialSimulation, K, p_stage, q, binding) {
  looks <- make_looks(K)
  if (binding) {
    futility_thr <- list(size = 0, type = "absolute.risk")
    fut_prob_mat <- matrix(1 - q, nrow = K, ncol = 1)
    fut_prob_mat[K, 1L] <- 0   # no futility at the final look (Section 2.1)
  } else {
    futility_thr <- NULL
    fut_prob_mat <- NULL
  }
  trialDesign <- list(
    numberOfSubjects = looks,
    effectThresholds = list(efficacy = list(size = 0, type = "absolute.risk"),
                            equivalence = NULL, futility = futility_thr),
    probabilityThresholds = list(efficacy = matrix(p_stage, nrow = K, ncol = 1),
                                 equivalence = NULL, futility = fut_prob_mat)
  )
  oc <- getOperatingCharacteristics(trialSimulation, trialDesign)
  err <- if (is.finite(oc$type1ErrorRate)) oc$type1ErrorRate else oc$type2ErrorRate
  list(error_rate = err, EN = oc$expSampleSize)
}

# ----------------------------------------------------------------------------
# Per-seed checkpointed sweep
# ----------------------------------------------------------------------------

cat("=== Calibration-design seed sensitivity ===\n")
for (d in designs)
  cat(sprintf("  K = %d : p_stage = (%s), q = %.2f\n", d$K,
              paste(formatC(d$p_stage, digits = 3, format = "f"), collapse = ", "), d$q))
cat(sprintf("Seeds: %s ; R = %s per cell\n\n",
            paste(SEED_GRID, collapse = ", "), formatC(R_TRIALS, big.mark = ",", format = "d")))

sig <- list(designs = designs, seeds = SEED_GRID, R = R_TRIALS,
            P = P_CONTROL, D = DELTA_ALT)
seed_results <- list()
if (file.exists(CKPT)) {
  e <- new.env(); load(CKPT, envir = e)
  if (identical(e$sig, sig)) {
    seed_results <- e$seed_results
    cat(sprintf("Resumed checkpoint: %d/%d seeds done\n", length(seed_results), length(SEED_GRID)))
  } else cat("Checkpoint signature changed; starting fresh.\n")
}

for (seed in SEED_GRID) {
  key <- as.character(seed)
  if (!is.null(seed_results[[key]])) { cat(sprintf("--- Seed %d cached ---\n", seed)); next }
  cat(sprintf("--- Seed %d ---\n", seed))
  tS <- proc.time()
  sim_H0 <- simulate_at_seed("H0", seed)
  sim_H1 <- simulate_at_seed("H1", seed)
  res <- list()
  for (d in designs) {
    bH0  <- eval_design(sim_H0, d$K, d$p_stage, d$q, binding = TRUE)
    bH1  <- eval_design(sim_H1, d$K, d$p_stage, d$q, binding = TRUE)
    nbH0 <- eval_design(sim_H0, d$K, d$p_stage, d$q, binding = FALSE)
    res[[as.character(d$K)]] <- list(type1_bind = bH0$error_rate, type1_nonbnd = nbH0$error_rate,
                                     power = 1 - bH1$error_rate, EN_H0 = bH0$EN, EN_H1 = bH1$EN)
    cat(sprintf("  K=%d  type1_bind = %.4f%%, type1_nonbnd = %.4f%%, power = %.4f%%\n",
                d$K, bH0$error_rate * 100, nbH0$error_rate * 100, (1 - bH1$error_rate) * 100))
  }
  seed_results[[key]] <- res
  save(seed_results, sig, file = CKPT)
  cat(sprintf("  seed %d done (%.0f s)\n", seed, (proc.time() - tS)[["elapsed"]]))
}

# ----------------------------------------------------------------------------
# Assemble result table from the cached seeds
# ----------------------------------------------------------------------------

result <- expand.grid(seed = SEED_GRID, K = design_K, stringsAsFactors = FALSE)
result$type1_bind <- NA_real_; result$type1_nonbnd <- NA_real_
result$power <- NA_real_; result$EN_H0_bind <- NA_real_; result$EN_H1_bind <- NA_real_
for (i in seq_len(nrow(result))) {
  r <- seed_results[[as.character(result$seed[i])]][[as.character(result$K[i])]]
  result$type1_bind[i] <- r$type1_bind; result$type1_nonbnd[i] <- r$type1_nonbnd
  result$power[i] <- r$power; result$EN_H0_bind[i] <- r$EN_H0; result$EN_H1_bind[i] <- r$EN_H1
}

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------

summarise <- function(df, K_pick) {
  sub <- df[df$K == K_pick, ]
  list(
    n = nrow(sub),
    type1_bind = list(mean = mean(sub$type1_bind), sd = sd(sub$type1_bind),
                      min = min(sub$type1_bind), max = max(sub$type1_bind),
                      n_over_25 = sum(sub$type1_bind > ALPHA_TARGET)),
    type1_nb = list(mean = mean(sub$type1_nonbnd), sd = sd(sub$type1_nonbnd),
                    min = min(sub$type1_nonbnd), max = max(sub$type1_nonbnd),
                    n_over_25 = sum(sub$type1_nonbnd > ALPHA_TARGET)),
    power = list(mean = mean(sub$power), sd = sd(sub$power),
                 min = min(sub$power), max = max(sub$power),
                 n_under_90 = sum(sub$power < 0.90))
  )
}

cat("\n=== Summary across", length(SEED_GRID), "seeds (alpha target = 0.025, power target = 0.90) ===\n")
for (d in designs) {
  s <- summarise(result, d$K)
  cat(sprintf("\nK = %d design (p = %s, q = %.2f):\n", d$K,
              paste(formatC(d$p_stage, digits = 3, format = "f"), collapse = ", "), d$q))
  cat(sprintf("  binding type I       mean = %.4f%%, sd = %.4f pp, range [%.4f%%, %.4f%%], cells over 2.5%% = %d/%d\n",
              s$type1_bind$mean * 100, s$type1_bind$sd * 100, s$type1_bind$min * 100, s$type1_bind$max * 100,
              s$type1_bind$n_over_25, s$n))
  cat(sprintf("  non-binding type I   mean = %.4f%%, sd = %.4f pp, range [%.4f%%, %.4f%%], cells over 2.5%% = %d/%d\n",
              s$type1_nb$mean * 100, s$type1_nb$sd * 100, s$type1_nb$min * 100, s$type1_nb$max * 100,
              s$type1_nb$n_over_25, s$n))
  cat(sprintf("  power                mean = %.4f%%, sd = %.4f pp, range [%.4f%%, %.4f%%], cells under 90%% = %d/%d\n",
              s$power$mean * 100, s$power$sd * 100, s$power$min * 100, s$power$max * 100,
              s$power$n_under_90, s$n))
}

save(result, designs, SEED_GRID, R_TRIALS, P_CONTROL, DELTA_ALT, ALPHA_TARGET, union_looks,
     file = RDA_PATH)
cat(sprintf("\nWrote %s\n", RDA_PATH))
if (!nzchar(Sys.getenv("KEEP_CKPT"))) unlink(CKPT)

bayesgsd_save_session("ADRENAL_Calibration_SeedSensitivity")
