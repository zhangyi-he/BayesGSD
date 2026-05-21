#' @title Find calibrated Bayesian GSDs for the ADRENAL re-design (Section 3.5)
#' @description Loads the cached OC-behaviour simulations (R = 1,000,000 trials,
#'   union of look times for K = 2..10) and sweeps the precomputed posterior
#'   probabilities over two complementary grids:
#'
#'     (1) Constant per-stage threshold sweep, in which a single efficacy
#'         threshold p is applied at every stage. This is the fixed-threshold
#'         grid of Section 3.4 (Figures 1-4) extended to a finer p mesh.
#'
#'     (2) Haybittle-Peto-style stage-specific sweep, in which the efficacy
#'         threshold takes a stringent value p_interim at each interim look
#'         (k = 1, ..., K - 1) and a separate near-nominal value p_final at
#'         the final analysis (k = K). This is the sweep described in
#'         Section 3.5 that delivers the calibrated three- and five-look
#'         designs.
#'
#'   The script uses a fast in-memory path that pre-extracts the per-stage
#'   posterior tail probability vectors from the cached simulations, then for
#'   each (K, p_stage, q) design computes type I error, power and expected
#'   sample size directly by elementwise comparison + cumulative-stop logic,
#'   bypassing the overhead of building the full OC data frame per design.
#'
#'   The grid covers K in {1, ..., 9} with single-criterion efficacy at the
#'   p values quoted in Section 3.4, plus dual-criterion for the constant
#'   sweep, and binding futility at q in {0.85, 0.90, 0.95, 1.00}. For each
#'   sweep the script reports the configuration that meets alpha* = 2.5%
#'   (one-sided, binding) and 1 - beta* = 90% with the smallest expected
#'   sample size under H_1 (and as a sanity check, the one minimising
#'   E[N | H_0]).
#'
#'   For each manuscript-quoted Haybittle-Peto design (Section 3.5), the
#'   script additionally evaluates the type I error rate under the
#'   non-binding convention of Section 3.2, by re-running fast_oc() with the
#'   futility stop disabled (q -> 1.00 inside fast_oc, which short-circuits
#'   to fut_hit = FALSE). The non-binding Type I is reported alongside the
#'   binding Type I so the difference is explicit.
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
  library(RBesT); library(parallel)
})
source("./Code/Code v1.0/bayseqSim_bern.R")

OUTPUT_DIR <- "./Output/Output v1.0"

# OC caches must exist before calibration. Fail fast with a clear message if
# the user hasn't run ADRENAL_Evaluation_OperatingCharacteristics.R yet.
oc_h0 <- file.path(OUTPUT_DIR, "ADRENAL_OperatingCharacteristics_H0.rda")
oc_h1 <- file.path(OUTPUT_DIR, "ADRENAL_OperatingCharacteristics_H1.rda")
missing <- c(if (!file.exists(oc_h0)) oc_h0,
             if (!file.exists(oc_h1)) oc_h1)
if (length(missing) > 0L) {
  stop(sprintf(
    "Required OC cache(s) missing:\n  %s\nRun ADRENAL_Evaluation_OperatingCharacteristics.R first.",
    paste(missing, collapse = "\n  ")
  ))
}

load(oc_h0)
sim_H0 <- trialSimulation
union_looks <- simulationSettings$numberOfSubjects
N_MAX <- max(union_looks)

load(oc_h1)
sim_H1 <- trialSimulation
# Cross-validate that the H0 and H1 caches share the same look schedule;
# otherwise downstream sweeps would be comparing mismatched designs.
stopifnot(identical(union_looks, simulationSettings$numberOfSubjects))

ALPHA_TARGET <- 0.025
POWER_TARGET <- 0.90

K_grid    <- 1:9
Q_grid    <- c(0.85, 0.90, 0.95, 1.00)

make_looks <- function(K) {
  as.integer(ceiling(seq(from = 0, to = N_MAX, length.out = K + 1L))[-1])
}

# ----------------------------------------------------------------------------
# Pre-extract the per-stage efficacy and futility tail probability vectors
# for each (K, sim) combination. For each K we identify the active stages,
# extract the corresponding posterior probability matrices, and store them as
# numeric matrices keyed by (K, hypothesis). The OC cache stores rows of
# `efficacy` as P(Delta < threshold) at thresholds (-0.05, 0); we keep both
# row indices so that the dual-criterion sweep can use both.
# ----------------------------------------------------------------------------

extract_pp_by_K <- function(sim) {
  cache_subj <- sim$virtualTrials$numberOfSubjects$control +
                sim$virtualTrials$numberOfSubjects$treatment
  pp_list <- sim$virtualTrials$posteriorProbabilities
  R       <- sim$virtualTrials$numberOfTrials
  out <- list()
  for (K in K_grid) {
    looks <- make_looks(K)
    active_stages <- which(cache_subj %in% looks)
    stopifnot(length(active_stages) == K)
    # Look up rows by numeric proximity against the cached `thresholds`
    # attribute, not by character rowname (which is fragile to float drift).
    pick_row <- function(stage_idx, kind, target) {
      mat <- pp_list[[stage_idx]][[kind]]
      stored <- attr(mat, "thresholds")
      j <- which(abs(stored - target) <= 1e-8)
      if (length(j) == 0L) {
        stop(sprintf("Threshold %g not found in stored %s cache (stage %d). Available: %s",
                     target, kind, stage_idx, paste(stored, collapse = ", ")))
      }
      mat[j[1L], ]
    }
    eff_lt_m05 <- vapply(active_stages, function(i) pick_row(i, "efficacy", -0.05), numeric(R))
    eff_lt_0   <- vapply(active_stages, function(i) pick_row(i, "efficacy",  0),    numeric(R))
    fut_lt_0   <- vapply(active_stages, function(i) pick_row(i, "futility",  0),    numeric(R))
    # vapply with vector FUN returns R x K with stages as columns
    out[[as.character(K)]] <- list(
      eff_lt_m05 = eff_lt_m05,    # R x K, P(Delta < -0.05)
      eff_lt_0   = eff_lt_0,      # R x K, P(Delta < 0)
      fut_lt_0   = fut_lt_0,      # R x K, P(Delta < 0) (also used for futility direction)
      R          = R,
      K          = K,
      looks      = looks
    )
  }
  out
}

cat("Pre-extracting per-stage posterior probability vectors ...\n")
t0 <- proc.time()
pp_H0 <- extract_pp_by_K(sim_H0)
pp_H1 <- extract_pp_by_K(sim_H1)
cat(sprintf("  elapsed = %.1f s\n", (proc.time() - t0)[["elapsed"]]))

# Free the big sim objects -- the extracted matrices hold what we need.
rm(sim_H0, sim_H1); gc(verbose = FALSE)

# ----------------------------------------------------------------------------
# Fast in-memory OC evaluation for one (K, p_stage, q, crit) design against
# a pre-extracted pp object. Returns list(type1_or_power = ..., EN = ...).
# The function does not distinguish H0 vs H1 -- the caller picks pp_H0 or
# pp_H1 and interprets the resulting "efficacy rate at final stage" as type I
# (under H0) or power (under H1).
# ----------------------------------------------------------------------------

fast_oc <- function(pp, p_stage, q, crit) {
  K <- pp$K
  R <- pp$R
  stopifnot(length(p_stage) == K)

  # Per-stage efficacy hits: R x K logical
  if (crit == "single") {
    # eff at stage k: P(Delta < 0) > p_stage[k]
    eff_hit <- sweep(pp$eff_lt_0, 2L, p_stage, FUN = ">")
  } else {
    # dual: eff at stage k iff P(Delta < 0) > p_stage[k] AND P(Delta < -0.05) > 0.5
    eff_hit_p   <- sweep(pp$eff_lt_0,   2L, p_stage, FUN = ">")
    eff_hit_inn <- pp$eff_lt_m05 > 0.5
    eff_hit     <- eff_hit_p & eff_hit_inn
  }
  # Per-stage futility hits (interim only): P(Delta > 0) > q
  #   <=> P(Delta < 0) < 1 - q
  # q = 1.00 deactivates futility (compared to 0 -- never true since P >= 0)
  fut_thr <- 1 - q
  if (K == 1L) {
    fut_hit <- matrix(FALSE, nrow = R, ncol = 1)
  } else {
    fut_hit <- pp$fut_lt_0 < fut_thr
    fut_hit[, K] <- FALSE   # futility only at interim looks 1..K-1
    if (q >= 1.0) fut_hit[, ] <- FALSE
  }

  # Cumulative-stop logic, vectorised. At each stage k a trial stops if either
  # eff_hit[, k] or fut_hit[, k] is TRUE; the trial's stop stage is the first
  # such k. Trials with no stop reach the final stage K. On ties (both eff and
  # fut at the same k) efficacy wins, because eff_hit[, k] is what is used for
  # the efficacy-rate calculation.
  combined <- eff_hit | fut_hit                  # R x K logical
  any_stop <- .rowSums(combined, R, K) > 0
  # max.col with ties.method = "first" returns the first TRUE column. For rows
  # with no TRUE, the returned value is arbitrary, so we overwrite with K.
  stop_stage <- max.col(combined, ties.method = "first")
  stop_stage[!any_stop] <- K
  # Efficacy declared iff eff_hit at the stop stage is TRUE (covers both
  # interim efficacy stops and final-look efficacy declarations).
  efficacy_declared <- eff_hit[cbind(seq_len(R), stop_stage)]
  efficacy_rate_at_K <- sum(efficacy_declared) / R
  # Expected sample size: average of `looks` at the trial's stop stage. Cast
  # to double to avoid integer overflow when R x max(looks) exceeds 2^31.
  EN <- mean(as.numeric(pp$looks[stop_stage]))

  list(efficacy_rate = efficacy_rate_at_K, EN = EN)
}

evaluate_design_fast <- function(K, p_stage, q, crit) {
  oc0 <- fast_oc(pp_H0[[as.character(K)]], p_stage, q, crit)
  oc1 <- fast_oc(pp_H1[[as.character(K)]], p_stage, q, crit)
  list(type1 = oc0$efficacy_rate,
       power = oc1$efficacy_rate,
       EN_H0 = oc0$EN,
       EN_H1 = oc1$EN)
}

# ----------------------------------------------------------------------------
# Sanity check: reproduce one cell from the head-to-head benchmark (K=9,
# constant p = 0.99, q = 0.90, single criterion). Expected from Table 1:
# type I error ~ 4.15%, power ~ 88.42%, E[N|H1] ~ 1,956.
# ----------------------------------------------------------------------------

san <- evaluate_design_fast(K = 9, p_stage = rep(0.99, 9), q = 0.90, crit = "single")
cat(sprintf("Sanity check (K=9, p=0.99, q=0.90): type1=%.4f power=%.4f EN_H1=%.0f\n",
            san$type1, san$power, san$EN_H1))

# ============================================================================
# Sweep 1: constant per-stage threshold (fine p mesh)
# ============================================================================

P_CONST_GRID <- c(0.90, 0.95, 0.97, 0.98, 0.99, 0.995, 0.998, 0.999)
CRIT_LEVELS  <- c("single", "dual")

cat("\nSweep 1: constant per-stage threshold p (",
    length(K_grid)*length(P_CONST_GRID)*length(Q_grid)*length(CRIT_LEVELS),
    "designs)\n")
t0 <- proc.time()
results_const <- vector("list", 0)
for (K in K_grid) for (p in P_CONST_GRID) for (q in Q_grid) for (crit in CRIT_LEVELS) {
  oc <- evaluate_design_fast(K, rep(p, K), q, crit)
  results_const[[length(results_const) + 1L]] <- c(
    list(K = K, p = p, q = q, criterion = crit, schedule = "constant"),
    oc
  )
}
df_const <- do.call(rbind, lapply(results_const, as.data.frame))
elapsed_const <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("  designs evaluated: %d  elapsed: %.1f s\n", nrow(df_const), elapsed_const))

# ============================================================================
# Sweep 2: Haybittle-Peto-style stage-specific threshold (p_interim, p_final)
# ============================================================================

# Stringent interim threshold + near-nominal final threshold. Grids are
# chosen so that the manuscript-quoted designs at (p_interim, p_final) =
# (0.995, 0.980) for K = 3 and (0.997, 0.980) for K = 5, both at q = 0.85,
# and the script-historical alternatives at (0.998, 0.9775) at q = 0.90,
# all lie on the mesh.
P_INTERIM_GRID <- c(0.990, 0.995, 0.997, 0.998, 0.999)
P_FINAL_GRID   <- c(0.950, 0.960, 0.970, 0.975, 0.9775, 0.980, 0.985,
                    0.990, 0.995, 0.998, 0.999)

cat("\nSweep 2: Haybittle-Peto-style (p_interim, p_final)\n")
cat(sprintf("  p_interim grid (%d values): %s\n",
            length(P_INTERIM_GRID),
            paste(sprintf("%.4f", P_INTERIM_GRID), collapse = ",")))
cat(sprintf("  p_final   grid (%d values): %s\n",
            length(P_FINAL_GRID),
            paste(sprintf("%.4f", P_FINAL_GRID), collapse = ",")))
t0 <- proc.time()
results_hp <- vector("list", 0)
for (K in K_grid) {
  if (K == 1L) {
    # Fixed-sample design: only the final-look threshold applies; q has no
    # effect because there are no interim looks. Loop over p_final only and
    # single criterion only (consistent with manuscript HP designs).
    for (p_fin in P_FINAL_GRID) {
      oc <- evaluate_design_fast(K, p_fin, q = 1.00, crit = "single")
      results_hp[[length(results_hp) + 1L]] <- c(
        list(K = K, p_interim = NA_real_, p_final = p_fin,
             q = NA_real_, criterion = "single", schedule = "haybittle-peto"),
        oc
      )
    }
  } else {
    for (p_int in P_INTERIM_GRID) for (p_fin in P_FINAL_GRID) for (q in Q_grid) {
      p_stage <- c(rep(p_int, K - 1L), p_fin)
      oc <- evaluate_design_fast(K, p_stage, q, crit = "single")
      results_hp[[length(results_hp) + 1L]] <- c(
        list(K = K, p_interim = p_int, p_final = p_fin,
             q = q, criterion = "single", schedule = "haybittle-peto"),
        oc
      )
    }
  }
}
df_hp <- do.call(rbind, lapply(results_hp, as.data.frame))
elapsed_hp <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("  designs evaluated: %d  elapsed: %.1f s\n", nrow(df_hp), elapsed_hp))

# ============================================================================
# Report best designs from each sweep
# ============================================================================

report_best <- function(df, schedule_name) {
  elig <- df[df$type1 <= ALPHA_TARGET & df$power >= POWER_TARGET, ]
  cat(sprintf("\n=== %s sweep: designs meeting alpha <= %.3f and power >= %.2f: %d ===\n",
              schedule_name, ALPHA_TARGET, POWER_TARGET, nrow(elig)))
  if (nrow(elig) == 0) return(invisible(NULL))
  cat("\n-- best (min E[N|H1]) overall --\n")
  print(elig[order(elig$EN_H1, elig$EN_H0), ][1, ])
  cat("\n-- best (min E[N|H1]) by K --\n")
  best_by_K <- do.call(rbind, lapply(split(elig, elig$K), function(g) {
    g[order(g$EN_H1, g$EN_H0), ][1, ]
  }))
  print(best_by_K, row.names = FALSE)
  invisible(list(elig = elig, best_by_K = best_by_K))
}

best_const <- report_best(df_const, "Constant per-stage threshold")
best_hp    <- report_best(df_hp,    "Haybittle-Peto stage-specific")

# Also report the manuscript-quoted 3-look and 5-look HP designs explicitly,
# so the reader can verify the numbers without re-deriving them. For each
# design, report both the binding type I error rate (the operating
# characteristic used to calibrate the design against alpha* = 2.5%) and the
# non-binding type I error rate (the regulator-facing rate computed without
# enforcing the futility stop, per the convention of Section 3.2 of the main
# paper).
#
# Two candidate threshold combinations are evaluated for each K, reflecting a
# minor inconsistency between the manuscript text and an earlier version of
# this script: the manuscript text at Section 3.5 reports
#   K = 3: p_{1:3} = (0.995, 0.995, 0.980), q = 0.85
#   K = 5: p_{1:5} = (0.997, 0.997, 0.997, 0.997, 0.980), q = 0.85
# and an older version of this script reported the alternative
#   K = 3: p_{1:3} = (0.998, 0.998, 0.9775), q = 0.90
#   K = 5: p_{1:5} = (0.998, 0.998, 0.998, 0.998, 0.9775), q = 0.90.
# Both sets are on the (P_INTERIM_GRID, P_FINAL_GRID) mesh and so are
# directly comparable. The reader can verify which set the manuscript
# numbers actually correspond to from the output below.
nonbinding_type1 <- function(K, p_stage, crit = "single") {
  # Non-binding type I: futility removed entirely. Use q = 1.00, which
  # short-circuits fut_hit to FALSE inside fast_oc(); under H0 the
  # cumulative-stop logic then reduces to "any efficacy hit across the K
  # looks", matching the regulator's non-binding type I definition.
  fast_oc(pp_H0[[as.character(K)]], p_stage, q = 1.00, crit = crit)$efficacy_rate
}

hp_designs <- list(
  list(label = "manuscript (q = 0.85)", q = 0.85,
       p_3 = c(0.995, 0.995, 0.980),
       p_5 = c(0.997, 0.997, 0.997, 0.997, 0.980)),
  list(label = "script-historical (q = 0.90)", q = 0.90,
       p_3 = c(0.998, 0.998, 0.9775),
       p_5 = c(0.998, 0.998, 0.998, 0.998, 0.9775))
)

hp_report <- vector("list", 0)
cat("\n=== Manuscript-quoted HP designs (single criterion) ===\n")
for (hp in hp_designs) {
  cat(sprintf("\n-- %s --\n", hp$label))
  for (K in c(3, 5)) {
    p_stage  <- if (K == 3L) hp$p_3 else hp$p_5
    oc       <- evaluate_design_fast(K, p_stage, q = hp$q, crit = "single")
    type1_nb <- nonbinding_type1(K, p_stage)
    cat(sprintf("K = %d, p_{1:%d} = (%s), q = %.2f:\n",
                K, K, paste(sprintf("%.4f", p_stage), collapse = ", "), hp$q))
    cat(sprintf("  Type I (binding,     q = %.2f) = %.4f\n",  hp$q,  oc$type1))
    cat(sprintf("  Type I (non-binding, q -> 1.00) = %.4f\n", type1_nb))
    cat(sprintf("  Power = %.4f, E[N|H0] = %.0f, E[N|H1] = %.0f\n",
                oc$power, oc$EN_H0, oc$EN_H1))
    hp_report[[length(hp_report) + 1L]] <- data.frame(
      label        = hp$label,
      K            = K,
      q            = hp$q,
      p_stage      = paste(sprintf("%.4f", p_stage), collapse = ","),
      type1_bind   = oc$type1,
      type1_nonbnd = type1_nb,
      power        = oc$power,
      EN_H0        = oc$EN_H0,
      EN_H1        = oc$EN_H1,
      stringsAsFactors = FALSE
    )
  }
}
df_hp_calibrated <- do.call(rbind, hp_report)
cat("\n-- summary table --\n")
print(df_hp_calibrated, row.names = FALSE)

save(df_const, df_hp, best_const, best_hp, df_hp_calibrated,
     P_CONST_GRID, P_INTERIM_GRID, P_FINAL_GRID, Q_grid, K_grid, CRIT_LEVELS,
     ALPHA_TARGET, POWER_TARGET,
     elapsed_const, elapsed_hp,
     file = file.path(OUTPUT_DIR, "ADRENAL_CalibratedDesign.rda"))
cat("\nSaved to", file.path(OUTPUT_DIR, "ADRENAL_CalibratedDesign.rda"), "\n")

# Persist sessionInfo and package versions alongside the cache.
bayesgsd_save_session("ADRENAL_Calibration_GridSearch")
