#' @title Find calibrated Bayesian GSDs for the ADRENAL re-design (Section 3.5)
#' @description Loads the cached OC-behaviour simulations (R = 1,000,000 trials,
#'   union of look times for K = 2..10) and sweeps the precomputed posterior
#'   probabilities over two complementary grids:
#'
#'     (1) Constant per-stage threshold sweep, in which a single efficacy
#'         threshold p is applied at every stage. This is the fixed-threshold
#'         grid of Section 3.4 (whose operating-characteristic surfaces are
#'         supplement Figures S1-S5) extended to a uniform 0.001 p mesh.
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
#'   sweep the script lists the designs meeting alpha* = 2.5% (one-sided,
#'   binding) and 1 - beta* = 90%, reporting the minimum-E[N | H_1] member
#'   only as a reference corner of the feasible family. The Section 3.5
#'   calibration selects for maximum power (the most stringent interim
#'   threshold, 0.999), not minimum expected sample size.
#'
#'   For each manuscript-quoted Haybittle-Peto design (Section 3.5), the
#'   script additionally evaluates the type I error rate under the
#'   non-binding convention of Section 3.2, by re-running fast_oc() with the
#'   futility stop disabled (q -> 1.00 inside fast_oc, which short-circuits
#'   to fut_hit = FALSE). The non-binding Type I is reported alongside the
#'   binding Type I so the difference is explicit.
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
  library(RBesT); library(parallel)
})
source("./Code/Code v1.0/adabay_bern.R")

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

# Load each cache into its own environment so the H1 load() does not silently
# overwrite the H0 simulationSettings / trialSimulation we depend on below.
env_H0 <- new.env(parent = emptyenv())
env_H1 <- new.env(parent = emptyenv())
load(oc_h0, envir = env_H0)
load(oc_h1, envir = env_H1)

sim_H0      <- env_H0$trialSimulation
sim_H1      <- env_H1$trialSimulation
union_looks <- env_H0$simulationSettings$numberOfSubjects

# N_MAX and make_looks() are the single source of truth in adabay_bern_setup.R
# (sourced via the bootstrap above). Cross-validate that the loaded cache uses
# the same maximum sample size and the same look schedule across H0 and H1;
# otherwise the downstream sweeps would compare mismatched designs.
stopifnot(
  identical(union_looks, env_H1$simulationSettings$numberOfSubjects),
  max(union_looks) == N_MAX
)

ALPHA_TARGET <- 0.025
POWER_TARGET <- 0.90

K_grid    <- 1:9
Q_grid    <- c(0.85, 0.90, 0.95, 1.00)

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
# type I error ~ 4.15% (Table 1), power ~ 88.43% (Table 2), E[N|H1] ~ 1,956.
# ----------------------------------------------------------------------------

san <- evaluate_design_fast(K = 9, p_stage = rep(0.99, 9), q = 0.90, crit = "single")
cat(sprintf("Sanity check (K=9, p=0.99, q=0.90): type1=%.4f power=%.4f EN_H1=%.0f\n",
            san$type1, san$power, san$EN_H1))
# Hard assertion (not just a printout): if the cached posteriors or the kernel
# drift so the K=9/p=0.99/q=0.90 cell no longer matches the published Table 1/2
# values, fail the run loudly rather than passing silently with new numbers.
stopifnot(
  abs(san$type1 - 0.0415) < 2e-3,
  abs(san$power - 0.8843) < 2e-3,
  abs(san$EN_H1 - 1956)   < 50
)

# ============================================================================
# Sweep 1: constant per-stage threshold (fine p mesh)
# ============================================================================

# Efficacy-threshold resolution. Default to a uniform 0.001 mesh (Section 3.5);
# MESH_STEP may be set coarser via the environment for a quick smoke test.
MESH_STEP    <- as.numeric(Sys.getenv("MESH_STEP", unset = "0.001"))
P_CONST_GRID <- round(seq(0.900, 0.999, by = MESH_STEP), 3L)
CRIT_LEVELS  <- c("single", "dual")

# Parallelise across designs and checkpoint per (sweep, K) batch, so a session
# interruption loses at most one K-batch and a re-launch resumes from the cache.
n_cores <- if (exists("N_CORES")) max(1L, as.integer(N_CORES)) else max(1L, parallel::detectCores() - 2L)
options(BayesGSD.cores = 1L)

# Haybittle-Peto grids, on the same uniform MESH_STEP mesh as the constant sweep.
P_INTERIM_GRID <- round(seq(0.990, 0.999, by = MESH_STEP), 3L)
P_FINAL_GRID   <- round(seq(0.950, 0.999, by = MESH_STEP), 3L)

CKPT <- file.path(OUTPUT_DIR, "ADRENAL_GridSearch_ckpt.rda")
sig  <- list(P_CONST = P_CONST_GRID, P_INT = P_INTERIM_GRID, P_FIN = P_FINAL_GRID,
             Q = Q_grid, K = K_grid, CRIT = CRIT_LEVELS)
batches <- list(); compute_seconds <- 0
if (file.exists(CKPT)) {
  e <- new.env(); load(CKPT, envir = e)
  if (identical(e$sig, sig)) {
    batches <- e$batches; compute_seconds <- e$compute_seconds
    cat(sprintf("Resumed checkpoint: %d/%d K-batches done (%.1f min compute so far)\n",
                length(batches), 2L * length(K_grid), compute_seconds / 60))
  } else cat("Checkpoint grid signature changed; starting fresh.\n")
}
run_batch <- function(key, grid_df, eval_row) {
  if (!is.null(batches[[key]])) { cat(sprintf("  %s cached\n", key)); return(invisible()) }
  t0  <- proc.time()
  res <- mclapply(seq_len(nrow(grid_df)), eval_row, mc.cores = n_cores)
  bad <- !vapply(res, is.data.frame, logical(1))
  if (any(bad)) stop(sprintf("%s: %d design(s) failed under mclapply", key, sum(bad)))
  batches[[key]]  <<- do.call(rbind, res)
  compute_seconds <<- compute_seconds + (proc.time() - t0)[["elapsed"]]
  save(batches, sig, compute_seconds, file = CKPT)
  cat(sprintf("  %s done: %d designs (%.1f min cumulative compute)\n",
              key, nrow(batches[[key]]), compute_seconds / 60))
}

cat(sprintf("\nSweep 1: constant threshold (%d designs over a %d-value p mesh)\n",
            length(K_grid) * length(P_CONST_GRID) * length(Q_grid) * length(CRIT_LEVELS),
            length(P_CONST_GRID)))
for (K in K_grid) {
  g <- expand.grid(p = P_CONST_GRID, q = Q_grid, criterion = CRIT_LEVELS,
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  run_batch(sprintf("const_K%d", K), g, function(i) {
    oc <- evaluate_design_fast(K, rep(g$p[i], K), g$q[i], g$criterion[i])
    data.frame(K = K, p = g$p[i], q = g$q[i], criterion = g$criterion[i],
               schedule = "constant", type1 = oc$type1, power = oc$power,
               EN_H0 = oc$EN_H0, EN_H1 = oc$EN_H1, stringsAsFactors = FALSE)
  })
}

cat(sprintf("\nSweep 2: Haybittle-Peto-style (p_interim %d values, p_final %d values)\n",
            length(P_INTERIM_GRID), length(P_FINAL_GRID)))
for (K in K_grid) {
  if (K == 1L)
    g <- data.frame(p_interim = NA_real_, p_final = P_FINAL_GRID, q = NA_real_,
                    stringsAsFactors = FALSE)
  else
    g <- expand.grid(p_interim = P_INTERIM_GRID, p_final = P_FINAL_GRID, q = Q_grid,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  run_batch(sprintf("hp_K%d", K), g, function(i) {
    if (K == 1L) { p_stage <- g$p_final[i]; q_use <- 1.00 }
    else { p_stage <- c(rep(g$p_interim[i], K - 1L), g$p_final[i]); q_use <- g$q[i] }
    oc <- evaluate_design_fast(K, p_stage, q_use, "single")
    data.frame(K = K, p_interim = g$p_interim[i], p_final = g$p_final[i], q = g$q[i],
               criterion = "single", schedule = "haybittle-peto",
               type1 = oc$type1, power = oc$power,
               EN_H0 = oc$EN_H0, EN_H1 = oc$EN_H1, stringsAsFactors = FALSE)
  })
}

df_const <- do.call(rbind, batches[sprintf("const_K%d", K_grid)])
df_hp    <- do.call(rbind, batches[sprintf("hp_K%d", K_grid)])
cat(sprintf("\nAll sweeps complete: %d constant + %d HP designs, %.1f min total compute\n",
            nrow(df_const), nrow(df_hp), compute_seconds / 60))

# ============================================================================
# Report best designs from each sweep
# ============================================================================

report_best <- function(df, schedule_name) {
  elig <- df[df$type1 <= ALPHA_TARGET & df$power >= POWER_TARGET, ]
  cat(sprintf("\n=== %s sweep: designs meeting alpha <= %.3f and power >= %.2f: %d ===\n",
              schedule_name, ALPHA_TARGET, POWER_TARGET, nrow(elig)))
  if (nrow(elig) == 0) return(invisible(NULL))
  cat("\n-- min E[N|H1] reference (feasible family) overall --\n")
  print(elig[order(elig$EN_H1, elig$EN_H0), ][1, ])
  cat("\n-- min E[N|H1] reference by K --\n")
  best_by_K <- do.call(rbind, lapply(split(elig, elig$K), function(g) {
    g[order(g$EN_H1, g$EN_H0), ][1, ]
  }))
  print(best_by_K, row.names = FALSE)
  invisible(list(elig = elig, best_by_K = best_by_K))
}

best_const <- report_best(df_const, "Constant per-stage threshold")
best_hp    <- report_best(df_hp,    "Haybittle-Peto stage-specific")

# Also report the headline 3-look and 5-look designs explicitly, so the reader
# can verify the numbers without re-deriving them. For each design, report both
# the binding type I error rate (the operating characteristic used to calibrate
# the design against alpha* = 2.5%) and the non-binding type I error rate (the
# regulator-facing rate computed without enforcing the futility stop, per the
# convention of Section 3.2 of the main paper).
#
# The headline designs are the maximum-power frontier points defined once in
# adabay_bern_setup.R (BAYESGSD_DESIGNS), which is THE single source of truth
# for the Section 3.5 calibrated designs:
#   K = 3: p_{1:3} = (0.999, 0.999, 0.977), q = 0.90
#   K = 5: p_{1:5} = (0.999, 0.999, 0.999, 0.999, 0.977), q = 0.90
# Both are on the (P_INTERIM_GRID, P_FINAL_GRID) mesh of this sweep. The band
# verification (final threshold raised one grid step if the band supremum
# exceeds 2.5%) is performed by ADRENAL_Calibration_FrontierVerification.R.
nonbinding_type1 <- function(K, p_stage, crit = "single") {
  # Non-binding type I: futility removed entirely. Use q = 1.00, which
  # short-circuits fut_hit to FALSE inside fast_oc(); under H0 the
  # cumulative-stop logic then reduces to "any efficacy hit across the K
  # looks", matching the regulator's non-binding type I definition.
  fast_oc(pp_H0[[as.character(K)]], p_stage, q = 1.00, crit = crit)$efficacy_rate
}

hp_designs <- list(
  list(label = "calibrated (maximum power, q = 0.90)", q = BAYESGSD_DESIGNS$q_futility,
       p_3 = bayesgsd_p_stage(3L),
       p_5 = bayesgsd_p_stage(5L))
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

# Persist the per-K Haybittle-Peto batches into the durable cache too, so that
# ADRENAL_Calibration_FeasibleFamily.R can read them from a stable artefact
# rather than from the resume checkpoint that is unlinked on success below.
hp_batches <- batches[sprintf("hp_K%d", K_grid)]

save(df_const, df_hp, best_const, best_hp, df_hp_calibrated, hp_batches,
     P_CONST_GRID, P_INTERIM_GRID, P_FINAL_GRID, Q_grid, K_grid, CRIT_LEVELS,
     ALPHA_TARGET, POWER_TARGET, compute_seconds,
     file = file.path(OUTPUT_DIR, "ADRENAL_CalibratedDesign.rda"))
cat("\nSaved to", file.path(OUTPUT_DIR, "ADRENAL_CalibratedDesign.rda"), "\n")
if (!nzchar(Sys.getenv("KEEP_CKPT"))) unlink(CKPT)  # drop the resume checkpoint on success

# Persist sessionInfo and package versions alongside the cache.
bayesgsd_save_session("ADRENAL_Calibration_GridSearch")
