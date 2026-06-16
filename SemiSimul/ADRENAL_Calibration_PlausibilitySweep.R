#' @title Plausibility-set operating characteristics for the calibrated designs
#'   (Section 3.5)
#' @description The headline calibration evaluates the type I error rate at the
#'   single null point (vartheta, delta) = (0.33, 0) and the power at the single
#'   alternative (0.33, -0.05). A regulatory assessment expects the frequentist
#'   operating characteristics over a plausibility set rather than a single
#'   point, and in particular the supremum of the type I error rate over the
#'   nuisance control rate. This script re-evaluates the two Section 3.5
#'   calibrated designs (three-look and five-look, both with common futility
#'   threshold q = 0.90, on the 0.001-mesh grid) over a wide band of control
#'   event rates at the null
#'   and a band of treatment effects at the alternative, at R = 10^6 virtual
#'   trials per cell. The type I error rate is reported under both the
#'   binding-futility convention and the non-binding (DSMC-always-overrides,
#'   futility-disabled) upper bound against which the designs are calibrated;
#'   both are read from a single simulation per cell against the same cache
#'   (Lemma 1).
#'
#'   Type I sweep: rates = (vartheta, vartheta) for vartheta in
#'   {0.15, 0.20, 0.25, 0.30, 0.33, 0.36, 0.40, 0.45, 0.50} (delta = 0).
#'   Power sweep: rates = (0.33, 0.33 + delta) for delta in
#'   {-0.03, -0.04, -0.05, -0.06, -0.07} (vartheta = 0.33).
#'
#'   Each (design, rate) cell is an independent R = 10^6 simulation and is
#'   checkpointed as soon as it completes, so a session interruption loses at
#'   most one cell and a re-launch resumes from the cache. Set MESH-style env
#'   overrides R_TRIALS (smaller R for a smoke test) and KEEP_CKPT (retain the
#'   resume checkpoint on success) as needed.
#'
#'   Output: prints both sweeps and the type I supremum over the nuisance band;
#'   writes `ADRENAL_Calibration_PlausibilitySweep.rda` plus a sessionInfo dump.
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

RDA_PATH <- file.path(OUTPUT_DIR, "ADRENAL_Calibration_PlausibilitySweep.rda")
CKPT     <- file.path(OUTPUT_DIR, "ADRENAL_PlausibilitySweep_ckpt.rda")

# N_MAX, make_looks(), SEED_DEFAULT and the calibrated maximum-power designs
# (bayesgsd_p_stage(), bayesgsd_q()) are the single source of truth defined in
# adabay_bern_setup.R (sourced via the bootstrap above).
SEED       <- SEED_DEFAULT
R_TRIALS   <- as.integer(Sys.getenv("R_TRIALS", unset = "1000000"))

# The two headline Section 3.5 calibrated designs (source of truth:
# BAYESGSD_DESIGNS in adabay_bern_setup.R), each the maximum-power design subject
# to the type I error rate <= 2.5%: interim 0.999, final 0.977, q = 0.90.
DESIGNS <- list(
  K3 = list(K = 3L, looks = make_looks(3L),
            p_stage = bayesgsd_p_stage(3L), q = BAYESGSD_DESIGNS$q_futility),
  K5 = list(K = 5L, looks = make_looks(5L),
            p_stage = bayesgsd_p_stage(5L), q = BAYESGSD_DESIGNS$q_futility)
)

# Evaluate one (design, rates) cell. A single simulation pass feeds two
# threshold applications against the same cache (Lemma 1): the binding-futility
# convention used to calibrate the design, and the futility-disabled non-binding
# upper bound. Under the null both returned errors are type I rates; under the
# alternative `binding` is the type II rate (so 1 - binding is power under
# binding futility, matching the main-paper convention).
eval_one <- function(design, rates) {
  simulationSettings <- list(
    numberOfTrials   = R_TRIALS,
    numberOfSubjects = design$looks,
    allocationRatio  = 1,
    rates            = rates,
    hypothesis       = if (rates[1] == rates[2]) "null" else "alternative",
    seed             = SEED
  )
  posteriorSettings <- list(
    effectThresholds = list(efficacy = list(size = 0, type = "absolute.risk"),
                            equivalence = NULL,
                            futility = list(size = 0, type = "absolute.risk")),
    alternative = "less"
  )
  ts <- runTrialMonitoring(simulationSettings, posteriorSettings,
                           priorSettings = initialisePriorSettings())
  K <- design$K
  make_design <- function(fut) list(
    numberOfSubjects = design$looks,
    effectThresholds = list(efficacy = list(size = 0, type = "absolute.risk"),
                            equivalence = NULL,
                            futility = list(size = 0, type = "absolute.risk")),
    probabilityThresholds = list(efficacy = matrix(design$p_stage, nrow = K, ncol = 1),
                                 equivalence = NULL, futility = fut))
  fut_bind <- matrix(1 - design$q, nrow = K, ncol = 1); fut_bind[K, 1L] <- 0  # no futility at final look
  fut_off  <- matrix(0, nrow = K, ncol = 1)   # futility never fires => non-binding upper bound
  oc_b <- getOperatingCharacteristics(ts, make_design(fut_bind))
  oc_n <- getOperatingCharacteristics(ts, make_design(fut_off))
  pick <- function(oc) if (is.finite(oc$type1ErrorRate)) oc$type1ErrorRate else oc$type2ErrorRate
  list(binding = pick(oc_b), nonbinding = pick(oc_n), EN = oc_b$expSampleSize)
}

VARTHETA_GRID <- c(0.15, 0.20, 0.25, 0.30, 0.33, 0.36, 0.40, 0.45, 0.50)
DELTA_GRID    <- c(-0.03, -0.04, -0.05, -0.06, -0.07)

# Cell list: type I over the control rate (delta = 0), then power over the
# treatment effect (vartheta = 0.33).
cells <- list()
for (dn in names(DESIGNS)) for (v in VARTHETA_GRID)
  cells[[length(cells) + 1L]] <- list(id = sprintf("t1_%s_v%.2f", dn, v), kind = "type1", design = dn, x = v)
for (dn in names(DESIGNS)) for (d in DELTA_GRID)
  cells[[length(cells) + 1L]] <- list(id = sprintf("pw_%s_d%+.2f", dn, d), kind = "power", design = dn, x = d)

# Per-cell checkpoint: a session interruption loses at most one R = 10^6 cell.
sig <- list(designs = DESIGNS, V = VARTHETA_GRID, D = DELTA_GRID, R = R_TRIALS, seed = SEED)
results <- list()
if (file.exists(CKPT)) {
  e <- new.env(); load(CKPT, envir = e)
  if (identical(e$sig, sig)) {
    results <- e$results
    cat(sprintf("Resumed checkpoint: %d/%d cells done\n", length(results), length(cells)))
  } else cat("Checkpoint signature changed; starting fresh.\n")
}
cat(sprintf("Plausibility sweep: %d cells at R = %s\n", length(cells),
            formatC(R_TRIALS, big.mark = ",", format = "d")))
for (cell in cells) {
  if (!is.null(results[[cell$id]])) { cat(sprintf("  %s cached\n", cell$id)); next }
  design <- DESIGNS[[cell$design]]
  t0 <- proc.time()
  if (cell$kind == "type1") {
    r <- eval_one(design, c(cell$x, cell$x))
    results[[cell$id]] <- list(binding = r$binding, nonbinding = r$nonbinding)
  } else {
    r <- eval_one(design, c(0.33, 0.33 + cell$x))
    results[[cell$id]] <- list(power = 1 - r$binding)
  }
  save(results, sig, file = CKPT)
  cat(sprintf("  %s done (%.0f s)\n", cell$id, (proc.time() - t0)[["elapsed"]]))
}

# ---- Assemble the result tables from the cached cells ----
type1 <- data.frame(vartheta = VARTHETA_GRID,
                    K3_binding = NA_real_, K3_nonbinding = NA_real_,
                    K5_binding = NA_real_, K5_nonbinding = NA_real_)
for (i in seq_along(VARTHETA_GRID)) {
  v <- VARTHETA_GRID[i]
  a <- results[[sprintf("t1_K3_v%.2f", v)]]; b <- results[[sprintf("t1_K5_v%.2f", v)]]
  type1$K3_binding[i] <- a$binding; type1$K3_nonbinding[i] <- a$nonbinding
  type1$K5_binding[i] <- b$binding; type1$K5_nonbinding[i] <- b$nonbinding
}
power <- data.frame(delta = DELTA_GRID, K3 = NA_real_, K5 = NA_real_)
for (i in seq_along(DELTA_GRID)) {
  d <- DELTA_GRID[i]
  power$K3[i] <- results[[sprintf("pw_K3_d%+.2f", d)]]$power
  power$K5[i] <- results[[sprintf("pw_K5_d%+.2f", d)]]$power
}

# ---- Print ----
cat("\n=== Type I error across the null plausibility band (delta = 0) ===\n")
cat(sprintf("%-9s %12s %12s %12s %12s\n",
            "vartheta", "K3 binding", "K3 nonbind", "K5 binding", "K5 nonbind"))
for (i in seq_along(VARTHETA_GRID))
  cat(sprintf("%-9.2f %11.4f%% %11.4f%% %11.4f%% %11.4f%%\n", VARTHETA_GRID[i],
              type1$K3_binding[i] * 100, type1$K3_nonbinding[i] * 100,
              type1$K5_binding[i] * 100, type1$K5_nonbinding[i] * 100))

cat("\n=== Power across the alternative plausibility band (vartheta = 0.33) ===\n")
for (i in seq_along(DELTA_GRID))
  cat(sprintf("  delta=%.2f: K3 power = %.4f%%, K5 power = %.4f%%\n",
              DELTA_GRID[i], power$K3[i] * 100, power$K5[i] * 100))

sup_at <- function(col) {
  i <- which.max(type1[[col]])
  sprintf("%.4f%% @ vartheta=%.2f", type1[[col]][i] * 100, type1$vartheta[i])
}
cat("\n=== Type I supremum over the nuisance band [0.15, 0.50] ===\n")
cat(sprintf("  K3 binding sup = %s ; K3 non-binding sup = %s\n",
            sup_at("K3_binding"), sup_at("K3_nonbinding")))
cat(sprintf("  K5 binding sup = %s ; K5 non-binding sup = %s\n",
            sup_at("K5_binding"), sup_at("K5_nonbinding")))
cat(sprintf("Type I at the design null (vartheta=0.33): K3 binding %.4f%% non-binding %.4f%% ; K5 binding %.4f%% non-binding %.4f%%\n",
            type1$K3_binding[type1$vartheta == 0.33] * 100, type1$K3_nonbinding[type1$vartheta == 0.33] * 100,
            type1$K5_binding[type1$vartheta == 0.33] * 100, type1$K5_nonbinding[type1$vartheta == 0.33] * 100))
cat(sprintf("Power at the design alternative (delta=-0.05): K3 %.4f%%, K5 %.4f%%\n",
            power$K3[power$delta == -0.05] * 100, power$K5[power$delta == -0.05] * 100))

save(type1, power, VARTHETA_GRID, DELTA_GRID, DESIGNS, R_TRIALS, SEED, file = RDA_PATH)
cat(sprintf("\nWrote %s\n", RDA_PATH))
if (!nzchar(Sys.getenv("KEEP_CKPT"))) unlink(CKPT)
bayesgsd_save_session("ADRENAL_Calibration_PlausibilitySweep")
