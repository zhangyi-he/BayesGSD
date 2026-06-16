#' @title Band verification of the maximum-power frontier (Section 3.5 / Suppl. S6)
#' @description The maximum-power calibration of Section 3.5 selects, at every
#'   number of looks K, the single-criterion Haybittle--Peto design that
#'   maximises power subject to the type I error rate not exceeding 2.5%. The
#'   selection is made on the grid-search cache at the single design null
#'   (vartheta = 0.33), where the cache type I estimate is optimistically biased;
#'   the operative constraint is the non-binding type I SUPREMUM over the
#'   control-rate band. This script re-evaluates each selected per-K design on
#'   its own analysis schedule over the band vartheta in
#'   {0.15, 0.20, 0.25, 0.30, 0.33, 0.36, 0.40, 0.45, 0.50} (delta = 0), under
#'   both the binding-futility convention and the non-binding
#'   (DSMC-always-overrides) upper bound, plus the power at the design
#'   alternative delta = -0.05, at R = 10^6 virtual trials per cell. The
#'   resulting band suprema and power are the frontier reported in Supplemental
#'   Material Section S6; any design whose non-binding band supremum exceeds 2.5%
#'   would have its final efficacy threshold raised one grid step (0.001).
#'
#'   Designs (interim 0.999 for K >= 2; K = 1 is the fixed-sample design;
#'   common futility threshold q = 0.90 with no futility rule at K = 1):
#'   final threshold 0.976 (K = 1, 2), 0.977 (K = 3..6), 0.978 (K = 7..10).
#'
#'   Each (design, rate) cell is an independent R = 10^6 simulation and is
#'   checkpointed as soon as it completes. Env overrides: R_TRIALS (smaller R for
#'   a smoke test) and KEEP_CKPT (retain the resume checkpoint on success).
#'
#'   Output: prints the per-K band suprema (binding and non-binding) and power;
#'   writes `ADRENAL_Calibration_FrontierVerification.rda` plus a sessionInfo dump.
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

RDA_PATH <- file.path(OUTPUT_DIR, "ADRENAL_Calibration_FrontierVerification.rda")
CKPT     <- file.path(OUTPUT_DIR, "ADRENAL_FrontierVerification_ckpt.rda")

# N_MAX, make_looks(), SEED_DEFAULT and the maximum-power frontier
# (BAYESGSD_DESIGNS, bayesgsd_p_stage(), bayesgsd_q()) are the single source of
# truth defined in adabay_bern_setup.R (sourced via the bootstrap above).
SEED       <- SEED_DEFAULT
R_TRIALS   <- as.integer(Sys.getenv("R_TRIALS", unset = "1000000"))

# Selected maximum-power final thresholds by K (interim 0.999 for K >= 2).
FINAL_BY_K <- BAYESGSD_DESIGNS$final_by_K
# q is carried as the common futility threshold for every K; eval_one() disables
# the futility rule at the final look (and there is no interim at K = 1), so the
# K = 1 q value never enters any operating characteristic.
mk <- function(K) {
  list(K = K, looks = make_looks(K), p_stage = bayesgsd_p_stage(K),
       q = BAYESGSD_DESIGNS$q_futility)
}
DESIGNS <- lapply(1:10, mk)

# One simulation pass feeds two threshold applications against the same cache
# (Lemma 1): binding futility (q) and futility-disabled (non-binding upper bound).
eval_one <- function(design, rates) {
  ss <- list(numberOfTrials = R_TRIALS, numberOfSubjects = design$looks, allocationRatio = 1,
             rates = rates, hypothesis = if (rates[1] == rates[2]) "null" else "alternative",
             seed = SEED)
  ps <- list(effectThresholds = list(efficacy = list(size = 0, type = "absolute.risk"),
             equivalence = NULL, futility = list(size = 0, type = "absolute.risk")),
             alternative = "less")
  ts <- runTrialMonitoring(ss, ps, priorSettings = initialisePriorSettings())
  K <- design$K
  md <- function(fut) list(numberOfSubjects = design$looks,
    effectThresholds = list(efficacy = list(size = 0, type = "absolute.risk"), equivalence = NULL,
      futility = list(size = 0, type = "absolute.risk")),
    probabilityThresholds = list(efficacy = matrix(design$p_stage, nrow = K, ncol = 1),
      equivalence = NULL, futility = fut))
  fb <- matrix(1 - design$q, nrow = K, ncol = 1); fb[K, 1L] <- 0
  fo <- matrix(0, nrow = K, ncol = 1)
  ob <- getOperatingCharacteristics(ts, md(fb)); on <- getOperatingCharacteristics(ts, md(fo))
  pk <- function(oc) if (is.finite(oc$type1ErrorRate)) oc$type1ErrorRate else oc$type2ErrorRate
  list(binding = pk(ob), nonbinding = pk(on), EN = ob$expSampleSize)
}

VARTHETA_GRID <- c(0.15, 0.20, 0.25, 0.30, 0.33, 0.36, 0.40, 0.45, 0.50)
DELTA_ALT     <- -0.05

cells <- list()
for (K in 1:10) {
  for (v in VARTHETA_GRID)
    cells[[length(cells) + 1L]] <- list(id = sprintf("t1_K%d_v%.2f", K, v), kind = "type1", K = K, x = v)
  cells[[length(cells) + 1L]] <- list(id = sprintf("pw_K%d", K), kind = "power", K = K, x = DELTA_ALT)
}

sig <- list(finals = FINAL_BY_K, V = VARTHETA_GRID, D = DELTA_ALT, R = R_TRIALS, seed = SEED)
results <- list()
if (file.exists(CKPT)) {
  e <- new.env(); load(CKPT, envir = e)
  if (identical(e$sig, sig)) { results <- e$results
    cat(sprintf("Resumed checkpoint: %d/%d cells done\n", length(results), length(cells)))
  } else cat("Checkpoint signature changed; starting fresh.\n")
}
cat(sprintf("Frontier band verification: %d cells at R = %s\n", length(cells),
            formatC(R_TRIALS, big.mark = ",", format = "d")))
for (cell in cells) {
  if (!is.null(results[[cell$id]])) { cat(sprintf("  %s cached\n", cell$id)); next }
  d <- DESIGNS[[cell$K]]; t0 <- proc.time()
  if (cell$kind == "type1") {
    r <- eval_one(d, c(cell$x, cell$x))
    results[[cell$id]] <- list(binding = r$binding, nonbinding = r$nonbinding, EN0 = r$EN)
  } else {
    r <- eval_one(d, c(0.33, 0.33 + cell$x))
    results[[cell$id]] <- list(power = 1 - r$binding, EN1 = r$EN)
  }
  save(results, sig, file = CKPT)
  cat(sprintf("  %s done (%.0f s)\n", cell$id, (proc.time() - t0)[["elapsed"]]))
}

# Assemble per-K band suprema and power.
frontier <- data.frame(K = 1:10, final = unname(FINAL_BY_K),
                       sup_bind = NA_real_, sup_nonbind = NA_real_, argsup = NA_real_,
                       power = NA_real_, EN_H0 = NA_real_, EN_H1 = NA_real_)
for (K in 1:10) {
  sb <- -Inf; sn <- -Inf; arg <- NA
  for (v in VARTHETA_GRID) {
    r <- results[[sprintf("t1_K%d_v%.2f", K, v)]]
    if (r$binding > sb) sb <- r$binding
    if (r$nonbinding > sn) { sn <- r$nonbinding; arg <- v }
  }
  en0 <- results[[sprintf("t1_K%d_v%.2f", K, 0.33)]]$EN0
  p <- results[[sprintf("pw_K%d", K)]]
  frontier[K, c("sup_bind", "sup_nonbind", "argsup", "power", "EN_H0", "EN_H1")] <-
    c(sb, sn, arg, p$power, en0, p$EN1)
}

cat("\n=== Maximum-power frontier: non-binding type I supremum over the band ===\n")
for (K in 1:10) {
  f <- frontier[K, ]
  cat(sprintf("K=%d final=%.3f : band-sup binding=%.3f%% non-binding=%.3f%% (@vartheta=%.2f) | power=%.3f%% | E(N|H0)=%.0f | E(N|H1)=%.0f %s\n",
              f$K, f$final, f$sup_bind * 100, f$sup_nonbind * 100, f$argsup,
              f$power * 100, f$EN_H0, f$EN_H1, if (f$sup_nonbind <= 0.025) "<=2.5% OK" else "OVER 2.5%"))
}

save(frontier, FINAL_BY_K, VARTHETA_GRID, DELTA_ALT, R_TRIALS, SEED, file = RDA_PATH)
cat(sprintf("\nWrote %s\n", RDA_PATH))
if (!nzchar(Sys.getenv("KEEP_CKPT"))) unlink(CKPT)
bayesgsd_save_session("ADRENAL_Calibration_FrontierVerification")
