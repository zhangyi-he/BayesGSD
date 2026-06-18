#' @title Shared environment setup for the BayesGSD reproduction scripts
#' @description Resolves the project root, configures parallelism, pins INLA
#'   threads if INLA is available, and provides a `bayesgsd_save_session()`
#'   helper that each script calls to persist its `sessionInfo()` and full
#'   package version inventory alongside the run-specific output cache.
#'
#'   Usage: each ADRENAL_*.R driver should `source()` this file early, before
#'   any setwd or library() calls of its own. After sourcing, the following
#'   are in scope:
#'     PROJECT_ROOT, OUTPUT_DIR, MANUSCRIPT_DIR, SUPPLEMENT_DIR,
#'     CODE_DIR, N_CORES, SEED_DEFAULT, N_MAX, make_looks(), and the
#'     calibrated Section 3.5 design constants (BAYESGSD_DESIGNS) described
#'     in the "Shared trial constants" block below, plus
#'     bayesgsd_save_session().
#' @author Zhangyi He, Feng Yu, Suzie Cro, Laurent Billot

# --- Project root resolution ------------------------------------------------
# Priority: (1) BAYESGSD_ROOT env var if set; (2) auto-detect by walking up
# from getwd() looking for the Code/ and Article/ siblings; (3) stop loudly.
.resolve_project_root <- function() {
  root <- Sys.getenv("BAYESGSD_ROOT", unset = "")
  if (nzchar(root) && dir.exists(root)) return(root)
  cur <- getwd()
  while (cur != dirname(cur)) {
    if (dir.exists(file.path(cur, "Code")) &&
        dir.exists(file.path(cur, "Article"))) return(cur)
    cur <- dirname(cur)
  }
  stop(
    "BAYESGSD_ROOT environment variable is not set and the project root ",
    "could not be auto-detected from the current working directory (",
    getwd(), "). Set it explicitly, e.g.:\n\n",
    "  R> Sys.setenv(BAYESGSD_ROOT = \"/absolute/path/to/repo\")\n",
    "  $ export BAYESGSD_ROOT=/absolute/path/to/repo  # before Rscript\n\n",
    "or run this script from inside the repository directory tree."
  )
}

PROJECT_ROOT   <- .resolve_project_root()
setwd(PROJECT_ROOT)

OUTPUT_DIR     <- file.path(".", "Output",   "Output v1.0")
MANUSCRIPT_DIR <- file.path(".", "Article",  "Manuscript", "Manuscript v1.0")
SUPPLEMENT_DIR <- file.path(MANUSCRIPT_DIR, "Supplemental Material")
CODE_DIR       <- file.path(".", "Code",    "Code v1.0")

# --- Parallelism ------------------------------------------------------------
N_CORES <- as.integer(Sys.getenv("BAYESGSD_CORES", unset = "8"))
options(BayesGSD.cores = N_CORES, mc.cores = N_CORES)

# --- INLA thread pinning (no-op if INLA is not loaded) ---------------------
if (requireNamespace("INLA", quietly = TRUE)) {
  INLA::inla.setOption(num.threads = "1:1")
}

# --- Shared trial constants (single source of truth) ------------------------
# These are defined once here so every ADRENAL_*.R driver reads the same
# values. Do NOT re-declare N_MAX, make_looks(), the project seed, or the
# Section 3.5 calibrated-design thresholds inside individual drivers: change
# them here and they propagate everywhere, which is what keeps the cache and
# the calibration consistent.

# Project canonical RNG seed. Mirrors the default in runTrialMonitoring()
# (adabay_bern.R); every driver uses this unless it deliberately sweeps seeds.
SEED_DEFAULT <- 21L

# Maximum total sample size of the ADRENAL re-design (Section 3.4 skeleton),
# and the equally-spaced look schedule used throughout: K analyses placed at
# ceil(seq(0, N_MAX, length.out = K + 1))[-1], i.e. K evenly spaced looks with
# the final look at N_MAX. make_looks(K) returns the per-look cumulative sizes.
N_MAX <- 3800L
make_looks <- function(K) {
  as.integer(ceiling(seq(from = 0, to = N_MAX, length.out = K + 1L))[-1])
}

# ---- Section 3.5 calibrated designs: the AUTHORITATIVE maximum-power frontier.
# With N_MAX and the futility boundary fixed, the calibration spends the type I
# budget and MAXIMISES power. This single object is the source of truth for the
# per-K frontier and for the two headline designs (K = 3 and K = 5); every
# downstream driver (FrontierSelection, FrontierVerification, FeasibleFamily, PlausibilitySweep,
# SeedSensitivity, OCSensitivity, Summary) reads from it instead of re-typing
# threshold literals.
#
#   * Common interim efficacy threshold 0.999 for K >= 2 (no interim at K = 1).
#   * Final efficacy threshold by K (band-verified maximum-power frontier):
#       0.976 (K = 1, 2), 0.977 (K = 3..6), 0.978 (K = 7..10).
#   * Common futility threshold q = 0.90 (no futility rule at K = 1, q = NA).
#   * Headline designs: K = 3 and K = 5 are 0.999 / 0.977, q = 0.90.
BAYESGSD_DESIGNS <- local({
  interim_eff <- 0.999          # common interim efficacy threshold for K >= 2
  q_futility  <- 0.90           # common futility threshold (NA at K = 1)
  # Final efficacy threshold by number of looks K (1..10).
  final_by_K  <- c(`1` = 0.976, `2` = 0.976, `3` = 0.977, `4` = 0.977,
                   `5` = 0.977, `6` = 0.977, `7` = 0.978, `8` = 0.978,
                   `9` = 0.978, `10` = 0.978)
  list(interim_eff = interim_eff, q_futility = q_futility, final_by_K = final_by_K)
})

# Per-stage efficacy threshold vector p_{1:K} for the maximum-power design at K
# looks: c(rep(interim_eff, K - 1), final_by_K[K]); at K = 1 it is just the
# final threshold (the fixed-sample design has no interim).
bayesgsd_p_stage <- function(K) {
  fin <- BAYESGSD_DESIGNS$final_by_K[[as.character(K)]]
  if (is.null(fin))
    stop(sprintf("No calibrated final threshold defined for K = %d", K))
  int <- if (K == 1L) numeric(0) else rep(BAYESGSD_DESIGNS$interim_eff, K - 1L)
  c(int, fin)
}

# Common futility threshold q for the calibrated design at K looks (NA at K = 1,
# which has no futility rule).
bayesgsd_q <- function(K) if (K == 1L) NA_real_ else BAYESGSD_DESIGNS$q_futility

# --- sessionInfo / version inventory dump ----------------------------------
# Each driver calls bayesgsd_save_session("ScriptName") at the end of its
# successful run; the dump is written under OUTPUT_DIR with the script name
# as a tag, so it's easy to correlate downstream caches with versions.
bayesgsd_save_session <- function(tag = "session") {
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  si <- sessionInfo()
  out <- list(
    tag           = tag,
    timestamp     = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    R_version     = R.version.string,
    platform      = R.version$platform,
    project_root  = PROJECT_ROOT,
    sessionInfo   = si,
    package_versions = setNames(
      vapply(installed.packages()[, "Package"], function(p) {
        tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
      }, character(1L)),
      installed.packages()[, "Package"]
    )
  )
  rds_path <- file.path(OUTPUT_DIR, sprintf("%s_sessionInfo.rds", tag))
  saveRDS(out, file = rds_path)
  txt_path <- file.path(OUTPUT_DIR, sprintf("%s_sessionInfo.txt", tag))
  con <- file(txt_path, open = "w")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    sprintf("# %s  (run %s)", tag, out$timestamp),
    sprintf("# %s, %s", out$R_version, out$platform),
    sprintf("# PROJECT_ROOT = %s", PROJECT_ROOT),
    "",
    capture.output(print(si))
  ), con)
  invisible(rds_path)
}

cat(sprintf("BayesGSD setup: PROJECT_ROOT = %s, N_CORES = %d\n",
            PROJECT_ROOT, N_CORES))
