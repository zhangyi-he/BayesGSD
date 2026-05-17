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
#'     CODE_DIR, N_CORES, bayesgsd_save_session().
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
