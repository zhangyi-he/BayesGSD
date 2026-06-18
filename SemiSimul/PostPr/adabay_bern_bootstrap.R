#' @title BayesGSD driver bootstrap
#' @description Resolves the project root and sources `adabay_bern_setup.R`.
#'   Every `ADRENAL_*.R` driver `source()`s this file early instead of
#'   duplicating the resolution + setup-source logic in each script.
#'
#'   Resolution order: (1) `BAYESGSD_ROOT` env var if it points to an existing
#'   directory; (2) auto-detect by walking up from `getwd()` looking for the
#'   sibling `Code/` and `Article/` directories that mark the repository root;
#'   (3) stop loudly with an actionable message.
#'
#'   On success, the corresponding `adabay_bern_setup.R` is sourced into
#'   the global environment, exporting `PROJECT_ROOT`, `OUTPUT_DIR`,
#'   `MANUSCRIPT_DIR`, `SUPPLEMENT_DIR`, `CODE_DIR`, `N_CORES` and
#'   `bayesgsd_save_session()` for driver use.
#' @author Zhangyi He, Feng Yu, Suzie Cro, Laurent Billot

.adabay_bootstrap <- function() {
  bootstrap_root <- Sys.getenv("BAYESGSD_ROOT", unset = "")
  if (!nzchar(bootstrap_root) || !dir.exists(bootstrap_root)) {
    cur <- getwd()
    while (cur != dirname(cur)) {
      if (dir.exists(file.path(cur, "Code")) &&
          dir.exists(file.path(cur, "Article"))) {
        bootstrap_root <- cur
        break
      }
      cur <- dirname(cur)
    }
  }
  if (!nzchar(bootstrap_root) || !dir.exists(bootstrap_root)) {
    stop("Could not resolve BAYESGSD_ROOT. Set it before running, e.g.:\n",
         "  R> Sys.setenv(BAYESGSD_ROOT = \"/absolute/path/to/repo\")\n",
         "  $ export BAYESGSD_ROOT=/absolute/path/to/repo")
  }
  source(file.path(bootstrap_root, "Code", "Code v1.0", "adabay_bern_setup.R"),
         local = FALSE)
}
.adabay_bootstrap()
# Tidy up the helper so it doesn't linger in globalenv after the bootstrap.
if (exists(".adabay_bootstrap", envir = globalenv(), inherits = FALSE)) {
  rm(".adabay_bootstrap", envir = globalenv())
}
