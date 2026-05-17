#' @title Splice the four benchmark tables into the manuscript .tex file
#' @description Renders the bench_log caches into LaTeX tabular code and
#'   replaces the block between matching `% BEGIN: <key>` and `% END: <key>`
#'   markers inside ZH2023_Manuscript.tex. Four tables are produced per
#'   prior, one for the type I error rate (H0 cells) and one for the power
#'   (H1 cells). Each table contains, per K, one row each for BATSS at
#'   R = 5,000, adaptr at R = 5,000, Proposed at R = 5,000 and Proposed at
#'   R = 1,000,000.
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


OUTPUT_DIR     <- "./Output/Output v1.0"
MANUSCRIPT_DIR <- "./Article/Manuscript/Manuscript v1.0"
TEX_PATH       <- file.path(MANUSCRIPT_DIR, "ZH2023_Manuscript.tex")

# --------------------------------------------------------------------------
# Formatting helpers
# --------------------------------------------------------------------------

mc_se <- function(p, R) sqrt(p * (1 - p) / R)
fmt_pct <- function(x, d = 2) {
  if (is.na(x)) "\\texttt{NA}" else paste0(formatC(100 * x, digits = d, format = "f"), "\\,\\%")
}
fmt_pp <- function(x, d = 3) {
  # MC SE column reported to 3 decimal places (pp) so that the matched-budget
  # comparisons in the body text remain interpretable when the SE is near the
  # rounding limit of the rate estimate itself.
  if (is.na(x)) "---" else paste0(formatC(100 * x, digits = d, format = "f"), "\\,pp")
}
fmt_int <- function(x) {
  if (is.na(x)) "\\texttt{NA}" else format(round(x), big.mark = ",")
}
fmt_sec <- function(x) {
  if (is.na(x)) "\\texttt{NA}"
  else if (x >= 100) format(round(x), big.mark = ",")
  else formatC(x, digits = 1, format = "f")
}
# Order in which methods appear inside each K block. adaptr is benchmarked
# at the matched-budget head-to-head R=5,000 alongside BATSS and the
# proposed framework; the proposed framework is also reported at R=1,000,000.
METHOD_ORDER <- list(
  list(slot = "batss",          label = "\\texttt{BATSS}", R_label = "$5{,}000$"),
  list(slot = "adaptr_small",   label = "\\texttt{adaptr}", R_label = "$5{,}000$"),
  list(slot = "proposed_small", label = "Proposed",         R_label = "$5{,}000$"),
  list(slot = "proposed_big",   label = "Proposed",         R_label = "$1{,}000{,}000$")
)

build_rows_filtered <- function(bench_log, hypothesis_keep) {
  keys <- names(bench_log)[order(sapply(bench_log, function(r) r$design$K),
                                  sapply(bench_log, function(r) r$hypothesis))]
  keys <- Filter(function(k) bench_log[[k]]$hypothesis == hypothesis_keep, keys)
  out <- character(0)
  for (key in keys) {
    r       <- bench_log[[key]]
    K_str   <- sprintf("$%d$", r$design$K)
    for (m in METHOD_ORDER) {
      cell <- r[[m$slot]]
      if (is.null(cell)) next
      out <- c(out, sprintf(
        "    %s & %s & %s & %s & %s & %s & %s \\\\",
        K_str, m$label, m$R_label,
        fmt_pct(cell$efficacy_prob),
        fmt_pp(mc_se(cell$efficacy_prob, cell$R)),
        fmt_int(cell$expected_n),
        fmt_sec(cell$elapsed)
      ))
    }
    out <- c(out, "    \\midrule")
  }
  if (length(out) && grepl("midrule", out[length(out)])) out <- out[-length(out)]
  out
}

build_table_block <- function(rda_path, marker, label, caption,
                              hypothesis, metric_header) {
  e <- new.env(); load(rda_path, envir = e)
  bench_log <- e$bench_log
  c(
    sprintf("%% BEGIN: %s", marker),
    "\\begin{table}[!ht]",
    "  \\centering",
    "  \\begin{tabular}{rllrrrr}",
    "    \\toprule",
    paste0("    $K$ & Method & $R$ & ", metric_header,
           " & MC SE & $\\mathbb{E}(N)$ & Time (s) \\\\"),
    "    \\midrule",
    build_rows_filtered(bench_log, hypothesis),
    "    \\bottomrule",
    "  \\end{tabular}",
    paste0("  \\caption{", caption, "}"),
    paste0("  \\label{", label, "}"),
    "\\end{table}",
    sprintf("%% END: %s", marker)
  )
}

# Splice a block of lines into the manuscript between
#   `% BEGIN: <marker>` and `% END: <marker>`. The marker lines themselves
# are emitted by `build_table_block`, so the splice replaces *including*
# those marker lines.
splice_block <- function(tex_lines, marker, new_block) {
  begin_re <- sprintf("^%% BEGIN: %s\\s*$", marker)
  end_re   <- sprintf("^%% END: %s\\s*$",   marker)
  i1 <- grep(begin_re, tex_lines)
  i2 <- grep(end_re,   tex_lines)
  if (length(i1) != 1L || length(i2) != 1L) {
    stop("Marker ", marker, " not found uniquely in manuscript (begin: ",
         length(i1), ", end: ", length(i2), ").")
  }
  c(tex_lines[seq_len(i1 - 1L)], new_block, tex_lines[seq.int(i2 + 1L, length(tex_lines))])
}

# --------------------------------------------------------------------------
# Captions
# --------------------------------------------------------------------------

main_common <- paste(
  "Estimates and wall-clock cost for the proposed semi-simulation framework, \\texttt{BATSS} and \\texttt{adaptr} on the fixed-sample design ($K=1$) and group sequential designs with three, five, seven and nine equally spaced looks ($K\\in\\{3,5,7,9\\}$) for the ADRENAL re-design ($\\vartheta=0.33$, $\\delta=-0.05$ under $H_{1}$, efficacy threshold $p=0.99$, futility threshold $q=0.90$, total sample size $N=3{,}800$).",
  "``MC SE'' denotes the Monte Carlo standard error of the corresponding estimate, in percentage points (pp).",
  "$\\mathbb{E}(N)$ is the estimated expected total sample size.",
  "All three methods use independent $\\operatorname{Beta}(1,1)$ priors on the per-arm rates, except that \\texttt{BATSS} uses \\texttt{INLA} with weakly informative Gaussian priors on the logit-scale coefficients."
)

supp_common <- paste(
  "As Table~\\ref{tab:331}, but replacing the independent $\\operatorname{Beta}(1,1)$ priors with an informative independent logit-normal prior on each arm, $\\operatorname{logit}(\\vartheta)\\sim\\mathcal{N}(\\mu,\\sigma^{2})$ with $\\mu=\\operatorname{logit}(0.33)$ and $\\sigma=0.5$.",
  "The proposed framework approximates this prior by a finite beta mixture via the prior-approximation step of Section~\\ref{sec:23}.",
  "\\texttt{adaptr} is supplied with the logit-normal prior directly through inverse-CDF sampling on a fine grid of $\\operatorname{logit}(\\vartheta)$, with no beta-mixture intermediate.",
  "\\texttt{BATSS} is supplied with the closest \\texttt{INLA}-representable approximation, namely independent normal priors on the regression coefficients of the model $y\\sim$\\,group with marginal-matched moments $\\beta_{0}\\sim\\mathcal{N}(\\mu,\\sigma^{2})$ and $\\beta_{1}\\sim\\mathcal{N}(0,2\\sigma^{2})$. The cross-covariance $\\operatorname{Cov}(\\beta_{0},\\beta_{1})=-\\sigma^{2}$ implied by the independent logit-normal on the arms cannot be expressed in \\texttt{INLA}'s fixed-effects framework."
)

# --------------------------------------------------------------------------
# Build the four table blocks
# --------------------------------------------------------------------------

block_type1 <- build_table_block(
  rda_path     = file.path(OUTPUT_DIR, "ADRENAL_Benchmark_UniformPrior.rda"),
  marker        = "benchmark_table_type1",
  label         = "tab:331",
  caption       = paste("Type I error rate estimates.", main_common),
  hypothesis    = "H0",
  metric_header = "Type I error"
)
block_power <- build_table_block(
  rda_path     = file.path(OUTPUT_DIR, "ADRENAL_Benchmark_UniformPrior.rda"),
  marker        = "benchmark_table_power",
  label         = "tab:332",
  caption       = "Power estimates. As Table~\\ref{tab:331}, but reporting the estimated power rather than the type I error rate.",
  hypothesis    = "H1",
  metric_header = "Power"
)
block_lt1 <- build_table_block(
  rda_path     = file.path(OUTPUT_DIR, "ADRENAL_Benchmark_LogitNormalPrior.rda"),
  marker        = "benchmark_table_logitnormal_type1",
  label         = "tab:333",
  caption       = paste("Type I error rate estimates under the logit-normal prior.", supp_common),
  hypothesis    = "H0",
  metric_header = "Type I error"
)
block_lp <- build_table_block(
  rda_path     = file.path(OUTPUT_DIR, "ADRENAL_Benchmark_LogitNormalPrior.rda"),
  marker        = "benchmark_table_logitnormal_power",
  label         = "tab:334",
  caption       = "Power estimates under the logit-normal prior. As Table~\\ref{tab:333}, but reporting the estimated power rather than the type I error rate.",
  hypothesis    = "H1",
  metric_header = "Power"
)

# --------------------------------------------------------------------------
# Splice the four blocks into the manuscript .tex file
# --------------------------------------------------------------------------

tex_lines <- readLines(TEX_PATH, warn = FALSE)
tex_lines <- splice_block(tex_lines, "benchmark_table_type1",              block_type1)
tex_lines <- splice_block(tex_lines, "benchmark_table_power",              block_power)
tex_lines <- splice_block(tex_lines, "benchmark_table_logitnormal_type1",  block_lt1)
tex_lines <- splice_block(tex_lines, "benchmark_table_logitnormal_power",  block_lp)
writeLines(tex_lines, TEX_PATH)
cat("Spliced 4 tables into", TEX_PATH, "\n")

# Persist sessionInfo and package versions alongside the cache.
bayesgsd_save_session("ADRENAL_Benchmark_GenerateTables")
