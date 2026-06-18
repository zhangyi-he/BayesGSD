#' @title Per-look stopping probabilities for the K=5 ADRENAL re-design
#'   (Supplement Table S2 / `tab:S2`)
#' @description Derives Supplement Table S2 directly from the cached posterior
#'   probabilities in `ADRENAL_OperatingCharacteristics_H[01].rda`. The table
#'   gives per-look efficacy- and futility-stopping probabilities and the
#'   cumulative expected sample size for the binary-endpoint re-design with
#'   K = 5 total analyses, single-criterion efficacy threshold p = 0.99, and
#'   binding futility threshold q = 0.90.
#'
#'   The cached caches are simulated at the union of looks across K = 2..10
#'   on R = 10^6 virtual trials per hypothesis. This script picks the five
#'   look times that correspond to the K = 5 design (cumulative both-arm
#'   N = 760, 1520, 2280, 3040, 3800) and applies the manuscript's
#'   binding-futility convention from Section 2.1 (decision criteria): efficacy is evaluated
#'   at every look, futility only at interim looks k = 1, ..., K - 1.
#'
#'   Output: prints the LaTeX-ready Table S2 rows and a per-cell comparison
#'   against the manuscript-supplement values; writes the structured result
#'   to `ADRENAL_Evaluation_StoppingProbabilities.rda` for downstream consumers, alongside
#'   `_sessionInfo.{rds,txt}` via `bayesgsd_save_session()`.
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

# ----------------------------------------------------------------------------
# Inputs
# ----------------------------------------------------------------------------

OC_H0_PATH <- file.path(OUTPUT_DIR, "ADRENAL_OperatingCharacteristics_H0.rda")
OC_H1_PATH <- file.path(OUTPUT_DIR, "ADRENAL_OperatingCharacteristics_H1.rda")
RDA_PATH   <- file.path(OUTPUT_DIR, "ADRENAL_Evaluation_StoppingProbabilities.rda")

missing <- c(if (!file.exists(OC_H0_PATH)) OC_H0_PATH,
             if (!file.exists(OC_H1_PATH)) OC_H1_PATH)
if (length(missing) > 0L) {
  stop(sprintf(
    "Required OC cache(s) missing:\n  %s\nRun ADRENAL_Evaluation_OperatingCharacteristics.R first.",
    paste(missing, collapse = "\n  ")
  ))
}

# K = 5 design specification.
K          <- 5L
N_TOT      <- c(760L, 1520L, 2280L, 3040L, 3800L)   # both-arm cumulative
PER_ARM    <- N_TOT / 2L                            # n_{c,k}
P_EFF      <- 0.99                                  # single-criterion efficacy
P_FUT_LOW  <- 1 - 0.90                              # PostPr(Delta < 0) < 0.10

# ----------------------------------------------------------------------------
# Compute per-look stopping probabilities for one cached cache
# ----------------------------------------------------------------------------

compute_table <- function(rda_path) {
  e <- new.env(parent = emptyenv())
  load(rda_path, envir = e)
  pp <- e$trialSimulation$virtualTrials$posteriorProbabilities
  union_looks <- e$simulationSettings$numberOfSubjects
  R <- e$trialSimulation$virtualTrials$numberOfTrials

  # Locate the K=5 looks within the simulated union.
  stage_idx <- vapply(N_TOT, function(n) {
    i <- which(union_looks == n)
    if (length(i) != 1L)
      stop(sprintf("Look N=%d not present (exactly once) in union_looks", n))
    as.integer(i)
  }, integer(1))

  # Section 2.1 decision rule: futility evaluated only at interim looks k = 1, ..., K - 1.
  eff_hit <- matrix(FALSE, nrow = R, ncol = K)
  fut_hit <- matrix(FALSE, nrow = R, ncol = K)
  for (k in seq_len(K)) {
    eff_hit[, k] <- pp[[stage_idx[k]]]$efficacy["0", ] > P_EFF
    if (k < K) {
      fut_hit[, k] <- pp[[stage_idx[k]]]$futility["0", ] < P_FUT_LOW
    }
  }

  # First-stop semantics: any trial not stopped at an interim is forced to
  # stop at the final look. There, the only outcomes are "efficacy at k=K"
  # (eff_hit[K] = TRUE) or "no efficacy" (counted as the final-stage beta cell
  # of Table S2). The ties.method = "first" rule mirrors fast_oc().
  combined <- eff_hit | fut_hit
  combined[, K] <- TRUE
  stage_stopped <- max.col(combined, ties.method = "first")
  is_eff <- eff_hit[cbind(seq_len(R), stage_stopped)]

  alpha_k <- numeric(K)
  beta_k  <- numeric(K)
  for (k in seq_len(K)) {
    alpha_k[k] <- mean(stage_stopped == k &  is_eff) * 100
    beta_k[k]  <- mean(stage_stopped == k & !is_eff) * 100
  }
  E_N <- mean(N_TOT[stage_stopped])

  list(R = R, K = K, per_arm = PER_ARM, N_tot = N_TOT,
       alpha_k = alpha_k, beta_k = beta_k,
       alpha_total = sum(alpha_k), beta_total = sum(beta_k),
       E_N = E_N)
}

cat(sprintf("Computing K=%d per-look stopping probabilities from cached posteriors\n", K))
cat(sprintf("Design: p = %.2f efficacy, q = %.2f binding futility\n", P_EFF, 1 - P_FUT_LOW))
cat("Cache inputs:\n")
cat(sprintf("  H0: %s\n", OC_H0_PATH))
cat(sprintf("  H1: %s\n\n", OC_H1_PATH))

tab_H0 <- compute_table(OC_H0_PATH)
tab_H1 <- compute_table(OC_H1_PATH)

# ----------------------------------------------------------------------------
# Print human-readable summary
# ----------------------------------------------------------------------------

print_summary <- function(tag, tab) {
  cat(sprintf("=== %s ===\n", tag))
  cat("   k  n_c,k  N_tot  alpha_k(%)  beta_k(%)\n")
  for (k in seq_len(tab$K)) {
    cat(sprintf("   %d  %5d  %5d   %8.4f   %8.4f\n",
                k, tab$per_arm[k], tab$N_tot[k],
                tab$alpha_k[k], tab$beta_k[k]))
  }
  cat(sprintf("   Total: alpha = %.4f%%, beta = %.4f%%, E[N] = %.2f\n\n",
              tab$alpha_total, tab$beta_total, tab$E_N))
}

print_summary("H0 (rates = 0.33, 0.33; Delta = 0)", tab_H0)
print_summary("H1 (rates = 0.33, 0.28; Delta = -0.05)", tab_H1)

# ----------------------------------------------------------------------------
# Emit LaTeX-ready Table S2 rows (display precision: 2 decimal places)
# ----------------------------------------------------------------------------

fmt_num <- function(x) {
  # Use comma thousands separator a la \texttt{1{,}520}.
  s <- formatC(x, big.mark = ",", format = "d")
  gsub(",", "{,}", s, fixed = TRUE)
}
fmt_pct <- function(x) sprintf("%.2f\\%%", round(x, 2))

cat("=== LaTeX-ready Table S2 body (paste between \\midrule and \\bottomrule) ===\n")
for (k in seq_len(K)) {
  cat(sprintf("%d & %s & %s & %s & %s & %s & %s \\\\\n",
              k, fmt_num(PER_ARM[k]),
              fmt_pct(tab_H0$alpha_k[k]), fmt_pct(tab_H0$beta_k[k]),
              fmt_pct(tab_H1$alpha_k[k]), fmt_pct(tab_H1$beta_k[k]),
              fmt_num(N_TOT[k])))
}
cat("\\midrule\n")
cat(sprintf("Total & --- & %s & %s & %s & %s & %s ($H_1$) \\\\\n",
            fmt_pct(tab_H0$alpha_total), fmt_pct(tab_H0$beta_total),
            fmt_pct(tab_H1$alpha_total), fmt_pct(tab_H1$beta_total),
            fmt_num(round(tab_H1$E_N))))
cat("\n")

# ----------------------------------------------------------------------------
# Persist for downstream consumers + sessionInfo
# ----------------------------------------------------------------------------

K_DESIGN <- K
save(tab_H0, tab_H1, K_DESIGN, N_TOT, PER_ARM, P_EFF, P_FUT_LOW,
     file = RDA_PATH)
cat(sprintf("Wrote %s\n", RDA_PATH))

bayesgsd_save_session("ADRENAL_Evaluation_StoppingProbabilities")
