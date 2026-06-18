#' @title Gauss-Legendre quadrature convergence study (Appendix A.2)
#' @description Sweeps the number of Gauss-Legendre quadrature nodes Q over
#'   {64, 128, 256, 512} for both the Beta(1,1) baseline and the
#'   §3.3.2 logit-normal beta-mixture prior. For each Q the script:
#'
#'     (1) Recomputes posterior tail probabilities P(Delta < tau | data) at
#'         the §3.4 efficacy threshold tau = 0 for a fixed set of (e, m, j, n)
#'         interim-look cells covering the support of the §3.4 simulation.
#'         The maximum absolute difference against the Q=512 reference is
#'         the quadrature error at each Q.
#'
#'     (2) Drives the same R = 50,000 design-level operating characteristics
#'         pipeline used by §3.3.1 at p = 0.99, q = 0.90, for K in {3, 5, 9},
#'         under both H0 and H1. Reports type I and type II error rate and
#'         expected sample size as a function of Q.
#'
#'     (3) Compares (1) and (2) against the Monte Carlo SE implied by the
#'         R-trial budget. As long as the Q-induced bias is at least an
#'         order of magnitude below the MC SE, the fixed-Q quadrature is
#'         adequate for the design-stage calibration.
#'
#'   The output (`ADRENAL_QuadratureConvergence.rda`) and the figure
#'   (`ADRENAL_QuadratureConvergence.jpeg`) support Appendix A.2 of the
#'   main paper.
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
  library(RBesT)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
})

# --------------------------------------------------------------------------
# Helpers: source the kernel with a given Q, run a small experiment, return.
# Sourcing the kernel re-initialises the Gauss-Legendre nodes/weights from
# the option, so swapping Q is just a matter of setting the option and
# re-sourcing.
# --------------------------------------------------------------------------

source_kernel_with_Q <- function(Q) {
  options(BayesGSD.quadrature.n = as.integer(Q))
  # Force a clean reload so .init_adabay() rebuilds nodes/weights.
  # Use the same explicit relative path every other ADRENAL_* driver uses,
  # rather than relying on CODE_DIR surviving the local({...}) bootstrap.
  source("./Code/Code v1.0/adabay_bern.R")
}

# Source the kernel once at the script's default Q so the public API
# (initialisePriorSettings, runTrialMonitoring, getOperatingCharacteristics,
# .posterior_tail_prob) is in scope before we build the prior mixture below.
# The per-Q sweep further down re-sources with each Q in Q_GRID.
source_kernel_with_Q(128L)

# --------------------------------------------------------------------------
# Priors used in (1) and (2).
# --------------------------------------------------------------------------

prior_beta11 <- list(
  list(weight = 1, alpha = 1, beta = 1),
  list(weight = 1, alpha = 1, beta = 1)
)

# Logit-normal mixture (§3.3.2). Built once with seed=21 to match the §3.3.2
# benchmark; we just store the mixture here to avoid bringing in RBesT's
# automixfit timing into the quadrature study.
set.seed(21L)
.logit <- function(p) log(p / (1 - p))
PRIOR_MU    <- .logit(0.33)
PRIOR_SIGMA <- 0.5
samples_c <- 1 / (1 + exp(-rnorm(50000L, PRIOR_MU, PRIOR_SIGMA)))
samples_t <- 1 / (1 + exp(-rnorm(50000L, PRIOR_MU, PRIOR_SIGMA)))
prior_logitnormal <- initialisePriorSettings(
  priorSamples       = list(samples_c, samples_t),
  numberOfComponents = 4L
)
cat("Logit-normal beta-mixture fit:\n")
print(prior_logitnormal[[1]])
print(prior_logitnormal[[2]])

# --------------------------------------------------------------------------
# (1) Per-trial posterior tail probability accuracy at fixed (e, m, j, n)
# cells. We choose cells spanning early-, mid-, and late-look event counts at
# the §3.4 design skeleton (N = 3,800, balanced allocation).
# --------------------------------------------------------------------------

# Mimic the §3.4 looks for K=9 (most stringent test of fixed-Q quadrature
# because the posterior is sharpest at the final stage). N_MAX and make_looks()
# are the single source of truth in adabay_bern_setup.R (sourced via bootstrap).
LOOKS_K9  <- make_looks(9L)  # 423, 845, ..., 3800
# Per-arm look sizes under 1:1 allocation.
NPC_K9    <- LOOKS_K9 - round(LOOKS_K9 / 2)
NPT_K9    <- round(LOOKS_K9 / 2)

# A small representative event-count grid per stage.
event_grid_per_stage <- function(m, n) {
  # Around the expected counts at theta = 0.33 (control) and 0.28 (treatment H1),
  # plus a span +/- 2 SD for binomial.
  exp_c <- round(m * 0.33); sd_c <- sqrt(m * 0.33 * 0.67)
  exp_t <- round(n * 0.28); sd_t <- sqrt(n * 0.28 * 0.72)
  ec <- as.integer(unique(round(c(exp_c - 2*sd_c, exp_c - sd_c, exp_c,
                                  exp_c + sd_c, exp_c + 2*sd_c))))
  et <- as.integer(unique(round(c(exp_t - 2*sd_t, exp_t - sd_t, exp_t,
                                  exp_t + sd_t, exp_t + 2*sd_t))))
  ec <- pmin(pmax(ec, 0L), as.integer(m))
  et <- pmin(pmax(et, 0L), as.integer(n))
  expand.grid(ec = ec, et = et, stringsAsFactors = FALSE)
}

# Sample 3 stages: stage 3 (early), stage 6 (mid), stage 9 (final).
STAGE_IDX <- c(3L, 6L, 9L)

Q_GRID    <- c(64L, 128L, 256L, 512L)
threshold <- 0    # the §3.4 efficacy threshold (P(Delta < 0))

# Storage: for each prior and each Q, the per-cell posterior tail prob.
quad_cells <- list()
for (prior_label in c("Beta(1,1)", "LogitNormal mixture")) {
  prior <- if (prior_label == "Beta(1,1)") prior_beta11 else prior_logitnormal
  for (Q in Q_GRID) {
    source_kernel_with_Q(Q)
    rows <- list()
    for (s in STAGE_IDX) {
      m  <- NPC_K9[s]; n <- NPT_K9[s]
      eg <- event_grid_per_stage(m, n)
      ec <- eg$ec; et <- eg$et
      tp <- .posterior_tail_prob(
        ec = ec, et = et, mc = m, nt = n,
        priors_c = prior[[1]], priors_t = prior[[2]],
        thresholds = threshold, type = "absolute.risk",
        alternative = "less",
        chunk_size = .Machine$integer.max, cores = 1L
      )
      rows[[length(rows) + 1L]] <- data.frame(
        prior = prior_label, Q = Q, stage = s, m = m, n = n,
        ec = ec, et = et, tail = as.numeric(tp)
      )
    }
    quad_cells[[paste(prior_label, Q, sep = "/")]] <- do.call(rbind, rows)
  }
}

# Q=512 is the reference. Tabulate max |dtail| vs Q for each prior.
ref_label <- "512"
diff_tbl <- do.call(rbind, lapply(c("Beta(1,1)", "LogitNormal mixture"), function(pl) {
  ref <- quad_cells[[paste(pl, ref_label, sep = "/")]]$tail
  out <- data.frame(prior = pl,
                    Q = Q_GRID,
                    max_abs_diff = NA_real_,
                    rms_diff     = NA_real_)
  for (i in seq_along(Q_GRID)) {
    Q <- Q_GRID[i]
    if (Q == 512L) { out$max_abs_diff[i] <- 0; out$rms_diff[i] <- 0; next }
    d <- quad_cells[[paste(pl, Q, sep = "/")]]$tail - ref
    out$max_abs_diff[i] <- max(abs(d))
    out$rms_diff[i]     <- sqrt(mean(d^2))
  }
  out
}))

cat("\n=== (1) Per-trial posterior tail prob accuracy vs Q=512 reference ===\n")
print(diff_tbl, row.names = FALSE)

# --------------------------------------------------------------------------
# (2) Design-level OC at R=50,000 under both priors, K in {3, 5, 9},
# H0 and H1. We use the §3.3 protocol (p = 0.99, q = 0.90, K = 1 fixed-sample
# baseline omitted for brevity). The proposed framework's MC SE at R = 50k
# is ~ sqrt(0.025 * 0.975 / 50000) ~ 0.07 pp for the type I error rate.
# --------------------------------------------------------------------------

R_QC          <- 50000L
P_EFF         <- 0.99
P_FUT_LOW     <- 1 - 0.90
P_CONTROL     <- 0.33
DELTA_ALT     <- -0.05
SEED          <- 21L
K_QC          <- c(3L, 5L, 9L)

# Equally-spaced looks at information fractions j/K for j = 1, ..., K, matching
# the §3.3 protocol; make_looks() is the single source of truth (setup helper).

run_oc_one <- function(prior, K, hypothesis, R, seed) {
  looks <- make_looks(K)
  rates <- if (hypothesis == "H0") c(P_CONTROL, P_CONTROL)
           else c(P_CONTROL, P_CONTROL + DELTA_ALT)
  sim <- runTrialMonitoring(
    simulationSettings = list(numberOfTrials = R, numberOfSubjects = looks,
                              allocationRatio = 1, rates = rates,
                              hypothesis = if (hypothesis == "H0") "null" else "alternative",
                              seed = seed),
    posteriorSettings  = list(effectThresholds =
                                list(efficacy = list(size = 0, type = "absolute.risk"),
                                     equivalence = NULL,
                                     futility = list(size = 0, type = "absolute.risk")),
                              alternative = "less"),
    priorSettings      = prior
  )
  fut_thr <- matrix(P_FUT_LOW, nrow = K, ncol = 1)
  fut_thr[K, 1L] <- 0
  td <- list(
    numberOfSubjects = looks,
    effectThresholds = list(
      efficacy = list(size = 0, type = "absolute.risk"),
      equivalence = NULL,
      futility = list(size = 0, type = "absolute.risk")
    ),
    probabilityThresholds = list(
      efficacy = matrix(P_EFF, K, 1L),
      equivalence = NULL,
      futility = fut_thr
    )
  )
  oc <- getOperatingCharacteristics(sim, td)
  list(
    eff_rate = if (hypothesis == "H0") oc$type1ErrorRate else 1 - oc$type2ErrorRate,
    EN       = oc$expSampleSize
  )
}

design_tbl <- list()
for (prior_label in c("Beta(1,1)", "LogitNormal mixture")) {
  prior <- if (prior_label == "Beta(1,1)") prior_beta11 else prior_logitnormal
  for (Q in Q_GRID) {
    source_kernel_with_Q(Q)
    for (K in K_QC) for (h in c("H0", "H1")) {
      cat(sprintf("  [%s, Q=%d, K=%d, %s] ... ", prior_label, Q, K, h))
      r <- run_oc_one(prior, K, h, R_QC, SEED)
      cat(sprintf("eff_rate=%.4f, E[N]=%.0f\n", r$eff_rate, r$EN))
      design_tbl[[length(design_tbl) + 1L]] <- data.frame(
        prior = prior_label, Q = Q, K = K, hypothesis = h,
        eff_rate = r$eff_rate, EN = r$EN
      )
    }
  }
}
design_tbl <- do.call(rbind, design_tbl)

# Difference vs Q=512 reference for each (prior, K, hypothesis).
design_diff <- do.call(rbind, by(design_tbl,
  list(design_tbl$prior, design_tbl$K, design_tbl$hypothesis), function(df) {
    ref_eff <- df$eff_rate[df$Q == 512L]
    ref_EN  <- df$EN[df$Q == 512L]
    mc_se   <- sqrt(ref_eff * (1 - ref_eff) / R_QC)
    df$abs_diff_eff <- abs(df$eff_rate - ref_eff)
    df$abs_diff_EN  <- abs(df$EN - ref_EN)
    df$mc_se_eff    <- mc_se
    df
  }))
rownames(design_diff) <- NULL

cat("\n=== (2) Design-level OC differences vs Q=512 reference (R = 50,000) ===\n")
print(design_diff, row.names = FALSE, digits = 4)

# --------------------------------------------------------------------------
# (3) Persist + plot
# --------------------------------------------------------------------------

OUT_PATH_RDA  <- file.path(OUTPUT_DIR, "ADRENAL_QuadratureConvergence.rda")
OUT_PATH_JPEG <- file.path(OUTPUT_DIR, "ADRENAL_QuadratureConvergence.jpeg")
save(quad_cells, diff_tbl, design_diff, file = OUT_PATH_RDA)
cat(sprintf("Wrote %s\n", OUT_PATH_RDA))

# Left panel: per-trial tail prob max |diff| vs Q (log scale on y).
# Style follows the other manuscript figures: theme_bw, ggsci BMJ palette,
# bottom legend, patchwork plot_annotation(tag_levels = "A") for the A/B
# subfigure tags.
p_left <- ggplot(diff_tbl[diff_tbl$Q < 512L, ],
                 aes(x = factor(Q), y = max_abs_diff, fill = prior)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.5) +
  scale_y_log10() +
  scale_fill_bmj() +
  xlab("Quadrature nodes Q") +
  ylab(expression("max |" * P[Q]("Δ<0" * "|" * data) - P[512]("Δ<0" * "|" * data) * "|")) +
  ggtitle("Per-trial posterior tail prob accuracy") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank())

# Right panel: design-level abs diff vs MC SE bound.
design_diff_q <- design_diff[design_diff$Q < 512L, ]
design_diff_q$panel <- sprintf("K=%d, %s", design_diff_q$K, design_diff_q$hypothesis)
mc_se_bound <- max(design_diff$mc_se_eff)

p_right <- ggplot(design_diff_q,
                  aes(x = factor(Q), y = abs_diff_eff, colour = prior,
                      shape = panel, group = interaction(prior, panel))) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = mc_se_bound, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 0.6, y = mc_se_bound, vjust = -0.4,
           label = sprintf("MC SE at R=%s", format(R_QC, big.mark = ",")),
           colour = "grey40", hjust = 0, size = 3.5) +
  scale_color_bmj() +
  xlab("Quadrature nodes Q") +
  ylab("|design type-I/II rate at Q - rate at Q=512|") +
  ggtitle("Design-level OC accuracy vs MC SE") +
  labs(colour = NULL, shape = "Design cell") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "vertical") +
  guides(colour = guide_legend(order = 1),
         shape  = guide_legend(order = 2, ncol = 5))

combined <- (p_left + p_right) +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(tag_levels = "A")
ggsave(OUT_PATH_JPEG, combined, width = 12, height = 6, dpi = 200)
cat(sprintf("Wrote %s\n", OUT_PATH_JPEG))


# Persist sessionInfo and package versions alongside the cache.
bayesgsd_save_session("ADRENAL_Diagnostics_QuadratureConvergence")
