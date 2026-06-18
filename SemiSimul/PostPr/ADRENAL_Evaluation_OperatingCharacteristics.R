#' @title Operating-characteristic behaviour of Bayesian GSDs for the ADRENAL
#'   re-design across a (K, p, q) design grid
#' @description Runs 1,000,000 virtual trials at the union of look times across
#'   K = 2,...,10 total analyses under both H0 and H1 (one shared trial cache
#'   per hypothesis), then sweeps the cached posterior probabilities over the
#'   design grid of Section 3.4. For each combination of (binding or
#'   non-binding futility) and (single- or dual-criterion efficacy), reports
#'   the type I and type II error rates, the early-stopping probability and
#'   the expected sample size, and renders the 12-panel figures over the
#'   p in {0.90, 0.95, 0.99} and q in {0.85, 0.90, 0.95, 1.00} thresholds.
#'   Heavy lifting (the simulation) is cached on disk and resumed on
#'   subsequent runs.
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
  library(parallel)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
  library(scales)
})

source("./Code/Code v1.0/adabay_bern.R")

OUTPUT_DIR <- "./Output/Output v1.0"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# Trial parameters (Section 3.4)
# ----------------------------------------------------------------------------

P_CONTROL  <- 0.33
DELTA_ALT  <- -0.05
P_TREAT_H1 <- P_CONTROL + DELTA_ALT
# N_MAX and make_looks() are the single source of truth in adabay_bern_setup.R
# (sourced via the bootstrap above); do not redefine them here.
SEED       <- 21L
R_TRIALS   <- 1000000L

# Posterior-probability thresholds for efficacy and (binding) futility, on the
# framework's lower-tail convention: P_FUT_LOW = 1 - q so that q in {0.85,
# 0.90, 0.95, 1.00} corresponds to P_FUT_LOW in {0.15, 0.10, 0.05, 0.00};
# q = 1.00 means no futility rule (P_FUT_LOW = 0.00 never fires).
P_EFF_GRID <- c(0.90, 0.95, 0.99)
P_FUT_GRID <- c(0.15, 0.10, 0.05, 0.00)
Q_FUT_GRID <- 1 - P_FUT_GRID

# Total analyses K ranges from 1 (fixed sample) to 10 (9 interim + final).
K_GRID <- 1:10

# Union of look times across K = 2,...,10 (i.e. 1 to 9 interim analyses).
# make_looks() is provided by adabay_bern_setup.R (single source of truth).
union_looks <- sort(unique(unlist(lapply(2:10, make_looks))))
cat("Union of look times (n =", length(union_looks), "):",
    paste(union_looks, collapse = ","), "\n")

# ----------------------------------------------------------------------------
# Phase 1: simulate virtual trials under H0 and H1 (cache on disk)
# ----------------------------------------------------------------------------

simulate_trials <- function(hypothesis) {
  stopifnot(hypothesis %in% c("H0", "H1"))
  rates <- if (hypothesis == "H0") c(P_CONTROL, P_CONTROL) else c(P_CONTROL, P_TREAT_H1)
  simulationSettings <- list(
    numberOfTrials   = R_TRIALS,
    numberOfSubjects = union_looks,
    allocationRatio  = 1,
    rates            = rates,
    hypothesis       = if (hypothesis == "H0") "null" else "alternative",
    seed             = SEED
  )
  posteriorSettings <- list(
    effectThresholds = list(
      efficacy    = list(size = c(-0.05, 0.00), type = "absolute.risk"),
      equivalence = NULL,
      futility    = list(size = c(0.00),        type = "absolute.risk")
    ),
    alternative = "less"
  )
  priorSettings <- initialisePriorSettings()

  t0 <- proc.time()
  trialSimulation <- runTrialMonitoring(
    simulationSettings = simulationSettings,
    posteriorSettings  = posteriorSettings,
    priorSettings      = priorSettings
  )
  t_elapsed <- (proc.time() - t0)[["elapsed"]]
  list(
    simulationSettings = simulationSettings,
    posteriorSettings  = posteriorSettings,
    priorSettings      = priorSettings,
    trialSimulation    = trialSimulation,
    elapsed            = t_elapsed
  )
}

for (hyp in c("H0", "H1")) {
  fp <- file.path(OUTPUT_DIR, sprintf("ADRENAL_OperatingCharacteristics_%s.rda", hyp))
  cache_valid <- FALSE
  if (file.exists(fp)) {
    e <- new.env(); load(fp, envir = e)
    cached_R <- e$simulationSettings$numberOfTrials
    if (identical(as.integer(cached_R), as.integer(R_TRIALS))) {
      cat(sprintf("Cache exists with matching R = %d: %s (skipping simulation)\n",
                  cached_R, fp))
      cache_valid <- TRUE
    } else {
      cat(sprintf("Cache at %s has R = %d, expected R = %d -- re-simulating.\n",
                  fp, cached_R, R_TRIALS))
    }
    rm(e)
  }
  if (!cache_valid) {
    cat(sprintf("\n=== Simulating virtual trials under %s ===\n", hyp))
    res <- simulate_trials(hyp)
    cat(sprintf("  elapsed = %.1f s\n", res$elapsed))
    simulationSettings <- res$simulationSettings
    posteriorSettings  <- res$posteriorSettings
    priorSettings      <- res$priorSettings
    trialSimulation    <- res$trialSimulation
    save(simulationSettings, posteriorSettings, priorSettings, trialSimulation,
         file = fp)
    cat(sprintf("  saved %s\n", fp))
  }
}

# ----------------------------------------------------------------------------
# Phase 2: sweep (K, p, q, criterion) grid for one (hypothesis, binding) cell
# ----------------------------------------------------------------------------

build_design <- function(K, p_eff, p_fut_low, crit, binding) {
  looks <- make_looks(K)

  if (binding) {
    futility_thr <- list(size = 0, type = "absolute.risk")
    fut_prob_mat <- matrix(p_fut_low, nrow = K, ncol = 1)
    # Manuscript decision criteria (Section 2.1): futility evaluated only at interim looks
    # (k = 1, ..., K - 1); deactivate the final-stage row by setting its
    # threshold to 0 (pp < 0 is never true). At K = 1, this disables
    # futility entirely, consistent with the fixed-sample design.
    fut_prob_mat[K, 1L] <- 0
  } else {
    futility_thr <- NULL
    fut_prob_mat <- NULL
  }

  if (crit == "single") {
    eff_thr      <- list(size = 0, type = "absolute.risk")
    eff_prob_mat <- matrix(p_eff, nrow = K, ncol = 1)
  } else {
    eff_thr      <- list(size = c(-0.05, 0.00), type = "absolute.risk")
    eff_prob_mat <- matrix(rep(c(0.50, p_eff), each = K), nrow = K, ncol = 2)
  }

  list(
    numberOfSubjects = looks,
    effectThresholds = list(
      efficacy    = eff_thr,
      equivalence = NULL,
      futility    = futility_thr
    ),
    probabilityThresholds = list(
      efficacy    = eff_prob_mat,
      equivalence = NULL,
      futility    = fut_prob_mat
    )
  )
}

CRIT_LEVELS <- c("single", "dual")

sweep_grid <- function(trialSimulation, hypothesis, binding) {
  dims <- c(length(CRIT_LEVELS), length(K_GRID), length(P_EFF_GRID), length(P_FUT_GRID))
  errorRate     <- array(NA_real_, dim = dims)
  earlyStop     <- array(NA_real_, dim = dims)
  expSampleSize <- array(NA_real_, dim = dims)

  for (cc in seq_along(CRIT_LEVELS)) for (kk in seq_along(K_GRID)) for (pp in seq_along(P_EFF_GRID)) for (qq in seq_along(P_FUT_GRID)) {
    design <- build_design(
      K          = K_GRID[kk],
      p_eff      = P_EFF_GRID[pp],
      p_fut_low  = P_FUT_GRID[qq],
      crit       = CRIT_LEVELS[cc],
      binding    = binding
    )
    oc <- getOperatingCharacteristics(
      trialSimulation = trialSimulation,
      trialDesign     = design
    )
    errorRate[cc, kk, pp, qq]     <- if (hypothesis == "H0") oc$type1ErrorRate else oc$type2ErrorRate
    earlyStop[cc, kk, pp, qq]     <- oc$earlyStoppingProb
    expSampleSize[cc, kk, pp, qq] <- oc$expSampleSize
  }

  list(errorRate = errorRate, earlyStop = earlyStop, expSampleSize = expSampleSize)
}

# ----------------------------------------------------------------------------
# Phase 3: 12-panel plotting helpers
# ----------------------------------------------------------------------------

build_panel <- function(arr_cc, kk_labels, p_eff, q_eff_thr,
                        ylab_text, ylim_range, hline_y = NULL,
                        percent_axis = FALSE) {
  df <- data.frame(
    K                 = rep(kk_labels, times = length(CRIT_LEVELS)),
    yval              = c(arr_cc[1, ], arr_cc[2, ]),
    efficacyCriterion = factor(rep(CRIT_LEVELS, each = length(kk_labels)),
                                levels = CRIT_LEVELS,
                                labels = c("single-criterion", "dual-criterion"))
  )
  # Colourblind-safe (Okabe-Ito) palette, additionally distinguished by
  # linetype and point shape so the two series remain separable in greyscale
  # and for readers with colour-vision deficiency. The colour, linetype and
  # shape scales share a label so they collapse into a single legend.
  g <- ggplot(df, aes(x = K, y = yval, color = efficacyCriterion,
                      linetype = efficacyCriterion, shape = efficacyCriterion)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("single-criterion" = "#0072B2",
                                  "dual-criterion"   = "#D55E00")) +
    scale_linetype_manual(values = c("single-criterion" = "solid",
                                     "dual-criterion"    = "longdash")) +
    scale_shape_manual(values = c("single-criterion" = 16,
                                  "dual-criterion"    = 17)) +
    scale_x_continuous(breaks = 1:10, limits = c(1, 10)) +
    ggtitle(sprintf("p = %.2f, q = %.2f", p_eff, q_eff_thr)) +
    xlab("K") + ylab(ylab_text) +
    labs(color = "efficacy rule", linetype = "efficacy rule",
         shape = "efficacy rule") +
    theme_bw() +
    theme(legend.position = "bottom")
  if (isTRUE(percent_axis)) {
    g <- g + scale_y_continuous(
      labels = scales::label_percent(accuracy = 0.1),
      limits = ylim_range
    )
  } else {
    g <- g + scale_y_continuous(limits = ylim_range)
  }
  if (!is.null(hline_y)) {
    g <- g + geom_hline(yintercept = hline_y, linetype = "dashed",
                        color = "black", linewidth = 1)
  }
  g
}

build_12panel <- function(arr, ylab_text, ylim_range, hline_y = NULL,
                          percent_axis = FALSE) {
  panels <- vector("list", length(P_EFF_GRID) * length(P_FUT_GRID))
  for (i in seq_along(P_EFF_GRID)) {
    for (j in seq_along(P_FUT_GRID)) {
      idx <- (i - 1L) * length(P_FUT_GRID) + j
      arr_cc <- arr[, , i, j]                       # dims: 2 (criterion) x length(K_GRID)
      panels[[idx]] <- build_panel(
        arr_cc       = arr_cc,
        kk_labels    = K_GRID,
        p_eff        = P_EFF_GRID[i],
        q_eff_thr    = Q_FUT_GRID[j],
        ylab_text    = ylab_text,
        ylim_range   = ylim_range,
        hline_y      = hline_y,
        percent_axis = percent_axis
      )
    }
  }
  p <- panels[[1]]
  for (k in 2:length(panels)) p <- p + panels[[k]]
  p <- p + plot_layout(ncol = 3, nrow = 4, byrow = FALSE, guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom")
  p
}

# ----------------------------------------------------------------------------
# Phase 4: drive the plots for each (hypothesis, binding) configuration
# ----------------------------------------------------------------------------

# The early-stopping probability under the alternative and the null is
# computed by sweep_grid() and stored in the swp_futON$earlyStop array for
# downstream inspection, but the corresponding 12-panel figures are not
# rendered. They are not referenced in the manuscript or in the supplement,
# and they are byte-identical between the binding and non-binding
# configurations (by the convention of Section 3.4 of the main paper). If
# rendering is needed for ad-hoc inspection, add the corresponding entries
# below.
PLOT_DEFINITIONS <- list(
  H0 = list(
    list(stat = "errorRate",     ylab = "Type I error rate",          ylim = c(0,    0.32), hline = 0.025,
         basename = "type1Error",       percent_axis = TRUE),
    list(stat = "expSampleSize", ylab = "Expected sample size E(N)",  ylim = c(2000, 3800), hline = 3800,
         basename = "expSampleSizeH0",  percent_axis = FALSE)
  ),
  H1 = list(
    list(stat = "errorRate",     ylab = "Type II error rate",         ylim = c(0,    0.55), hline = 0.10,
         basename = "type2Error",       percent_axis = TRUE),
    list(stat = "expSampleSize", ylab = "Expected sample size E(N)",  ylim = c(900,  3800), hline = 3800,
         basename = "expSampleSizeH1",  percent_axis = FALSE)
  )
)

for (hyp in names(PLOT_DEFINITIONS)) {
  cat(sprintf("\n=== Operating characteristics under %s ===\n", hyp))
  rda_path <- file.path(OUTPUT_DIR, sprintf("ADRENAL_OperatingCharacteristics_%s.rda", hyp))
  e <- new.env(); load(rda_path, envir = e)

  # Under standard (FDA/EMA) non-binding semantics the futility boundary is
  # still followed in practice, so type II, expected sample size and early
  # stopping probability are computed with the futility rule active. Only the
  # type I calculation drops the futility rule, to give a conservative bound
  # on the false-positive rate if the DSMB chose to override futility under
  # H0. We therefore need at most two sweeps per hypothesis: one with futility
  # active (used for every metric in the binding figures, and for all
  # non-binding figures except type I), and one with futility removed (used
  # only for the non-binding type I figure under H0).
  cat("  sweep with futility active ...\n")
  t0 <- proc.time()
  swp_futON <- sweep_grid(e$trialSimulation, hypothesis = hyp, binding = TRUE)
  cat(sprintf("     sweep elapsed = %.1f s\n", (proc.time() - t0)[["elapsed"]]))

  swp_futOFF <- NULL
  if (hyp == "H0") {
    cat("  sweep with futility removed (for non-binding type I) ...\n")
    t0 <- proc.time()
    swp_futOFF <- sweep_grid(e$trialSimulation, hypothesis = hyp, binding = FALSE)
    cat(sprintf("     sweep elapsed = %.1f s\n", (proc.time() - t0)[["elapsed"]]))
  }

  for (binding in c(TRUE, FALSE)) {
    binding_tag <- if (binding) "Binding" else "nonBinding"
    cat(sprintf("  -- %s futility figures --\n", if (binding) "binding" else "non-binding"))

    for (cfg in PLOT_DEFINITIONS[[hyp]]) {
      # Only non-binding type I under H0 drops the futility rule. Every other
      # (hypothesis, metric, binding) combination uses the futility-active
      # sweep, so type II / E[N|H0] / E[N|H1] / early stopping probability
      # match exactly between the binding and non-binding figures.
      use_futOFF <- (!binding) && hyp == "H0" && cfg$stat == "errorRate"
      data_arr   <- if (use_futOFF) swp_futOFF[[cfg$stat]] else swp_futON[[cfg$stat]]

      file_path <- file.path(
        OUTPUT_DIR,
        sprintf("ADRENAL_OperatingCharacteristics_%s_%s.jpeg", cfg$basename, binding_tag)
      )
      p <- build_12panel(
        arr          = data_arr,
        ylab_text    = cfg$ylab,
        ylim_range   = cfg$ylim,
        hline_y      = cfg$hline,
        percent_axis = isTRUE(cfg$percent_axis)
      )
      ggsave(filename = file_path, plot = p, width = 15, height = 12, dpi = 300)
      cat(sprintf("     saved %s\n", file_path))
    }
  }
}

cat("\nDone. Operating-characteristic behaviour evaluation complete.\n")

# Persist sessionInfo and package versions alongside the cache.
bayesgsd_save_session("ADRENAL_Evaluation_OperatingCharacteristics")
