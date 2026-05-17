#' @title Monte Carlo precision of the operating-characteristic estimates for
#'   the ADRENAL re-design (Section 3.5)
#' @description Estimates the sampling distribution of the type I and type II
#'   error rate estimators for Bayesian GSDs with five design skeletons
#'   K in {1, 3, 5, 7, 9} (equivalently 0, 2, 4, 6 and 8 equally spaced
#'   interim analyses), by repeating `runTrialMonitoring` 100 times at each
#'   of four virtual-trial budgets R in {1,000, 10,000, 100,000, 1,000,000}
#'   under both H0 and H1. Replicates within a single (R, K, hypothesis)
#'   cell use distinct seeds, and seed segments are disjoint across R values,
#'   so the 100 replicates at each R are independent of one another and of
#'   those at any other R. The decision rule fixes a single-criterion
#'   efficacy threshold p = 0.99 and a binding futility threshold q = 0.90,
#'   with an independent Beta(1,1) prior on each arm. The cache underlies
#'   the numerical summary table and the boxplot figure of Section 3.5.
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



suppressPackageStartupMessages({
  library(RBesT)
  library(parallel)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
})

source("./Code/Code v1.0/bayseqSim_bern.R")

OUTPUT_DIR <- "./Output/Output v1.0"
RDA_PATH   <- file.path(OUTPUT_DIR, "ADRENAL_MonteCarloError.rda")
FIG_PATH   <- file.path(OUTPUT_DIR, "ADRENAL_MonteCarloError.jpeg")

# ----------------------------------------------------------------------------
# Trial parameters
# ----------------------------------------------------------------------------

P_CONTROL  <- 0.33
DELTA_ALT  <- -0.05
P_TREAT_H1 <- P_CONTROL + DELTA_ALT
N_MAX      <- 3800L

P_EFF      <- 0.99
Q_FUT_MS   <- 0.90
P_FUT_LOW  <- 1 - Q_FUT_MS

# ----------------------------------------------------------------------------
# Simulation grid: virtual-trial budgets, design skeletons, replicates
# ----------------------------------------------------------------------------

R_GRID   <- c(1000L, 10000L, 100000L, 1000000L)
R_LABELS <- c("1,000", "10,000", "100,000", "1,000,000")
N_REPS   <- 100L

# Disjoint seed segments across R values: replicate j at R = R_GRID[r_ix]
# uses seed = SEED_OFFSET[r_ix] + j. With N_REPS = 100, this gives seed
# ranges 1..100 (R=1,000), 101..200 (R=10,000), 201..300 (R=100,000) and
# 301..400 (R=1,000,000), so the 100 replicates at each R are mutually
# independent and independent of those at any other R.
SEED_OFFSET <- as.integer(seq(0L, by = N_REPS, length.out = length(R_GRID)))

K_GRID <- c(1L, 3L, 5L, 7L, 9L)   # total analyses K; interim analyses = K - 1

make_looks <- function(K) {
  as.integer(ceiling(seq(from = 0, to = N_MAX, length.out = K + 1L))[-1])
}

DESIGNS <- list()
for (K in K_GRID) {
  DESIGNS[[sprintf("K%d", K)]] <- list(K = K, looks = make_looks(K))
}
cat("Designs:\n")
for (d in DESIGNS) cat(sprintf("  K=%d (interim analyses = %d) looks=%s\n",
                                d$K, d$K - 1L, paste(d$looks, collapse = ",")))
cat("R budgets:    ", paste(R_LABELS, collapse = ", "), "\n")
cat("Replicates:   ", N_REPS, "per (R, K, hypothesis) cell\n")
cat("Seed offsets: ", paste(SEED_OFFSET, collapse = ", "), "\n")

# ----------------------------------------------------------------------------
# Runner: one replicate of a (R, K, hypothesis) cell
# ----------------------------------------------------------------------------

run_one_replicate <- function(design, hypothesis, R, seed) {
  rates <- if (hypothesis == "H0") c(P_CONTROL, P_CONTROL) else c(P_CONTROL, P_TREAT_H1)
  simulationSettings <- list(
    numberOfTrials   = R,
    numberOfSubjects = design$looks,
    allocationRatio  = 1,
    rates            = rates,
    hypothesis       = if (hypothesis == "H0") "null" else "alternative",
    seed             = seed
  )
  posteriorSettings <- list(
    effectThresholds = list(
      efficacy    = list(size = 0, type = "absolute.risk"),
      equivalence = NULL,
      futility    = list(size = 0, type = "absolute.risk")
    ),
    alternative = "less"
  )
  priorSettings <- initialisePriorSettings()
  # Manuscript Eq. 3.2.3: futility evaluated only at interim looks
  # (k = 1, ..., K - 1); deactivate the final-stage row. Final-stage
  # outcomes are recorded as "no efficacy" rather than "futility", but
  # for the type I / type II error rate reported here the result is
  # identical to the legacy "futility-on-at-K" path.
  fut_thr <- matrix(P_FUT_LOW, nrow = design$K, ncol = 1)
  fut_thr[design$K, 1L] <- 0
  trialDesign <- list(
    numberOfSubjects = design$looks,
    effectThresholds = list(
      efficacy    = list(size = 0, type = "absolute.risk"),
      equivalence = NULL,
      futility    = list(size = 0, type = "absolute.risk")
    ),
    probabilityThresholds = list(
      efficacy    = matrix(P_EFF, nrow = design$K, ncol = 1),
      equivalence = NULL,
      futility    = fut_thr
    )
  )

  trialSimulation <- runTrialMonitoring(
    simulationSettings = simulationSettings,
    posteriorSettings  = posteriorSettings,
    priorSettings      = priorSettings
  )
  oc <- getOperatingCharacteristics(
    trialSimulation = trialSimulation,
    trialDesign     = trialDesign
  )
  if (hypothesis == "H0") oc$type1ErrorRate else oc$type2ErrorRate
}

# ----------------------------------------------------------------------------
# Cache: lists of (length(K_GRID) x N_REPS) matrices keyed by R_LABEL
# ----------------------------------------------------------------------------

init_matrix_list <- function() {
  out <- vector("list", length(R_LABELS))
  names(out) <- R_LABELS
  for (lab in R_LABELS) {
    out[[lab]] <- matrix(NA_real_, nrow = length(K_GRID), ncol = N_REPS,
                         dimnames = list(sprintf("K%d", K_GRID), NULL))
  }
  out
}

if (file.exists(RDA_PATH)) {
  load(RDA_PATH)
  if (!is.list(type1ErrorRate) ||
      !identical(names(type1ErrorRate), R_LABELS) ||
      !identical(dim(type1ErrorRate[[1]]), c(length(K_GRID), N_REPS))) {
    stop("Cache at ", RDA_PATH, " is not in the direct-simulation format ",
         "(expected lists indexed by R = ", paste(R_LABELS, collapse = ", "),
         " with ", length(K_GRID), " x ", N_REPS, " matrices). ",
         "Delete it before running this script.")
  }
  cat(sprintf("\nResuming from %s with %d filled cells (max %d).\n",
              RDA_PATH,
              sum(!is.na(unlist(type1ErrorRate)) & !is.na(unlist(type2ErrorRate))),
              length(R_LABELS) * length(K_GRID) * N_REPS))
} else {
  type1ErrorRate     <- init_matrix_list()
  type2ErrorRate     <- init_matrix_list()
  numberOfAnalyses   <- K_GRID
  numberOfReplicates <- N_REPS
}

save_cache <- function() {
  save(numberOfAnalyses, numberOfReplicates, R_GRID, R_LABELS, SEED_OFFSET,
       type1ErrorRate, type2ErrorRate,
       file = RDA_PATH)
}

# ----------------------------------------------------------------------------
# Main loop: for each R, K, hypothesis, run N_REPS replicates with
# disjoint seeds. Save after every replicate to make resume cheap.
# ----------------------------------------------------------------------------

for (r_ix in seq_along(R_GRID)) {
  R       <- R_GRID[r_ix]
  R_label <- R_LABELS[r_ix]
  cat(sprintf("\n========== R = %s ==========\n", R_label))
  for (i in seq_along(K_GRID)) {
    design <- DESIGNS[[i]]
    cat(sprintf("\n=== R=%s, K=%d (interim analyses = %d) ===\n",
                R_label, design$K, design$K - 1L))
    for (j in seq_len(N_REPS)) {
      done_h0 <- !is.na(type1ErrorRate[[R_label]][i, j])
      done_h1 <- !is.na(type2ErrorRate[[R_label]][i, j])
      if (done_h0 && done_h1) next

      seed_j <- SEED_OFFSET[r_ix] + j

      if (!done_h0) {
        type1ErrorRate[[R_label]][i, j] <-
          run_one_replicate(design, "H0", R = R, seed = seed_j)
      }
      if (!done_h1) {
        type2ErrorRate[[R_label]][i, j] <-
          run_one_replicate(design, "H1", R = R, seed = seed_j)
      }
      save_cache()

      if (j %% 10L == 0L || j == N_REPS) {
        cat(sprintf("  replicate %3d / %3d  (alpha=%.4f, beta=%.4f)\n",
                    j, N_REPS,
                    type1ErrorRate[[R_label]][i, j],
                    type2ErrorRate[[R_label]][i, j]))
      }
    }
  }
}

save_cache()
cat(sprintf("\nSaved %s.\n", RDA_PATH))

# ----------------------------------------------------------------------------
# Boxplot figure (Figure 3.5.1)
# ----------------------------------------------------------------------------

K_LABELS <- as.character(K_GRID)

build_boxplot_df <- function(errorList) {
  df <- NULL
  for (i in seq_along(K_LABELS)) {
    for (R_label in R_LABELS) {
      rates <- errorList[[R_label]][i, ]
      df <- rbind(df,
                  data.frame(K              = K_LABELS[i],
                             numberOfTrials = R_label,
                             errorRate      = rates))
    }
  }
  df$numberOfTrials <- factor(df$numberOfTrials, levels = R_LABELS)
  df
}

df_t1 <- build_boxplot_df(type1ErrorRate)
df_t2 <- build_boxplot_df(type2ErrorRate)

p1 <- ggplot(df_t1, aes(x = K, y = errorRate, fill = numberOfTrials)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  scale_fill_bmj() +
  xlab("K") + ylab("type I error rate") +
  labs(fill = "virtual trials") +
  theme_bw() +
  theme(legend.position = "bottom")

p2 <- ggplot(df_t2, aes(x = K, y = errorRate, fill = numberOfTrials)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  scale_fill_bmj() +
  xlab("K") + ylab("type II error rate") +
  labs(fill = "virtual trials") +
  theme_bw() +
  theme(legend.position = "bottom")

p <- (p1 + p2) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

ggsave(filename = FIG_PATH, plot = p, width = 12, height = 6, dpi = 300)
cat(sprintf("Saved %s.\n", FIG_PATH))

cat("\nDone. Monte Carlo precision evaluation complete.\n")

# Persist sessionInfo and package versions alongside the cache.
bayesgsd_save_session("ADRENAL_Evaluation_MonteCarloError")
