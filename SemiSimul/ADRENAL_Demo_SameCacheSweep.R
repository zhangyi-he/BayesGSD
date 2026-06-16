#' @title Same-cache demonstration of the precomputation strategy of Section 2.4 (Figure 1, Section 3.4)
#' @description Illustrates Lemma~1 (cache invariance) operationally: from a
#'   single cached Monte Carlo pass at the union of all candidate look times,
#'   the 438-design grid of Section 3.4 is evaluated one design at a time and
#'   the cumulative wall-clock is recorded after each design. The resulting
#'   curve shows the architectural claim directly --- a one-time simulation
#'   cost (the y-intercept) followed by a shallow, near-linear per-design
#'   sweep cost --- against the per-design re-simulation that BATSS and adaptr
#'   would require.
#'
#'   The designs sweep the same axes as
#'   ADRENAL_Evaluation_OperatingCharacteristics.R: K in 1..10, efficacy
#'   thresholds p in {0.90, 0.95, 0.99}, futility thresholds q in {0.85, 0.90,
#'   0.95, 1.00}, single-/dual-criterion efficacy, and binding/non-binding
#'   futility, for a total of 438 configurations, all evaluated against the
#'   one cached simulation per hypothesis.
#'
#'   Output: `ADRENAL_Demo_SameCacheSweep.{rda,jpeg}` (the figure is Figure 1
#'   of Section 3.4, the same-cache demonstration of Lemma~1) plus a
#'   sessionInfo dump.
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

suppressPackageStartupMessages({
  library(RBesT); library(parallel)
  library(ggplot2); library(ggsci)
})
source("./Code/Code v1.0/adabay_bern.R")

RDA_PATH <- file.path(OUTPUT_DIR, "ADRENAL_Demo_SameCacheSweep.rda")
FIG_PATH <- file.path(OUTPUT_DIR, "ADRENAL_Demo_SameCacheSweep.jpeg")

OC_H0 <- file.path(OUTPUT_DIR, "ADRENAL_OperatingCharacteristics_H0.rda")
OC_H1 <- file.path(OUTPUT_DIR, "ADRENAL_OperatingCharacteristics_H1.rda")
FGT   <- file.path(OUTPUT_DIR, "ADRENAL_Benchmark_FullGridTiming.rda")
for (f in c(OC_H0, OC_H1)) if (!file.exists(f))
  stop("Required OC cache missing: ", f,
       "\nRun ADRENAL_Evaluation_OperatingCharacteristics.R first.")

# One-time simulation cost (the y-intercept) from the timing cache, if present.
sim_cost <- c(H0 = NA_real_, H1 = NA_real_)
if (file.exists(FGT)) {
  e <- new.env(parent = emptyenv()); load(FGT, envir = e)
  sim_cost["H0"] <- e$timing$H0$elapsed
  sim_cost["H1"] <- e$timing$H1$elapsed
}

# ----------------------------------------------------------------------------
# Design grid (identical axes to ADRENAL_Evaluation_OperatingCharacteristics.R)
# ----------------------------------------------------------------------------

N_MAX       <- 3800L
make_looks  <- function(K) as.integer(ceiling(seq(0, N_MAX, length.out = K + 1L))[-1])
P_EFF_GRID  <- c(0.90, 0.95, 0.99)
P_FUT_LOW   <- c(0.15, 0.10, 0.05, 0.00)   # 1 - q for q in {0.85,0.90,0.95,1.00}
K_GRID      <- 1:10
CRIT_LEVELS <- c("single", "dual")

build_design <- function(K, p_eff, p_fut_low, crit, binding) {
  looks <- make_looks(K)
  if (binding) {
    futility_thr <- list(size = 0, type = "absolute.risk")
    fut_prob_mat <- matrix(p_fut_low, nrow = K, ncol = 1); fut_prob_mat[K, 1L] <- 0
  } else {
    futility_thr <- NULL; fut_prob_mat <- NULL
  }
  if (crit == "single") {
    eff_thr <- list(size = 0, type = "absolute.risk")
    eff_prob_mat <- matrix(p_eff, nrow = K, ncol = 1)
  } else {
    eff_thr <- list(size = c(-0.05, 0.00), type = "absolute.risk")
    eff_prob_mat <- matrix(rep(c(0.50, p_eff), each = K), nrow = K, ncol = 2)
  }
  list(numberOfSubjects = looks,
       effectThresholds = list(efficacy = eff_thr, equivalence = NULL,
                               futility = futility_thr),
       probabilityThresholds = list(efficacy = eff_prob_mat, equivalence = NULL,
                                    futility = fut_prob_mat))
}

# Enumerate the 438 designs (K=1 contributes 6; K=2..10 contribute 432).
designs <- list()
for (K in K_GRID) for (p_eff in P_EFF_GRID) for (crit in CRIT_LEVELS) {
  if (K == 1L) {
    designs[[length(designs) + 1L]] <- build_design(K, p_eff, 0, crit, TRUE)
  } else {
    for (p_fut in P_FUT_LOW) for (binding in c(TRUE, FALSE))
      designs[[length(designs) + 1L]] <- build_design(K, p_eff, p_fut, crit, binding)
  }
}
n_designs <- length(designs)
stopifnot(n_designs == 438L)
cat(sprintf("Enumerated %d designs.\n", n_designs))

# ----------------------------------------------------------------------------
# Sweep: evaluate each design against the single cached simulation, recording
# cumulative wall-clock after each design.
# ----------------------------------------------------------------------------

sweep_cumulative <- function(rda_path, hyp) {
  e <- new.env(parent = emptyenv()); load(rda_path, envir = e)
  ts <- e$trialSimulation
  cum <- numeric(n_designs)
  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(n_designs)) {
    invisible(getOperatingCharacteristics(ts, designs[[i]]))
    cum[i] <- proc.time()[["elapsed"]] - t0
  }
  data.frame(hypothesis = hyp, design_index = seq_len(n_designs),
             sweep_seconds = cum)
}

cat("Sweeping 438 designs against the H0 cache ...\n")
df_H0 <- sweep_cumulative(OC_H0, "H0")
cat(sprintf("  H0 sweep total: %.1f s\n", tail(df_H0$sweep_seconds, 1)))
cat("Sweeping 438 designs against the H1 cache ...\n")
df_H1 <- sweep_cumulative(OC_H1, "H1")
cat(sprintf("  H1 sweep total: %.1f s\n", tail(df_H1$sweep_seconds, 1)))

# Cumulative wall-clock including the one-time simulation cost (the intercept).
df <- rbind(df_H0, df_H1)
df$total_seconds <- df$sweep_seconds +
  ifelse(df$hypothesis == "H0", sim_cost["H0"], sim_cost["H1"])

# ----------------------------------------------------------------------------
# Figure: cumulative wall-clock vs number of designs evaluated
# ----------------------------------------------------------------------------

per_design_H0 <- tail(df_H0$sweep_seconds, 1) / n_designs
cat(sprintf("Per-design sweep cost: H0 %.3f s, H1 %.3f s.\n",
            per_design_H0, tail(df_H1$sweep_seconds, 1) / n_designs))

p <- ggplot(df, aes(x = design_index, y = total_seconds, color = hypothesis)) +
  geom_line(linewidth = 1) +
  scale_color_bmj() +
  scale_x_continuous(breaks = c(1, seq(50, n_designs, by = 50), n_designs)) +
  xlab("Number of candidate designs evaluated against the single cache") +
  ylab("Cumulative wall-clock (seconds)") +
  labs(color = "Hypothesis") +
  ggtitle("Same-cache evaluation of the 438-design grid") +
  theme_bw() + theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 300))
if (all(is.finite(sim_cost))) {
  p <- p + geom_hline(yintercept = sim_cost[["H0"]], linetype = "dotted",
                      color = "grey40") +
    annotate("text", x = n_designs * 0.55, y = sim_cost[["H0"]],
             label = "one-time simulation pass", hjust = 0.5, vjust = -0.6,
             size = 3, color = "grey30")
}
ggsave(FIG_PATH, p, width = 9, height = 5.5, dpi = 200)
cat(sprintf("Saved %s\n", FIG_PATH))

save(df, df_H0, df_H1, sim_cost, n_designs, per_design_H0, file = RDA_PATH)
cat(sprintf("Wrote %s\n", RDA_PATH))

bayesgsd_save_session("ADRENAL_Demo_SameCacheSweep")
