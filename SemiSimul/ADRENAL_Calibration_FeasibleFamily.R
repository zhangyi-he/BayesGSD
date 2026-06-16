#' @title Feasible interim/final threshold family for the headline designs
#'   (Section 3.5)
#' @description Section 3.5 frames calibration as selecting a point on the
#'   feasible family of efficacy-threshold pairs. This script makes that family
#'   explicit for the headline three- and five-look designs.
#'
#'   The family is the set of single-criterion Haybittle--Peto (common interim,
#'   final) threshold pairs at the fixed futility threshold q = 0.90 that meet
#'   the design targets: power >= 90% at the design alternative and a type I
#'   error SUPREMUM at most 2.5% over the plausible control-rate band
#'   vartheta_c in {0.15,...,0.50} (the non-binding, DSMC-always-overrides
#'   convention). Design-null operating characteristics for every pair are read
#'   from the grid-search cache produced by `ADRENAL_Calibration_GridSearch.R`;
#'   the band supremum --- the operative type I constraint --- is then obtained
#'   by one simulation pass per control rate (Lemma 1: the data differ with
#'   vartheta_c, so each rate is simulated afresh, but every pair is scored
#'   against that single pass). Pairs already failing the design-null target
#'   have a band supremum above 2.5% a fortiori and need no simulation.
#'
#'   The figure marks the maximum-power and minimum-E(N|H1) feasible designs and
#'   the selected headline design (interim SEL_INT, final SEL_FIN), which are the
#'   maximum-power frontier points defined once in adabay_bern_setup.R
#'   (BAYESGSD_DESIGNS) --- the single source of truth for the Section 3.5
#'   calibrated designs. For the headline three- and five-look designs the
#'   selected point coincides with the maximum-power design; the plotting code
#'   handles either case (coincident or distinct) automatically.
#'
#'   Design-null operating characteristics for every pair are read from the
#'   durable grid-search cache `ADRENAL_CalibratedDesign.rda` (object
#'   `hp_batches`) produced by `ADRENAL_Calibration_GridSearch.R`.
#'
#'   Output: prints the band-feasible family size per K; writes
#'   `ADRENAL_Calibration_FeasibleFamily.rda` and the figure
#'   `ADRENAL_Calibration_FeasibleFamily.jpeg` (the Section 3.5 figure).
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

suppressPackageStartupMessages({
  library(RBesT); library(parallel); library(ggplot2)
})
# ggrepel is used below for geom_text_repel; require it explicitly with an
# actionable message rather than failing opaquely at plot time.
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  stop("Package 'ggrepel' is required for the feasible-family figure ",
       "(geom_text_repel). Install it with install.packages(\"ggrepel\").")
}
source("./Code/Code v1.0/adabay_bern.R")

# Read the per-K Haybittle-Peto batches from the DURABLE grid-search cache
# (ADRENAL_CalibratedDesign.rda, object hp_batches) rather than from the resume
# checkpoint, which ADRENAL_Calibration_GridSearch.R unlinks on success.
CACHE    <- file.path(OUTPUT_DIR, "ADRENAL_CalibratedDesign.rda")
RDA_PATH <- file.path(OUTPUT_DIR, "ADRENAL_Calibration_FeasibleFamily.rda")
FIG_PATH <- file.path(OUTPUT_DIR, "ADRENAL_Calibration_FeasibleFamily.jpeg")

# N_MAX, make_looks(), SEED_DEFAULT and the selected design thresholds are the
# single source of truth in adabay_bern_setup.R (sourced via the bootstrap).
KS             <- c(3L, 5L)
ALPHA          <- 0.025                          # type I target (as a band supremum)
POWER          <- 0.90                           # power target at the design alternative
SEL_INT        <- BAYESGSD_DESIGNS$interim_eff   # selected design: interim threshold (0.999)
SEL_FIN        <- BAYESGSD_DESIGNS$final_by_K[["3"]]  # selected design: final threshold (0.977)
ALPHA_PLOT_MAX <- 0.028      # plot window: show the type I > 2.5% boundary up to here
POWER_PLOT_MIN <- 0.89       # plot window: show the power < 90% boundary down to here
SEED           <- SEED_DEFAULT
R_TRIALS       <- as.integer(Sys.getenv("R_TRIALS", unset = "1000000"))
VARTHETA       <- c(0.15, 0.20, 0.25, 0.30, 0.33, 0.36, 0.40, 0.45, 0.50)

# ---- design-null operating characteristics for all pairs (from the cache) ----
if (!file.exists(CACHE))
  stop("Run ADRENAL_Calibration_GridSearch.R first to produce ", CACHE,
       " (object hp_batches).")
e <- new.env(); load(CACHE, envir = e)
if (is.null(e$hp_batches))
  stop(CACHE, " does not contain hp_batches; re-run ADRENAL_Calibration_GridSearch.R ",
       "to refresh the durable cache.")
get_pairs <- function(K) {
  df <- e$hp_batches[[paste0("hp_K", K)]]
  df <- df[df$criterion == "single" & abs(df$q - 0.90) < 1e-9, ]
  data.frame(K = K, p_interim = df$p_interim, p_final = df$p_final,
             dp_type1 = df$type1, power = df$power, EN_H1 = df$EN_H1)
}
allp <- do.call(rbind, lapply(KS, get_pairs))

# ---- band-sup (non-binding) type I for the design-null-feasible pairs ----
# Operative feasibility is the type I SUPREMUM over the control-rate band, not the
# design-null value. The data differ with vartheta_c, so each control rate is a
# fresh simulation (Lemma 1: one pass per rate, every pair scored against it).
simulate_null <- function(K, v) {
  ss <- list(numberOfTrials = R_TRIALS, numberOfSubjects = make_looks(K), allocationRatio = 1,
             rates = c(v, v), hypothesis = "null", seed = SEED)
  ps <- list(effectThresholds = list(efficacy = list(size = 0, type = "absolute.risk"),
             equivalence = NULL, futility = list(size = 0, type = "absolute.risk")),
             alternative = "less")
  runTrialMonitoring(ss, ps, priorSettings = initialisePriorSettings())
}
score_nb <- function(sim, K, pin, pfn) {   # non-binding (futility-disabled) type I rate
  pv <- c(rep(pin, K - 1L), pfn)
  td <- list(numberOfSubjects = make_looks(K),
    effectThresholds = list(efficacy = list(size = 0, type = "absolute.risk"), equivalence = NULL,
      futility = list(size = 0, type = "absolute.risk")),
    probabilityThresholds = list(efficacy = matrix(pv, nrow = K, ncol = 1), equivalence = NULL,
      futility = matrix(0, nrow = K, ncol = 1)))
  getOperatingCharacteristics(sim, td)$type1ErrorRate
}
allp$sup_nb <- NA_real_
for (K in KS) {
  idx  <- which(allp$K == K & allp$dp_type1 <= ALPHA & allp$power >= POWER)
  sims <- lapply(VARTHETA, function(v) simulate_null(K, v))
  for (i in idx)
    allp$sup_nb[i] <- max(vapply(sims, function(sm)
      score_nb(sm, K, allp$p_interim[i], allp$p_final[i]), numeric(1)))
}
allp$bandfeas <- !is.na(allp$sup_nb) & allp$sup_nb <= ALPHA & allp$power >= POWER

for (K in KS)
  cat(sprintf("K=%d: %d band-feasible pairs (band-sup type I <= 2.5%%, power >= 90%%)\n",
              K, sum(allp$K == K & allp$bandfeas)))

# ---- classify for the plot: filled = band-feasible; x = band-sup > 2.5%;
#      open triangle = power < 90% ----
plt <- allp[allp$power >= POWER_PLOT_MIN & allp$dp_type1 <= ALPHA_PLOT_MAX, ]
plt$Status <- with(plt, ifelse(power < POWER, "power < 90%",
                        ifelse(bandfeas, "Feasible", "type I > 2.5%")))
plt$Status <- factor(plt$Status, levels = c("Feasible", "type I > 2.5%", "power < 90%"))
plt$Klab   <- factor(paste0(plt$K, " looks"))

# ---- objective corners among band-feasible designs: max power, min E(N|H1),
#      and the selected design (supplement Table S3; coincides with max power at K=5) ----
corners <- do.call(rbind, lapply(KS, function(k) {
  bf  <- allp[allp$K == k & allp$bandfeas, ]
  mp  <- bf[which.max(bf$power), ]
  me  <- bf[which.min(bf$EN_H1), ]
  sel <- bf[abs(bf$p_interim - SEL_INT) < 1e-9 & abs(bf$p_final - SEL_FIN) < 1e-9, ]
  same <- nrow(sel) > 0 && abs(mp$p_interim - sel$p_interim) < 1e-9 &&
          abs(mp$p_final - sel$p_final) < 1e-9
  if (same) {
    r <- data.frame(p_interim = c(sel$p_interim, me$p_interim),
                    p_final   = c(sel$p_final,   me$p_final),
                    corner = c("max power (selected)", "min E(N)"),
                    nx = c(0.0004, -0.0004), ny = c(-0.0019, 0.0018))
  } else {
    r <- data.frame(p_interim = c(mp$p_interim, sel$p_interim, me$p_interim),
                    p_final   = c(mp$p_final,   sel$p_final,   me$p_final),
                    corner = c("max power", "selected", "min E(N)"),
                    nx = c(-0.0006, 0.0006, -0.0004), ny = c(-0.0021, 0.0010, 0.0018))
  }
  cbind(Klab = factor(paste0(k, " looks")), r)
}))

save(allp, plt, corners, KS, ALPHA, POWER, SEL_INT, SEL_FIN, VARTHETA, R_TRIALS, SEED,
     file = RDA_PATH)
cat(sprintf("Wrote %s\n", RDA_PATH))

p <- ggplot(plt, aes(p_interim, p_final)) +
  geom_point(aes(colour = power * 100, shape = Status), size = 3.3) +
  geom_point(data = corners, shape = 21, size = 6.5, stroke = 1.2, colour = "black", fill = NA) +
  ggrepel::geom_text_repel(data = corners, aes(label = corner), size = 3,
    nudge_x = corners$nx, nudge_y = corners$ny, min.segment.length = 0,
    segment.size = 0.3, segment.colour = "grey45", box.padding = 0.3, seed = 1) +
  facet_wrap(~ Klab) +
  coord_cartesian(xlim = c(0.993, 0.9997), ylim = c(0.975, 0.983)) +
  scale_colour_viridis_c(name = "Power (%)") +
  scale_shape_manual(values = c("Feasible" = 16, "type I > 2.5%" = 4, "power < 90%" = 6), name = NULL) +
  labs(x = "Common interim efficacy threshold", y = "Final efficacy threshold") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom", legend.box = "vertical") +
  guides(colour = guide_colourbar(barwidth = 16, order = 1), shape = guide_legend(order = 2))
ggsave(FIG_PATH, p, width = 9, height = 5, dpi = 200)
cat(sprintf("Wrote %s\n", FIG_PATH))

bayesgsd_save_session("ADRENAL_Calibration_FeasibleFamily")
