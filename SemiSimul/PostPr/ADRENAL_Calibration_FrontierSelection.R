#' @title Maximum-power frontier across the number of looks (Section 3.5)
#' @description Post-processes the grid-search cache produced by
#'   `ADRENAL_Calibration_GridSearch.R` (no new simulation). With the maximum
#'   sample size and the futility boundary fixed, the expected sample size is a
#'   downstream consequence of the efficacy thresholds rather than a design
#'   objective; the calibration spends the type I budget and maximises power.
#'   For each number of looks K, among single-criterion designs at the common
#'   futility threshold q = 0.90, it selects the design that MAXIMISES power
#'   subject to the design-point type I error <= ALPHA_SELECT; maximising power
#'   drives the common interim threshold to its most stringent grid value (0.999).
#'   The cache type I estimate is at the single design null and is optimistically
#'   biased by the selection (winner's curse), and the operative constraint is
#'   the non-binding type I SUPREMUM over the control-rate band, so the per-K
#'   designs are CANDIDATES: each is confirmed <= 2.5% on its own analysis
#'   schedule over the band by ADRENAL_Calibration_PlausibilitySweep.R, with the
#'   final threshold raised one grid step if its band supremum exceeds the
#'   target. The band-verified per-K frontier is the table in Supplemental
#'   Material Section S6.
#'
#'   Input : `ADRENAL_CalibratedDesign.rda` (objects df_const, df_hp,
#'           ALPHA_TARGET, POWER_TARGET).
#'   Output: prints the frontier; writes `ADRENAL_Calibration_FrontierSelection.rda`
#'           (object df_frontier) plus a sessionInfo dump.
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

IN_PATH  <- file.path(OUTPUT_DIR, "ADRENAL_CalibratedDesign.rda")
RDA_PATH <- file.path(OUTPUT_DIR, "ADRENAL_Calibration_FrontierSelection.rda")
if (!file.exists(IN_PATH))
  stop("Run ADRENAL_Calibration_GridSearch.R first to produce ", IN_PATH)

e <- new.env(); load(IN_PATH, envir = e)
ALPHA_TARGET <- e$ALPHA_TARGET
POWER_TARGET <- e$POWER_TARGET
# Selection thresholds, optionally tightened below the reported targets to
# leave headroom for Monte Carlo noise and selection bias ("winner's curse"),
# so the chosen designs stay inside the targets on an independent own-looks
# check. Env overrides: ALPHA_SELECT, POWER_SELECT.
# Default to a small margin below the reported targets, chosen so the selected
# designs remain inside 2.5%/90% on an independent own-looks re-evaluation
# (Section 3.5 plausibility/seed checks); override via env to reproduce the
# unmargined frontier.
ALPHA_SELECT <- as.numeric(Sys.getenv("ALPHA_SELECT", unset = "0.0245"))
POWER_SELECT <- as.numeric(Sys.getenv("POWER_SELECT", unset = "0.902"))

# Unify the constant-threshold and Haybittle--Peto design families into one
# table with a common efficacy-threshold description.
dc <- e$df_const; dc$family <- "constant"
dc$efficacy <- sprintf("%.3f", dc$p)
dh <- e$df_hp;    dh$family <- "HaybittlePeto"
dh$efficacy <- sprintf("%.3f/%.3f", dh$p_interim, dh$p_final)
cols <- c("K", "family", "efficacy", "q", "criterion",
          "type1", "power", "EN_H0", "EN_H1")
all_designs <- rbind(dc[, cols], dh[, cols])

# Section 3.5 setting: single-criterion designs at the common futility threshold
# q = 0.90 (q is irrelevant to power; the K=1 fixed-sample design has no futility
# rule, q = NA). Per-K, take the MAXIMUM-power candidate with design-point type I
# <= ALPHA_SELECT; its band supremum is confirmed <= 2.5% by the plausibility
# sweep (final threshold raised one grid step if the band supremum exceeds it).
sel <- subset(all_designs, criterion == "single" & (is.na(q) | abs(q - 0.90) < 1e-9))
feasible <- subset(sel, type1 <= ALPHA_SELECT & power >= POWER_SELECT)
df_frontier <- do.call(rbind, lapply(sort(unique(sel$K)), function(k) {
  s <- feasible[feasible$K == k, ]
  if (!nrow(s)) return(NULL)
  s[which.max(s$power), ]
}))
df_frontier <- df_frontier[order(df_frontier$K), ]
rownames(df_frontier) <- NULL

cat(sprintf("Frontier (selection type I <= %.2f%%, power >= %.1f%%; targets %.1f%%, %.0f%%)\n\n",
            ALPHA_SELECT * 100, POWER_SELECT * 100, ALPHA_TARGET * 100, POWER_TARGET * 100))
disp <- within(df_frontier, {
  type1 <- sprintf("%.2f%%", type1 * 100)
  power <- sprintf("%.2f%%", power * 100)
  EN_H0 <- round(EN_H0); EN_H1 <- round(EN_H1)
})
print(disp[, c("K", "efficacy", "q", "criterion", "type1", "power", "EN_H0", "EN_H1")],
      row.names = FALSE)
cat(sprintf("\nE(N|H1) range over the frontier: %.0f (K=%d) to %.0f (K=%d)\n",
            max(df_frontier$EN_H1), df_frontier$K[which.max(df_frontier$EN_H1)],
            min(df_frontier$EN_H1), df_frontier$K[which.min(df_frontier$EN_H1)]))

save(df_frontier, ALPHA_TARGET, POWER_TARGET, ALPHA_SELECT, POWER_SELECT, file = RDA_PATH)
cat(sprintf("\nWrote %s\n", RDA_PATH))
bayesgsd_save_session("ADRENAL_Calibration_FrontierSelection")
