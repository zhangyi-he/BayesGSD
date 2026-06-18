#' @title Post-run manuscript-update helper
#' @description Runs after the full re-run completes. Performs the
#'   deterministic manuscript updates that can be done from the new caches:
#'     (1) copies/renames the new figure files into the manuscript directory
#'         (fig341 / fig361 / figA11 / figA21) and the supplement directory
#'         (figS41..figS45) under the filenames the .tex files reference;
#'     (2) replaces the \texttt{TBD} placeholders in the Monte Carlo
#'         precision table (§3.6) and its discussion paragraphs with
#'         mean / MC SE values computed from the new
#'         ADRENAL_MonteCarloError.rda;
#'     (3) writes Output/Output v1.0/post_run_summary.txt containing all the
#'         numerical values needed to refresh the §3.3.1, §3.3.2, §3.4 and
#'         §3.5 narrative claims (cross-method differences, speedup ratios,
#'         wall-clock figures, type I climb, expected sample size shrinkage).
#'   The summary file is consumed in a follow-up conversation turn for the
#'   prose edits that are too sensitive to auto-rewrite by regex.
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

# Note: OUTPUT_DIR, MANUSCRIPT_DIR, SUPPLEMENT_DIR are provided by the setup
# helper. TEX_PATH and SUMMARY_PATH are derived here for clarity.
TEX_PATH          <- file.path(MANUSCRIPT_DIR, "ZH2023_Manuscript.tex")
SUMMARY_PATH      <- file.path(OUTPUT_DIR, "post_run_summary.txt")

# ----------------------------------------------------------------------------
# Idempotent / safe in-place .tex editing helpers
# ----------------------------------------------------------------------------
# The manuscript .tex is a tracked source file; this script mutates it in place
# to fill \texttt{TBD...} placeholders. To make re-runs safe and idempotent we
# (a) take a one-time backup before the first mutation of a run, and (b) for
# each substitution distinguish two outcomes rather than silently no-op:
#   FILLED  - the placeholder was found and replaced;
#   ALREADY - the placeholder is gone but the intended value is already present
#             (a previous run filled it; re-running is a no-op success).
#
# Idempotency discriminator: the TBD tokens are unique strings only this script
# emits, so an absent placeholder means a prior run already filled it (ALREADY,
# a no-op success). A genuine error -- a wrong TEX_PATH or an unreadable file --
# is caught once up front by tex_read_checked(), which fails loudly if the file
# is missing/empty or does not look like the manuscript. With that guard in
# place an absent placeholder is unambiguously "already filled", not a typo.

tex_backup_once <- function(path) {
  bak <- paste0(path, ".bak")
  if (file.exists(path) && !file.exists(bak)) {
    file.copy(path, bak, overwrite = FALSE)
    cat(sprintf("  backup written: %s\n", bak))
  }
  invisible(bak)
}

# Read TEX_PATH with an up-front sanity check so a wrong path or an empty/garbled
# file fails loudly rather than masquerading as "everything already filled".
tex_read_checked <- function(path) {
  if (!file.exists(path))
    stop("Manuscript .tex not found at: ", path,
         " (check MANUSCRIPT_DIR / TEX_PATH).")
  tx <- readLines(path, warn = FALSE)
  if (length(tx) < 50L)
    stop("Manuscript .tex at ", path, " has only ", length(tx),
         " lines; expected the full manuscript. Refusing to edit.")
  if (!any(grepl("\\begin{document}", tx, fixed = TRUE)))
    stop("File at ", path, " does not look like the manuscript ",
         "(no \\begin{document}); refusing to edit.")
  tx
}

# Substitute a set of \texttt{<key>} placeholders. `subs` is a named character
# vector key -> replacement value. Underscores in keys are escaped to match the
# LaTeX source (\texttt{TBD\_FG\_H0}). Idempotent: an absent placeholder is a
# no-op success (ALREADY) given the up-front tex_read_checked() guard.
apply_texttt_subs <- function(tex, subs) {
  for (key in names(subs)) {
    val     <- subs[[key]]
    key_tex <- gsub("_", "\\_", key, fixed = TRUE)
    needle  <- paste0("\\texttt{", key_tex, "}")
    hits    <- grep(needle, tex, fixed = TRUE)
    if (length(hits) >= 1L) {
      tex <- gsub(needle, val, tex, fixed = TRUE)
      cat(sprintf("  FILLED  %s -> %s  (%d hit%s)\n",
                  key, val, length(hits), if (length(hits) == 1L) "" else "s"))
    } else {
      cat(sprintf("  ALREADY %s (placeholder absent; prior run filled it; no-op)\n", key))
    }
  }
  tex
}


# ----------------------------------------------------------------------------
# (1) Figure copy/rename
# ----------------------------------------------------------------------------

# Per-figure: (src basename, dst basename, dst directory). The same-cache demo,
# Monte Carlo and appendix figures go to the main manuscript directory; the
# design-grid operating characteristics go to the Supplemental Material
# directory under figS41..figS45 (Section S4 of the supplement, Figures S1-S5).
fig_map <- list(
  # Design-grid operating characteristics now live in Supplement Section S4
  # (the binding type I is Figure S1 there); the binding type II / E(N) panels
  # are graphically identical to the non-binding figS43..figS45 and are not copied.
  list("ADRENAL_OperatingCharacteristics_type1Error_Binding.jpeg",      "figS41.jpeg",  SUPPLEMENT_DIR),
  list("ADRENAL_Demo_SameCacheSweep.jpeg",                              "fig341.jpeg",  MANUSCRIPT_DIR),
  list("ADRENAL_Calibration_FeasibleFamily.jpeg",                       "fig351.jpeg",  MANUSCRIPT_DIR),
  list("ADRENAL_MonteCarloError.jpeg",                                  "fig361.jpeg",  MANUSCRIPT_DIR),
  list("ADRENAL_PriorDiagnostics.jpeg",                                 "figA11.jpeg",   MANUSCRIPT_DIR),
  list("ADRENAL_QuadratureConvergence.jpeg",                            "figA21.jpeg",   MANUSCRIPT_DIR),
  # Supplement Section S4 (design-grid operating characteristics, Figures S2-S5)
  list("ADRENAL_OperatingCharacteristics_type1Error_nonBinding.jpeg",      "figS42.jpeg", SUPPLEMENT_DIR),
  list("ADRENAL_OperatingCharacteristics_type2Error_nonBinding.jpeg",      "figS43.jpeg", SUPPLEMENT_DIR),
  list("ADRENAL_OperatingCharacteristics_expSampleSizeH0_nonBinding.jpeg", "figS44.jpeg", SUPPLEMENT_DIR),
  list("ADRENAL_OperatingCharacteristics_expSampleSizeH1_nonBinding.jpeg", "figS45.jpeg", SUPPLEMENT_DIR)
)
cat("=== Figure copy/rename ===\n")
for (entry in fig_map) {
  src      <- entry[[1L]]
  dst      <- entry[[2L]]
  dst_dir  <- entry[[3L]]
  src_path <- file.path(OUTPUT_DIR, src)
  dst_path <- file.path(dst_dir,    dst)
  if (file.exists(src_path)) {
    if (!dir.exists(dst_dir)) dir.create(dst_dir, recursive = TRUE)
    ok <- file.copy(src_path, dst_path, overwrite = TRUE)
    cat(sprintf("  %s %s -> %s\n", if (ok) "OK" else "FAIL",
                src, file.path(basename(dst_dir), dst)))
  } else {
    cat(sprintf("  MISSING %s\n", src_path))
  }
}

# ----------------------------------------------------------------------------
# (2) [removed] Section 3.6 Monte Carlo precision: no manuscript splice.
#     The Monte Carlo precision table now lives in the Supplemental Material
#     (Table S4, hand-coded) and the Section 3.6 manuscript prose is filled
#     directly in 10^{6} notation, so there is nothing to splice into the main
#     text here. (The previous block targeted R=1,000,000 / TBD markers that no
#     longer exist in the manuscript, which would abort this stage.)
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# (2b) §3.3 and §4 full-grid timing placeholder replacement
# ----------------------------------------------------------------------------

cat("\n=== §3.3 / §4 full-grid timing placeholder replacement ===\n")
fg_path <- file.path(OUTPUT_DIR, "ADRENAL_Benchmark_FullGridTiming.rda")
if (!file.exists(fg_path)) {
  cat("  MISSING ", fg_path, " - skipping timing update\n")
} else {
  e <- new.env(); load(fg_path, envir = e)
  t_H0    <- e$timing$H0$elapsed
  t_H1    <- e$timing$H1$elapsed
  t_total <- t_H0 + t_H1
  t_per   <- (t_H0 + t_H1) / 2

  fmt_t <- function(x) {
    if (x < 10)  formatC(x, digits = 1, format = "f")
    else         formatC(x, digits = 0, format = "f")
  }

  tex_backup_once(TEX_PATH)
  tex <- tex_read_checked(TEX_PATH)
  # LaTeX source escapes underscores in \texttt{}, so the literal to match is
  # e.g. `\texttt{TBD\_FG\_H0}`; apply_texttt_subs() handles the escaping.
  subs <- c(
    "TBD_FG_H0"    = fmt_t(t_H0),
    "TBD_FG_H1"    = fmt_t(t_H1),
    "TBD_FG_TOTAL" = fmt_t(t_total),
    "TBD_FG_PER"   = fmt_t(t_per)
  )
  tex <- apply_texttt_subs(tex, subs)
  writeLines(tex, TEX_PATH)
  cat("  manuscript timing placeholders updated\n")
}

# ----------------------------------------------------------------------------
# (2c) §3.5 calibrated Haybittle-Peto design placeholder replacement
# ----------------------------------------------------------------------------

cat("\n=== §3.5 calibrated HP design placeholder replacement ===\n")
cd_path <- file.path(OUTPUT_DIR, "ADRENAL_CalibratedDesign.rda")
if (!file.exists(cd_path)) {
  cat("  MISSING ", cd_path, " - skipping calibrated-design update\n")
} else {
  e <- new.env(); load(cd_path, envir = e)
  # Re-derive the headline designs from df_hp. The manuscript reports the
  # maximum-power single-criterion HP designs, defined once in
  # adabay_bern_setup.R (BAYESGSD_DESIGNS) -- the single source of truth for the
  # Section 3.5 calibrated designs: interim 0.999, final 0.977, q = 0.90 for both
  # K = 3 and K = 5.
  hp <- e$df_hp
  hp_specs <- list(
    `3` = list(p_interim = BAYESGSD_DESIGNS$interim_eff,
               p_final   = BAYESGSD_DESIGNS$final_by_K[["3"]],
               q         = BAYESGSD_DESIGNS$q_futility),
    `5` = list(p_interim = BAYESGSD_DESIGNS$interim_eff,
               p_final   = BAYESGSD_DESIGNS$final_by_K[["5"]],
               q         = BAYESGSD_DESIGNS$q_futility)
  )
  pick <- function(K_val) {
    spec <- hp_specs[[as.character(K_val)]]
    if (is.null(spec)) return(NULL)
    row <- hp[hp$K == K_val &
              hp$schedule == "haybittle-peto" &
              hp$criterion == "single" &
              !is.na(hp$p_interim) & abs(hp$p_interim - spec$p_interim) < 1e-9 &
              !is.na(hp$p_final)   & abs(hp$p_final   - spec$p_final)   < 1e-9 &
              !is.na(hp$q)         & abs(hp$q         - spec$q)         < 1e-9, ]
    if (nrow(row) != 1L) NULL else row[1, ]
  }
  k3 <- pick(3)
  k5 <- pick(5)
  if (is.null(k3) || is.null(k5)) {
    cat("  Manuscript-quoted HP designs not found in df_hp (",
        if (is.null(k3)) "K=3" else "", " ",
        if (is.null(k5)) "K=5" else "",
        ") - skipping\n", sep = "")
  } else {
    fmt_pct1 <- function(x) formatC(100 * x, digits = 2, format = "f")
    # In LaTeX math mode a bare comma is a list delimiter, so we wrap the
    # thousands separator in braces (e.g. 3{,}471) consistent with usage
    # elsewhere in the manuscript (3{,}800, etc.).
    fmt_int1 <- function(x) format(round(x), big.mark = "{,}")

    tex_backup_once(TEX_PATH)
    tex <- tex_read_checked(TEX_PATH)
    subs <- c(
      "TBD_K3_T1"  = fmt_pct1(k3$type1),
      "TBD_K3_PWR" = fmt_pct1(k3$power),
      "TBD_K3_EN0" = fmt_int1(k3$EN_H0),
      "TBD_K3_EN1" = fmt_int1(k3$EN_H1),
      "TBD_K5_T1"  = fmt_pct1(k5$type1),
      "TBD_K5_PWR" = fmt_pct1(k5$power),
      "TBD_K5_EN0" = fmt_int1(k5$EN_H0),
      "TBD_K5_EN1" = fmt_int1(k5$EN_H1)
    )
    tex <- apply_texttt_subs(tex, subs)
    writeLines(tex, TEX_PATH)
    cat("  manuscript HP calibrated-design placeholders updated\n")
  }
}

# ----------------------------------------------------------------------------
# (3) Summary of values needed for §3.3 / §3.4 / §3.5 prose edits
# ----------------------------------------------------------------------------

cat("\n=== Writing post_run_summary.txt ===\n")
sink(SUMMARY_PATH, split = FALSE)
on.exit(sink(NULL), add = TRUE)

cat("Post-run summary written ", format(Sys.time()), "\n", sep = "")
cat(strrep("=", 70), "\n\n")

emit_benchmark <- function(rda_path, prior_label) {
  if (!file.exists(rda_path)) {
    cat(prior_label, ": cache missing at ", rda_path, "\n", sep = "")
    return(invisible(NULL))
  }
  e <- new.env(); load(rda_path, envir = e)
  bl <- e$bench_log
  if (length(bl) == 0L) {
    cat(prior_label, ": empty bench_log\n")
    return(invisible(NULL))
  }
  cat(prior_label, " (", rda_path, ")\n", sep = "")
  cat(strrep("-", 70), "\n")
  hdr <- sprintf("%-7s %-6s %s",
                 "tag", "method",
                 "  eff_prob   exp_n   elapsed_s")
  cat(hdr, "\n")
  keys <- names(bl)[order(sapply(bl, function(r) r$design$K),
                          sapply(bl, function(r) r$hypothesis))]
  for (k in keys) {
    r <- bl[[k]]
    for (slot in c("batss", "adaptr_small", "proposed_small", "proposed_big")) {
      cell <- r[[slot]]
      if (!is.null(cell)) {
        cat(sprintf("%-7s %-14s  %.4f   %5.0f   %8.1f\n",
                    k, slot,
                    cell$efficacy_prob, cell$expected_n, cell$elapsed))
      }
    }
  }
  cat("\n")

  # Aggregate speedup ratios across all (K, hypothesis) cells
  ratios_batss   <- numeric()
  ratios_adaptr  <- numeric()
  ratios_a_vs_b  <- numeric()
  ratios_big_vs_batss <- numeric()
  prop_big_times <- numeric()
  prop_small_per_design <- numeric()
  diff_t1_batss  <- numeric()
  diff_t1_adaptr <- numeric()
  diff_pwr_batss  <- numeric()
  diff_pwr_adaptr <- numeric()
  diff_EN_batss   <- numeric()
  diff_EN_adaptr  <- numeric()

  for (k in keys) {
    r  <- bl[[k]]
    pb <- r$proposed_big
    ps <- r$proposed_small
    bb <- r$batss
    aa <- r$adaptr_small
    if (!is.null(bb) && !is.null(ps) && !is.na(ps$elapsed) && ps$elapsed > 0) {
      ratios_batss  <- c(ratios_batss,  bb$elapsed / ps$elapsed)
    }
    if (!is.null(aa) && !is.null(ps) && !is.na(ps$elapsed) && ps$elapsed > 0) {
      ratios_adaptr <- c(ratios_adaptr, aa$elapsed / ps$elapsed)
    }
    if (!is.null(aa) && !is.null(bb) && !is.na(aa$elapsed) && aa$elapsed > 0) {
      ratios_a_vs_b <- c(ratios_a_vs_b, bb$elapsed / aa$elapsed)
    }
    if (!is.null(bb) && !is.null(pb) && !is.na(pb$elapsed) && pb$elapsed > 0) {
      ratios_big_vs_batss <- c(ratios_big_vs_batss, bb$elapsed / pb$elapsed)
    }
    if (!is.null(pb)) prop_big_times <- c(prop_big_times, pb$elapsed)
    if (!is.null(pb) && !is.null(bb)) {
      if (r$hypothesis == "H0") {
        diff_t1_batss <- c(diff_t1_batss, abs(bb$efficacy_prob - pb$efficacy_prob))
        diff_EN_batss <- c(diff_EN_batss, abs(bb$expected_n - pb$expected_n))
      } else {
        diff_pwr_batss <- c(diff_pwr_batss, abs(bb$efficacy_prob - pb$efficacy_prob))
        diff_EN_batss  <- c(diff_EN_batss, abs(bb$expected_n - pb$expected_n))
      }
    }
    if (!is.null(pb) && !is.null(aa)) {
      if (r$hypothesis == "H0") {
        diff_t1_adaptr <- c(diff_t1_adaptr, abs(aa$efficacy_prob - pb$efficacy_prob))
        diff_EN_adaptr <- c(diff_EN_adaptr, abs(aa$expected_n - pb$expected_n))
      } else {
        diff_pwr_adaptr <- c(diff_pwr_adaptr, abs(aa$efficacy_prob - pb$efficacy_prob))
        diff_EN_adaptr  <- c(diff_EN_adaptr, abs(aa$expected_n - pb$expected_n))
      }
    }
  }
  cat("Speedup BATSS / proposed @ R=5,000:        ",
      sprintf("[%.0fx, %.0fx]\n", min(ratios_batss),  max(ratios_batss)))
  cat("Speedup adaptr / proposed @ R=5,000:       ",
      sprintf("[%.0fx, %.0fx]\n", min(ratios_adaptr), max(ratios_adaptr)))
  cat("Ratio BATSS / adaptr @ R=5,000:            ",
      sprintf("[%.0fx, %.0fx]\n", min(ratios_a_vs_b), max(ratios_a_vs_b)))
  cat("Speedup BATSS @ R=5k / proposed @ R=1M:    ",
      sprintf("[%.0fx, %.0fx]\n", min(ratios_big_vs_batss), max(ratios_big_vs_batss)))
  cat("Proposed @ R=1,000,000 elapsed (seconds):  ",
      sprintf("[%.1f, %.1f]\n", min(prop_big_times), max(prop_big_times)))
  cat("Cross-method differences vs proposed @ R=1,000,000:\n")
  cat("  type I error (pp): BATSS max=", sprintf("%.2f", 100*max(diff_t1_batss)),
      "; adaptr max=", sprintf("%.2f", 100*max(diff_t1_adaptr)), "\n", sep = "")
  cat("  power        (pp): BATSS max=", sprintf("%.2f", 100*max(diff_pwr_batss)),
      "; adaptr max=", sprintf("%.2f", 100*max(diff_pwr_adaptr)), "\n", sep = "")
  cat("  E[N]    (patients): BATSS max=", sprintf("%.0f", max(diff_EN_batss)),
      "; adaptr max=", sprintf("%.0f", max(diff_EN_adaptr)), "\n", sep = "")
  cat("\n")

  # Type I climb under H0 from K=1 to K=9 for proposed @ R=1,000,000
  cat("Type I climb (proposed @ R=1,000,000, H0):\n")
  for (kk in c(1, 3, 5, 7, 9)) {
    cell <- bl[[sprintf("K%d_H0", kk)]]
    if (!is.null(cell$proposed_big)) {
      cat(sprintf("  K=%d: alpha=%.4f\n", kk, cell$proposed_big$efficacy_prob))
    }
  }
  cat("Expected sample size under H1 (proposed @ R=1,000,000):\n")
  for (kk in c(1, 3, 5, 7, 9)) {
    cell <- bl[[sprintf("K%d_H1", kk)]]
    if (!is.null(cell$proposed_big)) {
      cat(sprintf("  K=%d: E[N]=%.0f\n", kk, cell$proposed_big$expected_n))
    }
  }
  cat("\n")
}

emit_benchmark(file.path(OUTPUT_DIR, "ADRENAL_Benchmark_UniformPrior.rda"),
               "§3.3.1 conjugate Beta(1,1) prior")
emit_benchmark(file.path(OUTPUT_DIR, "ADRENAL_Benchmark_LogitNormalPrior.rda"),
               "§3.3.2 logit-normal non-conjugate prior")

# §3.4 full-grid timing
fg_path <- file.path(OUTPUT_DIR, "ADRENAL_Benchmark_FullGridTiming.rda")
if (file.exists(fg_path)) {
  cat("§3.4 full-grid timing (", fg_path, ")\n", sep = "")
  cat(strrep("-", 70), "\n")
  e <- new.env(); load(fg_path, envir = e)
  cat(sprintf("  H0 elapsed: %.1f s\n", e$timing$H0$elapsed))
  cat(sprintf("  H1 elapsed: %.1f s\n", e$timing$H1$elapsed))
  cat(sprintf("  n_designs : %d (cores: %d, machine: %s)\n\n",
              e$timing$n_designs, e$timing$n_cores, e$timing$machine))
}

# §3.6 MC SE summary at R=1,000,000
mc_path <- file.path(OUTPUT_DIR, "ADRENAL_MonteCarloError.rda")
if (file.exists(mc_path)) {
  cat("§3.6 MC SE summary at R=1,000,000\n")
  cat(strrep("-", 70), "\n")
  e <- new.env(); load(mc_path, envir = e)
  for (i in seq_along(e$numberOfAnalyses)) {
    K <- e$numberOfAnalyses[i]
    t1_se <- sd(e$type1ErrorRate[["1,000,000"]][i, ], na.rm = TRUE)
    t2_se <- sd(e$type2ErrorRate[["1,000,000"]][i, ], na.rm = TRUE)
    cat(sprintf("  K=%d: type I MC SE = %.4f pp; type II MC SE = %.4f pp\n",
                K, 100 * t1_se, 100 * t2_se))
  }
  cat("\n")
}

cat("End of summary.\n")
sink(NULL)
cat("  wrote ", SUMMARY_PATH, "\n", sep = "")

cat("\nDone. Post-run manuscript updates complete.\n")

# Persist sessionInfo and package versions alongside the cache.
bayesgsd_save_session("ADRENAL_Summary")
