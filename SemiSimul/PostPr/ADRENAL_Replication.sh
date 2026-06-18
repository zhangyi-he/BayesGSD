#!/usr/bin/env bash
#
# ADRENAL_Replication.sh -- end-to-end reproduction of all numerical results, tables,
# and figures in the manuscript and supplement.
#
# Usage:
#
#   export BAYESGSD_ROOT=/absolute/path/to/repo
#   bash "${BAYESGSD_ROOT}/Code/Code v1.0/ADRENAL_Replication.sh"
#
# or, from anywhere inside the repository tree (auto-detects BAYESGSD_ROOT):
#
#   bash "Code/Code v1.0/ADRENAL_Replication.sh"
#
# Optional environment variables:
#
#   BAYESGSD_ROOT   absolute path to the project root. Auto-detected from this
#                   script's location if unset.
#   BAYESGSD_CORES  outer mclapply parallelism per stage (default 8).
#   BAYESGSD_OPENMP if "1", build the Rcpp kernel with inner-loop OpenMP
#                   (requires a toolchain that supports it).
#
# Stage order matches Code/Code v1.0/README.md. Each stage exits non-zero on
# failure; the script aborts on the first error so you can inspect the log
# before the next (potentially long) stage starts.
#
# Wall-clock cost on an 8-core commodity machine (R 4.5.2, x86_64-apple-darwin):
#
#   Stage  0 -- tests                                        ~ 30 s
#   Stage  1 -- OC simulation (LONG)                         ~ tens of CPU-hours
#   Stage  2 -- MC error study                               ~ minutes
#   Stage  3 -- Appendix A.1 prior approximation             ~ seconds
#   Stage  3b -- Appendix A.1 OC sensitivity to L            ~ 1-2 minutes
#   Stage  4 -- Appendix A.2 quadrature convergence          ~ minutes
#   Stage  5 -- Section 3.5 calibration grid search          ~ minutes
#   Stage  5b -- Section 3.5 seed-sensitivity sweep          ~ 2 minutes
#   Stage  5c -- Section 3.5 plausibility-set sweep          ~ 4 minutes
#   Stage  5d -- Section 3.5 maximum-power frontier table         < 1 minute
#   Stage  6 -- Section 3.3.1 benchmark vs BATSS+adaptr      ~ tens of CPU-hours
#   Stage  6b -- Section 3.3.1 K=9 matched-budget seed sweep ~ 3 minutes
#   Stage  7 -- Section 3.3.2 benchmark (logit-normal)       ~ tens of CPU-hours
#   Stage  8 -- Section 3.4 / Section 4 full-grid timing     ~ seconds
#   Stage  8b -- Section 2.4 same-cache demo figure          ~ minutes
#   Stage  9 -- Tables 1-4 generation                        ~ seconds
#   Stage  9b -- Supplement Table S2 (per-look stop probs)   ~ seconds
#   Stage 10 -- Figure copy/rename + summary                 ~ seconds
#   Stage 11 -- Rebuild manuscript + supplement PDFs         ~ minute
#
# Total: of order a CPU-day on 8 cores. Stages 1, 6, 7 dominate.

set -euo pipefail

# --- Resolve project root ---------------------------------------------------
if [ -z "${BAYESGSD_ROOT:-}" ]; then
  HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  CAND="$(cd "$HERE/../.." && pwd)"
  if [ -d "$CAND/Code" ] && [ -d "$CAND/Article" ]; then
    export BAYESGSD_ROOT="$CAND"
  else
    echo "ERROR: BAYESGSD_ROOT is not set and could not be auto-detected." >&2
    echo "Set it explicitly:" >&2
    echo "  export BAYESGSD_ROOT=/absolute/path/to/repo" >&2
    exit 1
  fi
fi

cd "$BAYESGSD_ROOT"
CODE_DIR="Code/Code v1.0"
OUTPUT_DIR="Output/Output v1.0"
MS_DIR="Article/Manuscript/Manuscript v1.0"
SUPP_DIR="$MS_DIR/Supplemental Material"
LOG_DIR="$OUTPUT_DIR/logs"
mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

echo "=== BayesGSD reproduction ==="
echo "    PROJECT_ROOT = $BAYESGSD_ROOT"
echo "    OUTPUT_DIR   = $OUTPUT_DIR"
echo "    BAYESGSD_CORES = ${BAYESGSD_CORES:-8}"
echo "    BAYESGSD_OPENMP = ${BAYESGSD_OPENMP:-0}"
echo

# --- Stage runner -----------------------------------------------------------
run_stage() {
  local script="$1"
  local label="$2"
  local ts; ts="$(date +%Y%m%d-%H%M%S)"
  local log="$LOG_DIR/${label}_${ts}.log"
  echo "=== [$ts] $label  ($script)"
  echo "    log: $log"
  if Rscript "$CODE_DIR/$script" >"$log" 2>&1; then
    echo "    OK   $(tail -n 1 "$log")"
  else
    echo "    FAIL (see $log)"
    exit 1
  fi
}

# --- 0. Sanity-check the kernel before any heavy stage ---------------------
run_stage "adabay_bern_tests.R"                          "00_tests"

# --- 1. Operating-characteristics simulation (LONG) ------------------------
#     Section 3.4 reference OC at R = 1,000,000 under H0 and H1.
run_stage "ADRENAL_Evaluation_OperatingCharacteristics.R" "01_oc"

# --- 2. Section 3.6 Monte Carlo error study --------------------------------
run_stage "ADRENAL_Evaluation_MonteCarloError.R"          "02_mc"

# --- 3. Appendix A.1 prior approximation diagnostics -----------------------
run_stage "ADRENAL_Diagnostics_PriorApproximation.R"      "03_priordiag"

# --- 3b. Appendix A.1 OC sensitivity to mixture order L --------------------
run_stage "ADRENAL_Diagnostics_OCSensitivity.R"           "03b_ocsens"

# --- 4. Appendix A.2 Gauss-Legendre quadrature convergence -----------------
run_stage "ADRENAL_Diagnostics_QuadratureConvergence.R"   "04_quadconv"

# --- 5. Section 3.5 calibration grid search --------------------------------
run_stage "ADRENAL_Calibration_GridSearch.R"              "05_calib"

# --- 5b. Section 3.5 calibration-design seed sensitivity -------------------
run_stage "ADRENAL_Calibration_SeedSensitivity.R"         "05b_calib_seed"

# --- 5c. Section 3.5 plausibility-set sweep --------------------------------
run_stage "ADRENAL_Calibration_PlausibilitySweep.R"       "05c_calib_plaus"
run_stage "ADRENAL_Calibration_FrontierSelection.R"       "05d_calib_frontier"

# --- 6. Section 3.3.1 benchmark vs BATSS + adaptr (LONG) -------------------
run_stage "ADRENAL_Benchmark_UniformPrior.R"              "06_bench_unif"

# --- 6b. Section 3.3.1 K=9 matched-budget seed sensitivity -----------------
run_stage "ADRENAL_Benchmark_SeedSensitivity.R"        "06b_bench_seed"

# --- 7. Section 3.3.2 benchmark under logit-normal prior (LONG) ------------
run_stage "ADRENAL_Benchmark_LogitNormalPrior.R"          "07_bench_logn"

# --- 8. Section 3.4 / Section 4 full-grid timing ---------------------------
run_stage "ADRENAL_Benchmark_FullGridTiming.R"            "08_fg_timing"

# --- 8b. Section 2.4 same-cache demonstration figure -----------------------
run_stage "ADRENAL_Demo_SameCacheSweep.R"                 "08b_samecache"

# --- 9. Generate Tables 1-4 and splice into the manuscript .tex ------------
run_stage "ADRENAL_Benchmark_GenerateTables.R"            "09_tables"

# --- 9b. Supplement Table S2 (per-look stopping probabilities, K = 5) ------
run_stage "ADRENAL_Evaluation_StoppingProbabilities.R"    "09b_stopprob"

# --- 10. Copy/rename figures + TBD_* placeholder splice + summary ----------
run_stage "ADRENAL_Summary.R"                             "10_summary"

# --- 11. Rebuild the manuscript and supplement PDFs ------------------------
echo
echo "=== [11_pdfs] Rebuilding manuscript and supplement PDFs"

# Capture pdflatex/bibtex output in per-document log files so failures can be
# inspected; matches the per-stage logging pattern used for stages 0-10.
PDF_LOG_MAIN="$LOG_DIR/11_pdfs_manuscript.log"
PDF_LOG_SUPP="$LOG_DIR/11_pdfs_supplement.log"

(
  cd "$MS_DIR"
  pdflatex -interaction=nonstopmode ZH2023_Manuscript.tex
  bibtex   ZH2023_Manuscript
  pdflatex -interaction=nonstopmode ZH2023_Manuscript.tex
  pdflatex -interaction=nonstopmode ZH2023_Manuscript.tex
) > "$PDF_LOG_MAIN" 2>&1
echo "    OK   $MS_DIR/ZH2023_Manuscript.pdf  (log: $PDF_LOG_MAIN)"

(
  cd "$SUPP_DIR"
  pdflatex -interaction=nonstopmode ZH2023_Supplemental_Material.tex
  bibtex   ZH2023_Supplemental_Material
  pdflatex -interaction=nonstopmode ZH2023_Supplemental_Material.tex
  pdflatex -interaction=nonstopmode ZH2023_Supplemental_Material.tex
) > "$PDF_LOG_SUPP" 2>&1
echo "    OK   $SUPP_DIR/ZH2023_Supplemental_Material.pdf  (log: $PDF_LOG_SUPP)"

echo
echo "=== All stages completed successfully ==="
echo "    Manuscript: $MS_DIR/ZH2023_Manuscript.pdf"
echo "    Supplement: $SUPP_DIR/ZH2023_Supplemental_Material.pdf"
