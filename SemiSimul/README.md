# BayesGSD reproduction code

R code that reproduces all numerical results, tables, and figures in the
manuscript **"Rapid evaluation and calibration of Bayesian group sequential
designs via conjugate-mixture semi-simulation"** (He, Yu, Cro, Billot — submitted
to *Statistics in Medicine*).

The proposed framework is implemented in a single file (`bayseqSim_bern.R`).
Everything else is driver code for the ADRENAL re-design used as the empirical
example throughout the paper.

## Repository layout

This repository is the code component of the project (the `Code/Code v1.0/`
directory of the full manuscript tree). Two usage modes are supported:

- **Standalone (this repository alone).** The core simulator
  (`bayseqSim_bern.R`) and the regression tests (`bayseqSim_bern_tests.R`)
  run directly from the repository root with no further setup. This is
  sufficient to use and verify the framework itself.
- **Full paper reproduction.** The `ADRENAL_*` drivers and
  `ADRENAL_Replication.sh` write caches to `Output/Output v1.0/` and splice
  tables and figures into `Article/Manuscript/Manuscript v1.0/`. They
  therefore require the full project tree (this directory placed at
  `<root>/Code/Code v1.0/` alongside `Article/` and `Output/`), with the
  root supplied via `BAYESGSD_ROOT` (see Configuration below). The manuscript
  sources are not part of this code repository.

## Prerequisites

- R ≥ 4.4 (last verified on R 4.5.2 / `x86_64-apple-darwin`).
- A C++ toolchain that `Rcpp::sourceCpp` can invoke (Xcode CLT on macOS,
  `r-base-dev` on Linux). If the toolchain is missing the framework falls back
  to a pure-R kernel that is approximately one order of magnitude slower; the
  numerical output is identical.
- R packages: `Rcpp`, `parallel`, `RBesT`, `ggplot2`, `ggsci`, `patchwork`.
  The benchmark scripts additionally require `BATSS` (and `INLA`) and
  `adaptr`. INLA is not on CRAN — install per
  <https://www.r-inla.org/download-install>.

## Configuration

The `ADRENAL_*` driver scripts resolve the project root in this order:

1. Environment variable `BAYESGSD_ROOT` (set before launching R).
2. Auto-detect by walking up from `getwd()` looking for the `Code/` and
   `Article/` sibling directories.

If neither succeeds, the driver errors out with a clear message telling you
to set `BAYESGSD_ROOT`. There is **no hard-coded default path** anywhere in
this directory. The regression tests fall back to the current working
directory, so `bayseqSim_bern_tests.R` runs standalone from the repository
root without `BAYESGSD_ROOT`.

Additional overrides:

| Env var                   | Default | Effect                                                                                                                |
|---------------------------|---------|-----------------------------------------------------------------------------------------------------------------------|
| `BAYESGSD_ROOT`  | unset | Required if auto-detection fails.                                                                          |
| `BAYESGSD_CORES` | `8`   | Outer `mclapply` parallelism per script.                                                                   |
| `BAYESGSD_OPENMP`| `0`   | If `1`, build the Rcpp kernel with inner-loop OpenMP enabled (requires a toolchain that supports OpenMP).  |

Quadrature-node count: `options(BayesGSD.quadrature.n = 128L)` (default).
Chunk size for `mclapply`: `options(BayesGSD.chunk_size = 50000L)` (default).
Set either before sourcing `bayseqSim_bern.R`.

## Quick start

Regression tests run standalone from the repository root, with no
environment variable required:

```bash
Rscript bayseqSim_bern_tests.R      # regression tests (~30 s)
```

If the tests pass, the kernel is healthy. Full paper reproduction additionally
requires the project tree (see Repository layout). With `BAYESGSD_ROOT`
pointing at a root that contains `Code/Code v1.0/`, `Article/` and `Output/`,
reproduce all numerical results, tables, figures, and PDFs end-to-end with:

```bash
export BAYESGSD_ROOT=/absolute/path/to/project-root
bash "$BAYESGSD_ROOT/Code/Code v1.0/ADRENAL_Replication.sh"
```

This runs every stage in the order listed below and rebuilds the manuscript
and supplement PDFs at the end. Total wall-clock cost is of order a CPU-day
on 8 cores. Stages 1, 6, and 7 dominate. Per-stage logs are written under
`Output/Output v1.0/logs/`.

## Run order

The scripts are written so that each subsequent stage consumes the cached
outputs of earlier stages. The recommended order is:

1. `ADRENAL_Evaluation_OperatingCharacteristics.R`
   - One-time simulation pass at the union of look times for `K ∈ {2, …, 10}`,
     `R = 1,000,000` virtual trials, under both `H₀` and `H₁`. Writes
     `ADRENAL_OperatingCharacteristics_H0.rda` and `…_H1.rda`
     (each ~ 450 MB) plus eight 12-panel jpegs covering the four
     operating-characteristic metrics (type I, type II, `E(N | H₀)`,
     `E(N | H₁)`) under both `Binding` and `nonBinding` futility. The
     `Binding` jpegs become Figures 1–4 of §3.4 of the main paper; the
     `nonBinding` jpegs become Figures S1–S4 of the supplement. Only the
     type-I-error panels differ numerically between the two binding states
     (see the "Futility convention" bullet under "Notes on the kernel" below).
   - **This is the long-running step (~ tens of CPU-hours).** Everything
     downstream is fast (seconds to minutes) because it reuses these caches.
2. `ADRENAL_Evaluation_MonteCarloError.R`
   - §3.6 simulation-budget study, 100 replicates × R ∈ {1k, 10k, 100k, 1M}.
   - Produces `ADRENAL_MonteCarloError.rda` and the boxplot Figure 5
     (`fig361.jpeg` after Stage 10).
3. `ADRENAL_Diagnostics_PriorApproximation.R`
   - Appendix A.1 prior-approximation diagnostics for the logit-normal example.
   - Produces `ADRENAL_PriorDiagnostics.{rda,jpeg}` (becomes Figure A.1).
4. `ADRENAL_Diagnostics_QuadratureConvergence.R`
   - Appendix A.2 Gauss–Legendre quadrature convergence study at
     `Q ∈ {64, 128, 256, 512}` under both priors.
   - Produces `ADRENAL_QuadratureConvergence.{rda,jpeg}` (becomes Figure A.2).
5. `ADRENAL_Calibration_GridSearch.R`
   - Sweeps the calibrated Haybittle–Peto designs of §3.5 over the cached
     OC. For each of the manuscript-quoted three- and five-look designs,
     reports both the **binding** type I error rate (calibration target)
     and the **non-binding** type I error rate (regulator-facing).
   - Produces `ADRENAL_CalibratedDesign.rda` with `df_const`, `df_hp`,
     `best_const`, `best_hp` (full sweep), and `df_hp_calibrated` (the
     four reported designs with both Type I values side-by-side).
6. `ADRENAL_Benchmark_UniformPrior.R` (§3.3.1)
   - Head-to-head against `BATSS` + `adaptr` under Beta(1,1). **This is the
     other long-running step (~ tens of CPU-hours) because of `BATSS`+INLA.**
   - Produces `ADRENAL_Benchmark_UniformPrior.rda`.
7. `ADRENAL_Benchmark_LogitNormalPrior.R` (§3.3.2)
   - Same benchmark under a logit-normal non-conjugate prior. Same scale of
     wall-clock cost as step 6.
8. `ADRENAL_Benchmark_FullGridTiming.R`
   - End-to-end timing of the 438-design grid sweep (~7 minutes on 8 CPU
     cores, per §3.4 and §4).
9. `ADRENAL_Benchmark_GenerateTables.R`
   - Splices the four benchmark tables (Tables 1–4) into the manuscript .tex.
10. `ADRENAL_Summary.R`
    - Copies/renames the OC jpegs into the manuscript and supplement
      directories under their final names (`fig341.jpeg`–`fig344.jpeg`,
      `fig361.jpeg`, `figA11.jpeg`, `figA12.jpeg`, `figS1.jpeg`–`figS4.jpeg`);
      fills the live `TBD_*` placeholders in the .tex with cached numbers;
      writes `post_run_summary.txt` with the numerical deltas referenced
      throughout §3.3, §3.4, §3.5 and §3.6.
11. Stage 11 of `ADRENAL_Replication.sh` rebuilds both the manuscript and
    the supplement PDFs (`pdflatex` × 3 plus `bibtex` for each document).

Each script writes a `<scriptname>_sessionInfo.rds` and `.txt` alongside its
output, capturing the exact R / package versions used. Together with the
single seed (`SEED = 21L`) this is sufficient for byte-level reproduction.

## Files

| File                                              | Purpose                                                            |
|---------------------------------------------------|--------------------------------------------------------------------|
| `bayseqSim_bern.R`                                   | Core simulator and public API.                                     |
| `bayseqSim_bern_setup.R`                               | Shared setup (project root, parallelism, sessionInfo helper).      |
| `bayseqSim_bern_tests.R`                             | Regression tests; run after any change to `bayseqSim_bern.R`.         |
| `ADRENAL_Evaluation_OperatingCharacteristics.R`   | §3.4 OC simulation (long-running).                                 |
| `ADRENAL_Evaluation_MonteCarloError.R`            | §3.6 MC error study.                                               |
| `ADRENAL_Diagnostics_PriorApproximation.R`        | Appendix A.1 prior-approximation diagnostics.                      |
| `ADRENAL_Diagnostics_QuadratureConvergence.R`     | Appendix A.2 Gauss–Legendre quadrature convergence study.          |
| `ADRENAL_Calibration_GridSearch.R`                | §3.5 calibrated three- and five-look designs.                      |
| `ADRENAL_Benchmark_UniformPrior.R`                | §3.3.1 head-to-head benchmark (long-running).                      |
| `ADRENAL_Benchmark_LogitNormalPrior.R`            | §3.3.2 head-to-head benchmark (long-running).                      |
| `ADRENAL_Benchmark_FullGridTiming.R`              | §3.4 / §4 full 438-design grid timing.                             |
| `ADRENAL_Benchmark_GenerateTables.R`              | Generates Tables 1–4 and splices them into the manuscript .tex.    |
| `ADRENAL_Summary.R`                               | Figure copy/rename + `TBD_*` placeholder splice + summary report.  |

## Notes on the kernel

- **Mixture-prior posterior weights.** As of the May 2026 fix, the posterior
  mixture weights are correctly renormalised by the per-component marginal
  likelihood `B(a + e, b + m − e) / B(a, b)` (computed in log space via
  `lbeta` / `lgamma`). The previous "legacy convention" of using prior
  weights directly is no longer applied. Single-component priors (the
  Beta(1,1) baseline) are unaffected.
- **Futility convention.** The framework supports both binding and
  non-binding futility, per §3.2 of the main paper. Under binding futility
  all stopping decisions enforce the futility rule. Under non-binding
  futility the type I error is computed *without* the futility stop (the
  rule is advisory for the false-positive accounting); the expected sample
  size, the early-stopping probability and the type II error are reported
  from the futility-active sweep on the assumption that the data-monitoring
  committee follows the rule in practice. As a consequence, only the
  type-I-error panel of the §3.4 figures differs between binding and
  non-binding (Figure 1 of the main paper vs Figure S1 of the supplement);
  Figures S2–S4 are graphically identical to Figures 2–4.
- **Gauss–Legendre quadrature.** Fixed 128 nodes on (0, 1). For sensitivity
  see the convergence study in Appendix A.2.
- **Parallel RNG.** The kernel installs `RNGkind("L'Ecuyer-CMRG")` at load
  time so that any future `mclapply`-based kernel that draws inside a chunk
  is reproducible.
- **Boundary handling.** For relative-risk and odds-ratio transforms the
  kernel zeroes out quadrature nodes whose transformed value falls outside
  (0, 1). This is correct: those `(θ_c, θ_t)` pairs are outside the support
  of the joint posterior on the unit square, not "lost probability mass."
- **Regression tests.** `bayseqSim_bern_tests.R` covers six checks: (1)
  two-component mixture posterior weight renormalisation, (2) Beta(1,1)
  delta-tail convolution against direct integration, (3) numeric-proximity
  threshold lookup, (4) `options(warn)` restored on exit, (5) the
  L'Ecuyer-CMRG RNG kind, and (6) the elbow rule for mixture-component
  selection. All six should pass in under 30 seconds.
