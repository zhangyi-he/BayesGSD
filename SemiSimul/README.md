# BayesGSD reproduction code

R code that reproduces all numerical results, tables, and figures in the
manuscript **"Rapid evaluation and calibration of Bayesian group sequential
designs via conjugate-mixture semi-simulation"** (He, Yu, Cro, Billot — submitted
to *Statistics in Medicine*).

The proposed framework is implemented in a single file (`adabay_bern.R`).
Everything else is driver code for the ADRENAL re-design used as the empirical
example throughout the paper.

## Contents

- **Framework.** `adabay_bern.R` (the simulator and public API) and
  `adabay_bern_tests.R` (its regression tests) run standalone from the
  repository root with no configuration — enough to use and verify the
  framework itself (see Quick start).
- **Drivers.** The `ADRENAL_*` scripts reproduce the analyses, tables, and
  figures for the ADRENAL re-design and are documented under Run order and
  Files below. Full end-to-end reproduction — running the drivers, caching
  their outputs, and splicing the tables and figures into the manuscript — is
  carried out within the complete project tree, which also holds the manuscript
  sources and the output caches and is not part of this code repository.

## Prerequisites

- R ≥ 4.4 (last verified on R 4.5.2 / `x86_64-apple-darwin`).
- A C++ toolchain that `Rcpp::sourceCpp` can invoke (Xcode CLT on macOS,
  `r-base-dev` on Linux). If the toolchain is missing the framework falls back
  to a pure-R kernel that is approximately one order of magnitude slower; the
  numerical output is identical.
- R packages: `Rcpp`, `parallel`, `RBesT`, `ggplot2`, `ggsci`, `ggrepel`,
  `patchwork`, `scales`. The benchmark scripts additionally require `BATSS`
  (and `INLA`) and `adaptr`. INLA is not on CRAN — install per
  <https://www.r-inla.org/download-install>.
- LaTeX (`pdflatex` and `bibtex`) is needed only for the final manuscript and
  supplement build; the R results and deposited caches do not depend on it.

### Pinned versions (`renv.lock`)

The exact package versions used to produce every result in the manuscript
are pinned in `renv.lock` (assembled from the per-script `*_sessionInfo.txt`
files that each driver writes). To install them into a project-local library
without touching your system library:

```r
install.packages("renv")
renv::restore()       # run from the repository root
```

See [`DEPENDENCIES.md`](DEPENDENCIES.md) for provenance, INLA repository
notes, and caveats.

## Quick start

The framework and its regression tests run standalone from the repository
root, with no configuration:

```bash
Rscript adabay_bern_tests.R      # regression tests (~30 s)
```

If the tests pass, the kernel is healthy and `adabay_bern.R` can be sourced
directly to use the public API. The `ADRENAL_*` drivers that reproduce the
paper's results are documented below.

## Run order

The scripts are written so that each subsequent stage consumes the cached
outputs of earlier stages. The dependency order is:

1. `ADRENAL_Evaluation_OperatingCharacteristics.R`
   - One-time simulation pass at the union of look times for `K ∈ {2, …, 10}`,
     `R = 1,000,000` virtual trials, under both `H₀` and `H₁`. Writes
     `ADRENAL_OperatingCharacteristics_H0.rda` and `…_H1.rda`
     (each ~ 450 MB) plus eight 12-panel jpegs covering the four
     operating-characteristic metrics (type I, type II, `E(N | H₀)`,
     `E(N | H₁)`) under both `Binding` and `nonBinding` futility. Five of
     these are copied into supplement Section S4 as Figures S1–S5: the
     `Binding` type-I jpeg becomes Figure S1 (`figS41.jpeg`) and the four
     `nonBinding` jpegs become Figures S2–S5 (`figS42`–`figS45`). Only the
     type-I-error panel differs between the two conventions, so the `Binding`
     type-II and `E(N)` jpegs are graphically identical to their `nonBinding`
     counterparts and are not copied. None of these are main-paper figures;
     the §3.4 figure (Figure 1) is the same-cache demonstration from Stage 8b
     (see the "Futility convention" note under "Notes on the kernel" below).
   - **This is the long-running step (~ tens of CPU-hours).** Everything
     downstream is fast (seconds to minutes) because it reuses these caches.
2. `ADRENAL_Evaluation_MonteCarloError.R`
   - §3.6 simulation-budget study, 100 replicates × R ∈ {1k, 10k, 100k, 1M}.
   - Produces `ADRENAL_MonteCarloError.rda` and the boxplot Figure 3
     (`fig361.jpeg` after Stage 10).
3. `ADRENAL_Diagnostics_PriorApproximation.R`
   - Appendix A.1 prior-approximation diagnostics for the logit-normal example.
   - Produces `ADRENAL_PriorDiagnostics.{rda,jpeg}` (becomes Figure A.1).
3b. `ADRENAL_Diagnostics_OCSensitivity.R`
   - Appendix A.1 diagnostic (iv): operating-characteristic sensitivity to the
     beta-mixture order L. Re-evaluates the §3.5 calibrated three-look design
     under the logit-normal prior at L = 1, 2, 3 (R = 10⁶, common seed) and
     reports type I error, power and E(N) (becomes Table A.2).
   - Produces `ADRENAL_Diagnostics_OCSensitivity.rda`. Runs in ~1-2 minutes.
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
5b. `ADRENAL_Calibration_SeedSensitivity.R` (§3.5 seed-robustness paragraph)
    - Re-evaluates the two manuscript-cited calibrated designs (K = 3 and
      K = 5) at 10 independent seeds (R = 10⁶ each) to pin seed-to-seed
      variability of the calibrated type I and power.
    - Produces `ADRENAL_Calibration_SeedSensitivity.rda`. Runs in ~2 minutes.
5c. `ADRENAL_Calibration_PlausibilitySweep.R` (§3.5 plausibility paragraph)
    - Re-evaluates the two calibrated designs' type I across a wide band of
      control rates (ϑ ∈ {0.15,0.20,…,0.50}, δ=0), reporting both the binding
      and the non-binding (futility-disabled) upper bound from one simulation
      per cell, plus power across a band of effects (δ ∈ {−0.03,…,−0.07},
      ϑ=0.33), at R = 10⁶ per cell. Also reports the type I supremum over the
      nuisance band under each convention.
    - Produces `ADRENAL_Calibration_PlausibilitySweep.rda`. Runs in ~4 minutes.
5d. `ADRENAL_Calibration_FrontierSelection.R` (§3.5 maximum-power frontier)
    - Post-processes the grid-search cache (no new simulation): for each K,
      selects the design maximising power subject to the design-point type I
      ≤ 2.5% and power ≥ 90%, i.e. the maximum-power frontier across looks.
    - Reads `ADRENAL_CalibratedDesign.rda`; produces
      `ADRENAL_Calibration_FrontierSelection.rda`. Runs in < 1 minute.
6. `ADRENAL_Benchmark_UniformPrior.R` (§3.3.1)
   - Head-to-head against `BATSS` + `adaptr` under Beta(1,1). **This is the
     other long-running step (~ tens of CPU-hours) because of `BATSS`+INLA.**
   - Produces `ADRENAL_Benchmark_UniformPrior.rda`.
6b. `ADRENAL_Benchmark_SeedSensitivity.R` (§3.3.1 K=9 multi-seed paragraph)
    - Re-runs the K = 9 matched-budget cell (R = 5,000) across 1000 independent
      seeds (configurable via `N_SEEDS`; seed 21 reproduces the cached value)
      to put empirical bars on the cross-method comparison of Table 1.
    - Produces `ADRENAL_Benchmark_SeedSensitivity.rda`. Runs in under a minute.
7. `ADRENAL_Benchmark_LogitNormalPrior.R` (§3.3.2)
   - Same benchmark under a logit-normal non-conjugate prior. Same scale of
     wall-clock cost as step 6.
8. `ADRENAL_Benchmark_FullGridTiming.R`
   - End-to-end timing of the 438-design grid sweep (per §3.4 and §4).
8b. `ADRENAL_Demo_SameCacheSweep.R`
   - §2.4 same-cache demonstration: evaluates the 438-design grid one design
     at a time against a single cached simulation pass, recording cumulative
     wall-clock, to illustrate Lemma 1 (cache invariance).
   - Produces `ADRENAL_Demo_SameCacheSweep.{rda,jpeg}` (the figure for §2.4).
9. `ADRENAL_Benchmark_GenerateTables.R`
   - Splices the four benchmark tables (Tables 1–4) into the manuscript .tex.
9b. `ADRENAL_Evaluation_StoppingProbabilities.R`
   - Re-derives Supplement Table S2 (per-look stopping probabilities and
     expected sample size for the K = 5, p = 0.99, q = 0.90 binding-futility
     design) directly from the deposited posterior-probability caches. Writes
     `ADRENAL_Evaluation_StoppingProbabilities.rda` and prints LaTeX-ready table
     rows so the supplement table can be regenerated from cached data rather
     than hand-typed. Runs in seconds.
10. `ADRENAL_Summary.R`
    - Copies/renames the OC jpegs into the manuscript and supplement
      directories under their final names (manuscript: `fig341.jpeg`,
      `fig351.jpeg`, `fig361.jpeg`, `figA11.jpeg`, `figA21.jpeg`; supplement Section S4:
      `figS41.jpeg`–`figS45.jpeg`);
      fills the live `TBD_*` placeholders in the .tex with cached numbers;
      writes `post_run_summary.txt` with the numerical deltas referenced
      throughout §3.3, §3.4, §3.5 and §3.6.
11. Stage 11 of `ADRENAL_Replication.sh` rebuilds both the manuscript and
    the supplement PDFs (`pdflatex` × 3 plus `bibtex` for each document).

Two further §3.5 analyses are provided as standalone scripts, not wired into
`ADRENAL_Replication.sh`:

- `ADRENAL_Calibration_FrontierVerification.R` — re-evaluates each per-look
  maximum-power design on its own analysis schedule over the control-rate band
  (binding and non-binding type I suprema, plus power at δ = −0.05) at
  R = 10⁶ per cell, producing the Supplement Section S6 frontier in
  `ADRENAL_Calibration_FrontierVerification.rda`. Long-running.
- `ADRENAL_Calibration_FeasibleFamily.R` — builds the feasible family of
  Haybittle–Peto efficacy-threshold pairs (fixed futility q = 0.90) meeting
  power ≥ 90% and a control-rate-band type I supremum ≤ 2.5% for the headline
  three- and five-look designs, and draws the §3.5 figure. Reads
  `ADRENAL_CalibratedDesign.rda`; produces
  `ADRENAL_Calibration_FeasibleFamily.{rda,jpeg}` (the source of `fig351.jpeg`,
  Figure 2). Requires `ggrepel`.

Each script writes a `<scriptname>_sessionInfo.rds` and `.txt` alongside its
output, capturing the exact R / package versions used. Together with the
single seed (`SEED = 21L`) this is sufficient for byte-level reproduction.

## Files

| File                                              | Purpose                                                            |
|---------------------------------------------------|--------------------------------------------------------------------|
| `adabay_bern.R`                                   | Core simulator and public API.                                     |
| `adabay_bern_bootstrap.R`                         | Find-and-source bootstrap each driver runs to locate the repo root and load the setup. |
| `adabay_bern_setup.R`                             | Shared setup (paths, parallelism, sessionInfo helper).             |
| `adabay_bern_tests.R`                             | Regression tests; run after any change to `adabay_bern.R`.         |
| `ADRENAL_Evaluation_OperatingCharacteristics.R`   | §3.4 OC simulation (long-running).                                 |
| `ADRENAL_Evaluation_MonteCarloError.R`            | §3.6 MC error study.                                               |
| `ADRENAL_Diagnostics_PriorApproximation.R`        | Appendix A.1 prior-approximation diagnostics.                      |
| `ADRENAL_Diagnostics_QuadratureConvergence.R`     | Appendix A.2 Gauss–Legendre quadrature convergence study.          |
| `ADRENAL_Diagnostics_OCSensitivity.R`             | Appendix A.1 diagnostic (iv): OC sensitivity to mixture order L.    |
| `ADRENAL_Calibration_GridSearch.R`                | §3.5 calibrated three- and five-look designs.                      |
| `ADRENAL_Calibration_SeedSensitivity.R`           | §3.5 seed-sensitivity sweep of the calibrated designs.             |
| `ADRENAL_Calibration_PlausibilitySweep.R`         | §3.5 plausibility-set type I / power over control & effect bands.   |
| `ADRENAL_Calibration_FrontierSelection.R`         | §3.5 maximum-power frontier across looks (cache post-processing).   |
| `ADRENAL_Calibration_FrontierVerification.R`      | §3.5 / Suppl. S6 band-verified per-K frontier (standalone).         |
| `ADRENAL_Calibration_FeasibleFamily.R`            | §3.5 feasible-family table + figure (fig351); standalone.           |
| `ADRENAL_Benchmark_UniformPrior.R`                | §3.3.1 head-to-head benchmark (long-running).                      |
| `ADRENAL_Benchmark_SeedSensitivity.R`             | §3.3.1 K=9 matched-budget 1000-seed type-I sweep.                  |
| `ADRENAL_Benchmark_LogitNormalPrior.R`            | §3.3.2 head-to-head benchmark (long-running).                      |
| `ADRENAL_Benchmark_FullGridTiming.R`              | §3.4 / §4 full 438-design grid timing.                             |
| `ADRENAL_Benchmark_GenerateTables.R`              | Generates Tables 1–4 and splices them into the manuscript .tex.    |
| `ADRENAL_Demo_SameCacheSweep.R`                   | §2.4 same-cache demonstration (Lemma 1) figure.                    |
| `ADRENAL_Evaluation_StoppingProbabilities.R`      | Supplement Table S2: per-look stopping probabilities (K = 5).      |
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
  type-I-error panel differs between binding and non-binding (supplement
  Figure S1 vs Figure S2); the type-II and expected-sample-size figures
  (Figures S3–S5) are graphically identical under both conventions.
- **Gauss–Legendre quadrature.** Fixed 128 nodes on (0, 1). For sensitivity
  see the convergence study in Appendix A.2.
- **Parallel RNG.** The kernel installs `RNGkind("L'Ecuyer-CMRG")` at load
  time so that any future `mclapply`-based kernel that draws inside a chunk
  is reproducible.
- **Boundary handling.** For relative-risk and odds-ratio transforms the
  kernel zeroes out quadrature nodes whose transformed value falls outside
  (0, 1). This is correct: those `(θ_c, θ_t)` pairs are outside the support
  of the joint posterior on the unit square, not "lost probability mass."
- **Regression tests.** `adabay_bern_tests.R` covers six checks: (1)
  two-component mixture posterior weight renormalisation, (2) Beta(1,1)
  delta-tail convolution against direct integration, (3) numeric-proximity
  threshold lookup, (4) `options(warn)` restored on exit, (5) the
  L'Ecuyer-CMRG RNG kind, and (6) the within-tolerance rule for mixture-component
  selection. All six should pass in under 30 seconds.
