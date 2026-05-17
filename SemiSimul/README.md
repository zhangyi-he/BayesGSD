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
     (each ~ 450 MB), plus the `Binding` / `nonBinding` Figure 1–4 jpegs.
   - **This is the long-running step (~ tens of CPU-hours).** Everything
     downstream is fast (seconds to minutes) because it reuses these caches.
2. `ADRENAL_Evaluation_MonteCarloError.R`
   - Section 3.5 simulation-budget study, 100 replicates × R ∈ {1k, 10k, 100k, 1M}.
   - Produces `ADRENAL_MonteCarloError.rda` and the boxplot Figure 5.
3. `ADRENAL_Diagnostics_PriorApproximation.R`
   - Appendix B prior-approximation diagnostics for the logit-normal example.
   - Produces `ADRENAL_PriorDiagnostics.jpeg`.
4. `ADRENAL_Diagnostics_QuadratureConvergence.R`
   - Appendix C Gauss–Legendre quadrature convergence study at
     `Q ∈ {64, 128, 256, 512}` under both priors.
   - Produces `ADRENAL_QuadratureConvergence.jpeg` and `_QuadratureConvergence.rda`.
5. `ADRENAL_Calibration_GridSearch.R`
   - Sweeps the calibrated Haybittle–Peto designs of §3.4 over the cached OC.
   - Produces `ADRENAL_CalibratedDesign.rda`.
6. `ADRENAL_Benchmark_UniformPrior.R` (§3.3.1)
   - Head-to-head against `BATSS` + `adaptr` under Beta(1,1). **This is the
     other long-running step (~ tens of CPU-hours) because of `BATSS`+INLA.**
   - Produces `ADRENAL_Benchmark_UniformPrior.rda`.
7. `ADRENAL_Benchmark_LogitNormalPrior.R` (§3.3.2)
   - Same benchmark under a logit-normal non-conjugate prior. Same scale of
     wall-clock cost as step 6.
8. `ADRENAL_Benchmark_FullGridTiming.R`
   - Sub-minute end-to-end timing of the 438-design grid sweep (§3.4 / §4).
9. `ADRENAL_Benchmark_GenerateTables.R`
   - Splices the four benchmark tables (Tables 1–4) into the manuscript .tex.
10. `ADRENAL_Summary.R`
   - Copies/renames figures into the manuscript directory; fills the live
     `TBD_*` placeholders; writes `post_run_summary.txt` with the numerical
     deltas needed for the §3.3 / §3.4 prose.

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
| `ADRENAL_Evaluation_MonteCarloError.R`            | §3.5 MC error study.                                               |
| `ADRENAL_Diagnostics_PriorApproximation.R`        | Appendix B prior approximation diagnostics.                        |
| `ADRENAL_Diagnostics_QuadratureConvergence.R`     | Appendix C Gauss–Legendre quadrature convergence study.            |
| `ADRENAL_Calibration_GridSearch.R`                | §3.4 calibrated three- and five-look designs.                      |
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
- **Gauss–Legendre quadrature.** Fixed 128 nodes on (0, 1). For sensitivity
  see the convergence study in Appendix B.
- **Parallel RNG.** The kernel installs `RNGkind("L'Ecuyer-CMRG")` at load
  time so that any future `mclapply`-based kernel that draws inside a chunk
  is reproducible.
- **Boundary handling.** For relative-risk and odds-ratio transforms the
  kernel zeroes out quadrature nodes whose transformed value falls outside
  (0, 1). This is correct: those `(θ_c, θ_t)` pairs are outside the support
  of the joint posterior on the unit square, not "lost probability mass."

