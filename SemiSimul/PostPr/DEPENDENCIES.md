# Dependencies and reproducibility pinning

This directory ships an `renv.lock` file pinning the exact package versions
used to generate every numerical result, table, and figure in the manuscript
and supplement. `renv.lock` is the standard R format consumed by the
[`renv`](https://rstudio.github.io/renv/) package.

## How to use the lockfile

```r
# Install renv if you don't have it already.
install.packages("renv")

# From the Code/Code v1.0/ directory:
renv::restore()
```

`renv::restore()` will install each package at the exact version listed,
into a project-local library, without touching your system library.

## How the lockfile was produced

The lockfile was assembled by hand from the `_sessionInfo.txt` files in
`Output/Output v1.0/`, which are written by each ADRENAL_*.R driver at the
end of a successful run (via `bayesgsd_save_session()` in
`adabay_bern_setup.R`). The versions are exactly those used to produce
the deposited `.rda` caches and figures.

The lockfile does **not** carry SHA-256 hashes; CRAN versions are resolved
by version string. If your local environment needs fully cryptographic
pinning, run a fresh `renv::snapshot()` after `renv::restore()` to add
hashes from your installed copies.

## Sources

- CRAN — all packages except `INLA`. Default repo: `https://cloud.r-project.org`.
- INLA — installed from `https://inla.r-inla-download.org/R/stable`
  (see <https://www.r-inla.org/download-install>). INLA is required by
  `BATSS` (one of the two comparators in §3.3 of the manuscript).

## Versions captured

The lockfile is anchored to the run that produced `Output v1.0/` and the
deposited figures:

- R `4.5.2 (2025-10-31)`, platform `x86_64-apple-darwin20` (macOS Ventura 13.7.8).
- Headline packages: `BATSS 1.2.0`, `INLA 25.10.19`, `RBesT 1.9-0`,
  `adaptr 1.5.0`, `Rcpp 1.1.1-1.1`, `ggplot2 4.0.3`, `ggsci 5.0.0`,
  `ggrepel 0.9.8`, `scales 1.4.0`, `patchwork 1.3.2`. (`ggrepel` is used by
  `ADRENAL_Calibration_FeasibleFamily.R` for the labelled corner annotations
  in the Section 3.5 feasible-family figure.)
- Full version inventory in `renv.lock`; per-script sessionInfo at
  `Output/Output v1.0/<script>_sessionInfo.txt`.

## Caveats

1. **Base packages are excluded** from the lockfile (`parallel`, `splines`,
   `stats4`, `grid`, `tools`, `compiler`, `methods`, `stats`, `utils`,
   `datasets`, `graphics`, `grDevices`, `base`). These ship with R itself
   at the R version listed; matching the R version is sufficient.
2. **INLA is not on CRAN.** Even with `renv::restore()`, the INLA repository
   must be reachable. The `renv.lock` includes the INLA repo URL so this is
   automatic in most environments, but corporate proxies that block
   non-CRAN repos may require manual install.
3. **`QuickJSR` appears at two versions** across the historical sessionInfo
   files (1.9.2 in earlier runs, 1.10.0 in later runs). The lockfile pins
   1.10.0. Either version reproduces the headline operating-characteristic
   numbers; `QuickJSR` is a transitive dependency of `rstan` and does not
   affect the simulation kernel.
