#' @title Rapid evaluation of the operational and statistical properties of
#'   Bayesian group sequential designs with a Bernoulli endpoint
#' @description Vectorised, Rcpp-accelerated implementation of the
#'   semi-simulation framework. Three changes relative to the legacy
#'   per-trial `integrate()` path: (a) fixed Gauss--Legendre quadrature on
#'   $(0,1)$ replaces adaptive per-trial quadrature so the integrand can be
#'   evaluated for all $R$ virtual trials at once via matrix arithmetic;
#'   (b) the hot inner kernel (mixture density times treatment-arm beta
#'   tail/density, weighted-summed over quadrature nodes) is shipped as an
#'   inline Rcpp routine compiled at module load (with a pure-R fallback if
#'   Rcpp is not available); (c) cross-trial parallelism is delivered by
#'   `parallel::mclapply` over chunks rather than per-stage `makeCluster`
#'   teardown/restart with `foreach %dopar%`. The public API
#'   (`initialisePriorSettings`, `getPosteriorProbabilities`,
#'   `getPosteriorDensities`, `runTrialMonitoring`,
#'   `getOperatingCharacteristics`, `runStageMonitoring`) and the on-disk
#'   `trialSimulation` schema are preserved.
#' @author Zhangyi He, Feng Yu, Suzie Cro, Laurent Billot

suppressPackageStartupMessages({
  library(parallel)
  library(RBesT)
})

# ============================================================================
# Module-level state, quadrature setup, and Rcpp kernel compilation
# ============================================================================

.adabay_env <- new.env(parent = baseenv())

# Inline C++ kernel. Compiled at load time; falls back to pure R if Rcpp is
# unavailable or compilation fails.
.adabay_env$cpp_src <- '
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
// Inner-loop parallelism via OpenMP is opt-in. By default the inner loop
// runs single-threaded, with cross-trial parallelism handled by mclapply
// in the front-end. On platforms where _OPENMP is defined at compile time
// (e.g. Linux with libgomp), this matters because the inner #pragma omp
// parallel for would otherwise stack underneath mclapply and oversubscribe
// the CPU. To enable inner-loop OpenMP, build with
//   BAYESGSD_OPENMP=1 Rscript ...
// (see the .init_adabay() R-side toggle below) and ensure your toolchain
// supports OpenMP. On macOS stock Apple Clang lacks omp.h; install libomp
// via Homebrew and add to ~/.R/Makevars:
//   CXXFLAGS += -Xclang -fopenmp
//   LDFLAGS  += -lomp
#ifndef BAYESGSD_INNER_OMP
#define BAYESGSD_INNER_OMP 0
#endif

// --------------------------------------------------------------------------
// Posterior tail-probability matrix (R trials x T thresholds) via Gauss-
// Legendre quadrature on (0, 1).
//
// Cache strategy: for a fixed stage the per-arm posterior beta parameters
// depend only on the integer event count, which takes values in {0, ..., m}
// for the control arm and {0, ..., n} for the treatment arm. We precompute a
// Q x (m+1) matrix of mixture-density values at quadrature nodes and, for
// each threshold, a Q x (n+1) matrix of mixture-tail values. The per-trial
// loop then collapses to a length-Q weighted dot product with no further
// dbeta/pbeta calls. With R >> m + n this is ~ (R / max(m+1, n+1)) times
// fewer transcendental calls than per-trial computation.
//
// transform(theta, tau) =
//   type == 0 (absolute risk): theta + tau
//   type == 1 (relative risk): theta * tau
//   type == 2 (odds ratio):    (theta * tau) / (1 - theta + theta * tau)
//
// alternative: 0 = "less" (lower tail), 1 = "greater" (upper tail).
// Posterior mixture weights renormalise by the per-component marginal
// likelihood B(a_l + e, b_l + m - e) / B(a_l, b_l), computed in log-space
// (lgamma differences) for numerical stability.
// [[Rcpp::export]]
NumericMatrix adabay_tail_prob(
    NumericVector nodes,
    NumericVector weights,
    IntegerVector ec,
    IntegerVector et,
    int mc,
    int nt,
    NumericVector pac, NumericVector pbc, NumericVector pwc,
    NumericVector pat, NumericVector pbt, NumericVector pwt,
    NumericVector thresholds,
    int type,
    int alternative
) {
  const int R = ec.size();
  const int T = thresholds.size();
  const int Q = nodes.size();
  const int Lc = pwc.size();
  const int Lt = pwt.size();
  NumericMatrix out(R, T);

  // --- Control-arm posterior cache ------------------------------------------
  // For each event count i in [0, mc], compute the posterior mixture
  // weights w_l_post(i) ∝ w_l * B(a_l+i, b_l+mc-i) / B(a_l, b_l), normalised
  // to sum to 1. Then form D[q + Q*i] = ∑_l w_l_post(i) * dbeta(nodes[q]; a_l+i, b_l+mc-i).
  // Cache the prior log-normalising constants once.
  std::vector<double> lB_prior_c((std::size_t)Lc);
  for (int lc = 0; lc < Lc; ++lc) {
    lB_prior_c[lc] = R::lbeta(pac[lc], pbc[lc]);
  }
  std::vector<double> D((std::size_t)Q * ((std::size_t)mc + 1), 0.0);
  std::vector<double> wpost_c((std::size_t)Lc);
  std::vector<double> lwm_c((std::size_t)Lc);
  for (int i = 0; i <= mc; ++i) {
    const double e = (double)i;
    double lmax = R_NegInf;
    for (int lc = 0; lc < Lc; ++lc) {
      const double a = pac[lc] + e;
      const double b = pbc[lc] + (double)mc - e;
      lwm_c[lc] = std::log(pwc[lc]) + R::lbeta(a, b) - lB_prior_c[lc];
      if (lwm_c[lc] > lmax) lmax = lwm_c[lc];
    }
    double zsum = 0.0;
    for (int lc = 0; lc < Lc; ++lc) {
      wpost_c[lc] = std::exp(lwm_c[lc] - lmax);
      zsum += wpost_c[lc];
    }
    for (int lc = 0; lc < Lc; ++lc) wpost_c[lc] /= zsum;
    for (int lc = 0; lc < Lc; ++lc) {
      const double a = pac[lc] + e;
      const double b = pbc[lc] + (double)mc - e;
      const double w = wpost_c[lc];
      const std::size_t base = (std::size_t)Q * (std::size_t)i;
      for (int q = 0; q < Q; ++q) {
        D[base + q] += w * R::dbeta(nodes[q], a, b, 0);
      }
    }
  }
  // Treatment-arm prior log-normalising constants (used for posterior weights
  // inside the per-threshold loop below).
  std::vector<double> lB_prior_t((std::size_t)Lt);
  for (int lt = 0; lt < Lt; ++lt) {
    lB_prior_t[lt] = R::lbeta(pat[lt], pbt[lt]);
  }

  // Per-threshold treatment-arm tail cache: T_t[q + Q*j] for q, j in [0, nt].
  std::vector<double> Tt((std::size_t)Q * ((std::size_t)nt + 1), 0.0);
  const int lower_tail = (alternative == 0) ? 1 : 0;
  std::vector<double> xt((std::size_t)Q);
  std::vector<unsigned char> valid((std::size_t)Q);
  std::vector<double> bound_tail((std::size_t)Q);

  for (int t = 0; t < T; ++t) {
    const double tau = thresholds[t];
    // Compute transformed node values once per threshold.
    for (int q = 0; q < Q; ++q) {
      const double theta = nodes[q];
      double x;
      bool ok = true;
      if (type == 0) {
        x = theta + tau;
      } else if (type == 1) {
        x = theta * tau;
      } else {
        const double denom = 1.0 - theta + theta * tau;
        if (denom <= 0.0) { ok = false; x = 0.0; }
        else x = (theta * tau) / denom;
      }
      xt[q] = x;
      valid[q] = ok ? 1u : 0u;
      // Resolve out-of-(0,1) boundary cases analytically (constant in j).
      if (!ok) bound_tail[q] = 0.0;
      else if (x <= 0.0) bound_tail[q] = (lower_tail ? 0.0 : 1.0);
      else if (x >= 1.0) bound_tail[q] = (lower_tail ? 1.0 : 0.0);
      else bound_tail[q] = -1.0;   // sentinel: needs pbeta evaluation per j
    }

    // Reset Tt for this threshold, then for each event count j compute the
    // posterior mixture weights and accumulate the per-component pbeta tail.
    std::fill(Tt.begin(), Tt.end(), 0.0);
    std::vector<double> wpost_t((std::size_t)Lt);
    std::vector<double> lwm_t((std::size_t)Lt);
    for (int j = 0; j <= nt; ++j) {
      const double e = (double)j;
      double lmax = R_NegInf;
      for (int lt = 0; lt < Lt; ++lt) {
        const double a = pat[lt] + e;
        const double b = pbt[lt] + (double)nt - e;
        lwm_t[lt] = std::log(pwt[lt]) + R::lbeta(a, b) - lB_prior_t[lt];
        if (lwm_t[lt] > lmax) lmax = lwm_t[lt];
      }
      double zsum = 0.0;
      for (int lt = 0; lt < Lt; ++lt) {
        wpost_t[lt] = std::exp(lwm_t[lt] - lmax);
        zsum += wpost_t[lt];
      }
      for (int lt = 0; lt < Lt; ++lt) wpost_t[lt] /= zsum;
      const std::size_t base = (std::size_t)Q * (std::size_t)j;
      for (int lt = 0; lt < Lt; ++lt) {
        const double a = pat[lt] + e;
        const double b = pbt[lt] + (double)nt - e;
        const double w = wpost_t[lt];
        for (int q = 0; q < Q; ++q) {
          double tail;
          if (bound_tail[q] >= 0.0) {
            tail = bound_tail[q];
          } else {
            tail = R::pbeta(xt[q], a, b, lower_tail, 0);
          }
          Tt[base + q] += w * tail;
        }
      }
    }

    // Per-trial weighted dot product (no transcendental work).
    double* outcol = &out(0, t);
#if defined(_OPENMP) && BAYESGSD_INNER_OMP
    #pragma omp parallel for schedule(static)
#endif
    for (int r = 0; r < R; ++r) {
      const int ic = ec[r];
      const int jt = et[r];
      const std::size_t Dbase = (std::size_t)Q * (std::size_t)ic;
      const std::size_t Tbase = (std::size_t)Q * (std::size_t)jt;
      double acc = 0.0;
      for (int q = 0; q < Q; ++q) {
        acc += weights[q] * D[Dbase + q] * Tt[Tbase + q];
      }
      outcol[r] = acc;
    }
  }
  return out;
}

// --------------------------------------------------------------------------
// Posterior density matrix (R x T) via Gauss-Legendre quadrature on (0, 1).
// Caches the per-event-count mixture density values at quadrature nodes for
// both arms (control once, treatment per threshold) and reduces the per-trial
// inner loop to a length-Q dot product. |J| depends on type:
//   type == 0: 1
//   type == 1: theta
//   type == 2: theta(1-theta) / (1 - theta + theta*delta)^2
// [[Rcpp::export]]
NumericMatrix adabay_density(
    NumericVector nodes,
    NumericVector weights,
    IntegerVector ec,
    IntegerVector et,
    int mc,
    int nt,
    NumericVector pac, NumericVector pbc, NumericVector pwc,
    NumericVector pat, NumericVector pbt, NumericVector pwt,
    NumericVector deltas,
    int type
) {
  const int R = ec.size();
  const int T = deltas.size();
  const int Q = nodes.size();
  const int Lc = pwc.size();
  const int Lt = pwt.size();
  NumericMatrix out(R, T);

  // Posterior mixture weights renormalised by per-component marginal
  // likelihood; see comment on `adabay_tail_prob` for details.
  std::vector<double> lB_prior_c((std::size_t)Lc);
  for (int lc = 0; lc < Lc; ++lc) lB_prior_c[lc] = R::lbeta(pac[lc], pbc[lc]);
  std::vector<double> lB_prior_t((std::size_t)Lt);
  for (int lt = 0; lt < Lt; ++lt) lB_prior_t[lt] = R::lbeta(pat[lt], pbt[lt]);

  std::vector<double> D((std::size_t)Q * ((std::size_t)mc + 1), 0.0);
  std::vector<double> wpost_c((std::size_t)Lc);
  std::vector<double> lwm_c((std::size_t)Lc);
  for (int i = 0; i <= mc; ++i) {
    const double e = (double)i;
    double lmax = R_NegInf;
    for (int lc = 0; lc < Lc; ++lc) {
      const double a = pac[lc] + e;
      const double b = pbc[lc] + (double)mc - e;
      lwm_c[lc] = std::log(pwc[lc]) + R::lbeta(a, b) - lB_prior_c[lc];
      if (lwm_c[lc] > lmax) lmax = lwm_c[lc];
    }
    double zsum = 0.0;
    for (int lc = 0; lc < Lc; ++lc) {
      wpost_c[lc] = std::exp(lwm_c[lc] - lmax);
      zsum += wpost_c[lc];
    }
    for (int lc = 0; lc < Lc; ++lc) wpost_c[lc] /= zsum;
    for (int lc = 0; lc < Lc; ++lc) {
      const double a = pac[lc] + e;
      const double b = pbc[lc] + (double)mc - e;
      const double w = wpost_c[lc];
      const std::size_t base = (std::size_t)Q * (std::size_t)i;
      for (int q = 0; q < Q; ++q) D[base + q] += w * R::dbeta(nodes[q], a, b, 0);
    }
  }

  std::vector<double> Dt((std::size_t)Q * ((std::size_t)nt + 1), 0.0);
  std::vector<double> xt((std::size_t)Q);
  std::vector<double> jac((std::size_t)Q);
  std::vector<unsigned char> valid((std::size_t)Q);
  std::vector<double> wpost_t((std::size_t)Lt);
  std::vector<double> lwm_t((std::size_t)Lt);

  for (int t = 0; t < T; ++t) {
    const double delta = deltas[t];
    for (int q = 0; q < Q; ++q) {
      const double theta = nodes[q];
      double x = 0.0, j = 1.0;
      bool ok = true;
      if (type == 0) {
        x = theta + delta;
      } else if (type == 1) {
        x = theta * delta;
        j = theta;
      } else {
        const double denom = 1.0 - theta + theta * delta;
        if (denom <= 0.0) { ok = false; }
        else { x = (theta * delta) / denom; j = theta * (1.0 - theta) / (denom * denom); }
      }
      if (!ok || x < 0.0 || x > 1.0) {
        valid[q] = 0u;
        xt[q]    = 0.0;
        jac[q]   = 0.0;
      } else {
        valid[q] = 1u;
        xt[q]    = x;
        jac[q]   = j;
      }
    }

    std::fill(Dt.begin(), Dt.end(), 0.0);
    for (int j = 0; j <= nt; ++j) {
      const double e = (double)j;
      double lmax = R_NegInf;
      for (int lt = 0; lt < Lt; ++lt) {
        const double a = pat[lt] + e;
        const double b = pbt[lt] + (double)nt - e;
        lwm_t[lt] = std::log(pwt[lt]) + R::lbeta(a, b) - lB_prior_t[lt];
        if (lwm_t[lt] > lmax) lmax = lwm_t[lt];
      }
      double zsum = 0.0;
      for (int lt = 0; lt < Lt; ++lt) {
        wpost_t[lt] = std::exp(lwm_t[lt] - lmax);
        zsum += wpost_t[lt];
      }
      for (int lt = 0; lt < Lt; ++lt) wpost_t[lt] /= zsum;
      const std::size_t base = (std::size_t)Q * (std::size_t)j;
      for (int lt = 0; lt < Lt; ++lt) {
        const double a = pat[lt] + e;
        const double b = pbt[lt] + (double)nt - e;
        const double w = wpost_t[lt];
        for (int q = 0; q < Q; ++q) {
          if (!valid[q]) continue;
          Dt[base + q] += w * R::dbeta(xt[q], a, b, 0);
        }
      }
    }

    double* outcol = &out(0, t);
#if defined(_OPENMP) && BAYESGSD_INNER_OMP
    #pragma omp parallel for schedule(static)
#endif
    for (int r = 0; r < R; ++r) {
      const int ic = ec[r];
      const int jt = et[r];
      const std::size_t Dbase  = (std::size_t)Q * (std::size_t)ic;
      const std::size_t Dtbase = (std::size_t)Q * (std::size_t)jt;
      double acc = 0.0;
      for (int q = 0; q < Q; ++q) {
        acc += weights[q] * D[Dbase + q] * Dt[Dtbase + q] * jac[q];
      }
      outcol[r] = acc;
    }
  }
  return out;
}
'

.init_adabay <- function() {
  # Use an RNG that supports reproducible parallel streams. mclapply forks
  # workers from the parent, so under the default Mersenne-Twister the child
  # state is unpredictable. L'Ecuyer-CMRG combined with `mc.set.seed = TRUE`
  # (the mclapply default) gives each child a deterministic substream.
  # The current kernel only samples in the parent, so this is forward-looking
  # protection; any future kernel that draws inside an mclapply chunk will
  # be reproducible without further changes.
  if (!identical(RNGkind()[1L], "L'Ecuyer-CMRG")) {
    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  }
  # Precompute Gauss-Legendre quadrature on (0, 1) via Golub-Welsch.
  n <- as.integer(getOption("BayesGSD.quadrature.n", 128L))
  k <- seq_len(n - 1L)
  beta_k <- k / sqrt(4 * k^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(2:n, seq_len(n - 1L))] <- beta_k
  J[cbind(seq_len(n - 1L), 2:n)] <- beta_k
  eig <- eigen(J, symmetric = TRUE)
  ord <- order(eig$values)
  nodes_pm1 <- eig$values[ord]
  wts <- 2 * eig$vectors[1L, ord]^2
  .adabay_env$gl_nodes   <- (1 + nodes_pm1) / 2
  .adabay_env$gl_weights <- wts / 2

  # Try to compile the Rcpp kernel; fall back to pure R on failure.
  # `BAYESGSD_OPENMP=1` enables the inner-loop OpenMP pragmas; otherwise the
  # inner kernel runs single-threaded (cross-trial parallelism is handled by
  # mclapply in the front-end).
  use_inner_omp <- identical(Sys.getenv("BAYESGSD_OPENMP", unset = "0"), "1")
  cxx_flags <- if (use_inner_omp) "-DBAYESGSD_INNER_OMP=1" else NULL
  .adabay_env$use_rcpp <- FALSE
  if (requireNamespace("Rcpp", quietly = TRUE)) {
    tryCatch({
      env_args <- list(code     = .adabay_env$cpp_src,
                       env      = .adabay_env,
                       verbose  = FALSE,
                       cacheDir = tempdir(),
                       rebuild  = FALSE)
      if (!is.null(cxx_flags)) {
        old <- Sys.getenv("PKG_CXXFLAGS", unset = NA_character_)
        Sys.setenv(PKG_CXXFLAGS = paste(Sys.getenv("PKG_CXXFLAGS"), cxx_flags))
        on.exit(if (is.na(old)) Sys.unsetenv("PKG_CXXFLAGS")
                else Sys.setenv(PKG_CXXFLAGS = old), add = TRUE)
      }
      do.call(Rcpp::sourceCpp, env_args)
      .adabay_env$use_rcpp <- TRUE
    }, error = function(e) {
      message("adabay: Rcpp kernel unavailable; using pure-R quadrature. (",
              conditionMessage(e), ")")
      .adabay_env$use_rcpp <- FALSE
    })
  }
}
.init_adabay()

# ============================================================================
# Internal coding helpers and pure-R fallback kernels
# ============================================================================

.type_code <- function(type) {
  switch(type,
         "absolute.risk" = 0L,
         "relative.risk" = 1L,
         "odds.ratio"    = 2L,
         stop("Unknown effect type: ", type))
}

.alternative_code <- function(alternative) {
  if (identical(alternative, "greater")) {
    1L
  } else if (identical(alternative, "less")) {
    0L
  } else {
    stop("Unknown alternative: ", alternative)
  }
}

# Pure-R Gauss-Legendre quadrature for the tail-probability matrix. Mirrors
# `adabay_tail_prob` and is used only when the Rcpp kernel is unavailable.
# Posterior mixture weights are renormalised by per-component marginal
# likelihood: w_l_post(e, m) ∝ w_l * B(a_l + e, b_l + m - e) / B(a_l, b_l).
.posterior_mixture_weights <- function(e_vec, m, alpha, beta, weight) {
  L <- length(weight)
  if (L == 1L) return(matrix(1, nrow = length(e_vec), ncol = 1L))
  lB_prior <- lbeta(alpha, beta)  # length L
  # lwm[r, l] = log w_l + lbeta(a_l + e_r, b_l + m - e_r) - lbeta(a_l, b_l)
  lwm <- matrix(NA_real_, nrow = length(e_vec), ncol = L)
  for (l in seq_len(L)) {
    lwm[, l] <- log(weight[l]) + lbeta(alpha[l] + e_vec, beta[l] + (m - e_vec)) - lB_prior[l]
  }
  lmax <- apply(lwm, 1L, max)
  w    <- exp(lwm - lmax)
  w / rowSums(w)
}

# Inner kernel: same algorithm as `.tail_prob_R` but on a single trial chunk.
# Allocates R*Q intermediate matrices, which is fine for the chunked path
# below (chunk_R is bounded by getOption("BayesGSD.fallback_chunk")) but
# would blow up memory on large R; never call directly.
.tail_prob_R_inner <- function(nodes, weights, ec, et, mc, nt,
                               pac, pbc, pwc, pat, pbt, pwt,
                               thresholds, type, alternative) {
  R <- length(ec)
  T <- length(thresholds)
  Q <- length(nodes)
  Lc <- length(pwc)
  Lt <- length(pwt)
  lower_tail <- (alternative == 0L)

  # Renormalised posterior mixture weights for each arm, per trial (R x L).
  wpost_c <- .posterior_mixture_weights(ec, mc, pac, pbc, pwc)
  wpost_t <- .posterior_mixture_weights(et, nt, pat, pbt, pwt)

  # Control-arm posterior density at each node, summed over mixture components
  # with per-trial posterior weights. Vectorised over R x Q grids.
  node_grid <- rep(nodes, each = R)
  dens_c <- matrix(0, nrow = R, ncol = Q)
  for (lc in seq_len(Lc)) {
    a <- rep(pac[lc] + ec,           Q)
    b <- rep(pbc[lc] + (mc - ec),    Q)
    dens_c <- dens_c + wpost_c[, lc] *
      matrix(dbeta(node_grid, shape1 = a, shape2 = b), R, Q)
  }

  out <- matrix(0, R, T)
  for (t in seq_len(T)) {
    tau <- thresholds[t]
    if (type == 0L) {
      xt <- nodes + tau
    } else if (type == 1L) {
      xt <- nodes * tau
    } else {
      denom <- 1 - nodes + nodes * tau
      xt <- ifelse(denom <= 0, NA_real_, (nodes * tau) / denom)
    }
    # Tail at each node (clamped to [0, 1]).
    xt_clamped <- pmin(pmax(xt, 0), 1)
    xt_grid <- rep(xt_clamped, each = R)
    tail_t <- matrix(0, R, Q)
    for (lt in seq_len(Lt)) {
      a <- rep(pat[lt] + et,        Q)
      b <- rep(pbt[lt] + (nt - et), Q)
      tail_t <- tail_t + wpost_t[, lt] *
        matrix(pbeta(xt_grid, a, b, lower.tail = lower_tail), R, Q)
    }
    # Boundary cases: xt <= 0 / xt >= 1 / xt = NA (odds-ratio numerator <= 0).
    # `pbeta` already returns 0 or 1 at the bounds; NA propagates and we zero
    # those columns explicitly so out(r,t) is well-defined.
    bad <- is.na(xt) | xt <= 0 | xt >= 1
    if (any(bad)) {
      # Replace tail entries for `bad` nodes with the analytic boundary value.
      bound_val <- ifelse(is.na(xt), 0,
                          ifelse(xt <= 0, ifelse(lower_tail, 0, 1),
                                 ifelse(lower_tail, 1, 0)))
      bad_cols <- which(bad)
      tail_t[, bad_cols] <- rep(bound_val[bad_cols], each = R)
    }
    out[, t] <- (dens_c * tail_t) %*% weights
  }
  out
}

# Pure-R fallback front-end: chunks the R trials so peak memory stays bounded
# at chunk_R * Q * 16 bytes regardless of total R. The default
# BayesGSD.fallback_chunk = 10,000 keeps the intermediate at ~ 10 MB for
# Q = 128, suitable for any commodity machine.
.tail_prob_R <- function(nodes, weights, ec, et, mc, nt,
                         pac, pbc, pwc, pat, pbt, pwt,
                         thresholds, type, alternative) {
  R <- length(ec)
  T <- length(thresholds)
  chunk_R <- as.integer(getOption("BayesGSD.fallback_chunk", 10000L))
  if (R <= chunk_R) {
    return(.tail_prob_R_inner(nodes, weights, ec, et, mc, nt,
                              pac, pbc, pwc, pat, pbt, pwt,
                              thresholds, type, alternative))
  }
  out <- matrix(0, R, T)
  starts <- seq.int(1L, R, by = chunk_R)
  for (s in starts) {
    idx <- seq.int(s, min(s + chunk_R - 1L, R))
    out[idx, ] <- .tail_prob_R_inner(nodes, weights, ec[idx], et[idx], mc, nt,
                                     pac, pbc, pwc, pat, pbt, pwt,
                                     thresholds, type, alternative)
  }
  out
}

.density_R_inner <- function(nodes, weights, ec, et, mc, nt,
                             pac, pbc, pwc, pat, pbt, pwt,
                             deltas, type) {
  R <- length(ec)
  T <- length(deltas)
  Q <- length(nodes)
  Lc <- length(pwc)
  Lt <- length(pwt)

  wpost_c <- .posterior_mixture_weights(ec, mc, pac, pbc, pwc)
  wpost_t <- .posterior_mixture_weights(et, nt, pat, pbt, pwt)

  node_grid <- rep(nodes, each = R)
  dens_c <- matrix(0, R, Q)
  for (lc in seq_len(Lc)) {
    a <- rep(pac[lc] + ec,        Q)
    b <- rep(pbc[lc] + (mc - ec), Q)
    dens_c <- dens_c + wpost_c[, lc] *
      matrix(dbeta(node_grid, a, b), R, Q)
  }

  out <- matrix(0, R, T)
  for (t in seq_len(T)) {
    delta <- deltas[t]
    if (type == 0L) {
      xt <- nodes + delta
      jac <- rep(1, Q)
    } else if (type == 1L) {
      xt <- nodes * delta
      jac <- nodes
    } else {
      denom <- 1 - nodes + nodes * delta
      xt <- ifelse(denom <= 0, NA_real_, (nodes * delta) / denom)
      jac <- ifelse(denom <= 0, 0, nodes * (1 - nodes) / denom^2)
    }
    valid <- !is.na(xt) & xt >= 0 & xt <= 1
    if (!any(valid)) next
    xt_v <- xt
    xt_v[!valid] <- 0   # placeholder; multiplied by 0 below
    xt_grid <- rep(xt_v, each = R)
    dens_t <- matrix(0, R, Q)
    for (lt in seq_len(Lt)) {
      a <- rep(pat[lt] + et,        Q)
      b <- rep(pbt[lt] + (nt - et), Q)
      dens_t <- dens_t + wpost_t[, lt] *
        matrix(dbeta(xt_grid, a, b), R, Q)
    }
    # Zero out invalid columns.
    if (any(!valid)) dens_t[, !valid] <- 0
    out[, t] <- (dens_c * dens_t) %*% (weights * jac)
  }
  out
}

.density_R <- function(nodes, weights, ec, et, mc, nt,
                       pac, pbc, pwc, pat, pbt, pwt,
                       deltas, type) {
  R <- length(ec)
  T <- length(deltas)
  chunk_R <- as.integer(getOption("BayesGSD.fallback_chunk", 10000L))
  if (R <= chunk_R) {
    return(.density_R_inner(nodes, weights, ec, et, mc, nt,
                            pac, pbc, pwc, pat, pbt, pwt,
                            deltas, type))
  }
  out <- matrix(0, R, T)
  starts <- seq.int(1L, R, by = chunk_R)
  for (s in starts) {
    idx <- seq.int(s, min(s + chunk_R - 1L, R))
    out[idx, ] <- .density_R_inner(nodes, weights, ec[idx], et[idx], mc, nt,
                                   pac, pbc, pwc, pat, pbt, pwt,
                                   deltas, type)
  }
  out
}

# --------------------------------------------------------------------------
# Front-end dispatch: chunk over R if requested, route to Rcpp or R kernel.
# --------------------------------------------------------------------------

.posterior_tail_prob <- function(ec, et, mc, nt, priors_c, priors_t,
                                  thresholds, type, alternative,
                                  ..., chunk_size = NULL, cores = NULL) {
  nodes   <- .adabay_env$gl_nodes
  weights <- .adabay_env$gl_weights
  pac <- as.numeric(priors_c$alpha);  pbc <- as.numeric(priors_c$beta);  pwc <- as.numeric(priors_c$weight)
  pat <- as.numeric(priors_t$alpha);  pbt <- as.numeric(priors_t$beta);  pwt <- as.numeric(priors_t$weight)
  type_i <- .type_code(type)
  alt_i  <- .alternative_code(alternative)
  ec_i <- as.integer(ec);  et_i <- as.integer(et)

  fn <- if (isTRUE(.adabay_env$use_rcpp)) {
    .adabay_env$adabay_tail_prob
  } else {
    .tail_prob_R
  }

  R <- length(ec_i)
  cs <- if (is.null(chunk_size)) {
    getOption("BayesGSD.chunk_size", 50000L)
  } else as.integer(chunk_size)
  cr <- if (is.null(cores)) {
    getOption("BayesGSD.cores", max(1L, parallel::detectCores() - 1L))
  } else as.integer(cores)

  if (R <= cs || cr <= 1L) {
    return(fn(nodes, weights, ec_i, et_i, mc, nt,
              pac, pbc, pwc, pat, pbt, pwt,
              as.numeric(thresholds), type_i, alt_i))
  }
  # Chunked path: split trials across `cr` cores via mclapply.
  starts <- seq.int(1L, R, by = cs)
  chunks <- lapply(starts, function(s) seq.int(s, min(s + cs - 1L, R)))
  parts <- parallel::mclapply(chunks, function(idx) {
    fn(nodes, weights, ec_i[idx], et_i[idx], mc, nt,
       pac, pbc, pwc, pat, pbt, pwt,
       as.numeric(thresholds), type_i, alt_i)
  }, mc.cores = cr, mc.preschedule = TRUE)
  do.call(rbind, parts)
}

.posterior_density <- function(ec, et, mc, nt, priors_c, priors_t,
                                deltas, type,
                                ..., chunk_size = NULL, cores = NULL) {
  nodes   <- .adabay_env$gl_nodes
  weights <- .adabay_env$gl_weights
  pac <- as.numeric(priors_c$alpha);  pbc <- as.numeric(priors_c$beta);  pwc <- as.numeric(priors_c$weight)
  pat <- as.numeric(priors_t$alpha);  pbt <- as.numeric(priors_t$beta);  pwt <- as.numeric(priors_t$weight)
  type_i <- .type_code(type)
  ec_i <- as.integer(ec);  et_i <- as.integer(et)

  fn <- if (isTRUE(.adabay_env$use_rcpp)) {
    .adabay_env$adabay_density
  } else {
    .density_R
  }

  R <- length(ec_i)
  cs <- if (is.null(chunk_size)) getOption("BayesGSD.chunk_size", 50000L) else as.integer(chunk_size)
  cr <- if (is.null(cores))      getOption("BayesGSD.cores",      max(1L, parallel::detectCores() - 1L)) else as.integer(cores)

  if (R <= cs || cr <= 1L) {
    return(fn(nodes, weights, ec_i, et_i, mc, nt,
              pac, pbc, pwc, pat, pbt, pwt,
              as.numeric(deltas), type_i))
  }
  starts <- seq.int(1L, R, by = cs)
  chunks <- lapply(starts, function(s) seq.int(s, min(s + cs - 1L, R)))
  parts <- parallel::mclapply(chunks, function(idx) {
    fn(nodes, weights, ec_i[idx], et_i[idx], mc, nt,
       pac, pbc, pwc, pat, pbt, pwt,
       as.numeric(deltas), type_i)
  }, mc.cores = cr, mc.preschedule = TRUE)
  do.call(rbind, parts)
}

# ============================================================================
# Public API
# ============================================================================

#' Initialise prior settings
#' @param priorSamples list of two vectors of prior samples (control, treatment); NULL gives the non-informative Beta(1,1) prior on each arm
#' @param numberOfComponents maximum number of beta-mixture components fitted by RBesT::automixfit when priorSamples is supplied
#' @param selection mixture-component selection rule. One of:
#'   * "within_tolerance" (default; matches Section 2.3 of the manuscript): the smallest
#'     L whose empirical forward KL divergence is within `tolerance` of the
#'     minimum achieved over `L = 1, ..., numberOfComponents`. This is
#'     well-defined even when the empirical KL is non-monotone in L (EM
#'     stochasticity) and never selects a mixture materially worse than the
#'     best available.
#'   * "aic": AIC-best mixture returned by RBesT::automixfit (its built-in default).
#'   * "max": legacy behaviour — pick the L = `numberOfComponents` fit.
#'     Provided for backward compatibility with pre-2026-05 results; not
#'     recommended.
#' @param tolerance selection tolerance for the within-tolerance rule (the value
#'   of epsilon in Section 2.3 of the manuscript). Default 1e-3 as adopted
#'   throughout the paper.
initialisePriorSettings <- function(priorSamples = NULL,
                                    numberOfComponents = 4L,
                                    selection = c("within_tolerance", "aic", "max"),
                                    tolerance = 1e-3) {
  selection <- match.arg(selection)
  if (is.null(priorSamples)) {
    return(list(
      list(weight = 1, alpha = 1, beta = 1),
      list(weight = 1, alpha = 1, beta = 1)
    ))
  }
  priors <- vector(mode = "list", length = length(priorSamples))
  for (i in seq_along(priorSamples)) {
    if (all(is.na(priorSamples[[i]]))) {
      priors[[i]] <- list(weight = 1, alpha = 1, beta = 1)
    } else {
      # Fit one mixture per candidate component count L = 1, ..., Lmax, so
      # that the empirical forward KL is a genuine function of L (the object
      # the Section 2.3 selection rule operates on).
      fits <- lapply(seq_len(numberOfComponents), function(L)
        automixfit(sample = priorSamples[[i]], Nc = L, type = "beta"))
      bm <- switch(selection,
                   "aic"   = automixfit(sample = priorSamples[[i]],
                                        Nc = 1:numberOfComponents, type = "beta"),
                   "max"   = fits[[numberOfComponents]],
                   "within_tolerance" = .select_within_tolerance(fits, priorSamples[[i]], tolerance))
      # ESS is the weighted sum of per-component prior ESS values
      # (a_l + b_l for Beta(a_l, b_l)), i.e. the mixture-weight-weighted mean
      # of the component ESS. This is a diagnostic-only summary used in the
      # prior-approximation report; it is NOT the moment-matching mixture ESS
      # of RBesT::ess(). For the RBesT definition, call RBesT::ess(bm) on the
      # returned mixture object directly. The diagnostic-only slot here is
      # not consumed by the downstream operating-characteristic pipeline.
      priors[[i]] <- list(
        weight = bm[1, ], alpha = bm[2, ], beta = bm[3, ],
        ESS = sum(bm[1, ] * (bm[2, ] + bm[3, ]))
      )
    }
  }
  priors
}

# Mixture-size selection rule implementing Section 2.3 of the manuscript:
# the smallest L whose empirical forward KL divergence is within `epsilon`
# of the minimum achieved over the candidate sizes. Since
# D_KL(p || p_hat_L) = -E_p[log p_hat_L] + E_p[log p] and E_p[log p] is
# constant in L, the rule is equivalent to one on the empirical cross-entropy
# CE_L = -(1/N) sum_i log p_hat_L(theta_i), evaluated on the prior sample.
# `fits` is a list of fitted beta mixtures, one per candidate L (L = index).
# The rule is well-defined even when CE_L is non-monotone in L and never
# selects a mixture materially worse than the best available.
.select_within_tolerance <- function(fits, samples, epsilon) {
  ce <- vapply(fits, function(bm) {
    lw <- log(bm[1, ])
    a  <- bm[2, ]
    b  <- bm[3, ]
    # log p_hat_L(theta_i) = log sum_l w_l * Beta(theta_i; a_l, b_l)
    log_p_per_sample <- vapply(samples, function(s) {
      lc <- lw + dbeta(s, a, b, log = TRUE)
      m  <- max(lc)
      m + log(sum(exp(lc - m)))
    }, numeric(1L))
    -mean(log_p_per_sample)
  }, numeric(1L))
  # Within-epsilon-of-best: smallest L with CE_L - min_L' CE_L' < epsilon.
  # which(...)[1] is always defined (the global argmin satisfies it), so the
  # rule is well-defined regardless of monotonicity.
  sel <- which(ce - min(ce) < epsilon)[1L]
  fits[[sel]]
}

#' Get posterior tail probabilities for a single trial across one or more thresholds.
#' Returns a numeric vector of length `length(treatmentEffects$size)`.
getPosteriorProbabilities <- function(treatmentEffects,
                                      numberOfEvents,
                                      numberOfSubjects,
                                      alternative,
                                      ...,
                                      priors = NULL) {
  if (is.null(priors)) priors <- initialisePriorSettings()
  out <- .posterior_tail_prob(
    ec = numberOfEvents[1], et = numberOfEvents[2],
    mc = numberOfSubjects[1], nt = numberOfSubjects[2],
    priors_c = priors[[1]], priors_t = priors[[2]],
    thresholds = treatmentEffects$size, type = treatmentEffects$type,
    alternative = alternative,
    chunk_size = .Machine$integer.max, cores = 1L
  )
  as.numeric(out[1, ])
}

#' Get posterior densities for a single trial across one or more effect sizes.
#' Returns a numeric vector of length `length(treatmentEffects$size)`.
getPosteriorDensities <- function(treatmentEffects,
                                  numberOfEvents,
                                  numberOfSubjects,
                                  ...,
                                  priors = NULL) {
  if (is.null(priors)) priors <- initialisePriorSettings()
  out <- .posterior_density(
    ec = numberOfEvents[1], et = numberOfEvents[2],
    mc = numberOfSubjects[1], nt = numberOfSubjects[2],
    priors_c = priors[[1]], priors_t = priors[[2]],
    deltas = treatmentEffects$size, type = treatmentEffects$type,
    chunk_size = .Machine$integer.max, cores = 1L
  )
  as.numeric(out[1, ])
}

# --------------------------------------------------------------------------
# Cross-trial helper: posterior probabilities at a single stage for a vector
# of trials. Returns the legacy list-of-matrices structure expected by
# getOperatingCharacteristics.
# --------------------------------------------------------------------------

.per_stage_posterior_probabilities <- function(numberOfEvents,
                                                numberOfSubjects,
                                                posteriorSettings,
                                                priorSettings) {
  R <- nrow(numberOfEvents)
  ec <- numberOfEvents[, 1L]
  et <- numberOfEvents[, 2L]
  mc <- numberOfSubjects[1L]
  nt <- numberOfSubjects[2L]

  out <- list(efficacy = NULL, equivalence = NULL, futility = NULL)
  alternative <- posteriorSettings$alternative
  eff <- posteriorSettings$effectThresholds$efficacy
  equ <- posteriorSettings$effectThresholds$equivalence
  fut <- posteriorSettings$effectThresholds$futility

  if (!is.null(eff)) {
    m <- .posterior_tail_prob(
      ec, et, mc, nt, priorSettings[[1]], priorSettings[[2]],
      thresholds = eff$size, type = eff$type, alternative = alternative
    )
    # Rows = thresholds, columns = trials (preserves legacy schema).
    m <- t(m)
    rownames(m) <- as.character(eff$size)
    # Robust numeric key alongside the rowname for index-based lookup in
    # `.trigger_stopping`, immune to as.character() precision drift.
    attr(m, "thresholds") <- as.numeric(eff$size)
    out$efficacy <- m
  }
  if (!is.null(equ)) {
    # Legacy equivalence representation: two side-by-side thresholds
    # (lower bound and upper bound) per row.
    lo <- equ$size[, 1L]; hi <- equ$size[, 2L]
    m_lo <- t(.posterior_tail_prob(
      ec, et, mc, nt, priorSettings[[1]], priorSettings[[2]],
      thresholds = lo, type = equ$type, alternative = alternative
    ))
    m_hi <- t(.posterior_tail_prob(
      ec, et, mc, nt, priorSettings[[1]], priorSettings[[2]],
      thresholds = hi, type = equ$type, alternative = alternative
    ))
    m <- m_hi - m_lo
    rownames(m) <- paste(lo, hi, sep = ",")
    attr(m, "thresholds") <- cbind(as.numeric(lo), as.numeric(hi))
    out$equivalence <- m
  }
  if (!is.null(fut)) {
    m <- t(.posterior_tail_prob(
      ec, et, mc, nt, priorSettings[[1]], priorSettings[[2]],
      thresholds = fut$size, type = fut$type, alternative = alternative
    ))
    rownames(m) <- as.character(fut$size)
    attr(m, "thresholds") <- as.numeric(fut$size)
    out$futility <- m
  }
  out
}

# Match a vector of numeric query thresholds against a numeric stored vector
# by proximity. Returns integer indices into `stored`. Fails loudly if any
# query has no match within `tol`; `tol` is generous (1e-8) because the
# stored thresholds were typed in by the user and are not at machine epsilon.
.match_thresholds <- function(query, stored, tol = 1e-8, what = "threshold") {
  query  <- as.numeric(query)
  stored <- as.numeric(stored)
  idx <- vapply(query, function(q) {
    j <- which(abs(stored - q) <= tol)
    if (length(j) == 0L) NA_integer_
    else if (length(j) > 1L) j[which.min(abs(stored[j] - q))]
    else j
  }, integer(1L))
  if (anyNA(idx)) {
    bad <- query[is.na(idx)]
    stop(sprintf(
      "%s value(s) not found in stored cache (tol = %g): %s. Available: %s",
      what, tol, paste(bad, collapse = ", "),
      paste(stored, collapse = ", ")
    ))
  }
  idx
}

#' Run trial monitoring (semi-simulation).
#' @param simulationSettings list with numberOfTrials, numberOfSubjects, allocationRatio, rates, hypothesis, seed
#' @param posteriorSettings list with effectThresholds (efficacy/equivalence/futility) and alternative. The `equivalence` slot is supported by the API (P(lo < Delta < hi) > p_equ rule) but is not exercised by the ADRENAL re-design in the manuscript; users running a two-sided / equivalence-style design can pass it as a list(size = cbind(lo, hi), type = ...).
#' @param priorSettings two-component list returned by initialisePriorSettings (defaults to non-informative)
runTrialMonitoring <- function(simulationSettings,
                               posteriorSettings,
                               ...,
                               priorSettings = NULL) {
  # Suppress benign convergence warnings from RBesT/automixfit during the
  # inner kernel; restore the caller's warn level on exit (incl. error paths).
  old_warn <- getOption("warn")
  on.exit(options(warn = old_warn), add = TRUE)
  options(warn = -1)
  set.seed(seed = ifelse(is.null(simulationSettings$seed), 21,
                         simulationSettings$seed))
  if (is.null(priorSettings)) priorSettings <- initialisePriorSettings()

  numberOfStages <- length(simulationSettings$numberOfSubjects)
  numberOfSubjects <- list(
    control   = simulationSettings$numberOfSubjects -
      round(simulationSettings$numberOfSubjects / (1 + simulationSettings$allocationRatio)),
    treatment = round(simulationSettings$numberOfSubjects / (1 + simulationSettings$allocationRatio))
  )
  numberOfEvents <- list(control = NULL, treatment = NULL)
  posteriorProbabilities <- vector("list", length = numberOfStages)
  R <- simulationSettings$numberOfTrials

  # Per-stage progress messages are gated behind options(BayesGSD.verbose).
  # Default is silent because they fire on every stage for every call and
  # become noisy at routine R = 10^6 simulation budgets used in calibration.
  # Set options(BayesGSD.verbose = TRUE) for diagnostic runs.
  .verbose <- isTRUE(getOption("BayesGSD.verbose", FALSE))

  # Store cumulative events as R x stage matrices; vectorised rbinom across trials.
  ec_mat <- matrix(NA_integer_, nrow = R, ncol = numberOfStages)
  et_mat <- matrix(NA_integer_, nrow = R, ncol = numberOfStages)

  if (.verbose) message("calculating posterior probabilities at stage 1 ... ")
  ec_mat[, 1L] <- rbinom(
    n    = R,
    size = numberOfSubjects$control[1L],
    prob = simulationSettings$rates[1L]
  )
  et_mat[, 1L] <- rbinom(
    n    = R,
    size = numberOfSubjects$treatment[1L],
    prob = simulationSettings$rates[2L]
  )
  posteriorProbabilities[[1L]] <- .per_stage_posterior_probabilities(
    numberOfEvents   = cbind(ec_mat[, 1L], et_mat[, 1L]),
    numberOfSubjects = c(numberOfSubjects$control[1L], numberOfSubjects$treatment[1L]),
    posteriorSettings = posteriorSettings,
    priorSettings     = priorSettings
  )
  if (numberOfStages > 1L) {
    for (i in 2L:numberOfStages) {
      if (.verbose)
        message(sprintf("calculating posterior probabilities at stage %d ... ", i))
      ec_mat[, i] <- ec_mat[, i - 1L] + rbinom(
        n    = R,
        size = numberOfSubjects$control[i] - numberOfSubjects$control[i - 1L],
        prob = simulationSettings$rates[1L]
      )
      et_mat[, i] <- et_mat[, i - 1L] + rbinom(
        n    = R,
        size = numberOfSubjects$treatment[i] - numberOfSubjects$treatment[i - 1L],
        prob = simulationSettings$rates[2L]
      )
      posteriorProbabilities[[i]] <- .per_stage_posterior_probabilities(
        numberOfEvents   = cbind(ec_mat[, i], et_mat[, i]),
        numberOfSubjects = c(numberOfSubjects$control[i], numberOfSubjects$treatment[i]),
        posteriorSettings = posteriorSettings,
        priorSettings     = priorSettings
      )
    }
  }
  numberOfEvents$control   <- ec_mat
  numberOfEvents$treatment <- et_mat
  if (.verbose) message("Done.")

  list(
    virtualTrials = list(
      numberOfTrials         = R,
      numberOfStages         = numberOfStages,
      numberOfSubjects       = numberOfSubjects,
      numberOfEvents         = numberOfEvents,
      posteriorProbabilities = posteriorProbabilities,
      allocationRatio        = simulationSettings$allocationRatio,
      rates                  = simulationSettings$rates,
      hypothesis             = simulationSettings$hypothesis,
      seed                   = simulationSettings$seed
    ),
    priorSettings     = priorSettings,
    posteriorSettings = posteriorSettings
  )
}

# ============================================================================
# Operating characteristics
# ============================================================================

# --- Internal helper: stage-by-trial stop indicator matrices ---------------
.trigger_stopping <- function(posteriorProbabilities,
                              effectThresholds,
                              probabilityThresholds,
                              numberOfStages,
                              numberOfTrials) {
  out_eff <- NULL; out_equ <- NULL; out_fut <- NULL

  if (!is.null(effectThresholds$efficacy)) {
    stored <- attr(posteriorProbabilities[[1L]]$efficacy, "thresholds")
    idx <- .match_thresholds(effectThresholds$efficacy$size, stored,
                             what = "efficacy threshold")
    out_eff <- matrix(TRUE, numberOfStages, numberOfTrials)
    for (i in seq_len(numberOfStages)) {
      pp <- posteriorProbabilities[[i]]$efficacy
      for (j in seq_along(idx)) {
        out_eff[i, ] <- out_eff[i, ] &
          (pp[idx[j], ] > probabilityThresholds$efficacy[i, j])
      }
    }
  }
  if (!is.null(effectThresholds$equivalence)) {
    stored <- attr(posteriorProbabilities[[1L]]$equivalence, "thresholds")
    query <- cbind(as.numeric(effectThresholds$equivalence$size[, 1L]),
                   as.numeric(effectThresholds$equivalence$size[, 2L]))
    idx <- vapply(seq_len(nrow(query)), function(r) {
      j <- which(abs(stored[, 1L] - query[r, 1L]) <= 1e-8 &
                 abs(stored[, 2L] - query[r, 2L]) <= 1e-8)
      if (length(j) == 0L) NA_integer_ else j[1L]
    }, integer(1L))
    if (anyNA(idx)) stop("equivalence threshold pair not found in stored cache.")
    out_equ <- matrix(TRUE, numberOfStages, numberOfTrials)
    for (i in seq_len(numberOfStages)) {
      pp <- posteriorProbabilities[[i]]$equivalence
      for (j in seq_along(idx)) {
        out_equ[i, ] <- out_equ[i, ] &
          (pp[idx[j], ] > probabilityThresholds$equivalence[i, j])
      }
    }
  }
  if (!is.null(effectThresholds$futility)) {
    stored <- attr(posteriorProbabilities[[1L]]$futility, "thresholds")
    idx <- .match_thresholds(effectThresholds$futility$size, stored,
                             what = "futility threshold")
    out_fut <- matrix(TRUE, numberOfStages, numberOfTrials)
    for (i in seq_len(numberOfStages)) {
      pp <- posteriorProbabilities[[i]]$futility
      for (j in seq_along(idx)) {
        out_fut[i, ] <- out_fut[i, ] &
          (pp[idx[j], ] < probabilityThresholds$futility[i, j])
      }
    }
  }
  list(efficacyStop = out_eff, equivalenceStop = out_equ, futilityStop = out_fut)
}

# --- Internal helper: collapse per-stage stop matrices to per-trial outcome -
# Returns a list with `efficacyRates`, `equivalenceRates`, `futilityRates`,
# `inconclusiveRates` (all length-K cumulative counts), and
# `stoppingProbabilities` (length-K per-stage stop probs).
.collapse_stops <- function(effStop, equStop, futStop, numberOfStages, numberOfTrials) {
  has_eff <- !is.null(effStop)
  has_equ <- !is.null(equStop)
  has_fut <- !is.null(futStop)
  results <- rep(0L, numberOfTrials)
  stoppingProbabilities <- numeric(numberOfStages)
  inconclusiveRates <- integer(numberOfStages)
  efficacyRates     <- integer(numberOfStages)
  equivalenceRates  <- integer(numberOfStages)
  futilityRates     <- integer(numberOfStages)

  for (i in seq_len(numberOfStages)) {
    cur_eff <- if (has_eff) effStop[i, ] else logical(numberOfTrials)
    cur_equ <- if (has_equ) equStop[i, ] else logical(numberOfTrials)
    cur_fut <- if (has_fut) futStop[i, ] else logical(numberOfTrials)
    any_stop <- (results == 0L) & (cur_eff | cur_equ | cur_fut)
    stoppingProbabilities[i] <- sum(any_stop)
    if (has_eff) results[(results == 0L) & cur_eff] <- 1L
    if (has_equ) results[(results == 0L) & cur_equ] <- 2L
    if (has_fut) results[(results == 0L) & cur_fut] <- 3L
    inconclusiveRates[i] <- sum(results == 0L)
    efficacyRates[i]     <- sum(results == 1L)
    equivalenceRates[i]  <- sum(results == 2L)
    futilityRates[i]     <- sum(results == 3L)
  }
  stoppingProbabilities <- stoppingProbabilities / numberOfTrials
  # Match legacy NA-vs-0 conventions per active-rule combination.
  if (!has_eff) efficacyRates    <- rep(0L, numberOfStages)
  if (!has_equ) equivalenceRates <- rep(0L, numberOfStages)
  if (!has_fut) futilityRates    <- rep(0L, numberOfStages)
  list(
    stoppingProbabilities = stoppingProbabilities,
    inconclusiveRates     = inconclusiveRates,
    efficacyRates         = efficacyRates,
    equivalenceRates      = equivalenceRates,
    futilityRates         = futilityRates
  )
}

.expected_sample_size <- function(stoppingProbabilities, numberOfSubjects) {
  K <- length(numberOfSubjects)
  if (K == 1L) return(list(early = 0, expN = numberOfSubjects[1L]))
  # Cumulative early stopping over interim looks (stages 1 .. K-1).
  early <- sum(stoppingProbabilities[seq_len(K - 1L)])
  expN  <- sum(numberOfSubjects[seq_len(K - 1L)] *
                 stoppingProbabilities[seq_len(K - 1L)]) +
    numberOfSubjects[K] * (1 - early)
  list(early = early, expN = expN)
}

#' Summarise operating characteristics for a chosen trialDesign against a cached trialSimulation.
getOperatingCharacteristics <- function(trialSimulation,
                                        trialDesign) {
  numberOfStages <- length(trialDesign$numberOfSubjects)
  numberOfTrials <- trialSimulation$virtualTrials$numberOfTrials

  # Align the cached posteriorProbabilities to the design's chosen look times.
  cache_subj <- trialSimulation$virtualTrials$numberOfSubjects$control +
    trialSimulation$virtualTrials$numberOfSubjects$treatment
  active_stages <- which(cache_subj %in% trialDesign$numberOfSubjects)
  # Fail fast if the design's look times are not all represented in the cache:
  # otherwise active_stages would silently mis-align with the threshold matrices.
  # This holds for every design in this paper's grids (each design's looks are a
  # subset of the simulated union, per the Lemma 1 cache-invariance condition).
  if (length(active_stages) != length(trialDesign$numberOfSubjects)) {
    stop("getOperatingCharacteristics: the design's look times are not all ",
         "present in the cached simulation. Re-simulate at the union of all ",
         "candidate look times before evaluating this design.")
  }
  trialStatus <- .trigger_stopping(
    posteriorProbabilities =
      trialSimulation$virtualTrials$posteriorProbabilities[active_stages],
    effectThresholds      = trialDesign$effectThresholds,
    probabilityThresholds = trialDesign$probabilityThresholds,
    numberOfStages        = numberOfStages,
    numberOfTrials        = numberOfTrials
  )
  effStop <- trialStatus$efficacyStop
  equStop <- trialStatus$equivalenceStop
  futStop <- trialStatus$futilityStop

  if (is.null(effStop) && is.null(equStop) && is.null(futStop)) {
    return(NULL)
  }
  agg <- .collapse_stops(effStop, equStop, futStop, numberOfStages, numberOfTrials)
  esz <- .expected_sample_size(agg$stoppingProbabilities, trialDesign$numberOfSubjects)

  oc <- cbind(
    cumInfo              = trialDesign$numberOfSubjects,
    cumEfficacyRate      = agg$efficacyRates / numberOfTrials,
    cumEquivalenceRate   = agg$equivalenceRates / numberOfTrials,
    cumFutilityRate      = agg$futilityRates / numberOfTrials,
    cumInconclusiveRate  = agg$inconclusiveRates / numberOfTrials,
    stoppingProbPerStage = agg$stoppingProbabilities
  )
  out <- list(
    operatingCharacteristics = as.data.frame(oc),
    earlyStoppingProb        = esz$early,
    expSampleSize            = esz$expN,
    maxSampleSize            = max(trialDesign$numberOfSubjects)
  )
  # Type-stable schema: always return both type1ErrorRate and type2ErrorRate
  # slots, with the field that doesn't apply set to NA. This avoids forcing
  # every caller to branch on `hypothesis` before destructuring the result.
  final_eff_rate <- tail(agg$efficacyRates / numberOfTrials, n = 1L)
  if (trialSimulation$virtualTrials$hypothesis == "null") {
    out$type1ErrorRate <- final_eff_rate
    out$type2ErrorRate <- NA_real_
  } else {
    out$type1ErrorRate <- NA_real_
    out$type2ErrorRate <- 1 - final_eff_rate
  }
  out
}

# ============================================================================
# Single-stage interim wrapper (legacy API, used by older drivers)
# ============================================================================

runStageMonitoring <- function(trialData,
                               trialDesign,
                               ...,
                               priorSettings = NULL) {
  if (is.null(priorSettings)) priorSettings <- initialisePriorSettings()
  numberOfSubjects <- c(trialData$numberOfSubjects$control,
                        trialData$numberOfSubjects$treatment)
  numberOfEvents   <- c(trialData$numberOfEvents$control,
                        trialData$numberOfEvents$treatment)
  trialResult <- list(posteriorPDF = NULL, posteriorCDF = NULL,
                       posteriorProbabilities = NULL, earlyStop = NULL)

  # Domain for the posterior CDF/PDF support.
  is_absolute <-
    isTRUE(trialDesign$effectThresholds$efficacy$type == "absolute.risk") ||
    isTRUE(trialDesign$effectThresholds$futility$type == "absolute.risk")
  pdf_type <- trialDesign$effectThresholds$efficacy$type
  cdf_sizes <- if (is_absolute) seq(-1, 1, length.out = 1e4)
               else seq(0, 1e2, length.out = 1e4)

  cdf_probs <- getPosteriorProbabilities(
    treatmentEffects = list(size = cdf_sizes, type = pdf_type),
    numberOfEvents   = numberOfEvents,
    numberOfSubjects = numberOfSubjects,
    alternative      = trialDesign$alternative,
    priors           = priorSettings
  )
  cdf_sizes <- seq(
    from = max(cdf_sizes[cdf_probs < 1e-3]),
    to   = min(cdf_sizes[cdf_probs > 1 - 1e-3]),
    length.out = 1e4
  )
  cdf_probs <- getPosteriorProbabilities(
    treatmentEffects = list(size = cdf_sizes, type = pdf_type),
    numberOfEvents   = numberOfEvents,
    numberOfSubjects = numberOfSubjects,
    alternative      = trialDesign$alternative,
    priors           = priorSettings
  )
  trialResult$posteriorCDF <- list(type = pdf_type, sizes = cdf_sizes,
                                     probabilities = cdf_probs)
  trialResult$posteriorPDF <- list(
    type  = pdf_type,
    sizes = cdf_sizes,
    densities = getPosteriorDensities(
      treatmentEffects = list(size = cdf_sizes, type = pdf_type),
      numberOfEvents   = numberOfEvents,
      numberOfSubjects = numberOfSubjects,
      priors           = priorSettings
    )
  )

  for (kind in c("efficacy", "equivalence", "futility")) {
    thr <- trialDesign$effectThresholds[[kind]]
    pthr <- trialDesign$probabilityThresholds[[kind]]
    if (is.null(thr)) {
      trialResult$posteriorProbabilities[[kind]] <- NA
      trialResult$earlyStop[[kind]] <- FALSE
      next
    }
    p <- getPosteriorProbabilities(
      treatmentEffects = thr,
      numberOfEvents   = numberOfEvents,
      numberOfSubjects = numberOfSubjects,
      alternative      = trialDesign$alternative,
      priors           = priorSettings
    )
    trialResult$posteriorProbabilities[[kind]] <- p
    threshold_row <- pthr[trialData$stageNumber, ]
    trialResult$earlyStop[[kind]] <- if (kind == "futility") {
      all(p < threshold_row)
    } else {
      all(p > threshold_row)
    }
  }
  trialResult
}
