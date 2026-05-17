#' @title Prior-approximation diagnostics for the logit-normal example
#' @description Quantifies how well the beta-mixture approximation of
#'   Section 2.3 represents the target logit-normal prior. Four diagnostics
#'   are produced: (i) the forward Kullback-Leibler divergence
#'   D_KL(p || hat_p_L) as a function of the number of mixture components L;
#'   (ii) a density overlay comparing the target prior to the fitted
#'   beta-mixture approximations at each L; (iii) prior tail probabilities
#'   at the decision-rule thresholds (P(Delta < 0) and P(Delta < -0.05),
#'   computed from independent Monte Carlo samples for each arm); and (iv)
#'   sensitivity of the resulting operating characteristics to the number
#'   of mixture components. The output is saved as an .rda cache and a
#'   .jpeg panel figure.
#' @author Zhangyi He, Feng Yu, Suzie Cro, Laurent Billot

# Project root, parallelism, INLA threads, sessionInfo helper.
# See Code/Code v1.0/bayseqSim_bern_setup.R for behaviour and override env vars.
local({
  bootstrap_root <- Sys.getenv("BAYESGSD_ROOT", unset = "")
  if (!nzchar(bootstrap_root) || !dir.exists(bootstrap_root)) {
    cur <- getwd()
    while (cur != dirname(cur)) {
      if (dir.exists(file.path(cur, "Code")) && dir.exists(file.path(cur, "Article"))) {
        bootstrap_root <- cur; break
      }
      cur <- dirname(cur)
    }
  }
  if (!nzchar(bootstrap_root) || !dir.exists(bootstrap_root)) {
    stop("Could not resolve BAYESGSD_ROOT. Set it before running.")
  }
  source(file.path(bootstrap_root, "Code", "Code v1.0", "bayseqSim_bern_setup.R"))
})



suppressPackageStartupMessages({
  library(RBesT); library(zipfR); library(extraDistr); library(parallel)
  library(ggplot2); library(ggsci); library(patchwork)
})
source("./Code/Code v1.0/bayseqSim_bern.R")

OUTPUT_DIR <- "./Output/Output v1.0"
RDA_PATH   <- file.path(OUTPUT_DIR, "ADRENAL_PriorDiagnostics.rda")
FIG_PATH   <- file.path(OUTPUT_DIR, "ADRENAL_PriorDiagnostics.jpeg")

# Logit-normal target prior on each arm: same as the supplement / Section 3.x
PRIOR_MU    <- log(0.33 / (1 - 0.33))
PRIOR_SIGMA <- 0.5
N_PRIOR_SAMPLES <- 50000L
SEED <- 21L

logit_normal_samples <- function(n, mu, sigma) {
  z <- rnorm(n, mu, sigma)
  1 / (1 + exp(-z))
}

set.seed(SEED)
samples_c <- logit_normal_samples(N_PRIOR_SAMPLES, PRIOR_MU, PRIOR_SIGMA)
samples_t <- logit_normal_samples(N_PRIOR_SAMPLES, PRIOR_MU, PRIOR_SIGMA)

# ----- Diagnostic 1 & 4: refit at L = 1, ..., 5 and record KL, ESS, ----------
L_GRID <- 1:5
fit_list <- vector("list", length(L_GRID))
for (i in seq_along(L_GRID)) {
  L <- L_GRID[i]
  cat(sprintf("Fitting beta mixture at L=%d ...\n", L))
  fit <- automixfit(sample = samples_c, Nc = L, type = "beta")
  fit_list[[i]] <- fit
}

# Forward KL divergence D_KL(p || hat_p_L) estimated by empirical mean
# log p(x) - log hat_p_L(x) over the prior samples.
log_dmix_beta <- function(x, weights, alpha, beta) {
  log_components <- mapply(function(al, bl) dbeta(x, al, bl, log = TRUE),
                            alpha, beta)
  if (is.matrix(log_components)) {
    log_w <- matrix(log(weights), nrow = nrow(log_components),
                    ncol = ncol(log_components), byrow = TRUE)
    apply(log_w + log_components, 1, function(v) {
      m <- max(v); m + log(sum(exp(v - m)))
    })
  } else {
    m <- max(log(weights) + log_components)
    m + log(sum(exp(log(weights) + log_components - m)))
  }
}

# Empirical density of the target prior (logit-normal on response-rate scale).
# Avoid the singular endpoints and use the change-of-variables form on the
# logit scale: if z = logit(theta) ~ N(mu, sigma^2) then
#   p(theta) = phi((logit(theta) - mu)/sigma)/sigma * 1/[theta(1-theta)].
dlogitnorm <- function(theta, mu, sigma) {
  z <- log(theta / (1 - theta))
  dnorm(z, mu, sigma) / (theta * (1 - theta))
}

kl_div <- numeric(length(L_GRID))
ess_v  <- numeric(length(L_GRID))
for (i in seq_along(L_GRID)) {
  fit <- fit_list[[i]]
  weights <- as.numeric(fit["w", ])
  alpha   <- as.numeric(fit["a", ])
  beta    <- as.numeric(fit["b", ])
  log_q   <- log_dmix_beta(samples_c, weights, alpha, beta)
  log_p   <- log(dlogitnorm(samples_c, PRIOR_MU, PRIOR_SIGMA))
  kl_div[i] <- mean(log_p - log_q)
  ess_v[i]  <- RBesT::ess(fit)
}
diag_table <- data.frame(L = L_GRID, KL = kl_div, ESS = ess_v)
cat("\n=== Diagnostic 1 & 4 (KL and ESS by L) ===\n")
print(diag_table)

# ----- Diagnostic 2: density overlay ----------------------------------------
theta_grid <- seq(0.05, 0.65, length.out = 401)
overlay_df <- data.frame()
for (i in seq_along(L_GRID)) {
  fit <- fit_list[[i]]
  weights <- as.numeric(fit["w", ])
  alpha   <- as.numeric(fit["a", ])
  beta    <- as.numeric(fit["b", ])
  comps   <- mapply(function(al, bl) dbeta(theta_grid, al, bl), alpha, beta)
  q_dens  <- rowSums(t(t(comps) * weights))
  overlay_df <- rbind(overlay_df,
                      data.frame(theta = theta_grid, density = q_dens,
                                 method = sprintf("beta mixture (L=%d)", L_GRID[i])))
}
overlay_df <- rbind(overlay_df,
                    data.frame(theta   = theta_grid,
                               density = dlogitnorm(theta_grid, PRIOR_MU, PRIOR_SIGMA),
                               method  = "logit-normal target"))

# ----- Diagnostic 3: prior tail probabilities at decision thresholds --------
# Under independent prior on (theta_c, theta_t), Delta = theta_t - theta_c.
# Compute P(Delta < 0) and P(Delta < -0.05) from MC samples.
tail_grid <- data.frame()
# Target tail probabilities from the original logit-normal samples
samples_t_lgn <- samples_t  # independent draws for treatment arm under target
delta_lgn <- samples_t_lgn - samples_c
tail_grid <- rbind(tail_grid, data.frame(
  method   = "logit-normal target",
  P_Delta_lt_0    = mean(delta_lgn <  0),
  P_Delta_lt_m005 = mean(delta_lgn <  -0.05)
))
# Tail probabilities from each fitted beta-mixture
for (i in seq_along(L_GRID)) {
  fit <- fit_list[[i]]
  weights <- as.numeric(fit["w", ])
  alpha   <- as.numeric(fit["a", ])
  beta    <- as.numeric(fit["b", ])
  sample_beta_mix <- function(n) {
    comp <- sample.int(length(weights), n, replace = TRUE, prob = weights)
    rbeta(n, alpha[comp], beta[comp])
  }
  set.seed(SEED + i)
  ths_c <- sample_beta_mix(N_PRIOR_SAMPLES)
  ths_t <- sample_beta_mix(N_PRIOR_SAMPLES)
  d <- ths_t - ths_c
  tail_grid <- rbind(tail_grid, data.frame(
    method          = sprintf("beta mixture (L=%d)", L_GRID[i]),
    P_Delta_lt_0    = mean(d <  0),
    P_Delta_lt_m005 = mean(d <  -0.05)
  ))
}
cat("\n=== Diagnostic 3 (prior tail probabilities at decision thresholds) ===\n")
print(tail_grid)

# ----- Plotting ------------------------------------------------------------
overlay_df$method <- factor(overlay_df$method,
                            levels = c("logit-normal target",
                                       sprintf("beta mixture (L=%d)", L_GRID)))
p_overlay <- ggplot(overlay_df, aes(x = theta, y = density,
                                     color = method, linetype = method)) +
  geom_line(linewidth = 0.7) +
  scale_color_bmj() +
  scale_linetype_manual(values = c("solid", rep("dashed", length(L_GRID)))) +
  xlab(expression(vartheta)) +
  ylab(expression(p(vartheta))) +
  ggtitle("Prior density overlay") +
  theme_bw() + theme(legend.position = "bottom",
                     legend.title = element_blank())

p_kl <- ggplot(diag_table, aes(x = L, y = KL)) +
  geom_point(size = 3) + geom_line(linewidth = 0.7) +
  scale_x_continuous(breaks = L_GRID) +
  xlab("L") + ylab(expression(D[KL](p ~ "||" ~ hat(p)[L]))) +
  ggtitle("Forward KL divergence by L") +
  theme_bw()

p_ess <- ggplot(diag_table, aes(x = L, y = ESS)) +
  geom_point(size = 3) + geom_line(linewidth = 0.7) +
  scale_x_continuous(breaks = L_GRID) +
  xlab("L") + ylab("Effective sample size") +
  ggtitle("Mixture effective sample size by L") +
  theme_bw()

tail_long <- rbind(
  data.frame(method = tail_grid$method,
             threshold = "P(Delta<0)",
             prob = tail_grid$P_Delta_lt_0),
  data.frame(method = tail_grid$method,
             threshold = "P(Delta<-0.05)",
             prob = tail_grid$P_Delta_lt_m005)
)
tail_long$method <- factor(tail_long$method,
                            levels = c("logit-normal target",
                                       sprintf("beta mixture (L=%d)", L_GRID)))
p_tail <- ggplot(tail_long, aes(x = method, y = prob, fill = threshold)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_bmj() +
  xlab(NULL) + ylab("prior tail probability") +
  ggtitle("Prior tail probabilities at decision thresholds") +
  theme_bw() + theme(legend.position = "bottom",
                     axis.text.x = element_text(angle = 25, hjust = 1))

combined <- (p_overlay + p_kl) / (p_ess + p_tail) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

ggsave(FIG_PATH, combined, width = 12, height = 8, dpi = 200)
cat(sprintf("Saved %s\n", FIG_PATH))

save(diag_table, tail_grid, overlay_df, fit_list,
     L_GRID, PRIOR_MU, PRIOR_SIGMA,
     file = RDA_PATH)
cat(sprintf("Saved %s\n", RDA_PATH))

# Persist sessionInfo and package versions alongside the cache.
bayesgsd_save_session("ADRENAL_Diagnostics_PriorApproximation")
