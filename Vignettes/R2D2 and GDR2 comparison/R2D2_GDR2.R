if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("yandorazhang/R2D2",force=T)
devtools::install_github("AugustArnstad/BayesianVariableImportance")
library(R2D2)
library(BayesianVariableImportance)
library(mvtnorm)
library(MCMCpack)
library(INLA)
library(relaimpo)
library(rstan)
library(MASS)
library(ggplot2)
library(dplyr)
library(tidyr)


## Please note that this example was created with the help of the R2D2 package from the Github repository of Yandora Zhang.
## Please note that this example was created with the help of the Stan files https://osf.io/ns2cv/, by Paul Bürkner and Javier Romero.
## The Stan files are used in this example and must be downloaded.


rinvgauss_new <- function(n, mean = 1, shape = NULL, dispersion = 1) {
  if (length(n) > 1)
    n <- length(n)
  if (n < 0)
    stop("n can't be negative")
  n <- as.integer(n)
  if (n == 0)
    return(numeric(0))
  if (!is.null(shape))
    dispersion <- 1/shape
  mu <- rep_len(mean, n)
  phi <- rep_len(dispersion, n)
  r <- rep_len(0, n)
  i <- (mu > 0 & phi > 0)
  if (!all(i)) {
    r[!i] <- NA
    n <- sum(i)
  }
  phi[i] <- phi[i] * mu[i]

  repeat {
    Y <- stats::rchisq(n, df = 1)
    X1 <- 1 + phi[i]/2 * (Y - sqrt(4 * Y/phi[i] + Y^2))
    X1[which(abs(X1) < 1e-06)] = 0  # 12/05/2014
    te = 1/(1 + X1)
    te1 = unlist(lapply(te, function(x) min(x, 1)))
    # te2 = unlist(lapply(te1, function(x) max(x,0)))
    firstroot <- as.logical(stats::rbinom(n, size = 1L, prob = te1))
    if (all(!is.na(firstroot)))
      break
  }

  r[i][firstroot] <- X1[firstroot]
  r[i][!firstroot] <- 1/X1[!firstroot]
  mu * r
}

r2d2marg <- function(x, y, hyper, mcmc.n = 10000, eps = 1e-07, thin = 1, print = TRUE) {

  #------------------------------------------
  EPS = 1e+20
  # 1. n, p, max(n,p)
  p <- ncol(x)
  n <- nrow(x)
  max.np <- max(n, p)

  # 2. hyperparameter values (a_pi, b, a1, b1)
  if (missing(hyper)) {
    b <- 0.5
    a_pi <- 1/(p^(b/2) * n^(b/2) * log(n))
    a1 <- 0.001
    b1 <- 0.001
  } else {
    a_pi <- hyper$a_pi
    b <- hyper$b
    a1 <- hyper$a1
    b1 <- hyper$b1
  }

  a <- a_pi * p

  # 3. Define variables to store posterior samples
  keep.beta = keep.phi = keep.psi = matrix(0, nrow = mcmc.n, ncol = p)
  keep.w = keep.xi = keep.sigma2 = rep(0, mcmc.n)

  # 4. Initial values.
  ## beta
  tem <- stats::coef(stats::lm(y ~ x - 1))
  tem[which(is.na(tem))] <- eps
  beta <- tem
  ## sigma
  sigma2 <- MCMCpack::rinvgamma(1, shape = a1 + n/2, scale = b1 + 0.5 * crossprod(y - x %*% beta))
  sigma <- sqrt(sigma2)
  ## phi
  phi <- rep(a_pi, p)
  ## xi
  xi <- a/(2 * b)
  ## w
  w <- b
  ## psi
  psi <- stats::rexp(p, rate = 0.5)
  #------------------------------------------

  XTX <- crossprod(x)  #t(X)%*%X
  XTY <- crossprod(x, y)  #t(X)%*%y

  keep.beta[1, ] <- beta
  keep.sigma2[1] <- sigma2
  keep.phi[1, ] <- phi
  keep.w[1] <- w
  keep.xi[1] <- xi
  keep.psi[1, ] <- psi

  #------------------------------------------
  #------------------------------------------
  # Now MCMC runnings!
  for (i in 2:mcmc.n) {
    for (j in 1:thin) {
      a <- a_pi * p

      # draw a number from N(0,1)
      Z <- stats::rnorm(p, mean = 0, sd = 1)

      #------------------------------------------
      # (i) Sample beta| phi, w, sigma2, y
      d.inv <- 1/(psi * phi * w)
      ad.inv <- abs(d.inv)
      sd.inv <- sign(d.inv)

      inx.e <- which(ad.inv < eps)
      inx.E <- which(ad.inv > EPS)
      d.inv[inx.e] <- eps * sd.inv[inx.e]
      d.inv[which(is.infinite((d.inv)))] <- EPS
      d.inv[inx.E] <- EPS * sd.inv[inx.E]

      if (length(d.inv) == 1) {
        Dinv <- d.inv
      } else {
        Dinv <- diag(d.inv)
      }

      Vinv <- XTX + Dinv

      # ********************************#********************************
      temQ <- chol(Vinv, pivot = T, tol = 1e-100000)
      pivot <- attr(temQ, "pivot")
      temc <- temQ[, order(pivot)]
      B <- XTY + t(temc) %*% Z * sigma


      # s = svd(Vinv) tes = s$u %*% diag(sqrt(s$d)) %*% t(s$v) b = XTY + t(tes)%*%Z*sqrt(sigma2[k-1])

      # b = XTY + t(chol(Vinv))%*%Z*sqrt(sigma2[k-1])
      # ********************************#********************************

      beta <- solve(Vinv, B, tol = 1e-10000)
      beta[which(is.na(beta))] <- 0

      #------------------------------------------
      # (ii) Sample sigma2| beta, phi, w, y
      sigma2 <- MCMCpack::rinvgamma(1, shape = a1 + n/2 + p/2, scale = b1 + 0.5 * crossprod(y - x %*%
                                                                                              beta) + sum(beta * d.inv * beta)/2)
      sigma <- sqrt(sigma2)

      #------------------------------------------
      # (iii) Sample psi| beta, phi, w, sigma2
      #------------------------------------------
      mu.te <- sigma * sqrt(phi * w)/abs(beta)
      amu.te <- abs(mu.te)
      smu.te <- sign(mu.te)
      inx.e <- which(amu.te < eps)
      inx.E <- which(amu.te > EPS)
      mu.te[inx.E] <- EPS * smu.te[inx.E]
      mu.te[inx.e] <- eps * smu.te[inx.e]
      ## rinvgauss_new: self-defined function
      psi <- 1/rinvgauss_new(p, mean = mu.te, shape = 1)


      #------------------------------------------
      # (iv) Sample w|beta, phi, xi, sigma2
      chi.te <- sum(beta^2/(psi * phi))/sigma2
      achi.te <- abs(chi.te)
      schi.te <- sign(chi.te)
      inx.e <- which(achi.te < eps)
      inx.E <- which(achi.te > EPS)
      chi.te[inx.e] <- eps * schi.te[inx.e]
      chi.te[inx.E] <- EPS * schi.te[inx.E]
      chi.te[which(chi.te == 0)] <- eps
      w <- GIGrvg::rgig(n = 1, lambda = a - p/2, chi = chi.te, psi = 4 * xi)



      #------------------------------------------
      # (v) Sample xi|w
      xi <- stats::rgamma(1, shape = a + b, rate = 1 + 2 * w)



      #------------------------------------------
      # (vi) Sample phi|beta, xi, sigma2, y
      TT <- rep(0, p)
      tem <- (beta^2/sigma2)/psi
      atem <- abs(tem)
      stem <- sign(tem)
      inx.e <- which(atem < eps)
      inx.E <- which(atem > EPS)
      tem[inx.e] <- eps * stem[inx.e]
      tem[inx.E] <- EPS * stem[inx.E]
      tem[which(tem == 0)] <- eps
      TT <- apply(matrix(tem, ncol = 1), 1, function(xx) return(GIGrvg::rgig(n = 1, lambda = a_pi -
                                                                               0.5, chi = xx, psi = 4 * xi)))

      Ts <- sum(TT)
      phi <- TT/Ts

    }

    if (print) {
      if (i%%500 == 0) {
        print(paste(c("The ", i, "th sample."), collapse = ""))
      }
    }

    keep.beta[i, ] <- beta
    keep.sigma2[i] <- sigma2
    keep.phi[i, ] <- phi
    keep.w[i] <- w
    keep.xi[i] <- xi
    keep.psi[i, ] <- psi

    if (sqrt(sum((keep.beta[i, ] - keep.beta[i - 1, ])^2)) < eps & abs(keep.xi[i] - keep.xi[i - 1]) <
        eps & sqrt(sum((keep.phi[i, ] - keep.phi[i - 1, ])^2)) < eps & abs(keep.w[i] - keep.w[i -
                                                                                              1]) < eps & abs(keep.sigma2[i] - keep.sigma2[i - 1]) < eps & sqrt(sum((keep.psi[i, ] - keep.psi[i -
                                                                                                                                                                                              1, ])^2)) < eps)
      break

  }

  return(list(beta = keep.beta, sigma2 = keep.sigma2, psi = keep.psi, w = keep.w, xi = keep.xi, phi = keep.phi))
}


# Data generation -----------------------------------------------------
generate_data <- function(rho, beta_true, n, p) {
  V <- matrix(c(1, rho, rho,
                rho, 1, rho,
                rho, rho, 1), nrow = p, ncol = p)
  x <- MASS::mvrnorm(n, rep(0, p), V)
  y <- x %*% beta_true + stats::rnorm(n)
  data <- data.frame(y, X1 = x[, 1], X2 = x[, 2], X3 = x[, 3])
  return(data)
}

# Define parameters
set.seed(1)
rho_values <- c(-0.4, -0.1, 0, 0.1, 0.4)
beta_true <- c(1, sqrt(2), sqrt(3))
p <- 3
n <- 1000
b <- 0.0575
a_pi <- 1/(p^(b/2) * n^(b/2) * log(n))
mu_R2 <- a_pi/(a_pi + b)
mu_R2
phi_R2 <- a_pi+b
phi_R2
hyper = list(b = b, a_pi = a_pi, a1 = 0.001, b1 = 0.001)

# R2D2 Fit -----------------------------------------------------------

r2d2_fits <- function(data, hyper){
  x <- as.matrix(data[, c("X1", "X2", "X3")])
  y <- data$y
  mcmc.n <- 10000
  burnIn <- 9000

  r2d2_fit <- r2d2marg(x = x, y = y, hyper=hyper, mcmc.n = mcmc.n, print = FALSE)

  beta_samples <- r2d2_fit$beta[burnIn:10000, ]


  psi_samples <- r2d2_fit$psi[burnIn:10000, ]


  phi_samples <- r2d2_fit$phi[burnIn:10000, ]


  r2_samples <- r2d2_fit$w[burnIn:10000] / (1 + r2d2_fit$w[burnIn:10000])

  #phi_samples <- phi_samples*r2_samples

  return(list(fit = r2d2_fit, beta = beta_samples, phi = phi_samples, psi = psi_samples, r2 = r2_samples))
}

data_high_neg <- generate_data(rho_values[1], beta_true, n, p)
data_low_neg <- generate_data(rho_values[2], beta_true, n, p)
data_no <- generate_data(rho_values[3], beta_true, n, p)
data_low_pos <- generate_data(rho_values[4], beta_true, n, p)
data_high_pos <- generate_data(rho_values[5], beta_true, n, p)

r2d2_high_neg <- r2d2_fits(data_high_neg, hyper)
r2d2_low_neg <- r2d2_fits(data_low_neg, hyper)
r2d2_no <- r2d2_fits(data_no, hyper)
r2d2_low_pos <- r2d2_fits(data_low_pos, hyper)
r2d2_high_pos <- r2d2_fits(data_high_pos, hyper)



# BVI ----------------------------------------------------------------
bvi_high_neg <- BayesianVariableImportance::perform_inla_analysis(data_high_neg, y ~ X1 + X2 + X3, family = "gaussian")
bvi_low_neg <- BayesianVariableImportance::perform_inla_analysis(data_low_neg, y ~ X1 + X2 + X3, family = "gaussian")
bvi_no <- BayesianVariableImportance::perform_inla_analysis(data_no, y ~ X1 + X2 + X3, family = "gaussian")
bvi_low_pos <- BayesianVariableImportance::perform_inla_analysis(data_low_pos, y ~ X1 + X2 + X3, family = "gaussian")
bvi_high_pos <- BayesianVariableImportance::perform_inla_analysis(data_high_pos, y ~ X1 + X2 + X3, family = "gaussian")


bvi_samples_high_neg <- BayesianVariableImportance::sample_posterior_gaussian(bvi_high_neg, y ~ X1 + X2 + X3, data_high_neg, n_samp=1000, additive_param = NULL)
bvi_samples_low_neg <- BayesianVariableImportance::sample_posterior_gaussian(bvi_low_neg, y ~ X1 + X2 + X3, data_low_neg, n_samp=1000, additive_param = NULL)
bvi_samples_no <- BayesianVariableImportance::sample_posterior_gaussian(bvi_no, y ~ X1 + X2 + X3, data_no, n_samp=1000, additive_param = NULL)
bvi_samples_low_pos <- BayesianVariableImportance::sample_posterior_gaussian(bvi_low_pos, y ~ X1 + X2 + X3, data_low_pos, n_samp=1000, additive_param = NULL)
bvi_samples_high_pos <- BayesianVariableImportance::sample_posterior_gaussian(bvi_high_pos, y ~ X1 + X2 + X3, data_high_pos, n_samp=1000, additive_param = NULL)

# GDR2 models -----------------------------------------------------------

# Define the function to fit models with GDR2 priors
fit_gdr2_models <- function(rho_values, beta_true, n, p, a_pi, iter = 10000, chains = 4, warmup=9000) {
  results <- list()

  # The Stan files are found from https://osf.io/ns2cv/, see article Bürkner et al. (2020) - Generalized Decomposition Priors on R2

  # Compile the Stan model once,
  stan_model_gdr2 <- stan_model(file = "path_to_Stan_file/logitR2D2.stan")

  # Define the mean vector and covariance matrix for logitphi

  mean <- digamma(a_pi) - digamma(a_pi)

  sigma_ii <- trigamma(a_pi) + trigamma(a_pi)
  sigma_ij <- trigamma(a_pi)
  mu_logitphi <- c(mean, mean)  # mean vector
  sqcov_logiphi <- matrix(c(sigma_ii, sigma_ij,
                            sigma_ij, sigma_ii), ncol= (p - 1), nrow=(p-1))  # identity matrix

  for (rho in rho_values) {
    # Generate data
    data <- generate_data(rho, beta_true, (n + 100), p)
    test_y <- data$y[1001:1100]
    test_x <- model.matrix(~ ., data = data[1001:1100, c("X1", "X2", "X3")])

    # Prepare Stan data
    stan_data_gdr2 <- list(
      N = 1000,
      y = data$y[1:1000],
      p = p + 1,  # including intercept
      X = model.matrix(~ ., data = data[1:1000, c("X1", "X2", "X3")]),
      mu_logitphi = mu_logitphi,  # mean vector for logitphi
      sqcov_logiphi = sqcov_logiphi,  # square root of covariance matrix for logitphi
      Ntest = 100,
      ytest = test_y,
      Xtest = test_x,
      R2D2_mean_R2 = 6/7,  # example mean of the R2 prior
      R2D2_prec_R2 = 0.4024984,  # example precision of the R2 prior
      prior_only = 0
    )

    # Fit the model
    fit <- sampling(stan_model_gdr2, data = stan_data_gdr2, iter = iter, chains = chains, warmup = warmup)

    extracted_samples <- rstan::extract(fit, pars = c("R2D2_phi", "R2D2_R2"))

    # Store the extracted samples in the results list
    results[[as.character(rho)]] <- extracted_samples. # Only store what is needed to avoid large objects

  }

  return(results)
}

# Define parameters
n <- 1000
p <- 3
rho_values <- c(-0.4, -0.1, 0, 0.1, 0.4)
beta_true <- c(1, sqrt(2), sqrt(3))

a_pi = 0.7

# Fit the models
gdr2_results <- fit_gdr2_models(rho_values, beta_true, n, p, a_pi)
# We recommend one fit and then storing the results. It may be time consuming if one must fit the models each time.
setwd("path_to_store/R2D2_GDR2_results")
saveRDS(gdr2_results, "GDR2_priors_results_15.05.rds")
gdr2_results <- readRDS("GDR2_priors_results_15.05.rds")

posterior_samples_high_neg <- gdr2_results$`-0.4`
posterior_samples_low_neg <- gdr2_results$`-0.1`
posterior_samples_no <- gdr2_results$`0`
posterior_samples_low_pos <- gdr2_results$`0.1`
posterior_samples_high_pos <- gdr2_results$`0.4`

posterior_samples_high_neg <- gdr2_results$`-0.4`
posterior_samples_low_neg <- gdr2_results$`-0.1`
posterior_samples_no <- gdr2_results$`0`
posterior_samples_low_pos <- gdr2_results$`0.1`
posterior_samples_high_pos <- gdr2_results$`0.4`


# Relaimpo ----------------------------------------------------------------
lm_high_neg <- lm(y ~ X1 + X2 + X3, data = data_high_neg)
lm_low_neg <- lm(y ~ X1 + X2 + X3, data = data_low_neg)
lm_no <- lm(y ~ X1 + X2 + X3, data = data_no)
lm_low_pos <- lm(y ~ X1 + X2 + X3, data = data_low_pos)
lm_high_pos <- lm(y ~ X1 + X2 + X3, data = data_high_pos)

relaimpo_high_neg <- boot.relimp(lm_high_neg, b=1000, type="lmg")
relaimpo_low_neg <- boot.relimp(lm_low_neg, b=1000, type="lmg")
relaimpo_no <- boot.relimp(lm_no, b=1000, type="lmg")
relaimpo_low_pos <- boot.relimp(lm_low_pos, b=1000, type="lmg")
relaimpo_high_pos <- boot.relimp(lm_high_pos, b=1000, type="lmg")

relaimpo_high_neg <- booteval.relimp(relaimpo_high_neg)
relaimpo_low_neg <- booteval.relimp(relaimpo_low_neg)
relaimpo_no <- booteval.relimp(relaimpo_no)
relaimpo_low_pos <- booteval.relimp(relaimpo_low_pos)
relaimpo_high_pos <- booteval.relimp(relaimpo_high_pos)


# Importance plots ----------------------------------------------------------------


# Define a function to extract data from the fits. Note that we scale the R2D2 and GDR2 results by the true R2 for comparison.
extract_r2d2_results <- function(fit, rho) {
  phi = fit$phi
  phi_samples <- as.data.frame(fit$phi)
  colnames(phi_samples) <- c("X1", "X2", "X3")
  phi_samples <- phi_samples %>%
    pivot_longer(cols = starts_with("X"), names_to = "Variable", values_to = "Importance") %>%
    mutate(Method = "R2D2", Correlation = rho)
  return(phi_samples)
}

extract_bvi_results <- function(bvi_sample, rho) {
  bvi_samples <- as.data.frame(bvi_sample$scaled_importance_samples)
  colnames(bvi_samples) <- c("X1", "X2", "X3")
  bvi_samples <- bvi_samples %>%
    pivot_longer(cols = starts_with("X"), names_to = "Variable", values_to = "Importance") %>%
    mutate(Method = "BVI", Correlation = rho)
  return(bvi_samples)
}

extract_gdr_results <- function(gdr_sample, rho) {
  gdr_sample <- gdr_sample
  gdr_samples <- as.data.frame(gdr_sample)
  colnames(gdr_samples) <- c("X1", "X2", "X3")
  gdr_samples <- gdr_samples %>%
    pivot_longer(cols = starts_with("X"), names_to = "Variable", values_to = "Importance") %>%
    mutate(Method = "GDR2", Correlation = rho)
  return(gdr_samples)
}

# Extract results from R2D2 fits and BVI samples
violin_data <- bind_rows(
  extract_r2d2_results(r2d2_high_neg, rho_values[1]),
  extract_r2d2_results(r2d2_low_neg, rho_values[2]),
  extract_r2d2_results(r2d2_no, rho_values[3]),
  extract_r2d2_results(r2d2_low_pos, rho_values[4]),
  extract_r2d2_results(r2d2_high_pos, rho_values[5]),
  extract_bvi_results(bvi_samples_high_neg, rho_values[1]),
  extract_bvi_results(bvi_samples_low_neg, rho_values[2]),
  extract_bvi_results(bvi_samples_no, rho_values[3]),
  extract_bvi_results(bvi_samples_low_pos, rho_values[4]),
  extract_bvi_results(bvi_samples_high_pos, rho_values[5]),
  extract_gdr_results(posterior_samples_high_neg$R2D2_phi, rho_values[1]),
  extract_gdr_results(posterior_samples_low_neg$R2D2_phi, rho_values[2]),
  extract_gdr_results(posterior_samples_no$R2D2_phi, rho_values[3]),
  extract_gdr_results(posterior_samples_low_pos$R2D2_phi, rho_values[4]),
  extract_gdr_results(posterior_samples_high_pos$R2D2_phi, rho_values[5])
)


# Extract LMG results and confidence intervals
extract_lmg_results <- function(lmg_boot, rho) {
  lmg_values <- data.frame(
    Variable = c("X1", "X2", "X3"),
    LMG = lmg_boot@lmg,
    Lower = c(lmg_boot@lmg.lower[1], lmg_boot@lmg.lower[2], lmg_boot@lmg.lower[3]),
    Upper = c(lmg_boot@lmg.upper[1], lmg_boot@lmg.upper[2], lmg_boot@lmg.upper[3]),
    Method = "LMG",
    Correlation = rho
  )
  return(lmg_values)
}

lmg_data <- bind_rows(
  extract_lmg_results(relaimpo_high_neg, rho_values[1]),
  extract_lmg_results(relaimpo_low_neg, rho_values[2]),
  extract_lmg_results(relaimpo_no, rho_values[3]),
  extract_lmg_results(relaimpo_low_pos, rho_values[4]),
  extract_lmg_results(relaimpo_high_pos, rho_values[5])
)

# Plotting with boxplots
r2d2_bvi_boxplot <- ggplot(violin_data, aes(x = factor(Correlation), y = Importance, fill = Method)) +
  geom_boxplot(alpha = 0.6, outlier.color = NA) +
  coord_cartesian(ylim=c(0, 0.7)) +
  geom_point(data = lmg_data, aes(x = factor(Correlation), y = LMG, shape = "LMG", color = "LMG"), size = 4) +
  #geom_point(data = lmg_data, aes(x = factor(Correlation), y = Lower, shape = "LMG Lower Bound", color = "LMG Lower Bound"), size = 2) +
  #geom_point(data = lmg_data, aes(x = factor(Correlation), y = Upper, shape = "LMG Upper Bound", color = "LMG Upper Bound"), size = 2) +
  facet_wrap(~Variable, scales = "free_y") +
  theme_minimal() +
  theme(
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size=12)
  )  +
  labs(#title = "Importance Distribution by Method and Correlation Level",
       x = "Correlation Level",
       y = "Importance") +
  scale_fill_manual(values = c("R2D2" = "#C6CDF7", "BVI" = "#E6A0C4", "GDR2"="#C6F7CD")) +
  #scale_shape_manual(name = "LMG", values = c("LMG Estimate" = 18, "LMG Lower Bound" = 24, "LMG Upper Bound" = 25)) +
  #scale_color_manual(name = "LMG", values = c("LMG Estimate" = "#32CD32", "LMG Lower Bound" = "grey", "LMG Upper Bound" = "grey")) +
  scale_shape_manual(name = "", values = c("LMG" = 18)) +
  scale_color_manual(name = "", values = c("LMG" = "#EF5350")) +
  guides(fill = guide_legend(override.aes = list(shape = NA)), shape = guide_legend(override.aes = list(fill = "grey")))

r2d2_bvi_boxplot


# R2 plots ---------------------------------------------------------------

# Extract R2 values from LMG results
extract_lmg_r2 <- function(lmg_boot, method_name, correlation) {
  data.frame(
    Method = method_name,
    Correlation = correlation,
    R2 = lmg_boot@R2,
    Lower = quantile(lmg_boot@R2.boot, 0.025),
    Upper = quantile(lmg_boot@R2.boot, 0.975)
  )
}

# Extract R2 values and confidence intervals for each correlation level
lmg_r2_data <- bind_rows(
  extract_lmg_r2(relaimpo_high_neg, "LMG", rho_values[1]),
  extract_lmg_r2(relaimpo_low_neg, "LMG", rho_values[2]),
  extract_lmg_r2(relaimpo_no, "LMG", rho_values[3]),
  extract_lmg_r2(relaimpo_low_pos, "LMG", rho_values[4]),
  extract_lmg_r2(relaimpo_high_pos, "LMG", rho_values[5])
)

# Combine R2D2 and BVI R2 samples
r2d2_r2_data <- bind_rows(
  data.frame(Method = "R2D2", Correlation = rho_values[1], R2 = r2d2_high_neg$r2),
  data.frame(Method = "R2D2", Correlation = rho_values[2], R2 = r2d2_low_neg$r2),
  data.frame(Method = "R2D2", Correlation = rho_values[3], R2 = r2d2_no$r2),
  data.frame(Method = "R2D2", Correlation = rho_values[4], R2 = r2d2_low_pos$r2),
  data.frame(Method = "R2D2", Correlation = rho_values[5], R2 = r2d2_high_pos$r2)
)

bvi_r2_data <- bind_rows(
  data.frame(Method = "BVI", Correlation = rho_values[1], R2 = unlist(bvi_samples_high_neg$R2_marginal, use.names = FALSE)),
  data.frame(Method = "BVI", Correlation = rho_values[2], R2 = unlist(bvi_samples_low_neg$R2_marginal, use.names = FALSE)),
  data.frame(Method = "BVI", Correlation = rho_values[3], R2 = unlist(bvi_samples_no$R2_marginal, use.names = FALSE)),
  data.frame(Method = "BVI", Correlation = rho_values[4], R2 = unlist(bvi_samples_low_pos$R2_marginal, use.names = FALSE)),
  data.frame(Method = "BVI", Correlation = rho_values[5], R2 = unlist(bvi_samples_high_pos$R2_marginal, use.names = FALSE))
)

gdr_r2_data <- bind_rows(
  data.frame(Method = "GDR2", Correlation = rho_values[1], R2 = posterior_samples_high_neg$R2D2_R2),
  data.frame(Method = "GDR2", Correlation = rho_values[2], R2 = posterior_samples_low_neg$R2D2_R2),
  data.frame(Method = "GDR2", Correlation = rho_values[3], R2 = posterior_samples_no$R2D2_R2),
  data.frame(Method = "GDR2", Correlation = rho_values[4], R2 = posterior_samples_low_pos$R2D2_R2),
  data.frame(Method = "GDR2", Correlation = rho_values[5], R2 = posterior_samples_high_pos$R2D2_R2)
)


# Combine all R2 data into one data frame
r2_data <- bind_rows(r2d2_r2_data, bvi_r2_data, gdr_r2_data)

# Define shapes for LMG points
shape <- c("LMG" = 18)
colors <- c("LMG" = "#EF5350")

# Plot R2 values
r2_plot <- ggplot(r2_data, aes(x = factor(Correlation), y = R2, fill = Method)) +
  geom_boxplot(alpha = 0.6, outlier.color = NA) +
  coord_cartesian(ylim=c(0.5, 1)) +
  geom_point(data = lmg_r2_data, aes(x = factor(Correlation), y = R2, shape = "LMG", color = "LMG"), size = 4) +
  theme_minimal() +
  theme(
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size=12)
  ) +
  labs(x = "Correlation Level", y = "R2 Value") +
  scale_shape_manual(name = "", values = shape) +
  scale_color_manual(name = "", values = colors) +
  scale_fill_manual(values = c("R2D2" = "#C6CDF7", "BVI" = "#E6A0C4", "GDR2" = "#C6F7CD"))

# Display the plot
r2_plot



