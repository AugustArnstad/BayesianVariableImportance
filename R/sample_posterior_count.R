#' Sample Posterior Distributions for Non-Gaussian Model Parameters
#'
#' Samples from the posterior distributions of model parameters and computes derived quantities for non-Gaussian models.
#'
#' @param model The fitted model object from INLA.
#' @param formula The model formula.
#' @param data The data used for fitting the model.
#' @param n_samp Number of samples to draw from the posterior distributions.
#' @param additive_param Optional; specifies an additive parameter for which to compute repeatability
#' @param repeatability Optional; calculates the repeatability (variance of additive_param divided by sum of random variances) according to Stoffel et. al - rptR: repeatability estimation and variance decomposition by generalized linear mixed-effects models (2017).
#' @return A list containing matrices and data frames of sampled values and derived quantities.
#' @examples
#' # Assuming 'result' is an INLA model object and 'data_binomial' is available
#' samples <- sample_posterior_count(result, formula, data_binomial, n_samp = 100)
#' @export
sample_posterior_count <- function(model, formula, data, n_samp=1000, additive_param=NULL, repeatability = FALSE) {


  response <- all.vars(formula)[1]
  scaled_response <- scale(data[, response])
  scale_const <- attr(scaled_response, 'scaled:scale')

  effects <- extract_effects(formula)
  fixed <- effects$fixed_effects

  fam <- model$.args$family

  link <- model$.args$control.family[[1]]$link

  distribution = paste0("Distributional variance: ", link)

  response <- all.vars(formula)[1]

  variance_marginals_list <- lapply(model$marginals.hyperpar, function(x) inla.tmarginal(function(t) 1/t, x))
  random <- names(variance_marginals_list)
  random <- gsub("Precision for the ", "", random)
  random <- gsub("Precision for ", "", random)

  random <- c(distribution, random)

  beta_mat <- matrix(NA, nrow=n_samp, ncol=length(fixed))
  scaled_beta_mat <- matrix(NA, nrow=n_samp, ncol=length(fixed))
  importance_mat <- matrix(NA, nrow=n_samp, ncol=length(fixed))
  scaled_importance_mat <- matrix(NA, nrow=n_samp, ncol=length(fixed))
  names(importance_mat) <- fixed
  R2_mat <- matrix(NA, nrow=n_samp, ncol=1)
  R2_cond_mat <- matrix(NA, nrow=n_samp, ncol=1)
  var_pred_mat <- matrix(NA, nrow=n_samp, ncol=1)
  h2_mat <- matrix(NA, nrow=n_samp, ncol=1)
  repeat_mat <- matrix(NA, nrow=n_samp, ncol=1)
  random_mat <- matrix(NA, nrow=n_samp, ncol=length(random))
  scaled_random_mat <- matrix(NA, nrow=n_samp, ncol=length(random))

  if(!is.null(fixed)){
    SVD <- VariableImportanceINLA::SVD_decomp(data[fixed])

    lambda <- SVD$lambda
  }else {
    lambda <- diag(length(fixed))
  }

  samps_Z <- inla.posterior.sample(model, n = n_samp)

  latent_row_names <- rownames(samps_Z[[1]]$latent)

  output_length=length(samps_Z[[1]]$latent)

  not_na = which(!is.na(data[response]))

  for (i in 1:n_samp){
    # Extract all sampled values, separate them by covariate/predictor, and assign them to the sampled matrix
    samples_length <- 0
    predictor <- paste0("^Predictor:")
    predictor_names <- grep(predictor, latent_row_names, value = TRUE)
    predictor_samples <- samps_Z[[i]]$latent[predictor_names, , drop = FALSE]
    samples_tot <- length(predictor_samples)


    total_latent_var <- 0
    if (length(random)>1){
      for (j in 2:length(random)){
        pattern <- paste0("^", random[j], ":")
        random_names <- grep(pattern, latent_row_names, value = TRUE)
        random_samples <- samps_Z[[i]]$latent[random_names, , drop = FALSE]
        samples_tot <- samples_tot + length(random_samples)
        random_mat[i, j] <- var(random_samples)

        total_latent_var <- total_latent_var + var(random_samples)
      }
    }else{
      if (i==n_samp){
        print("No random effects, only resiudal variance")
      }

    }

    if (fam == "binomial"){
      if (link == "probit"){
        distribution_var <- 1
      } else if (link == "logit"){
        distribution_var <- pi^2/3

        # I think this could be difficult to implement. Ask Steffi.

        #intercept <- samps_Z[[i]]$latent[output_length-length(fixed)]
        #fixed_contribution <- as.matrix(data[, fixed]) %*% samps_Z[[i]]$latent[(samples_tot+2):output_length]

        #distribution_var <- inv.logit(intercept + fixed_contribution - 0.5*total_latent_var * tanh(((intercept+)*(1+2exp(-0.5*total_latent_var)))/6))
      }
    }else if (fam == "poisson"){
      if (link == "log"){
        intercept <- samps_Z[[i]]$latent[output_length-length(fixed)]
        lambda_pois <- exp(intercept + 0.5*total_latent_var)
        distribution_var <- log(1 + 1/lambda_pois)
      }else if (link == "root"){
        distribution_var <- 0.25
      }
    }

    random_mat[i, 1] <- distribution_var

    beta <- samps_Z[[i]]$latent[(samples_tot+2):output_length]  #Skip intercept
    beta_mat[i, ] <- beta
    importance_mat[i, ] <- lambda^2 %*% beta^2
  }

  rowsum <- rowSums(random_mat) + rowSums(importance_mat)
  scaled_random_mat <- random_mat/rowsum
  scaled_beta_mat <- beta_mat/rowsum
  scaled_importance_mat <- importance_mat/rowsum

  R2_mat <- rowSums(importance_mat) / (rowSums(importance_mat) + rowSums(random_mat))

  if (length(random)>2){
    R2_cond <- (rowSums(importance_mat) + rowSums(random_mat[, -1])) / (rowSums(importance_mat) + rowSums(random_mat))
  }else if (length(random)==2){
    R2_cond <- (rowSums(importance_mat) + random_mat[, -1] ) / (rowSums(importance_mat) + rowSums(random_mat))
  }else{
    R2_cond <- rowSums(importance_mat)  / (rowSums(importance_mat) + rowSums(random_mat))
  }

  beta_mat <- as.data.frame(beta_mat)
  names(beta_mat) <- fixed
  importance_mat <- as.data.frame(importance_mat)
  names(importance_mat) <- fixed
  scaled_beta_mat <- as.data.frame(scaled_beta_mat)
  names(scaled_beta_mat) <- fixed
  scaled_importance_mat <- as.data.frame(scaled_importance_mat)
  names(scaled_importance_mat) <- fixed

  random_mat <- as.data.frame(random_mat)
  names(random_mat) <- random
  scaled_random_mat <- as.data.frame(scaled_random_mat)
  names(scaled_random_mat) <- random

  R2_mat <- as.data.frame(R2_mat)
  names(R2_mat) <- "Marginal R2"
  R2_cond <- as.data.frame(R2_cond)
  names(R2_cond) <- "Conditional R2"

  if (!is.null(additive_param) && repeatability){
    repeat_mat <- random_mat[, additive_param]/(rowSums(random_mat))
    repeat_mat <- as.data.frame(repeat_mat)
    names(repeat_mat) <- paste0("Repeatability of: ", additive_param)
  }

  return(list(beta_samples = beta_mat,
              importance_samples = importance_mat,
              scaled_beta_samples = scaled_beta_mat,
              scaled_importance_samples = scaled_importance_mat,
              random_samples = random_mat,
              scaled_random_samples = scaled_random_mat,
              R2_marginal = R2_mat,
              R2_conditional = R2_cond,
              var_y = var_pred_mat,
              repeatability = repeat_mat))
}
