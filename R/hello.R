# The method can handle factor variables with more than 2 levels, as long as these have distinct names,
# and are encoded as dummy variables.


# Ensure the INLA package is loaded
if (!requireNamespace("INLA", quietly = TRUE)) {
  install.packages("INLA", repos = c(INLA = "https://inla.r-inla-download.org/R/stable"), type = "source")
}
library(INLA)



#' Singular Value Decomposition (SVD) and Scaling
#'
#' This function performs a singular value decomposition (SVD) on a numeric matrix and scales the data.
#' Non-numeric columns are converted to numeric values by treating factors and characters appropriately.
#'
#' @param X A data frame or matrix with numeric and possibly non-numeric columns.
#' @return A list containing the SVD-transformed matrix `Z`, the rotation matrix `R`,
#'         and the scaling factors `lambda`.
#' @examples
#' data(mtcars)
#' result <- SVD_decomp(mtcars)
#' @export
SVD_decomp <- function(X){

  for (col in colnames(X)) {
    if (!is.numeric(X[[col]])) {
      # If column is factor, convert to numeric by converting to integer codes
      # If column is character, first convert to factor then to numeric
      X[[col]] <- as.numeric(as.factor(X[[col]]))
    }
  }


  original_colnames <- colnames(X)
  X <- scale(X)

  # Calculate eigenvalues and eigenvectors
  e <- eigen(t(X) %*% X)
  Q <- e$vectors
  D <- diag(sqrt(e$values), nrow = nrow(Q))
  Dinv <- diag(1/sqrt(e$values), nrow = nrow(Q))

  # Calculate R_xx^(-0.5)
  R <- sqrt(nrow(X) - 1) * Q %*% Dinv %*% t(Q)

  # Calculate the transformed numerical fixed effects
  Z <- X %*% R
  colnames(Z) <- original_colnames <- colnames(X)


  lambda <- 1 / sqrt(nrow(X) - 1) * Q %*% D %*% t(Q)
  return(list(Z = Z, R = R, lambda = lambda))
}


#' Extract Fixed and Random Effects from a Model Formula
#'
#' This function identifies and separates fixed and random effects specified in a model formula.
#'
#' @param formula A model formula.
#' @return A list with two elements: `fixed_effects` and `random_effects`,
#'         each containing the names of the respective effects extracted from the formula.
#' @examples
#' formula <- y ~ x + f(group, model = "iid")
#' effects <- extract_effects(formula)
#' @export
extract_effects <- function(formula) {
  # Get the terms object from the formula
  formula_terms <- terms(formula)

  # Extract term labels
  term_labels <- attr(formula_terms, "term.labels")

  # Initialize lists to store fixed and random effects
  fixed_effects <- character()
  random_effects <- character()

  # Iterate through term labels to categorize them
  for (label in term_labels) {
    if (grepl("^f\\(", label)) {
      # If label starts with 'f(', it's a random effect
      random_effects <- c(random_effects, label)
    } else {
      # Otherwise, it's a fixed effect
      fixed_effects <- c(fixed_effects, label)
    }
  }

  if (length(fixed_effects) == 0){
    fixed_effects=NULL
  }
  if (length(random_effects) == 0){
    random_effects=NULL
  }

  # Return a list containing both fixed and random effects
  return(list(fixed_effects = fixed_effects, random_effects = random_effects))
}


#' Perform INLA Analysis on Specified Data and Formula
#'
#' Fits a model using the INLA approach, accommodating both fixed and random effects as specified.
#' The function handles preprocessing such as SVD decomposition and scaling.
#'
#' @param data The dataset to be analyzed.
#' @param formula The formula specifying the model to be fitted.
#' @param family The model family (e.g., "binomial", "poisson").
#' @param priors Optional list of prior specifications for model parameters.
#' @param link_func A link function to be used in the model (e.g., "identity", "logit").
#' @return An object of class `inla` representing the fitted model.
#' @examples
#' data <- data_binomial # Assuming data_binomial is available
#' formula <- y ~ x + f(group, model = "iid")
#' result <- perform_inla_analysis(data, formula, family = "binomial")
#' @export
perform_inla_analysis <- function(data, formula, family, link_func="identity", priors = NULL) {

  data_copy <- data

  response <- all.vars(formula)[1]


  # Set default priors if none are specified
  if (is.null(priors) && !family %in% c("binomial", "poisson")) {
    priors <- list(
      prec = list(
        prior = "pc.prec",
        param = c(1, 0.01),
        initial = log(1)
      )
    )
  }

  #formula_str <- deparse(formula)

  #simplified_formula <- gsub("f\\([^()]*\\)", "", formula_str)
  #simplified_vars <- strsplit(simplified_formula, "[ ~+]+")[[1]]
  #simplified_vars <- simplified_vars[2:(length(simplified_vars)-1)]

  effects <- extract_effects(formula)
  fixed_effects <- effects$fixed_effects

  #   if(!is.null(fixed_effects)){
  if(length(fixed_effects)>1){
    X <- data_copy[, c(fixed_effects)]

    SVD <- SVD_decomp(X)

    data_copy[, c(fixed_effects)] <- SVD$Z
  }

  scaled <- scale(data_copy[, response])
  #data_copy[, response] <- scaled

  # Process priors to include them in the formula if needed
  # This is a simplified example; actual implementation may need to adjust
  # based on how priors are specified and used in your formula

  # Fit the model using INLA
  inla_result <- inla(formula,
                      family = family,
                      data = data_copy,
                      control.family = list(hyper = priors, link = link_func),
                      control.compute = list(dic = FALSE, return.marginals=TRUE, config=TRUE, waic = TRUE))

  #inla_result$summary.fitted.values$descaled <- inla_result$summary.fitted.values$mean*attr(scaled, 'scaled:scale') + attr(scaled, 'scaled:center')

  # Return the INLA result object
  return(inla_result)
}

#' Extract Importance of Fixed and Random Effects
#'
#' Calculates the relative importance of both fixed and random effects from a fitted INLA model.
#'
#' @param model The fitted model object.
#' @param data The data used for fitting the model.
#' @param random_names A vector of names of the random effects.
#' @param fixed_names A vector of names of the fixed effects.
#' @param dist_factor Distributional variance factor, used to calculate residual variance. If NULL, the family and link will be checked. If the family is binomial with logit link, the dist_factor will be set to pi^2 / 3. If the family is binomial with probit link, the dist_factor will be set to 1. If the family is poisson, the dist_factor will be calculated as described in Nakagawa et. al - A general and simple method for obtaining R2 from generalized linear mixed-effects models (2013) and Nakagawa et. al - The coefficient of determination R2 and intra-class correlation coefficient from generalized linear mixed-effects models revisited and expanded (2017). If the family does not match any of these, the dist_factor will be set to 0, corresponding to Gaussian responses (LMM)
#' @return A list containing normalized importance values for random and fixed effects, marginal and conditional R2 values, and expected importances.
#' @examples
#' # Assuming `model` is an INLA model object and `data` contains appropriate predictors
#' importance <- extract_importances(model, data, dist_factor = pi^2 / 3, random_names = "Z1", fixed_names = c("X1", "X2", "X3"))
#' @export
extract_importances_test <- function(model, data, random_names, fixed_names, dist_factor=NULL) {
  # Decompose the fixed effects matrix using SVD
  SVD <- BayesianImpGLMM::SVD_decomp(data[, fixed_names])

  # Sample from the posterior distribution
  sample_post <- inla.posterior.sample(model, n = 1)

  # Extract latent effects names
  latent_names <- rownames(sample_post[[1]]$latent)

  # Initialize a named vector to store variances for each random effect
  imp_random <- setNames(numeric(length(random_names)), random_names)

  # Iterate through each random effect name and extract corresponding samples
  for (i in seq_along(random_names)) {
    pattern_random <- paste0("^", random_names[i], ":")
    random_indices <- grep(pattern_random, latent_names, value = TRUE)
    random_samples <- sample_post[[1]]$latent[random_indices, , drop = FALSE]

    # Variance of the current random effect
    imp_random[i] <- var(random_samples)
  }

  # Calculate residual variance
  residual_var <- 0
  if ("Precision for the Gaussian observations" %in% rownames(model$summary.hyperpar)) {
    residual_var <- 1 / model$summary.hyperpar["Precision for the Gaussian observations", "mean"]
  }

  # Initialize a named vector to store the fixed effects importances
  imp_fixed <- setNames(numeric(length(fixed_names)), fixed_names)

  # Extract and calculate importance for each fixed effect
  fixed_effects <- numeric(length(fixed_names))
  for (i in seq_along(fixed_names)) {
    fixed_indices <- grep(paste0("^", fixed_names[i], ":"), latent_names, value = TRUE)
    fixed_samples <- sample_post[[1]]$latent[fixed_indices, , drop = FALSE]

    # Store the extracted fixed effect sample
    fixed_effects[i] <- fixed_samples
  }

  # Calculate fixed effect importance using SVD
  imp_fixed <- (SVD$lambda^2 %*% (fixed_effects^2))

  fam <- model$.args$family

  link <- model$.args$control.family[[1]]$link

  if (is.null(dist_factor)){
    if (fam == "binomial"){
      if (link == "probit"){
        dist_factor <- 1
      } else if (link == "logit"){
        dist_factor <- pi^2/3

        # I think this could be difficult to implement. Ask Steffi.

        #intercept <- samps_Z[[i]]$latent[output_length-length(fixed)]
        #fixed_contribution <- as.matrix(data[, fixed]) %*% samps_Z[[i]]$latent[(samples_tot+2):output_length]

        #distribution_var <- inv.logit(intercept + fixed_contribution - 0.5*total_latent_var * tanh(((intercept+)*(1+2exp(-0.5*total_latent_var)))/6))
      }
    }else if (fam == "poisson"){
      if (link == "log"){
        intercept <- model$summary.fixed["(Intercept)", "mean"]
        lambda_pois <- exp(intercept + 0.5*sum(imp_random))
        dist_factor <- log(1 + 1/lambda_pois)
      }else if (link == "root"){
        dist_factor <- 0.25
      }
    }else{
      dist_factor <- 0
    }
  }


  # Calculate total variance
  total_var <- as.double(dist_factor + sum(imp_random) + sum(imp_fixed) + residual_var)

  # Calculate marginal and conditional R^2
  r2m <- sum(imp_fixed) / total_var
  r2c <- (sum(imp_fixed) + sum(imp_random)) / total_var

  residual_imp <- 0
  if ("Precision for the Gaussian observations" %in% rownames(model$summary.hyperpar)) {
    residual_imp <- 1 - r2c
  }else{
    residual_imp <- NA
  }

  # Normalize importances
  imp_random <- imp_random / total_var
  imp_fixed <- imp_fixed / total_var

  rownames(imp_fixed) <- fixed_names


  # Return the results as a list with named importances
  return(list(
    random_importance = imp_random,
    fixed_importance = imp_fixed,
    residual_importance = residual_imp,
    r2m = r2m,
    r2c = r2c
    #expected_importance = 1 / (dist_factor + sum(1:length(random_names)) + sum(1:length(fixed_names)) + residual_var + 1),
    #expected_r2m = sum(1:length(fixed_names)) / (dist_factor + sum(1:length(random_names)) + sum(1:length(fixed_names)) + residual_var + 1),
    #expected_r2c = (sum(1:length(fixed_names)) + sum(1:length(random_names)) + residual_var) / (dist_factor + sum(1:length(random_names)) + sum(1:length(fixed_names)) + residual_var + 1)
  ))
}

#' Sample Posterior Distributions for Gaussian Model Parameters
#'
#' Samples from the posterior distributions of model parameters and computes derived quantities for Gaussian models.
#'
#' @param model The fitted model object from INLA.
#' @param formula The model formula.
#' @param data The data used for fitting the model.
#' @param n_samp Number of samples to draw from the posterior distributions.
#' @param additive_param Optional; specifies an additive parameter for which to compute heritability.
#' @param param_of_interest Optional; specifies parameters of interest for detailed sampling.
#' @return A list containing matrices and data frames of sampled values and derived quantities.
#' @examples
#' # Assuming 'result' is an INLA model object and 'data_binomial' is available
#' samples <- sample_posterior_gaussian(result, formula, data_binomial, n_samp = 100)
#' @export
sample_posterior_gaussian <- function(model, formula, data, n_samp=1000, additive_param=NULL, param_of_interest=NULL) {

  # Make the param_of_interest a general input object for all functions. That was a nice way to put it


  response <- all.vars(formula)[1]
  scaled_response <- scale(data[, response])
  scale_const <- attr(scaled_response, 'scaled:scale')

  effects <- extract_effects(formula)
  fixed <- effects$fixed_effects

  response <- all.vars(formula)[1]

  variance_marginals_list <- lapply(model$marginals.hyperpar, function(x) inla.tmarginal(function(t) 1/t, x))

  hyperparam_names <- names(variance_marginals_list)

  random <- c()
  correlations <- c()

  for (name in hyperparam_names) {
    if (grepl("Precision", name)) {
      # If the name contains "Precision", it's a random effect
      random <- c(random, name)
    }else {
      # Otherwise, consider it a correlation parameter
      correlations <- c(correlations, name)
    }
  }


  #random <- names(variance_marginals_list)
  random <- gsub("Precision for the ", "", random)
  random <- gsub("Precision for ", "", random)
  # random <- names(variance_marginals_list)
  # random <- gsub("Precision for the ", "", random)
  # random <- gsub("Precision for ", "", random)

  sum <- summary(model)
  iidkd_indices <- which(sum$random.model == "IIDKD model")

  if (length(iidkd_indices) != 0) {
    iidkd_names <- sum$random.names[iidkd_indices]

    qq_matrices <- list()

    # Loop over each IIDKD model effect to perform sampling and calculate qq
    for (iidkd_effect in iidkd_names) {
      xx <- inla.iidkd.sample(500, model, iidkd_effect)
      order <- dim(xx[[1]])[1]
      qq <- matrix(rowMeans(matrix(unlist(xx), nrow = order^2)), order, order)
      qq_matrices[[iidkd_effect]] <- qq
    }
  }else{
    qq_matrices <- NULL
  }


  beta_mat <- matrix(NA, nrow=n_samp, ncol=length(fixed))
  scaled_beta_mat <- matrix(NA, nrow=n_samp, ncol=length(fixed))
  importance_mat <- matrix(NA, nrow=n_samp, ncol=length(fixed))
  scaled_importance_mat <- matrix(NA, nrow=n_samp, ncol=length(fixed))
  names(importance_mat) <- fixed
  R2_mat <- matrix(NA, nrow=n_samp, ncol=1)
  R2_cond_mat <- matrix(NA, nrow=n_samp, ncol=1)
  var_pred_mat <- matrix(NA, nrow=n_samp, ncol=1)
  h2_mat <- matrix(NA, nrow=n_samp, ncol=1)
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

    gaussian <- data[[response]][not_na] - predictor_samples[not_na]
    var_pred_mat[i] <- var(predictor_samples)
    random_mat[i, 1] <- var(gaussian)


    if (length(random)>1){
      for (j in 2:length(random)){
        pattern <- paste0("^", random[j], ":")
        random_names <- grep(pattern, latent_row_names, value = TRUE)
        random_samples <- samps_Z[[i]]$latent[random_names, , drop = FALSE]
        samples_tot <- samples_tot + length(random_samples)
        random_mat[i, j] <- var(random_samples)
      }
    }else{
      if (i==n_samp){
        print("No random effects, only resiudal variance")
      }

    }

    beta <- samps_Z[[i]]$latent[(output_length-(length(fixed)-1)):output_length]  #Skip intercept
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
  } else if (length(random)==2){
    R2_cond <- (rowSums(importance_mat) + random_mat[, -1]) / (rowSums(importance_mat) + rowSums(random_mat))
  }else{
    R2_cond <- rowSums(importance_mat) / (rowSums(importance_mat) + rowSums(random_mat))
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

  if (!is.null(additive_param)){
    #h2_mat <- random_mat[, additive_param]/(rowSums(importance_mat) + rowSums(random_mat))
    h2_mat <- random_mat[, additive_param]/(rowSums(random_mat))
    h2_mat <- as.data.frame(h2_mat)
    names(h2_mat) <- paste0("Heritability of: ", additive_param)
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
              heritability = h2_mat,
              iidkd_effects = qq_matrices))
}


#' Sample Posterior Distributions for Non-Gaussian Model Parameters
#'
#' Samples from the posterior distributions of model parameters and computes derived quantities for non-Gaussian models.
#'
#' @param model The fitted model object from INLA.
#' @param formula The model formula.
#' @param data The data used for fitting the model.
#' @param n_samp Number of samples to draw from the posterior distributions.
#' @param additive_param Optional; specifies an additive parameter for which to compute heritability.
#' @param param_of_interest Optional; specifies parameters of interest for detailed sampling.
#' @return A list containing matrices and data frames of sampled values and derived quantities.
#' @examples
#' # Assuming 'result' is an INLA model object and 'data_binomial' is available
#' samples <- sample_posterior_count(result, formula, data_binomial, n_samp = 100)
#' @export
sample_posterior_count <- function(model, formula, data, n_samp=1000, additive_param=NULL, param_of_interest=NULL) {

  # Make the param_of_interest a general input object for all functions. That was a nice way to put it

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
        # intercept <- samps_Z[[i]]$latent[output_length-length(fixed)]
        # distribution_var <- log((1/exp(intercept)) + 1)
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

  if (!is.null(additive_param)){
    #h2_mat <- random_mat[, additive_param]/(rowSums(importance_mat) + rowSums(random_mat) + distribution_var)
    h2_mat <- random_mat[, additive_param]/(rowSums(random_mat))
    h2_mat <- as.data.frame(h2_mat)
    names(h2_mat) <- paste0("Heritability of: ", additive_param)
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
              heritability = h2_mat))
}


#' Plot Samples from Posterior Distributions
#'
#' Generates plots for the scaled importance of fixed and random effects,
#' as well as marginal and conditional R2 values from sampled posterior distributions.
#'
#' @param samples A list object containing sampled values and derived quantities from `sample_posterior`.
#' @return A list of ggplot objects for each of the plotted quantities.
#' @examples
#' # Assuming 'samples' contains results from `sample_posterior`
#' plots <- plot_samples(samples)
#' @export
plot_samples <- function(samples) {
  plots <- list() # Initialize an empty list to store plots

  # Check and plot scaled importance for fixed effects
  if (!is.null(samples$scaled_importance_samples) && dim(samples$scaled_importance_samples)[2] > 0) {
    melted_scaled_importance <- melt(as.data.frame(samples$scaled_importance_samples))
    fixed_effects_plot <- ggplot(melted_scaled_importance, aes(x = value)) +
      geom_histogram(aes(y = ..density..), fill = "blue", alpha = 0.5) +
      geom_density(colour = "red", adjust = 1.5) +
      labs(title = "Relative Importance of Fixed Effects", x = "Value", y = "Density") +
      theme_minimal() +
      facet_wrap(~ variable, scales = "free_x")
    plots$fixed_effects <- fixed_effects_plot
  }

  # Check and plot scaled importance for random effects
  if (!is.null(samples$scaled_random_samples)) {
    melted_scaled_random <- melt(as.data.frame(samples$scaled_random_samples))
    random_effects_plot <- ggplot(melted_scaled_random, aes(x = value)) +
      geom_histogram(aes(y = ..density..), fill = "green", alpha = 0.5) +
      geom_density(colour = "darkgreen", adjust = 1.5) +
      labs(title = "Relative Importance of Random Effects", x = "Value", y = "Density") +
      theme_minimal() +
      facet_wrap(~ variable, scales = "free_x")
    plots$random_effects <- random_effects_plot
  }

  # Adjusted handling for heritability if it has values and not all are NA
  if (!is.null(samples$heritability) && any(!is.na(samples$heritability))) {
    # Extract the actual column name for heritability
    heritability_colname <- names(samples$heritability)

    # Dynamically specify the column name in aes() using rlang's sym() and !! for tidy evaluation
    heritability_plot <- ggplot(samples$heritability, aes(x = !!sym(heritability_colname))) +
      geom_histogram(aes(y = ..density..), fill = "purple", alpha = 0.5) +
      geom_density(color = "purple", adjust = 1.5) +
      labs(title = paste("Heritability of:", heritability_colname), x = heritability_colname, y = "Density") +
      theme_minimal()

    plots$heritability <- heritability_plot
  }

  if (!is.null(samples$R2_marginal) && !is.null(samples$R2_conditional)) {
    # Combine R2_marginal and R2_conditional into a single data frame
    R2_data <- data.frame(R2_marginal = samples$R2_marginal$`Marginal R2`,
                          R2_conditional = samples$R2_conditional$`Conditional R2`)

    # Convert from wide to long format
    R2_long <- pivot_longer(R2_data, cols = c(R2_marginal, R2_conditional), names_to = "Type", values_to = "Value")

    R2_long$Type <- factor(R2_long$Type, levels = c("R2_marginal", "R2_conditional"))

    # Plot
    R2_plot <- ggplot(R2_long, aes(x = Value, fill = Type)) +
      geom_histogram(aes(y = ..density..), alpha = 0.5, position = "identity") +
      geom_density(alpha = 0.75, adjust = 1.5) +
      labs(title = "Marginal and Conditional R2", x = "R2 Value", y = "Density") +
      scale_fill_manual(values = c("R2_marginal" = "blue", "R2_conditional" = "green")) +
      theme_minimal() +
      facet_wrap(~ Type, scales = "free_x")

    plots$R2 <- R2_plot
  }

  return(plots)
}


#' Plot Posterior Distributions of Variance Components and Heritability
#'
#' Generates plots for the posterior distributions of variance components and, optionally, heritability.
#'
#' @param model An `inla` model object.
#' @param random_effect_name Optional name of the random effect for which heritability is to be plotted.
#' @return A list of plot objects: `variance_plot` for variance components and `heritability_plot` for heritability.
#' @examples
#' # Assuming 'result' is an inla model object
#' plots <- plot_posteriors_and_heritability(result, random_effect_name = "group")
#' @export
plot_posteriors_and_heritability <- function(model, random_effect_name = NULL) {
  # Calculate marginals for variance components
  variance_marginals_list <- lapply(model$marginals.hyperpar, function(x) inla.tmarginal(function(t) 1/t, x))

  # Create data frames for plotting
  df_list <- lapply(names(variance_marginals_list), function(effect_name) {
    data.frame(
      x = variance_marginals_list[[effect_name]][, 1],
      y = variance_marginals_list[[effect_name]][, 2],
      Effect = effect_name
    )
  })
  df_combined <- do.call(rbind, df_list)

  df_combined$Effect <- gsub("Precision for ", "", df_combined$Effect, fixed = TRUE)

  # Plot posterior distributions
  p1 <- ggplot(df_combined, aes(x = x, y = y, color = Effect)) +
    geom_line() +
    labs(title = "Posterior Distributions of Variance Components",
         x = "Variance", y = "Density") +
    theme_minimal() +
    theme(legend.position = "right")

  # Plot heritability if random effect name is provided and correctly calculate heritability
  if (!is.null(random_effect_name) && random_effect_name %in% names(variance_marginals_list)) {
    # Calculate heritability distribution
    names <- names(model$marginals.hyperpar)
    total <- rep(0, length(variance_marginals_list[[random_effect_name]]))
    for (name in names){
      total <- total + variance_marginals_list[[name]][, 1]
    }

    heritability <- variance_marginals_list[[random_effect_name]][, 1]/total

    heritability_density <- variance_marginals_list[[random_effect_name]][ ,2]

    df_heritability <- data.frame(Heritability = heritability, Density = heritability_density)

    p2 <- ggplot(df_heritability, aes(x = Heritability, y = Density)) +
      geom_line() +
      labs(title = paste("Heritability Distribution for", random_effect_name),
           x = "Heritability", y = "Density") +
      theme_minimal()
  } else {
    p2 <- NULL
  }


  return(list(variance_plot = p1, heritability_plot = p2))
}


