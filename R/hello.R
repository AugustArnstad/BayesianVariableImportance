# The method cannot handle factor variables with more than 2 levels, as this creates two or more dummy variables,
# which causes a problem with the relative weights method. I could try to bypass this by handling each factor variable
# if it is worth the work.


# Ensure the INLA package is loaded
if (!requireNamespace("INLA", quietly = TRUE)) {
  install.packages("INLA", repos = c(INLA = "https://inla.r-inla-download.org/R/stable"), type = "source")
}
library(INLA)


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

# Function to perform analysis with INLA
perform_inla_analysis <- function(data, formula, family, priors = NULL) {

  data_copy <- data

  response <- all.vars(formula)[1]


  # Set default priors if none are specified
  if (is.null(priors)) {
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

  if(!is.null(fixed_effects)){
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
                      control.family = list(hyper = priors),
                      control.compute = list(dic = FALSE, return.marginals=TRUE, config=TRUE, waic = TRUE))

  #inla_result$summary.fitted.values$descaled <- inla_result$summary.fitted.values$mean*attr(scaled, 'scaled:scale') + attr(scaled, 'scaled:center')

  # Return the INLA result object
  return(inla_result)
}

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

sample_posterior <- function(model, formula, data, n_samp=1000, additive_param=NULL, param_of_interest=NULL, distribution_var=0, kuk=NULL) {

  # Make sure its correct with distributional variance
  # Make the param_of_interest a general input object for all functions. That was a nice way to put it
  # Scaling not complete


  response <- all.vars(formula)[1]
  scaled_response <- scale(data[, response])
  scale_const <- attr(scaled_response, 'scaled:scale')

  #fixed <- model$names.fixed[2:length(model$names.fixed)]
  effects <- extract_effects(formula)
  fixed <- effects$fixed_effects



  response <- all.vars(formula)[1]

  variance_marginals_list <- lapply(model$marginals.hyperpar, function(x) inla.tmarginal(function(t) 1/t, x))
  random <- names(variance_marginals_list)
  random <- gsub("Precision for the ", "", random)
  random <- gsub("Precision for ", "", random)

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
    #print(predictor_samples)
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

    beta <- samps_Z[[i]]$latent[(samples_tot+2):output_length]  #Skip intercept
    beta_mat[i, ] <- beta
    scaled_beta <- beta/scale_const
    scaled_beta_mat[i, ] <- scaled_beta
    importance_mat[i, ] <- lambda^2 %*% beta^2
    #scaled_importance_mat[i, ] <- lambda^2 %*% scaled_beta^2


  }

  rowsum <- rowSums(random_mat) + rowSums(importance_mat)
  scaled_random_mat <- random_mat/rowsum
  scaled_importance_mat <- importance_mat/rowsum

  # Do not think these are correct!!!!! They do not contain the guassian observations for example.
  # Could possibly use the variance of the predictor as the measure

  R2_mat <- rowSums(importance_mat) / (rowSums(importance_mat) + rowSums(random_mat) + distribution_var)

  if (length(random)>2){
    R2_cond <- (rowSums(importance_mat) + rowSums(random_mat[, -1]) + distribution_var) / (rowSums(importance_mat) + rowSums(random_mat) + distribution_var)
  } else if (length(random)==2){
    R2_cond <- (rowSums(importance_mat) + random_mat[, -1] + distribution_var) / (rowSums(importance_mat) + rowSums(random_mat) + distribution_var)
  }else{
    R2_cond <- rowSums(importance_mat) / (rowSums(importance_mat) + rowSums(random_mat) + distribution_var)
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
    h2_mat <- random_mat[, additive_param]/(rowSums(importance_mat) + rowSums(random_mat) + distribution_var)
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

plot_samples <- function(samples) {
  plots <- list() # Initialize an empty list to store plots

  # Check and plot scaled importance for fixed effects
  if (!is.null(samples$scaled_importance_samples)) {
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




