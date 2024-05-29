#' Extract Importance of Fixed and Random Effects
#'
#' Calculates the relative importance of both fixed and random effects from a fitted INLA model.
#'
#' @param model The fitted model object.
#' @param data The data used for fitting the model.
#' @param dist_factor Distributional variance factor, used to calculate residual variance. If NULL, the family and link will be checked. If the family is binomial with logit link, the dist_factor will be set to pi^2 / 3. If the family is binomial with probit link, the dist_factor will be set to 1. If the family is poisson, the dist_factor will be calculated as described in Nakagawa et. al - A general and simple method for obtaining R2 from generalized linear mixed-effects models (2013) and Nakagawa et. al - The coefficient of determination R2 and intra-class correlation coefficient from generalized linear mixed-effects models revisited and expanded (2017). If the family does not match any of these, the dist_factor will be set to 0, corresponding to Gaussian responses (LMM)
#' @return A list containing normalized importance values for random and fixed effects, marginal and conditional R2 values, and expected importances.
#' @examples
#' # Assuming `model` is an INLA model object with log, logit or probit link and `data` contains appropriate predictors
#' importance <- extract_importances(model, data, dist_factor = pi^2 / 3)
#' @export
extract_importances <- function(model, data, dist_factor=NULL) {

  # Calculate residual variance
  residual_var <- 0

  if ("Precision for the Gaussian observations" %in% rownames(model$summary.hyperpar)) {
    residual_var <- 1 / model$summary.hyperpar["Precision for the Gaussian observations", "mean"]
    randoms <- summary(model)$hyperpar[, c("mode")]
    random_effects <- randoms[-1]
    random_names <- rownames(model$summary.hyperpar)[-1]
  }else{
    randoms <- summary(model)$hyperpar[, c("mode")]
    random_names <- rownames(model$summary.hyperpar)
  }

  imp_random <- setNames(numeric(length(random_names)), random_names)
  for (i in seq_along(random_names)) {
    prec_random <- summary(model)$hyperpar[c(i), c("mode")]
    imp_random[i] <- 1/prec_random
  }

  random_names <- gsub("Precision for ", "", random_names)

  fixed_means <- summary(model)$fixed[, c("mode")]
  fixed_effects <- fixed_means[-1]
  fixed_names <- names(fixed_effects)

  SVD <- BayesianVariableImportance::SVD_decomp(data[, fixed_names])

  imp_fixed <- setNames(numeric(length(fixed_effects)), fixed_names)

  imp_fixed <- (SVD$lambda^2 %*% (fixed_effects^2))

  fam <- model$.args$family

  link <- model$.args$control.family[[1]]$link

  if (is.null(dist_factor)){
    if (fam == "binomial"){
      if (link == "probit"){
        dist_factor <- 1
      } else if (link == "logit"){
        dist_factor <- pi^2/3

      }
    }else if (fam == "poisson"){
      if (link == "log"){
        intercept <- model$summary.fixed["(Intercept)", "mean"]
        lambda_pois <- exp(intercept + 0.5*(sum(imp_random) + sum(imp_fixed)))
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

  #rownames(imp_fixed) <- fixed_names


  # Return the results as a list with named importances
  return(list(
    random_importance = imp_random,
    fixed_importance = imp_fixed,
    residual_importance = residual_imp,
    r2m = r2m,
    r2c = r2c
  ))
}
