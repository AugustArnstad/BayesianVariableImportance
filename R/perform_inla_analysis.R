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
#' @param inla_strat Decides the strategy to approximate marginals by. Simplified Laplace is default
#' @param int_strat Decides the integration strategy for numerical intergrations with INLA. Default is auto, grid for hyperparameters of dimension 2 or less, CCD otherwise
#' @return An object of class `inla` representing the fitted model.
#' @examples
#' data <- data_binomial # Assuming data_binomial is available
#' formula <- y ~ x + f(group, model = "iid")
#' result <- perform_inla_analysis(data, formula, family = "binomial")
#' @export
perform_inla_analysis <- function(data, formula, family, link_func="identity", inla_strat="simplified.laplace", int_strat = "auto", priors = NULL) {

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


  effects <- extract_effects(formula)
  fixed_effects <- effects$fixed_effects

  #   if(!is.null(fixed_effects)){
  if(length(fixed_effects)>1){
    X <- data_copy[, c(fixed_effects)]

    SVD <- SVD_decomp(X)

    data_copy[, c(fixed_effects)] <- SVD$Z
  }

  scaled <- scale(data_copy[, response])


  # Fit the model using INLA
  inla_result <- inla(formula,
                      family = family,
                      data = data_copy,
                      control.family = list(hyper = priors, link = link_func),
                      control.inla = list(strategy = inla_strat, int.strategy=int_strat),
                      control.compute = list(dic = FALSE, return.marginals=TRUE, config=TRUE, waic = TRUE))


  # Return the INLA result object
  return(inla_result)
}
