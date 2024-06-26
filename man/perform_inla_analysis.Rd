% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hello.R
\name{perform_inla_analysis}
\alias{perform_inla_analysis}
\title{Perform INLA Analysis on Specified Data and Formula}
\usage{
perform_inla_analysis(
  data,
  formula,
  family,
  link_func = "identity",
  inla_strat = "simplified.laplace",
  int_strat = "auto",
  priors = NULL
)
}
\arguments{
\item{data}{The dataset to be analyzed.}

\item{formula}{The formula specifying the model to be fitted.}

\item{family}{The model family (e.g., "binomial", "poisson").}

\item{link_func}{A link function to be used in the model (e.g., "identity", "logit").}

\item{inla_strat}{Decides the strategy to approximate marginals by. Simplified Laplace is default}

\item{int_strat}{Decides the integration strategy for numerical intergrations with INLA. Default is auto, grid for hyperparameters of dimension 2 or less, CCD otherwise}

\item{priors}{Optional list of prior specifications for model parameters.}
}
\value{
An object of class `inla` representing the fitted model.
}
\description{
Fits a model using the INLA approach, accommodating both fixed and random effects as specified.
The function handles preprocessing such as SVD decomposition and scaling.
}
\examples{
data <- data_binomial # Assuming data_binomial is available
formula <- y ~ x + f(group, model = "iid")
result <- perform_inla_analysis(data, formula, family = "binomial")
}
