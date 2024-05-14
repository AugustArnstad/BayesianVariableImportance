# BayesianVariableImportance
**Bayesian Variable Importance for GLMMs using INLA**
\\
This is a package developed for my masters thesis at NTNU. It is a further developement and extension of the previous package BayesianImportance, available in full at https://github.com/AugustArnstad/BayesianImportance.

`BayesianVariableImportance` is an R package designed to compute Bayesian variable importance metrics for Generalized Linear Mixed Models (GLMMs) utilizing the Integrated Nested Laplace Approximation (INLA) methodology.

## Features
- **Bayesian Variable Importance Computation**: Allows for the quantification of the importance/statistical evidence of predictors in GLMMs in a Bayesian framework. Currently it handles gaussian (identity link), binomial (probit or logit link) and Poisson (log link) data, but more extension are desirable in the near future.
- **INLA Integration**: Leverages the computational advantages of INLA, a popular method for Bayesian inference for latent Gaussian models.
- **Support for Various GLMMs**: Compatible with a wide range of generalized linear mixed models.
- **Extensible**: Designed with the modern R user in mind, offering a range of utilities to further expand upon the base functionality.
- **Priors**: INLAs default penalize complexity priors are used as defaults, but others can be specified if desirable

## Installation
To install the latest version of `BayesianVariableImportance` from GitHub you need INLA. You can therefore use the following command:
```R
install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)

# If not already installed, install the 'devtools' package
if(!require(devtools)) install.packages("devtools")
# Install BayesianVariableImportance from GitHub
devtools::install_github("AugustArnstad/BayesianVariableImportance")
``` 

## Usage
```R
# Create a model using the INLA formula syntax, deciding on a prior if necessary. Specify correlation structure of random effects in the formula
# as is standard for INLA models.
glmm_pois <- y_pois ~ X1 + X2 + X3 + f(Z1, model="iid", hyper=list(prec = list(
        prior = "pc.prec",
        param = c(1, 0.01),
        initial = log(1)
      ))
)
model_pois <- BayesianVariableImportance::perform_inla_analysis(datasets$poisson, glmm_pois, family = "poisson", link_func = "log")

# Extract some summary statistics from the fitted model.
imp_pois <- BayesianVariableImportance::extract_importances(model_pois, 
                                datasets$poisson,
                                random_names=c("Z1"), 
                                fixed_names=c("X1", "X2", "X3"))


# Compute the variable importance metrics by sampling from the joint posterior. Specify the additive parameter, which represents the additive 
# variance component of an animal model (see Kruuk - Estimating genetic parameters in natural populations using the ‘animal model’ (2004))
samples_pois <- BayesianVariableImportance::sample_posterior_count(model_pois, glmm_pois, datasets$poisson, n_samp=5000, additive_param = "Z1")

# Plot the variable importance metrics
plots_pois <- BayesianVariableImportance::plot_samples(samples_pois)
``` 

A worked example can be found in the BVI Usage.Rmd file.
## Simulation study
A full simulation study, which has been described, reported and discussed in my master thesis, can be found in the file Simulation study.Rmd

## Documentation
Further documentation and function references can be found within the package. Use the standard R help and documentation commands to access detailed information about each function.

## Contributing
Contributions to `BayesianVariableImportance` are welcome. Please ensure that you adhere to standard coding practices and include tests for any new features. Open a pull request with details about your changes.

## License
This project is licensed under the MIT License - see the LICENSE.txt file for details

## Acknowledgements
INLA team for their outstanding work on the INLA methodology and R package.
My counsellor Stefanie Muff, Associate Professor, Department of Mathematical Sciences, NTNU Trondheim
