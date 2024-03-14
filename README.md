# VariableImportanceINLA
**Bayesian Variable Importance for GLMMs using INLA**
This is a package developed for my masters thesis at NTNU. It is a further developement of a previous package BayesianImportance, available in full at https://github.com/AugustArnstad/BayesianImportance.

`VariableImportanceINLA` is an R package designed to compute Bayesian variable importance metrics for Generalized Linear Mixed Models (GLMMs) utilizing the Integrated Nested Laplace Approximation (INLA) methodology.

## Features
- **Bayesian Variable Importance Computation**: Allows for the quantification of the importance of predictors in GLMMs in a Bayesian framework. Currently it only works with LMM's but this is thought to be extended in the near future.
- **INLA Integration**: Leverages the computational advantages of INLA, a popular method for Bayesian inference for latent Gaussian models.
- **Support for Various GLMMs**: Compatible with a wide range of generalized linear mixed models.
- **Extensible**: Designed with the modern R user in mind, offering a range of utilities to further expand upon the base functionality.
- **Priors**: As of right now it uses the default priors that INLA provides. Adding user specified priors is perhaps the most desirable extension.

## Installation
To install the latest version of `VariableImportanceINLA` from GitHub you need INLA. You can therefore use the following command:
```R
install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)

# If not already installed, install the 'devtools' package
if(!require(devtools)) install.packages("devtools")
# Install BayesianImportance
devtools::install_github("AugustArnstad/VariableImportanceINLA")
``` 

## Usage


## Simulation study


## Documentation
Further documentation and function references can be found within the package. Use the standard R help and documentation commands to access detailed information about each function.

## Contributing
Contributions to `VariableImportanceINLA` are welcome. Please ensure that you adhere to standard coding practices and include tests for any new features. Open a pull request with details about your changes.

## License
This project is licensed under the MIT License - see the LICENSE.txt file for details

## Acknowledgements
INLA team for their outstanding work on the INLA methodology and R package.
My counsellor Stefanie Muff, Associate Professor, Department of Mathematical Sciences, NTNU Trondheim
