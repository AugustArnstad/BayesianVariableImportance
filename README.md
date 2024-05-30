# BayesianVariableImportance
**Bayesian Variable Importance for GLMMs using INLA** is a package developed for my master's thesis at NTNU. It is a further development and extension of the previous package BayesianImportance, available in full at https://github.com/AugustArnstad/BayesianImportance.

`BayesianVariableImportance` is an R package designed to compute Bayesian variable importance metrics for Generalized Linear Mixed Models (GLMMs) utilizing the Integrated Nested Laplace Approximation (INLA) methodology.

## Features
- **Bayesian Variable Importance Computation**: Allows for the quantification of the importance/statistical evidence of predictors in GLMMs in a Bayesian framework. Currently it handles gaussian (identity link), binomial (probit or logit link) and Poisson (log link) data, but more extension are desirable in the near future.
- **INLA Integration**: Leverages the computational advantages of INLA, a popular method for Bayesian inference for latent Gaussian models.
- **Priors**: INLAs penalizing complexity priors are used as defaults, but others can be specified if desirable.
- **Extensible**: Designed with the modern R user in mind, offering a range of utilities to further expand upon the base functionality.
- **Usage examples**: In the package, we showcase a usage example, simulation study, analysis on real data gathered from house sparrows on Helgelandskysten, Norway, a case study comparing it to the `rptR` package and a comparison with Bayesian shrinkage priors methods.

## Installation
To install the latest version of `BayesianVariableImportance` from GitHub you need INLA. You can therefore use the following command:
```R
install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)

# If not already installed, install the 'devtools' package
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# Install BayesianVariableImportance from GitHub
devtools::install_github("AugustArnstad/BayesianVariableImportance")
``` 
To load all other packages needed to use the BVI method to its full extent, use the following command:

```R
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
    Matrix, ggplot2, reshape2, MASS, dplyr, tidyr,
    gridExtra, latex2exp, scales, patchwork, ggtext, mnormt,
    formatR, knitr, devtools, remotes, microbenchmark, 
    rstan, MCMCpack, relaimpo, rptR
)
``` 

## Usage
```R
# Create a model using the INLA formula syntax, deciding on a prior if necessary. 
# Specify correlation structure of random effects in the formula as is standard for INLA models.

glmm_pois <- y_pois ~ X1 + X2 + X3 + f(Z1, model="iid", hyper=list(prec = list(
        prior = "pc.prec",
        param = c(1, 0.01),
        initial = log(1)
      )
   )
)

model_pois <- BayesianVariableImportance::perform_inla_analysis(
    data = datasets$poisson, 
    formula = glmm_pois, 
    family = "poisson", 
    link_func = "log"
)

# Extract the posterior mean importances for each variable. Note that sampling is adviced!

imp_pois <- BayesianVariableImportance::extract_importances(
    model = model_pois, 
    data = datasets$poisson
)


# Compute the variable importance metrics by sampling from the joint posterior. 
# Specify the additive parameter, which represents the additive variance component for computing ICC/Readbility/heritability
# and what scale to compute heritability/repeatability on (see Kruuk - Estimating genetic parameters in natural populations using the ‘animal model’ (2004) and
# Stoffel - rptR: repeatability estimation and variance decomposition by generalized linear mixed-effects models for explanations).
# Note that the additive_param is optional, and not necessary for the computation of the variable importance metrics.

samples_pois <- BayesianVariableImportance::sample_posterior_count(
    model = model_pois, 
    formula = glmm_pois, 
    data = datasets$poisson, 
    n_samp=5000, 
    additive_param = "Z1", 
    repeatability = FALSE
)

# Plot the variable importance metrics
plots_pois <- BayesianVariableImportance::plot_samples(samples_pois)
``` 

## Vignettes

### Usage example
A worked example on relative variable importance calculations for simulated data can be found in the BVI Usage.Rmd file.

### Real data analysis
A study on the heritability of phenotypic traits in house sparrows is conducted in the file Animal_model.Rmd. Please read the method, results and discussion part of my masters thesis for a full description of the study.

### Simulation study
A full simulation study, which has been described, reported and discussed in my master thesis, can be found in the file Simulation study.Rmd
See also `BayesianImportance` for a simulation study on Gaussian LMMs.

### Case studies
A case study, following the vignette of the `rptR` package, which has been described, reported and discussed in my master thesis, can be found in the file Stoffel_comparison.Rmd

### Comparing the BVI method to the Dirichlet and Generalized Decomposition priors on $R^2$
A comparison of the BVI method to the Dirichlet and Generalized Decomposition priors on $R^2$ can be found in the R2D2_GDR2.R file. The R2D2 and GDR2 methods decompose the $R^2$ value, and in that sense can be seen as variable importance methods. However, they were not originally intended to be. Nonetheless, we have used the methods in such a way that the output is analogous to the BVI method. Please read the method, results and discussion part of my masters thesis for more details and references to the original papers.

## Documentation
Further documentation and function references can be found within the package. Use the standard R help and documentation commands to access detailed information about each function.

## Contributing
Contributions to `BayesianVariableImportance` are welcome. Please ensure that you adhere to standard coding practices and include tests for any new features. Open a pull request with details about your changes.

## License
This project is licensed under the MIT License - see the LICENSE.txt file for details

## Acknowledgements
INLA team for their outstanding work on the INLA methodology and R package.
My counsellor Stefanie Muff, Associate Professor, Department of Mathematical Sciences, NTNU Trondheim.

## Contact
For questions, suggestions, or any other inquiries, please contact me at augustarnstad@gmail.com
