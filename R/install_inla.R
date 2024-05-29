#' Ensure the INLA package is installed
#'
#' This function checks if INLA is installed and installs it if necessary.
#' @export
install_inla <- function() {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    install.packages("INLA", repos = c(INLA = "https://inla.r-inla-download.org/R/stable"), type = "source")
  }
  library(INLA)
}
