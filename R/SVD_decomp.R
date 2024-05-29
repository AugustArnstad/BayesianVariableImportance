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

