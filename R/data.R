#' d.morph data
#'
#' A dataset containing morphological measurements of a house sparrow population on the coast of Helgeland, Norway.
#'
#' @format A data frame with 4625 rows and 24 variables:
#' \describe{
#'   \item{ID}{int: ID of the individual}
#'   \item{indnr}{int: Individual number}
#'   \item{sex}{num: Sex of the individual}
#'   \item{island_current}{Factor: Current island}
#'   \item{adultSNPisland}{int: Adult SNP island}
#'   \item{hatchisland}{int: Hatch island}
#'   \item{hatchyear}{Factor: Hatch year}
#'   \item{year}{Factor: Year}
#'   \item{month}{num: Month}
#'   \item{age}{num: Age}
#'   \item{tarsus}{num: Tarsus length}
#'   \item{wing}{num: Wing length}
#'   \item{mass}{num: Mass}
#'   \item{island}{int: Island}
#'   \item{FGRM}{num: FGRM value}
#'   \item{id}{int: ID}
#'   \item{billD}{num: Bill depth}
#'   \item{inner}{num: Inner measurement}
#'   \item{outer}{num: Outer measurement}
#'   \item{other}{num: Other measurement}
#'   \item{IDC}{int: IDC value}
#'   \item{IDC2}{int: IDC2 value}
#'   \item{IDC3}{int: IDC3 value}
#'   \item{IDC4}{int: IDC4 value}
#' }
#' @source Data from sparrows on the Helgeland coast, Norway
"d.morph_no_ringnr"

#' Cmatrix: Inverse of the Relatedness Matrix A
#'
#' The `Cmatrix` object is a sparse matrix of class `dgCMatrix` representing the inverse of the relatedness matrix \( A \).
#' This matrix is used in computations involving genetic relatedness.
#'
#' @format A sparse matrix of class `dgCMatrix` with 3116 rows and 3116 columns:
#' \describe{
#'   \item{i}{An integer vector indicating the row indices of non-zero elements.}
#'   \item{p}{An integer vector of pointers to the column indices of non-zero elements.}
#'   \item{Dim}{An integer vector of length 2 indicating the dimensions of the matrix (3116 x 3116).}
#'   \item{Dimnames}{A list of length 2 containing the row names (individual IDs) and NULL for column names.}
#'   \item{x}{A numeric vector of non-zero elements in the matrix.}
#'   \item{factors}{A list containing factors of the matrix (usually empty).}
#' }
#' @source Generated as the inverse of the relatedness matrix \( A \).
"Cmatrix"
