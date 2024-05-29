#' Extract Fixed and Random Effects from a Model Formula
#'
#' This function identifies and separates fixed and random effects specified in a model formula.
#'
#' @param formula A model formula.
#' @return A list with two elements: `fixed_effects` and `random_effects`,
#'         each containing the names of the respective effects extracted from the formula.
#' @examples
#' formula <- y ~ x + f(group, model = "iid")
#' effects <- extract_effects(formula)
#' @export
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
