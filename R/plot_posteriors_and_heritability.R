#' Plot Posterior Distributions of Variance Components and Heritability
#'
#' Generates plots for the posterior distributions of variance components and, optionally, heritability.
#'
#' @param model An `inla` model object.
#' @param random_effect_name Optional name of the random effect for which heritability is to be plotted.
#' @return A list of plot objects: `variance_plot` for variance components and `heritability_plot` for heritability.
#' @examples
#' # Assuming 'result' is an inla model object
#' plots <- plot_posteriors_and_heritability(result, random_effect_name = "group")
#' @export
plot_posteriors_and_heritability <- function(model, random_effect_name = NULL) {
  # Calculate marginals for variance components
  variance_marginals_list <- lapply(model$marginals.hyperpar, function(x) inla.tmarginal(function(t) 1/t, x))

  # Create data frames for plotting
  df_list <- lapply(names(variance_marginals_list), function(effect_name) {
    data.frame(
      x = variance_marginals_list[[effect_name]][, 1],
      y = variance_marginals_list[[effect_name]][, 2],
      Effect = effect_name
    )
  })
  df_combined <- do.call(rbind, df_list)

  df_combined$Effect <- gsub("Precision for ", "", df_combined$Effect, fixed = TRUE)

  # Plot posterior distributions
  p1 <- ggplot(df_combined, aes(x = x, y = y, color = Effect)) +
    geom_line() +
    labs(title = "Posterior Distributions of Variance Components",
         x = "Variance", y = "Density") +
    theme_minimal() +
    theme(legend.position = "right")

  # Plot heritability if random effect name is provided and correctly calculate heritability
  if (!is.null(random_effect_name) && random_effect_name %in% names(variance_marginals_list)) {
    # Calculate heritability distribution
    names <- names(model$marginals.hyperpar)
    total <- rep(0, length(variance_marginals_list[[random_effect_name]]))
    for (name in names){
      total <- total + variance_marginals_list[[name]][, 1]
    }

    heritability <- variance_marginals_list[[random_effect_name]][, 1]/total

    heritability_density <- variance_marginals_list[[random_effect_name]][ ,2]

    df_heritability <- data.frame(Heritability = heritability, Density = heritability_density)

    p2 <- ggplot(df_heritability, aes(x = Heritability, y = Density)) +
      geom_line() +
      labs(title = paste("Heritability Distribution for", random_effect_name),
           x = "Heritability", y = "Density") +
      theme_minimal()
  } else {
    p2 <- NULL
  }


  return(list(variance_plot = p1, heritability_plot = p2))
}
