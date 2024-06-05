#' Plot Samples from Posterior Distributions
#'
#' Generates plots for the scaled importance of fixed and random effects,
#' as well as marginal and conditional R2 values from sampled posterior distributions.
#'
#' @param samples A list object containing sampled values and derived quantities from `sample_posterior`.
#' @return A list of ggplot objects for each of the plotted quantities.
#' @examples
#' # Assuming 'samples' contains results from `sample_posterior`
#' plots <- plot_samples(samples)
#' @export
plot_samples <- function(samples) {
  plots <- list() # Initialize an empty list to store plots

  # Check and plot scaled importance for fixed effects
  if (!is.null(samples$scaled_importance_samples) && dim(samples$scaled_importance_samples)[2] > 0) {
    melted_scaled_importance <- melt(as.data.frame(samples$scaled_importance_samples))
    fixed_effects_plot <- ggplot(melted_scaled_importance, aes(x = value)) +
      #geom_histogram(aes(y = ..density..), fill = "#C6CDF7", alpha = 0.5) +
      geom_histogram(aes(y = after_stat(density)), fill = "#C6CDF7", alpha = 0.5) +
      #stat_density(geom = 'line', colour = "#E6C6DF", adjust=1.5, linewidth =1.5) +
      geom_density(colour = "#E6C6DF", adjust = 1.5, linewidth=1.5) +
      labs(title = "Fixed Effects", x = "Relative Importance", y = "Density") +
      theme_minimal() +
      # theme(text = element_text(family="LM Roman 10"),
      #       plot.title = element_text(size = 18, hjust = 0.5, face="bold")) +
      facet_wrap(~ variable, scales = "free")
    plots$fixed_effects <- fixed_effects_plot
  }

  # Check and plot scaled importance for random effects
  if (!is.null(samples$scaled_random_samples)) {
    melted_scaled_random <- melt(as.data.frame(samples$scaled_random_samples))
    random_effects_plot <- ggplot(melted_scaled_random, aes(x = value)) +
      #geom_histogram(aes(y = ..density..), fill = "#C6F7CD", alpha = 0.5) +
      geom_histogram(aes(y = after_stat(density)), fill = "#C6F7CD", alpha = 0.5) +
      #stat_density(geom = 'line', colour = "#E6C6DF", adjust=1.5, linewidth =1.5) +
      geom_density(colour = "#E6C6DF", adjust = 1.5, linewidth=1.5) +
      labs(title = "Random Effects", x = "Relative Importance", y = "Density") +
      theme_minimal() +
      # theme(text = element_text(family="LM Roman 10"),
      #       plot.title = element_text(size = 18, hjust = 0.5, face="bold")) +
      facet_wrap(~ variable, scales = "free")
    plots$random_effects <- random_effects_plot
  }

  # Adjusted handling for heritability if it has values and not all are NA
  if (!is.null(samples$heritability) && any(!is.na(samples$heritability))) {
    # Extract the actual column name for heritability
    heritability_colname <- names(samples$heritability)
    heritability_plot <- ggplot(samples$heritability, aes(x = !!sym(heritability_colname))) +
      #geom_histogram(aes(y = ..density..), fill = "purple", alpha = 0.5) +
      geom_histogram(aes(y = after_stat(density)), fill = "purple", alpha = 0.5) +
      #stat_density(geom = 'line', colour = "orange", adjust=1.5, linewidth =1.5) +
      geom_density(color = "orange", adjust = 1.5, linewidth=1.5) +
      labs(title = paste("Heritability of:", heritability_colname), x = heritability_colname, y = "Frequency") +
      theme_minimal() +
      # theme(text = element_text(family="LM Roman 10"),
      #       plot.title = element_text(size = 18, hjust = 0.5, face="bold")) +
      plots$heritability <- heritability_plot
  }

  if (!is.null(samples$R2_marginal) && !is.null(samples$R2_conditional)) {
    # Combine R2_marginal and R2_conditional into a single data frame
    R2_data <- data.frame(R2_marginal = samples$R2_marginal$`Marginal R2`,
                          R2_conditional = samples$R2_conditional$`Conditional R2`)

    # Convert from wide to long format
    R2_long <- pivot_longer(R2_data, cols = c(R2_marginal, R2_conditional), names_to = "Type", values_to = "Value")

    R2_long$Type <- factor(R2_long$Type, levels = c("R2_marginal", "R2_conditional"))

    # Plot
    R2_plot <- ggplot(R2_long, aes(x = Value, fill = Type)) +
      #geom_histogram(aes(y = ..density..), alpha = 0.5, position = "identity") +
      geom_histogram(aes(y = after_stat(density)), alpha = 0.5, position = "identity") +
      #stat_density(geom = 'line', colour = "#E6C6DF", adjust=1.5, linewidth =1.5) +
      geom_density(colour = "#E6C6DF", alpha = 0.75, adjust = 1.5, linewidth=1.5) +
      labs(title = "Marginal and Conditional R2", x = "R2 estimate", y = "Density") +
      scale_fill_manual(values = c("R2_marginal" = "#C6CDF7", "R2_conditional" = "#C6F7CD")) +
      theme_minimal() +
      # theme(text = element_text(family="LM Roman 10"),
      #       plot.title = element_text(size = 18, hjust = 0.5, face="bold")) +
      facet_wrap(~ Type, scales = "free")

    plots$R2 <- R2_plot
  }

  return(plots)
}
