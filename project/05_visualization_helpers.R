
library(tidyverse)
library(patchwork)

plot_density_comparison <- function(grid) {

  grid_long <- grid |>
    pivot_longer(
      c(true_dens, arf_dens),
      names_to = "source",
      values_to = "density"
    )

  ggplot(
    grid_long,
    aes(Fare, density, color = source)
  ) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ sex + Pclass, scales = "free_y") +
    scale_x_log10() +
    theme_minimal()
}

combine_plots <- function(plot_list) {

  wrap_plots(plotlist = plot_list) +
    plot_annotation(tag_levels = "A")
}
