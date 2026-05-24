
library(tidyverse)
library(arf)

source("01_config_and_dgp.R")

dir.create("results/grids", recursive = TRUE, showWarnings = FALSE)
dir.create("results/metrics", recursive = TRUE, showWarnings = FALSE)

gen_ll_grid <- function(data, arf_fit, psi, n_grid = 300, verbose = 0) {

  groups <- data |> distinct(sex, Pclass, male)
  fares  <- exp(seq(log(1), log(400), length.out = n_grid))

  grid <- expand_grid(groups, Fare = fares)

  grid$true_dens <- purrr::pmap_dbl(
    list(grid$Fare, grid$sex, as.character(grid$Pclass)),
    get_cond_density
  )
  if (verbose>1) browser()
  # q_df <- as.data.frame(grid[, "Fare", drop = FALSE])
  # 
  # e_df <- as.data.frame(
  #   grid[, c("male", "Pclass"), drop = FALSE]
  # )
  # 
  # grid$arf_ll <- lik(
  #   psi,
  #   query = q_df,
  #   evidence = e_df[1,],
  #   arf = arf_fit,
  #   log = TRUE
  # )
  #grid$arf_dens <- exp(grid$arf_ll)
  
  
  # ============================================================
  # Compute conditional likelihoods groupwise
  # ============================================================
  
  # We split by conditioning variables because lik()
  # does not reliably preserve rowwise alignment between
  # query and evidence.
  
  group_grids <- grid |>
    group_split(sex, Pclass)
  
  group_results <- purrr::map_dfr(group_grids, function(g) {
    
    # ----------------------------------------------------------
    # Evidence: one row per conditioning group
    # ----------------------------------------------------------
    
    e_df <- g |>
      slice(1) |>
      select(male, Pclass) |>
      as.data.frame()
    
    # ----------------------------------------------------------
    # Query: varying Fare values
    # ----------------------------------------------------------
    
    q_df <- g |>
      select(Fare) |>
      as.data.frame()
    
    # ----------------------------------------------------------
    # Compute conditional log likelihoods
    # ----------------------------------------------------------
    
    g$arf_ll <- lik(
      psi,
      query = q_df,
      evidence = e_df,
      arf = arf_fit,
      log = TRUE
    )
    
    g$arf_dens <- exp(g$arf_ll)
    
    g
  })
  
  # Restore original ordering
  grid <- group_results |>
    arrange(sex, Pclass, Fare)

  grid
}

#chatGPT:
# A more “tidyverse-native” version would be:
# grid_nested <- grid |>
#   nest(data = -c(sex, Pclass))
# grid_nested <- grid_nested |>
#   mutate(
#     data = map(data, ...)
#   )
# grid <- unnest(grid_nested, data)


compute_metrics <- function(grid) {

  grid |>
    group_by(sex, Pclass) |>
    arrange(Fare) |>
    mutate(
      delta = Fare - lag(Fare),
      error = abs(true_dens - arf_dens)
    ) |>
    summarize(
      IAE = sum(error * delta, na.rm = TRUE),
      .groups = "drop"
    )
}

train_data <- simulate_titanic_dgp(n = 3000)

psi_registry <- readRDS("results/psi_registry.rds")

for(i in seq_len(nrow(psi_registry))) {

  row <- psi_registry[i, ]

  grid_file <- glue::glue(
    "results/grids/grid_mtry{row$mtry}_seed{row$seed}_{row$family}_{row$finite_bounds}.rds"
  )

  metric_file <- glue::glue(
    "results/metrics/metric_mtry{row$mtry}_seed{row$seed}_{row$family}_{row$finite_bounds}.rds"
  )

  if(!file.exists(grid_file)) {

    cat("Grid:", grid_file, "\\n")

    arf_obj <- readRDS(row$arf_file)
    psi <- readRDS(row$psi_file)

    grid <- gen_ll_grid(
      train_data,
      arf_fit = arf_obj$arf,
      psi = psi,
      verbose = 2
    )

    saveRDS(grid, grid_file)

    metrics <- compute_metrics(grid)

    saveRDS(metrics, metric_file)
  }
}
