
library(tidyverse)
library(arf)

source("01_config_and_dgp.R")

dir.create("results/arf", recursive = TRUE, showWarnings = FALSE)

fit_single_arf <- function(data, mtry, seed, num_trees = 100) {

  set.seed(seed)

  X <- data |> select(male, Pclass, Fare)

  arf_fit <- adversarial_rf(
    X,
    mtry = mtry,
    num_trees = num_trees
  )

  list(
    arf = arf_fit,
    X = X
  )
}

train_data <- simulate_titanic_dgp(n = 3000)

arf_registry <- EXPERIMENT_GRID |>
  distinct(mtry, seed) |>
  mutate(
    file = glue::glue("results/arf/arf_mtry{mtry}_seed{seed}.rds")
  )

for(i in seq_len(nrow(arf_registry))) {

  row <- arf_registry[i, ]

  if(!file.exists(row$file)) {

    cat("Fitting:", row$file, "\\n")

    obj <- fit_single_arf(
      train_data,
      mtry = row$mtry,
      seed = row$seed
    )

    saveRDS(obj, row$file)
  }
}

saveRDS(arf_registry, "results/arf_registry.rds")
