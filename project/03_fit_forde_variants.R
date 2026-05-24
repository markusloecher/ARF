
library(tidyverse)
library(arf)

source("01_config_and_dgp.R")

dir.create("results/psi", recursive = TRUE, showWarnings = FALSE)

arf_registry <- readRDS("results/arf_registry.rds")

psi_registry <- EXPERIMENT_GRID |>
  mutate(
    arf_file = glue::glue(
      "results/arf/arf_mtry{mtry}_seed{seed}.rds"
    ),

    psi_file = glue::glue(
      "results/psi/psi_mtry{mtry}_seed{seed}_{family}_{finite_bounds}.rds"
    )
  )

for(i in seq_len(nrow(psi_registry))) {

  row <- psi_registry[i, ]

  if(!file.exists(row$psi_file)) {

    cat("Computing:");print(as.data.frame(row))#print(as.data.frame(row)[1:3])

    arf_obj <- readRDS(row$arf_file)

    psi <- forde(
      arf_obj$arf,
      arf_obj$X,
      family = row$family,
      finite_bounds = row$finite_bounds
    )

    saveRDS(psi, row$psi_file)
  }
}

saveRDS(psi_registry, "results/psi_registry.rds")
