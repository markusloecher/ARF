library(arf)
library(tidyverse)
library(patchwork)

# ============================================================
# 1. CENTRAL PARAMETER REGISTRY (Single Source of Truth)
# ============================================================
# We define everything here so that simulation, likelihood, 
# and grid generation are always in sync.

DGP_CONFIG <- list(
  # Marginal Probabilities
  p_sex = c(male = 0.65, female = 0.35),
  p_pclass = c("1" = 0.25, "2" = 0.30, "3" = 0.45),
  
  # Conditional Fare Logic: f(Fare | Sex, Pclass)
  rules = list(
    "female 1" = list(type = "bimodal_norm",  m1=50,  s1=8,    m2=200, s2=20,  w1=0.5),
    "male 1" = list(type = "unimodal_norm",  m=50,  s=10),
    #"male 1"   = list(type = "bimodal_lnorm", m1=log(50), s1=0.2, m2=log(250), s2=0.2, w1=0.5),
    "female 2" = list(type = "power_law",     x_min=5, x_max=500),
    "male 2"   = list(type = "uniform",       min=5,   max=150),
    "female 3" = list(type = "gamma",         shape=2, rate=2/20),
    "male 3"   = list(type = "spike_gamma",   m_n=8,   s_n=1.5,  sh_g=1.2, ra_g=1.2/40, w1=0.7)
  )
)

# Helper: Power Law Density (Normalized)
power_dens <- function(x, x_min = 5, x_max = 500) {
  C <- 1 / (1/x_min - 1/x_max)
  ifelse(x >= x_min & x <= x_max, C / (x^2), 0)
}

# Unified Helper for Conditional Density
get_cond_density <- function(f, sex, pclass) {
  key <- paste(as.character(sex), as.character(pclass))
  p <- DGP_CONFIG$rules[[key]]
  if (is.null(p)) return(rep(0, length(f)))
  
  switch(p$type,
         "bimodal_norm"  = p$w1 * dnorm(f, p$m1, p$s1) + (1-p$w1) * dnorm(f, p$m2, p$s2),
         "bimodal_lnorm" = p$w1 * dlnorm(f, p$m1, p$s1) + (1-p$w1) * dlnorm(f, p$m2, p$s2),
         "unimodal_norm" = dnorm(f, p$m, p$s),
         "power_law"     = power_dens(f, p$x_min, p$x_max),
         "uniform"       = dunif(f, p$min, p$max),
         "gamma"         = dgamma(f, shape = p$shape, rate = p$rate),
         "spike_gamma"   = p$w1 * dnorm(f, p$m_n, p$s_n) + (1-p$w1) * dgamma(f, shape = p$sh_g, rate = p$ra_g)
  )
}

# ============================================================
# 2. SIMULATION & LIKELIHOOD
# ============================================================

simulate_titanic_dgp <- function(n = 2000, seed = 42) {
  set.seed(seed)
  
  # Sample sex and pclass as characters to prevent indexing bugs
  sex_vec <- sample(names(DGP_CONFIG$p_sex), n, replace = TRUE, prob = DGP_CONFIG$p_sex)
  pc_vec  <- sample(names(DGP_CONFIG$p_pclass), n, replace = TRUE, prob = DGP_CONFIG$p_pclass)
  
  fare <- numeric(n)
  for (i in seq_len(n)) {
    key <- paste(sex_vec[i], pc_vec[i])
    p   <- DGP_CONFIG$rules[[key]]
    
    fare[i] <- switch(p$type,
                      "bimodal_norm"  = if(runif(1) < p$w1) rnorm(1, p$m1, p$s1) else rnorm(1, p$m2, p$s2),
                      "bimodal_lnorm" = if(runif(1) < p$w1) rlnorm(1, p$m1, p$s1) else rlnorm(1, p$m2, p$s2),
                      "unimodal_norm" = rnorm(1, p$m, p$s),
                      "power_law"     = { u <- runif(1); p$x_min / (1 - u * (1 - p$x_min / p$x_max)) },
                      "uniform"       = runif(1, p$min, p$max),
                      "gamma"         = rgamma(1, shape = p$shape, rate = p$rate),
                      "spike_gamma"   = if(runif(1) < p$w1) rnorm(1, p$m_n, p$s_n) else rgamma(1, shape = p$sh_g, rate = p$ra_g)
    )
  }
  
  tibble(
    sex = sex_vec, 
    Pclass = factor(pc_vec, levels = c("1", "2", "3")), 
    Fare = pmax(fare, 0.5),
    male = as.integer(sex_vec == "male")
  )
}

true_log_lik <- function(data) {
  map_dbl(1:nrow(data), function(i) {
    f  <- data$Fare[i]
    s  <- as.character(data$sex[i])
    pc <- as.character(data$Pclass[i])
    
    # Joint = P(Sex) * P(Pclass) * f(Fare | Sex, Pclass)
    p_mrg  <- DGP_CONFIG$p_sex[s] * DGP_CONFIG$p_pclass[pc]
    f_cond <- get_cond_density(f, s, pc)
    
    log(p_mrg) + log(f_cond)
  })
}

# ============================================================
# 3. FITTING & GRID GENERATION
# ============================================================

fit_arf <- function(data, num.trees = 100, min.node.size = 10, seed = 42) {
  set.seed(seed)
  X   <- data[, c("male", "Pclass", "Fare")]
  arf  <- adversarial_rf(X, num_trees = num.trees, min_node_size = min.node.size)
  psi  <- forde(arf, X)
  list(arf = arf, psi = psi, X = X)
}

gen_ll_grid <- function(data, arf_fit, n_grid = 300, use_arf_conditional = TRUE) {
  # 1. Create the grid combinations
  groups <- data |> distinct(sex, Pclass, male)
  fares  <- exp(seq(log(1), log(500), length.out = n_grid))
  grid   <- expand_grid(groups, Fare = fares)
  
  # 2. Calculate True Conditional Density (Red Line)
  grid$true_dens <- pmap_dbl(list(grid$Fare, grid$sex, as.character(grid$Pclass)), 
                             get_cond_density)
  
  # 3. Calculate ARF Density (Blue Line)
  if (use_arf_conditional) {
    # This route uses the 'evidence' argument but ensures 1-to-1 row alignment
    # by stripping tibble metadata which often triggers the Cartesian join error.
    q_df <- as.data.frame(grid[, "Fare", drop = FALSE])
    e_df <- as.data.frame(grid[, c("male", "Pclass"), drop = FALSE])
    
    grid$arf_ll <- lik(
      arf_fit$psi, 
      query = q_df, 
      evidence = e_df,
      arf = arf_fit$arf, 
      log = TRUE
    )
    grid$arf_dens <- exp(grid$arf_ll)
    
  } else {
    # ROUTE B: Manual Joint / Marginal
    # To avoid the 'obs' column error, we compute likelihoods separately
    # and combine them using standard vector math rather than data.table joins.
    
    # A. Joint Likelihood f(Fare, Sex, Pclass)
    joint_input <- as.data.frame(grid[, c("male", "Pclass", "Fare")])
    joint_ll <- lik(arf_fit$psi, joint_input, arf = arf_fit$arf, log = TRUE)
    
    # B. Marginal Likelihood f(Sex, Pclass)
    # We evaluate the marginal at the same categorical points
    marg_input <- as.data.frame(grid[, c("male", "Pclass")])
    marg_ll <- lik(arf_fit$psi, marg_input, arf = arf_fit$arf, log = TRUE)
    
    print("ARF log marginals", marg_ll)
    # C. Conditional = exp(Joint_LL - Marg_LL)
    grid$arf_dens <- exp(joint_ll - marg_ll)
  }
  
  return(grid)
}

arf_marg_ll <- function(data, arf_fit, n_grid = 300) {
  # 1. Create the grid combinations
  groups <- data |> distinct(sex, Pclass, male)
  fares  <- exp(seq(log(1), log(500), length.out = n_grid))
  grid   <- expand_grid(groups, Fare = fares)
  
  #marg_input <- as.data.frame(grid[, c("male", "Pclass")])
  # 1. Get unique group combinations
  unique_groups <- grid |> distinct(male, sex, Pclass)
  
  # 2. Compute marginals once per group
  unique_groups$arf_marginal <- exp(lik(
    arf_fit$psi, 
    as.data.frame(unique_groups[, c("male", "Pclass")]), 
    arf = arf_fit$arf, 
    log = TRUE
  ))
  
  actual_marginals <- data |>
    count(sex, Pclass) |>
    mutate(p_mrg_actual = n / nrow(data)) |>
    select(-n)
  
  print(actual_marginals)
  
  unique_groups$arf_marginal <- actual_marginals
    
  #browser()
  return(unique_groups)
  print("ARF log marginals", unique_groups$arf_marginal)
  
  
  
  # 3. Join back to the full grid
  grid <- grid |>
    left_join(unique_groups, by = c("male", "sex", "Pclass")) |>
    mutate(arf_dens = exp(arf_ll) / arf_marginal)
  
  
  marg_ll <- lik(arf_fit$psi, marg_input, arf = arf_fit$arf, log = TRUE)
  
  
}


gen_ll_grid_try2 <- function(data, arf_fit, n_grid = 300, use_arf_conditional = TRUE, data_table_bug = TRUE) {
  # 1. Create the grid
  groups <- data |> distinct(sex, Pclass, male)
  fares  <- exp(seq(log(1), log(500), length.out = n_grid))
  grid   <- expand_grid(groups, Fare = fares)
  
  # 2. Calculate True Conditional Density
  grid$true_dens <- pmap_dbl(list(grid$Fare, grid$sex, as.character(grid$Pclass)), 
                             get_cond_density)
  
  # 3. Calculate ARF Density
  if (use_arf_conditional) {
    if (data_table_bug){
      # Alternative if 'evidence' remains buggy:
      joint_lik <- lik(arf_fit$psi, grid |> select(male, Pclass, Fare), arf = arf_fit$arf)
      marginal_lik <- lik(arf_fit$psi, grid |> select(male, Pclass), arf = arf_fit$arf)
      grid$arf_dens <- joint_lik / marginal_lik
    } else {
      # FIX: Ensure query and evidence are treated as aligned rows.
      # We pass them as simple data frames and ensure no row names/extra attributes 
      # interfere with the internal data.table join.
      grid$arf_ll <- lik(
        arf_fit$psi, 
        query = as.data.frame(grid[, "Fare", drop = FALSE]), 
        evidence = as.data.frame(grid[, c("male", "Pclass"), drop = FALSE]),
        arf = arf_fit$arf, 
        log = TRUE
      )
      grid$arf_dens <- exp(grid$arf_ll)
    }
    
  } else {
    # Manual route (unchanged)
    actual_marginals <- data |>
      count(sex, Pclass) |>
      mutate(p_mrg_actual = n / nrow(data))
    
    grid <- grid |>
      left_join(actual_marginals, by = c("sex", "Pclass")) |>
      mutate(
        joint_ll = lik(arf_fit$psi, pick(male, Pclass, Fare), arf = arf_fit$arf, log = TRUE),
        arf_dens = exp(joint_ll) / p_mrg_actual
      )
  }
  
  return(grid)
}

gen_ll_grid_try1 <- function(data, arf_fit, n_grid = 600, use_arf_conditional = TRUE) {
  # 1. Create the grid
  groups <- data |> distinct(sex, Pclass, male)
  fares  <- exp(seq(log(1), log(500), length.out = n_grid))
  grid   <- expand_grid(groups, Fare = fares)
  
  # 2. Calculate True Conditional Density (for the red line)
  grid$true_dens <- pmap_dbl(list(grid$Fare, grid$sex, as.character(grid$Pclass)), 
                             get_cond_density)
  
  # 3. Calculate ARF Density (for the blue line)
  if (use_arf_conditional) {
    # ROUTE A: Direct Conditional Likelihood via 'evidence'
    # The ARF calculates f(Fare | male, Pclass) directly.
    grid$arf_ll <- lik(
      arf_fit$psi, 
      query = grid |> select(Fare), 
      evidence = grid |> select(male, Pclass),
      arf = arf_fit$arf, 
      log = TRUE
    )
    grid$arf_dens <- exp(grid$arf_ll)
    
  } else {
    # ROUTE B: Manual Normalization
    # Calculate Joint Likelihood and divide by empirical marginals.
    actual_marginals <- data |>
      count(sex, Pclass) |>
      mutate(p_mrg_actual = n / nrow(data))
    
    grid <- grid |>
      left_join(actual_marginals, by = c("sex", "Pclass")) |>
      mutate(
        joint_ll = lik(arf_fit$psi, pick(male, Pclass, Fare), arf = arf_fit$arf, log = TRUE),
        arf_dens = exp(joint_ll) / p_mrg_actual
      )
  }
  
  return(grid)
}

gen_ll_grid_arf_margs <- function(data, arf_fit, n_grid = 750) {
  # 1. Create a clean grid of all categories and fares
  # This ensures we have a predictable row order for the likelihood function
  groups <- data |> distinct(sex, Pclass, male)
  fares  <- exp(seq(log(1), log(500), length.out = n_grid))
  grid   <- expand_grid(groups, Fare = fares)
  
  # 2. Calculate actual empirical marginals from the simulation data
  # Using actual proportions (e.g., 0.64 instead of theoretical 0.65) 
  # prevents the ARF density from being scaled incorrectly in the plot.
  actual_marginals <- data |>
    count(sex, Pclass) |>
    mutate(p_mrg_actual = n / nrow(data)) |>
    select(-n)
  print(paste("empirical marginals", actual_marginals))
  # 3. Calculate densities and join with actual marginals
  grid <- grid |>
    left_join(actual_marginals, by = c("sex", "Pclass")) |>
    mutate(
      # True Conditional Density from your DGP_CONFIG
      true_dens = pmap_dbl(list(Fare, sex, as.character(Pclass)), get_cond_density),
      
      # ARF Joint Likelihood (f(Fare, Sex, Pclass))
      # We use pick() to ensure the columns match the ARF's expected input
      arf_ll = lik(arf_fit$psi, pick(male, Pclass, Fare), arf = arf_fit$arf, log = TRUE),
      
      # Convert Joint Density -> Conditional Density: 
      # f(Fare | Sex, Pclass) = f(Fare, Sex, Pclass) / P(Sex, Pclass)
      arf_dens = exp(arf_ll) / p_mrg_actual
    )
  
  return(grid)
}

gen_ll_grid_th_margs <- function(data, arf_fit, n_grid = 750) {
  # 1. Create a clean grid of all categories and fares
  groups <- data |> distinct(sex, Pclass, male)
  fares  <- exp(seq(log(1), log(500), length.out = n_grid))
  grid   <- expand_grid(groups, Fare = fares)
  
  # 2. Calculate densities
  grid <- grid |>
    mutate(
      # True Conditional Density
      true_dens = pmap_dbl(list(Fare, sex, as.character(Pclass)), get_cond_density),
      
      # ARF Joint Likelihood
      arf_ll = lik(arf_fit$psi, pick(male, Pclass, Fare), arf = arf_fit$arf, log = TRUE),
      
      # Marginal probabilities for normalization (Joint -> Conditional)
      p_mrg = DGP_CONFIG$p_sex[sex] * DGP_CONFIG$p_pclass[as.character(Pclass)],
      arf_dens = exp(arf_ll) / p_mrg
    )
  
  return(grid)
}

# ============================================================
# 4. VISUALIZATION
# ============================================================

plot_density_comparison <- function(data, grid) {
  grid_long <- grid |>
    pivot_longer(c(true_dens, arf_dens), names_to = "source", values_to = "density") |>
    mutate(source = recode(source, true_dens = "True DGP", arf_dens = "ARF estimate"))
  
  ggplot(grid_long, aes(x = Fare, y = density, color = source, linetype = source)) +
    geom_line(linewidth = 0.9) +
    # Rug plot of actual data observations to verify peaks
    geom_rug(data = data, aes(x = Fare), sides = "b", inherit.aes = FALSE, alpha = 0.2) +
    facet_wrap(~ sex + Pclass, scales = "free_y", ncol = 3) +
    scale_x_log10() +
    scale_color_manual(values = c("True DGP" = "#E41A1C", "ARF estimate" = "#377EB8")) +
    theme_minimal() +
    labs(
      title = "Conditional Density: f(Fare | Sex, Pclass)",
      subtitle = "ARF joint density normalized by theoretical marginals for fair comparison",
      x = "Fare (log scale)", 
      y = "Density"
    ) +
    theme(legend.position = "bottom")
}

# ============================================================
# 5. Difference between densities
# ============================================================

IAE = function(grid, ndigits = 2){
  
  #1. Integrated Absolute Error (IAE)
  #This is the most intuitive "single number" for your plot. It represents the total area between the red and blue curves.
  
  iae_results <- grid |>
    group_by(sex, Pclass) |>
    arrange(Fare) |>
    mutate(
      # Width of the interval on the x-axis
      delta = Fare - lag(Fare),
      # Absolute vertical difference
      error = abs(true_dens - arf_dens)
    ) |>
    summarize(
      IAE = round(sum(error * delta, na.rm = TRUE), ndigits),
      .groups = "drop"
    ) |>
    arrange(desc(IAE))
  
  
  tvd_results <- grid |>
    group_by(sex, Pclass) |>
    arrange(Fare) |>
    mutate(
      # 1. Calculate the integration width
      delta = Fare - lag(Fare, default = first(Fare)),
      
      # 2. Normalize the ARF density so the area equals 1.0
      # This removes the calibration error of the marginals
      current_area = sum(arf_dens * delta, na.rm = TRUE),
      arf_dens_norm = arf_dens / current_area
    ) |>
    summarize(
      # 3. TVD = 0.5 * Integral |f_true - f_est_norm|
      # This is always between 0 (perfect) and 1 (no overlap)
      TVD = round(0.5 * sum(abs(true_dens - arf_dens_norm) * delta, na.rm = TRUE), ndigits),
      .groups = "drop"
    ) |>
    arrange(desc(TVD))
  
  #print(tvd_results)
  return(list(IAE =  as.data.frame(iae_results), TVD = as.data.frame(tvd_results)))
}
