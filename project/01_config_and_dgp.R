
library(tidyverse)

# ============================================================
# CENTRAL PARAMETER REGISTRY
# ============================================================

DGP_CONFIG <- list(
  p_sex = c(male = 0.65, female = 0.35),
  p_pclass = c("1" = 0.25, "2" = 0.30, "3" = 0.45),

  rules = list(
    "female 1" = list(type = "bimodal_norm",  m1=50,  s1=8,    m2=200, s2=20,  w1=0.5),
    "male 1"   = list(type = "bimodal_lnorm", m1=log(50), s1=0.2, m2=log(200), s2=0.2, w1=0.5),
    "female 2" = list(type = "power_law",     x_min=5, x_max=500),
    "male 2"   = list(type = "uniform",       min=5,   max=150),
    "female 3" = list(type = "gamma",         shape=2, rate=2/20),
    "male 3"   = list(type = "spike_gamma",   m_n=8,   s_n=1.5,  sh_g=1.2, ra_g=1.2/40, w1=0.7)
  )
)

EXPERIMENT_GRID <- crossing(
  mtry = c(1,2,3),
  finite_bounds = c("local", "no"),
  family = c("truncnorm", "unif"),
  seed = 1
)

power_dens <- function(x, x_min = 5, x_max = 500) {
  C <- 1 / (1/x_min - 1/x_max)
  ifelse(x >= x_min & x <= x_max, C / (x^2), 0)
}

get_cond_density <- function(f, sex, pclass) {
  key <- paste(as.character(sex), as.character(pclass))
  p <- DGP_CONFIG$rules[[key]]

  if (is.null(p)) return(rep(0, length(f)))

  switch(
    p$type,
    "bimodal_norm"  = p$w1 * dnorm(f, p$m1, p$s1) + (1-p$w1) * dnorm(f, p$m2, p$s2),
    "bimodal_lnorm" = p$w1 * dlnorm(f, p$m1, p$s1) + (1-p$w1) * dlnorm(f, p$m2, p$s2),
    "uniform"       = dunif(f, p$min, p$max),
    "gamma"         = dgamma(f, shape = p$shape, rate = p$rate),
    "power_law"     = power_dens(f, p$x_min, p$x_max),
    "spike_gamma"   = p$w1 * dnorm(f, p$m_n, p$s_n) +
                       (1-p$w1) * dgamma(f, shape = p$sh_g, rate = p$ra_g)
  )
}

simulate_titanic_dgp <- function(n = 2000, seed = 42) {

  set.seed(seed)

  sex_vec <- sample(
    names(DGP_CONFIG$p_sex),
    n,
    replace = TRUE,
    prob = DGP_CONFIG$p_sex
  )

  pc_vec <- sample(
    names(DGP_CONFIG$p_pclass),
    n,
    replace = TRUE,
    prob = DGP_CONFIG$p_pclass
  )

  fare <- numeric(n)

  for(i in seq_len(n)) {

    key <- paste(sex_vec[i], pc_vec[i])
    p <- DGP_CONFIG$rules[[key]]

    fare[i] <- switch(
      p$type,
      "bimodal_norm"  = if(runif(1) < p$w1) rnorm(1, p$m1, p$s1)
                         else rnorm(1, p$m2, p$s2),

      "bimodal_lnorm" = if(runif(1) < p$w1) rlnorm(1, p$m1, p$s1)
                         else rlnorm(1, p$m2, p$s2),

      "uniform"       = runif(1, p$min, p$max),

      "gamma"         = rgamma(1, shape = p$shape, rate = p$rate),

      "power_law"     = {
        u <- runif(1)
        p$x_min / (1 - u * (1 - p$x_min / p$x_max))
      },

      "spike_gamma"   = if(runif(1) < p$w1)
                           rnorm(1, p$m_n, p$s_n)
                         else
                           rgamma(1, shape = p$sh_g, rate = p$ra_g)
    )
  }

  tibble(
    sex = sex_vec,
    Pclass = factor(pc_vec, levels = c("1", "2", "3")),
    Fare = pmax(fare, 0.5),
    male = as.integer(sex_vec == "male")
  )
}
