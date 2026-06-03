library(dplyr)
library(tidyr)
library(purrr)

## ------------------------------------------------------------
## Simulation size
## ------------------------------------------------------------

n_sim_final <- 100   # pilot run; change to 8000 for final run

## ------------------------------------------------------------
## Fixed design parameters
## ------------------------------------------------------------

methods_main <- c("2SPS", "2SRI", "GMM", "IV-MVB")

n_iv_main <- 30

alpha0 <- 0
beta0 <- 0
sigma_u <- 0.5
ci_alpha <- 0.05

pz <- c(
  rep(0.1, 10),
  rep(0.2, 10),
  rep(0.3, 10)
)

r_min <- 0.10
r_max <- 0.20

## ------------------------------------------------------------
## Scenario design values
## ------------------------------------------------------------

b1_null <- 0
b1_alt <- c(1, 2, 3)

conf_grid <- c(0.5, 1.0, 1.5)

n_fixed <- 1000
n_grid <- c(500, 1500, 2500, 5000, 7500, 10000)

iv_strength_levels <- c(
  "weak",
  "moderate",
  "strong",
  "very_strong"
)

## ------------------------------------------------------------
## IV strength definitions
## ------------------------------------------------------------

make_alpha_z <- function(iv_strength, n_iv = 30) {
  switch(
    iv_strength,
    
    weak =
      seq(0.001, 0.005, length.out = n_iv),
    
    moderate =
      seq(0.4, 0.8, length.out = n_iv),
    
    strong =
      rep(1.2, times = n_iv),
    
    very_strong =
      rep(3, times = n_iv),
    
    stop("Unknown IV strength: ", iv_strength)
  )
}

## ------------------------------------------------------------
## Scenario grid from revised Table 1
## ------------------------------------------------------------

scenario_grid <- bind_rows(
  
  ## Scenario A
  crossing(
    scenario = "A",
    b1 = b1_null,
    n = n_fixed,
    iv_strength = iv_strength_levels,
    conf = conf_grid
  ),
  
  ## Scenario B
  crossing(
    scenario = "B",
    b1 = b1_alt,
    n = n_fixed,
    iv_strength = iv_strength_levels,
    conf = conf_grid
  ),
  
  ## Scenario C
  crossing(
    scenario = "C",
    b1 = b1_null,
    n = n_grid,
    iv_strength = "weak",
    conf = conf_grid
  ),
  
  ## Scenario D
  crossing(
    scenario = "D",
    b1 = b1_alt,
    n = n_grid,
    iv_strength = "weak",
    conf = conf_grid
  )
) |>
  mutate(
    setting_id = row_number(),
    
    c_x = conf,
    c_y = conf,
    
    n_iv = n_iv_main,
    
    alpha0 = alpha0,
    beta0 = beta0,
    sigma_u = sigma_u,
    ci_alpha = ci_alpha,
    
    r_min = r_min,
    r_max = r_max,
    
    alpha_min = map_dbl(
      iv_strength,
      ~ min(make_alpha_z(.x, n_iv = n_iv_main))
    ),
    
    alpha_max = map_dbl(
      iv_strength,
      ~ max(make_alpha_z(.x, n_iv = n_iv_main))
    )
  ) |>
  select(
    setting_id,
    scenario,
    n,
    n_iv,
    b1,
    iv_strength,
    alpha_min,
    alpha_max,
    c_x,
    c_y,
    alpha0,
    beta0,
    sigma_u,
    ci_alpha,
    r_min,
    r_max
  )

## ------------------------------------------------------------
## Checks
## ------------------------------------------------------------

stopifnot(nrow(scenario_grid) == 120)
stopifnot(length(pz) == n_iv_main)

scenario_counts <- scenario_grid |>
  count(scenario, name = "n_settings")

print(scenario_counts)

n_jobs_final <- nrow(scenario_grid) * n_sim_final

cat("Total scenario settings:", nrow(scenario_grid), "\n")
cat("Simulation replicates per setting:", n_sim_final, "\n")
cat("Total jobs/replicates:", n_jobs_final, "\n")