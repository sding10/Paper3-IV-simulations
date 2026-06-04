source("scripts/Paper3_methods_revised_clean.R")
source("scripts/Paper3_simulation_revised_clean.R")
source("scripts/Paper3_scenarios_revised_clean.R")

args <- commandArgs(trailingOnly = TRUE)

cell_id <- as.integer(args[1])

ivnum_grid <- c(1, 5, 10, 15, 20, 25, 30)
n_iv_cell <- ivnum_grid[cell_id]

if (is.na(n_iv_cell)) {
  stop("cell_id must be 1, 2, 3, 4, 5, 6, or 7.")
}

dir.create("results_ivnum", showWarnings = FALSE)

res_ivnum <- run_sim(
  n_sim = n_sim_final,
  n = 1000,
  b1 = 1,
  iv_strength = "strong",
  c_x = 1.5,
  c_y = 1.5,
  n_iv = n_iv_cell,
  pz = rep(0.3, n_iv_cell),
  sigma_u = sigma_u,
  alpha0 = alpha0,
  beta0 = beta0,
  r_min = r_min,
  r_max = r_max,
  scenario = "NumIV",
  setting_id = cell_id,
  seed_base = 9000 + 100000 * cell_id,
  ci_alpha = ci_alpha
)

res_ivnum <- res_ivnum %>%
  dplyr::mutate(
    pz_pattern = "pz_0.3_all",
    method_version = "revised_with_F"
  )

out_file <- sprintf(
  "results_ivnum/Paper3_NumIV_%02dIV_%dsim.rds",
  n_iv_cell,
  n_sim_final
)

saveRDS(res_ivnum, out_file)

message("Saved: ", out_file)