source("Paper3_methods_noF.R")
source("Paper3_simulation_noF.R")
source("Paper3_scenarios_noF.R")

args <- commandArgs(trailingOnly = TRUE)

cell_id <- as.integer(args[1])

num_iv_grid <- c(1, 5, 10, 15, 20, 25, 30)
num_iv <- num_iv_grid[cell_id]

if (is.na(num_iv)) {
  stop("cell_id must be 1, 2, 3, 4, 5, 6, or 7.")
}

tmp <- run_sim(
  n_sim = n_sim_final,
  n = 1000,
  a1 = rep(0.05, num_iv),
  b1 = 1,
  c_x = 1.5,
  c_y = 1.5,
  p_z = rep(0.3, num_iv),
  sigma_u = 0.5,
  scen = "NumIV",
  strength_id = "moderate",
  strength_label = "Moderate IVs",
  a1_value = 0.05,
  seed = 9000 + 100000 * cell_id
)

tmp$num_iv <- num_iv
tmp$n_iv <- num_iv
tmp$p_z_pattern <- "pZ_0.3_all"
tmp$method_version <- "noF"

out_file <- sprintf(
  "results_noF_NumIV_%02dIV_100sim.rds",
  num_iv
)

saveRDS(tmp, out_file)

message("Saved: ", out_file)