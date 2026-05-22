source("Paper3_methods_noF.R")
source("Paper3_simulation_noF.R")
source("Paper3_scenarios_noF.R")

args <- commandArgs(trailingOnly = TRUE)

scen <- args[1]
n_i <- as.numeric(args[2])
cell_id <- as.integer(args[3])

if (!exists("n_cells_per_group")) {
  n_cells_per_group <- 9
}


if (scen %in% c("A", "B")) {
  
  strength_ids <- names(a1_strength_grid)
  
  strength_index <- ceiling(cell_id / n_cells_per_group)
  cell_within_group <- ((cell_id - 1) %% n_cells_per_group) + 1
  
  if (strength_index < 1 || strength_index > length(strength_ids)) {
    stop("cell_id is out of range for scenario ", scen)
  }
  
  cxcy_grid <- expand.grid(
    c_x = c_x_grid,
    c_y = c_y_grid,
    stringsAsFactors = FALSE
  )
  
  pars <- cxcy_grid[cell_within_group, ]
  
  strength_id_use <- strength_ids[strength_index]
  a1_use <- a1_strength_grid[[strength_id_use]]
  strength_label_use <- a1_strength_labels[[strength_id_use]]
  a1_value_use <- unique(a1_use)
  
  b1_use <- ifelse(scen == "A", 0, 1)
  
  group_id <- strength_index
  group_label <- strength_label_use
  
} else if (scen %in% c("C", "D")) {
  
 
  cxcy_grid <- expand.grid(
    c_x = c_x_grid,
    c_y = c_y_grid,
    stringsAsFactors = FALSE
  )
  
  cell_within_group <- cell_id
  
  if (cell_within_group < 1 || cell_within_group > nrow(cxcy_grid)) {
    stop("cell_id is out of range for scenario ", scen)
  }
  
  pars <- cxcy_grid[cell_within_group, ]
  
  a1_use <- a1_cd
  strength_id_use <- "weak"
  strength_label_use <- a1_cd_label
  a1_value_use <- unique(a1_use)
  
  b1_use <- ifelse(scen == "C", 0, 1)
  
  group_id <- 1
  group_label <- paste0("n=", n_i)
  
} else {
  stop("scen must be one of A, B, C, or D.")
}


sim_ids <- split(
  seq_len(n_sim_final),
  rep(seq_len(n_cells_per_group), length.out = n_sim_final)
)

sim_ids_this_cell <- sim_ids[[cell_within_group]]
n_sim_this_cell <- length(sim_ids_this_cell)

if (n_sim_this_cell == 0) {
  stop("No simulations assigned to this cell.")
}


seed_base <- switch(
  scen,
  A = 1000,
  B = 2000,
  C = 3000,
  D = 4000
)

seed_use <- seed_base +
  1000000 * group_id +
  10000 * cell_within_group +
  as.integer(n_i)


tmp <- run_sim(
  n_sim = n_sim_this_cell,
  n = n_i,
  a1 = a1_use,
  b1 = b1_use,
  c_x = pars$c_x,
  c_y = pars$c_y,
  p_z = p_z_30,
  sigma_u = 0.5,
  scen = scen,
  strength_id = strength_id_use,
  strength_label = strength_label_use,
  a1_value = a1_value_use,
  seed = seed_use
)

tmp$cell_id <- cell_id
tmp$group_id <- group_id
tmp$group_label <- group_label
tmp$cell_within_group <- cell_within_group
tmp$n_sim_this_cell <- n_sim_this_cell
tmp$n_sim_group_target <- n_sim_final
tmp$n_cells_per_group <- n_cells_per_group
tmp$n_iv <- length(a1_use)
tmp$p_z_pattern <- "pZ_0.1x10_0.2x10_0.3x10"
tmp$method_version <- "noF"

out_file <- sprintf(
  "results_noF_scenario_%s_n%s_cell%02d.rds",
  scen,
  n_i,
  cell_id
)

saveRDS(tmp, out_file)

message("Saved: ", out_file)
message("Scenario: ", scen)
message("Group ID: ", group_id)
message("Group label: ", group_label)
message("Cell within group: ", cell_within_group)
message("n_sim_this_cell: ", n_sim_this_cell)
message("n_sim_group_target: ", n_sim_final)