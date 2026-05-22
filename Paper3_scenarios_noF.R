n_sim_final <- 100
n_cells_per_group <- 9

c_x_grid <- c(0.5, 1, 1.5)
c_y_grid <- c(0.5, 1, 1.5)

J_main <- 30

p_z_30 <- c(
  rep(0.1, 10),
  rep(0.2, 10),
  rep(0.3, 10)
)

a1_strength_grid <- list(
  very_weak = rep(0.001, J_main),
  weak = rep(0.005, J_main),
  moderately_weak = rep(0.010, J_main),
  moderate = rep(0.050, J_main),
  moderately_strong = rep(0.070, J_main)
)


a1_strength_labels <- c(
  very_weak = "Very Weak IVs",
  weak = "Weak IVs",
  moderately_weak = "Moderately Weak IVs",
  moderate = "Moderate IVs",
  moderately_strong = "Moderately Strong IVs"
)

a1_cd <- a1_strength_grid$weak
a1_cd_label <- "Weak IVs"

n_ab <- 1000
n_grid_cd <- c(500, 1500, 2500, 5000, 7500, 10000)