source("scripts/Paper3_methods_2.R")
source("scripts/Paper3_simulation_revised_clean.R")
source("scripts/Paper3_scenarios_revised_clean.R")

args <- commandArgs(trailingOnly = TRUE)

job_id <- as.integer(args[1])

if (is.na(job_id)) {
  stop("Please provide a job ID, e.g., Rscript Paper3_run_cell_revised_clean.R 1")
}

if (job_id < 1 || job_id > n_jobs_final) {
  stop("job_id must be between 1 and ", n_jobs_final)
}

dir.create("results", showWarnings = FALSE)

setting_id <- ceiling(job_id / n_sim_final)
sim_id <- ((job_id - 1) %% n_sim_final) + 1

setting <- scenario_grid[setting_id, ]
if (nrow(setting) != 1) {
  stop("Could not identify exactly one scenario setting for setting_id = ", setting_id)
}

seed_use <- 100000 + job_id

tmp <- sim_once(
  n = setting$n,
  b1 = setting$b1,
  iv_strength = setting$iv_strength,
  c_x = setting$c_x,
  c_y = setting$c_y,
  n_iv = setting$n_iv,
  pz = pz,
  sigma_u = setting$sigma_u,
  alpha0 = setting$alpha0,
  beta0 = setting$beta0,
  r_min = setting$r_min,
  r_max = setting$r_max,
  seed = seed_use,
  ci_alpha = setting$ci_alpha
)

tmp <- tmp %>%
  dplyr::mutate(
    job_id = job_id,
    setting_id = setting_id,
    sim_id = sim_id,
    
    scenario = setting$scenario,
    n = setting$n,
    n_iv = setting$n_iv,
    b1 = setting$b1,
    
    iv_strength = setting$iv_strength,
    
    c_x = setting$c_x,
    c_y = setting$c_y,
    sigma_u = setting$sigma_u,
    alpha0 = setting$alpha0,
    beta0 = setting$beta0,
    
    ci_alpha = setting$ci_alpha,
    pz_pattern = "pz_0.1x10_0.2x10_0.3x10",
    iv_corr = "random_r_0.10_to_0.20",
    method_version = "constant_rho_ivmvb_F_inside_2SPS"
  )

out_file <- sprintf(
  "results/paper3_setting%03d_sim%04d.rds",
  setting_id,
  sim_id
)

saveRDS(tmp, out_file)

message("Saved: ", out_file)
message("Job ID: ", job_id)
message("Setting ID: ", setting_id)
message("Simulation ID within setting: ", sim_id)
message("Scenario: ", setting$scenario)
message("n: ", setting$n)
message("b1: ", setting$b1)
message("IV strength: ", setting$iv_strength)
message("c_x = c_y: ", setting$c_x)
message("n_iv: ", setting$n_iv)
message("Methods: 2SPS, 2SRI, GMM, IV-MVB")