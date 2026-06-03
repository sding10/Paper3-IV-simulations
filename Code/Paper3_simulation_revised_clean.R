library(MASS)
library(Matrix)
library(dplyr)
library(purrr)
library(tidyr)
library(readr)

source("scripts/Paper3_methods_revised_clean.R")


sim_binary_iv <- function(n, pz, r_min = 0.10, r_max = 0.20) {
  n_iv <- length(pz)

  R <- diag(1, n_iv)

  if (n_iv > 1) {
    idx <- upper.tri(R)
    R[idx] <- runif(sum(idx), min = r_min, max = r_max)
    R[lower.tri(R)] <- t(R)[lower.tri(R)]
  }

  R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)

  Z_latent <- MASS::mvrnorm(
    n = n,
    mu = rep(0, n_iv),
    Sigma = R
  )

  cuts <- qnorm(1 - pz)
  Z_mat <- sweep(Z_latent, 2, cuts, FUN = ">") * 1
  Z_mat <- matrix(Z_mat, nrow = n, ncol = n_iv)
  colnames(Z_mat) <- paste0("Z", seq_len(n_iv))

  Z_mat
}

## ------------------------------------------------------------
## Table 1 design parameters
## ------------------------------------------------------------

make_alpha_z <- function(iv_strength, n_iv = 30) {
  if (!iv_strength %in% c("weak", "moderate", "strong", "very_strong")) {
    stop("iv_strength must be weak, moderate, strong, or very_strong.")
  }

  switch(
    iv_strength,
    weak        = seq(0.001, 0.005, length.out = n_iv),
    moderate    = seq(0.4, 0.8, length.out = n_iv),
    strong      = rep(1.2, n_iv),
    very_strong = rep(3.0, n_iv)
  )
}

make_pz <- function(n_iv = 30) {
  if (n_iv != 30) {
    stop("make_pz() is only for the main revised Table 1 design with exactly 30 IVs. For NumIV sensitivity analyses, pass pz manually, e.g., rep(0.3, n_iv).")
  }

  c(rep(0.1, 10), rep(0.2, 10), rep(0.3, 10))
}

make_scenario_grid <- function() {
  iv_strength_levels <- c("weak", "moderate", "strong", "very_strong")
  conf_levels <- c(0.5, 1.0, 1.5)
  n_grid <- c(500, 1500, 2500, 5000, 7500, 10000)

  bind_rows(
    crossing(
      scenario = "A",
      b1 = 0,
      n = 1000,
      iv_strength = iv_strength_levels,
      c_x = conf_levels,
      c_y = conf_levels
    ) %>% filter(c_x == c_y),

    crossing(
      scenario = "B",
      b1 = c(1, 2, 3),
      n = 1000,
      iv_strength = iv_strength_levels,
      c_x = conf_levels,
      c_y = conf_levels
    ) %>% filter(c_x == c_y),

    crossing(
      scenario = "C",
      b1 = 0,
      n = n_grid,
      iv_strength = "weak",
      c_x = conf_levels,
      c_y = conf_levels
    ) %>% filter(c_x == c_y),

    crossing(
      scenario = "D",
      b1 = c(1, 2, 3),
      n = n_grid,
      iv_strength = "weak",
      c_x = conf_levels,
      c_y = conf_levels
    ) %>% filter(c_x == c_y)
  ) %>%
    mutate(
      setting_id = row_number(),
      n_iv = 30,
      alpha0 = 0,
      beta0 = 0,
      sigma_u = 0.5,
      alpha_min = map_dbl(iv_strength, ~ min(make_alpha_z(.x, n_iv = 30))),
      alpha_max = map_dbl(iv_strength, ~ max(make_alpha_z(.x, n_iv = 30)))
    ) %>%
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
      sigma_u
    )
}


sim_once <- function(n,
                     b1,
                     iv_strength,
                     c_x,
                     c_y,
                     n_iv = 30,
                     pz = make_pz(n_iv),
                     sigma_u = 0.5,
                     alpha0 = 0,
                     beta0 = 0,
                     r_min = 0.10,
                     r_max = 0.20,
                     seed = NULL,
                     ci_alpha = 0.05,
                     maxit = 1000) {

  if (!is.null(seed)) set.seed(seed)

  alpha_z <- make_alpha_z(iv_strength = iv_strength, n_iv = n_iv)

  if (length(pz) == 1) pz <- rep(pz, n_iv)
  if (length(pz) != n_iv) stop("Length of pz must be 1 or equal to n_iv.")
  if (length(alpha_z) != n_iv) stop("Length of alpha_z must equal n_iv.")

  Z_mat <- sim_binary_iv(
    n = n,
    pz = pz,
    r_min = r_min,
    r_max = r_max
  )

  U <- rnorm(n, mean = 0, sd = sigma_u)

  lin_x <- alpha0 + as.numeric(Z_mat %*% alpha_z) + c_x * U
  p_x <- trim_prob(expit(lin_x))
  X <- rbinom(n, size = 1, prob = p_x)

  lin_y <- beta0 + b1 * X + c_y * U
  p_y <- trim_prob(expit(lin_y))
  Y <- rbinom(n, size = 1, prob = p_y)

  dat <- data.frame(
    X = X,
    Y = Y,
    U = U,
    Z_mat,
    check.names = FALSE
  )

  z_names <- colnames(Z_mat)
  mr_dat <- dat[, c(z_names, "X", "Y"), drop = FALSE]

  fit <- tryCatch(
    MRbinary(
      dat = mr_dat,
      z = z_names,
      x = "X",
      y = "Y",
      alpha = ci_alpha,
      maxit = maxit
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    out <- tibble(
      method = c("2SPS", "2SRI", "GMM", "IV-MVB"),
      est = NA_real_,
      se = NA_real_,
      lcl = NA_real_,
      ucl = NA_real_,
      z = NA_real_,
      p = NA_real_,
      conv = 0,
      R2_mcf = NA_real_,
      R2_cs = NA_real_,
      R2_nk = NA_real_,
      F_mcf = NA_real_,
      F_cs = NA_real_,
      F_nk = NA_real_,
      F_anova = NA_real_
    )
  } else {
    out <- as_tibble(fit$summary)
  }

  out %>%
    mutate(
      n = n,
      n_iv = n_iv,
      b1 = b1,
      iv_strength = iv_strength,
      alpha_min = min(alpha_z),
      alpha_max = max(alpha_z),
      c_x = c_x,
      c_y = c_y,
      alpha0 = alpha0,
      beta0 = beta0,
      sigma_u = sigma_u,
      r_min = r_min,
      r_max = r_max,
      ci_alpha = ci_alpha,
      mean_X = mean(X),
      mean_Y = mean(Y),
      sd_X = sd(X),
      sd_Y = sd(Y),
      bias = est - b1,
      cover = ifelse(
        conv == 1 & is.finite(lcl) & is.finite(ucl),
        as.numeric(lcl <= b1 & ucl >= b1),
        NA_real_
      ),
      reject = ifelse(
        conv == 1 & is.finite(lcl) & is.finite(ucl),
        as.numeric(ucl < 0 | lcl > 0),
        NA_real_
      )
    )
}


run_one_setting <- function(setting_row,
                            n_sim = 8000,
                            seed_base = 2026,
                            out_dir = NULL,
                            ci_alpha = 0.05,
                            maxit = 1000) {

  setting_row <- as_tibble(setting_row)
  if (nrow(setting_row) != 1) stop("setting_row must contain exactly one row.")

  setting_id <- setting_row$setting_id

  res <- map_dfr(seq_len(n_sim), function(sim_id) {
    seed <- seed_base + setting_id * 100000 + sim_id

    sim_once(
      n = setting_row$n,
      b1 = setting_row$b1,
      iv_strength = setting_row$iv_strength,
      c_x = setting_row$c_x,
      c_y = setting_row$c_y,
      n_iv = setting_row$n_iv,
      pz = make_pz(setting_row$n_iv),
      sigma_u = setting_row$sigma_u,
      alpha0 = setting_row$alpha0,
      beta0 = setting_row$beta0,
      seed = seed,
      ci_alpha = ci_alpha,
      maxit = maxit
    ) %>%
      mutate(
        setting_id = setting_id,
        scenario = setting_row$scenario,
        sim_id = sim_id,
        seed = seed,
        .before = 1
      )
  })

  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    out_file <- file.path(
      out_dir,
      sprintf("sim_setting%03d_scenario%s.rds", setting_id, setting_row$scenario)
    )
    saveRDS(res, out_file)
  }

  res
}

run_sim <- function(n_sim,
                    n,
                    b1,
                    iv_strength,
                    c_x,
                    c_y,
                    n_iv = 30,
                    pz = make_pz(n_iv),
                    sigma_u = 0.5,
                    alpha0 = 0,
                    beta0 = 0,
                    r_min = 0.10,
                    r_max = 0.20,
                    scenario = NA_character_,
                    setting_id = NA_integer_,
                    seed_base = 2026,
                    ci_alpha = 0.05,
                    maxit = 1000) {

  map_dfr(seq_len(n_sim), function(sim_id) {
    seed <- seed_base + sim_id

    sim_once(
      n = n,
      b1 = b1,
      iv_strength = iv_strength,
      c_x = c_x,
      c_y = c_y,
      n_iv = n_iv,
      pz = pz,
      sigma_u = sigma_u,
      alpha0 = alpha0,
      beta0 = beta0,
      r_min = r_min,
      r_max = r_max,
      seed = seed,
      ci_alpha = ci_alpha,
      maxit = maxit
    ) %>%
      mutate(
        setting_id = setting_id,
        scenario = scenario,
        sim_id = sim_id,
        seed = seed,
        .before = 1
      )
  })
}

run_all_settings <- function(n_sim = 8000,
                             settings = make_scenario_grid(),
                             seed_base = 2026,
                             out_dir = "results",
                             ci_alpha = 0.05,
                             maxit = 1000) {

  if (nrow(settings) != 120) {
    stop("Expected 120 scenario settings from Table 1, but found ", nrow(settings), ".")
  }

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  map_dfr(seq_len(nrow(settings)), function(i) {
    message(
      "Running setting ", i, " of ", nrow(settings),
      " | Scenario ", settings$scenario[i],
      " | b1=", settings$b1[i],
      " | n=", settings$n[i],
      " | IV=", settings$iv_strength[i],
      " | c_x=c_y=", settings$c_x[i]
    )

    run_one_setting(
      setting_row = settings[i, ],
      n_sim = n_sim,
      seed_base = seed_base,
      out_dir = out_dir,
      ci_alpha = ci_alpha,
      maxit = maxit
    )
  })
}


summarize_sim_results <- function(sim_res) {
  sim_res %>%
    group_by(
      scenario,
      setting_id,
      method,
      n,
      n_iv,
      b1,
      iv_strength,
      alpha_min,
      alpha_max,
      c_x,
      c_y
    ) %>%
    summarise(
      n_sim = n_distinct(sim_id),
      n_conv = sum(conv == 1, na.rm = TRUE),

      median_bias = median(bias[conv == 1], na.rm = TRUE),
      q1_bias = quantile(bias[conv == 1], 0.25, na.rm = TRUE),
      q3_bias = quantile(bias[conv == 1], 0.75, na.rm = TRUE),

      coverage = mean(cover[conv == 1], na.rm = TRUE),
      rejection = mean(reject[conv == 1], na.rm = TRUE),
      type1 = ifelse(first(b1) == 0, rejection, NA_real_),
      power = ifelse(first(b1) != 0, rejection, NA_real_),

      median_R2_mcf = median(R2_mcf, na.rm = TRUE),
      median_R2_cs = median(R2_cs, na.rm = TRUE),
      median_R2_nk = median(R2_nk, na.rm = TRUE),

      median_F_mcf = median(F_mcf, na.rm = TRUE),
      q1_F_mcf = quantile(F_mcf, 0.25, na.rm = TRUE),
      q3_F_mcf = quantile(F_mcf, 0.75, na.rm = TRUE),

      median_F_cs = median(F_cs, na.rm = TRUE),
      median_F_nk = median(F_nk, na.rm = TRUE),
      median_F_anova = median(F_anova, na.rm = TRUE),

      mean_X = mean(mean_X, na.rm = TRUE),
      mean_Y = mean(mean_Y, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      bias_text = sprintf("%.2f (%.2f, %.2f)", median_bias, q1_bias, q3_bias),
      coverage_pct = 100 * coverage,
      type1_pct = 100 * type1,
      power_pct = 100 * power,
      F_mcf_text = sprintf("%.2f (%.2f, %.2f)", median_F_mcf, q1_F_mcf, q3_F_mcf),
      F_cat_mcf = case_when(
        median_F_mcf < 10 ~ "Weak",
        median_F_mcf >= 10 & median_F_mcf <= 30 ~ "Moderate",
        median_F_mcf > 30 & median_F_mcf <= 100 ~ "Strong",
        median_F_mcf > 100 ~ "Very strong",
        TRUE ~ NA_character_
      )
    )
}

save_sim_outputs <- function(sim_res,
                             out_dir = "results",
                             prefix = "Paper3_sim") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  sim_sum <- summarize_sim_results(sim_res)

  write_csv(sim_res, file.path(out_dir, paste0(prefix, "_raw.csv")))
  write_csv(sim_sum, file.path(out_dir, paste0(prefix, "_summary.csv")))
  saveRDS(sim_res, file.path(out_dir, paste0(prefix, "_raw.rds")))
  saveRDS(sim_sum, file.path(out_dir, paste0(prefix, "_summary.rds")))

  invisible(list(raw = sim_res, summary = sim_sum))
}
