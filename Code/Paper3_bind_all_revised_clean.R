library(dplyr)
library(purrr)
library(readr)
library(tibble)

source("scripts/Paper3_scenarios_revised_clean.R")

dir.create("results_combined", showWarnings = FALSE)

all_files <- list.files(
  path = "results",
  pattern = "^paper3_setting[0-9]+_sim[0-9]+\\.rds$",
  full.names = TRUE
)

cat("Total simulation result files found:", length(all_files), "\n")

if (length(all_files) != n_jobs_final) {
  stop("Expected ", n_jobs_final, " result files, but found ", length(all_files))
}

sim_res_all <- map_dfr(all_files, readRDS)

sim_res_all <- sim_res_all %>%
  mutate(
    method = factor(
      as.character(method),
      levels = methods_main
    ),
    bias = est - b1,
    cover = ifelse(
      conv == 1 & is.finite(lcl) & is.finite(ucl),
      as.numeric(lcl <= b1 & ucl >= b1),
      NA_real_
    ),
    reject = ifelse(
      conv == 1 & is.finite(p),
      as.numeric(p < ci_alpha),
      NA_real_
    )
  )

cat("\nRows by scenario, method, and convergence:\n")
print(
  sim_res_all %>%
    count(scenario, method, conv, name = "n_rows") %>%
    arrange(scenario, method, conv)
)

cat("\nRows by setting:\n")
setting_counts <- sim_res_all %>%
  count(setting_id, scenario, b1, n, iv_strength, c_x, c_y, name = "n_rows") %>%
  arrange(setting_id)

print(setting_counts)

bad_settings <- setting_counts %>%
  filter(n_rows != n_sim_final * length(methods_main))

if (nrow(bad_settings) > 0) {
  warning("Some settings do not have the expected number of rows:")
  print(bad_settings)
}

analysis_dat <- sim_res_all %>%
  filter(conv == 1, is.finite(bias))

summary_main <- analysis_dat %>%
  group_by(
    setting_id,
    scenario,
    b1,
    n,
    n_iv,
    iv_strength,
    c_x,
    c_y,
    method
  ) %>%
  summarise(
    n_total = n_sim_final,
    n_converged = n(),
    
    median_bias = median(bias, na.rm = TRUE),
    q1_bias = quantile(bias, 0.25, na.rm = TRUE),
    q3_bias = quantile(bias, 0.75, na.rm = TRUE),
    
    coverage = mean(cover, na.rm = TRUE),
    rejection = mean(reject, na.rm = TRUE),
    
    type1 = ifelse(first(b1) == 0, rejection, NA_real_),
    power = ifelse(first(b1) != 0, rejection, NA_real_),
    
    n_finite_R2_mcf = sum(is.finite(R2_mcf)),
    median_R2_mcf = median(R2_mcf, na.rm = TRUE),
    q1_R2_mcf = quantile(R2_mcf, 0.25, na.rm = TRUE),
    q3_R2_mcf = quantile(R2_mcf, 0.75, na.rm = TRUE),
    
    n_finite_R2_cs = sum(is.finite(R2_cs)),
    median_R2_cs = median(R2_cs, na.rm = TRUE),
    q1_R2_cs = quantile(R2_cs, 0.25, na.rm = TRUE),
    q3_R2_cs = quantile(R2_cs, 0.75, na.rm = TRUE),
    
    n_finite_R2_nk = sum(is.finite(R2_nk)),
    median_R2_nk = median(R2_nk, na.rm = TRUE),
    q1_R2_nk = quantile(R2_nk, 0.25, na.rm = TRUE),
    q3_R2_nk = quantile(R2_nk, 0.75, na.rm = TRUE),
    
    n_finite_F_mcf = sum(is.finite(F_mcf)),
    median_F_mcf = median(F_mcf, na.rm = TRUE),
    q1_F_mcf = quantile(F_mcf, 0.25, na.rm = TRUE),
    q3_F_mcf = quantile(F_mcf, 0.75, na.rm = TRUE),
    
    n_finite_F_cs = sum(is.finite(F_cs)),
    median_F_cs = median(F_cs, na.rm = TRUE),
    q1_F_cs = quantile(F_cs, 0.25, na.rm = TRUE),
    q3_F_cs = quantile(F_cs, 0.75, na.rm = TRUE),
    
    n_finite_F_nk = sum(is.finite(F_nk)),
    median_F_nk = median(F_nk, na.rm = TRUE),
    q1_F_nk = quantile(F_nk, 0.25, na.rm = TRUE),
    q3_F_nk = quantile(F_nk, 0.75, na.rm = TRUE),
    
    n_finite_F_anova = sum(is.finite(F_anova)),
    median_F_anova = median(F_anova, na.rm = TRUE),
    q1_F_anova = quantile(F_anova, 0.25, na.rm = TRUE),
    q3_F_anova = quantile(F_anova, 0.75, na.rm = TRUE),
    
    mean_X = mean(mean_X, na.rm = TRUE),
    sd_X = mean(sd_X, na.rm = TRUE),
    mean_Y = mean(mean_Y, na.rm = TRUE),
    sd_Y = mean(sd_Y, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    bias_text = sprintf("%.2f (%.2f, %.2f)", median_bias, q1_bias, q3_bias),
    coverage_text = sprintf("%.1f%%", 100 * coverage),
    type1_text = ifelse(is.finite(type1), sprintf("%.1f%%", 100 * type1), NA_character_),
    power_text = ifelse(is.finite(power), sprintf("%.1f%%", 100 * power), NA_character_),
    
    F_mcf_text = sprintf("%.2f (%.2f, %.2f)", median_F_mcf, q1_F_mcf, q3_F_mcf),
    F_cs_text = sprintf("%.2f (%.2f, %.2f)", median_F_cs, q1_F_cs, q3_F_cs),
    F_nk_text = sprintf("%.2f (%.2f, %.2f)", median_F_nk, q1_F_nk, q3_F_nk),
    F_anova_text = sprintf("%.2f (%.2f, %.2f)", median_F_anova, q1_F_anova, q3_F_anova),
    
    R2_mcf_text = sprintf("%.4f (%.4f, %.4f)", median_R2_mcf, q1_R2_mcf, q3_R2_mcf),
    R2_cs_text = sprintf("%.4f (%.4f, %.4f)", median_R2_cs, q1_R2_cs, q3_R2_cs),
    R2_nk_text = sprintf("%.4f (%.4f, %.4f)", median_R2_nk, q1_R2_nk, q3_R2_nk),
    
    F_cat_mcf = case_when(
      median_F_mcf < 10 ~ "Weak",
      median_F_mcf >= 10 & median_F_mcf <= 30 ~ "Moderate",
      median_F_mcf > 30 & median_F_mcf <= 100 ~ "Strong",
      median_F_mcf > 100 ~ "Very strong",
      TRUE ~ NA_character_
    ),
    
    F_cat_cs = case_when(
      median_F_cs < 10 ~ "Weak",
      median_F_cs >= 10 & median_F_cs <= 30 ~ "Moderate",
      median_F_cs > 30 & median_F_cs <= 100 ~ "Strong",
      median_F_cs > 100 ~ "Very strong",
      TRUE ~ NA_character_
    ),
    
    F_cat_nk = case_when(
      median_F_nk < 10 ~ "Weak",
      median_F_nk >= 10 & median_F_nk <= 30 ~ "Moderate",
      median_F_nk > 30 & median_F_nk <= 100 ~ "Strong",
      median_F_nk > 100 ~ "Very strong",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(setting_id, method)

cat("\nMain simulation summary:\n")
print(summary_main)

write_csv(sim_res_all, "results_combined/Paper3_allmethods_30IV_raw.csv")
write_csv(summary_main, "results_combined/Paper3_allmethods_30IV_summary.csv")

saveRDS(sim_res_all, "results_combined/Paper3_allmethods_30IV_raw.rds")
saveRDS(summary_main, "results_combined/Paper3_allmethods_30IV_summary.rds")

cat("\nSaved:\n")
cat("results_combined/Paper3_allmethods_30IV_raw.csv\n")
cat("results_combined/Paper3_allmethods_30IV_summary.csv\n")
cat("results_combined/Paper3_allmethods_30IV_raw.rds\n")
cat("results_combined/Paper3_allmethods_30IV_summary.rds\n")