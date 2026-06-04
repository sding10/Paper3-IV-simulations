library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

options(bitmapType = "cairo")

dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)


raw_rds <- "results_combined/Paper3_allmethods_30IV_raw.rds"
raw_csv <- "results_combined/Paper3_allmethods_30IV_raw.csv"

if (file.exists(raw_rds)) {
  sim_res_all <- readRDS(raw_rds)
} else if (file.exists(raw_csv)) {
  sim_res_all <- read_csv(raw_csv, show_col_types = FALSE)
} else {
  stop("Cannot find combined raw simulation results. Expected: ", raw_rds,
       " or ", raw_csv)
}

method_levels <- c("2SPS", "2SRI", "GMM", "IV-MVB")

iv_strength_levels <- c(
  "weak",
  "moderate",
  "strong",
  "very_strong"
)

iv_strength_labels <- c(
  weak = "Weak IVs",
  moderate = "Moderate IVs",
  strong = "Strong IVs",
  very_strong = "Very Strong IVs"
)

method_cols <- c(
  "2SPS" = "grey90",
  "2SRI" = "grey75",
  "GMM" = "grey50",
  "IV-MVB" = "grey15"
)

if (!("ci_alpha" %in% names(sim_res_all))) {
  sim_res_all$ci_alpha <- 0.05
}

sim_res_all <- sim_res_all %>%
  mutate(
    scenario = factor(as.character(scenario), levels = c("A", "B", "C", "D")),
    method = factor(as.character(method), levels = method_levels),
    iv_strength = factor(as.character(iv_strength), levels = iv_strength_levels),
    iv_strength_label = factor(
      iv_strength_labels[as.character(iv_strength)],
      levels = unname(iv_strength_labels)
    ),
    bias = ifelse(is.finite(est), est - b1, NA_real_),
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


cat("\nScenario/sample-size check:\n")
print(sim_res_all %>% count(scenario, n, name = "n_rows"))

cat("\nA/B IV-strength check:\n")
print(
  sim_res_all %>%
    filter(scenario %in% c("A", "B")) %>%
    count(scenario, n, iv_strength_label, c_x, c_y, name = "n_rows")
)

cat("\nC/D sample-size check:\n")
print(
  sim_res_all %>%
    filter(scenario %in% c("C", "D")) %>%
    count(scenario, n, c_x, c_y, name = "n_rows")
)

sim_group_check <- sim_res_all %>%
  group_by(setting_id, scenario, b1, n, iv_strength_label, c_x, c_y, method) %>%
  summarise(
    n_rows = n(),
    n_sim = n_distinct(sim_id),
    n_converged = sum(conv == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(setting_id, method)

write_csv(sim_group_check, "tables/check_sim_group_counts_30IV.csv")

## ------------------------------------------------------------
## Figure 1: Bias distributions
## ------------------------------------------------------------

fig1_dat <- sim_res_all %>%
  filter(conv == 1, is.finite(bias)) %>%
  mutate(
    panel_x = case_when(
      scenario %in% c("A", "B") ~ as.character(iv_strength_label),
      scenario %in% c("C", "D") ~ as.character(n),
      TRUE ~ NA_character_
    ),
    panel_x = factor(
      panel_x,
      levels = c(unname(iv_strength_labels), as.character(sort(unique(n))))
    ),
    conf_label = paste0("c_x=c_y=", c_x),
    b1_label = paste0("b1=", b1)
  )

fig1 <- ggplot(fig1_dat, aes(x = panel_x, y = bias, fill = method)) +
  geom_hline(yintercept = 0, color = "darkgreen", linewidth = 0.30) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.42,
    linewidth = 0.25,
    outlier.size = 0.10,
    outlier.alpha = 0.12,
    color = "black"
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.30,
    linewidth = 0.22,
    color = "red",
    position = position_dodge(width = 0.75),
    show.legend = FALSE
  ) +
  facet_grid(conf_label + b1_label ~ scenario, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = method_cols, drop = FALSE) +
  scale_x_discrete(drop = TRUE) +
  coord_cartesian(ylim = c(-10, 10)) +
  labs(
    title = "Figure 1. Full distribution of bias estimates for each MR method from the revised simulation scenarios.",
    subtitle = "Panels separate confounding strength and true causal effect to avoid pooling across scenario settings.",
    x = NULL,
    y = "Estimated Bias in Causal Parameter",
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15, family = "serif"),
    plot.subtitle = element_text(size = 11, family = "serif"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 7, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 9),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  "figures/Figure1_bias_distribution_30IV.png",
  fig1,
  width = 16,
  height = 11,
  dpi = 300,
  device = "png"
)

ggsave(
  "figures/Figure1_bias_distribution_30IV.pdf",
  fig1,
  width = 16,
  height = 11,
  device = cairo_pdf
)

## ------------------------------------------------------------
## Figure 2: Type I error and power
## ------------------------------------------------------------

fig2_dat <- sim_res_all %>%
  filter(conv == 1, is.finite(reject)) %>%
  group_by(scenario, b1, n, iv_strength_label, c_x, c_y, method) %>%
  summarise(
    rejection_rate = mean(reject, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    estimand = ifelse(b1 == 0, "Type I error", "Power"),
    panel_x = case_when(
      scenario %in% c("A", "B") ~ as.character(iv_strength_label),
      scenario %in% c("C", "D") ~ as.character(n),
      TRUE ~ NA_character_
    ),
    panel_x = factor(
      panel_x,
      levels = c(unname(iv_strength_labels), as.character(sort(unique(n))))
    ),
    conf_label = paste0("c_x=c_y=", c_x),
    b1_label = paste0("b1=", b1)
  )

fig2 <- ggplot(fig2_dat, aes(x = panel_x, y = rejection_rate, fill = method)) +
  geom_hline(yintercept = 0.05, color = "darkgreen", linewidth = 0.30, linetype = "dashed") +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black",
    linewidth = 0.20
  ) +
  facet_grid(conf_label + estimand + b1_label ~ scenario, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = method_cols, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accurac_y = 1)) +
  labs(
    title = "Figure 2. Rejection rates for each MR method across the revised simulation scenarios.",
    subtitle = "For null scenarios, rejection rate is type I error; for positive-effect scenarios, rejection rate is power.",
    x = NULL,
    y = "Rejection Rate",
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15, family = "serif"),
    plot.subtitle = element_text(size = 11, family = "serif"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 8),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 7, angle = 30, hjust = 1),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  "figures/Figure2_typeI_power_30IV.png",
  fig2,
  width = 16,
  height = 12,
  dpi = 300,
  device = "png"
)

ggsave(
  "figures/Figure2_typeI_power_30IV.pdf",
  fig2,
  width = 16,
  height = 12,
  device = cairo_pdf
)

## ------------------------------------------------------------
## Supplementary Table 1: Bias summary
## ------------------------------------------------------------

bias_summary <- sim_res_all %>%
  filter(conv == 1, is.finite(bias)) %>%
  group_by(scenario, b1, n, iv_strength_label, c_x, c_y, method) %>%
  summarise(
    num_converged = n(),
    median_bias = median(bias, na.rm = TRUE),
    q1_bias = quantile(bias, 0.25, na.rm = TRUE),
    q3_bias = quantile(bias, 0.75, na.rm = TRUE),
    bias_text = sprintf("%.2f (%.2f, %.2f)", median_bias, q1_bias, q3_bias),
    .groups = "drop"
  )

bias_ab <- bias_summary %>%
  filter(scenario %in% c("A", "B")) %>%
  arrange(scenario, c_x, b1, iv_strength_label, method)

bias_cd <- bias_summary %>%
  filter(scenario %in% c("C", "D")) %>%
  arrange(scenario, c_x, b1, n, method)

write_csv(bias_ab, "tables/Supplementary_Table1A_Bias_AB_30IV.csv")
write_csv(bias_cd, "tables/Supplementary_Table1B_Bias_CD_30IV.csv")

bias_ab_wide <- bias_summary %>%
  filter(scenario %in% c("A", "B")) %>%
  select(scenario, b1, c_x, iv_strength_label, method, num_converged, bias_text) %>%
  pivot_wider(
    names_from = scenario,
    values_from = c(num_converged, bias_text),
    names_glue = "Scenario_{scenario}_{.value}"
  ) %>%
  arrange(c_x, b1, iv_strength_label, method)

bias_cd_wide <- bias_summary %>%
  filter(scenario %in% c("C", "D")) %>%
  select(scenario, b1, c_x, n, method, num_converged, bias_text) %>%
  pivot_wider(
    names_from = scenario,
    values_from = c(num_converged, bias_text),
    names_glue = "Scenario_{scenario}_{.value}"
  ) %>%
  arrange(c_x, b1, n, method)

write_csv(bias_ab_wide, "tables/Supplementary_Table1A_Bias_AB_30IV_wide.csv")
write_csv(bias_cd_wide, "tables/Supplementary_Table1B_Bias_CD_30IV_wide.csv")

## ------------------------------------------------------------
## Supplementary Table 2: Coverage, type I error, and power
## ------------------------------------------------------------

cov_summary <- sim_res_all %>%
  group_by(scenario, b1, n, iv_strength_label, c_x, c_y, method) %>%
  summarise(
    num_converged = sum(conv == 1, na.rm = TRUE),
    coverage = mean(cover[conv == 1], na.rm = TRUE),
    rejection_rate = mean(reject[conv == 1], na.rm = TRUE),
    type1 = ifelse(first(b1) == 0, rejection_rate, NA_real_),
    power = ifelse(first(b1) != 0, rejection_rate, NA_real_),
    coverage_text = sprintf("%.1f%%", 100 * coverage),
    type1_text = ifelse(is.finite(type1), sprintf("%.1f%%", 100 * type1), NA_character_),
    power_text = ifelse(is.finite(power), sprintf("%.1f%%", 100 * power), NA_character_),
    .groups = "drop"
  )

cov_ab <- cov_summary %>%
  filter(scenario %in% c("A", "B")) %>%
  arrange(scenario, c_x, b1, iv_strength_label, method)

cov_cd <- cov_summary %>%
  filter(scenario %in% c("C", "D")) %>%
  arrange(scenario, c_x, b1, n, method)

write_csv(cov_ab, "tables/Supplementary_Table2A_Coverage_TypeI_Power_AB_30IV.csv")
write_csv(cov_cd, "tables/Supplementary_Table2B_Coverage_TypeI_Power_CD_30IV.csv")

cov_ab_wide <- cov_summary %>%
  filter(scenario %in% c("A", "B")) %>%
  select(scenario, b1, c_x, iv_strength_label, method,
         num_converged, coverage_text, type1_text, power_text) %>%
  pivot_wider(
    names_from = scenario,
    values_from = c(num_converged, coverage_text, type1_text, power_text),
    names_glue = "Scenario_{scenario}_{.value}"
  ) %>%
  arrange(c_x, b1, iv_strength_label, method)

cov_cd_wide <- cov_summary %>%
  filter(scenario %in% c("C", "D")) %>%
  select(scenario, b1, c_x, n, method,
         num_converged, coverage_text, type1_text, power_text) %>%
  pivot_wider(
    names_from = scenario,
    values_from = c(num_converged, coverage_text, type1_text, power_text),
    names_glue = "Scenario_{scenario}_{.value}"
  ) %>%
  arrange(c_x, b1, n, method)

write_csv(cov_ab_wide, "tables/Supplementary_Table2A_Coverage_TypeI_Power_AB_30IV_wide.csv")
write_csv(cov_cd_wide, "tables/Supplementary_Table2B_Coverage_TypeI_Power_CD_30IV_wide.csv")

## ------------------------------------------------------------
## Supplementary Table 3: First-stage F and pseudo-R2 summaries
## ------------------------------------------------------------

f_dat <- sim_res_all %>%
  distinct(
    setting_id, sim_id, scenario, b1, n, n_iv,
    iv_strength_label, c_x, c_y,
    .keep_all = TRUE
  )

f_summary <- f_dat %>%
  group_by(setting_id, scenario, b1, n, n_iv, iv_strength_label, c_x, c_y) %>%
  summarise(
    n_sim = n_distinct(sim_id),

    median_R2_mcf = median(R2_mcf, na.rm = TRUE),
    q1_R2_mcf = quantile(R2_mcf, 0.25, na.rm = TRUE),
    q3_R2_mcf = quantile(R2_mcf, 0.75, na.rm = TRUE),

    median_R2_cs = median(R2_cs, na.rm = TRUE),
    q1_R2_cs = quantile(R2_cs, 0.25, na.rm = TRUE),
    q3_R2_cs = quantile(R2_cs, 0.75, na.rm = TRUE),

    median_R2_nk = median(R2_nk, na.rm = TRUE),
    q1_R2_nk = quantile(R2_nk, 0.25, na.rm = TRUE),
    q3_R2_nk = quantile(R2_nk, 0.75, na.rm = TRUE),

    median_F_mcf = median(F_mcf, na.rm = TRUE),
    q1_F_mcf = quantile(F_mcf, 0.25, na.rm = TRUE),
    q3_F_mcf = quantile(F_mcf, 0.75, na.rm = TRUE),

    median_F_cs = median(F_cs, na.rm = TRUE),
    q1_F_cs = quantile(F_cs, 0.25, na.rm = TRUE),
    q3_F_cs = quantile(F_cs, 0.75, na.rm = TRUE),

    median_F_nk = median(F_nk, na.rm = TRUE),
    q1_F_nk = quantile(F_nk, 0.25, na.rm = TRUE),
    q3_F_nk = quantile(F_nk, 0.75, na.rm = TRUE),

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
    R2_mcf_text = sprintf("%.4f (%.4f, %.4f)", median_R2_mcf, q1_R2_mcf, q3_R2_mcf),
    R2_cs_text = sprintf("%.4f (%.4f, %.4f)", median_R2_cs, q1_R2_cs, q3_R2_cs),
    R2_nk_text = sprintf("%.4f (%.4f, %.4f)", median_R2_nk, q1_R2_nk, q3_R2_nk),
    F_mcf_text = sprintf("%.2f (%.2f, %.2f)", median_F_mcf, q1_F_mcf, q3_F_mcf),
    F_cs_text = sprintf("%.2f (%.2f, %.2f)", median_F_cs, q1_F_cs, q3_F_cs),
    F_nk_text = sprintf("%.2f (%.2f, %.2f)", median_F_nk, q1_F_nk, q3_F_nk),
    F_anova_text = sprintf("%.2f (%.2f, %.2f)", median_F_anova, q1_F_anova, q3_F_anova),
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
  arrange(setting_id)

write_csv(f_summary, "tables/Supplementary_Table3_FirstStage_F_R2_30IV.csv")


check_counts <- sim_res_all %>%
  count(scenario, b1, n, iv_strength_label, c_x, c_y, method, conv, name = "n_rows") %>%
  arrange(scenario, b1, n, iv_strength_label, c_x, method, conv)

write_csv(check_counts, "tables/check_counts_30IV.csv")

message("Done. Figures saved as PNG/PDF and tables saved as CSV.")
