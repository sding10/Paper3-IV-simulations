library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

options(bitmapType = "cairo")

dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

sim_res_all <- readRDS("results_withF_newStrength_all_scenarios_30IV_100sim.rds")

method_levels <- c("2SPS", "2SRI", "GMM", "IV-MVB")

strength_levels <- c(
  "Very Weak IVs",
  "Weak IVs",
  "Moderately Weak IVs",
  "Moderate IVs",
  "Moderately Strong IVs"
)

method_cols <- c(
  "2SPS" = "grey90",
  "2SRI" = "grey75",
  "GMM" = "grey50",
  "IV-MVB" = "grey15"
)

sim_res_all <- sim_res_all %>%
  mutate(
    scen = factor(as.character(scen), levels = c("A", "B", "C", "D")),
    method = factor(as.character(method), levels = method_levels),
    strength_label = factor(as.character(strength_label), levels = strength_levels),
    bias = est - b1,
    cover = ifelse(
      conv == 1 & is.finite(lcl) & is.finite(ucl),
      as.numeric(lcl <= b1 & ucl >= b1),
      NA_real_
    ),
    reject = ifelse(
      conv == 1 & is.finite(p),
      as.numeric(p < 0.05),
      NA_real_
    )
  )

## ============================================================
## First-stage F summaries
## ============================================================

f_summary_ab <- sim_res_all %>%
  filter(scen %in% c("A", "B")) %>%
  distinct(scen, strength_label, cell_id, cell_within_group, sim, first_stage_F) %>%
  group_by(scen, strength_label) %>%
  summarise(
    median_F = median(first_stage_F, na.rm = TRUE),
    q1_F = quantile(first_stage_F, 0.25, na.rm = TRUE),
    q3_F = quantile(first_stage_F, 0.75, na.rm = TRUE),
    min_F = min(first_stage_F, na.rm = TRUE),
    max_F = max(first_stage_F, na.rm = TRUE),
    F_text = sprintf("%.2f (%.2f, %.2f)", median_F, q1_F, q3_F),
    .groups = "drop"
  ) %>%
  arrange(scen, strength_label)

f_summary_cd <- sim_res_all %>%
  filter(scen %in% c("C", "D")) %>%
  distinct(scen, n, cell_id, cell_within_group, sim, first_stage_F) %>%
  group_by(scen, n) %>%
  summarise(
    median_F = median(first_stage_F, na.rm = TRUE),
    q1_F = quantile(first_stage_F, 0.25, na.rm = TRUE),
    q3_F = quantile(first_stage_F, 0.75, na.rm = TRUE),
    min_F = min(first_stage_F, na.rm = TRUE),
    max_F = max(first_stage_F, na.rm = TRUE),
    F_text = sprintf("%.2f (%.2f, %.2f)", median_F, q1_F, q3_F),
    .groups = "drop"
  ) %>%
  arrange(scen, n)

write_csv(f_summary_ab, "tables/FirstStageF_Summary_AB_30IV_withF_newStrength.csv")
write_csv(f_summary_cd, "tables/FirstStageF_Summary_CD_30IV_withF_newStrength.csv")

## ============================================================
## Figure 1
## ============================================================

fig1_dat <- sim_res_all %>%
  filter(conv == 1, is.finite(bias)) %>%
  mutate(
    panel_x = case_when(
      scen %in% c("A", "B") ~ as.character(strength_label),
      scen %in% c("C", "D") ~ as.character(n),
      TRUE ~ NA_character_
    ),
    panel_x = factor(
      panel_x,
      levels = c(strength_levels, "500")
    )
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
  facet_wrap(~ scen, nrow = 2, scales = "free_x") +
  scale_fill_manual(values = method_cols, drop = FALSE) +
  scale_x_discrete(drop = TRUE) +
  coord_cartesian(ylim = c(-10, 10)) +
  labs(
    title = "Figure 1. Bias estimates for each MR method with fifteen instruments.",
    subtitle = "withF new-strength test: A/B vary IV strength at n = 1000; C/D are evaluated at n = 500.",
    x = NULL,
    y = "Estimated Bias in Causal Parameter",
    fill = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16, family = "serif"),
    plot.subtitle = element_text(size = 12, family = "serif"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  "figures/Figure1_bias_distribution_30IV_withF_newStrength.png",
  fig1,
  width = 12,
  height = 8.5,
  dpi = 300,
  device = "png"
)

ggsave(
  "figures/Figure1_bias_distribution_30IV_withF_newStrength.pdf",
  fig1,
  width = 12,
  height = 8.5,
  device = cairo_pdf
)

## ============================================================
## Supplementary Table 1: Bias + F diagnostic
## ============================================================

bias_summary <- sim_res_all %>%
  filter(conv == 1, is.finite(bias)) %>%
  group_by(scen, strength_label, n, method) %>%
  summarise(
    num_converged = n(),
    median_bias = median(bias, na.rm = TRUE),
    q1_bias = quantile(bias, 0.25, na.rm = TRUE),
    q3_bias = quantile(bias, 0.75, na.rm = TRUE),
    bias_text = sprintf("%.2f (%.2f, %.2f)", median_bias, q1_bias, q3_bias),
    .groups = "drop"
  )

table1_ab <- bias_summary %>%
  filter(scen %in% c("A", "B")) %>%
  left_join(
    f_summary_ab %>% select(scen, strength_label, F_text),
    by = c("scen", "strength_label")
  ) %>%
  select(scen, strength_label, F_text, method, num_converged, bias_text) %>%
  pivot_wider(
    names_from = scen,
    values_from = c(F_text, num_converged, bias_text),
    names_glue = "{scen}_{.value}"
  ) %>%
  arrange(strength_label, method) %>%
  rename(
    `Instrument Strengths` = strength_label,
    `MR Method` = method,
    `Scenario A F Median (Q1, Q3)` = A_F_text,
    `Scenario A Num Converged` = A_num_converged,
    `Scenario A Median (Q1, Q3) Bias` = A_bias_text,
    `Scenario B F Median (Q1, Q3)` = B_F_text,
    `Scenario B Num Converged` = B_num_converged,
    `Scenario B Median (Q1, Q3) Bias` = B_bias_text
  )

table1_cd <- bias_summary %>%
  filter(scen %in% c("C", "D")) %>%
  left_join(
    f_summary_cd %>% select(scen, n, F_text),
    by = c("scen", "n")
  ) %>%
  select(scen, n, F_text, method, num_converged, bias_text) %>%
  pivot_wider(
    names_from = scen,
    values_from = c(F_text, num_converged, bias_text),
    names_glue = "{scen}_{.value}"
  ) %>%
  arrange(n, method) %>%
  rename(
    `Sample Size` = n,
    `MR Method` = method,
    `Scenario C F Median (Q1, Q3)` = C_F_text,
    `Scenario C Num Converged` = C_num_converged,
    `Scenario C Median (Q1, Q3) Bias` = C_bias_text,
    `Scenario D F Median (Q1, Q3)` = D_F_text,
    `Scenario D Num Converged` = D_num_converged,
    `Scenario D Median (Q1, Q3) Bias` = D_bias_text
  )

write_csv(table1_ab, "tables/Supplementary_Table1A_Bias_F_AB_30IV_withF_newStrength.csv")
write_csv(table1_cd, "tables/Supplementary_Table1B_Bias_F_CD_30IV_withF_newStrength.csv")

## ============================================================
## Supplementary Table 2: Coverage + Type I error / Power + F
## ============================================================

cov_summary <- sim_res_all %>%
  group_by(scen, strength_label, n, method) %>%
  summarise(
    num_converged = sum(conv == 1, na.rm = TRUE),
    coverage = mean(cover, na.rm = TRUE),
    rejection_rate = mean(reject, na.rm = TRUE),
    coverage_text = sprintf("%.1f%%", 100 * coverage),
    rejection_text = sprintf("%.1f%%", 100 * rejection_rate),
    .groups = "drop"
  )

table2_ab <- cov_summary %>%
  filter(scen %in% c("A", "B")) %>%
  left_join(
    f_summary_ab %>% select(scen, strength_label, F_text),
    by = c("scen", "strength_label")
  ) %>%
  select(scen, strength_label, F_text, method, num_converged, coverage_text, rejection_text) %>%
  pivot_wider(
    names_from = scen,
    values_from = c(F_text, num_converged, coverage_text, rejection_text),
    names_glue = "{scen}_{.value}"
  ) %>%
  arrange(strength_label, method) %>%
  rename(
    `Instrument Strengths` = strength_label,
    `MR Method` = method,
    `Scenario A F Median (Q1, Q3)` = A_F_text,
    `Scenario A Num Converged` = A_num_converged,
    `Scenario A Coverage` = A_coverage_text,
    `Scenario A Type I Error` = A_rejection_text,
    `Scenario B F Median (Q1, Q3)` = B_F_text,
    `Scenario B Num Converged` = B_num_converged,
    `Scenario B Coverage` = B_coverage_text,
    `Scenario B Power` = B_rejection_text
  )

table2_cd <- cov_summary %>%
  filter(scen %in% c("C", "D")) %>%
  left_join(
    f_summary_cd %>% select(scen, n, F_text),
    by = c("scen", "n")
  ) %>%
  select(scen, n, F_text, method, num_converged, coverage_text, rejection_text) %>%
  pivot_wider(
    names_from = scen,
    values_from = c(F_text, num_converged, coverage_text, rejection_text),
    names_glue = "{scen}_{.value}"
  ) %>%
  arrange(n, method) %>%
  rename(
    `Sample Size` = n,
    `MR Method` = method,
    `Scenario C F Median (Q1, Q3)` = C_F_text,
    `Scenario C Num Converged` = C_num_converged,
    `Scenario C Coverage` = C_coverage_text,
    `Scenario C Type I Error` = C_rejection_text,
    `Scenario D F Median (Q1, Q3)` = D_F_text,
    `Scenario D Num Converged` = D_num_converged,
    `Scenario D Coverage` = D_coverage_text,
    `Scenario D Power` = D_rejection_text
  )

write_csv(table2_ab, "tables/Supplementary_Table2A_Coverage_TypeI_Power_F_AB_30IV_withF_newStrength.csv")
write_csv(table2_cd, "tables/Supplementary_Table2B_Coverage_TypeI_Power_F_CD_30IV_withF_newStrength.csv")

check_counts <- sim_res_all %>%
  count(scen, n, strength_label, method, conv)

write_csv(check_counts, "tables/check_counts_30IV_withF_newStrength.csv")

message("Done. Figures and tables saved for 30IV withF new-strength test.")