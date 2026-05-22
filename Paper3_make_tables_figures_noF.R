library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

options(bitmapType = "cairo")

dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

sim_res_all <- readRDS("results_noF_all_scenarios_30IV_100sim.rds")

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
## Sanity checks
## ============================================================

cat("\nScenario/sample-size check:\n")
print(sim_res_all %>% count(scen, n, name = "n_rows"))

cat("\nA/B strength check:\n")
print(
  sim_res_all %>%
    filter(scen %in% c("A", "B")) %>%
    count(scen, n, strength_label, name = "n_rows")
)

cat("\nC/D sample-size check:\n")
print(
  sim_res_all %>%
    filter(scen %in% c("C", "D")) %>%
    count(scen, n, name = "n_rows")
)

sim_group_check <- sim_res_all %>%
  group_by(scen, n, strength_label, method) %>%
  summarise(
    n_sim_group = n_distinct(sim, cell_within_group),
    .groups = "drop"
  )

print(sim_group_check)

write_csv(sim_group_check, "tables/check_sim_group_counts_30IV_noF.csv")

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
      levels = c(
        strength_levels,
        "500"
      )
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
    title = "Figure 1. Full distribution of bias estimates for each MR method from the four simulation scenarios with thirty instruments.",
    subtitle = "Test run: Scenarios A and B vary instrument strengths at n = 1000; Scenarios C and D are evaluated at n = 500.",
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
  "figures/Figure1_bias_distribution_30IV_noF.png",
  fig1,
  width = 12,
  height = 8.5,
  dpi = 300,
  device = "png"
)

ggsave(
  "figures/Figure1_bias_distribution_30IV_noF.pdf",
  fig1,
  width = 12,
  height = 8.5,
  device = cairo_pdf
)

## ============================================================
## Supplementary Table 1: Bias summary
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
  select(scen, strength_label, method, num_converged, bias_text) %>%
  pivot_wider(
    names_from = scen,
    values_from = c(num_converged, bias_text),
    names_glue = "{scen}_{.value}"
  ) %>%
  arrange(strength_label, method) %>%
  rename(
    `Instrument Strengths` = strength_label,
    `MR Method` = method,
    `Scenario A Num Converged` = A_num_converged,
    `Scenario A Median (Q1, Q3) Bias` = A_bias_text,
    `Scenario B Num Converged` = B_num_converged,
    `Scenario B Median (Q1, Q3) Bias` = B_bias_text
  )

table1_cd <- bias_summary %>%
  filter(scen %in% c("C", "D")) %>%
  select(scen, n, method, num_converged, bias_text) %>%
  pivot_wider(
    names_from = scen,
    values_from = c(num_converged, bias_text),
    names_glue = "{scen}_{.value}"
  ) %>%
  arrange(n, method) %>%
  rename(
    `Sample Size` = n,
    `MR Method` = method,
    `Scenario C Num Converged` = C_num_converged,
    `Scenario C Median (Q1, Q3) Bias` = C_bias_text,
    `Scenario D Num Converged` = D_num_converged,
    `Scenario D Median (Q1, Q3) Bias` = D_bias_text
  )

write_csv(table1_ab, "tables/Supplementary_Table1A_Bias_AB_30IV_noF.csv")
write_csv(table1_cd, "tables/Supplementary_Table1B_Bias_CD_30IV_noF.csv")

## ============================================================
## Supplementary Table 2: Coverage + Type I error / Power
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
  select(scen, strength_label, method, num_converged, coverage_text, rejection_text) %>%
  pivot_wider(
    names_from = scen,
    values_from = c(num_converged, coverage_text, rejection_text),
    names_glue = "{scen}_{.value}"
  ) %>%
  arrange(strength_label, method) %>%
  rename(
    `Instrument Strengths` = strength_label,
    `MR Method` = method,
    `Scenario A Num Converged` = A_num_converged,
    `Scenario A Coverage` = A_coverage_text,
    `Scenario A Type I Error` = A_rejection_text,
    `Scenario B Num Converged` = B_num_converged,
    `Scenario B Coverage` = B_coverage_text,
    `Scenario B Power` = B_rejection_text
  )

table2_cd <- cov_summary %>%
  filter(scen %in% c("C", "D")) %>%
  select(scen, n, method, num_converged, coverage_text, rejection_text) %>%
  pivot_wider(
    names_from = scen,
    values_from = c(num_converged, coverage_text, rejection_text),
    names_glue = "{scen}_{.value}"
  ) %>%
  arrange(n, method) %>%
  rename(
    `Sample Size` = n,
    `MR Method` = method,
    `Scenario C Num Converged` = C_num_converged,
    `Scenario C Coverage` = C_coverage_text,
    `Scenario C Type I Error` = C_rejection_text,
    `Scenario D Num Converged` = D_num_converged,
    `Scenario D Coverage` = D_coverage_text,
    `Scenario D Power` = D_rejection_text
  )

write_csv(table2_ab, "tables/Supplementary_Table2A_Coverage_TypeI_Power_AB_30IV_noF.csv")
write_csv(table2_cd, "tables/Supplementary_Table2B_Coverage_TypeI_Power_CD_30IV_noF.csv")

## ============================================================
## Extra check file
## ============================================================

check_counts <- sim_res_all %>%
  count(scen, n, strength_label, method, conv, name = "n_rows")

write_csv(check_counts, "tables/check_counts_30IV_noF.csv")

message("Done. Figures saved as PNG/PDF and tables saved as CSV.")