library(dplyr)
library(purrr)
library(ggplot2)
library(readr)

options(bitmapType = "cairo")

dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

numiv_files <- list.files(
  pattern = "^results_withF_NumIV_[0-9]+IV_100sim\\.rds$",
  full.names = TRUE
)

numiv_files <- numiv_files[!grepl("NumIV_03IV", numiv_files)]

cat("Number of NumIV files found:", length(numiv_files), "\n")
print(basename(numiv_files))

sim_numiv <- map_dfr(numiv_files, readRDS)

method_levels <- c("2SPS", "2SRI", "GMM", "IV-MVB")

method_cols <- c(
  "2SPS" = "grey90",
  "2SRI" = "grey75",
  "GMM" = "grey50",
  "IV-MVB" = "grey15"
)

fig2_dat <- sim_numiv %>%
  mutate(
    bias = est - b1,
    method = factor(as.character(method), levels = method_levels),
    
    F_group = case_when(
      first_stage_F < 2 ~ "Weak IVs: F < 2",
      first_stage_F >= 2 & first_stage_F < 5 ~ "Moderate IVs: 2 ≤ F < 5",
      first_stage_F >= 5 ~ "Strong IVs: F ≥ 5",
      TRUE ~ NA_character_
    ),
    
    F_group = factor(
      F_group,
      levels = c(
        "Weak IVs: F < 2",
        "Moderate IVs: 2 ≤ F < 5",
        "Strong IVs: F ≥ 5"
      )
    )
  ) %>%
  filter(conv == 1, is.finite(bias), !is.na(F_group))

## Check how many simulations fall into each F group
F_count_check <- fig2_dat %>%
  distinct(num_iv, sim, first_stage_F, F_group) %>%
  count(F_group, num_iv)

print(F_count_check)

write_csv(
  F_count_check,
  "tables/Figure2_Fgroup_NumIV_Count_Check_withF.csv"
)

F_summary <- fig2_dat %>%
  distinct(num_iv, sim, first_stage_F, F_group) %>%
  group_by(F_group) %>%
  summarise(
    n_replicates = n(),
    median_F = median(first_stage_F, na.rm = TRUE),
    q1_F = quantile(first_stage_F, 0.25, na.rm = TRUE),
    q3_F = quantile(first_stage_F, 0.75, na.rm = TRUE),
    min_F = min(first_stage_F, na.rm = TRUE),
    max_F = max(first_stage_F, na.rm = TRUE),
    .groups = "drop"
  )

print(F_summary)

write_csv(
  F_summary,
  "tables/Figure2_Fgroup_Summary_withF.csv"
)

fig2 <- ggplot(fig2_dat, aes(x = F_group, y = bias, fill = method)) +
  geom_hline(
    yintercept = 0,
    color = "darkgreen",
    linewidth = 0.30
  ) +
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
  scale_fill_manual(values = method_cols, drop = FALSE) +
  coord_cartesian(ylim = c(-10, 10)) +
  labs(
    title = "Figure 2. Full distribution of bias estimates by empirical first-stage F-statistic group.",
    subtitle = "Instrument strength groups are defined using the realized first-stage F-statistic.",
    x = "Empirical First-Stage F-Statistic Group",
    y = "Estimated Bias in Causal Parameter",
    fill = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16, family = "serif"),
    plot.subtitle = element_text(size = 11, family = "serif"),
    axis.title = element_text(face = "bold", size = 13),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  "figures/Figure2_bias_by_Fgroup_30IV_withF.png",
  fig2,
  width = 10,
  height = 7,
  dpi = 300,
  device = "png"
)

ggsave(
  "figures/Figure2_bias_by_Fgroup_30IV_withF.pdf",
  fig2,
  width = 10,
  height = 7,
  device = cairo_pdf
)

saveRDS(sim_numiv, "results_withF_NumIV_all_100sim.rds")

message("Figure 2 based fully on F-statistic groups saved.")