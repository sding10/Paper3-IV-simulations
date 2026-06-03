library(dplyr)
library(purrr)
library(ggplot2)
library(readr)

options(bitmapType = "cairo")

dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)
dir.create("results_ivnum_combined", showWarnings = FALSE)

ivnum_files <- list.files(
  path = "results_ivnum",
  pattern = "^Paper3_NumIV_[0-9]+IV_[0-9]+sim\\.rds$",
  full.names = TRUE
)

cat("Number of NumIV files found:", length(ivnum_files), "\n")
print(basename(ivnum_files))

if (length(ivnum_files) != 7) {
  warning("Expected 7 NumIV files: 1, 5, 10, 15, 20, 25, and 30 IVs.")
}

sim_ivnum <- map_dfr(ivnum_files, readRDS)

method_levels <- c("2SPS", "2SRI", "GMM", "IV-MVB")
ivnum_levels <- c(1, 5, 10, 15, 20, 25, 30)

method_cols <- c(
  "2SPS" = "grey90",
  "2SRI" = "grey75",
  "GMM" = "grey50",
  "IV-MVB" = "grey15"
)

if (!("ci_alpha" %in% names(sim_ivnum))) {
  sim_ivnum$ci_alpha <- 0.05
}

sim_ivnum <- sim_ivnum %>%
  mutate(
    method = factor(as.character(method), levels = method_levels),
    n_iv = as.integer(n_iv),
    n_iv_label = factor(as.character(n_iv), levels = as.character(ivnum_levels)),
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

cat("\nFigure 2 count check:\n")
print(sim_ivnum %>% count(n_iv, method, conv, name = "n_rows"))

fig2_dat <- sim_ivnum %>%
  filter(conv == 1, is.finite(bias))

fig2 <- ggplot(fig2_dat, aes(x = n_iv_label, y = bias, fill = method)) +
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
  scale_fill_manual(values = method_cols, drop = FALSE) +
  coord_cartesian(ylim = c(-10, 10)) +
  labs(
    title = "Figure 2. Full distribution of bias estimates for each MR method with varying number of instruments.",
    x = "Number of Instruments",
    y = "Estimated Bias in Causal Parameter",
    fill = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16, family = "serif"),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  "figures/Figure2_bias_by_numIV_revised.png",
  fig2,
  width = 9,
  height = 7,
  dpi = 300,
  device = "png"
)

ggsave(
  "figures/Figure2_bias_by_numIV_revised.pdf",
  fig2,
  width = 9,
  height = 7,
  device = cairo_pdf
)

ivnum_summary <- sim_ivnum %>%
  filter(conv == 1, is.finite(bias)) %>%
  group_by(n_iv, method) %>%
  summarise(
    n_converged = n(),
    median_bias = median(bias, na.rm = TRUE),
    q1_bias = quantile(bias, 0.25, na.rm = TRUE),
    q3_bias = quantile(bias, 0.75, na.rm = TRUE),
    coverage = mean(cover, na.rm = TRUE),
    power = mean(reject, na.rm = TRUE),
    median_F_mcf = median(F_mcf, na.rm = TRUE),
    q1_F_mcf = quantile(F_mcf, 0.25, na.rm = TRUE),
    q3_F_mcf = quantile(F_mcf, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    bias_text = sprintf("%.2f (%.2f, %.2f)", median_bias, q1_bias, q3_bias),
    coverage_text = sprintf("%.1f%%", 100 * coverage),
    power_text = sprintf("%.1f%%", 100 * power),
    F_mcf_text = sprintf("%.2f (%.2f, %.2f)", median_F_mcf, q1_F_mcf, q3_F_mcf)
  ) %>%
  arrange(n_iv, method)

write_csv(sim_ivnum, "results_ivnum_combined/Paper3_NumIV_raw_revised.csv")
write_csv(ivnum_summary, "tables/Paper3_NumIV_summary_revised.csv")

saveRDS(sim_ivnum, "results_ivnum_combined/Paper3_NumIV_raw_revised.rds")
saveRDS(ivnum_summary, "results_ivnum_combined/Paper3_NumIV_summary_revised.rds")

message("Figure 2 and NumIV summary saved.")