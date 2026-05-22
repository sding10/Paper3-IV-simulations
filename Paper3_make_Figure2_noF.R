library(dplyr)
library(purrr)
library(ggplot2)

options(bitmapType = "cairo")

dir.create("figures", showWarnings = FALSE)

numiv_files <- list.files(
  pattern = "^results_noF_NumIV_[0-9]+IV_100sim\\.rds$",
  full.names = TRUE
)

## Exclude old 3-IV file if it still exists
numiv_files <- numiv_files[!grepl("NumIV_03IV", numiv_files)]

cat("Number of NumIV files found:", length(numiv_files), "\n")
print(basename(numiv_files))

if (length(numiv_files) != 7) {
  warning("Expected 7 NumIV files: 1, 5, 10, 15, 20, 25, and 30 IVs.")
}

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
    num_iv = factor(
      as.character(num_iv),
      levels = c("1", "5", "10", "15", "20", "25", "30")
    )
  ) %>%
  filter(conv == 1, is.finite(bias))

cat("\nFigure 2 count check:\n")
print(fig2_dat %>% count(num_iv, method))

fig2 <- ggplot(fig2_dat, aes(x = num_iv, y = bias, fill = method)) +
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
  "figures/Figure2_bias_by_numIV_30IV_noF.png",
  fig2,
  width = 9,
  height = 7,
  dpi = 300,
  device = "png"
)

ggsave(
  "figures/Figure2_bias_by_numIV_30IV_noF.pdf",
  fig2,
  width = 9,
  height = 7,
  device = cairo_pdf
)

saveRDS(sim_numiv, "results_noF_NumIV_all_100sim.rds")

message("Figure 2 saved.")