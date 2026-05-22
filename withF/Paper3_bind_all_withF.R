library(dplyr)
library(purrr)
library(readr)
library(tibble)

all_files <- list.files(
  pattern = "^results_withF_newStrength_scenario_[ABCD]_n[0-9]+_cell[0-9]+\\.rds$",
  full.names = TRUE
)

files <- all_files

cat("Total withF new-strength scenario result files found:", length(all_files), "\n")
cat("Files used:", length(files), "\n\n")

file_check <- tibble(file = basename(files)) %>%
  mutate(
    scen = sub(".*scenario_([ABCD])_n.*", "\\1", file),
    n = as.numeric(sub(".*_n([0-9]+)_cell.*", "\\1", file)),
    cell = as.integer(sub(".*_cell([0-9]+)\\.rds", "\\1", file))
  ) %>%
  count(scen, n, name = "n_files")

cat("File count check:\n")
print(file_check)

expected_check <- bind_rows(
  tibble(scen = "A", n = 1000, expected = 45),
  tibble(scen = "B", n = 1000, expected = 45),
  tibble(scen = "C", n = 500, expected = 9),
  tibble(scen = "D", n = 500, expected = 9)
)

check_compare <- expected_check %>%
  left_join(file_check, by = c("scen", "n")) %>%
  mutate(
    n_files = ifelse(is.na(n_files), 0L, n_files),
    ok = n_files == expected
  )

cat("\nExpected vs observed files:\n")
print(check_compare)

if (!all(check_compare$ok)) {
  stop("Some scenario/sample-size cells have missing or extra files. Check before binding.")
}

if (length(files) != 108) {
  stop("Expected 108 files for test run, but found ", length(files))
}

sim_res_all <- map_dfr(files, readRDS)

if (!"first_stage_F" %in% names(sim_res_all)) {
  stop("first_stage_F column is missing. Check withF methods/simulation files.")
}

sim_res_all <- sim_res_all %>%
  mutate(
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
    ),
    scen = factor(as.character(scen), levels = c("A", "B", "C", "D")),
    method = factor(as.character(method), levels = c("2SPS", "2SRI", "GMM", "IV-MVB")),
    strength_label = factor(
      as.character(strength_label),
      levels = c(
        "Very Weak IVs",
        "Weak IVs",
        "Moderately Weak IVs",
        "Moderate IVs",
        "Moderately Strong IVs"
      )
    )
  )

cat("\nRows by scenario and sample size:\n")
print(sim_res_all %>% count(scen, n, name = "n_rows"))

cat("\nRows by scenario, method, convergence:\n")
print(sim_res_all %>% count(scen, method, conv, name = "n_rows"))

cat("\nA/B strength checks:\n")
print(
  sim_res_all %>%
    filter(scen %in% c("A", "B")) %>%
    count(scen, n, strength_label, name = "n_rows")
)

cat("\nC/D sample-size checks:\n")
print(
  sim_res_all %>%
    filter(scen %in% c("C", "D")) %>%
    count(scen, n, name = "n_rows")
)

cat("\nNumber of IVs check:\n")
print(sim_res_all %>% count(n_iv, name = "n_rows"))

cat("\nSimulation-count check by displayed group:\n")
group_check <- sim_res_all %>%
  group_by(scen, n, strength_label, method) %>%
  summarise(
    n_sim_group = n_distinct(sim, cell_within_group),
    .groups = "drop"
  )

print(group_check)

bad_groups <- group_check %>%
  filter(n_sim_group != 100)

if (nrow(bad_groups) > 0) {
  warning("Some groups do not have exactly 100 simulations. Check below:")
  print(bad_groups)
}

cat("\nFirst-stage F summary by scenario and sample size:\n")
print(
  sim_res_all %>%
    distinct(scen, n, strength_label, cell_id, cell_within_group, sim, first_stage_F) %>%
    group_by(scen, n) %>%
    summarise(
      median_F = median(first_stage_F, na.rm = TRUE),
      q1_F = quantile(first_stage_F, 0.25, na.rm = TRUE),
      q3_F = quantile(first_stage_F, 0.75, na.rm = TRUE),
      min_F = min(first_stage_F, na.rm = TRUE),
      max_F = max(first_stage_F, na.rm = TRUE),
      .groups = "drop"
    )
)

cat("\nFirst-stage F summary by A/B IV-strength group:\n")
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
    .groups = "drop"
  )

print(f_summary_ab)

write_csv(group_check, "tables/check_sim_group_counts_30IV_withF_newStrength.csv")
write_csv(f_summary_ab, "tables/FirstStageF_Summary_AB_30IV_withF_newStrength.csv")

saveRDS(sim_res_all, "results_withF_newStrength_all_scenarios_30IV_100sim.rds")
write_csv(sim_res_all, "results_withF_newStrength_all_scenarios_30IV_100sim.csv")

cat("\nSaved:\n")
cat("results_withF_newStrength_all_scenarios_30IV_100sim.rds\n")
cat("results_withF_newStrength_all_scenarios_30IV_100sim.csv\n")