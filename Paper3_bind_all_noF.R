library(dplyr)
library(purrr)
library(readr)
library(tibble)

all_files <- list.files(
  pattern = "^results_noF_scenario_[ABCD]_n[0-9]+_cell[0-9]+\\.rds$",
  full.names = TRUE
)

files <- all_files

cat("Total scenario result files found:", length(all_files), "\n")
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
  expand.grid(
    scen = c("C", "D"),
    n = c(500)
  ) %>%
    as_tibble() %>%
    mutate(expected = 9)
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
  stop("Expected 108 files, but found ", length(files))
}

sim_res_all <- map_dfr(files, readRDS)

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

saveRDS(sim_res_all, "results_noF_all_scenarios_30IV_100sim.rds")
write_csv(sim_res_all, "results_noF_all_scenarios_30IV_100sim.csv")

cat("\nSaved:\n")
cat("results_noF_all_scenarios_30IV_100sim.rds\n")
cat("results_noF_all_scenarios_30IV_100sim.csv\n")