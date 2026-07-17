Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())
gc()

setwd("/Volumes/FlynnBio/XinyeLi/GSE299207_RAW/")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

input_file <- "cellchat_t_plasma_tables/sample_level_t_plasma_interactions.csv"
output_dir <- "cellchat_t_plasma_tables"

sample_level_df <- read_csv(input_file, show_col_types = FALSE)

sample_totals_df <- sample_level_df %>%
  group_by(
    sample, disease, direction, interaction_name, interaction_name_2,
    ligand, receptor, pathway_name, annotation
  ) %>%
  summarise(prob = sum(prob, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample, disease) %>%
  summarise(
    count = n(),
    weight = sum(prob, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(disease = factor(disease, levels = c("Ctrl", "MGUS", "SM", "MM"))) %>%
  arrange(disease, sample)

group_summary_df <- sample_totals_df %>%
  group_by(disease) %>%
  summarise(
    n_samples = n(),
    mean_count = mean(count),
    median_count = median(count),
    sd_count = sd(count),
    mean_weight = mean(weight),
    median_weight = median(weight),
    sd_weight = sd(weight),
    .groups = "drop"
  )

metric_long_df <- sample_totals_df %>%
  pivot_longer(cols = c(count, weight), names_to = "metric", values_to = "value")

kruskal_df <- metric_long_df %>%
  group_by(metric) %>%
  group_modify(~{
    test <- kruskal.test(value ~ disease, data = .x)
    tibble(
      statistic = unname(test$statistic),
      df = unname(test$parameter),
      p_value = test$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(p_adj_bh = p.adjust(p_value, method = "BH"))

pairwise_compare <- function(dat, metric_name) {
  disease_levels <- levels(dat$disease)
  pairs <- combn(disease_levels, 2, simplify = FALSE)

  pair_rows <- map_dfr(pairs, function(pair) {
    sub_dat <- dat %>%
      filter(disease %in% pair) %>%
      droplevels()

    test <- suppressWarnings(
      wilcox.test(
        value ~ disease,
        data = sub_dat,
        exact = FALSE,
        conf.int = TRUE
      )
    )

    tibble(
      metric = metric_name,
      group1 = pair[[1]],
      group2 = pair[[2]],
      n_group1 = sum(sub_dat$disease == pair[[1]]),
      n_group2 = sum(sub_dat$disease == pair[[2]]),
      group1_mean = mean(sub_dat$value[sub_dat$disease == pair[[1]]]),
      group2_mean = mean(sub_dat$value[sub_dat$disease == pair[[2]]]),
      group1_median = median(sub_dat$value[sub_dat$disease == pair[[1]]]),
      group2_median = median(sub_dat$value[sub_dat$disease == pair[[2]]]),
      statistic = unname(test$statistic),
      p_value = test$p.value,
      conf_low = if (!is.null(test$conf.int)) test$conf.int[[1]] else NA_real_,
      conf_high = if (!is.null(test$conf.int)) test$conf.int[[2]] else NA_real_,
      estimate = if (!is.null(test$estimate)) unname(test$estimate) else NA_real_
    )
  })

  pair_rows %>%
    mutate(
      p_adj_bh = p.adjust(p_value, method = "BH"),
      significant_bh = p_adj_bh < 0.05
    ) %>%
    arrange(p_value, p_adj_bh)
}

pairwise_df <- metric_long_df %>%
  group_split(metric) %>%
  map_dfr(~pairwise_compare(.x, unique(.x$metric)))

focus_df <- pairwise_df %>%
  filter((group1 == "SM" & group2 == "MM") | (group1 == "Ctrl" & group2 == "MM") | (group1 == "MGUS" & group2 == "MM")) %>%
  arrange(metric, group1, group2)

write_csv(sample_totals_df, file.path(output_dir, "average_per_sample_totals.csv"))
write_csv(group_summary_df, file.path(output_dir, "average_per_sample_group_summary.csv"))
write_csv(kruskal_df, file.path(output_dir, "average_per_sample_kruskal_wallis.csv"))
write_csv(pairwise_df, file.path(output_dir, "average_per_sample_pairwise_wilcoxon.csv"))
write_csv(focus_df, file.path(output_dir, "average_per_sample_pairwise_wilcoxon_focus_MM.csv"))

message("Done. Statistical test tables written to: ", output_dir)
