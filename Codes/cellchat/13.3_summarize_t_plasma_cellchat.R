Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())
gc()

setwd("/Volumes/FlynnBio/XinyeLi/GSE299207_RAW/")

suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(readr)
  library(stringr)
  library(tibble)
  library(qs)
  library(SpatialCellChat)
})

input_dir <- "SpatialCellChat_results"
output_dir <- "cellchat_t_plasma_tables"
dir.create(output_dir, showWarnings = FALSE)

if (!dir.exists(input_dir)) {
  stop("Input directory not found: ", input_dir)
}

sample_to_disease <- function(sample_name) {
  dplyr::case_when(
    str_starts(sample_name, "Ctrl-") ~ "Ctrl",
    str_starts(sample_name, "MGUS-") ~ "MGUS",
    str_starts(sample_name, "SM-") ~ "SM",
    str_starts(sample_name, "MM-") ~ "MM",
    TRUE ~ NA_character_
  )
}

is_t_label <- function(x) {
  str_detect(x, "(^CD4\\+ T|^CD8\\+ T|^T\\b)")
}

is_plasma_label <- function(x) {
  x %in% c("Plasma", "Malignant_Plasma")
}

collect_evidence <- function(x) {
  vals <- sort(unique(stats::na.omit(x)))
  if (length(vals) == 0) {
    return(NA_character_)
  }
  paste(vals, collapse = " | ")
}

safe_min_numeric <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  min(x, na.rm = TRUE)
}

extract_sample_name <- function(path) {
  basename(path) %>%
    str_remove("^cellchat_") %>%
    str_remove("[.]qs$")
}

extract_lr_metadata <- function(cellchat_obj) {
  lr_df <- cellchat_obj@LR$LRsig %>%
    as.data.frame(stringsAsFactors = FALSE)

  if (!"interaction_name" %in% colnames(lr_df)) {
    lr_df$interaction_name <- rownames(lr_df)
  }

  as_tibble(lr_df) %>%
    transmute(
      interaction_name = as.character(interaction_name),
      interaction_name_2 = as.character(interaction_name_2),
      pathway_name = as.character(pathway_name),
      ligand = as.character(ligand),
      receptor = as.character(receptor),
      annotation = as.character(annotation),
      evidence = as.character(evidence)
    )
}

finalize_t_plasma_table <- function(comm_df, sample_name, disease, mode_label) {
  if (nrow(comm_df) == 0) {
    return(tibble())
  }

  comm_df %>%
    mutate(
      sample = sample_name,
      disease = disease,
      sample_extraction_mode = mode_label,
      source_raw = as.character(source_raw),
      target_raw = as.character(target_raw),
      source_is_t = is_t_label(source_raw),
      target_is_t = is_t_label(target_raw),
      source_is_plasma = is_plasma_label(source_raw),
      target_is_plasma = is_plasma_label(target_raw)
    ) %>%
    filter(
      (source_is_t & target_is_plasma) |
        (source_is_plasma & target_is_t)
    ) %>%
    mutate(
      source_std = if_else(source_is_t, "T", "Plasma"),
      target_std = if_else(target_is_t, "T", "Plasma"),
      direction = if_else(source_is_t, "T_to_Plasma", "Plasma_to_T"),
      t_subtype = if_else(source_is_t, source_raw, target_raw),
      plasma_type = if_else(source_is_plasma, source_raw, target_raw)
    ) %>%
    transmute(
      sample,
      disease,
      sample_extraction_mode,
      source_raw,
      target_raw,
      source_std,
      target_std,
      direction,
      t_subtype,
      plasma_type,
      ligand,
      receptor,
      interaction_name,
      interaction_name_2,
      pathway_name,
      prob = as.numeric(prob),
      pval = as.numeric(pval),
      annotation,
      evidence
    )
}

extract_from_group_prob <- function(cellchat_obj, sample_name, disease, thresh = 0.05) {
  lr_meta <- extract_lr_metadata(cellchat_obj)
  prob_arr <- cellchat_obj@net$prob
  pval_arr <- cellchat_obj@net$pval

  if (is.null(prob_arr)) {
    return(tibble())
  }

  dims <- dim(prob_arr)
  dimn <- dimnames(prob_arr)

  comm_df <- expand.grid(
    source_raw = dimn[[1]],
    target_raw = dimn[[2]],
    interaction_name = dimn[[3]],
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    mutate(
      prob = as.numeric(prob_arr),
      pval = if (!is.null(pval_arr)) as.numeric(pval_arr) else NA_real_
    ) %>%
    filter(prob > 0) %>%
    filter(if (!is.null(pval_arr)) pval < thresh else TRUE) %>%
    left_join(lr_meta, by = "interaction_name")

  finalize_t_plasma_table(comm_df, sample_name, disease, "group_prob")
}

extract_from_cell_prob <- function(cellchat_obj, sample_name, disease) {
  prob_list <- cellchat_obj@net$tmp$prob.cell
  lr_meta <- extract_lr_metadata(cellchat_obj)

  if (is.null(prob_list) || length(prob_list) == 0) {
    return(tibble())
  }

  group_lookup <- setNames(
    as.character(cellchat_obj@meta$cellchat_group),
    rownames(cellchat_obj@meta)
  )
  cell_names <- rownames(cellchat_obj@meta)
  group_sizes <- table(group_lookup)

  comm_list <- vector("list", length(prob_list))
  keep_idx <- 0L

  for (lr_name in names(prob_list)) {
    prob_mat <- prob_list[[lr_name]]

    if (!inherits(prob_mat, "dgCMatrix") || length(prob_mat@x) == 0) {
      next
    }

    triplet_df <- Matrix::summary(prob_mat) %>%
      as_tibble() %>%
      transmute(
        source_raw = group_lookup[cell_names[i]],
        target_raw = group_lookup[cell_names[j]],
        prob_value = x
      ) %>%
      filter(
        (is_t_label(source_raw) & is_plasma_label(target_raw)) |
          (is_plasma_label(source_raw) & is_t_label(target_raw))
      )

    if (nrow(triplet_df) == 0) {
      next
    }

    keep_idx <- keep_idx + 1L
    comm_list[[keep_idx]] <- triplet_df %>%
      group_by(source_raw, target_raw) %>%
      summarise(
        prob = sum(prob_value, na.rm = TRUE) /
          (as.numeric(group_sizes[first(source_raw)]) *
             as.numeric(group_sizes[first(target_raw)])),
        pval = NA_real_,
        .groups = "drop"
      ) %>%
      mutate(interaction_name = lr_name) %>%
      left_join(lr_meta, by = "interaction_name")
  }

  if (keep_idx == 0L) {
    return(tibble())
  }

  comm_df <- bind_rows(comm_list[seq_len(keep_idx)])
  finalize_t_plasma_table(comm_df, sample_name, disease, "cell_pair_mean")
}

extract_t_plasma_for_sample <- function(file_path, thresh = 0.05) {
  sample_name <- extract_sample_name(file_path)
  disease <- sample_to_disease(sample_name)

  if (is.na(disease)) {
    stop("Cannot map sample to disease group: ", sample_name)
  }

  message("Processing ", sample_name)
  cellchat_obj <- qread(file_path)

  has_group_prob <- "prob" %in% names(cellchat_obj@net)
  has_cell_prob <- "prob.cell" %in% names(cellchat_obj@net) ||
    (!is.null(cellchat_obj@net$tmp$prob.cell))

  comm_df <- if (has_group_prob) {
    extract_from_group_prob(cellchat_obj, sample_name, disease, thresh = thresh)
  } else if (has_cell_prob) {
    extract_from_cell_prob(cellchat_obj, sample_name, disease)
  } else {
    stop("Unsupported CellChat net structure for sample: ", sample_name)
  }

  rm(cellchat_obj)
  gc()
  comm_df
}

sample_files <- list.files(
  input_dir,
  pattern = "^cellchat_.*[.]qs$",
  full.names = TRUE
) %>%
  sort()

if (length(sample_files) == 0) {
  stop("No cellchat result files found in ", input_dir)
}

sample_manifest <- tibble(
  file_path = sample_files,
  sample = vapply(sample_files, extract_sample_name, character(1)),
  disease = vapply(
    vapply(sample_files, extract_sample_name, character(1)),
    sample_to_disease,
    character(1)
  )
) %>%
  arrange(factor(disease, levels = c("Ctrl", "MGUS", "SM", "MM")), sample)

if (any(is.na(sample_manifest$disease))) {
  stop("Unknown disease group detected in sample manifest.")
}

sample_counts <- sample_manifest %>%
  count(disease, name = "n_samples_group")

sample_level_df <- bind_rows(lapply(sample_manifest$file_path, extract_t_plasma_for_sample))

if (nrow(sample_level_df) == 0) {
  stop("No significant T-Plasma interactions were found.")
}

sample_level_df <- sample_level_df %>%
  arrange(
    factor(disease, levels = c("Ctrl", "MGUS", "SM", "MM")),
    sample,
    direction,
    desc(prob),
    interaction_name
  )

sample_overall_lr <- sample_level_df %>%
  group_by(
    sample, disease, direction, interaction_name, interaction_name_2,
    ligand, receptor, pathway_name, annotation
  ) %>%
  summarise(
    prob = sum(prob, na.rm = TRUE),
    pval = safe_min_numeric(pval),
    n_component_rows = n(),
    t_subtypes = paste(sort(unique(t_subtype)), collapse = " | "),
    plasma_types = paste(sort(unique(plasma_type)), collapse = " | "),
    evidence = collect_evidence(evidence),
    .groups = "drop"
  )

group_lr_summary <- sample_overall_lr %>%
  group_by(
    disease, direction, interaction_name, interaction_name_2,
    ligand, receptor, pathway_name, annotation
  ) %>%
  summarise(
    n_samples_detected = n_distinct(sample),
    mean_prob = mean(prob, na.rm = TRUE),
    median_prob = median(prob, na.rm = TRUE),
    sd_prob = sd(prob, na.rm = TRUE),
    min_prob = min(prob, na.rm = TRUE),
    max_prob = max(prob, na.rm = TRUE),
    mean_component_rows = mean(n_component_rows, na.rm = TRUE),
    samples = paste(sort(unique(sample)), collapse = " | "),
    t_subtypes_seen = paste(sort(unique(t_subtypes)), collapse = " || "),
    plasma_types_seen = paste(sort(unique(plasma_types)), collapse = " || "),
    evidence = collect_evidence(evidence),
    .groups = "drop"
  ) %>%
  left_join(sample_counts, by = "disease") %>%
  mutate(prop_samples_detected = n_samples_detected / n_samples_group) %>%
  select(
    disease, n_samples_group, direction, interaction_name, interaction_name_2,
    ligand, receptor, pathway_name, annotation, n_samples_detected,
    prop_samples_detected, mean_prob, median_prob, sd_prob, min_prob, max_prob,
    mean_component_rows, samples, t_subtypes_seen, plasma_types_seen, evidence
  ) %>%
  arrange(
    factor(disease, levels = c("Ctrl", "MGUS", "SM", "MM")),
    direction,
    desc(n_samples_detected),
    desc(median_prob),
    desc(mean_prob),
    interaction_name
  )

sample_overall_pathway <- sample_level_df %>%
  group_by(sample, disease, direction, pathway_name, annotation) %>%
  summarise(
    prob_sum = sum(prob, na.rm = TRUE),
    n_unique_lr_pairs = n_distinct(interaction_name),
    t_subtypes = paste(sort(unique(t_subtype)), collapse = " | "),
    plasma_types = paste(sort(unique(plasma_type)), collapse = " | "),
    evidence = collect_evidence(evidence),
    .groups = "drop"
  )

group_pathway_summary <- sample_overall_pathway %>%
  group_by(disease, direction, pathway_name, annotation) %>%
  summarise(
    n_samples_detected = n_distinct(sample),
    mean_prob_sum = mean(prob_sum, na.rm = TRUE),
    median_prob_sum = median(prob_sum, na.rm = TRUE),
    sd_prob_sum = sd(prob_sum, na.rm = TRUE),
    min_prob_sum = min(prob_sum, na.rm = TRUE),
    max_prob_sum = max(prob_sum, na.rm = TRUE),
    mean_n_unique_lr_pairs = mean(n_unique_lr_pairs, na.rm = TRUE),
    samples = paste(sort(unique(sample)), collapse = " | "),
    t_subtypes_seen = paste(sort(unique(t_subtypes)), collapse = " || "),
    plasma_types_seen = paste(sort(unique(plasma_types)), collapse = " || "),
    evidence = collect_evidence(evidence),
    .groups = "drop"
  ) %>%
  left_join(sample_counts, by = "disease") %>%
  mutate(prop_samples_detected = n_samples_detected / n_samples_group) %>%
  select(
    disease, n_samples_group, direction, pathway_name, annotation,
    n_samples_detected, prop_samples_detected, mean_prob_sum, median_prob_sum,
    sd_prob_sum, min_prob_sum, max_prob_sum, mean_n_unique_lr_pairs, samples,
    t_subtypes_seen, plasma_types_seen, evidence
  ) %>%
  arrange(
    factor(disease, levels = c("Ctrl", "MGUS", "SM", "MM")),
    direction,
    desc(n_samples_detected),
    desc(median_prob_sum),
    desc(mean_prob_sum),
    pathway_name
  )

sample_subtype_lr <- sample_level_df %>%
  group_by(
    sample, disease, t_subtype, plasma_type, direction, interaction_name,
    interaction_name_2, ligand, receptor, pathway_name, annotation
  ) %>%
  summarise(
    prob = sum(prob, na.rm = TRUE),
    pval = safe_min_numeric(pval),
    evidence = collect_evidence(evidence),
    .groups = "drop"
  )

group_subtype_summary <- sample_subtype_lr %>%
  group_by(
    disease, t_subtype, plasma_type, direction, interaction_name,
    interaction_name_2, ligand, receptor, pathway_name, annotation
  ) %>%
  summarise(
    n_samples_detected = n_distinct(sample),
    mean_prob = mean(prob, na.rm = TRUE),
    median_prob = median(prob, na.rm = TRUE),
    sd_prob = sd(prob, na.rm = TRUE),
    min_prob = min(prob, na.rm = TRUE),
    max_prob = max(prob, na.rm = TRUE),
    samples = paste(sort(unique(sample)), collapse = " | "),
    evidence = collect_evidence(evidence),
    .groups = "drop"
  ) %>%
  left_join(sample_counts, by = "disease") %>%
  mutate(prop_samples_detected = n_samples_detected / n_samples_group) %>%
  select(
    disease, n_samples_group, t_subtype, plasma_type, direction,
    interaction_name, interaction_name_2, ligand, receptor, pathway_name,
    annotation, n_samples_detected, prop_samples_detected, mean_prob,
    median_prob, sd_prob, min_prob, max_prob, samples, evidence
  ) %>%
  arrange(
    factor(disease, levels = c("Ctrl", "MGUS", "SM", "MM")),
    t_subtype,
    direction,
    desc(n_samples_detected),
    desc(median_prob),
    desc(mean_prob),
    interaction_name
  )

write_csv(sample_level_df, file.path(output_dir, "sample_level_t_plasma_interactions.csv"))
write_csv(group_lr_summary, file.path(output_dir, "group_level_t_plasma_lr_summary.csv"))
write_csv(group_pathway_summary, file.path(output_dir, "group_level_t_plasma_pathway_summary.csv"))
write_csv(group_subtype_summary, file.path(output_dir, "group_level_t_subtype_plasma_summary.csv"))

message("Done. Tables written to: ", output_dir)
