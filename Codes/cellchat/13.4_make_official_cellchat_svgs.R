Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())
gc()

setwd("/Volumes/FlynnBio/XinyeLi/GSE299207_RAW/")

font_cache_dir <- file.path(getwd(), ".fontconfig-cache")
dir.create(font_cache_dir, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = font_cache_dir)

suppressPackageStartupMessages({
  library(CellChat)
  library(ComplexHeatmap)
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(Matrix)
  library(qs)
})

output_root <- "cellchat_svg_figures_official"
compare_dir <- file.path(output_root, "01_compare_coarse")
subtype_dir <- file.path(output_root, "02_compare_subtype")
spatial_dir <- file.path(output_root, "03_spatial")
dir.create(output_root, showWarnings = FALSE)
dir.create(compare_dir, showWarnings = FALSE)
dir.create(subtype_dir, showWarnings = FALSE)
dir.create(spatial_dir, showWarnings = FALSE)

sample_level_file <- "cellchat_t_plasma_tables/sample_level_t_plasma_interactions.csv"
group_lr_file <- "cellchat_t_plasma_tables/group_level_t_plasma_lr_summary.csv"
group_pathway_file <- "cellchat_t_plasma_tables/group_level_t_plasma_pathway_summary.csv"

sample_level_df <- read_csv(sample_level_file, show_col_types = FALSE)
group_lr_df <- read_csv(group_lr_file, show_col_types = FALSE)
group_pathway_df <- read_csv(group_pathway_file, show_col_types = FALSE)

sample_coarse_summary_df <- sample_level_df %>%
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
  mutate(disease = factor(disease, levels = c("Ctrl", "MGUS", "SM", "MM")))

template_obj <- qread("SpatialCellChat_results/cellchat_SM-1.qs")
DB_use <- template_obj@DB
rm(template_obj)
gc()

coarse_groups <- c("T", "Plasma")
shared_t_subtypes <- c("CD4+ T Naive", "CD8+ T Effector", "T Proliferating")
subtype_groups <- c(shared_t_subtypes, "Plasma")

disease_sets <- list(
  Ctrl = "Ctrl",
  MGUS = "MGUS",
  SM = "SM",
  MM = "MM",
  preMM = c("Ctrl", "MGUS", "SM")
)

plot_log <- list()

log_plot <- function(file, status, note = "") {
  plot_log[[length(plot_log) + 1L]] <<- tibble(
    file = file,
    status = status,
    note = note
  )
}

save_svg_plot <- function(file, width, height, expr) {
  tryCatch({
    svg(file, width = width, height = height)
    on.exit(dev.off(), add = TRUE)
    force(expr)
    log_plot(file, "ok")
  }, error = function(e) {
    if (dev.cur() > 1) {
      try(dev.off(), silent = TRUE)
    }
    log_plot(file, "error", conditionMessage(e))
  })
}

make_dummy_cellchat <- function(groups, object_name) {
  cell_names <- paste0(object_name, "_", make.names(groups))
  meta <- data.frame(
    group = factor(groups, levels = groups),
    samples = factor("sample1"),
    row.names = cell_names
  )
  mat <- Matrix(
    matrix(
      1,
      nrow = 2,
      ncol = length(groups),
      dimnames = list(c("gene1", "gene2"), cell_names)
    ),
    sparse = TRUE
  )
  createCellChat(mat, meta = meta, group.by = "group", datatype = "RNA")
}

build_coarse_object <- function(object_name, disease_vec) {
  lr_df <- sample_level_df %>%
    filter(disease %in% disease_vec) %>%
    group_by(
      sample, direction, interaction_name, interaction_name_2,
      ligand, receptor, pathway_name, annotation, evidence
    ) %>%
    summarise(prob = sum(prob, na.rm = TRUE), .groups = "drop") %>%
    group_by(
      direction, interaction_name, interaction_name_2,
      ligand, receptor, pathway_name, annotation, evidence
    ) %>%
    summarise(
      prob = median(prob, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      source = if_else(direction == "T_to_Plasma", "T", "Plasma"),
      target = if_else(direction == "T_to_Plasma", "Plasma", "T")
    )

  pathway_df <- sample_level_df %>%
    filter(disease %in% disease_vec) %>%
    group_by(sample, direction, pathway_name, annotation) %>%
    summarise(prob = sum(prob, na.rm = TRUE), .groups = "drop") %>%
    group_by(direction, pathway_name, annotation) %>%
    summarise(prob = median(prob, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      source = if_else(direction == "T_to_Plasma", "T", "Plasma"),
      target = if_else(direction == "T_to_Plasma", "Plasma", "T")
    )

  lr_names <- sort(unique(lr_df$interaction_name))
  pathway_names <- sort(unique(pathway_df$pathway_name))

  prob_arr <- array(
    0,
    dim = c(length(coarse_groups), length(coarse_groups), length(lr_names)),
    dimnames = list(coarse_groups, coarse_groups, lr_names)
  )
  pval_arr <- array(
    1,
    dim = c(length(coarse_groups), length(coarse_groups), length(lr_names)),
    dimnames = list(coarse_groups, coarse_groups, lr_names)
  )

  for (i in seq_len(nrow(lr_df))) {
    prob_arr[lr_df$source[i], lr_df$target[i], lr_df$interaction_name[i]] <- lr_df$prob[i]
    pval_arr[lr_df$source[i], lr_df$target[i], lr_df$interaction_name[i]] <- 0
  }

  pathway_arr <- array(
    0,
    dim = c(length(coarse_groups), length(coarse_groups), length(pathway_names)),
    dimnames = list(coarse_groups, coarse_groups, pathway_names)
  )

  for (i in seq_len(nrow(pathway_df))) {
    pathway_arr[pathway_df$source[i], pathway_df$target[i], pathway_df$pathway_name[i]] <- pathway_df$prob[i]
  }

  obj <- make_dummy_cellchat(coarse_groups, object_name)
  obj@DB <- DB_use
  obj@LR$LRsig <- lr_df %>%
    distinct(
      interaction_name, pathway_name, ligand, receptor,
      annotation, interaction_name_2, evidence
    ) %>%
    as.data.frame(stringsAsFactors = FALSE)
  rownames(obj@LR$LRsig) <- obj@LR$LRsig$interaction_name
  obj@net <- list(
    prob = prob_arr,
    pval = pval_arr,
    count = apply(prob_arr > 0, c(1, 2), sum),
    weight = apply(prob_arr, c(1, 2), sum),
    LR.sig = dimnames(prob_arr)[[3]]
  )
  obj@netP <- list(
    pathways = pathway_names,
    prob = pathway_arr
  )
  obj <- netAnalysis_computeCentrality(obj, slot.name = "netP")
  obj
}

build_subtype_object <- function(object_name, disease_vec) {
  lr_df <- sample_level_df %>%
    filter(disease %in% disease_vec, t_subtype %in% shared_t_subtypes) %>%
    group_by(
      sample, direction, t_subtype, interaction_name, interaction_name_2,
      ligand, receptor, pathway_name, annotation, evidence
    ) %>%
    summarise(prob = sum(prob, na.rm = TRUE), .groups = "drop") %>%
    group_by(
      direction, t_subtype, interaction_name, interaction_name_2,
      ligand, receptor, pathway_name, annotation, evidence
    ) %>%
    summarise(prob = median(prob, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      source = if_else(direction == "T_to_Plasma", t_subtype, "Plasma"),
      target = if_else(direction == "T_to_Plasma", "Plasma", t_subtype)
    )

  pathway_df <- sample_level_df %>%
    filter(disease %in% disease_vec, t_subtype %in% shared_t_subtypes) %>%
    group_by(sample, direction, t_subtype, pathway_name, annotation) %>%
    summarise(prob = sum(prob, na.rm = TRUE), .groups = "drop") %>%
    group_by(direction, t_subtype, pathway_name, annotation) %>%
    summarise(prob = median(prob, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      source = if_else(direction == "T_to_Plasma", t_subtype, "Plasma"),
      target = if_else(direction == "T_to_Plasma", "Plasma", t_subtype)
    )

  lr_names <- sort(unique(lr_df$interaction_name))
  pathway_names <- sort(unique(pathway_df$pathway_name))

  prob_arr <- array(
    0,
    dim = c(length(subtype_groups), length(subtype_groups), length(lr_names)),
    dimnames = list(subtype_groups, subtype_groups, lr_names)
  )
  pval_arr <- array(
    1,
    dim = c(length(subtype_groups), length(subtype_groups), length(lr_names)),
    dimnames = list(subtype_groups, subtype_groups, lr_names)
  )

  for (i in seq_len(nrow(lr_df))) {
    prob_arr[lr_df$source[i], lr_df$target[i], lr_df$interaction_name[i]] <- lr_df$prob[i]
    pval_arr[lr_df$source[i], lr_df$target[i], lr_df$interaction_name[i]] <- 0
  }

  pathway_arr <- array(
    0,
    dim = c(length(subtype_groups), length(subtype_groups), length(pathway_names)),
    dimnames = list(subtype_groups, subtype_groups, pathway_names)
  )

  for (i in seq_len(nrow(pathway_df))) {
    pathway_arr[pathway_df$source[i], pathway_df$target[i], pathway_df$pathway_name[i]] <- pathway_df$prob[i]
  }

  obj <- make_dummy_cellchat(subtype_groups, object_name)
  obj@DB <- DB_use
  obj@LR$LRsig <- lr_df %>%
    distinct(
      interaction_name, pathway_name, ligand, receptor,
      annotation, interaction_name_2, evidence
    ) %>%
    as.data.frame(stringsAsFactors = FALSE)
  rownames(obj@LR$LRsig) <- obj@LR$LRsig$interaction_name
  obj@net <- list(
    prob = prob_arr,
    pval = pval_arr,
    count = apply(prob_arr > 0, c(1, 2), sum),
    weight = apply(prob_arr, c(1, 2), sum),
    LR.sig = dimnames(prob_arr)[[3]]
  )
  obj@netP <- list(
    pathways = pathway_names,
    prob = pathway_arr
  )
  obj <- netAnalysis_computeCentrality(obj, slot.name = "netP")
  obj
}

coarse_objects <- list(
  Ctrl = build_coarse_object("Ctrl", disease_sets$Ctrl),
  MGUS = build_coarse_object("MGUS", disease_sets$MGUS),
  SM = build_coarse_object("SM", disease_sets$SM),
  MM = build_coarse_object("MM", disease_sets$MM),
  preMM = build_coarse_object("preMM", disease_sets$preMM)
)

subtype_objects <- list(
  SM = build_subtype_object("SM_subtype", disease_sets$SM),
  MM = build_subtype_object("MM_subtype", disease_sets$MM)
)

merged_4group <- mergeCellChat(
  object.list = coarse_objects[c("Ctrl", "MGUS", "SM", "MM")],
  add.names = c("Ctrl", "MGUS", "SM", "MM")
)

merged_preMM_MM <- mergeCellChat(
  object.list = coarse_objects[c("preMM", "MM")],
  add.names = c("preMM", "MM")
)

merged_SM_MM_coarse <- mergeCellChat(
  object.list = coarse_objects[c("SM", "MM")],
  add.names = c("SM", "MM")
)

merged_SM_MM_subtype <- mergeCellChat(
  object.list = subtype_objects[c("SM", "MM")],
  add.names = c("SM", "MM")
)

## Compare plots
save_svg_plot(
  file.path(compare_dir, "average_per_sample_interactions_4groups_count.svg"),
  7, 5,
  print(
    ggplot2::ggplot(
      sample_coarse_summary_df,
      ggplot2::aes(x = disease, y = count, fill = disease)
    ) +
      ggplot2::stat_summary(fun = mean, geom = "col", width = 0.65) +
      ggplot2::stat_summary(
        fun = mean,
        geom = "text",
        ggplot2::aes(label = round(after_stat(y), 1)),
        vjust = -0.4,
        size = 3
      ) +
      ggplot2::geom_point(
        position = ggplot2::position_jitter(width = 0.12, height = 0),
        size = 1.8,
        alpha = 0.8,
        color = "black"
      ) +
      ggplot2::scale_fill_manual(values = c(
        Ctrl = "#1b9e77",
        MGUS = "#d95f02",
        SM = "#7570b3",
        MM = "#e7298a"
      )) +
      ggplot2::labs(
        x = NULL,
        y = "Average inferred interactions per sample",
        title = "Average T-Plasma interaction count per sample"
      ) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
  )
)
save_svg_plot(
  file.path(compare_dir, "average_per_sample_interactions_4groups_weight.svg"),
  7, 5,
  print(
    ggplot2::ggplot(
      sample_coarse_summary_df,
      ggplot2::aes(x = disease, y = weight, fill = disease)
    ) +
      ggplot2::stat_summary(fun = mean, geom = "col", width = 0.65) +
      ggplot2::stat_summary(
        fun = mean,
        geom = "text",
        ggplot2::aes(label = signif(after_stat(y), 3)),
        vjust = -0.4,
        size = 3
      ) +
      ggplot2::geom_point(
        position = ggplot2::position_jitter(width = 0.12, height = 0),
        size = 1.8,
        alpha = 0.8,
        color = "black"
      ) +
      ggplot2::scale_fill_manual(values = c(
        Ctrl = "#1b9e77",
        MGUS = "#d95f02",
        SM = "#7570b3",
        MM = "#e7298a"
      )) +
      ggplot2::labs(
        x = NULL,
        y = "Average interaction weight per sample",
        title = "Average T-Plasma interaction weight per sample"
      ) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
  )
)

save_svg_plot(
  file.path(compare_dir, "compareInteractions_4groups_count.svg"),
  7, 5,
  print(compareInteractions(merged_4group, measure = "count", title.name = "T-Plasma interaction count"))
)
save_svg_plot(
  file.path(compare_dir, "compareInteractions_4groups_weight.svg"),
  7, 5,
  print(compareInteractions(merged_4group, measure = "weight", title.name = "T-Plasma interaction weight"))
)
save_svg_plot(
  file.path(compare_dir, "compareInteractions_preMM_vs_MM_count.svg"),
  6, 5,
  print(compareInteractions(merged_preMM_MM, measure = "count", title.name = "preMM vs MM interaction count"))
)
save_svg_plot(
  file.path(compare_dir, "compareInteractions_preMM_vs_MM_weight.svg"),
  6, 5,
  print(compareInteractions(merged_preMM_MM, measure = "weight", title.name = "preMM vs MM interaction weight"))
)

save_svg_plot(
  file.path(compare_dir, "netVisual_diffInteraction_SM_vs_MM_count.svg"),
  6, 6,
  netVisual_diffInteraction(merged_SM_MM_coarse, comparison = c(1, 2), measure = "count", weight.scale = TRUE)
)
save_svg_plot(
  file.path(compare_dir, "netVisual_diffInteraction_SM_vs_MM_weight.svg"),
  6, 6,
  netVisual_diffInteraction(merged_SM_MM_coarse, comparison = c(1, 2), measure = "weight", weight.scale = TRUE)
)

save_svg_plot(
  file.path(compare_dir, "netVisual_heatmap_SM_vs_MM_count.svg"),
  6, 6,
  draw(netVisual_heatmap(merged_SM_MM_coarse, comparison = c(1, 2), measure = "count", slot.name = "net"))
)
save_svg_plot(
  file.path(compare_dir, "netVisual_heatmap_SM_vs_MM_weight.svg"),
  6, 6,
  draw(netVisual_heatmap(merged_SM_MM_coarse, comparison = c(1, 2), measure = "weight", slot.name = "net"))
)

save_svg_plot(
  file.path(compare_dir, "rankNet_SM_vs_MM_T_to_Plasma.svg"),
  12, 16,
  print(rankNet(
    merged_SM_MM_coarse,
    mode = "comparison",
    comparison = c(1, 2),
    sources.use = "T",
    targets.use = "Plasma",
    title = "SM vs MM: T to Plasma",
    font.size = 6,
    bar.w = 0.6
  ))
)
save_svg_plot(
  file.path(compare_dir, "rankNet_SM_vs_MM_Plasma_to_T.svg"),
  12, 16,
  print(rankNet(
    merged_SM_MM_coarse,
    mode = "comparison",
    comparison = c(1, 2),
    sources.use = "Plasma",
    targets.use = "T",
    title = "SM vs MM: Plasma to T",
    font.size = 6,
    bar.w = 0.6
  ))
)

## Subtype plots
save_svg_plot(
  file.path(subtype_dir, "netVisual_bubble_SM_vs_MM_Tsubtype_to_Plasma_MM_high.svg"),
  12, 7,
  print(netVisual_bubble(
    merged_SM_MM_subtype,
    sources.use = shared_t_subtypes,
    targets.use = "Plasma",
    comparison = c(1, 2),
    max.dataset = 2,
    angle.x = 45,
    title.name = "SM vs MM: T subtype to Plasma, MM enriched"
  ))
)
save_svg_plot(
  file.path(subtype_dir, "netVisual_bubble_SM_vs_MM_Tsubtype_to_Plasma_SM_high.svg"),
  12, 7,
  print(netVisual_bubble(
    merged_SM_MM_subtype,
    sources.use = shared_t_subtypes,
    targets.use = "Plasma",
    comparison = c(1, 2),
    min.dataset = 1,
    angle.x = 45,
    title.name = "SM vs MM: T subtype to Plasma, SM enriched"
  ))
)
save_svg_plot(
  file.path(subtype_dir, "netVisual_bubble_SM_vs_MM_Plasma_to_Tsubtype_MM_high.svg"),
  12, 7,
  print(netVisual_bubble(
    merged_SM_MM_subtype,
    sources.use = "Plasma",
    targets.use = shared_t_subtypes,
    comparison = c(1, 2),
    max.dataset = 2,
    angle.x = 45,
    title.name = "SM vs MM: Plasma to T subtype, MM enriched"
  ))
)
save_svg_plot(
  file.path(subtype_dir, "netVisual_bubble_SM_vs_MM_Plasma_to_Tsubtype_SM_high.svg"),
  12, 7,
  print(netVisual_bubble(
    merged_SM_MM_subtype,
    sources.use = "Plasma",
    targets.use = shared_t_subtypes,
    comparison = c(1, 2),
    min.dataset = 1,
    angle.x = 45,
    title.name = "SM vs MM: Plasma to T subtype, SM enriched"
  ))
)

save_svg_plot(
  file.path(subtype_dir, "netAnalysis_signalingRole_heatmap_SM.svg"),
  12, 8,
  draw(netAnalysis_signalingRole_heatmap(subtype_objects$SM, pattern = "all", width = 10, height = 8))
)
save_svg_plot(
  file.path(subtype_dir, "netAnalysis_signalingRole_heatmap_MM.svg"),
  12, 8,
  draw(netAnalysis_signalingRole_heatmap(subtype_objects$MM, pattern = "all", width = 10, height = 8))
)
save_svg_plot(
  file.path(subtype_dir, "netAnalysis_signalingRole_scatter_SM.svg"),
  8, 6,
  print(netAnalysis_signalingRole_scatter(subtype_objects$SM, title = "SM signaling roles"))
)
save_svg_plot(
  file.path(subtype_dir, "netAnalysis_signalingRole_scatter_MM.svg"),
  8, 6,
  print(netAnalysis_signalingRole_scatter(subtype_objects$MM, title = "MM signaling roles"))
)
save_svg_plot(
  file.path(subtype_dir, "netAnalysis_diff_signalingRole_scatter_SM_vs_MM.svg"),
  8, 6,
  print(netAnalysis_diff_signalingRole_scatter(merged_SM_MM_subtype, comparison = c(1, 2)))
)

## Spatial plots
rep_samples <- sample_level_df %>%
  count(sample, disease, wt = prob, name = "total_prob") %>%
  filter(disease %in% c("SM", "MM")) %>%
  arrange(disease, desc(total_prob)) %>%
  group_by(disease) %>%
  slice_head(n = 1) %>%
  ungroup()

write_csv(rep_samples, file.path(spatial_dir, "representative_samples.csv"))

sample_objects <- list()
for (i in seq_len(nrow(rep_samples))) {
  disease_name <- rep_samples$disease[i]
  sample_name <- rep_samples$sample[i]
  sample_objects[[disease_name]] <- qread(file.path("SpatialCellChat_results", paste0("cellchat_", sample_name, ".qs")))
}

pathway_candidates <- list(
  SM = c("LAMININ", "CXCL", "GAS"),
  MM = c("LAMININ", "CXCL", "GAS", "TNF")
)

pair_candidates <- list(
  SM = c("LAMC1_CD44", "CXCL12_CXCR4"),
  MM = c("LAMC1_CD44", "CXCL9_CXCR3", "TNF_TNFRSF1A")
)

for (disease_name in names(sample_objects)) {
  obj <- sample_objects[[disease_name]]
  sample_name <- rep_samples$sample[rep_samples$disease == disease_name][1]
  available_pathways <- obj@netP$pathways
  available_pairs <- obj@LR$LRsig$interaction_name

  for (pathway_name in intersect(pathway_candidates[[disease_name]], available_pathways)) {
    save_svg_plot(
      file.path(spatial_dir, paste0(sample_name, "_netVisual_aggregate_", pathway_name, "_spatial.svg")),
      8, 7,
      netVisual_aggregate(obj, signaling = pathway_name, layout = "spatial")
    )
  }

  for (pair_name in intersect(pair_candidates[[disease_name]], available_pairs)) {
    save_svg_plot(
      file.path(spatial_dir, paste0(sample_name, "_spatialFeaturePlot_", pair_name, ".svg")),
      8, 6,
      print(spatialFeaturePlot(obj, pairLR.use = pair_name, do.binary = FALSE))
    )
  }
}

plot_log_df <- bind_rows(plot_log)
write_csv(plot_log_df, file.path(output_root, "plot_status.csv"))

message("Done. SVGs written to: ", output_root)
