# ==============================================================================
# 12.3.1  Niche 与临床关联分析（3 分类器版，合并 _neg/_pos → 真实患者级）
#
# 与 12.3 的区别：
#   去掉 sample_id 的 _neg/_pos 后缀，按真实患者聚合软概率均值
#   等价于把同一患者的所有细胞合并后取均值（细胞数加权）
#
# 输出（./niche_classifier_3class/clinical_correlation_merged/）：
#   {disease}_{ct}_soft_scores_per_patient_merged.csv
#   {disease}_{ct}_spearman_all3_merged.csv
#   {disease}_{ct}_prob_T_vs_grade_scatter_merged.pdf
#   {disease}_{ct}_prob_T_vs_grade_boxplot_merged.pdf
#   {disease}_{ct}_soft_composition_bar_merged.pdf
#   {disease}_{ct}_3niche_spearman_panel_merged.pdf
#   {disease}_all_summary_merged.pdf
# ==============================================================================

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls()); gc()

setwd("/Volumes/FlynnBio/XinyeLi/GSE299207_RAW/")

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

PRED_DIR <- "./niche_classifier_3class/predictions"
OUTDIR   <- "./niche_classifier_3class/clinical_correlation_merged"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

niche_colors <- c(
  "T_Dominant"      = "#E64B35",
  "Mixed"           = "#F39B7F",
  "Plasma_Dominant" = "#4DBBD5"
)

# ── 辅助函数 ──────────────────────────────────────────────────────────────────

strip_suffix <- function(x) sub("_(neg|pos)$", "", x)

# 按真实患者聚合软概率（细胞数加权 = 直接取所有细胞均值）
compute_soft_scores_merged <- function(pred_df) {
  pred_df %>%
    mutate(patient_id = strip_suffix(sample_id)) %>%
    group_by(patient_id, clinical_grade) %>%
    summarise(
      mean_prob_T      = mean(prob_T_Dominant,     na.rm = TRUE),
      mean_prob_Mixed  = mean(prob_Mixed,           na.rm = TRUE),
      mean_prob_Plasma = mean(prob_Plasma_Dominant, na.rm = TRUE),
      n_cells          = n(),
      n_subsamples     = n_distinct(sample_id),
      .groups          = "drop"
    )
}

# Spearman 相关
run_spearman <- function(summary_df, y_col, grade_levels) {
  grade_num <- setNames(seq_along(grade_levels), grade_levels)
  x <- grade_num[as.character(summary_df$clinical_grade)]
  y <- summary_df[[y_col]]
  ok <- !is.na(x) & !is.na(y)
  test <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  data.frame(
    niche   = y_col,
    rho     = round(unname(test$estimate), 4),
    p.value = round(test$p.value, 4),
    n       = sum(ok)
  )
}

# ── 主流程 ────────────────────────────────────────────────────────────────────

datasets <- list(
  list(disease     = "SM",
       grade_label = "SM Risk Grade",
       grade_order = c("low", "intermediate", "high")),
  list(disease     = "MM",
       grade_label = "R-ISS",
       grade_order = c("1", "2", "3"))
)

cell_types <- c("T_cell", "Plasma_cell")

for (ds in datasets) {
  disease     <- ds$disease
  grade_label <- ds$grade_label
  grade_order <- ds$grade_order

  cat("\n", strrep("=", 60), "\n")
  cat("Disease:", disease, "\n")
  cat(strrep("=", 60), "\n")

  scatter_T_list <- list()

  for (ct in cell_types) {
    pred_file <- file.path(PRED_DIR, paste0(disease, "_", ct, "_niche_pred.csv"))
    if (!file.exists(pred_file)) {
      cat("  File not found:", pred_file, "\n"); next
    }

    pred <- read.csv(pred_file)
    pred$clinical_grade <- as.character(pred$clinical_grade)
    if (disease == "MM") {
      pred$clinical_grade <- sub("\\.0$", "", pred$clinical_grade)
    }
    pred <- pred %>% filter(!is.na(clinical_grade), clinical_grade != "NA")
    pred$clinical_grade <- factor(pred$clinical_grade, levels = grade_order)

    # ── 合并 _neg/_pos → 真实患者 ─────────────────────────────────────────
    scores <- compute_soft_scores_merged(pred)
    scores$clinical_grade <- factor(scores$clinical_grade, levels = grade_order)

    cat(sprintf("\n  %s | %d subsamples → %d true patients\n",
                ct, n_distinct(pred$sample_id), nrow(scores)))

    out_scores <- file.path(OUTDIR,
                            paste0(disease, "_", ct, "_soft_scores_per_patient_merged.csv"))
    write.csv(scores, out_scores, row.names = FALSE)
    cat("  Saved:", out_scores, "\n")
    print(scores)

    # ── Spearman：3 个 niche 各跑一次 ────────────────────────────────────
    sp_list <- lapply(
      c("mean_prob_T", "mean_prob_Mixed", "mean_prob_Plasma"),
      function(col) run_spearman(scores, col, grade_order)
    )
    sp_all <- do.call(rbind, sp_list)
    sp_all$disease   <- disease
    sp_all$cell_type <- ct
    cat("\n  Spearman results (merged):\n")
    print(sp_all)

    out_sp <- file.path(OUTDIR,
                        paste0(disease, "_", ct, "_spearman_all3_merged.csv"))
    write.csv(sp_all, out_sp, row.names = FALSE)
    cat("  Saved:", out_sp, "\n")

    # ── 散点图（T_Dominant 主图）─────────────────────────────────────────
    grade_num_map <- setNames(seq_along(grade_order), grade_order)
    scores$grade_num <- grade_num_map[as.character(scores$clinical_grade)]

    sp_T     <- sp_all[sp_all$niche == "mean_prob_T", ]
    label_T  <- sprintf("Spearman ρ = %.3f\np = %.4f\nN = %d patients",
                        sp_T$rho, sp_T$p.value, sp_T$n)

    p_scatter <- ggplot(scores, aes(x = grade_num, y = mean_prob_T)) +
      geom_smooth(method = "lm", se = TRUE, color = "grey60", fill = "grey85",
                  linewidth = 0.8) +
      geom_jitter(aes(color = clinical_grade), width = 0.05, size = 3.5, alpha = 0.9) +
      annotate("text", x = min(scores$grade_num) + 0.1,
               y = max(scores$mean_prob_T, na.rm = TRUE),
               label = label_T, hjust = 0, vjust = 1, size = 3.5) +
      scale_x_continuous(breaks = seq_along(grade_order), labels = grade_order) +
      scale_color_brewer(palette = "Set1") +
      labs(
        title = paste0(disease, " — ", ct, " (3-class, merged)"),
        x     = grade_label,
        y     = "Mean prob(T_Dominant) per patient",
        color = grade_label
      ) +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))

    out_scatter <- file.path(OUTDIR,
                             paste0(disease, "_", ct, "_prob_T_vs_grade_scatter_merged.pdf"))
    ggsave(out_scatter, p_scatter, width = 5, height = 4, dpi = 300)
    cat("  Saved:", out_scatter, "\n")
    scatter_T_list[[ct]] <- p_scatter

    # ── Boxplot ───────────────────────────────────────────────────────────
    p_box <- ggboxplot(
      scores,
      x = "clinical_grade", y = "mean_prob_T",
      color = "clinical_grade", fill = "clinical_grade",
      alpha = 0.4, add = "jitter",
      add.params = list(size = 2.5, alpha = 0.9)
    ) +
      labs(
        title = paste0(disease, " — ", ct, " (3-class, merged)"),
        x     = grade_label,
        y     = "Mean prob(T_Dominant) per patient"
      ) +
      theme_classic(base_size = 12) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold", hjust = 0.5)) +
      scale_color_brewer(palette = "Set1") +
      scale_fill_brewer(palette = "Set1")

    out_box <- file.path(OUTDIR,
                         paste0(disease, "_", ct, "_prob_T_vs_grade_boxplot_merged.pdf"))
    ggsave(out_box, p_box, width = 5, height = 4, dpi = 300)
    cat("  Saved:", out_box, "\n")

    # ── 软成分堆积柱状图 ──────────────────────────────────────────────────
    bar_df <- scores %>%
      select(patient_id, clinical_grade,
             T_Dominant      = mean_prob_T,
             Mixed           = mean_prob_Mixed,
             Plasma_Dominant = mean_prob_Plasma) %>%
      pivot_longer(cols = c(T_Dominant, Mixed, Plasma_Dominant),
                   names_to = "niche", values_to = "mean_prob") %>%
      mutate(niche = factor(niche, levels = c("T_Dominant", "Mixed", "Plasma_Dominant")))

    order_df <- bar_df %>%
      filter(niche == "T_Dominant") %>%
      arrange(clinical_grade, desc(mean_prob))
    bar_df$patient_id <- factor(bar_df$patient_id, levels = order_df$patient_id)

    p_bar <- ggplot(bar_df, aes(x = patient_id, y = mean_prob, fill = niche)) +
      geom_col(width = 0.8) +
      scale_fill_manual(values = niche_colors, name = "Niche") +
      facet_grid(~ clinical_grade, scales = "free_x", space = "free_x") +
      labs(
        title = paste0(disease, " — ", ct, " soft niche composition (3-class, merged)"),
        x     = NULL,
        y     = "Mean probability"
      ) +
      theme_classic(base_size = 11) +
      theme(
        axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
        plot.title       = element_text(face = "bold", hjust = 0.5),
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text       = element_text(face = "bold")
      )

    out_bar <- file.path(OUTDIR,
                         paste0(disease, "_", ct, "_soft_composition_bar_merged.pdf"))
    ggsave(out_bar, p_bar,
           width  = max(6, n_distinct(bar_df$patient_id) * 0.7),
           height = 4, dpi = 300)
    cat("  Saved:", out_bar, "\n")

    # ── 3 niche Spearman 联合散点图 ───────────────────────────────────────
    make_scatter <- function(y_col, niche_name, color) {
      sp_row <- sp_all[sp_all$niche == y_col, ]
      lbl    <- sprintf("ρ = %.3f\np = %.4f", sp_row$rho, sp_row$p.value)
      y_max  <- max(scores[[y_col]], na.rm = TRUE)
      ggplot(scores, aes_string(x = "grade_num", y = y_col)) +
        geom_smooth(method = "lm", se = TRUE, color = "grey60", fill = "grey85",
                    linewidth = 0.7) +
        geom_jitter(width = 0.05, size = 2.5, alpha = 0.9, color = color) +
        annotate("text", x = min(scores$grade_num) + 0.1, y = y_max,
                 label = lbl, hjust = 0, vjust = 1, size = 3) +
        scale_x_continuous(breaks = seq_along(grade_order), labels = grade_order) +
        labs(title = niche_name, x = grade_label,
             y = paste0("Mean prob(", niche_name, ")")) +
        theme_classic(base_size = 11) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5, color = color))
    }

    panel <- make_scatter("mean_prob_T",     "T_Dominant",      niche_colors["T_Dominant"]) +
             make_scatter("mean_prob_Mixed",  "Mixed",           niche_colors["Mixed"]) +
             make_scatter("mean_prob_Plasma", "Plasma_Dominant", niche_colors["Plasma_Dominant"]) +
      plot_annotation(
        title = paste0(disease, " — ", ct, " soft niche scores vs ",
                       grade_label, " (merged)"),
        theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))
      )

    out_panel <- file.path(OUTDIR,
                           paste0(disease, "_", ct, "_3niche_spearman_panel_merged.pdf"))
    ggsave(out_panel, panel, width = 12, height = 4, dpi = 300)
    cat("  Saved:", out_panel, "\n")
  }

  # ── 合并大图 ───────────────────────────────────────────────────────────────
  if (length(scatter_T_list) == 2) {
    combined <- wrap_plots(scatter_T_list, ncol = 2) +
      plot_annotation(
        title = paste0(disease, " — Mean prob(T_Dominant) vs ", grade_label,
                       " (3-class, true patients)"),
        theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))
      )
    out_all <- file.path(OUTDIR, paste0(disease, "_all_summary_merged.pdf"))
    ggsave(out_all, combined, width = 10, height = 4, dpi = 300)
    cat("\n  Combined figure saved:", out_all, "\n")
  }
}

cat("\n", strrep("=", 60), "\n")
cat("All done. Outputs in:", OUTDIR, "\n")
cat(strrep("=", 60), "\n")
