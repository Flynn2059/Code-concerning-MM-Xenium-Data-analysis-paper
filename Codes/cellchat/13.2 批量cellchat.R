Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 16 * 1024^3)  
rm(list = ls())
gc()

setwd("/Volumes/FlynnBio/XinyeLi/GSE299207_RAW/")

# ---- packages ----
library(Seurat)
library(data.table)
library(dplyr)
library(tidyverse)
library(qs)
library(SpatialCellChat)
library(Matrix)
library(RANN)
library(spatstat.sparse)

# =========================================================
# patch SpatialCellChat internal future helpers to pure serial
# avoids:
# 1. "should not be called with handlers on the stack"
# 2. future.globals.maxSize overflow in centrality calculation
# =========================================================
patch_spatialcellchat_serial <- function() {
  my_future_lapply_serial <- function(
    X, FUN, ..., future.seed = TRUE, simplify = FALSE,
    hint.message = "Computing...", strategy.message = TRUE,
    future.label = "future_lapply-%d"
  ) {
    if (strategy.message) {
      cat("\u2139 Future strategy in use: `sequential`\n")
    }
    res <- lapply(X, FUN, ...)
    if (simplify) res <- unlist(res)
    res
  }
  
  my_future_sapply_serial <- function(
    X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE,
    future.envir = parent.frame(), future.seed = TRUE,
    future.label = "future_sapply-%d",
    hint.message = "Computing...", strategy.message = TRUE
  ) {
    if (strategy.message) {
      cat("\u2139 Future strategy in use: `sequential`\n")
    }
    ans <- lapply(X, FUN, ...)
    if (USE.NAMES && is.character(X) && is.null(names(ans))) {
      names(ans) <- X
    }
    if (!isFALSE(simplify) && length(ans)) {
      simplify2array(ans, higher = (identical(simplify, "array")))
    } else {
      ans
    }
  }
  
  ns <- asNamespace("SpatialCellChat")
  
  unlockBinding("my_future_lapply", ns)
  assign("my_future_lapply", my_future_lapply_serial, envir = ns)
  lockBinding("my_future_lapply", ns)
  
  unlockBinding("my_future_sapply", ns)
  assign("my_future_sapply", my_future_sapply_serial, envir = ns)
  lockBinding("my_future_sapply", ns)
}

patch_spatialcellchat_serial()

# ---- input ----
seurat_all <- qread("all_sample_plasma_and_t_cells.qs")
seurat_list <- SplitObject(seurat_all, split.by = "orig.ident")

# ---- output dir ----
dir.create("SpatialCellChat_results", showWarnings = FALSE)

# =========================================================
# helper 1: find raw Xenium directory for one sample
# =========================================================
find_sample_dir <- function(sample_name, root_dir = getwd()) {
  hits <- Sys.glob(file.path(root_dir, paste0("*_", sample_name)))
  hits <- hits[dir.exists(hits)]
  
  if (length(hits) != 1) {
    stop("Cannot uniquely find raw Xenium directory for sample: ", sample_name)
  }
  
  hits
}

# =========================================================
# helper 2: extract coordinates from raw cells.csv.gz
# Seurat cell names: GSM9035042_SM-4_aaaadojg-1
# raw cell_id:       aaaadojg-1
# =========================================================
get_sample_coordinates <- function(seurat_obj, sample_name, root_dir = getwd()) {
  sample_dir <- find_sample_dir(sample_name, root_dir)
  cell_file <- file.path(sample_dir, "cells.csv.gz")
  
  raw_cells <- fread(
    cell_file,
    select = c("cell_id", "x_centroid", "y_centroid")
  )
  
  cell_map <- data.frame(
    cell = colnames(seurat_obj),
    cell_id = sub("^[^_]+_[^_]+_", "", colnames(seurat_obj)),
    stringsAsFactors = FALSE
  )
  
  coords <- merge(
    cell_map,
    raw_cells,
    by = "cell_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  coords <- coords[match(colnames(seurat_obj), coords$cell), ]
  
  if (any(is.na(coords$x_centroid)) || any(is.na(coords$y_centroid))) {
    bad_n <- sum(is.na(coords$x_centroid) | is.na(coords$y_centroid))
    stop("Unmatched cells in cells.csv.gz for sample ", sample_name, ": ", bad_n)
  }
  
  coords_mat <- as.data.frame(coords[, c("x_centroid", "y_centroid")])
  colnames(coords_mat) <- c("x", "y")
  rownames(coords_mat) <- coords$cell
  
  stopifnot(identical(rownames(coords_mat), colnames(seurat_obj)))
  
  coords_mat
}

# =========================================================
# helper 3: estimate scale.distance
# target scaled minimum nearest-neighbor distance ~ 1.2
# =========================================================
estimate_scale_distance <- function(coords, target_min = 1.2) {
  nn <- RANN::nn2(as.matrix(coords), as.matrix(coords), k = 2)
  min_nn <- min(nn$nn.dists[, 2][nn$nn.dists[, 2] > 0])
  
  if (!is.finite(min_nn) || min_nn <= 0) {
    stop("Failed to estimate nearest-neighbor distance from coordinates")
  }
  
  target_min / min_nn
}

# =========================================================
# helper 4: fix SpatialCellChat aggregateNet bug on sparse3Darray
# =========================================================
aggregateNet_fixed <- function(object, thresh = 0.05) {
  net <- object@net
  
  prob <- net$prob
  pval <- net$pval
  
  prob_arr <- if (inherits(prob, "sparse3Darray")) as.array(prob) else prob
  pval_arr <- if (inherits(pval, "sparse3Darray")) as.array(pval) else pval
  
  pval_arr[prob_arr == 0] <- 1
  prob_arr[pval_arr >= thresh] <- 0
  
  net$count <- apply(prob_arr > 0, c(1, 2), sum)
  net$weight <- apply(prob_arr, c(1, 2), sum)
  
  net$count[is.na(net$count)] <- 0
  net$weight[is.na(net$weight)] <- 0
  net$LR.sig <- dimnames(prob_arr)[[3]][apply(prob_arr, 3, sum) > 0]
  
  if ("prob.cell" %in% names(net)) {
    prob.cell <- net$prob.cell
    prob.cell.positive <- prob.cell
    prob.cell.positive$x <- rep.int(1, length(prob.cell.positive$x))
    
    net$count.cell <- spatstat.sparse::marginSumsSparse(
      prob.cell.positive,
      MARGIN = c(1, 2)
    )
    
    net$weight.cell <- spatstat.sparse::marginSumsSparse(
      prob.cell,
      MARGIN = c(1, 2)
    )
    
    prob.cell.sum <- spatstat.sparse::marginSumsSparse(
      prob.cell,
      MARGIN = c(3)
    )
    
    net$LR.sig.cell <- dimnames(prob.cell)[[3]][prob.cell.sum@i]
  }
  
  object@net <- net
  object
}

# =========================================================
# helper 5: run one sample
# =========================================================
run_spatialcellchat_one <- function(
    seurat_obj,
    sample_name,
    root_dir = getwd(),
    assay.use = "Xenium",
    group.by = "major_annotation",
    species = "human",
    min.cells = 10,
    ratio = 1,
    tol = 5,
    interaction.range = 250,
    contact.range = 10,
    use.projected = FALSE
) {
  DefaultAssay(seurat_obj) <- assay.use
  
  # ---- normalize ----
  ncount_col <- paste0("nCount_", assay.use)
  seurat_obj <- NormalizeData(
    seurat_obj,
    assay = assay.use,
    scale.factor = median(seurat_obj@meta.data[[ncount_col]])
  )
  
  # ---- clean grouping variable ----
  if (!group.by %in% colnames(seurat_obj@meta.data)) {
    stop("group.by column not found in meta.data: ", group.by)
  }
  
  group_vec <- seurat_obj@meta.data[[group.by]]
  keep_groups <- names(which(table(group_vec) >= min.cells))
  keep_cells <- rownames(seurat_obj@meta.data)[
    !is.na(group_vec) & group_vec %in% keep_groups
  ]
  
  seurat_obj <- subset(seurat_obj, cells = keep_cells)
  seurat_obj$cellchat_group <- droplevels(factor(seurat_obj@meta.data[[group.by]]))
  
  # ---- coordinates from raw Xenium output ----
  coords <- get_sample_coordinates(
    seurat_obj = seurat_obj,
    sample_name = sample_name,
    root_dir = root_dir
  )
  
  # ---- Xenium spatial factors ----
  spatial.factors <- data.frame(ratio = ratio, tol = tol)
  
  # ---- auto scale.distance ----
  scale.distance <- estimate_scale_distance(coords, target_min = 1.2)
  message(sample_name, " scale.distance = ", signif(scale.distance, 4))
  
  # ---- create SpatialCellChat object ----
  cellchat <- createSpatialCellChat(
    object = seurat_obj,
    group.by = "cellchat_group",
    assay = assay.use,
    datatype = "spatial",
    coordinates = coords,
    spatial.factors = spatial.factors
  )
  
  # ---- database ----
  if (species == "human") {
    cellchat@DB <- CellChatDB.human
    PPI.use <- PPI.human
  } else if (species == "mouse") {
    cellchat@DB <- CellChatDB.mouse
    PPI.use <- PPI.mouse
  } else {
    stop("species must be 'human' or 'mouse'")
  }
  
  # ---- main pipeline ----
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  if (use.projected) {
    cellchat <- projectData(cellchat, PPI.use)
  }
  
  cellchat <- computeCommunProb(
    cellchat,
    raw.use = !use.projected,
    distance.use = TRUE,
    interaction.range = interaction.range,
    contact.range = contact.range,
    contact.dependent = TRUE,
    scale.distance = scale.distance
  )
  
  cellchat <- filterCommunication(cellchat, min.cells = min.cells)
  cellchat <- computeCommunProbPathway(cellchat)
  
  # replace buggy aggregateNet()
  cellchat <- aggregateNet_fixed(cellchat)
  
  # centrality now uses patched serial helper
  cellchat <- netAnalysis_computeCentrality(
    cellchat,
    slot.name = "netP"
  )
  
  return(cellchat)
}

# =========================================================
# batch run
# no tryCatch around the main loop
# =========================================================

sample_order <- names(seurat_list)

# 先试单样本时可改成：
# sample_order <- "SM-4"

for (nm in sample_order) {
  out_file <- file.path("SpatialCellChat_results", paste0("cellchat_", nm, ".qs"))
  
  if (file.exists(out_file)) {
    message("Skip existing result: ", nm)
    next
  }
  
  message("Running SpatialCellChat for: ", nm)
  
  cellchat_obj <- run_spatialcellchat_one(
    seurat_obj = seurat_list[[nm]],
    sample_name = nm,
    root_dir = "/Volumes/FlynnBio/XinyeLi/GSE299207_RAW/",
    assay.use = "Xenium",
    group.by = "major_annotation",
    species = "human",
    min.cells = 10,
    ratio = 1,
    tol = 5,
    interaction.range = 250,
    contact.range = 10,
    use.projected = FALSE
  )
  
  qsave(cellchat_obj, out_file)
  rm(cellchat_obj)
  gc()
}

# =========================================================
# combine all saved sample results into one list
# =========================================================

res_files <- list.files(
  "SpatialCellChat_results",
  pattern = "^cellchat_.*\\.qs$",
  full.names = TRUE
)

cellchat_list <- lapply(res_files, qread)
names(cellchat_list) <- gsub("^cellchat_|\\.qs$", "", basename(res_files))

qsave(cellchat_list, "cellchat_list_by_sample.qs")

# =========================================================
# export communication tables
# =========================================================

comm_df_list <- lapply(cellchat_list, subsetCommunication)
qsave(comm_df_list, "cellchat_communication_df_by_sample.qs")

comm_df_all <- bind_rows(
  lapply(names(comm_df_list), function(nm) {
    df <- comm_df_list[[nm]]
    df$sample <- nm
    df
  })
)

fwrite(comm_df_all, "cellchat_communication_all_samples.tsv.gz", sep = "\t")