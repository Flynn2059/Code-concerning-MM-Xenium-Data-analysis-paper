# ==== 调整自己喜欢的系统设置 ====
Sys.setenv(LANGUAGE="en")
options(StringAsFactors=FALSE)
rm(list=ls());gc()
setwd("/Volumes/FlynnDisk/XinyeLi/GSE299207_RAW/")
getwd()
list.files()
list.files("./RCTD/")
# ==== 加载我们要用到的包 ====
library(ggplot2)
library(data.table)
library(dplyr)
library(qs)
library(RColorBrewer)
library(spacexr)

# ==== 设置路径 ====
rctd_dir <- "./RCTD/"

# 获取所有 qs 文件
qs_files <- list.files(
  path = rctd_dir,
  pattern = "\\.qs$",
  full.names = TRUE
)

print(qs_files)

# ==== 遍历处理 ====
for (qs_file in qs_files) {
  
  cat("Processing:", qs_file, "\n")
  
  # 读取 RCTD 对象
  RCTD_obj <- qread(qs_file)
  
  # 提取结果
  result <- RCTD_obj@results$results_df
  
  # 过滤 reject
  result <- result[result$spot_class != "reject", ]
  
  # 构建输出 dataframe
  result_df <- data.frame(
    row.names = rownames(result),
    celltype = result$first_type
  )
  
  # 构建输出文件名
  csv_file <- file.path(
    rctd_dir,
    paste0(tools::file_path_sans_ext(basename(qs_file)), ".csv")
  )
  
  # 写出 csv
  write.csv(result_df, file = csv_file)
  
  cat("Saved:", csv_file, "\n")
}

cat("All files processed.\n")
