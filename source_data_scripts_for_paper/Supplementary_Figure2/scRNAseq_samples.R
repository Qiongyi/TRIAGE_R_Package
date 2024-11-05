library(data.table)

setwd("/home/uqqzhao/QY/others/202410_BiB_revision/PBMC_ifnb_data")

# Load the gene expression data and metadata
expr_data <- fread("pbmc_ifnb_Seurat_CCAint_allgenes.csv")
metadata <- fread("pbmc_ifnb_metadata_withUMAP.csv")

# sample sizes
#target_sizes <- c(20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000)
target_sizes <- c(5000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000)

original_cell_count <- ncol(expr_data) - 1
expanded_expr_data <- expr_data
expanded_metadata <- metadata

for (i in 1:8) {
  copy_expr_data <- copy(expr_data)
  copy_metadata <- copy(metadata)
  new_cell_names <- paste0(colnames(copy_expr_data)[-1], ".", i)
  setnames(copy_expr_data, old = colnames(copy_expr_data)[-1], new = new_cell_names)
  copy_metadata[, V1 := new_cell_names]
  noise_level <- 0.01
  copy_expr_data[, (new_cell_names) := lapply(.SD, function(x) x + rnorm(length(x), mean = 0, sd = noise_level * abs(x))), .SDcols = new_cell_names]
  expanded_expr_data <- cbind(expanded_expr_data, copy_expr_data[, -1, with = FALSE])
  expanded_metadata <- rbind(expanded_metadata, copy_metadata)
}

for (size in target_sizes) {
  cat("Processing sample size:", size, "...\n")
  sample_needed <- size - original_cell_count
  formatted_size <- format(size, scientific = FALSE, trim = TRUE)

  if (size <= original_cell_count) {
    sampled_cells <- sample(colnames(expr_data)[-1], size, replace = FALSE)
    final_expr_data <- expr_data[, c("V1", sampled_cells), with = FALSE]
    final_metadata <- metadata[V1 %in% sampled_cells]
  } else {
    sample_needed <- size - original_cell_count
    sampled_cells <- sample(colnames(expanded_expr_data)[-(1:(original_cell_count + 1))], sample_needed)
    final_cells <- c(colnames(expr_data)[-1], sampled_cells)
    final_expr_data <- expanded_expr_data[, c("V1", final_cells), with = FALSE]
    final_metadata <- expanded_metadata[V1 %in% final_cells]
  }  
  
  final_metadata$cell_name <- final_metadata$V1
  expr_filename <- paste0("sample_expr_data_", formatted_size, ".csv")
  metadata_filename <- paste0("sample_metadata_", formatted_size, ".csv")
  
  fwrite(final_expr_data, expr_filename)
  fwrite(final_metadata, metadata_filename)
  
  cat("Sample size", size, "data saved to", expr_filename, "and", metadata_filename, "\n")
}
