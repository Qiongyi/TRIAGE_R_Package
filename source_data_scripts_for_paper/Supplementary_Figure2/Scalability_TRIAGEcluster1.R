library(reticulate)
setwd("/home/uqqzhao/QY/others/202410_BiB_revision/PBMC_ifnb_data") 
results_file <- "TRIAGEcluster_runtime_results.txt"

if (!file.exists(results_file)) {
  cat("Cell_number\tRuntime (mins)\n", file = results_file)
}

# sample sizes
sample_sizes <- c(5000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000)
# Create R scripts for each sample size (number of cells)
for (num_cells in sample_sizes) {
  cat("Creating script for", num_cells, "samples...\n")
  formatted_num <- format(num_cells, scientific = FALSE, trim = TRUE)
  
  # Script for each sample size
  script_content <- paste0(
    "library(reticulate)\n",
    "use_condaenv(\"r-reticulate\", required = TRUE)\n",
    "library(TRIAGE)\n",
    "setwd(\"/home/uqqzhao/QY/others/202410_BiB_revision/PBMC_ifnb_data\")\n",
    
    "start_time <- Sys.time()\n",
    # Run TRIAGEcluster
    "TRIAGEcluster(expr = \"sample_expr_data_", formatted_num, ".csv\",\n",
    "              metadata = \"sample_metadata_", formatted_num, ".csv\",\n",
    "              outdir= \"TRIAGEcluster_", formatted_num, "_cells\",\n",
    "              output_prefix = \"TRIAGEcluster\",\n",
    "              cell_column = \"cell_name\",\n",
    "              umap_column = \"UMAP_\",\n",
    "              seed=37,\n",
    "              bw = 0.4)\n",
    "end_time <- Sys.time()\n",
    
    "run_time <- end_time - start_time\n",
    "run_time_mins <- as.numeric(run_time, units = 'mins')\n",
    
    "cat(\"Start time:\", start_time, \"\\n\")\n",
    "cat(\"End time:\", end_time, \"\\n\")\n",
    "cat(\"Runtime (TRIAGEcluster for ", formatted_num, " cells):\", run_time_mins, \"minutes\\n\")\n",
    
    "cat(\"", formatted_num, "\\t\", run_time_mins, \"\\n\", file=\"", results_file, "\", append=TRUE)\n"
  )
  
  script_file <- paste0("TRIAGEcluster_", formatted_num, "_cells.R")
  writeLines(script_content, script_file)
}

# conda activate r-env
# R --vanilla < /home/uqqzhao/QY/others/202410_BiB_revision/PBMC_ifnb_data/Scalability_TRIAGEcluster1.R > /dev/null 2>&1