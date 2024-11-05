# Set working directory
setwd("/home/uqqzhao/QY/others/202410_BiB_revision/GSE95755") 
results_file <- "TRIAGEgene_runtime_results.txt"

if (!file.exists(results_file)) {
  cat("Sample_size\tRuntime (mins)\n", file = results_file)
}

sample_sizes <- c(4, 10, 50, 100, 500, 1000, 2000, 3000, 4000, 5000)
# Create R scripts for each sample size
for (num_samples in sample_sizes) {
  cat("Creating script for", num_samples, "samples...\n")
  formatted_num_samples <- format(num_samples, scientific = FALSE, trim = TRUE)
  script_content <- paste0(
    "library(reticulate)\n",
    "use_condaenv(\"r-reticulate\", required = TRUE)\n",
    "library(TRIAGE)\n",
    "setwd(\"/home/uqqzhao/QY/others/202410_BiB_revision/GSE95755\")\n",
    "# Load original dataset\n",
    "file <- read.table(\"GSE95755_MultiCellularRNAseq_EdgeR_CPM.txt\", header=T, quote=\"\", sep=\"\\t\")\n",
    "cpm <- file[,c(\"Symbol\", \"MIP56_Myo_1\", \"MIP56_Myo_2\", \"MIP56_Myo_3\", \"MIP56_Myo_4\")]\n",
    "rownames(cpm) <- cpm$Symbol\n",
    "cpm <- cpm[,-1]\n",
    "simulate_samples <- function(data, num_samples) {\n",
    "  if (num_samples <= ncol(data)) {\n",
    "    return(data[, 1:num_samples])\n",
    "  } else {\n",
    "    additional_samples <- num_samples - ncol(data)\n",
    "    extended_data <- data.frame(data, replicate(additional_samples, sample(unlist(data), nrow(data), replace = TRUE)))\n",
    "    colnames(extended_data) <- paste0(\"Sample_\", 1:ncol(extended_data))\n",
    "    return(extended_data)\n",
    "  }\n",
    "}\n",
    "simulated_cpm <- simulate_samples(cpm, ", formatted_num_samples, ")\n",
    "# Measure runtime for TRIAGEgene\n",
    "start_time <- Sys.time()\n",
    "ds <- TRIAGEgene(simulated_cpm, species = \"Mouse\", pvalue = T)\n",
    "end_time <- Sys.time()\n",
    "run_time <- end_time - start_time\n",
    "run_time_mins <- as.numeric(run_time, units = 'mins')\n",
    "# write.table(ds, \"DS_table_", formatted_num_samples, "_samples.xls\", sep=\"\\t\", row.names=F, quote=F)\n",
    "cat(\"Start time:\", start_time, \"\\n\")\n",
    "cat(\"End time:\", end_time, \"\\n\")\n",
    "cat(\"Runtime (TRIAGEgene for ", formatted_num_samples, " samples):\", run_time_mins, \"\\n\")\n",
    "cat(\"", formatted_num_samples, "\\t\", run_time_mins, \"\\n\", file=\"", results_file, "\", append=TRUE)\n"
  )
  
  script_file <- paste0("TRIAGEgene_", formatted_num_samples, "_samples.R")
  writeLines(script_content, script_file)
}

# conda activate r-env
# R --vanilla < /home/uqqzhao/QY/others/202410_BiB_revision/GSE95755/Scalability_TRIAGEgene1.R > /dev/null 2>&1

