# scripts/run_example.R

# Load required libraries
library(readr)
library(dplyr)
library(ggplot2)

# Load BSArray functions
source("BSArray_v1.0.0.R")

# -----------------------------
# Input files
# -----------------------------
bulk_file <- "example_data/Rcs2_bulk_matrix.txt"
parent_file <- "example_data/Rcs2_parent_raw.txt"

# -----------------------------
# Run BSArray
# -----------------------------
result <- run_bsarray(
  bulk_file = bulk_file,
  parent_file = parent_file,
  window_size = 1e6,
  step_size = 1e5,
  min_snps = 20
)

# -----------------------------
# Save outputs
# -----------------------------
dir.create("example_output", showWarnings = FALSE)

write.csv(result$all_results,
          "example_output/all_results.csv",
          row.names = FALSE)

ggsave("example_output/landscape_plot.png",
       result$plot,
       width = 10, height = 6)

# -----------------------------
# Message
# -----------------------------
cat("✅ BSArray example completed successfully\n")