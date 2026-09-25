# Load required packages
library(pheatmap)
library(dplyr)
library(scico)

# Read data
r_data <- read.csv("R.csv", row.names = 1, check.names = FALSE)
significant_data <- read.csv("R-p.csv", row.names = 1, check.names = FALSE)

# No Urban/Rural ordering is required; keep the original row order
# Ensure that the row order of the significance data matches the main data
significant_data <- significant_data[rownames(r_data), ]

# Calculate the variation of each column using standard deviation
# and sort columns from low to high variation
col_variation <- apply(r_data, 2, function(x) sd(x, na.rm = TRUE))
ordered_indices <- order(col_variation)
r_data <- r_data[, ordered_indices]
significant_data <- significant_data[, ordered_indices]

# Convert P values to significance symbols
# *** p < 0.001, ** p < 0.01, * p < 0.05
sig_marks <- matrix("", nrow = nrow(significant_data), ncol = ncol(significant_data))
sig_marks[significant_data < 0.001] <- "***"
sig_marks[significant_data >= 0.001 & significant_data < 0.01] <- "**"
sig_marks[significant_data >= 0.01 & significant_data < 0.05] <- "*"

# Color palette: blue-white-pink
target_palette <- c("#8ECAE6", "#FFFFFF", "#F4ACB7")
color_count <- 101
my_color <- colorRampPalette(target_palette)(color_count)

# Define data range and color breaks
data_min <- min(r_data, na.rm = TRUE)
data_max <- max(r_data, na.rm = TRUE)
ph_breaks <- seq(data_min, data_max, length.out = color_count + 1)

# Draw the heatmap and export as PDF
pdf("heatmap_output_with_significance.pdf", width = 12, height = 8)

pheatmap(
  as.matrix(r_data),
  color = my_color,
  breaks = ph_breaks,
  cluster_rows = FALSE,        # Keep the original row order
  cluster_cols = FALSE,        # Keep columns ordered by variation
  cellwidth = 30,
  cellheight = 30,
  treeheight_row = 0,          # Hide row clustering tree
  treeheight_col = 0,          # Hide column clustering tree
  annotation_names_row = FALSE,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 90,              # Display column names vertically
  na_col = "grey80",           # Display NA values in gray
  legend_labels = "log2 OR",   # Legend label
  display_numbers = sig_marks, # Display significance symbols
  number_color = "black",      # Significance symbol color
  number_size = 8              # Significance symbol size
)

dev.off()

# Display the heatmap in R
pheatmap(
  as.matrix(r_data),
  color = my_color,
  breaks = ph_breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cellwidth = 30,
  cellheight = 30,
  treeheight_row = 0,
  treeheight_col = 0,
  annotation_names_row = FALSE,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 90,
  na_col = "grey80",
  legend_labels = "log2 OR",
  display_numbers = sig_marks,
  number_color = "black",
  number_size = 8
)
```
