
# Load required packages
library(pheatmap)
library(dplyr)
library(scico)

# Read data
r_data <- read.csv("Odds ratio log2.csv", row.names = 1, check.names = FALSE)
group_col <- read.csv("classification.csv", row.names = 1, check.names = FALSE)
significant_data <- read.csv("significant.csv", row.names = 1, check.names = FALSE)

# Set row order: Urban and Rural first, followed by the remaining ethnic groups
ordered_rows <- c("Urban", "Rural", setdiff(rownames(r_data), c("Urban", "Rural")))
r_data <- r_data[ordered_rows, ]

# Ensure the row order of the significance data matches the heatmap data
significant_data <- significant_data[ordered_rows, ]

# ----------------------
# Calculate variation for each column using standard deviation
# Larger values indicate greater variation
# ----------------------

col_variation <- apply(r_data, 2, function(x) sd(x, na.rm = TRUE))

# Sort columns from low to high variation
ordered_indices <- order(col_variation)
r_data <- r_data[, ordered_indices]
significant_data <- significant_data[, ordered_indices]

# Urban and Rural occupy the first two rows
gap_position <- 2

# ----------------------
# Convert P values to significance symbols
# ----------------------

# Significance levels:
# *** p < 0.001, ** p < 0.01, * p < 0.05
sig_marks <- matrix(
  "",
  nrow = nrow(significant_data),
  ncol = ncol(significant_data)
)

sig_marks[significant_data < 0.001] <- "***"
sig_marks[significant_data >= 0.001 & significant_data < 0.01] <- "**"
sig_marks[significant_data >= 0.01 & significant_data < 0.05] <- "*"

# ----------------------
# Color settings
# ----------------------

# Main blue-white-pink color palette
target_palette <- c("#8ECAE6", "#FFFFFF", "#F4ACB7")
color_count <- 101
my_color <- colorRampPalette(target_palette)(color_count)

# Define data range and color breaks
data_min <- min(r_data, na.rm = TRUE)
data_max <- max(r_data, na.rm = TRUE)
ph_breaks <- seq(data_min, data_max, length.out = color_count + 1)

# Group colors
group_colors <- list(
  group = setNames(
    c("#D885A3", "#FBD38D", "#90CDF4", "#FAF089", "#A7F3D0", "#F97316"),
    nm = group_col$group %>% unique() %>% sort()
  )
)

# ----------------------
# Draw heatmap and export as PDF
# ----------------------

pdf(
  "heatmap_output_with_significance.pdf",
  width = 12,
  height = 8
)

pheatmap(
  as.matrix(r_data),
  color = my_color,
  breaks = ph_breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = gap_position,
  gap_width = 2,
  cellwidth = 10,
  cellheight = 10,
  treeheight_row = 0,
  treeheight_col = 0,
  annotation_col = group_col,
  annotation_colors = group_colors,
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

dev.off()

# Display the heatmap in R
pheatmap(
  as.matrix(r_data),
  color = my_color,
  breaks = ph_breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = gap_position,
  gap_width = 2,
  cellwidth = 10,
  cellheight = 10,
  treeheight_row = 0,
  treeheight_col = 0,
  annotation_col = group_col,
  annotation_colors = group_colors,
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
