
# Load libraries
library(ggplot2)
library(dplyr)

# Read CSV file
df <- read.csv("permanova.csv")

# Define color palette
mycol <- c("#E2B0BA","#87AFD5","#6EB9C3","#C98B88",
           "#93C89A","#FFCC98","#E1D1BA","#EDAFA9",
           "#0AB0C8","#9781BB","#E8BD65","#E39844","#ADD1E5")

darken_colors <- function(colors, factor = 1.3) {
  sapply(colors, function(col) {
    rgb_col <- col2rgb(col)
    hsv_col <- rgb2hsv(rgb_col)
    # Increase saturation
    hsv_col[2] <- pmin(hsv_col[2] * factor, 1)
    hsv(hsv_col[1], hsv_col[2], hsv_col[3])
  })
}

# Darken the colors
mycol_dark <- darken_colors(mycol, factor = 1.5)

# Select required columns and sort by R2 in ascending order
plot_data <- df %>%
  select(Group = 1, R2 = 3, p_adj = 5) %>%
  # Sort by R2 in ascending order
  arrange(R2) %>%
  # Convert Group to a factor using the sorted order
  mutate(Group = factor(Group, levels = unique(Group))) %>%
  # Assign significance symbols based on adjusted P values
  mutate(
    significance = case_when(
      p_adj < 0.05 ~ "*"
    )
  ) %>%
  # Assign colors to each group
  mutate(row_num = row_number(),
         point_color = mycol_dark[(row_num - 1) %% length(mycol_dark) + 1])

# Draw horizontal lollipop plot ordered by R2
p <- ggplot(plot_data, aes(x = R2, y = Group)) +
  # Draw horizontal gray bars
  geom_col(fill = "gray", width = 0.1, alpha = 0.7) +
  # Add colored points
  geom_point(aes(color = Group), size = 12, show.legend = FALSE) +
  # Apply custom colors
  scale_color_manual(values = setNames(mycol_dark[1:nrow(plot_data)], levels(plot_data$Group))) +
  # Add significance symbols
  geom_text(aes(label = significance),
            nudge_x = max(plot_data$R2) * 0.05,
            nudge_y = -0.14,
            hjust = 0.2,
            size = 9, color = "black") +
  theme_minimal() +
  theme(
    # Remove background grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Axis text and titles
    axis.text.y = element_text(size = 13, color = "black"),
    axis.text.x = element_text(size = 13, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    # Remove plot title
    plot.title = element_blank(),
    # Set axis lines and ticks to black
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    # Adjust plot margins
    plot.margin = margin(15, 15, 5, 5)
  ) +
  labs(
    x = expression(Variance~explained~(R^2)),
    y = ""
  ) +
  # Adjust x-axis range to keep significance symbols visible
  xlim(0, max(plot_data$R2) * 1.15) +
  # Adjust y- and x-axis expansion
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)))

# Display plot
print(p)

# Export as PDF
ggsave("r2_sorted_ascending.pdf", plot = p, width = 8, height = 6, device = "pdf")
```
