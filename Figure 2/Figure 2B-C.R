
# Load libraries
library(ggplot2)
library(dplyr)

# Read CSV file
df <- read.csv("PERMANOVA.csv")

# Define color palette
mycol <- c("#E2B0BA","#87AFD5","#6EB9C3","#C98B88",
           "#93C89A","#FFCC98","#E1D1BA","#EDAFA9",
           "#0AB0C8","#9781BB","#E8BD65","#E39844","#ADD1E5")

darken_colors <- function(colors, factor = 1.3) {
  sapply(colors, function(col) {
    rgb_col <- col2rgb(col)
    hsv_col <- rgb2hsv(rgb_col)
    hsv_col[2] <- pmin(hsv_col[2] * factor, 1)
    hsv(hsv_col[1], hsv_col[2], hsv_col[3])
  })
}

# Darken the colors
mycol_dark <- darken_colors(mycol, factor = 1.5)

# Select required columns and sort by R2 in ascending order
# Column order: Group, Factor_Count, Group_Marginal_R2, p_value, p_adj
# Therefore: Group = column 1, R2 = column 3, p_adj = column 5
plot_data <- df %>%
  select(Group = 1, R2 = 3, p_adj = 5) %>%
  arrange(R2) %>%
  mutate(Group = factor(Group, levels = unique(Group))) %>%
  mutate(
    significance = case_when(
      p_adj < 0.05 ~ "*"
    )
  ) %>%
  mutate(
    row_num = row_number(),
    point_color = mycol_dark[(row_num - 1) %% length(mycol_dark) + 1]
  )

# Draw horizontal lollipop plot
p <- ggplot(plot_data, aes(x = R2, y = Group)) +
  geom_col(fill = "gray", width = 0.1, alpha = 0.7) +
  geom_point(aes(color = Group), size = 12, show.legend = FALSE) +
  scale_color_manual(values = setNames(mycol_dark[1:nrow(plot_data)], levels(plot_data$Group))) +
  geom_text(
    aes(label = significance),
    nudge_x = max(plot_data$R2) * 0.05,
    nudge_y = -0.14,
    hjust = 0.2,
    size = 9,
    color = "black"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.text.x = element_text(size = 13, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    plot.title = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.margin = margin(15, 15, 5, 5)
  ) +
  labs(
    x = expression(Variance~explained~(R^2)),
    y = ""
  ) +
  xlim(0, max(plot_data$R2) * 1.15) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)))

# Display plot
print(p)

# Export as PDF
ggsave("urbanization_dsv.pdf", plot = p, width = 8, height = 6, device = "pdf")
```
