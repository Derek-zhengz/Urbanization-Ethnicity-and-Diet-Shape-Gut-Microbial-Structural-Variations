
library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(extrafont)

df <- read.csv("counts.csv")

# Rename columns to match the Excel headers
colnames(df) <- c("Genus", "vsv_counts", "sample_dsgv_counts")

# Format genus names: abbreviate the part before "_" to the first letter
df$Genus <- str_replace_all(df$Genus, "([A-Za-z]+)_([A-Za-z_]+)", "\\1.\\2")
df$Genus <- str_replace(df$Genus, "^([A-Za-z])([A-Za-z]+)\\.", "\\1.")

# Convert vsv_counts to integer
df$vsv_counts <- as.integer(df$vsv_counts)

# Check and merge duplicated genus names
if (any(duplicated(df$Genus))) {
  df <- df %>%
    group_by(Genus) %>%
    summarise(
      vsv_counts = sum(vsv_counts, na.rm = TRUE),
      sample_dsgv_counts = sum(sample_dsgv_counts, na.rm = TRUE)
    ) %>%
    ungroup()
  cat("Duplicated genus names were merged by summing their counts.\n")
}

# Calculate the total number of SVs for each genus
df <- df %>%
  mutate(Total_SV = vsv_counts + sample_dsgv_counts)

# Select the top 50 genera ranked by Total_SV
top_70_df <- df %>%
  arrange(desc(Total_SV)) %>%
  slice_head(n = 50)

# Reshape the filtered data to long format
df_long <- top_70_df %>%
  pivot_longer(
    cols = c("vsv_counts", "sample_dsgv_counts"),
    names_to = "SV_Type",
    values_to = "Count"
  )

# Rename SV types
df_long$SV_Type <- ifelse(df_long$SV_Type == "vsv_counts", "Variable SVs", "Deletion SVs")

# Order genera by total SV count
df_long$Genus <- factor(
  df_long$Genus,
  levels = top_70_df %>%
    arrange(desc(Total_SV)) %>%
    pull(Genus)
)

# Draw stacked bar plot
bar_plot <- ggplot(df_long, aes(x = Genus, y = Count, fill = SV_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c(
    "Variable SVs" = "#8ECAE6",
    "Deletion SVs" = "#F4ACB7"
  )) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    
    # Legend style
    legend.title = element_text(
      size = 14,
      color = "black",
      family = "sans"
    ),
    legend.text = element_text(
      size = 12,
      color = "black",
      family = "sans"
    ),
    legend.key.size = unit(1, "cm"),
    
    # X-axis text style
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      face = "italic",
      size = 12,
      margin = margin(t = 2),
      color = "black",
      family = "sans"
    ),
    axis.title.x = element_text(
      size = 0,
      margin = margin(t = 20),
      color = "black",
      family = "sans"
    ),
    
    # Y-axis text style
    axis.text.y = element_text(
      size = 12,
      color = "black",
      family = "sans"
    ),
    axis.title.y = element_text(
      size = 15,
      margin = margin(r = 10),
      color = "black",
      family = "sans"
    ),
    axis.line = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.ticks.margin = unit(0.1, "cm"),
    axis.ticks = element_line(color = "black")
  ) +
  labs(y = "SV number", fill = "SV type", x = "")

# Display the plot
print(bar_plot)

# Export as a PDF using cairo_pdf for better font compatibility
ggsave(
  filename = "SV_barplot.pdf",
  plot = bar_plot,
  device = cairo_pdf,
  width = 17,
  height = 8,
  dpi = 300,
  bg = "transparent"
)
```
