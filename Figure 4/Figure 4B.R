
# ============================================================================
# Diet PERMANOVA R2 Bar Plot
# Significant dietary factors are colored by classification;
# non-significant factors are shown in gray
# ============================================================================

library(tidyverse)
library(stringr)
library(ggplot2)

# ============================================================================
# 1. Read data
# ============================================================================

df <- read.csv(
  "die_PERMANOVA.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

classification <- read.csv(
  "classification.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("PERMANOVA result columns:\n")
print(colnames(df))

cat("\nClassification columns:\n")
print(colnames(classification))

# ============================================================================
# 2. Check required columns
# ============================================================================

required_cols <- c(
  "Factor",
  "R2",
  "Combined_P_value",
  "Combined_p_adj_BH"
)

missing_cols <- setdiff(
  required_cols,
  colnames(df)
)

if (length(missing_cols) > 0) {
  stop(
    "The PERMANOVA results are missing the following columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

if (!all(c("Factor", "Group") %in% colnames(classification))) {
  stop(
    "classification.csv must contain the columns Factor and Group. Current columns: ",
    paste(colnames(classification), collapse = ", ")
  )
}

# ============================================================================
# 3. Name standardization functions
# Handle differences caused by underscores, extra spaces, and capitalization
# ============================================================================

clean_key <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "_", " ")
  x <- str_squish(x)
  x <- str_to_lower(x)
  x
}

clean_label <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "_", " ")
  x <- str_squish(x)
  x
}

# ============================================================================
# 4. Prepare classification data
# ============================================================================

classification_clean <- classification %>%
  transmute(
    Match_Key = clean_key(Factor),
    Diet_Group = as.character(Group)
  ) %>%
  filter(
    !is.na(Match_Key),
    Match_Key != "",
    !is.na(Diet_Group),
    Diet_Group != ""
  ) %>%
  distinct(
    Match_Key,
    .keep_all = TRUE
  )

cat("\nClassification groups:\n")
print(
  classification_clean %>%
    count(Diet_Group)
)

# ============================================================================
# 5. Merge PERMANOVA results with classification
# ============================================================================

plot_data <- df %>%
  transmute(
    Factor_original = as.character(Factor),
    Match_Key = clean_key(Factor),
    Factor = clean_label(Factor),
    R2 = as.numeric(R2),
    Combined_P_value = as.numeric(Combined_P_value),
    p_adj = as.numeric(Combined_p_adj_BH)
  ) %>%
  left_join(
    classification_clean,
    by = "Match_Key"
  )

# ============================================================================
# 6. Check factors without a matched classification
# ============================================================================

unmatched <- plot_data %>%
  filter(is.na(Diet_Group)) %>%
  select(Factor_original)

if (nrow(unmatched) > 0) {
  cat("\nThe following factors were not matched to classification.csv:\n")
  print(unmatched)
} else {
  cat("\nAll factors were successfully matched to a classification.\n")
}

# ============================================================================
# 7. Define colors
# ============================================================================

group_colors_raw <- c(
  "#E2B0BA",
  "#87AFD5",
  "#6EB9C3",
  "#C98B88",
  "#93C89A",
  "#FFCC98",
  "#E1D1BA",
  "#EDAFA9",
  "#0AB0C8",
  "#9781BB",
  "#E8BD65",
  "#E39844",
  "#ADD1E5"
)

diet_groups <- sort(
  unique(
    na.omit(
      plot_data$Diet_Group
    )
  )
)

if (length(diet_groups) > length(group_colors_raw)) {
  stop("The number of diet classifications exceeds the number of predefined colors.")
}

group_colors <- setNames(
  group_colors_raw[
    seq_along(diet_groups)
  ],
  diet_groups
)

non_sig_color <- "#CCCCCC"

# ============================================================================
# 8. Define significance and ordering
# ============================================================================

plot_data <- plot_data %>%
  mutate(
    Significant = ifelse(
      !is.na(p_adj) & p_adj < 0.05,
      "Significant",
      "Not significant"
    ),
    Fill_Group = ifelse(
      Significant == "Significant" & !is.na(Diet_Group),
      Diet_Group,
      "Not significant"
    ),
    Significance_Label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  arrange(
    desc(R2)
  ) %>%
  mutate(
    Factor = factor(
      Factor,
      levels = rev(unique(Factor))
    )
  )

cat(
  "\nNumber of significant factors:",
  sum(plot_data$p_adj < 0.05, na.rm = TRUE),
  "/", nrow(plot_data), "\n"
)

# ============================================================================
# 9. Export processed plotting data
# ============================================================================

write.csv(
  plot_data,
  "permanova_diet_R2_plot_data.csv",
  row.names = FALSE,
  quote = FALSE
)

# ============================================================================
# 10. Plot colors
# ============================================================================

plot_colors <- c(
  group_colors,
  "Not significant" = non_sig_color
)

# ============================================================================
# 11. Horizontal bar plot
# ============================================================================

p <- ggplot(
  plot_data,
  aes(
    x = Factor,
    y = R2,
    fill = Fill_Group
  )
) +
  geom_col(
    width = 0.72,
    color = NA
  ) +
  geom_text(
    aes(
      label = Significance_Label
    ),
    hjust = -0.25,
    size = 3.5,
    color = "black"
  ) +
  coord_flip(
    clip = "off"
  ) +
  scale_fill_manual(
    values = plot_colors,
    name = "Diet category"
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.12)
    )
  ) +
  labs(
    x = NULL,
    y = expression("marginal " * R^2)
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    axis.text.y = element_text(
      size = 8.5,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 10,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 12,
      color = "black"
    ),
    axis.line = element_line(
      color = "black",
      linewidth = 0.5
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.4
    ),
    legend.position = "right",
    legend.title = element_text(
      size = 11
    ),
    legend.text = element_text(
      size = 10
    ),
    plot.margin = margin(
      t = 10,
      r = 25,
      b = 10,
      l = 10
    )
  )

print(p)

# ============================================================================
# 12. Export PDF
# ============================================================================

ggsave(
  "permanova_diet_R2_barplot.pdf",
  p,
  width = 8,
  height = max(
    7,
    nrow(plot_data) * 0.20
  ),
  device = cairo_pdf,
  bg = "white"
)

# ============================================================================
# 13. Export PNG
# ============================================================================

ggsave(
  "permanova_diet_R2_barplot.png",
  p,
  width = 8,
  height = max(
    7,
    nrow(plot_data) * 0.20
  ),
  dpi = 600,
  bg = "white"
)

cat("\nPlotting completed.\n")
cat("PDF: permanova_diet_R2_barplot.pdf\n")
cat("PNG: permanova_diet_R2_barplot.png\n")
cat("Plot data: permanova_diet_R2_plot_data.csv\n")

