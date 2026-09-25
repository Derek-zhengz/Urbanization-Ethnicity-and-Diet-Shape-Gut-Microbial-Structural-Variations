library(ggplot2)
library(dplyr)

#==================================================
# 1. Color palette
#==================================================

mycol <- c(
  "#E2B0BA", "#87AFD5", "#6EB9C3", "#C98B88",
  "#93C89A", "#FFCC98", "#E1D1BA", "#EDAFA9",
  "#0AB0C8", "#9781BB", "#E8BD65", "#E39844",
  "#ADD1E5"
)

#==================================================
# 2. Read PERMANOVA results
#==================================================

data <- read.csv(
  "PERMANOVA.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Keep only the required columns
data <- data %>%
  select(
    SV_Type,
    Group,
    Marginal_R2,
    p_adj
  )

#==================================================
# 3. Calculate mean marginal R2 across dSV and vSV
#==================================================

group_order <- c(
  "General",
  "Ethnicity",
  "Diet",
  "Medication",
  "Urbanization",
  "BMI",
  "AGE",
  "SEX"
)

group_colors <- c(
  "General"      = mycol[1],
  "Ethnicity"    = mycol[2],
  "Diet"         = mycol[3],
  "Medication"   = mycol[4],
  "Urbanization" = mycol[5],
  "BMI"          = mycol[6],
  "AGE"          = mycol[7],
  "SEX"          = mycol[8]
)

data_mean <- data %>%
  filter(Group %in% group_order) %>%
  group_by(Group) %>%
  summarise(
    Mean_R2 = mean(Marginal_R2, na.rm = TRUE),
    # Add * only when both dSV and vSV are significant
    significant = ifelse(
      n_distinct(SV_Type) == 2 &&
        all(p_adj < 0.05, na.rm = TRUE),
      "*",
      ""
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Group = factor(Group, levels = group_order),
    Mean_R2_percent = Mean_R2 * 100
  )

print(data_mean)

#==================================================
# 4. Plot
#==================================================

p <- ggplot(
  data_mean,
  aes(
    x = Group,
    y = Mean_R2_percent,
    fill = Group
  )
) +
  geom_col(
    width = 0.7,
    color = "black",
    linewidth = 0.35
  ) +
  geom_text(
    aes(
      y = Mean_R2_percent + 0.12,
      label = significant
    ),
    size = 5,
    color = "black",
    vjust = 0
  ) +
  scale_fill_manual(
    values = group_colors,
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 5.2),
    breaks = seq(0, 5, 1),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = expression("Mean marginal " * R^2 * " (%)")
  ) +
  theme_classic() +
  theme(
    text = element_text(
      family = "sans",
      color = "black"
    ),
    axis.text.x = element_text(
      size = 11,
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 11,
      color = "black"
    ),
    axis.title.y = element_text(
      size = 14
    ),
    axis.line = element_line(
      color = "black",
      linewidth = 0.5
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    plot.margin = margin(
      t = 15,
      r = 15,
      b = 10,
      l = 15
    )
  )

print(p)

#==================================================
# 5. Save as PDF
#==================================================

ggsave(
  "Validation_cohort_mean_marginal_R2_barplot.pdf",
  plot = p,
  device = cairo_pdf,
  width = 7.5,
  height = 5,
  dpi = 300,
  bg = "transparent"
)