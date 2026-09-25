
library(tidyverse)
library(patchwork)

# =============================================================================
# 1. Read validation results
# =============================================================================

replication_df <- read.csv(
  "Ethnicity_SV_Discovery_Validation_Beta_STRICT.csv",
  stringsAsFactors = FALSE
)

replication_df <- replication_df %>%
  mutate(
    Same_Direction = sign(Discovery_Beta) == sign(Validation_Beta),
    SV_Type = factor(SV_Type, levels = c("dSV", "vSV"))
  )

dsv_df <- replication_df %>% filter(SV_Type == "dSV")
vsv_df <- replication_df %>% filter(SV_Type == "vSV")

# =============================================================================
# 2. Statistics
# =============================================================================

calc_stats <- function(dat) {
  n <- nrow(dat)
  concordance <- mean(dat$Same_Direction, na.rm = TRUE) * 100
  
  if (
    n < 3 ||
    length(unique(dat$Discovery_Beta)) < 2 ||
    length(unique(dat$Validation_Beta)) < 2
  ) {
    return(
      tibble(
        N = n,
        Spearman_Rho = NA_real_,
        P_value = NA_real_,
        Directional_Concordance = concordance
      )
    )
  }
  
  ct <- suppressWarnings(
    cor.test(
      dat$Discovery_Beta,
      dat$Validation_Beta,
      method = "spearman",
      exact = FALSE
    )
  )
  
  tibble(
    N = n,
    Spearman_Rho = as.numeric(ct$estimate),
    P_value = ct$p.value,
    Directional_Concordance = concordance
  )
}

dsv_stats <- calc_stats(dsv_df)
vsv_stats <- calc_stats(vsv_df)

print(dsv_stats)
print(vsv_stats)

# =============================================================================
# 3. Overall Discovery vs Validation correlation plot
# =============================================================================

make_overall_plot <- function(dat, sv_type, point_color) {
  st <- calc_stats(dat)
  
  label <- paste0(
    "Spearman \u03c1 = ",
    ifelse(is.na(st$Spearman_Rho), "NA", sprintf("%.3f", st$Spearman_Rho)),
    "\nP = ",
    ifelse(
      is.na(st$P_value),
      "NA",
      format.pval(st$P_value, digits = 3, eps = 1e-4)
    ),
    "\nDirectional concordance = ",
    sprintf("%.1f", st$Directional_Concordance),
    "%",
    "\nn = ",
    st$N
  )
  
  xr <- range(dat$Discovery_Beta, na.rm = TRUE)
  yr <- range(dat$Validation_Beta, na.rm = TRUE)
  
  ggplot(
    dat,
    aes(
      x = Discovery_Beta,
      y = Validation_Beta
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey65",
      linewidth = 0.45
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey65",
      linewidth = 0.45
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = "grey45",
      fill = "grey88",
      linewidth = 0.7
    ) +
    geom_point(
      size = 2.8,
      alpha = 0.78,
      color = point_color
    ) +
    annotate(
      "text",
      x = xr[1] + 0.03 * diff(xr),
      y = yr[2] - 0.03 * diff(yr),
      label = label,
      hjust = 0,
      vjust = 1,
      size = 4
    ) +
    labs(
      x = paste0("Discovery ", sv_type, " effect size (\u03b2)"),
      y = paste0("Validation ", sv_type, " effect size (\u03b2)")
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(size = 13)
    )
}

# =============================================================================
# 4. Overall dSV correlation plot
# =============================================================================

p_dsv <- make_overall_plot(
  dsv_df,
  "dSV",
  "#F4ACB7"
)

print(p_dsv)

ggsave(
  "dSV_Overall_Correlation.pdf",
  p_dsv,
  device = cairo_pdf,
  width = 6.5,
  height = 6,
  dpi = 300
)

# =============================================================================
# 5. Overall vSV correlation plot
# =============================================================================

p_vsv <- make_overall_plot(
  vsv_df,
  "vSV",
  "#8ECAE6"
)

print(p_vsv)

ggsave(
  "vSV_Overall_Correlation.pdf",
  p_vsv,
  device = cairo_pdf,
  width = 6.5,
  height = 6,
  dpi = 300
)

# =============================================================================
# 6. Overall directional concordance
# =============================================================================

concordance_df <- bind_rows(
  calc_stats(dsv_df) %>% mutate(SV_Type = "dSV"),
  calc_stats(vsv_df) %>% mutate(SV_Type = "vSV")
) %>%
  mutate(
    SV_Type = factor(SV_Type, levels = c("dSV", "vSV"))
  )

print(concordance_df)

write.csv(
  concordance_df,
  "Overall_Directional_Concordance.csv",
  row.names = FALSE
)

p_concordance <- ggplot(
  concordance_df,
  aes(
    x = SV_Type,
    y = Directional_Concordance,
    fill = SV_Type
  )
) +
  geom_col(
    width = 0.6,
    color = "black",
    linewidth = 0.35
  ) +
  geom_hline(
    yintercept = 50,
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.5
  ) +
  geom_text(
    aes(
      label = paste0(
        sprintf("%.1f", Directional_Concordance),
        "%\n(n=",
        N,
        ")"
      )
    ),
    vjust = -0.3,
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "dSV" = "#F4ACB7",
      "vSV" = "#8ECAE6"
    ),
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 110),
    breaks = seq(0, 100, 20),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Directional concordance (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(
      color = "black",
      size = 12
    ),
    axis.title.y = element_text(
      size = 13
    )
  )

print(p_concordance)

ggsave(
  "Overall_Directional_Concordance.pdf",
  p_concordance,
  device = cairo_pdf,
  width = 5,
  height = 5,
  dpi = 300
)

# =============================================================================
# 7. Optional combined figure
# =============================================================================

combined_plot <- p_dsv + p_vsv + p_concordance +
  plot_layout(
    widths = c(1, 1, 0.8)
  )

print(combined_plot)

ggsave(
  "Overall_Replication_Combined.pdf",
  combined_plot,
  device = cairo_pdf,
  width = 16,
  height = 5.5,
  dpi = 300
)
