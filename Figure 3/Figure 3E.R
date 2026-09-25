
# ============================================================================
# Yunnan Cohort PCoA, Envfit, and Combined Visualization
# ============================================================================

library(vegan)
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ---------- 1. Read and merge metadata ----------
cat("========== 1. Read and merge metadata ==========\n")

metadata_main <- read_excel("metadata.xlsx")
colnames(metadata_main)[1] <- "Sample_ID"

urbanization <- read_excel("Urbanization.xlsx")
colnames(urbanization)[1] <- "Sample_ID"

medication <- read_excel("Medication.xlsx")
colnames(medication)[1] <- "Sample_ID"

general <- read_excel("Gerneral metadata.xlsx")
colnames(general)[1] <- "Sample_ID"

diet <- read_excel("Diet.xlsx")
colnames(diet)[1] <- "Sample_ID"

# Diet prevalence filtering
cat("\n--- Diet filtering (ethnicity-specific prevalence >= 1%) ---\n")

eth_col <- grep(
  "ethnic|民族",
  colnames(metadata_main),
  value = TRUE,
  ignore.case = TRUE
)[1]

if (is.na(eth_col)) stop("Ethnicity column was not found.")

diet_with_ethnic <- diet %>%
  inner_join(
    metadata_main %>% select(Sample_ID, !!sym(eth_col)),
    by = "Sample_ID"
  ) %>%
  rename(ethnicity = !!sym(eth_col))

clean_diet <- function(data, th = 1) {
  non <- c("Sample_ID", "ethnicity")
  
  ep <- data %>%
    group_by(ethnicity) %>%
    summarise(
      across(
        -any_of(non),
        ~ sum(.x > 0 & !is.na(.x)) / n() * 100
      )
    ) %>%
    pivot_longer(
      -ethnicity,
      names_to = "item",
      values_to = "prev"
    )
  
  keep <- ep %>%
    group_by(item) %>%
    summarise(keep = any(prev >= th)) %>%
    filter(keep) %>%
    pull(item)
  
  data %>%
    select(Sample_ID, all_of(keep))
}

diet_filtered <- clean_diet(diet_with_ethnic, 1)

# Add prefixes
colnames(urbanization)[-1] <- paste0("Urban_", colnames(urbanization)[-1])
colnames(medication)[-1] <- paste0("Med_", colnames(medication)[-1])
colnames(general)[-1] <- paste0("Gen_", colnames(general)[-1])
colnames(diet_filtered)[-1] <- paste0("Diet_", colnames(diet_filtered)[-1])

metadata_combined <- metadata_main %>%
  inner_join(urbanization, by = "Sample_ID") %>%
  inner_join(medication, by = "Sample_ID") %>%
  inner_join(general, by = "Sample_ID") %>%
  inner_join(diet_filtered, by = "Sample_ID") %>%
  rename(sample_id = Sample_ID) %>%
  mutate(sample_id = trimws(as.character(sample_id)))

cat("Total samples after metadata merging:", nrow(metadata_combined), "\n")

# ---------- 2. Build group.csv mapping ----------
base_vars <- intersect(
  c("Ethnicity", "AGE", "SEX", "BMI"),
  colnames(metadata_combined)
)

urban_cols <- grep("^Urban_", colnames(metadata_combined), value = TRUE)
med_cols <- grep("^Med_", colnames(metadata_combined), value = TRUE)
gen_cols <- grep("^Gen_", colnames(metadata_combined), value = TRUE)
diet_cols <- grep("^Diet_", colnames(metadata_combined), value = TRUE)

group_list <- list(
  Ethnicity = base_vars[base_vars == "Ethnicity"],
  AGE = base_vars[base_vars == "AGE"],
  SEX = base_vars[base_vars == "SEX"],
  BMI = base_vars[base_vars == "BMI"],
  Urbanization = urban_cols,
  Medication = med_cols,
  General = gen_cols,
  Diet = diet_cols
)

group_info <- data.frame(
  Factor = unlist(group_list, use.names = FALSE),
  Group = rep(names(group_list), times = sapply(group_list, length)),
  stringsAsFactors = FALSE
)

write_csv(group_info, "group.csv")

# ---------- 3. Calculate species-abundance Bray-Curtis PCoA covariates ----------
cat("\n========== 3. Species-abundance Bray-Curtis PCoA ==========\n")

metaphlan_raw <- read_tsv(
  "yunnan_7_years_later_combined_metaphlan_profile.tsv",
  comment = ""
)

if (grepl("#", colnames(metaphlan_raw)[1])) {
  colnames(metaphlan_raw)[1] <- "clade_name"
}

mat <- as.data.frame(metaphlan_raw)
rownames(mat) <- mat$clade_name
mat <- mat[, -1]

tmat <- as.matrix(t(mat))
rownames(tmat) <- trimws(as.character(rownames(tmat)))

common_samples <- intersect(
  metadata_combined$sample_id,
  rownames(tmat)
)

tmat_sub <- tmat[common_samples, , drop = FALSE]
tmat_sub <- tmat_sub[
  rowSums(tmat_sub) > 0,
  colSums(tmat_sub) > 0
]

dist_abun <- vegdist(tmat_sub, method = "bray")
pcoa_abun <- cmdscale(dist_abun, k = 3, eig = TRUE)

pcoa_df <- data.frame(
  sample_id = rownames(tmat_sub),
  Abun_PCoA1 = pcoa_abun$points[, 1],
  Abun_PCoA2 = pcoa_abun$points[, 2],
  Abun_PCoA3 = pcoa_abun$points[, 3],
  stringsAsFactors = FALSE
)

# ---------- 4. Metadata cleaning and data type conversion ----------
all_factors <- intersect(
  unique(group_info$Factor),
  colnames(metadata_combined)
)

metadata_complete <- metadata_combined %>%
  inner_join(pcoa_df, by = "sample_id")

numeric_targets <- c(
  "AGE",
  "BMI",
  diet_cols
)

for (nv in intersect(numeric_targets, colnames(metadata_complete))) {
  metadata_complete[[nv]] <- suppressWarnings(
    as.numeric(as.character(metadata_complete[[nv]]))
  )
}

valid_factors <- c()

for (var in all_factors) {
  vec <- metadata_complete[[var]]
  
  if (is.numeric(vec)) {
    metadata_complete[[var]][is.na(vec)] <- median(vec, na.rm = TRUE)
  } else {
    vec <- as.character(vec)
    vec[is.na(vec)] <- "Unknown"
    metadata_complete[[var]] <- as.factor(vec)
  }
  
  if (
    length(
      unique(
        metadata_complete[[var]][!is.na(metadata_complete[[var]])]
      )
    ) > 1
  ) {
    valid_factors <- c(valid_factors, var)
  }
}

valid_ids <- metadata_complete$sample_id

# ---------- 5. PCoA and Envfit analysis function ----------
run_pcoa_analysis <- function(
    dist_matrix,
    data_df,
    valid_factors,
    group_info,
    sv_type = "dSV"
) {
  cat(
    sprintf(
      "\n>>> Running PCoA analysis [%s] <<<\n",
      sv_type
    )
  )
  
  pcoa_fit <- cmdscale(
    dist_matrix,
    k = 5,
    eig = TRUE
  )
  
  eig_vals <- pcoa_fit$eig
  var_explained <- eig_vals / sum(eig_vals[eig_vals > 0]) * 100
  
  pcoa_points <- as.data.frame(
    pcoa_fit$points[, 1:3]
  )
  
  colnames(pcoa_points) <- paste0(
    sv_type,
    "_PCoA",
    1:3
  )
  
  pcoa_points$sample_id <- rownames(pcoa_points)
  
  write_csv(
    pcoa_points,
    sprintf("%s_pcoa_coords.csv", sv_type)
  )
  
  env_df <- data_df %>%
    filter(sample_id %in% pcoa_points$sample_id) %>%
    select(all_of(valid_factors))
  
  fit_env <- envfit(
    pcoa_fit$points[, 1:2],
    env_df,
    permutations = 999,
    na.rm = TRUE
  )
  
  res_vectors <- data.frame()
  
  if (!is.null(fit_env$vectors)) {
    res_vectors <- data.frame(
      Factor = names(fit_env$vectors$r),
      Type = "Numeric_Vector",
      r2 = fit_env$vectors$r,
      p_value = fit_env$vectors$pvals,
      stringsAsFactors = FALSE
    )
  }
  
  res_factors <- data.frame()
  
  if (!is.null(fit_env$factors)) {
    res_factors <- data.frame(
      Factor = names(fit_env$factors$r),
      Type = "Categorical_Factor",
      r2 = fit_env$factors$r,
      p_value = fit_env$factors$pvals,
      stringsAsFactors = FALSE
    )
  }
  
  fit_combined <- bind_rows(
    res_vectors,
    res_factors
  ) %>%
    left_join(group_info, by = "Factor") %>%
    mutate(
      SV_Type = sv_type,
      p_adj = p.adjust(p_value, method = "BH")
    ) %>%
    arrange(desc(r2))
  
  if (!is.null(fit_env$factors)) {
    centroids_df <- as.data.frame(
      fit_env$factors$centroids
    )
    
    centroids_df$Factor_Level <- rownames(centroids_df)
    
    write_csv(
      centroids_df,
      sprintf(
        "%s_pcoa_factor_centroids.csv",
        sv_type
      )
    )
  }
  
  return(
    list(
      pcoa_points = pcoa_points,
      envfit_res = fit_combined,
      var_explained = var_explained[1:3]
    )
  )
}

# ---------- 6. Run dSV analysis ----------
cat("\n========== 6. Run dSV analysis (Jaccard PCoA) ==========\n")

dsv <- read_csv("filtered_dsv.csv")
colnames(dsv)[1] <- "sample_id"

dsv_cleaned <- dsv %>%
  mutate(
    across(
      -sample_id,
      ~ replace_na(as.numeric(.), 0)
    )
  )

dsv_filtered <- dsv_cleaned %>%
  filter(sample_id %in% valid_ids)

dsv_mat <- as.matrix(
  dsv_filtered %>%
    select(-sample_id)
)

rownames(dsv_mat) <- dsv_filtered$sample_id

dsv_mat_valid <- dsv_mat[
  rowSums(dsv_mat != 0) > 0,
  colSums(dsv_mat != 0) > 0,
  drop = FALSE
]

metadata_dsv <- metadata_complete %>%
  filter(sample_id %in% rownames(dsv_mat_valid))

dist_dsv <- vegdist(
  dsv_mat_valid,
  method = "jaccard",
  binary = TRUE
)

dsv_pcoa_results <- run_pcoa_analysis(
  dist_dsv,
  metadata_dsv,
  valid_factors,
  group_info,
  sv_type = "dSV"
)

# ---------- 7. Run vSV analysis ----------
cat("\n========== 7. Run vSV analysis (Canberra PCoA) ==========\n")

vsv <- read_csv("filtered_vsv.csv")
colnames(vsv)[1] <- "sample_id"

min_max_scale <- function(x) {
  if (all(is.na(x))) return(x)
  
  non_na <- x[!is.na(x)]
  min_val <- min(non_na)
  max_val <- max(non_na)
  
  if (min_val == max_val) {
    x[!is.na(x)] <- 0.5
  } else {
    x[!is.na(x)] <- (non_na - min_val) / (max_val - min_val)
  }
  
  x[is.na(x)] <- 0
  return(x)
}

vsv_processed <- bind_cols(
  sample_id = vsv$sample_id,
  vsv %>%
    select(-sample_id) %>%
    mutate(
      across(
        everything(),
        min_max_scale
      )
    )
)

vsv_filtered <- vsv_processed %>%
  filter(sample_id %in% valid_ids)

vsv_mat <- as.matrix(
  vsv_filtered %>%
    select(-sample_id)
)

rownames(vsv_mat) <- vsv_filtered$sample_id

vsv_mat_valid <- vsv_mat[
  rowSums(vsv_mat != 0) > 0,
  colSums(vsv_mat != 0) > 0,
  drop = FALSE
]

metadata_vsv <- metadata_complete %>%
  filter(sample_id %in% rownames(vsv_mat_valid))

dist_vsv <- vegdist(
  vsv_mat_valid,
  method = "canberra"
)

vsv_pcoa_results <- run_pcoa_analysis(
  dist_vsv,
  metadata_vsv,
  valid_factors,
  group_info,
  sv_type = "vSV"
)

bind_rows(
  dsv_pcoa_results$envfit_res,
  vsv_pcoa_results$envfit_res
) %>%
  write_csv(
    "pcoa_envfit_combined_results.csv"
  )

# ---------- 8. Combined PCoA visualization ----------
cat("\n========== 8. Generate combined PCoA visualization ==========\n")

# Match samples and calculate combined distance matrix
common_plot_ids <- Reduce(
  intersect,
  list(
    rownames(vsv_mat_valid),
    rownames(dsv_mat_valid),
    metadata_complete$sample_id
  )
)

d_vsv_sub <- vegdist(
  vsv_mat_valid[common_plot_ids, ],
  method = "canberra"
)

d_dsv_sub <- vegdist(
  dsv_mat_valid[common_plot_ids, ],
  method = "jaccard",
  binary = TRUE
)

d_comb <- (d_vsv_sub + d_dsv_sub) / 2

# PCoA
pcoa_cap <- capscale(
  d_comb ~ 1,
  add = TRUE
)

pcoa_pts <- scores(
  pcoa_cap,
  display = "sites"
)[, 1:2]

eig <- pcoa_cap$CA$eig
pc1 <- round(
  100 * eig[1] / sum(eig),
  2
)

pc2 <- round(
  100 * eig[2] / sum(eig),
  2
)

# Prepare metadata for visualization
eth_var_name <- eth_col

plotdata <- data.frame(
  sample = rownames(pcoa_pts),
  dim1 = pcoa_pts[, 1],
  dim2 = pcoa_pts[, 2],
  stringsAsFactors = FALSE
) %>%
  left_join(
    metadata_complete %>%
      select(
        sample_id,
        group = !!sym(eth_var_name)
      ),
    by = c(
      "sample" = "sample_id"
    )
  ) %>%
  mutate(
    group = as.character(group)
  )

# PERMANOVA
set.seed(123)

adonis_res <- adonis2(
  d_comb ~ group,
  data = plotdata,
  permutations = 999
)

R2_val <- adonis_res$R2[1]
pval_val <- adonis_res$`Pr(>F)`[1]

stat_text <- sprintf(
  "PERMANOVA:\nR² = %.4f\np = %.4f",
  R2_val,
  pval_val
)

# Colors and point shapes
unique_groups <- unique(plotdata$group)
length_group <- length(unique_groups)

mycol <- c(
  "#E2B0BA", "#87AFD5", "#6EB9C3", "#C98B88",
  "#93C89A", "#FFCC98", "#E1D1BA", "#EDAFA9",
  "#0AB0C8", "#9781BB", "#E8BD65", "#E39844",
  "#ADD1E5"
)

times_col <- length_group %/% length(mycol)
res_col <- length_group %% length(mycol)

col <- c(
  rep(mycol, times_col),
  mycol[1:res_col]
)

pich <- rep(
  16,
  length_group
)

names(pich) <- unique_groups

legend_ncol <- if (
  length_group > 8
) {
  2
} else {
  1
}

# Main PCoA plot
pcoa_plot <- ggplot(
  plotdata,
  aes(
    x = dim1,
    y = dim2
  )
) +
  geom_point(
    aes(
      colour = group,
      shape = group
    ),
    size = 3
  ) +
  stat_ellipse(
    aes(
      x = dim1,
      y = dim2,
      color = group
    ),
    level = 0.95
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted"
  ) +
  xlab(
    paste0(
      "PCoA1 (",
      pc1,
      "%)"
    )
  ) +
  ylab(
    paste0(
      "PCoA2 (",
      pc2,
      "%)"
    )
  ) +
  scale_colour_manual(
    values = col
  ) +
  scale_shape_manual(
    values = pich
  ) +
  theme_classic(
    base_size = 14
  ) +
  theme(
    axis.line = element_blank(),
    axis.line.x.bottom = element_line(
      colour = "black"
    ),
    axis.line.y.left = element_line(
      colour = "black"
    ),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.text.y.right = element_blank(),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.grid = element_blank(),
    axis.title = element_text(
      color = "black",
      size = 18
    ),
    axis.text = element_text(
      colour = "black",
      size = 16,
      margin = margin(
        0.6,
        0.6,
        0.6,
        0.6,
        "lines"
      )
    ),
    legend.title = element_blank(),
    legend.text = element_text(
      size = 10
    ),
    legend.key = element_blank(),
    legend.position = c(
      0.85,
      0.8
    ),
    legend.background = element_rect(
      fill = "transparent",
      colour = NA
    )
  ) +
  guides(
    col = guide_legend(
      ncol = legend_ncol
    ),
    shape = guide_legend(
      ncol = legend_ncol
    )
  )

# Top boxplot
box_top <- ggplot(
  plotdata,
  aes(
    x = group,
    y = dim1,
    fill = group
  )
) +
  geom_boxplot(
    show.legend = FALSE,
    outlier.colour = "gray50",
    outlier.size = 1,
    outlier.alpha = 0.5
  ) +
  stat_boxplot(
    geom = "errorbar",
    width = 0.1,
    size = 0.1
  ) +
  geom_jitter(
    show.legend = FALSE,
    color = "gray50",
    alpha = 0.5,
    size = 1
  ) +
  scale_fill_manual(
    values = col
  ) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.line = element_line(
      colour = "black"
    ),
    axis.ticks = element_line(
      color = "black"
    ),
    axis.text.x = element_blank(),
    axis.text.y = element_text(
      colour = "black",
      size = 12
    )
  )

# Right-side boxplot
box_right <- ggplot(
  plotdata,
  aes(
    x = group,
    y = dim2,
    fill = group
  )
) +
  geom_boxplot(
    show.legend = FALSE,
    outlier.colour = "gray50",
    outlier.size = 1,
    outlier.alpha = 0.5
  ) +
  stat_boxplot(
    geom = "errorbar",
    width = 0.1,
    size = 0.1
  ) +
  geom_jitter(
    show.legend = FALSE,
    color = "gray50",
    alpha = 0.5,
    size = 1,
    width = 0.1,
    height = 0
  ) +
  scale_fill_manual(
    values = col
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.line = element_line(
      colour = "black"
    ),
    axis.ticks = element_line(
      color = "black"
    ),
    axis.text.x = element_text(
      colour = "black",
      size = 12,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_blank()
  )

# PERMANOVA statistics panel
stat_box <- ggplot(
  plotdata,
  aes(
    dim1,
    dim2
  )
) +
  annotate(
    "text",
    x = mean(
      range(
        plotdata$dim1
      )
    ),
    y = mean(
      range(
        plotdata$dim2
      )
    ),
    label = stat_text,
    size = 4,
    hjust = 0.5,
    vjust = 0.5
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )

# Combine plots and export
final_pcoa_plot <- box_top + stat_box +
  pcoa_plot + box_right +
  plot_layout(
    heights = c(1, 4),
    widths = c(4, 1),
    ncol = 2,
    nrow = 2
  )

pdf(
  "PCoA_result.pdf",
  width = 9,
  height = 8
)

print(final_pcoa_plot)

dev.off()

cat(
  "\nPCoA combined figure exported to PCoA_result.pdf.\n"
)
```
