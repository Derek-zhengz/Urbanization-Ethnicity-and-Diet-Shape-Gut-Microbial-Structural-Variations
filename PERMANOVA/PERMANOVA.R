# =============================================================================
# PERMANOVA analysis of dSV and vSV
# =============================================================================

library(vegan)
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

set.seed(12345)

# =============================================================================
# 1. Input files
# =============================================================================

metadata_file <- "metadata_clean.xlsx"
dsv_file <- "final_filtered_dsv.csv"
vsv_file <- "final_filtered_vsv.csv"
group_file <- "group.csv"
classification_file <- "classification.csv"
species_map_file <- "sv_to_metaphlan_final_mapping.csv"
metaphlan_file <- "yunnan_7_years_ago_combined_metaphlan_profile.tsv"

# =============================================================================
# 2. Load data once
# =============================================================================

cat("========== Load input data ==========\n")

metadata <- as.data.frame(read_excel(metadata_file))
dsv <- read.csv(dsv_file, check.names = FALSE, stringsAsFactors = FALSE)
vsv <- read.csv(vsv_file, check.names = FALSE, stringsAsFactors = FALSE)
group_info <- read_csv(group_file, show_col_types = FALSE)
classification <- read_csv(classification_file, show_col_types = FALSE)
species_map <- read.csv(species_map_file, check.names = FALSE, stringsAsFactors = FALSE)
metaphlan_raw <- read_tsv(metaphlan_file, comment = "", show_col_types = FALSE)

if (grepl("#", colnames(metaphlan_raw)[1])) {
  colnames(metaphlan_raw)[1] <- "clade_name"
}

# =============================================================================
# 3. Helper functions
# =============================================================================

quote_var <- function(x) {
  paste0("`", gsub("`", "", x), "`")
}

make_formula <- function(response, variables) {
  as.formula(
    paste(response, "~", paste(quote_var(variables), collapse = " + "))
  )
}

normalize_name <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- str_squish(x)
  tolower(x)
}

min_max_scale <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  
  if (all(is.na(x))) {
    return(rep(0, length(x)))
  }
  
  non_na <- x[!is.na(x)]
  min_val <- min(non_na)
  max_val <- max(non_na)
  
  if (min_val == max_val) {
    x[!is.na(x)] <- 0.5
  } else {
    x[!is.na(x)] <- (x[!is.na(x)] - min_val) / (max_val - min_val)
  }
  
  x[is.na(x)] <- 0
  x
}

combine_two_p <- function(p1, p2) {
  p <- c(p1, p2)
  p <- p[!is.na(p)]
  
  if (length(p) == 0) {
    return(NA_real_)
  }
  
  p[p <= 0] <- 1 / 1000
  fisher_stat <- -2 * sum(log(p))
  
  pchisq(
    fisher_stat,
    df = 2 * length(p),
    lower.tail = FALSE
  )
}

get_model_R2 <- function(dist_mat, metadata_df, variables) {
  fit <- adonis2(
    make_formula("dist_mat", variables),
    data = metadata_df,
    permutations = 0
  )
  
  if ("Model" %in% rownames(fit)) {
    return(as.numeric(fit["Model", "R2"]))
  }
  
  sum(
    fit[!rownames(fit) %in% c("Residual", "Total"), "R2"],
    na.rm = TRUE
  )
}

# Calculate marginal R2 and p-value for one factor.
# The full/reduced model logic is unchanged.
calc_factor_marginal <- function(dist_mat, metadata_df, all_factors, target_factor) {
  base_covars <- c("Abun_PCoA1", "Abun_PCoA2", "Abun_PCoA3")
  
  fit_full <- adonis2(
    make_formula("dist_mat", c(base_covars, all_factors)),
    data = metadata_df,
    permutations = 999
  )
  
  remaining_factors <- setdiff(all_factors, target_factor)
  
  fit_reduced <- adonis2(
    make_formula("dist_mat", c(base_covars, remaining_factors)),
    data = metadata_df,
    permutations = 999
  )
  
  if ("Model" %in% rownames(fit_full)) {
    total_R2_full <- as.numeric(fit_full["Model", "R2"])
  } else {
    total_R2_full <- sum(
      fit_full[!rownames(fit_full) %in% c("Residual", "Total"), "R2"],
      na.rm = TRUE
    )
  }
  
  if ("Model" %in% rownames(fit_reduced)) {
    total_R2_reduced <- as.numeric(fit_reduced["Model", "R2"])
  } else {
    total_R2_reduced <- sum(
      fit_reduced[!rownames(fit_reduced) %in% c("Residual", "Total"), "R2"],
      na.rm = TRUE
    )
  }
  
  marginal_R2 <- max(0, total_R2_full - total_R2_reduced)
  
  fit_factor <- adonis2(
    make_formula("dist_mat", c(base_covars, target_factor)),
    data = metadata_df,
    permutations = 999,
    by = "margin"
  )
  
  rn_clean <- gsub("`", "", rownames(fit_factor))
  target_idx <- which(rn_clean == target_factor)
  
  p_val <- if (length(target_idx) > 0) {
    as.numeric(fit_factor[target_idx[1], "Pr(>F)"])
  } else {
    NA_real_
  }
  
  data.frame(
    Factor = target_factor,
    Marginal_R2 = marginal_R2,
    P_value = p_val,
    stringsAsFactors = FALSE
  )
}

run_factor_group <- function(
    group_name,
    dist_mat,
    metadata_df,
    all_factors,
    group_info,
    sv_type
) {
  factors <- group_info %>%
    filter(Group == group_name) %>%
    pull(Factor)
  
  factors <- factors[factors %in% all_factors]
  
  if (length(factors) == 0) {
    stop(paste0("No valid factors found for group: ", group_name))
  }
  
  cat("\n--- ", sv_type, " | ", group_name, " ---\n", sep = "")
  
  map_df(
    factors,
    ~ calc_factor_marginal(dist_mat, metadata_df, all_factors, .x)
  ) %>%
    mutate(
      SV_Type = sv_type,
      p_adj = p.adjust(P_value, method = "BH")
    ) %>%
    arrange(desc(Marginal_R2))
}

# Group-level PERMANOVA using the original leave-one-group-out logic.
run_permanova_group_analysis <- function(
    dist_matrix,
    data_df,
    valid_factors,
    group_info,
    sv_type
) {
  cat(sprintf("\n>>> Start group-level analysis [%s] <<<\n", sv_type))
  
  base_covars <- c("Abun_PCoA1", "Abun_PCoA2", "Abun_PCoA3")
  group_mapping <- group_info %>% filter(Factor %in% valid_factors)
  groups_list <- unique(group_mapping$Group)
  
  fit_base <- adonis2(
    make_formula("dist_matrix", base_covars),
    data = data_df,
    permutations = 999
  )
  
  base_species_R2 <- sum(fit_base[base_covars, "R2"], na.rm = TRUE)
  base_p_val <- fit_base$`Pr(>F)`[1]
  
  fit_full <- adonis2(
    make_formula("dist_matrix", c(base_covars, valid_factors)),
    data = data_df,
    permutations = 999
  )
  
  all_model_rows <- setdiff(rownames(fit_full), c("Residual", "Total"))
  total_model_terms_R2 <- sum(fit_full[all_model_rows, "R2"], na.rm = TRUE)
  full_p_val <- fit_full$`Pr(>F)`[1]
  
  all_factors_combined_pure_R2 <- max(
    0,
    total_model_terms_R2 - base_species_R2
  )
  
  group_results <- map_df(groups_list, function(grp) {
    grp_factors <- group_mapping %>%
      filter(Group == grp) %>%
      pull(Factor)
    
    remaining_factors <- setdiff(valid_factors, grp_factors)
    
    fit_reduced <- adonis2(
      make_formula("dist_matrix", c(base_covars, remaining_factors)),
      data = data_df,
      permutations = 999
    )
    
    reduced_rows <- setdiff(rownames(fit_reduced), c("Residual", "Total"))
    reduced_model_terms_R2 <- sum(
      fit_reduced[reduced_rows, "R2"],
      na.rm = TRUE
    )
    
    grp_marginal_R2 <- max(
      0,
      total_model_terms_R2 - reduced_model_terms_R2
    )
    
    fit_grp <- adonis2(
      make_formula("dist_matrix", c(base_covars, grp_factors)),
      data = data_df,
      permutations = 999,
      by = "margin"
    )
    
    env_rows <- setdiff(
      rownames(fit_grp),
      c(base_covars, "Residual", "Total")
    )
    
    p_vals <- fit_grp[env_rows, "Pr(>F)"]
    p_vals_clean <- p_vals[!is.na(p_vals) & is.finite(p_vals)]
    grp_p_val <- if (length(p_vals_clean) > 0) min(p_vals_clean) else 0.001
    
    data.frame(
      SV_Type = sv_type,
      Group = grp,
      Factor_Count = length(grp_factors),
      Group_Marginal_R2 = grp_marginal_R2,
      p_value = grp_p_val,
      stringsAsFactors = FALSE
    )
  }) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(desc(Group_Marginal_R2))
  
  overall_rows <- data.frame(
    SV_Type = sv_type,
    Group = c(
      "ALL_ENVIRONMENTAL_FACTORS_COMBINED",
      "SPECIES_ABUNDANCE_COVARIATES",
      "TOTAL_EXPLAINED_VARIANCE"
    ),
    Factor_Count = c(length(valid_factors), 3, length(valid_factors) + 3),
    Group_Marginal_R2 = c(
      all_factors_combined_pure_R2,
      base_species_R2,
      total_model_terms_R2
    ),
    p_value = c(full_p_val, base_p_val, full_p_val),
    p_adj = c(full_p_val, base_p_val, full_p_val),
    stringsAsFactors = FALSE
  )
  
  bind_rows(group_results, overall_rows)
}

# Classification-level block permutation PERMANOVA.
calc_classification_permanova <- function(
    dist_mat,
    metadata_df,
    all_factors,
    classification_factors,
    classification_name,
    sv_type,
    nperm = 999
) {
  base_covars <- c("Abun_PCoA1", "Abun_PCoA2", "Abun_PCoA3")
  
  classification_factors <- intersect(
    classification_factors,
    all_factors
  )
  
  if (length(classification_factors) == 0) {
    return(
      data.frame(
        Classification = classification_name,
        SV_Type = sv_type,
        N_Factors = 0,
        Marginal_R2 = NA_real_,
        P_value = NA_real_
      )
    )
  }
  
  cat(
    "\n", sv_type, " | ", classification_name,
    " | ", length(classification_factors), " factors\n",
    sep = ""
  )
  
  full_variables <- c(base_covars, all_factors)
  remaining_factors <- setdiff(all_factors, classification_factors)
  reduced_variables <- c(base_covars, remaining_factors)
  
  R2_full <- get_model_R2(dist_mat, metadata_df, full_variables)
  R2_reduced <- get_model_R2(dist_mat, metadata_df, reduced_variables)
  observed_R2 <- max(0, R2_full - R2_reduced)
  
  cat("Observed marginal R2 =", round(observed_R2, 6), "\n")
  
  perm_R2 <- numeric(nperm)
  
  for (b in seq_len(nperm)) {
    perm_metadata <- metadata_df
    perm_index <- sample(seq_len(nrow(metadata_df)))
    
    # Apply the same permutation index to every variable in the block.
    perm_metadata[, classification_factors] <- metadata_df[
      perm_index,
      classification_factors,
      drop = FALSE
    ]
    
    R2_perm_full <- get_model_R2(
      dist_mat,
      perm_metadata,
      full_variables
    )
    
    perm_R2[b] <- max(0, R2_perm_full - R2_reduced)
    
    if (b == 1 || b %% 100 == 0 || b == nperm) {
      cat(
        "\r", sv_type, " | ", classification_name,
        " | ", b, "/", nperm,
        " (", round(b / nperm * 100, 1), "%)",
        sep = ""
      )
      flush.console()
    }
  }
  
  cat("\n")
  
  p_value <- (sum(perm_R2 >= observed_R2) + 1) / (nperm + 1)
  
  data.frame(
    Classification = classification_name,
    SV_Type = sv_type,
    N_Factors = length(classification_factors),
    Marginal_R2 = observed_R2,
    P_value = p_value,
    stringsAsFactors = FALSE
  )
}

run_classification_analysis <- function(
    dist_mat,
    metadata_df,
    valid_factors,
    classification_for_model,
    sv_type,
    nperm = 999
) {
  diet_classes <- unique(classification_for_model$Classification)
  
  results <- map_df(seq_along(diet_classes), function(i) {
    class_name <- diet_classes[i]
    
    class_factors <- classification_for_model %>%
      filter(Classification == class_name) %>%
      pull(Factor)
    
    cat(
      "\n>>> ", sv_type, " ", i, "/", length(diet_classes),
      ": ", class_name, "\n",
      sep = ""
    )
    
    calc_classification_permanova(
      dist_mat = dist_mat,
      metadata_df = metadata_df,
      all_factors = valid_factors,
      classification_factors = class_factors,
      classification_name = class_name,
      sv_type = sv_type,
      nperm = nperm
    )
  })
  
  results %>%
    mutate(p_adj_BH = p.adjust(P_value, method = "BH")) %>%
    arrange(desc(Marginal_R2))
}

# =============================================================================
# 4. Species-abundance PCoA (calculated once)
# =============================================================================

cat("\n========== Build species-abundance PCoA ==========\n")

matched_clades <- unique(species_map$metaphlan_clade)

abundance_matrix <- metaphlan_raw %>%
  filter(clade_name %in% matched_clades) %>%
  left_join(species_map, by = c("clade_name" = "metaphlan_clade")) %>%
  select(-clade_name, -metaphlan_species) %>%
  pivot_longer(
    cols = -sv_species,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  mutate(
    sample_id = as.character(sample_id),
    abundance = suppressWarnings(as.numeric(abundance))
  ) %>%
  group_by(sample_id, sv_species) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    id_cols = sample_id,
    names_from = sv_species,
    values_from = abundance,
    values_fill = 0
  )

abun_numeric <- as.matrix(abundance_matrix %>% select(-sample_id))
rownames(abun_numeric) <- abundance_matrix$sample_id

abun_numeric <- abun_numeric[
  rowSums(abun_numeric) > 0,
  colSums(abun_numeric) > 0,
  drop = FALSE
]

dist_abun <- vegdist(abun_numeric, method = "bray")
pcoa_abun <- cmdscale(dist_abun, k = 3, eig = TRUE)

pcoa_df <- data.frame(
  sample_id = rownames(pcoa_abun$points),
  Abun_PCoA1 = pcoa_abun$points[, 1],
  Abun_PCoA2 = pcoa_abun$points[, 2],
  Abun_PCoA3 = pcoa_abun$points[, 3],
  stringsAsFactors = FALSE
)

cat("Abundance PCoA complete. Samples:", nrow(pcoa_df), "\n")

# =============================================================================
# 5. Metadata preparation (performed once)
# =============================================================================

all_factors <- unique(group_info$Factor)
model_vars <- c(
  "sample_id",
  "Abun_PCoA1",
  "Abun_PCoA2",
  "Abun_PCoA3",
  all_factors
)

metadata_complete <- metadata %>%
  left_join(pcoa_df, by = "sample_id") %>%
  select(any_of(model_vars))

for (var in intersect(all_factors, colnames(metadata_complete))) {
  if (is.character(metadata_complete[[var]])) {
    metadata_complete[[var]] <- factor(metadata_complete[[var]])
  }
}

metadata_complete <- metadata_complete %>% drop_na()

valid_factors <- all_factors[all_factors %in% colnames(metadata_complete)]
valid_factors <- valid_factors[
  sapply(metadata_complete[valid_factors], function(x) length(unique(x)) > 1)
]

valid_ids <- as.character(metadata_complete$sample_id)

cat("Valid samples:", length(valid_ids), "\n")
cat("Valid factors:", length(valid_factors), "\n")

# =============================================================================
# 6. Build dSV and vSV distance matrices once
# =============================================================================

cat("\n========== Build dSV distance matrix ==========\n")

feature_cols_dsv <- setdiff(colnames(dsv), "sample_id")
dsv_cleaned <- dsv

dsv_cleaned[feature_cols_dsv] <- lapply(
  dsv_cleaned[feature_cols_dsv],
  function(x) {
    x <- suppressWarnings(as.numeric(x))
    x[is.na(x)] <- 0
    x
  }
)

dsv_filtered <- dsv_cleaned[
  as.character(dsv_cleaned$sample_id) %in% valid_ids,
  ,
  drop = FALSE
]

dsv_filtered <- dsv_filtered[
  match(valid_ids, as.character(dsv_filtered$sample_id), nomatch = 0),
  ,
  drop = FALSE
]

dsv_mat <- as.matrix(dsv_filtered[, feature_cols_dsv, drop = FALSE])
storage.mode(dsv_mat) <- "numeric"
rownames(dsv_mat) <- as.character(dsv_filtered$sample_id)

dsv_mat <- dsv_mat[
  ,
  colSums(dsv_mat != 0, na.rm = TRUE) > 0,
  drop = FALSE
]

non_empty_samples_dsv <- rowSums(dsv_mat != 0, na.rm = TRUE) > 0

cat(
  sprintf(
    "dSV: original samples %d, all-zero samples removed %d, remaining %d\n",
    nrow(dsv_mat),
    sum(!non_empty_samples_dsv),
    sum(non_empty_samples_dsv)
  )
)

dsv_mat_valid <- dsv_mat[non_empty_samples_dsv, , drop = FALSE]

metadata_dsv_complete <- metadata_complete[
  as.character(metadata_complete$sample_id) %in% rownames(dsv_mat_valid),
  ,
  drop = FALSE
]

metadata_dsv_complete <- metadata_dsv_complete[
  match(rownames(dsv_mat_valid), as.character(metadata_dsv_complete$sample_id)),
  ,
  drop = FALSE
]

dsv_mat_valid <- dsv_mat_valid[
  as.character(metadata_dsv_complete$sample_id),
  ,
  drop = FALSE
]

dist_dsv <- vegdist(dsv_mat_valid, method = "jaccard", binary = TRUE)

cat("dSV distance matrix complete. Samples:", nrow(metadata_dsv_complete), "\n")

cat("\n========== Build vSV distance matrix ==========\n")

feature_cols_vsv <- setdiff(colnames(vsv), "sample_id")

vsv_numeric <- vsv[, feature_cols_vsv, drop = FALSE]
vsv_numeric[] <- lapply(vsv_numeric, min_max_scale)

vsv_processed <- data.frame(
  sample_id = vsv$sample_id,
  vsv_numeric,
  check.names = FALSE
)

vsv_filtered <- vsv_processed[
  as.character(vsv_processed$sample_id) %in% valid_ids,
  ,
  drop = FALSE
]

vsv_filtered <- vsv_filtered[
  match(valid_ids, as.character(vsv_filtered$sample_id), nomatch = 0),
  ,
  drop = FALSE
]

vsv_mat <- as.matrix(vsv_filtered[, feature_cols_vsv, drop = FALSE])
storage.mode(vsv_mat) <- "numeric"
rownames(vsv_mat) <- as.character(vsv_filtered$sample_id)

vsv_mat <- vsv_mat[
  ,
  colSums(vsv_mat != 0, na.rm = TRUE) > 0,
  drop = FALSE
]

non_empty_samples_vsv <- rowSums(vsv_mat != 0, na.rm = TRUE) > 0

cat(
  sprintf(
    "vSV: original samples %d, all-zero samples removed %d, remaining %d\n",
    nrow(vsv_mat),
    sum(!non_empty_samples_vsv),
    sum(non_empty_samples_vsv)
  )
)

vsv_mat_valid <- vsv_mat[non_empty_samples_vsv, , drop = FALSE]

metadata_vsv_complete <- metadata_complete[
  as.character(metadata_complete$sample_id) %in% rownames(vsv_mat_valid),
  ,
  drop = FALSE
]

metadata_vsv_complete <- metadata_vsv_complete[
  match(rownames(vsv_mat_valid), as.character(metadata_vsv_complete$sample_id)),
  ,
  drop = FALSE
]

vsv_mat_valid <- vsv_mat_valid[
  as.character(metadata_vsv_complete$sample_id),
  ,
  drop = FALSE
]

dist_vsv <- vegdist(vsv_mat_valid, method = "euclidean")

cat("vSV distance matrix complete. Samples:", nrow(metadata_vsv_complete), "\n")

# =============================================================================
# 7. Overall group-level PERMANOVA
# =============================================================================

cat("\n========== Overall group-level PERMANOVA ==========\n")

dsv_results <- run_permanova_group_analysis(
  dist_dsv,
  metadata_dsv_complete,
  valid_factors,
  group_info,
  "dSV"
)

vsv_results <- run_permanova_group_analysis(
  dist_vsv,
  metadata_vsv_complete,
  valid_factors,
  group_info,
  "vSV"
)

write_csv(dsv_results, "permanova_dsv_results.csv")
write_csv(vsv_results, "permanova_vsv_results.csv")
write_csv(
  bind_rows(dsv_results, vsv_results),
  "permanova_vsv_dsv_bray_results.csv"
)

# =============================================================================
# 8. urabnization factor-level PERMANOVA
# =============================================================================

cat("\n========== urabnization factor-level PERMANOVA ==========\n")

urabnization_dsv <- run_factor_group(
  "urabnization",
  dist_dsv,
  metadata_dsv_complete,
  valid_factors,
  group_info,
  "dSV"
)

urabnization_vsv <- run_factor_group(
  "urabnization",
  dist_vsv,
  metadata_vsv_complete,
  valid_factors,
  group_info,
  "vSV"
)

write_csv(
  urabnization_dsv,
  "permanova_urabnization_factors_dsv.csv"
)

write_csv(
  urabnization_vsv,
  "permanova_urabnization_factors_vsv.csv"
)

write_csv(
  bind_rows(urabnization_dsv, urabnization_vsv),
  "permanova_urabnization_factors_combined.csv"
)

# =============================================================================
# 9. Prepare Diet classification mapping
# =============================================================================

factor_candidates <- c(
  "Factor", "factor", "Variable", "variable",
  "Diet_factor", "diet_factor"
)

class_candidates <- c(
  "Classification", "classification",
  "Category", "category",
  "Class", "class",
  "Group", "group"
)

factor_col <- factor_candidates[
  factor_candidates %in% colnames(classification)
][1]

class_col <- class_candidates[
  class_candidates %in% colnames(classification)
][1]

if (is.na(factor_col) || is.na(class_col)) {
  stop(
    "Could not identify Factor or Classification column in classification.csv. Current columns: ",
    paste(colnames(classification), collapse = ", ")
  )
}

classification_clean <- classification %>%
  transmute(
    Factor_classification = as.character(.data[[factor_col]]),
    Classification = as.character(.data[[class_col]]),
    Factor_key = normalize_name(.data[[factor_col]])
  ) %>%
  filter(
    !is.na(Factor_key),
    Factor_key != "",
    !is.na(Classification),
    Classification != ""
  ) %>%
  distinct(Factor_key, .keep_all = TRUE)

diet_factor_df <- group_info %>%
  filter(Group == "Diet") %>%
  transmute(
    Factor = as.character(Factor),
    Factor_key = normalize_name(Factor)
  ) %>%
  distinct() %>%
  filter(Factor %in% valid_factors)

diet_factors <- diet_factor_df$Factor

if (length(diet_factors) == 0) {
  stop("No valid Diet factors found. Check group.csv and metadata column names.")
}

diet_classification_map <- diet_factor_df %>%
  left_join(
    classification_clean %>% select(Factor_key, Classification),
    by = "Factor_key"
  )

cat(
  "Diet factors matched to classification.csv:",
  sum(!is.na(diet_classification_map$Classification)),
  "/",
  nrow(diet_classification_map),
  "\n"
)

unclassified_factors <- diet_classification_map %>%
  filter(is.na(Classification))

if (nrow(unclassified_factors) > 0) {
  cat("\nDiet factors still unmatched after normalization:\n")
  print(unclassified_factors$Factor)
}

# =============================================================================
# 10. Diet factor-level PERMANOVA
# =============================================================================

cat("\n========== Diet factor-level PERMANOVA ==========\n")

diet_results_dsv <- run_factor_group(
  "Diet",
  dist_dsv,
  metadata_dsv_complete,
  valid_factors,
  group_info,
  "dSV"
)

diet_results_vsv <- run_factor_group(
  "Diet",
  dist_vsv,
  metadata_vsv_complete,
  valid_factors,
  group_info,
  "vSV"
)

write_csv(diet_results_dsv, "permanova_diet_factors_dsv.csv")
write_csv(diet_results_vsv, "permanova_diet_factors_vsv.csv")

diet_results_merged <- full_join(
  diet_results_dsv %>%
    select(
      Factor,
      dSV_Marginal_R2 = Marginal_R2,
      dSV_P_value = P_value,
      dSV_p_adj = p_adj
    ),
  diet_results_vsv %>%
    select(
      Factor,
      vSV_Marginal_R2 = Marginal_R2,
      vSV_P_value = P_value,
      vSV_p_adj = p_adj
    ),
  by = "Factor"
) %>%
  mutate(Factor_key = normalize_name(Factor)) %>%
  left_join(
    classification_clean %>% select(Factor_key, Classification),
    by = "Factor_key"
  ) %>%
  rowwise() %>%
  mutate(
    Mean_Marginal_R2 = mean(
      c(dSV_Marginal_R2, vSV_Marginal_R2),
      na.rm = TRUE
    ),
    Combined_P_value = combine_two_p(dSV_P_value, vSV_P_value)
  ) %>%
  ungroup() %>%
  mutate(
    Combined_p_adj_BH = p.adjust(Combined_P_value, method = "BH"),
    Classification = ifelse(
      is.na(Classification),
      "Unclassified",
      Classification
    )
  ) %>%
  select(
    Classification,
    Factor,
    dSV_Marginal_R2,
    vSV_Marginal_R2,
    Mean_Marginal_R2,
    dSV_P_value,
    vSV_P_value,
    dSV_p_adj,
    vSV_p_adj,
    Combined_P_value,
    Combined_p_adj_BH
  ) %>%
  arrange(Classification, desc(Mean_Marginal_R2))

write_csv(
  diet_results_merged,
  "permanova_diet_factors_dSV_vSV_with_classification.csv"
)

# =============================================================================
# 11. Diet classification-level PERMANOVA
# =============================================================================

cat("\n========== Diet classification-level PERMANOVA ==========\n")

classification_for_model <- classification_clean %>%
  inner_join(
    diet_factor_df %>% select(Factor, Factor_key),
    by = "Factor_key"
  ) %>%
  filter(Factor %in% valid_factors) %>%
  distinct(Classification, Factor)

cat("\nFactors entering classification-level PERMANOVA:\n")
print(
  classification_for_model %>%
    count(Classification, name = "N_Factors")
)

classification_results_dsv <- run_classification_analysis(
  dist_dsv,
  metadata_dsv_complete,
  valid_factors,
  classification_for_model,
  "dSV",
  nperm = 999
)

classification_results_vsv <- run_classification_analysis(
  dist_vsv,
  metadata_vsv_complete,
  valid_factors,
  classification_for_model,
  "vSV",
  nperm = 999
)

write_csv(
  classification_results_dsv,
  "PERMANOVA_Diet_Classification_dSV.csv"
)

write_csv(
  classification_results_vsv,
  "PERMANOVA_Diet_Classification_vSV.csv"
)

classification_results_merged <- full_join(
  classification_results_dsv %>%
    select(
      Classification,
      N_Factors_dSV = N_Factors,
      dSV_Marginal_R2 = Marginal_R2,
      dSV_P_value = P_value,
      dSV_p_adj_BH = p_adj_BH
    ),
  classification_results_vsv %>%
    select(
      Classification,
      N_Factors_vSV = N_Factors,
      vSV_Marginal_R2 = Marginal_R2,
      vSV_P_value = P_value,
      vSV_p_adj_BH = p_adj_BH
    ),
  by = "Classification"
) %>%
  rowwise() %>%
  mutate(
    Mean_Marginal_R2 = mean(
      c(dSV_Marginal_R2, vSV_Marginal_R2),
      na.rm = TRUE
    ),
    Combined_P_value = combine_two_p(dSV_P_value, vSV_P_value)
  ) %>%
  ungroup() %>%
  mutate(
    Combined_p_adj_BH = p.adjust(Combined_P_value, method = "BH")
  ) %>%
  arrange(desc(Mean_Marginal_R2))

write_csv(
  classification_results_merged,
  "PERMANOVA_Diet_Classification_dSV_vSV_merged.csv"
)

# =============================================================================
# 12. Final output
# =============================================================================

cat("\n============================================================\n")
cat("Diet classification-level PERMANOVA\n")
cat("============================================================\n")

print(
  classification_results_merged %>%
    select(
      Classification,
      N_Factors_dSV,
      N_Factors_vSV,
      dSV_Marginal_R2,
      vSV_Marginal_R2,
      Mean_Marginal_R2,
      dSV_P_value,
      vSV_P_value,
      Combined_P_value,
      Combined_p_adj_BH
    )
)

cat("\n========== PERMANOVA analysis complete ==========\n")
