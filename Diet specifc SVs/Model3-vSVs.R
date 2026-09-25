# ============================================================================
# Diet-associated vSV analysis
# Stage 1: raw linear screening for each dietary factor
# Stage 2: multivariable linear regression with species coverage as a covariate
# Only Version A is included: no coverage cutoff is applied
# ============================================================================

set.seed(12345)

library(tidyverse)
library(broom)
library(openxlsx)
library(caret)
library(stringr)

if (!dir.exists("results_vsv")) dir.create("results_vsv")

safe_write_xlsx <- function(df, path) {
  tryCatch(
    write.xlsx(df, path, rowNames = FALSE),
    error = function(e) {
      alt <- gsub("\\.xlsx$", paste0("_backup_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"), path)
      write.xlsx(df, alt, rowNames = FALSE)
    }
  )
}

normalize_species_name <- function(x) {
  species <- str_extract(x, "^[A-Za-z]+_[a-z]+")
  as.character(species)[1]
}

# ----------------------------------------------------------------------------
# 1. Read species abundance and coverage data
# ----------------------------------------------------------------------------

cat("========== 1. Read species abundance and coverage data ==========\n")

species_map <- read.csv("sv_to_metaphlan_final_mapping.csv", check.names = FALSE)
metaphlan_raw <- read_tsv(
  "yunnan_7_years_ago_combined_metaphlan_profile.tsv",
  comment = "",
  show_col_types = FALSE
)

if (grepl("#", colnames(metaphlan_raw)[1])) colnames(metaphlan_raw)[1] <- "clade_name"

matched_clades <- unique(species_map$metaphlan_clade)

abundance_matrix <- metaphlan_raw %>%
  filter(clade_name %in% matched_clades) %>%
  left_join(species_map, by = c("clade_name" = "metaphlan_clade")) %>%
  select(-clade_name, -metaphlan_species) %>%
  pivot_longer(cols = -sv_species, names_to = "sample_id", values_to = "abundance") %>%
  mutate(abundance = suppressWarnings(as.numeric(abundance))) %>%
  group_by(sample_id, sv_species) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(id_cols = sample_id, names_from = sv_species, values_from = abundance, values_fill = 0)

raw_cov_data <- read.csv("merged_species_mean_coverage_matrix.csv", check.names = FALSE)
colnames(raw_cov_data)[1] <- "sample_id"

coverage_matrix <- raw_cov_data %>%
  pivot_longer(cols = -sample_id, names_to = "sv_species", values_to = "coverage") %>%
  mutate(coverage = suppressWarnings(as.numeric(coverage))) %>%
  group_by(sample_id, sv_species) %>%
  summarise(coverage = mean(coverage, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(id_cols = sample_id, names_from = sv_species, values_from = coverage, values_fill = 0)

colnames(coverage_matrix)[-1] <- paste0(colnames(coverage_matrix)[-1], "_coverage")

cat("Abundance matrix:", nrow(abundance_matrix), "samples x",
    ncol(abundance_matrix) - 1, "species\n")
cat("Coverage matrix:", nrow(coverage_matrix), "samples x",
    ncol(coverage_matrix) - 1, "species\n")

# ----------------------------------------------------------------------------
# 2. Read vSV, metadata, diet, and covariates
# ----------------------------------------------------------------------------

cat("\n========== 2. Read vSV and metadata ==========\n")

vsv_data <- read.csv(
  "final_filtered_vsv.csv",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Ensure unique vSV feature names without removing any features
if (anyDuplicated(colnames(vsv_data)) > 0) {
  cat("Duplicate vSV feature names detected; making names unique.\n")
  colnames(vsv_data) <- make.unique(colnames(vsv_data), sep = "__dup")
}

metadata <- read.csv("metadata.csv", check.names = FALSE)
diet_data <- read.csv("final_filtered_diet.csv", check.names = FALSE)
urbanization_data <- read.csv("Urbanization.csv", check.names = FALSE)
medication_data <- read.csv("Medication.csv", check.names = FALSE)
general_metadata_data <- read.csv("General metadata.csv", check.names = FALSE)

required_meta_cols <- c("sample_id", "residency", "Age", "Sex", "BMI", "ethnicity")
missing_cols <- setdiff(required_meta_cols, names(metadata))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

diet_factors <- diet_data %>%
  select(-any_of(c("sample_id", "ethnicity", "residency", "region"))) %>%
  select(where(is.numeric)) %>%
  colnames()

cat("Number of dietary factors:", length(diet_factors), "\n")

# ----------------------------------------------------------------------------
# 3. Prepare covariates
# ----------------------------------------------------------------------------

process_covariate <- function(data, covariate_name) {
  numeric_cols <- data %>%
    select(-any_of(c("sample_id", "ethnicity", "residency", "region"))) %>%
    select(where(is.numeric))
  
  if (ncol(numeric_cols) == 0) {
    return(data %>% select(sample_id) %>% mutate(!!covariate_name := NA_real_))
  }
  
  data.frame(
    sample_id = data$sample_id,
    score = rowSums(numeric_cols, na.rm = TRUE)
  ) %>%
    rename(!!covariate_name := score)
}

reverse_ordinal <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  if (all(is.na(x_num))) return(x)
  rng <- range(x_num, na.rm = TRUE)
  rng[1] + rng[2] - x_num
}

# Reverse-code animal contact before constructing the Urbanization composite
animal_col <- names(urbanization_data)[
  tolower(trimws(names(urbanization_data))) == "contact with animals"
]

if (length(animal_col) == 1) {
  urbanization_data[[animal_col]] <- reverse_ordinal(urbanization_data[[animal_col]])
}

medication_cov <- process_covariate(medication_data, "medication")
urbanization_cov <- process_covariate(urbanization_data, "urbanization")

# Keep each general metadata variable as an individual covariate
general_metadata_cov <- general_metadata_data
general_cols_original <- setdiff(names(general_metadata_cov), "sample_id")
general_covars <- paste0("general__", make.names(general_cols_original, unique = TRUE))
names(general_metadata_cov)[match(general_cols_original, names(general_metadata_cov))] <- general_covars

general_metadata_cov <- general_metadata_cov %>%
  mutate(across(all_of(general_covars), function(x) {
    if (is.numeric(x)) return(x)
    
    x_char <- trimws(as.character(x))
    x_lower <- tolower(x_char)
    
    if (all(x_lower %in% c("yes", "no", "y", "n", "true", "false", "1", "0", "", NA))) {
      return(case_when(
        x_lower %in% c("yes", "y", "true", "1") ~ 1,
        x_lower %in% c("no", "n", "false", "0") ~ 0,
        TRUE ~ NA_real_
      ))
    }
    
    numeric_x <- suppressWarnings(as.numeric(x_char))
    non_missing_original <- !is.na(x_char) & x_char != ""
    
    if (sum(!is.na(numeric_x)) == sum(non_missing_original)) return(numeric_x)
    
    factor(x_char)
  }))

# Remove any dietary columns already present in metadata; diet variables are taken from final_filtered_diet.csv
metadata_base <- metadata %>%
  select(-any_of(intersect(diet_factors, colnames(metadata))))

diet_data_join <- diet_data %>%
  select(sample_id, all_of(diet_factors))

combined_metadata <- metadata_base %>%
  inner_join(diet_data_join, by = "sample_id") %>%
  inner_join(medication_cov, by = "sample_id") %>%
  inner_join(urbanization_cov, by = "sample_id") %>%
  inner_join(general_metadata_cov, by = "sample_id") %>%
  left_join(abundance_matrix, by = "sample_id") %>%
  left_join(coverage_matrix, by = "sample_id") %>%
  mutate(
    Sex = as.factor(Sex),
    ethnicity = as.factor(ethnicity),
    residency = as.factor(residency)
  )

cat("Samples after metadata merge:", nrow(combined_metadata), "\n")
cat("Dietary factors matched:", sum(diet_factors %in% colnames(combined_metadata)),
    "/", length(diet_factors), "\n")

# ----------------------------------------------------------------------------
# 4. vSV quality control
# ----------------------------------------------------------------------------

cat("\n========== 3. vSV quality control ==========\n")

valid_vsv <- colnames(vsv_data)[sapply(vsv_data, function(x) {
  x <- suppressWarnings(as.numeric(x))
  mean(!is.na(x)) >= 0.10 &&
    !is.na(sd(x, na.rm = TRUE)) &&
    sd(x, na.rm = TRUE) > 0
})]

vsv_filtered <- vsv_data[, valid_vsv, drop = FALSE]

cat("vSVs retained:", ncol(vsv_data), "->", ncol(vsv_filtered), "\n")

# ----------------------------------------------------------------------------
# 5. Stage 1: raw linear screening
# ----------------------------------------------------------------------------

cat("\n========== 4. Raw linear screening ==========\n")

run_raw_vsv <- function(diet_var, meta_df, vsv_df) {
  out <- vector("list", ncol(vsv_df))
  
  for (j in seq_len(ncol(vsv_df))) {
    vsv_id <- colnames(vsv_df)[j]
    
    df <- meta_df %>%
      mutate(
        vsv_value = suppressWarnings(as.numeric(vsv_df[[vsv_id]][match(sample_id, rownames(vsv_df))])),
        target_diet_factor = .data[[diet_var]]
      ) %>%
      drop_na(vsv_value, target_diet_factor)
    
    if (nrow(df) < 30 ||
        length(unique(df$vsv_value)) < 2 ||
        length(unique(df$target_diet_factor)) < 2) next
    
    fit <- tryCatch(
      lm(vsv_value ~ target_diet_factor, data = df),
      error = function(e) NULL
    )
    
    if (is.null(fit)) next
    
    tmp <- tidy(fit) %>%
      filter(term == "target_diet_factor")
    
    if (nrow(tmp) == 0) next
    
    out[[j]] <- tmp %>%
      transmute(
        diet_factor = diet_var,
        vsv = vsv_id,
        N_raw = nrow(df),
        raw_beta = estimate,
        raw_std_error = std.error,
        raw_statistic = statistic,
        raw_p = p.value
      )
  }
  
  bind_rows(out)
}

raw_vsv_list <- vector("list", length(diet_factors))

for (k in seq_along(diet_factors)) {
  cat(sprintf("[RAW %d/%d] %s\n", k, length(diet_factors), diet_factors[k]))
  raw_vsv_list[[k]] <- run_raw_vsv(
    diet_factors[k],
    combined_metadata,
    vsv_filtered
  )
}

raw_vsv_results <- bind_rows(raw_vsv_list)

raw_vsv_candidates <- raw_vsv_results %>%
  filter(!is.na(raw_p), raw_p < 0.05) %>%
  arrange(diet_factor, raw_p)

safe_write_xlsx(
  raw_vsv_results,
  "results_vsv/all_diet_vsv_raw_screening.xlsx"
)

safe_write_xlsx(
  raw_vsv_candidates,
  "results_vsv/diet_vsv_raw_p005_candidates.xlsx"
)

cat("Diet-vSV pairs passing raw P < 0.05:", nrow(raw_vsv_candidates), "\n")

# ----------------------------------------------------------------------------
# 6. Stage 2: adjusted Version A model
# ----------------------------------------------------------------------------

run_adjusted_vsv <- function(diet_var, vsv_id, meta_df, vsv_df) {
  if (!vsv_id %in% colnames(vsv_df)) return(NULL)
  
  species <- normalize_species_name(vsv_id)
  if (is.na(species) || species == "") return(NULL)
  
  spec_candidates <- c(species, paste0("s__", species), gsub("_", " ", species))
  spec_col <- colnames(meta_df)[colnames(meta_df) %in% spec_candidates][1]
  
  if (is.na(spec_col)) return(NULL)
  
  cov_col <- paste0(species, "_coverage")
  
  df <- meta_df %>%
    mutate(
      vsv_value = suppressWarnings(as.numeric(vsv_df[[vsv_id]][match(sample_id, rownames(vsv_df))])),
      target_diet_factor = .data[[diet_var]],
      species_abundance = .data[[spec_col]],
      species_coverage = if (cov_col %in% colnames(.)) .data[[cov_col]] else NA_real_
    ) %>%
    drop_na(
      vsv_value, target_diet_factor, Age, Sex, BMI, ethnicity, residency,
      medication, urbanization, all_of(general_covars),
      species_abundance, species_coverage
    ) %>%
    droplevels()
  
  if (nrow(df) < 30 ||
      length(unique(df$vsv_value)) < 2 ||
      length(unique(df$target_diet_factor)) < 2) return(NULL)
  
  model_covars <- c(
    "target_diet_factor", "species_coverage", "species_abundance",
    "Age", "Sex", "BMI", "ethnicity", "residency",
    "medication", "urbanization", general_covars
  )
  
  valid_covars <- model_covars[sapply(model_covars, function(v) {
    if (!v %in% colnames(df)) return(FALSE)
    x <- df[[v]]
    
    if (is.factor(x) || is.character(x)) {
      length(unique(na.omit(as.character(x)))) >= 2
    } else {
      length(unique(na.omit(x))) >= 2
    }
  })]
  
  if (!"target_diet_factor" %in% valid_covars) return(NULL)
  
  model_formula <- as.formula(
    paste("vsv_value ~", paste(sprintf("`%s`", valid_covars), collapse = " + "))
  )
  
  fit <- tryCatch(
    lm(model_formula, data = df),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NULL)
  
  tmp <- tidy(fit) %>%
    filter(term == "target_diet_factor")
  
  if (nrow(tmp) == 0) return(NULL)
  
  tmp %>%
    transmute(
      diet_factor = diet_var,
      vsv = vsv_id,
      species = species,
      N = nrow(df),
      Beta = estimate,
      std.error = std.error,
      statistic = statistic,
      p.value = p.value,
      CI_low = estimate - 1.96 * std.error,
      CI_high = estimate + 1.96 * std.error
    )
}

cat("\n========== 5. Adjusted model (Version A) ==========\n")

adjusted_A <- vector("list", nrow(raw_vsv_candidates))

if (nrow(raw_vsv_candidates) > 0) {
  pb <- txtProgressBar(min = 0, max = nrow(raw_vsv_candidates), style = 3)
  
  for (i in seq_len(nrow(raw_vsv_candidates))) {
    adjusted_A[[i]] <- run_adjusted_vsv(
      raw_vsv_candidates$diet_factor[i],
      raw_vsv_candidates$vsv[i],
      combined_metadata,
      vsv_filtered
    )
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
}

all_vsv_A <- bind_rows(adjusted_A)

if (nrow(all_vsv_A) > 0) {
  all_vsv_A <- all_vsv_A %>%
    left_join(
      raw_vsv_candidates %>%
        select(diet_factor, vsv, N_raw, raw_beta, raw_std_error, raw_statistic, raw_p),
      by = c("diet_factor", "vsv")
    ) %>%
    group_by(diet_factor) %>%
    mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
    ungroup() %>%
    arrange(diet_factor, p_adj)
}

sig_vsv_A <- all_vsv_A %>%
  filter(!is.na(p_adj), p_adj < 0.05)

safe_write_xlsx(
  all_vsv_A,
  "results_vsv/all_diet_factors_vsv_adjusted_verA.xlsx"
)

safe_write_xlsx(
  sig_vsv_A,
  "results_vsv/significant_diet_factors_vsv_verA.xlsx"
)

cat("Version A associations with FDR < 0.05:", nrow(sig_vsv_A), "\n")

# ----------------------------------------------------------------------------
# 7. Ten-fold cross-validation and marginal R2
# ----------------------------------------------------------------------------

run_cv_vsv <- function(diet_var, vsv_name, meta_df, vsv_df) {
  if (!vsv_name %in% colnames(vsv_df)) return(NULL)
  
  species <- normalize_species_name(vsv_name)
  if (is.na(species) || species == "") return(NULL)
  
  spec_candidates <- c(species, paste0("s__", species), gsub("_", " ", species))
  spec_col <- colnames(meta_df)[colnames(meta_df) %in% spec_candidates][1]
  
  if (is.na(spec_col)) return(NULL)
  
  cov_col <- paste0(species, "_coverage")
  
  df <- meta_df %>%
    mutate(
      vsv_value = suppressWarnings(as.numeric(vsv_df[[vsv_name]][match(sample_id, rownames(vsv_df))])),
      target_diet_factor = .data[[diet_var]],
      species_abundance = .data[[spec_col]],
      species_coverage = if (cov_col %in% colnames(.)) .data[[cov_col]] else NA_real_
    ) %>%
    drop_na(
      vsv_value, target_diet_factor, Age, Sex, BMI, ethnicity, residency,
      medication, urbanization, all_of(general_covars),
      species_abundance, species_coverage
    ) %>%
    droplevels()
  
  if (nrow(df) < 30 || length(unique(df$vsv_value)) < 2) return(NULL)
  
  model_covars <- c(
    "target_diet_factor", "species_coverage", "species_abundance",
    "Age", "Sex", "BMI", "ethnicity", "residency",
    "medication", "urbanization", general_covars
  )
  
  valid_covars <- model_covars[sapply(model_covars, function(v) {
    if (!v %in% colnames(df)) return(FALSE)
    x <- df[[v]]
    
    if (is.factor(x) || is.character(x)) {
      length(unique(na.omit(as.character(x)))) >= 2
    } else {
      length(unique(na.omit(x))) >= 2
    }
  })]
  
  if (!"target_diet_factor" %in% valid_covars) return(NULL)
  
  fml <- as.formula(
    paste("vsv_value ~", paste(sprintf("`%s`", valid_covars), collapse = " + "))
  )
  
  fit <- tryCatch(
    lm(fml, data = df),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NULL)
  
  full_r2 <- summary(fit)$r.squared
  
  delta_r2 <- function(v) {
    if (!v %in% valid_covars) return(0)
    
    reduced_covars <- setdiff(valid_covars, v)
    if (length(reduced_covars) == 0) return(full_r2)
    
    reduced_formula <- as.formula(
      paste("vsv_value ~", paste(sprintf("`%s`", reduced_covars), collapse = " + "))
    )
    
    reduced_model <- tryCatch(
      lm(reduced_formula, data = df),
      error = function(e) NULL
    )
    
    if (is.null(reduced_model)) return(NA_real_)
    
    max(0, full_r2 - summary(reduced_model)$r.squared)
  }
  
  ctrl <- trainControl(
    method = "cv",
    number = 10,
    savePredictions = TRUE
  )
  
  cv <- tryCatch(
    train(
      fml,
      data = df,
      method = "lm",
      trControl = ctrl,
      metric = "Rsquared"
    ),
    error = function(e) NULL
  )
  
  if (is.null(cv)) return(NULL)
  
  result <- data.frame(
    diet_factor = diet_var,
    vsv = vsv_name,
    species = species,
    N = nrow(df),
    Full_Model_R2 = round(full_r2, 4),
    CV_Mean_Rsquared = round(mean(cv$results$Rsquared, na.rm = TRUE), 4),
    CV_Mean_RMSE = round(mean(cv$results$RMSE, na.rm = TRUE), 4),
    CV_Mean_MAE = round(mean(cv$results$MAE, na.rm = TRUE), 4),
    Marginal_Delta_R2_TargetDietFactor = round(delta_r2("target_diet_factor"), 4),
    Marginal_Delta_R2_SpeciesCoverage = round(delta_r2("species_coverage"), 4),
    Marginal_Delta_R2_SpeciesAbundance = round(delta_r2("species_abundance"), 4),
    Marginal_Delta_R2_Ethnicity = round(delta_r2("ethnicity"), 4),
    Marginal_Delta_R2_Residency = round(delta_r2("residency"), 4),
    Marginal_Delta_R2_Urbanization = round(delta_r2("urbanization"), 4),
    Marginal_Delta_R2_Medication = round(delta_r2("medication"), 4),
    Marginal_Delta_R2_Age = round(delta_r2("Age"), 4),
    Marginal_Delta_R2_Sex = round(delta_r2("Sex"), 4),
    Marginal_Delta_R2_BMI = round(delta_r2("BMI"), 4),
    stringsAsFactors = FALSE
  )
  
  for (g in general_covars) {
    result[[paste0("Marginal_Delta_R2_", g)]] <- round(delta_r2(g), 4)
  }
  
  result
}

if (nrow(sig_vsv_A) > 0) {
  cat("\n========== 6. Ten-fold CV and marginal R2 ==========\n")
  
  cv_list_A <- vector("list", nrow(sig_vsv_A))
  pb <- txtProgressBar(min = 0, max = nrow(sig_vsv_A), style = 3)
  
  for (i in seq_len(nrow(sig_vsv_A))) {
    cv_list_A[[i]] <- run_cv_vsv(
      sig_vsv_A$diet_factor[i],
      sig_vsv_A$vsv[i],
      combined_metadata,
      vsv_filtered
    )
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  cv_A <- bind_rows(cv_list_A)
  
  if (nrow(cv_A) > 0) {
    safe_write_xlsx(
      cv_A,
      "results_vsv/cv_10fold_diet_vsv_verA.xlsx"
    )
  }
}

# ----------------------------------------------------------------------------
# 8. Summary of vSV counts for each dietary factor
# ----------------------------------------------------------------------------

summary_A <- data.frame(
  Diet = diet_factors,
  Raw_P_lt_0.05_vSV_N = sapply(diet_factors, function(x) {
    raw_vsv_candidates %>%
      filter(diet_factor == x) %>%
      distinct(vsv) %>%
      nrow()
  }),
  Adjusted_Significant_vSV_N = sapply(diet_factors, function(x) {
    sig_vsv_A %>%
      filter(diet_factor == x) %>%
      distinct(vsv) %>%
      nrow()
  })
) %>%
  arrange(desc(Raw_P_lt_0.05_vSV_N))

write.csv(
  summary_A,
  "results_vsv/diet_vsv_raw_and_adjusted_counts_verA.csv",
  row.names = FALSE,
  quote = FALSE
)

cat("\n========== Analysis completed ==========\n")
cat("Raw P < 0.05 candidate pairs:", nrow(raw_vsv_candidates), "\n")
cat("Version A associations with FDR < 0.05:", nrow(sig_vsv_A), "\n")
cat("Results saved in the 'results_vsv' directory.\n")
