# ============================================================================
# Diet-associated dSV analysis
# Stage 1: raw logistic screening for each dietary factor
# Stage 2: multivariable logistic regression with species coverage as a covariate
# Only Version A is included: no coverage cutoff is applied
# ============================================================================

set.seed(12345)

library(tidyverse)
library(broom)
library(openxlsx)
library(brglm2)
library(caret)
library(pROC)

if (!dir.exists("results")) dir.create("results")

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
# 2. Read dSV, metadata, diet, and covariates
# ----------------------------------------------------------------------------

cat("\n========== 2. Read dSV and metadata ==========\n")

dsv_data <- read.csv(
  "final_filtered_dsv.csv",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Ensure unique dSV feature names without removing any features
if (anyDuplicated(colnames(dsv_data)) > 0) {
  cat("Duplicate dSV feature names detected; making names unique.\n")
  colnames(dsv_data) <- make.unique(colnames(dsv_data), sep = "__dup")
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

combined_metadata <- metadata %>%
  inner_join(diet_data, by = "sample_id") %>%
  inner_join(medication_cov, by = "sample_id") %>%
  inner_join(urbanization_cov, by = "sample_id") %>%
  inner_join(general_metadata_cov, by = "sample_id") %>%
  left_join(abundance_matrix, by = "sample_id") %>%
  left_join(coverage_matrix, by = "sample_id") %>%
  mutate(
    Sex = as.factor(Sex),
    residency = as.factor(residency),
    ethnicity = as.factor(ethnicity)
  )

cat("Samples after metadata merge:", nrow(combined_metadata), "\n")

# ----------------------------------------------------------------------------
# 4. dSV prevalence filtering
# ----------------------------------------------------------------------------

cat("\n========== 3. dSV prevalence filtering ==========\n")

dsv_data_prevalence <- dsv_data
dsv_data_prevalence[] <- lapply(dsv_data_prevalence, function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 0
  x
})

dsv_prevalence <- data.frame(
  dsv = colnames(dsv_data_prevalence),
  count_1 = colSums(dsv_data_prevalence == 1, na.rm = TRUE),
  count_0 = colSums(dsv_data_prevalence == 0, na.rm = TRUE),
  prevalence = colMeans(dsv_data_prevalence == 1, na.rm = TRUE)
)

valid_dsv <- dsv_prevalence %>%
  filter(prevalence >= 0.1, count_1 >= 10, count_0 >= 10) %>%
  pull(dsv)

dsv_filtered <- dsv_data[, valid_dsv, drop = FALSE]

dsv_long <- dsv_filtered %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(cols = -sample_id, names_to = "dsv", values_to = "status") %>%
  mutate(status = as.numeric(status))

combined_data <- dsv_long %>%
  inner_join(combined_metadata, by = "sample_id")

cat("dSVs retained after prevalence filtering:", length(valid_dsv), "\n")

# ----------------------------------------------------------------------------
# 5. Stage 1: raw logistic screening
# ----------------------------------------------------------------------------

cat("\n========== 4. Raw logistic screening ==========\n")

run_raw_screen <- function(diet_var, data) {
  dsv_list <- unique(data$dsv)
  out <- vector("list", length(dsv_list))
  
  for (i in seq_along(dsv_list)) {
    dsv_id <- dsv_list[i]
    
    df <- data %>%
      filter(dsv == dsv_id) %>%
      mutate(target_diet_factor = .data[[diet_var]]) %>%
      drop_na(status, target_diet_factor) %>%
      filter(status %in% c(0, 1))
    
    if (nrow(df) < 30 ||
        length(unique(df$status)) < 2 ||
        length(unique(df$target_diet_factor)) < 2) next
    
    cc <- table(df$status)
    if (length(cc) < 2 || any(cc < 10)) next
    
    fit <- tryCatch(
      glm(status ~ target_diet_factor, data = df, family = binomial()),
      error = function(e) NULL
    )
    
    if (is.null(fit)) next
    
    tmp <- broom::tidy(fit) %>%
      filter(term == "target_diet_factor")
    
    if (nrow(tmp) == 0) next
    
    out[[i]] <- tmp %>%
      transmute(
        diet_factor = diet_var,
        dsv = dsv_id,
        N_raw = nrow(df),
        raw_estimate = estimate,
        raw_std_error = std.error,
        raw_statistic = statistic,
        raw_p = p.value
      )
  }
  
  bind_rows(out)
}

raw_results_list <- vector("list", length(diet_factors))

for (k in seq_along(diet_factors)) {
  factor_name <- diet_factors[k]
  cat(sprintf("[RAW %d/%d] %s\n", k, length(diet_factors), factor_name))
  raw_results_list[[k]] <- run_raw_screen(factor_name, combined_data)
}

raw_diet_dsv_results <- bind_rows(raw_results_list)

write.xlsx(
  raw_diet_dsv_results,
  "results/all_diet_dsv_raw_screening.xlsx",
  rowNames = FALSE
)

raw_candidates <- raw_diet_dsv_results %>%
  filter(!is.na(raw_p), raw_p < 0.05) %>%
  arrange(diet_factor, raw_p)

write.xlsx(
  raw_candidates,
  "results/diet_dsv_raw_p005_candidates.xlsx",
  rowNames = FALSE
)

cat("Diet-dSV pairs passing raw P < 0.05:", nrow(raw_candidates), "\n")

# ----------------------------------------------------------------------------
# 6. Stage 2: adjusted Version A model
# ----------------------------------------------------------------------------

run_adjusted_model <- function(diet_var, dsv_id, data) {
  dsv_species <- str_extract(dsv_id, "^[A-Za-z0-9_]+(?=\\.)")
  if (is.na(dsv_species)) {
    dsv_species <- str_extract(dsv_id, "^[A-Za-z0-9]+_[A-Za-z0-9]+")
  }
  
  dsv_species <- as.character(dsv_species)[1]
  
  if (is.na(dsv_species) || !dsv_species %in% colnames(combined_metadata)) {
    return(NULL)
  }
  
  cov_col_name <- paste0(dsv_species, "_coverage")
  
  df <- data %>%
    filter(dsv == dsv_id) %>%
    mutate(
      target_diet_factor = .data[[diet_var]],
      species_abundance = .data[[dsv_species]],
      species_coverage = if (cov_col_name %in% colnames(.)) .data[[cov_col_name]] else NA_real_
    ) %>%
    drop_na(
      status, target_diet_factor, Age, Sex, BMI, ethnicity, residency,
      medication, urbanization, all_of(general_covars),
      species_abundance, species_coverage
    ) %>%
    filter(status %in% c(0, 1)) %>%
    droplevels()
  
  if (nrow(df) < 30 || length(unique(df$status)) < 2) return(NULL)
  
  cc <- table(df$status)
  if (length(cc) < 2 || any(cc < 10)) return(NULL)
  
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
  
  full_formula <- as.formula(
    paste("status ~", paste(sprintf("`%s`", valid_covars), collapse = " + "))
  )
  
  fit <- tryCatch(
    glm(full_formula, data = df, family = binomial()),
    error = function(e) NULL
  )
  
  if (is.null(fit) || any(!is.finite(coef(fit)))) {
    fit <- tryCatch(
      glm(full_formula, data = df, family = binomial(), method = brglm2::brglmFit),
      error = function(e) NULL
    )
  }
  
  if (is.null(fit)) return(NULL)
  
  tmp <- broom::tidy(fit) %>%
    filter(term == "target_diet_factor")
  
  if (nrow(tmp) == 0) return(NULL)
  
  tmp %>%
    transmute(
      diet_factor = diet_var,
      dsv = dsv_id,
      species = dsv_species,
      N = nrow(df),
      estimate = estimate,
      std.error = std.error,
      statistic = statistic,
      p.value = p.value,
      OR = exp(estimate),
      OR_low = exp(estimate - 1.96 * std.error),
      OR_high = exp(estimate + 1.96 * std.error)
    )
}

cat("\n========== 5. Adjusted model (Version A) ==========\n")

adjusted_A_list <- vector("list", nrow(raw_candidates))

if (nrow(raw_candidates) > 0) {
  pb <- txtProgressBar(min = 0, max = nrow(raw_candidates), style = 3)
  
  for (i in seq_len(nrow(raw_candidates))) {
    adjusted_A_list[[i]] <- run_adjusted_model(
      raw_candidates$diet_factor[i],
      raw_candidates$dsv[i],
      combined_data
    )
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
}

all_diet_factor_df_A <- bind_rows(adjusted_A_list)

if (nrow(all_diet_factor_df_A) > 0) {
  all_diet_factor_df_A <- all_diet_factor_df_A %>%
    left_join(
      raw_candidates %>%
        select(diet_factor, dsv, N_raw, raw_estimate, raw_std_error, raw_statistic, raw_p),
      by = c("diet_factor", "dsv")
    ) %>%
    group_by(diet_factor) %>%
    mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
    ungroup() %>%
    arrange(diet_factor, p_adj)
}

write.xlsx(
  all_diet_factor_df_A,
  "results/all_diet_factors_dsv_adjusted_verA.xlsx",
  rowNames = FALSE
)

sig_factor_dsvs_A <- all_diet_factor_df_A %>%
  filter(!is.na(p_adj), p_adj < 0.05)

write.xlsx(
  sig_factor_dsvs_A,
  "results/significant_diet_factors_dsv_verA.xlsx",
  rowNames = FALSE
)

cat("Version A associations with FDR < 0.05:", nrow(sig_factor_dsvs_A), "\n")

# ----------------------------------------------------------------------------
# 7. Ten-fold cross-validation and marginal R2
# ----------------------------------------------------------------------------

run_cv_for_factor_dsv <- function(diet_var, dsv_name, data_meta, dsv_filtered_data) {
  dsv_df_temp <- dsv_filtered_data %>% rownames_to_column("sample_id")
  dsv_status_df <- dsv_df_temp %>%
    select(sample_id, status_val = all_of(dsv_name))
  
  dsv_species <- str_extract(dsv_name, "^[A-Za-z0-9_]+(?=\\.)")
  if (is.na(dsv_species)) {
    dsv_species <- str_extract(dsv_name, "^[A-Za-z0-9]+_[A-Za-z0-9]+")
  }
  
  dsv_species <- as.character(dsv_species)[1]
  
  if (is.na(dsv_species) || !dsv_species %in% colnames(data_meta)) {
    return(NULL)
  }
  
  cov_col_name <- paste0(dsv_species, "_coverage")
  
  df <- data_meta %>%
    inner_join(dsv_status_df, by = "sample_id") %>%
    mutate(
      status_raw = as.numeric(status_val),
      target_diet_factor = .data[[diet_var]],
      species_abundance = .data[[dsv_species]],
      species_coverage = if (cov_col_name %in% colnames(.)) .data[[cov_col_name]] else NA_real_
    ) %>%
    filter(status_raw %in% c(0, 1)) %>%
    drop_na(
      status_raw, target_diet_factor, Age, Sex, BMI, ethnicity, residency,
      medication, urbanization, all_of(general_covars),
      species_abundance, species_coverage
    ) %>%
    droplevels()
  
  cc <- table(df$status_raw)
  if (nrow(df) < 30 || length(cc) < 2 || any(cc < 5)) return(NULL)
  
  df <- df %>%
    mutate(status = factor(ifelse(status_raw == 1, "Pos", "Neg"), levels = c("Neg", "Pos")))
  
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
  
  full_formula <- as.formula(
    paste("status ~", paste(sprintf("`%s`", valid_covars), collapse = " + "))
  )
  
  full_glm <- tryCatch(
    glm(full_formula, data = df, family = binomial()),
    error = function(e) NULL
  )
  
  if (is.null(full_glm)) return(NULL)
  
  calc_r2 <- function(m) {
    if (is.null(m)) return(NA_real_)
    1 - m$deviance / m$null.deviance
  }
  
  full_r2 <- calc_r2(full_glm)
  
  get_delta_r2 <- function(drop_var) {
    if (!drop_var %in% valid_covars) return(0)
    
    keep_vars <- setdiff(valid_covars, drop_var)
    if (length(keep_vars) == 0) return(full_r2)
    
    reduced_formula <- as.formula(
      paste("status ~", paste(sprintf("`%s`", keep_vars), collapse = " + "))
    )
    
    reduced_model <- tryCatch(
      glm(reduced_formula, data = df, family = binomial()),
      error = function(e) NULL
    )
    
    if (is.null(reduced_model)) return(NA_real_)
    full_r2 - calc_r2(reduced_model)
  }
  
  ctrl <- trainControl(
    method = "cv",
    number = 10,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = TRUE
  )
  
  model_cv <- tryCatch(
    train(
      full_formula,
      data = df,
      method = "glm",
      family = binomial(),
      trControl = ctrl,
      metric = "ROC"
    ),
    error = function(e) NULL
  )
  
  if (is.null(model_cv)) return(NULL)
  
  result <- data.frame(
    diet_factor = diet_var,
    dsv = dsv_name,
    species = dsv_species,
    N = nrow(df),
    Full_Model_R2 = round(full_r2, 4),
    R2_TargetDietFactor = round(get_delta_r2("target_diet_factor"), 4),
    R2_SpeciesCoverage = round(get_delta_r2("species_coverage"), 4),
    R2_SpeciesAbundance = round(get_delta_r2("species_abundance"), 4),
    R2_Ethnicity = round(get_delta_r2("ethnicity"), 4),
    R2_Residency = round(get_delta_r2("residency"), 4),
    R2_Urbanization = round(get_delta_r2("urbanization"), 4),
    R2_Medication = round(get_delta_r2("medication"), 4),
    R2_Age = round(get_delta_r2("Age"), 4),
    R2_Sex = round(get_delta_r2("Sex"), 4),
    R2_BMI = round(get_delta_r2("BMI"), 4),
    CV_Mean_ROC = round(mean(model_cv$results$ROC, na.rm = TRUE), 4),
    CV_Mean_Sens = round(mean(model_cv$results$Sens, na.rm = TRUE), 4),
    CV_Mean_Spec = round(mean(model_cv$results$Spec, na.rm = TRUE), 4),
    CV_ROC_SD = round(sd(model_cv$resample$ROC, na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  )
  
  for (g in general_covars) {
    result[[paste0("R2_", g)]] <- round(get_delta_r2(g), 4)
  }
  
  result
}

if (nrow(sig_factor_dsvs_A) > 0) {
  cat("\n========== 6. Ten-fold CV and marginal R2 ==========\n")
  
  cv_list_A <- vector("list", nrow(sig_factor_dsvs_A))
  pb <- txtProgressBar(min = 0, max = nrow(sig_factor_dsvs_A), style = 3)
  
  for (i in seq_len(nrow(sig_factor_dsvs_A))) {
    cv_list_A[[i]] <- run_cv_for_factor_dsv(
      sig_factor_dsvs_A$diet_factor[i],
      sig_factor_dsvs_A$dsv[i],
      combined_metadata,
      dsv_filtered
    )
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  cv_df_A <- bind_rows(cv_list_A)
  
  if (nrow(cv_df_A) > 0) {
    write.xlsx(
      cv_df_A,
      "results/cv_10fold_diet_factors_verA_coverage_as_covariate.xlsx",
      rowNames = FALSE
    )
  }
}

# ----------------------------------------------------------------------------
# 8. Summary of dSV counts for each dietary factor
# ----------------------------------------------------------------------------

diet_dsv_count_summary <- data.frame(
  Diet = diet_factors,
  Raw_P_lt_0.05_dSV_N = sapply(diet_factors, function(x) {
    raw_candidates %>%
      filter(diet_factor == x) %>%
      distinct(dsv) %>%
      nrow()
  }),
  Adjusted_Significant_dSV_N = sapply(diet_factors, function(x) {
    sig_factor_dsvs_A %>%
      filter(diet_factor == x) %>%
      distinct(dsv) %>%
      nrow()
  })
) %>%
  arrange(desc(Raw_P_lt_0.05_dSV_N))

write.csv(
  diet_dsv_count_summary,
  "results/diet_dsv_raw_and_adjusted_counts_verA.csv",
  row.names = FALSE,
  quote = FALSE
)

cat("\n========== Analysis completed ==========\n")
cat("Raw P < 0.05 candidate pairs:", nrow(raw_candidates), "\n")
cat("Version A associations with FDR < 0.05:", nrow(sig_factor_dsvs_A), "\n")
cat("Results saved in the 'results' directory.\n")
