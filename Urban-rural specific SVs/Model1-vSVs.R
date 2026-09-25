# ============================================================================
# Residency-associated vSV analysis
# Version A: coverage included as a covariate without a 5X cutoff
# Version B: sensitivity analysis restricted to species coverage >= 5X
# ============================================================================

library(tidyverse)
library(broom)
library(openxlsx)
library(caret)

set.seed(12345)

# ============================================================================
# SHARED DATA PREPARATION
# ============================================================================

data_dir <- "."
results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir)

cat(sprintf("[%s] Reading input files...\n", Sys.time()))
vsv_data <- read.csv(file.path(data_dir, "final_filtered_vsv.csv"), row.names = 1, check.names = FALSE)
metadata <- read.csv("metadata.csv", check.names = FALSE)
diet_data <- read.csv("final_filtered_diet.csv", check.names = FALSE)
medication_data <- read.csv("Medication.csv", check.names = FALSE)
urbanization_data <- read.csv("Urbanization.csv", check.names = FALSE)
general_metadata_data <- read.csv("General metadata.csv", check.names = FALSE)

required_cols <- c("sample_id", "residency", "Age", "Sex", "BMI", "ethnicity")
missing_cols <- setdiff(required_cols, names(metadata))
if (length(missing_cols) > 0) stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))

# ----------------------------------------------------------------------------
# Composite covariates: diet, medication, and urbanization
# ----------------------------------------------------------------------------
process_covariate <- function(data, covariate_name) {
  vars_df <- data %>% select(-sample_id)
  numeric_df <- vars_df %>%
    mutate(across(everything(), function(col) {
      if (is.numeric(col)) return(col)
      col_char <- trimws(as.character(col))
      col_lower <- tolower(col_char)
      if (all(col_lower %in% c("yes", "no", "y", "n", "true", "false", "1", "0", "", NA))) {
        return(case_when(
          col_lower %in% c("yes", "y", "true", "1") ~ 1,
          col_lower %in% c("no", "n", "false", "0") ~ 0,
          TRUE ~ NA_real_
        ))
      }
      num_col <- suppressWarnings(as.numeric(col_char))
      if (all(is.na(num_col[!is.na(col_char) & col_char != ""]))) num_col <- as.numeric(as.factor(col_char))
      num_col
    }))
  
  data %>%
    select(sample_id) %>%
    mutate(!!covariate_name := rowSums(numeric_df, na.rm = TRUE))
}

# Reverse-code Contact with Animals before calculating the Urbanization composite.
reverse_ordinal <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  if (all(is.na(x_num))) return(x)
  rng <- range(x_num, na.rm = TRUE)
  rng[1] + rng[2] - x_num
}

animal_col <- names(urbanization_data)[tolower(trimws(names(urbanization_data))) == "contact with animals"]
if (length(animal_col) == 1) {
  urbanization_data[[animal_col]] <- reverse_ordinal(urbanization_data[[animal_col]])
}

diet_cov <- process_covariate(diet_data, "diet")
medication_cov <- process_covariate(medication_data, "medication")
urbanization_cov <- process_covariate(urbanization_data, "urbanization")

# ----------------------------------------------------------------------------
# General metadata: retain each item as an individual model covariate
# ----------------------------------------------------------------------------
general_metadata_cov <- general_metadata_data %>%
  mutate(across(-sample_id, function(col) {
    if (is.numeric(col)) return(col)
    col_char <- trimws(as.character(col))
    col_lower <- tolower(col_char)
    
    if (all(col_lower %in% c("yes", "no", "y", "n", "true", "false", "1", "0", "", NA))) {
      return(case_when(
        col_lower %in% c("yes", "y", "true", "1") ~ 1,
        col_lower %in% c("no", "n", "false", "0") ~ 0,
        TRUE ~ NA_real_
      ))
    }
    
    num_col <- suppressWarnings(as.numeric(col_char))
    nonmissing_original <- !is.na(col_char) & col_char != ""
    if (sum(!is.na(num_col)) == sum(nonmissing_original)) return(num_col)
    factor(col_char)
  }))

general_cols_original <- setdiff(names(general_metadata_cov), "sample_id")
general_covars <- paste0("general__", make.names(general_cols_original, unique = TRUE))
names(general_metadata_cov)[match(general_cols_original, names(general_metadata_cov))] <- general_covars

# ----------------------------------------------------------------------------
# Species abundance and coverage matrices
# ----------------------------------------------------------------------------
cat(sprintf("[%s] Reading species abundance and coverage data...\n", Sys.time()))
species_map <- read.csv("sv_to_metaphlan_final_mapping.csv", check.names = FALSE)
metaphlan_raw <- read_tsv("yunnan_7_years_ago_combined_metaphlan_profile.tsv", comment = "", show_col_types = FALSE)
if (grepl("#", colnames(metaphlan_raw)[1])) colnames(metaphlan_raw)[1] <- "clade_name"

abundance_matrix <- metaphlan_raw %>%
  filter(clade_name %in% unique(species_map$metaphlan_clade)) %>%
  left_join(species_map, by = c("clade_name" = "metaphlan_clade")) %>%
  select(-clade_name, -metaphlan_species) %>%
  pivot_longer(cols = -sv_species, names_to = "sample_id", values_to = "abundance") %>%
  mutate(abundance = as.numeric(abundance)) %>%
  group_by(sample_id, sv_species) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(id_cols = sample_id, names_from = sv_species, values_from = abundance, values_fill = 0)

raw_cov_data <- read.csv("merged_species_mean_coverage_matrix.csv", check.names = FALSE)
colnames(raw_cov_data)[1] <- "sample_id"
coverage_raw <- raw_cov_data %>%
  pivot_longer(cols = -sample_id, names_to = "sv_species", values_to = "coverage") %>%
  mutate(coverage = as.numeric(coverage)) %>%
  group_by(sample_id, sv_species) %>%
  summarise(coverage = mean(coverage, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(id_cols = sample_id, names_from = sv_species, values_from = coverage, values_fill = 0)
colnames(coverage_raw)[-1] <- paste0(colnames(coverage_raw)[-1], "_coverage")

combined_metadata <- metadata %>%
  inner_join(diet_cov, by = "sample_id") %>%
  inner_join(medication_cov, by = "sample_id") %>%
  inner_join(urbanization_cov, by = "sample_id") %>%
  inner_join(general_metadata_cov, by = "sample_id") %>%
  left_join(abundance_matrix, by = "sample_id") %>%
  left_join(coverage_raw, by = "sample_id") %>%
  mutate(Sex = as.factor(Sex), residency = as.factor(residency), ethnicity = as.factor(ethnicity))

cat("Samples after metadata merge:", nrow(combined_metadata), "\n")

# ----------------------------------------------------------------------------
# vSV preprocessing: 0 -> NA -> z-score; retain prevalence >= 10%
# ----------------------------------------------------------------------------
vsv_na <- vsv_data
vsv_na[vsv_na == 0] <- NA
vsv_scaled <- as.data.frame(scale(vsv_na, center = TRUE, scale = TRUE))
vsv_prevalence <- colMeans(vsv_data != 0, na.rm = TRUE)
vsv_filtered <- vsv_scaled[, names(vsv_prevalence[vsv_prevalence >= 0.1]), drop = FALSE]

combined_data <- vsv_filtered %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(cols = -sample_id, names_to = "vsv", values_to = "expression") %>%
  mutate(expression = as.numeric(expression)) %>%
  inner_join(combined_metadata, by = "sample_id")

extract_species <- function(vsv_name) {
  sp <- str_extract(vsv_name, "^[A-Za-z0-9.]+\\.[a-zA-Z0-9_.]+|^[A-Za-z0-9]+_[a-zA-Z0-9_]+")
  if (is.na(sp)) sp <- str_split(vsv_name, "_")[[1]][1]
  sp
}

quote_vars <- function(x) paste0("`", gsub("`", "", x), "`")
make_formula <- function(response, vars) as.formula(paste(response, "~", paste(quote_vars(vars), collapse = " + ")))

# ============================================================================
# VERSION A
# Coverage is included as a covariate; no 5X cutoff is applied.
# ============================================================================

cat("\n================ VERSION A ================\n")

screen_vsv_residency_verA <- function(data, residency_name) {
  data <- data %>% mutate(residency_binary = ifelse(residency == residency_name, 1, 0))
  results <- list()
  vsv_names <- unique(data$vsv)
  
  for (i in seq_along(vsv_names)) {
    vsv_item <- vsv_names[i]
    vsv_species <- extract_species(vsv_item)
    cov_col_name <- paste0(vsv_species, "_coverage")
    cov_vec <- if (!is.na(cov_col_name) && cov_col_name %in% colnames(data)) data[[cov_col_name]] else 0
    
    current_data <- data %>%
      filter(vsv == !!vsv_item) %>%
      mutate(species_coverage = cov_vec[data$vsv == !!vsv_item]) %>%
      drop_na(expression, residency_binary, species_coverage)
    
    if (nrow(current_data) < 30) {
      results[[vsv_item]] <- tibble(vsv = vsv_item, p.value = NA_real_, rho = NA_real_, N_cov_available = nrow(current_data))
      next
    }
    
    cor_test <- cor.test(current_data$expression, current_data$residency_binary, method = "spearman", exact = FALSE)
    results[[vsv_item]] <- tibble(vsv = vsv_item, p.value = cor_test$p.value, rho = unname(cor_test$estimate), N_cov_available = nrow(current_data))
  }
  
  bind_rows(results)
}

analyze_residency_verA <- function(residency_name, data) {
  data <- data %>% mutate(residency_binary = ifelse(residency == residency_name, 1, 0))
  results <- list()
  vsv_names <- unique(data$vsv)
  
  for (i in seq_along(vsv_names)) {
    vsv_item <- vsv_names[i]
    vsv_species <- extract_species(vsv_item)
    cov_col_name <- paste0(vsv_species, "_coverage")
    abun_vec <- if (!is.na(vsv_species) && vsv_species %in% colnames(data)) data[[vsv_species]] else 0
    cov_vec <- if (!is.na(cov_col_name) && cov_col_name %in% colnames(data)) data[[cov_col_name]] else abun_vec
    
    current_data <- data %>%
      filter(vsv == !!vsv_item) %>%
      mutate(species_abundance = abun_vec[data$vsv == !!vsv_item], species_coverage = cov_vec[data$vsv == !!vsv_item]) %>%
      drop_na(expression, residency_binary, Age, Sex, BMI, ethnicity, diet, medication, urbanization,
              all_of(general_covars), species_abundance, species_coverage)
    
    if (nrow(current_data) < 10) next
    
    model_vars <- c("residency_binary", "species_coverage", "species_abundance", "Age", "Sex", "BMI", "ethnicity",
                    "diet", "medication", "urbanization", general_covars)
    model <- tryCatch(lm(make_formula("expression", model_vars), data = current_data), error = function(e) NULL)
    if (is.null(model)) next
    
    results[[vsv_item]] <- tidy(model, conf.int = TRUE) %>%
      filter(term == "residency_binary") %>%
      mutate(vsv = vsv_item, species = ifelse(is.na(vsv_species), "Unknown", vsv_species), N_samples = nrow(current_data)) %>%
      select(vsv, species, N_samples, estimate, std.error, statistic, p.value, conf.low, conf.high)
  }
  
  if (length(results) == 0) return(NULL)
  bind_rows(results) %>% mutate(residency = residency_name) %>% select(residency, vsv, species, N_samples, everything())
}

all_results_verA <- list()
for (res in c("Urban", "Rural")) {
  if (!(res %in% levels(combined_metadata$residency))) next
  cat(sprintf("[%s] Version A: analyzing %s vs others\n", Sys.time(), res))
  screen_res <- screen_vsv_residency_verA(combined_data, res)
  passed_vsv <- screen_res %>% filter(!is.na(p.value), p.value < 0.05) %>% pull(vsv)
  if (length(passed_vsv) == 0) next
  
  current_result <- analyze_residency_verA(res, combined_data %>% filter(vsv %in% passed_vsv))
  if (!is.null(current_result) && nrow(current_result) > 0) {
    current_result <- current_result %>% mutate(p_adj = p.adjust(p.value, method = "BH"))
    all_results_verA[[res]] <- current_result
    write.xlsx(current_result, file.path(results_dir, paste0(res, "_vs_Others_vsv_results_verA_coverage_as_covariate.xlsx")), rowNames = FALSE)
  }
}

final_results_verA <- bind_rows(all_results_verA)
significant_results_verA <- final_results_verA %>% filter(!is.na(p_adj), p_adj < 0.05) %>% arrange(p_adj)
write.xlsx(final_results_verA, file.path(results_dir, "all_residencies_vsv_results_verA_coverage_as_covariate.xlsx"), rowNames = FALSE)
write.xlsx(significant_results_verA, file.path(results_dir, "all_residencies_significant_verA_coverage_as_covariate.xlsx"), rowNames = FALSE)

run_cv_verA <- function(vsv_name, data_meta, vsv_df) {
  vsv_expr_df <- vsv_df %>% select(sample_id, expression = all_of(vsv_name))
  vsv_species <- extract_species(vsv_name)
  cov_col_name <- paste0(vsv_species, "_coverage")
  
  df <- data_meta %>% inner_join(vsv_expr_df, by = "sample_id") %>% mutate(expression = as.numeric(expression))
  df$species_abundance <- if (!is.na(vsv_species) && vsv_species %in% colnames(df)) df[[vsv_species]] else 0
  df$species_coverage <- if (!is.na(cov_col_name) && cov_col_name %in% colnames(df)) df[[cov_col_name]] else df$species_abundance
  
  df <- df %>% drop_na(expression, residency, Age, Sex, BMI, ethnicity, diet, medication, urbanization,
                       all_of(general_covars), species_abundance, species_coverage)
  if (nrow(df) < 20 || length(unique(df$expression)) < 2) return(NULL)
  
  covars <- c("ethnicity", "residency", "Age", "Sex", "BMI", "diet", "medication", "urbanization",
              general_covars, "species_abundance", "species_coverage")
  full_lm <- tryCatch(lm(make_formula("expression", covars), data = df), error = function(e) NULL)
  if (is.null(full_lm)) return(NULL)
  full_r2 <- summary(full_lm)$r.squared
  
  get_delta_r2 <- function(drop_var) {
    keep_vars <- setdiff(covars, drop_var)
    fit_red <- tryCatch(lm(make_formula("expression", keep_vars), data = df), error = function(e) NULL)
    if (is.null(fit_red)) return(NA_real_)
    max(0, full_r2 - summary(fit_red)$r.squared)
  }
  
  r2_vals <- setNames(lapply(covars, get_delta_r2), covars)
  model_data <- df %>% select(expression, all_of(covars))
  ctrl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)
  model_lm <- tryCatch(train(expression ~ ., data = model_data, method = "lm", trControl = ctrl, metric = "Rsquared"), error = function(e) NULL)
  if (is.null(model_lm)) return(NULL)
  
  res_row <- data.frame(
    vsv = vsv_name, species = vsv_species, N = nrow(model_data), Full_Model_R2 = round(full_r2, 4),
    CV_Mean_Rsquared = round(model_lm$results$Rsquared, 4),
    CV_Rsquared_SD = round(sd(model_lm$resample$Rsquared, na.rm = TRUE), 4),
    CV_Mean_RMSE = round(model_lm$results$RMSE, 4),
    CV_RMSE_SD = round(sd(model_lm$resample$RMSE, na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  )
  
  for (nm in names(r2_vals)) res_row[[paste0("R2_", make.names(nm))]] <- round(as.numeric(r2_vals[[nm]]), 4)
  res_row
}

if (nrow(significant_results_verA) > 0) {
  vsv_scaled_df <- vsv_scaled %>% rownames_to_column("sample_id")
  cv_verA <- map_df(unique(significant_results_verA$vsv), ~ run_cv_verA(.x, combined_metadata, vsv_scaled_df))
  if (nrow(cv_verA) > 0) write.xlsx(cv_verA, file.path(results_dir, "vsv_verA_cv_marginal_r2.xlsx"), rowNames = FALSE)
}

# ============================================================================
# VERSION B
# Sensitivity analysis restricted to species coverage >= 5X.
# Coverage remains included as a covariate.
# ============================================================================

cat("\n================ VERSION B ================\n")

screen_vsv_residency_verB <- function(data, residency_name) {
  data <- data %>% mutate(residency_binary = ifelse(residency == residency_name, 1, 0))
  results <- list()
  vsv_names <- unique(data$vsv)
  
  for (i in seq_along(vsv_names)) {
    vsv_item <- vsv_names[i]
    vsv_species <- extract_species(vsv_item)
    cov_col_name <- paste0(vsv_species, "_coverage")
    cov_vec <- if (!is.na(cov_col_name) && cov_col_name %in% colnames(data)) data[[cov_col_name]] else 0
    
    current_data <- data %>%
      filter(vsv == !!vsv_item) %>%
      mutate(species_coverage = cov_vec[data$vsv == !!vsv_item]) %>%
      drop_na(expression, residency_binary, species_coverage) %>%
      filter(species_coverage >= 5)
    
    if (nrow(current_data) < 30) {
      results[[vsv_item]] <- tibble(vsv = vsv_item, p.value = NA_real_, rho = NA_real_, N_cov5x = nrow(current_data))
      next
    }
    
    cor_test <- cor.test(current_data$expression, current_data$residency_binary, method = "spearman", exact = FALSE)
    results[[vsv_item]] <- tibble(vsv = vsv_item, p.value = cor_test$p.value, rho = unname(cor_test$estimate), N_cov5x = nrow(current_data))
  }
  
  bind_rows(results)
}

analyze_residency_verB <- function(residency_name, data) {
  data <- data %>% mutate(residency_binary = ifelse(residency == residency_name, 1, 0))
  results <- list()
  vsv_names <- unique(data$vsv)
  
  for (i in seq_along(vsv_names)) {
    vsv_item <- vsv_names[i]
    vsv_species <- extract_species(vsv_item)
    cov_col_name <- paste0(vsv_species, "_coverage")
    abun_vec <- if (!is.na(vsv_species) && vsv_species %in% colnames(data)) data[[vsv_species]] else 0
    cov_vec <- if (!is.na(cov_col_name) && cov_col_name %in% colnames(data)) data[[cov_col_name]] else abun_vec
    
    current_data <- data %>%
      filter(vsv == !!vsv_item) %>%
      mutate(species_abundance = abun_vec[data$vsv == !!vsv_item], species_coverage = cov_vec[data$vsv == !!vsv_item]) %>%
      drop_na(expression, residency_binary, Age, Sex, BMI, ethnicity, diet, medication, urbanization,
              all_of(general_covars), species_abundance, species_coverage) %>%
      filter(species_coverage >= 5)
    
    if (nrow(current_data) < 10) next
    
    model_vars <- c("residency_binary", "species_coverage", "species_abundance", "Age", "Sex", "BMI", "ethnicity",
                    "diet", "medication", "urbanization", general_covars)
    model <- tryCatch(lm(make_formula("expression", model_vars), data = current_data), error = function(e) NULL)
    if (is.null(model)) next
    
    results[[vsv_item]] <- tidy(model, conf.int = TRUE) %>%
      filter(term == "residency_binary") %>%
      mutate(vsv = vsv_item, species = ifelse(is.na(vsv_species), "Unknown", vsv_species), N_samples_cov5x = nrow(current_data)) %>%
      select(vsv, species, N_samples_cov5x, estimate, std.error, statistic, p.value, conf.low, conf.high)
  }
  
  if (length(results) == 0) return(NULL)
  bind_rows(results) %>% mutate(residency = residency_name) %>% select(residency, vsv, species, N_samples_cov5x, everything())
}

all_results_verB <- list()
for (res in c("Urban", "Rural")) {
  if (!(res %in% levels(combined_metadata$residency))) next
  cat(sprintf("[%s] Version B: analyzing %s vs others\n", Sys.time(), res))
  screen_res <- screen_vsv_residency_verB(combined_data, res)
  passed_vsv <- screen_res %>% filter(!is.na(p.value), p.value < 0.05) %>% pull(vsv)
  if (length(passed_vsv) == 0) next
  
  current_result <- analyze_residency_verB(res, combined_data %>% filter(vsv %in% passed_vsv))
  if (!is.null(current_result) && nrow(current_result) > 0) {
    current_result <- current_result %>% mutate(p_adj = p.adjust(p.value, method = "BH"))
    all_results_verB[[res]] <- current_result
    write.xlsx(current_result, file.path(results_dir, paste0(res, "_vs_Others_vsv_results_verB_cov5x.xlsx")), rowNames = FALSE)
  }
}

final_results_verB <- bind_rows(all_results_verB)
significant_results_verB <- final_results_verB %>% filter(!is.na(p_adj), p_adj < 0.05) %>% arrange(p_adj)
write.xlsx(final_results_verB, file.path(results_dir, "all_residencies_vsv_results_verB_cov5x.xlsx"), rowNames = FALSE)
write.xlsx(significant_results_verB, file.path(results_dir, "all_residencies_significant_verB_cov5x.xlsx"), rowNames = FALSE)

run_cv_verB <- function(vsv_name, data_meta, vsv_df) {
  vsv_expr_df <- vsv_df %>% select(sample_id, expression = all_of(vsv_name))
  vsv_species <- extract_species(vsv_name)
  cov_col_name <- paste0(vsv_species, "_coverage")
  
  df <- data_meta %>% inner_join(vsv_expr_df, by = "sample_id") %>% mutate(expression = as.numeric(expression))
  df$species_abundance <- if (!is.na(vsv_species) && vsv_species %in% colnames(df)) df[[vsv_species]] else 0
  df$species_coverage <- if (!is.na(cov_col_name) && cov_col_name %in% colnames(df)) df[[cov_col_name]] else df$species_abundance
  
  df <- df %>%
    drop_na(expression, residency, Age, Sex, BMI, ethnicity, diet, medication, urbanization,
            all_of(general_covars), species_abundance, species_coverage) %>%
    filter(species_coverage >= 5)
  if (nrow(df) < 20 || length(unique(df$expression)) < 2) return(NULL)
  
  covars <- c("ethnicity", "residency", "Age", "Sex", "BMI", "diet", "medication", "urbanization",
              general_covars, "species_abundance", "species_coverage")
  full_lm <- tryCatch(lm(make_formula("expression", covars), data = df), error = function(e) NULL)
  if (is.null(full_lm)) return(NULL)
  full_r2 <- summary(full_lm)$r.squared
  
  get_delta_r2 <- function(drop_var) {
    keep_vars <- setdiff(covars, drop_var)
    fit_red <- tryCatch(lm(make_formula("expression", keep_vars), data = df), error = function(e) NULL)
    if (is.null(fit_red)) return(NA_real_)
    max(0, full_r2 - summary(fit_red)$r.squared)
  }
  
  r2_vals <- setNames(lapply(covars, get_delta_r2), covars)
  model_data <- df %>% select(expression, all_of(covars))
  ctrl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)
  model_lm <- tryCatch(train(expression ~ ., data = model_data, method = "lm", trControl = ctrl, metric = "Rsquared"), error = function(e) NULL)
  if (is.null(model_lm)) return(NULL)
  
  res_row <- data.frame(
    vsv = vsv_name, species = vsv_species, N_cov5x = nrow(model_data), Full_Model_R2 = round(full_r2, 4),
    CV_Mean_Rsquared = round(model_lm$results$Rsquared, 4),
    CV_Rsquared_SD = round(sd(model_lm$resample$Rsquared, na.rm = TRUE), 4),
    CV_Mean_RMSE = round(model_lm$results$RMSE, 4),
    CV_RMSE_SD = round(sd(model_lm$resample$RMSE, na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  )
  
  for (nm in names(r2_vals)) res_row[[paste0("R2_", make.names(nm))]] <- round(as.numeric(r2_vals[[nm]]), 4)
  res_row
}

if (nrow(significant_results_verB) > 0) {
  vsv_scaled_df <- vsv_scaled %>% rownames_to_column("sample_id")
  cv_verB <- map_df(unique(significant_results_verB$vsv), ~ run_cv_verB(.x, combined_metadata, vsv_scaled_df))
  if (nrow(cv_verB) > 0) write.xlsx(cv_verB, file.path(results_dir, "vsv_verB_cov5x_cv_marginal_r2.xlsx"), rowNames = FALSE)
}

cat("\nResidency-associated vSV Version A and Version B analyses completed.\n")
