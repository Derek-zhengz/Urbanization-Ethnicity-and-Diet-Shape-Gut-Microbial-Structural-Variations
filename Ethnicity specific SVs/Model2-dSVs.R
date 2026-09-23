# ============================================================================
# Reproducible dSV ethnicity association models
# Version A: sequencing coverage included as a covariate (no coverage cutoff)
# Version B: sensitivity analysis restricted to sequencing coverage >= 5X
#
# Statistical procedures, thresholds, covariates, model formulas, FDR correction,
# cross-validation, and R2 calculations are preserved from the original analysis.
# Terminology is standardized to "Urbanization" throughout.
# ============================================================================

library(tidyverse)
library(broom)
library(openxlsx)
library(lmtest)
library(caret)

set.seed(12345)

# ----------------------------------------------------------------------------
# 1. Data paths and input files
# ----------------------------------------------------------------------------
data_dir <- "."

cat(sprintf("[%s] Reading input files...\n", Sys.time()))
dsv_data <- read.csv(file.path(data_dir, "final_filtered_dsv.csv"), row.names = 1)
metadata <- read.csv("metadata.csv")
diet_data <- read.csv("final_filtered_diet.csv", check.names = FALSE)
medication_data <- read.csv("Medication.csv", check.names = FALSE)
urbanization_data <- read.csv("Urbanization.csv", check.names = FALSE)
general_metadata_data <- read.csv("General metadata.csv", check.names = FALSE)

# 2. Check required metadata columns
required_cols <- c("sample_id", "residency", "Age", "Sex", "BMI", "ethnicity")
missing_cols <- setdiff(required_cols, names(metadata))
if (length(missing_cols) > 0) stop("Missing required columns in metadata: ", paste(missing_cols, collapse = ", "))

# 3. Construct aggregate covariates
process_covariate <- function(data, covariate_name) {
  vars_df <- data %>% select(-sample_id)
  numeric_df <- vars_df %>%
    mutate(across(everything(), function(col) {
      if (is.numeric(col)) return(col)
      col_char <- as.character(col) %>% trimws()
      col_lower <- tolower(col_char)
      if (all(col_lower %in% c("yes", "no", "y", "n", "true", "false", "1", "0", "", NA))) {
        return(case_when(
          col_lower %in% c("yes", "y", "true", "1") ~ 1,
          col_lower %in% c("no", "n", "false", "0") ~ 0,
          TRUE ~ NA_real_
        ))
      }
      num_col <- suppressWarnings(as.numeric(col_char))
      if (all(is.na(num_col[!is.na(col_char) & col_char != ""]))) {
        num_col <- as.numeric(as.factor(col_char))
      }
      return(num_col)
    }))
  
  result <- data %>% select(sample_id) %>% mutate(!!covariate_name := rowSums(numeric_df, na.rm = TRUE))
  return(result)
}

cat(sprintf("[%s] Constructing aggregate covariates...\n", Sys.time()))
diet_total <- process_covariate(diet_data, "diet")
medication_total <- process_covariate(medication_data, "medication")
urbanization_total <- process_covariate(urbanization_data, "urbanization")
general_metadata_total <- process_covariate(general_metadata_data, "general_metadata")

# 4. Read species abundance and sequencing coverage matrices
cat(sprintf("[%s] Reading species abundance and sequencing coverage data...\n", Sys.time()))
species_map <- read.csv("sv_to_metaphlan_final_mapping.csv")
metaphlan_raw <- read_tsv("yunnan_7_years_ago_combined_metaphlan_profile.tsv", comment = "", show_col_types = FALSE)
if (grepl("#", colnames(metaphlan_raw)[1])) colnames(metaphlan_raw)[1] <- "clade_name"

abundance_matrix <- metaphlan_raw %>%
  filter(clade_name %in% unique(species_map$metaphlan_clade)) %>%
  left_join(species_map, by = c("clade_name" = "metaphlan_clade")) %>%
  select(-clade_name, -metaphlan_species) %>%
  pivot_longer(cols = -sv_species, names_to = "sample_id", values_to = "abundance") %>%
  mutate(abundance = as.numeric(abundance)) %>%
  pivot_wider(id_cols = sample_id, names_from = sv_species, values_from = abundance, values_fill = list(abundance = 0))

raw_cov_data <- read.csv("merged_species_mean_coverage_matrix.csv", check.names = FALSE)
colnames(raw_cov_data)[1] <- "sample_id"
coverage_raw <- raw_cov_data %>%
  pivot_longer(cols = -sample_id, names_to = "sv_species", values_to = "coverage") %>%
  mutate(coverage = as.numeric(coverage)) %>%
  group_by(sample_id, sv_species) %>% summarise(coverage = mean(coverage, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(id_cols = sample_id, names_from = sv_species, values_from = coverage, values_fill = list(coverage = 0))
colnames(coverage_raw)[-1] <- paste0(colnames(coverage_raw)[-1], "_coverage")

# 5. Merge metadata and covariates
combined_metadata <- metadata %>%
  inner_join(diet_total, by = "sample_id") %>% inner_join(medication_total, by = "sample_id") %>%
  inner_join(urbanization_total, by = "sample_id") %>% inner_join(general_metadata_total, by = "sample_id") %>%
  left_join(abundance_matrix, by = "sample_id") %>% left_join(coverage_raw, by = "sample_id") %>%
  mutate(Sex = as.factor(Sex), residency = as.factor(residency), ethnicity = as.factor(ethnicity))

cat("\n================ Cohort summary ================\n")
cat("Sample counts by ethnicity:\n")
print(table(combined_metadata$ethnicity))
cat("==============================================\n\n")

# 6. Preprocess dSV data (retain NA values)
cat(sprintf("[%s] Preprocessing dSV data (retaining NA values)...\n", Sys.time()))
dsv_clean <- dsv_data %>%
  mutate(across(everything(), as.numeric))

# Filter low-prevalence dSVs (>=10%; NA values excluded from prevalence calculation)
dsv_prevalence <- colMeans(dsv_clean != 0, na.rm = TRUE)
dsv_filtered <- dsv_clean[, names(dsv_prevalence[dsv_prevalence >= 0.1]), drop = FALSE]

# Convert to long format while retaining NA values
combined_data <- dsv_filtered %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(cols = -sample_id, names_to = "dsv", values_to = "status") %>%
  mutate(status = as.numeric(status)) %>%  # retain NA values
  inner_join(combined_metadata, by = "sample_id")

# ----------------------------------------------------------------------------
# 7. Screening function: no coverage cutoff, chi-squared test (p < 0.05)
# ----------------------------------------------------------------------------
screen_dsv_ethnicity_no_cov_filter <- function(data, target_eth) {
  data <- data %>% mutate(ethnicity_binary = ifelse(ethnicity == target_eth, 1, 0))
  results <- list()
  dsv_names <- unique(data$dsv)
  total_dsv <- length(dsv_names)
  
  cat(sprintf("  -> Chi-squared screening (no coverage cutoff): %d features\n", total_dsv))
  pb <- txtProgressBar(min = 0, max = total_dsv, style = 3)
  
  for (i in seq_along(dsv_names)) {
    setTxtProgressBar(pb, i)
    dsv_item <- dsv_names[i]
    
    # Extract species name for coverage matching
    dsv_species <- str_extract(dsv_item, "^[A-Za-z0-9.]+\\.[a-zA-Z0-9_.]+|^[A-Za-z0-9]+_[a-zA-Z0-9_]+")
    if (is.na(dsv_species)) dsv_species <- str_split(dsv_item, "_")[[1]][1]
    cov_col_name <- paste0(dsv_species, "_coverage")
    cov_vec <- if (!is.na(cov_col_name) && cov_col_name %in% colnames(data)) data[[cov_col_name]] else 0
    
    current_data <- data %>%
      filter(dsv == !!dsv_item) %>%
      mutate(species_coverage = cov_vec[data$dsv == !!dsv_item]) %>%
      drop_na(status, ethnicity_binary, species_coverage)  # require available coverage data
    
    if (nrow(current_data) < 30) {
      results[[dsv_item]] <- tibble(dsv = dsv_item, p.value = NA, N_cov_available = nrow(current_data), note = "Insufficient sample size (<30)")
      next
    }
    
    tbl <- table(current_data$status, current_data$ethnicity_binary)
    if (any(rowSums(tbl) == 0) || any(colSums(tbl) == 0)) {
      results[[dsv_item]] <- tibble(dsv = dsv_item, p.value = NA, N_cov_available = nrow(current_data), note = "Invalid contingency table")
      next
    }
    
    chi_test <- tryCatch({ chisq.test(tbl) }, error = function(e) NULL)
    if (is.null(chi_test)) {
      fisher_test <- fisher.test(tbl)
      p_val <- fisher_test$p.value
      note <- "Fisher exact test"
    } else {
      if (any(chi_test$expected < 5)) {
        fisher_test <- fisher.test(tbl)
        p_val <- fisher_test$p.value
        note <- "Fisher exact test (expected count <5)"
      } else {
        p_val <- chi_test$p.value
        note <- "Chi-squared test"
      }
    }
    
    results[[dsv_item]] <- tibble(
      dsv = dsv_item,
      p.value = p_val,
      N_cov_available = nrow(current_data),
      note = note
    )
  }
  close(pb)
  bind_rows(results)
}

# ----------------------------------------------------------------------------
# 8. Logistic regression: coverage as a covariate, no coverage cutoff
# ----------------------------------------------------------------------------
analyze_dsv_ethnicity_no_cov_filter <- function(target_eth, data) {
  data <- data %>% mutate(ethnicity_binary = ifelse(ethnicity == target_eth, 1, 0))
  results <- list()
  dsv_names <- unique(data$dsv)
  total_dsv <- length(dsv_names)
  
  cat(sprintf("\n  -> Multivariable logistic regression (coverage as a covariate, no cutoff): %d candidate features\n", total_dsv))
  pb <- txtProgressBar(min = 0, max = total_dsv, style = 3)
  
  for (i in seq_along(dsv_names)) {
    setTxtProgressBar(pb, i)
    dsv_item <- dsv_names[i]
    dsv_species <- str_extract(dsv_item, "^[A-Za-z0-9.]+\\.[a-zA-Z0-9_.]+|^[A-Za-z0-9]+_[a-zA-Z0-9_]+")
    if (is.na(dsv_species)) dsv_species <- str_split(dsv_item, "_")[[1]][1]
    
    cov_col_name <- paste0(dsv_species, "_coverage")
    abun_vec <- if (!is.na(dsv_species) && dsv_species %in% colnames(data)) data[[dsv_species]] else 0
    cov_vec  <- if (!is.na(cov_col_name) && cov_col_name %in% colnames(data)) data[[cov_col_name]] else abun_vec
    
    current_data <- data %>%
      filter(dsv == !!dsv_item) %>%
      mutate(species_abundance = abun_vec[data$dsv == !!dsv_item], 
             species_coverage = cov_vec[data$dsv == !!dsv_item]) %>%
      drop_na(status, ethnicity_binary, Age, Sex, BMI, residency,
              diet, medication, urbanization, general_metadata,
              species_abundance, species_coverage)
    
    if (nrow(current_data) < 10) next
    
    # Check outcome class balance
    status_table <- table(current_data$status)
    if (length(status_table) < 2 || any(status_table < 3)) next
    
    model_formula <- as.formula(
      "status ~ ethnicity_binary + species_coverage + species_abundance + Age + Sex + BMI + residency +
       diet + medication + urbanization + general_metadata"
    )
    
    model <- tryCatch({ glm(model_formula, data = current_data, family = binomial()) }, error = function(e) NULL)
    if (is.null(model)) next
    
    # Fall back to Wald confidence intervals if profile-likelihood CI fails
    tidy_res <- tryCatch({
      tidy(model, conf.int = TRUE)
    }, error = function(e) {
      tidy(model, conf.int = FALSE) %>%
        mutate(
          conf.low = estimate - 1.96 * std.error,
          conf.high = estimate + 1.96 * std.error
        )
    })
    if (is.null(tidy_res) || nrow(tidy_res) == 0) next
    
    results[[dsv_item]] <- tidy_res %>%
      filter(term == "ethnicity_binary") %>%
      mutate(
        dsv = dsv_item,
        species = ifelse(is.na(dsv_species), "Unknown", dsv_species),
        N_samples = nrow(current_data),
        N_positive = status_table["1"],
        N_negative = status_table["0"],
        or = exp(estimate)
      ) %>%
      select(dsv, species, N_samples, N_positive, N_negative, estimate, std.error, statistic, p.value, conf.low, conf.high, or)
  }
  close(pb)
  
  if(length(results) == 0) return(NULL)
  bind_rows(results) %>% mutate(ethnicity = target_eth) %>% 
    select(ethnicity, dsv, species, N_samples, N_positive, N_negative, everything())
}

# ----------------------------------------------------------------------------
# 9. Run the analysis for each ethnicity
# ----------------------------------------------------------------------------
if (!dir.exists("results_ethnicity")) dir.create("results_ethnicity")

ethnicities <- as.character(unique(combined_metadata$ethnicity))
all_ethnicity_results <- list()
all_screening_results <- list()

for (eth in ethnicities) {
  cat(sprintf("\n==================================================\n"))
  cat(sprintf("[%s] Analyzing %s (1 vs Others) - Version A (coverage as a covariate)\n", Sys.time(), eth))
  
  screen_res <- screen_dsv_ethnicity_no_cov_filter(combined_data, eth) %>% mutate(target_ethnicity = eth)
  all_screening_results[[eth]] <- screen_res
  
  passed_dsv <- screen_res %>% filter(!is.na(p.value), p.value < 0.05) %>% pull(dsv)
  cat(sprintf("\n  -> Number of dSVs passing screening: %d\n", length(passed_dsv)))
  
  if(length(passed_dsv) == 0) next
  
  eth_data_filtered <- combined_data %>% filter(dsv %in% passed_dsv)
  current_result <- analyze_dsv_ethnicity_no_cov_filter(eth, eth_data_filtered)
  
  if (!is.null(current_result) && nrow(current_result) > 0) {
    current_result <- current_result %>% 
      mutate(p_adj = p.adjust(p.value, method = "BH"),
             fdr_note = ifelse(p_adj < 0.05, "Significant after FDR correction", "Not significant after FDR correction"))
    all_ethnicity_results[[eth]] <- current_result
    write.xlsx(current_result, paste0("results_ethnicity/", eth, "_vs_Others_dsv_results_verA_coverage_as_covariate.xlsx"), rowNames = FALSE)
  }
}

final_results <- bind_rows(all_ethnicity_results)
significant_results <- final_results %>% filter(!is.na(p_adj), p_adj < 0.05) %>% arrange(p_adj)

write.xlsx(final_results, "results_ethnicity/all_ethnicities_dsv_results_verA_coverage_as_covariate.xlsx", rowNames = FALSE)
write.xlsx(significant_results, "results_ethnicity/all_ethnicities_significant_verA_coverage_as_covariate.xlsx", rowNames = FALSE)

# ============================================================================
# 10-fold cross-validation for significant dSVs without a coverage cutoff
# ============================================================================

# Load logistf for Firth logistic regression
library(logistf)

run_cv_with_varimp_dsv_no_cov_filter <- function(target_ethnicity, dsv_name, data_meta, dsv_data_clean) {
  
  # Ensure dsv_data_clean contains sample_id
  if (!"sample_id" %in% colnames(dsv_data_clean)) {
    dsv_data_clean <- dsv_data_clean %>% rownames_to_column("sample_id")
  }
  
  if (!dsv_name %in% colnames(dsv_data_clean)) return(NULL)
  
  current_dsv_df <- dsv_data_clean %>%
    select(sample_id, all_of(dsv_name)) %>%
    rename(status = !!sym(dsv_name))
  
  dsv_species <- str_extract(dsv_name, "^[A-Za-z0-9.]+\\.[a-zA-Z0-9_.]+|^[A-Za-z0-9]+_[a-zA-Z0-9_]+")
  if (is.na(dsv_species)) {
    dsv_species <- str_split(dsv_name, "_")[[1]][1]
  }
  cov_col_name <- paste0(dsv_species, "_coverage")
  
  df <- data_meta %>%
    inner_join(current_dsv_df, by = "sample_id") %>%
    mutate(
      ethnicity_binary = ifelse(ethnicity == target_ethnicity, 1, 0),
      species_abundance = if (!is.na(dsv_species) && dsv_species %in% colnames(.)) .[[dsv_species]] else 0,
      species_coverage = if (!is.na(cov_col_name) && cov_col_name %in% colnames(.)) .[[cov_col_name]] else species_abundance
    ) %>%
    drop_na(status, ethnicity_binary, Age, Sex, BMI, residency,
            diet, medication, 
            urbanization, general_metadata,
            species_abundance, species_coverage)
  
  if (nrow(df) < 20) return(NULL)
  
  status_table <- table(df$status)
  if (length(status_table) < 2 || any(status_table < 5)) return(NULL)
  
  covars_all <- c("ethnicity_binary", "species_coverage", "species_abundance", "residency",
                  "urbanization", "diet", "medication", 
                  "Age", "Sex", "BMI", "general_metadata")
  
  valid_covars <- covars_all[sapply(df[covars_all], function(x) length(unique(x)) > 1)]
  
  if (!"ethnicity_binary" %in% valid_covars) return(NULL)
  
  full_formula <- as.formula(paste("status ~", paste(valid_covars, collapse = " + ")))
  
  # Use Firth logistic regression to handle complete separation
  fit_full <- tryCatch({
    logistf(full_formula, data = df, family = binomial())
  }, error = function(e) NULL)
  if (is.null(fit_full)) return(NULL)
  
  # McFadden R² = 1 - (logLik_full / logLik_null)
  null_formula <- as.formula("status ~ 1")
  fit_null <- logistf(null_formula, data = df, family = binomial())
  r2_full <- 1 - (fit_full$loglik / fit_null$loglik)
  
  get_marginal_delta_r2 <- function(drop_var) {
    if (!drop_var %in% valid_covars) return(0)
    sub_covars <- setdiff(valid_covars, drop_var)
    if (length(sub_covars) == 0) return(r2_full)
    f_drop <- as.formula(paste("status ~", paste(sub_covars, collapse = " + ")))
    fit_drop <- tryCatch({
      logistf(f_drop, data = df, family = binomial())
    }, error = function(e) NULL)
    if (is.null(fit_drop)) return(NA_real_)
    r2_drop <- 1 - (fit_drop$loglik / fit_null$loglik)
    return(max(0, r2_full - r2_drop))
  }
  
  marginal_r2_ethnicity   <- get_marginal_delta_r2("ethnicity_binary")
  marginal_r2_coverage    <- get_marginal_delta_r2("species_coverage")
  marginal_r2_residency   <- get_marginal_delta_r2("residency")
  marginal_r2_urban       <- get_marginal_delta_r2("urbanization")
  marginal_r2_diet        <- get_marginal_delta_r2("diet")
  marginal_r2_spec        <- get_marginal_delta_r2("species_abundance")
  
  # 10-fold cross-validation with valid factor labels for status
  df_cv <- df
  df_cv$status <- factor(ifelse(df_cv$status == 1, "Pos", "Neg"), levels = c("Neg", "Pos"))
  
  ctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)
  
  model_cv <- tryCatch({
    suppressWarnings(train(full_formula, data = df_cv, method = "glm", family = binomial,
                           trControl = ctrl, metric = "ROC"))
  }, error = function(e) NULL)
  
  if (is.null(model_cv)) return(NULL)
  
  resample_stats <- model_cv$resample %>%
    summarise(
      ROC_SD = sd(ROC, na.rm = TRUE),
      Sens_SD = sd(Sens, na.rm = TRUE),
      Spec_SD = sd(Spec, na.rm = TRUE)
    )
  
  imp_df <- varImp(model_cv, scale = TRUE)$importance %>%
    rownames_to_column(var = "predictor") %>%
    rename(importance = Overall)
  
  get_imp_val <- function(p_name) {
    val <- imp_df %>% filter(grepl(p_name, predictor, fixed = TRUE)) %>% pull(importance)
    if (length(val) == 0) return(0) else return(round(max(val), 2))
  }
  
  return(data.frame(
    Target_Ethnicity = target_ethnicity,
    dsv = dsv_name,
    species = ifelse(is.na(dsv_species), "Unknown", dsv_species),
    N_no_filter = nrow(df),
    Full_Model_R2 = round(r2_full, 4),
    CV_Mean_ROC = round(mean(model_cv$results$ROC, na.rm = TRUE), 4),
    CV_Mean_Sens = round(mean(model_cv$results$Sens, na.rm = TRUE), 4),
    CV_Mean_Spec = round(mean(model_cv$results$Spec, na.rm = TRUE), 4),
    CV_ROC_SD = round(resample_stats$ROC_SD, 4),
    CV_Sens_SD = round(resample_stats$Sens_SD, 4),
    CV_Spec_SD = round(resample_stats$Spec_SD, 4),
    Marginal_Delta_R2_Ethnicity = round(marginal_r2_ethnicity, 4),
    Marginal_Delta_R2_SpeciesCoverage = round(marginal_r2_coverage, 4),
    Marginal_Delta_R2_Residency = round(marginal_r2_residency, 4),
    Marginal_Delta_R2_Urbanization = round(marginal_r2_urban, 4),
    Marginal_Delta_R2_Diet = round(marginal_r2_diet, 4),
    Marginal_Delta_R2_SpeciesAbundance = round(marginal_r2_spec, 4),
    Imp_Ethnicity = get_imp_val("ethnicity_binary"),
    Imp_Residency = get_imp_val("residency"),
    Imp_SpeciesCoverage = get_imp_val("species_coverage"),
    Imp_SpeciesAbundance = get_imp_val("species_abundance"),
    Imp_Urbanization = get_imp_val("urbanization"),
    Imp_Diet = get_imp_val("diet"),
    Imp_Medication = get_imp_val("medication"),
    Imp_Age = get_imp_val("Age"),
    Imp_BMI = get_imp_val("BMI"),
    stringsAsFactors = FALSE
  ))
}

# Re-run 10-fold cross-validation for significant ethnicity-dSV pairs
sig_dsv_pairs <- significant_results %>% select(ethnicity, dsv) %>% distinct()
if (nrow(sig_dsv_pairs) > 0) {
  cat(sprintf("\nRunning 10-fold CV for %d significant ethnicity-dSV pairs (no coverage cutoff)...\n", nrow(sig_dsv_pairs)))
  cv_list <- list()
  
  for (i in 1:nrow(sig_dsv_pairs)) {
    eth <- sig_dsv_pairs$ethnicity[i]
    dsv <- sig_dsv_pairs$dsv[i]
    
    res_cv <- run_cv_with_varimp_dsv_no_cov_filter(eth, dsv, combined_metadata, dsv_filtered)
    if (!is.null(res_cv)) cv_list[[i]] <- res_cv
  }
  
  cv_df <- bind_rows(cv_list)
  
  if (!dir.exists("results")) dir.create("results")
  
  openxlsx::write.xlsx(cv_df, "results/cv_10fold_ethnicity_dsv_sig_svs_no_filter.xlsx", rowNames = FALSE)
  cat(sprintf("Results saved successfully: %d records\n", nrow(cv_df)))
}


# ============================================================================
# dSV multi-ethnicity analysis - Version B: sequencing coverage >= 5X
# Sensitivity analysis with a strict sequencing coverage cutoff
# ============================================================================

library(tidyverse)
library(broom)
library(openxlsx)
library(lmtest)
library(caret)

set.seed(12345)

# ----------------------------------------------------------------------------
# 1. Data paths and input files
# ----------------------------------------------------------------------------
data_dir <- "."

cat(sprintf("[%s] Reading input files...\n", Sys.time()))
dsv_data <- read.csv(file.path(data_dir, "final_filtered_dsv.csv"), row.names = 1)
metadata <- read.csv("metadata.csv")
diet_data <- read.csv("final_filtered_diet.csv", check.names = FALSE)
medication_data <- read.csv("Medication.csv", check.names = FALSE)
urbanization_data <- read.csv("Urbanization.csv", check.names = FALSE)
general_metadata_data <- read.csv("General metadata.csv", check.names = FALSE)

# 2. Check required metadata columns
required_cols <- c("sample_id", "residency", "Age", "Sex", "BMI", "ethnicity")
missing_cols <- setdiff(required_cols, names(metadata))
if (length(missing_cols) > 0) stop("Missing required columns in metadata: ", paste(missing_cols, collapse = ", "))

# 3. Construct aggregate covariates
process_covariate <- function(data, covariate_name) {
  vars_df <- data %>% select(-sample_id)
  numeric_df <- vars_df %>%
    mutate(across(everything(), function(col) {
      if (is.numeric(col)) return(col)
      col_char <- as.character(col) %>% trimws()
      col_lower <- tolower(col_char)
      if (all(col_lower %in% c("yes", "no", "y", "n", "true", "false", "1", "0", "", NA))) {
        return(case_when(
          col_lower %in% c("yes", "y", "true", "1") ~ 1,
          col_lower %in% c("no", "n", "false", "0") ~ 0,
          TRUE ~ NA_real_
        ))
      }
      num_col <- suppressWarnings(as.numeric(col_char))
      if (all(is.na(num_col[!is.na(col_char) & col_char != ""]))) {
        num_col <- as.numeric(as.factor(col_char))
      }
      return(num_col)
    }))
  
  result <- data %>% select(sample_id) %>% mutate(!!covariate_name := rowSums(numeric_df, na.rm = TRUE))
  return(result)
}

cat(sprintf("[%s] Constructing aggregate covariates...\n", Sys.time()))
diet_total <- process_covariate(diet_data, "diet")
medication_total <- process_covariate(medication_data, "medication")
urbanization_total <- process_covariate(urbanization_data, "urbanization")
general_metadata_total <- process_covariate(general_metadata_data, "general_metadata")

# 4. Read species abundance and sequencing coverage matrices
cat(sprintf("[%s] Reading species abundance and sequencing coverage data...\n", Sys.time()))
species_map <- read.csv("sv_to_metaphlan_final_mapping.csv")
metaphlan_raw <- read_tsv("yunnan_7_years_ago_combined_metaphlan_profile.tsv", comment = "", show_col_types = FALSE)
if (grepl("#", colnames(metaphlan_raw)[1])) colnames(metaphlan_raw)[1] <- "clade_name"

abundance_matrix <- metaphlan_raw %>%
  filter(clade_name %in% unique(species_map$metaphlan_clade)) %>%
  left_join(species_map, by = c("clade_name" = "metaphlan_clade")) %>%
  select(-clade_name, -metaphlan_species) %>%
  pivot_longer(cols = -sv_species, names_to = "sample_id", values_to = "abundance") %>%
  mutate(abundance = as.numeric(abundance)) %>%
  pivot_wider(id_cols = sample_id, names_from = sv_species, values_from = abundance, values_fill = list(abundance = 0))

raw_cov_data <- read.csv("merged_species_mean_coverage_matrix.csv", check.names = FALSE)
colnames(raw_cov_data)[1] <- "sample_id"
coverage_raw <- raw_cov_data %>%
  pivot_longer(cols = -sample_id, names_to = "sv_species", values_to = "coverage") %>%
  mutate(coverage = as.numeric(coverage)) %>%
  group_by(sample_id, sv_species) %>% summarise(coverage = mean(coverage, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(id_cols = sample_id, names_from = sv_species, values_from = coverage, values_fill = list(coverage = 0))
colnames(coverage_raw)[-1] <- paste0(colnames(coverage_raw)[-1], "_coverage")

# 5. Merge metadata and covariates
combined_metadata <- metadata %>%
  inner_join(diet_total, by = "sample_id") %>% inner_join(medication_total, by = "sample_id") %>%
  inner_join(urbanization_total, by = "sample_id") %>% inner_join(general_metadata_total, by = "sample_id") %>%
  left_join(abundance_matrix, by = "sample_id") %>% left_join(coverage_raw, by = "sample_id") %>%
  mutate(Sex = as.factor(Sex), residency = as.factor(residency), ethnicity = as.factor(ethnicity))

cat("\n================ Cohort summary ================\n")
cat("Sample counts by ethnicity:\n")
print(table(combined_metadata$ethnicity))
cat("==============================================\n\n")

# 6. Preprocess dSV data (retain NA values)
cat(sprintf("[%s] Preprocessing dSV data (retaining NA values)...\n", Sys.time()))
dsv_clean <- dsv_data %>%
  mutate(across(everything(), as.numeric))

# Filter low-prevalence dSVs (>=10%; NA values excluded from prevalence calculation)
dsv_prevalence <- colMeans(dsv_clean != 0, na.rm = TRUE)
dsv_filtered <- dsv_clean[, names(dsv_prevalence[dsv_prevalence >= 0.1]), drop = FALSE]

# Convert to long format while retaining NA values
combined_data <- dsv_filtered %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(cols = -sample_id, names_to = "dsv", values_to = "status") %>%
  mutate(status = as.numeric(status)) %>%  # retain NA values
  inner_join(combined_metadata, by = "sample_id")

# ----------------------------------------------------------------------------
# 7. Screening function: coverage >= 5X, chi-squared test (p < 0.05)
# ----------------------------------------------------------------------------
screen_dsv_ethnicity_cov5x <- function(data, target_eth) {
  data <- data %>% mutate(ethnicity_binary = ifelse(ethnicity == target_eth, 1, 0))
  results <- list()
  dsv_names <- unique(data$dsv)
  total_dsv <- length(dsv_names)
  
  cat(sprintf("  -> Chi-squared screening (coverage >= 5X): %d features\n", total_dsv))
  pb <- txtProgressBar(min = 0, max = total_dsv, style = 3)
  
  for (i in seq_along(dsv_names)) {
    setTxtProgressBar(pb, i)
    dsv_item <- dsv_names[i]
    dsv_species <- str_extract(dsv_item, "^[A-Za-z0-9.]+\\.[a-zA-Z0-9_.]+|^[A-Za-z0-9]+_[a-zA-Z0-9_]+")
    if (is.na(dsv_species)) dsv_species <- str_split(dsv_item, "_")[[1]][1]
    
    cov_col_name <- paste0(dsv_species, "_coverage")
    cov_vec <- if (!is.na(cov_col_name) && cov_col_name %in% colnames(data)) data[[cov_col_name]] else 0
    
    current_data <- data %>%
      filter(dsv == !!dsv_item) %>%
      mutate(species_coverage = cov_vec[data$dsv == !!dsv_item]) %>%
      drop_na(status, ethnicity_binary, species_coverage) %>%
      filter(species_coverage >= 5)  # apply coverage cutoff
    
    if (nrow(current_data) < 30) {
      results[[dsv_item]] <- tibble(dsv = dsv_item, p.value = NA, N_cov5x = nrow(current_data), note = "Insufficient sample size (<30)")
      next
    }
    
    tbl <- table(current_data$status, current_data$ethnicity_binary)
    if (any(rowSums(tbl) == 0) || any(colSums(tbl) == 0)) {
      results[[dsv_item]] <- tibble(dsv = dsv_item, p.value = NA, N_cov5x = nrow(current_data), note = "Invalid contingency table")
      next
    }
    
    chi_test <- tryCatch({ chisq.test(tbl) }, error = function(e) NULL)
    if (is.null(chi_test)) {
      fisher_test <- fisher.test(tbl)
      p_val <- fisher_test$p.value
      note <- "Fisher exact test"
    } else {
      if (any(chi_test$expected < 5)) {
        fisher_test <- fisher.test(tbl)
        p_val <- fisher_test$p.value
        note <- "Fisher exact test (expected count <5)"
      } else {
        p_val <- chi_test$p.value
        note <- "Chi-squared test"
      }
    }
    
    results[[dsv_item]] <- tibble(
      dsv = dsv_item,
      p.value = p_val,
      N_cov5x = nrow(current_data),
      note = note
    )
  }
  close(pb)
  bind_rows(results)
}

# ----------------------------------------------------------------------------
# 8. Logistic regression restricted to coverage >= 5X
# ----------------------------------------------------------------------------
analyze_dsv_ethnicity_cov5x <- function(target_eth, data) {
  data <- data %>% mutate(ethnicity_binary = ifelse(ethnicity == target_eth, 1, 0))
  results <- list()
  dsv_names <- unique(data$dsv)
  total_dsv <- length(dsv_names)
  
  cat(sprintf("\n  -> Multivariable logistic regression (coverage >= 5X): %d candidate features\n", total_dsv))
  pb <- txtProgressBar(min = 0, max = total_dsv, style = 3)
  
  for (i in seq_along(dsv_names)) {
    setTxtProgressBar(pb, i)
    dsv_item <- dsv_names[i]
    dsv_species <- str_extract(dsv_item, "^[A-Za-z0-9.]+\\.[a-zA-Z0-9_.]+|^[A-Za-z0-9]+_[a-zA-Z0-9_]+")
    if (is.na(dsv_species)) dsv_species <- str_split(dsv_item, "_")[[1]][1]
    
    cov_col_name <- paste0(dsv_species, "_coverage")
    abun_vec <- if (!is.na(dsv_species) && dsv_species %in% colnames(data)) data[[dsv_species]] else 0
    cov_vec  <- if (!is.na(cov_col_name) && cov_col_name %in% colnames(data)) data[[cov_col_name]] else abun_vec
    
    current_data <- data %>%
      filter(dsv == !!dsv_item) %>%
      mutate(species_abundance = abun_vec[data$dsv == !!dsv_item], 
             species_coverage = cov_vec[data$dsv == !!dsv_item]) %>%
      drop_na(status, ethnicity_binary, Age, Sex, BMI, residency,
              diet, medication, urbanization, general_metadata,
              species_abundance, species_coverage) %>%
      filter(species_coverage >= 5)  # apply coverage cutoff
    
    if (nrow(current_data) < 10) next
    
    status_table <- table(current_data$status)
    if (length(status_table) < 2 || any(status_table < 3)) next
    
    model_formula <- as.formula(
      "status ~ ethnicity_binary + species_coverage + species_abundance + Age + Sex + BMI + residency +
       diet + medication + urbanization + general_metadata"
    )
    
    model <- tryCatch({ glm(model_formula, data = current_data, family = binomial()) }, error = function(e) NULL)
    if (is.null(model)) next
    
    tidy_res <- tryCatch({
      tidy(model, conf.int = TRUE)
    }, error = function(e) {
      tidy(model, conf.int = FALSE) %>%
        mutate(
          conf.low = estimate - 1.96 * std.error,
          conf.high = estimate + 1.96 * std.error
        )
    })
    if (is.null(tidy_res) || nrow(tidy_res) == 0) next
    
    results[[dsv_item]] <- tidy_res %>%
      filter(term == "ethnicity_binary") %>%
      mutate(
        dsv = dsv_item,
        species = ifelse(is.na(dsv_species), "Unknown", dsv_species),
        N_samples_cov5x = nrow(current_data),
        N_positive = status_table["1"],
        N_negative = status_table["0"],
        or = exp(estimate)
      ) %>%
      select(dsv, species, N_samples_cov5x, N_positive, N_negative, estimate, std.error, statistic, p.value, conf.low, conf.high, or)
  }
  close(pb)
  
  if(length(results) == 0) return(NULL)
  bind_rows(results) %>% mutate(ethnicity = target_eth) %>% 
    select(ethnicity, dsv, species, N_samples_cov5x, N_positive, N_negative, everything())
}

# ----------------------------------------------------------------------------
# 9. Run the analysis for each ethnicity
# ----------------------------------------------------------------------------
if (!dir.exists("results_ethnicity")) dir.create("results_ethnicity")

ethnicities <- as.character(unique(combined_metadata$ethnicity))
all_ethnicity_results <- list()
all_screening_results <- list()

for (eth in ethnicities) {
  cat(sprintf("\n==================================================\n"))
  cat(sprintf("[%s] Analyzing %s (1 vs Others) - Version B (coverage >= 5X)\n", Sys.time(), eth))
  
  screen_res <- screen_dsv_ethnicity_cov5x(combined_data, eth) %>% mutate(target_ethnicity = eth)
  all_screening_results[[eth]] <- screen_res
  
  passed_dsv <- screen_res %>% filter(!is.na(p.value), p.value < 0.05) %>% pull(dsv)
  cat(sprintf("\n  -> Number of dSVs passing screening: %d\n", length(passed_dsv)))
  
  if(length(passed_dsv) == 0) next
  
  eth_data_filtered <- combined_data %>% filter(dsv %in% passed_dsv)
  current_result <- analyze_dsv_ethnicity_cov5x(eth, eth_data_filtered)
  
  if (!is.null(current_result) && nrow(current_result) > 0) {
    current_result <- current_result %>% 
      mutate(p_adj = p.adjust(p.value, method = "BH"),
             fdr_note = ifelse(p_adj < 0.05, "Significant after FDR correction", "Not significant after FDR correction"))
    all_ethnicity_results[[eth]] <- current_result
    write.xlsx(current_result, paste0("results_ethnicity/", eth, "_vs_Others_dsv_results_verB_cov5x.xlsx"), rowNames = FALSE)
  }
}

final_results <- bind_rows(all_ethnicity_results)
significant_results <- final_results %>% filter(!is.na(p_adj), p_adj < 0.05) %>% arrange(p_adj)

write.xlsx(final_results, "results_ethnicity/all_ethnicities_dsv_results_verB_cov5x.xlsx", rowNames = FALSE)
write.xlsx(significant_results, "results_ethnicity/all_ethnicities_significant_verB_cov5x.xlsx", rowNames = FALSE)


# 10. 10-fold cross-validation (coverage >= 5X) with marginal R2 estimates
# Report marginal delta R2 and variable importance for all covariates

run_cv_verB_dsv <- function(target_eth, dsv_name, data_meta, dsv_data_clean) {
  
  # Ensure dsv_data_clean contains sample_id
  if (!"sample_id" %in% colnames(dsv_data_clean)) {
    dsv_data_clean <- dsv_data_clean %>% rownames_to_column("sample_id")
  }
  
  if (!dsv_name %in% colnames(dsv_data_clean)) return(NULL)
  
  current_dsv_df <- dsv_data_clean %>% 
    select(sample_id, all_of(dsv_name)) %>%
    rename(status = !!sym(dsv_name))
  
  dsv_species <- str_extract(dsv_name, "^[A-Za-z0-9.]+\\.[a-zA-Z0-9_.]+|^[A-Za-z0-9]+_[a-zA-Z0-9_]+")
  if (is.na(dsv_species)) dsv_species <- str_split(dsv_name, "_")[[1]][1]
  
  cov_col_name <- paste0(dsv_species, "_coverage")
  
  df <- data_meta %>%
    inner_join(current_dsv_df, by = "sample_id") %>%
    mutate(
      ethnicity_binary = ifelse(ethnicity == target_eth, 1, 0),
      species_abundance = if (!is.na(dsv_species) && dsv_species %in% colnames(.)) .[[dsv_species]] else 0,
      species_coverage = if (!is.na(cov_col_name) && cov_col_name %in% colnames(.)) .[[cov_col_name]] else species_abundance
    ) %>%
    drop_na(status, ethnicity_binary, Age, Sex, BMI, residency,
            diet, medication, urbanization, general_metadata,
            species_abundance, species_coverage) %>%
    filter(species_coverage >= 5)  # Version B coverage cutoff
  
  if (nrow(df) < 20) return(NULL)
  
  status_table <- table(df$status)
  if (length(status_table) < 2 || any(status_table < 5)) return(NULL)
  
  covars_all <- c("ethnicity_binary", "species_coverage", "species_abundance", "residency",
                  "urbanization", "diet", "medication", 
                  "Age", "Sex", "BMI", "general_metadata")
  valid_covars <- covars_all[sapply(df[covars_all], function(x) length(unique(x)) > 1)]
  if (!"ethnicity_binary" %in% valid_covars) return(NULL)
  
  full_formula <- as.formula(paste("status ~", paste(valid_covars, collapse = " + ")))
  
  # Full-model McFadden R2 using standard logistic regression
  fit_full <- tryCatch({ glm(full_formula, data = df, family = binomial()) }, error = function(e) NULL)
  if (is.null(fit_full)) return(NULL)
  r2_full <- 1 - (fit_full$deviance / fit_full$null.deviance)
  
  get_marginal_delta_r2 <- function(drop_var) {
    if (!drop_var %in% valid_covars) return(0)
    sub_covars <- setdiff(valid_covars, drop_var)
    if (length(sub_covars) == 0) return(r2_full)
    f_drop <- as.formula(paste("status ~", paste(sub_covars, collapse = " + ")))
    fit_drop <- tryCatch({ glm(f_drop, data = df, family = binomial()) }, error = function(e) NULL)
    if (is.null(fit_drop)) return(NA_real_)
    r2_drop <- 1 - (fit_drop$deviance / fit_drop$null.deviance)
    return(max(0, r2_full - r2_drop))
  }
  
  # Calculate marginal delta R2 for all covariates
  marginal_r2_ethnicity   <- get_marginal_delta_r2("ethnicity_binary")
  marginal_r2_coverage    <- get_marginal_delta_r2("species_coverage")
  marginal_r2_residency   <- get_marginal_delta_r2("residency")
  marginal_r2_urban       <- get_marginal_delta_r2("urbanization")
  marginal_r2_diet        <- get_marginal_delta_r2("diet")
  marginal_r2_spec        <- get_marginal_delta_r2("species_abundance")
  marginal_r2_medication  <- get_marginal_delta_r2("medication")
  marginal_r2_age         <- get_marginal_delta_r2("Age")
  marginal_r2_sex         <- get_marginal_delta_r2("Sex")
  marginal_r2_bmi         <- get_marginal_delta_r2("BMI")
  marginal_r2_general     <- get_marginal_delta_r2("general_metadata")
  
  # 10-fold cross-validation with valid factor labels for status
  df_cv <- df
  df_cv$status <- factor(ifelse(df_cv$status == 1, "Pos", "Neg"), levels = c("Neg", "Pos"))
  
  ctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)
  model_cv <- tryCatch({ 
    suppressWarnings(train(full_formula, data = df_cv, method = "glm", family = binomial,
                           trControl = ctrl, metric = "ROC")) 
  }, error = function(e) NULL)
  if (is.null(model_cv)) return(NULL)
  
  # Cross-validation performance statistics
  resample_stats <- model_cv$resample %>%
    summarise(
      ROC_SD = sd(ROC, na.rm = TRUE),
      Sens_SD = sd(Sens, na.rm = TRUE),
      Spec_SD = sd(Spec, na.rm = TRUE)
    )
  
  # Extract variable importance
  imp_df <- varImp(model_cv, scale = TRUE)$importance %>% 
    rownames_to_column(var = "predictor") %>% rename(importance = Overall)
  
  get_imp_val <- function(p_name) {
    val <- imp_df %>% filter(grepl(p_name, predictor, fixed = TRUE)) %>% pull(importance)
    if (length(val) == 0) return(0) else return(round(max(val), 2))
  }
  
  return(data.frame(
    Target_Ethnicity = target_eth,
    dsv = dsv_name,
    species = ifelse(is.na(dsv_species), "Unknown", dsv_species),
    N_cov5x = nrow(df),
    
    Full_Model_R2 = round(r2_full, 4),
    CV_Mean_ROC = round(mean(model_cv$results$ROC, na.rm = TRUE), 4),
    CV_Mean_Sens = round(mean(model_cv$results$Sens, na.rm = TRUE), 4),
    CV_Mean_Spec = round(mean(model_cv$results$Spec, na.rm = TRUE), 4),
    CV_ROC_SD = round(resample_stats$ROC_SD, 4),
    CV_Sens_SD = round(resample_stats$Sens_SD, 4),
    CV_Spec_SD = round(resample_stats$Spec_SD, 4),
    
    # Marginal delta R2 for all covariates
    Marginal_Delta_R2_Ethnicity = round(marginal_r2_ethnicity, 4),
    Marginal_Delta_R2_SpeciesCoverage = round(marginal_r2_coverage, 4),
    Marginal_Delta_R2_Residency = round(marginal_r2_residency, 4),
    Marginal_Delta_R2_Urbanization = round(marginal_r2_urban, 4),
    Marginal_Delta_R2_Diet = round(marginal_r2_diet, 4),
    Marginal_Delta_R2_SpeciesAbundance = round(marginal_r2_spec, 4),
    Marginal_Delta_R2_Medication = round(marginal_r2_medication, 4),
    Marginal_Delta_R2_Age = round(marginal_r2_age, 4),
    Marginal_Delta_R2_Sex = round(marginal_r2_sex, 4),
    Marginal_Delta_R2_BMI = round(marginal_r2_bmi, 4),
    Marginal_Delta_R2_GeneralMetadata = round(marginal_r2_general, 4),
    
    # Variable importance (0-100)
    Imp_Ethnicity = get_imp_val("ethnicity_binary"),
    Imp_SpeciesCoverage = get_imp_val("species_coverage"),
    Imp_Residency = get_imp_val("residency"),
    Imp_Urbanization = get_imp_val("urbanization"),
    Imp_Diet = get_imp_val("diet"),
    Imp_SpeciesAbundance = get_imp_val("species_abundance"),
    Imp_Medication = get_imp_val("medication"),
    Imp_Age = get_imp_val("Age"),
    Imp_Sex = get_imp_val("Sex"),
    Imp_BMI = get_imp_val("BMI"),
    Imp_GeneralMetadata = get_imp_val("general_metadata"),
    stringsAsFactors = FALSE
  ))
}

if (nrow(significant_results) > 0) {
  cat(sprintf("\n[%s] Running 10-fold CV for significant dSVs (Version B: coverage >= 5X)...\n", Sys.time()))
  cv_list <- list()
  for (i in 1:nrow(significant_results)) {
    eth <- significant_results$ethnicity[i]
    dsv <- significant_results$dsv[i]
    res_cv <- run_cv_verB_dsv(eth, dsv, combined_metadata, dsv_filtered)
    if (!is.null(res_cv)) cv_list[[i]] <- res_cv
  }
  if (length(cv_list) > 0) {
    cv_df <- bind_rows(cv_list)
    write.xlsx(cv_df, "results_ethnicity/cv_10fold_verB_cov5x_dsv.xlsx", rowNames = FALSE)
    cat(sprintf("Version B CV results saved: %d records\n", nrow(cv_df)))
  } else {
    cat("Version B produced no valid CV results\n")
  }
}

cat(sprintf("\n[%s] Version B (coverage >= 5X) analysis completed.\n", Sys.time()))
