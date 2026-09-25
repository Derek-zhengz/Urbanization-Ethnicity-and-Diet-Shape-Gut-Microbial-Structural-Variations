# ============================================================================
# Reproducible dSV ethnicity association models
# Version A: sequencing coverage included as a covariate (no coverage cutoff)
# Version B: sensitivity analysis restricted to sequencing coverage >= 5X
# ============================================================================

library(tidyverse)
library(broom)
library(openxlsx)
library(lmtest)
library(caret)
library(randomForest)
library(logistf)

set.seed(12345)

# ----------------------------------------------------------------------------
# 1. Data paths and input files
# ----------------------------------------------------------------------------
data_dir <- "."

cat(sprintf("[%s] Reading input files...\n", Sys.time()))
dsv_data <- read.csv(file.path(data_dir, "final_filtered_dsv.csv"),
                     row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read.csv("metadata.csv", check.names = FALSE)
diet_data <- read.csv("final_filtered_diet.csv", check.names = FALSE)
medication_data <- read.csv("Medication.csv", check.names = FALSE)
urbanization_data <- read.csv("Urbanization.csv", check.names = FALSE)
general_metadata_data <- read.csv("General metadata.csv", check.names = FALSE)

# Ensure unique dSV feature names
if (anyDuplicated(colnames(dsv_data)) > 0) {
  cat("Duplicate dSV feature names detected; making names unique.\n")
  colnames(dsv_data) <- make.unique(colnames(dsv_data), sep = "__dup")
}

# 2. Check required metadata columns
required_cols <- c("sample_id", "residency", "Age", "Sex", "BMI", "ethnicity")
missing_cols <- setdiff(required_cols, names(metadata))
if (length(missing_cols) > 0) {
  stop("Missing required columns in metadata: ", paste(missing_cols, collapse = ", "))
}

# ----------------------------------------------------------------------------
# 3. Process aggregate covariates (diet, medication, urbanization)
# ----------------------------------------------------------------------------
process_covariate_sum <- function(data, covariate_name) {
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
        # 对多分类字符列，这里仍按因子编码为数值，用于求行均值
        # 若某 aggregate 打分希望保留顺序，请自行指定
        num_col <- as.numeric(as.factor(col_char))
      }
      return(num_col)
    }))
  
  result <- data %>% select(sample_id) %>%
    mutate(!!covariate_name := rowMeans(numeric_df, na.rm = TRUE))
  return(result)
}

# ----------------------------------------------------------------------------
# 4. Expand General Metadata: each column independently, with prefix
#    - binary yes/no -> 0/1
#    - pure numeric   -> numeric
#    - multi-class    -> factor (NOT numeric)
#    - add prefix "general__" to avoid name collision
# ----------------------------------------------------------------------------
process_general_metadata_columns <- function(data, prefix = "general__") {
  sample_ids <- data %>% select(sample_id)
  vars_df <- data %>% select(-sample_id)
  
  original_names <- colnames(vars_df)
  safe_names <- paste0(prefix, make.names(original_names, unique = TRUE))
  colnames(vars_df) <- safe_names
  
  cleaned_vars <- vars_df %>%
    mutate(across(everything(), function(col) {
      if (is.numeric(col)) return(col)
      
      col_char <- as.character(col) %>% trimws()
      col_lower <- tolower(col_char)
      
      # yes/no/true/false -> 0/1
      if (all(col_lower %in% c("yes", "no", "y", "n", "true", "false", "1", "0", "", NA))) {
        return(case_when(
          col_lower %in% c("yes", "y", "true", "1") ~ 1,
          col_lower %in% c("no", "n", "false", "0") ~ 0,
          TRUE ~ NA_real_
        ))
      }
      
      # pure numeric string -> numeric
      num_col <- suppressWarnings(as.numeric(col_char))
      non_missing_original <- !is.na(col_char) & col_char != ""
      if (sum(!is.na(num_col)) == sum(non_missing_original)) {
        return(num_col)
      }
      
      # otherwise -> factor (do NOT coerce to integer)
      factor(col_char)
    }))
  
  bind_cols(sample_ids, cleaned_vars)
}

cat(sprintf("[%s] Processing covariates and expanding general metadata columns...\n", Sys.time()))
diet_total <- process_covariate_sum(diet_data, "diet")
medication_total <- process_covariate_sum(medication_data, "medication")
urbanization_total <- process_covariate_sum(urbanization_data, "urbanization")

general_metadata_expanded <- process_general_metadata_columns(general_metadata_data, "general__")
gen_meta_cols <- setdiff(colnames(general_metadata_expanded), "sample_id")

cat(sprintf("General metadata columns incorporated individually: %d\n", length(gen_meta_cols)))

# ----------------------------------------------------------------------------
# 5. Read species abundance and sequencing coverage matrices
# ----------------------------------------------------------------------------
cat(sprintf("[%s] Reading species abundance and sequencing coverage data...\n", Sys.time()))
species_map <- read.csv("sv_to_metaphlan_final_mapping.csv", check.names = FALSE)
metaphlan_raw <- read_tsv("yunnan_7_years_ago_combined_metaphlan_profile.tsv",
                          comment = "", show_col_types = FALSE)
if (grepl("#", colnames(metaphlan_raw)[1])) colnames(metaphlan_raw)[1] <- "clade_name"

abundance_matrix <- metaphlan_raw %>%
  filter(clade_name %in% unique(species_map$metaphlan_clade)) %>%
  left_join(species_map, by = c("clade_name" = "metaphlan_clade")) %>%
  select(-clade_name, -metaphlan_species) %>%
  pivot_longer(cols = -sv_species, names_to = "sample_id", values_to = "abundance") %>%
  mutate(abundance = suppressWarnings(as.numeric(abundance))) %>%
  group_by(sample_id, sv_species) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(id_cols = sample_id, names_from = sv_species,
              values_from = abundance, values_fill = 0)

raw_cov_data <- read.csv("merged_species_mean_coverage_matrix.csv", check.names = FALSE)
colnames(raw_cov_data)[1] <- "sample_id"
coverage_raw <- raw_cov_data %>%
  pivot_longer(cols = -sample_id, names_to = "sv_species", values_to = "coverage") %>%
  mutate(coverage = suppressWarnings(as.numeric(coverage))) %>%
  group_by(sample_id, sv_species) %>%
  summarise(coverage = mean(coverage, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(id_cols = sample_id, names_from = sv_species,
              values_from = coverage, values_fill = 0)
colnames(coverage_raw)[-1] <- paste0(colnames(coverage_raw)[-1], "_coverage")

# ----------------------------------------------------------------------------
# 6. Merge metadata and all covariates
# ----------------------------------------------------------------------------

combined_metadata <- metadata %>%
  inner_join(diet_total, by = "sample_id") %>%
  inner_join(medication_total, by = "sample_id") %>%
  inner_join(urbanization_total, by = "sample_id") %>%
  inner_join(general_metadata_expanded, by = "sample_id") %>%
  left_join(abundance_matrix, by = "sample_id") %>%
  left_join(coverage_raw, by = "sample_id") %>%
  mutate(
    Sex = as.factor(Sex),
    residency = as.factor(residency),
    ethnicity = as.factor(ethnicity)
  )

cat("\n================ Cohort summary ================\n")
cat("Sample counts by ethnicity:\n")
print(table(combined_metadata$ethnicity))
cat(sprintf("Number of General Metadata variables incorporated individually: %d\n",
            length(gen_meta_cols)))

# Check complete-case sample size for general metadata
complete_gen <- sum(complete.cases(combined_metadata[, gen_meta_cols, drop = FALSE]))
cat(sprintf("Complete cases across general metadata columns: %d / %d\n",
            complete_gen, nrow(combined_metadata)))
cat("==============================================\n\n")

# ----------------------------------------------------------------------------
# 7. Preprocess dSV data (no prevalence filtering, retaining NA values)
# ----------------------------------------------------------------------------
cat(sprintf("[%s] Preprocessing dSV data...\n", Sys.time()))
dsv_clean <- dsv_data %>% mutate(across(everything(), as.numeric))
dsv_filtered <- dsv_clean

combined_data <- dsv_filtered %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(cols = -sample_id, names_to = "dsv", values_to = "status") %>%
  mutate(status = as.numeric(status)) %>%
  inner_join(combined_metadata, by = "sample_id")

# Dynamic general metadata formula part (backtick-quoted)
if (length(gen_meta_cols) > 0) {
  gen_meta_formula_part <- paste(sprintf("`%s`", gen_meta_cols), collapse = " + ")
} else {
  gen_meta_formula_part <- NULL
}

# Helper to build the formula string
build_model_formula <- function(response) {
  base_terms <- "ethnicity_binary + species_coverage + species_abundance + Age + Sex + BMI + residency + diet + medication + urbanization"
  if (!is.null(gen_meta_formula_part) && nzchar(gen_meta_formula_part)) {
    as.formula(paste(response, "~", base_terms, "+", gen_meta_formula_part))
  } else {
    as.formula(paste(response, "~", base_terms))
  }
}

# ----------------------------------------------------------------------------
# 8. Screening function: chi-squared test (no coverage cutoff)
# ----------------------------------------------------------------------------
screen_dsv_ethnicity_no_cov_filter <- function(data, target_eth) {
  data <- data %>% mutate(ethnicity_binary = ifelse(ethnicity == target_eth, 1, 0))
  results <- list()
  dsv_names <- unique(data$dsv)
  total_dsv <- length(dsv_names)
  
  cat(sprintf("  -> Chi-squared screening: %d features\n", total_dsv))
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
      drop_na(status, ethnicity_binary, species_coverage)
    
    if (nrow(current_data) < 30) {
      results[[dsv_item]] <- tibble(dsv = dsv_item, p.value = NA,
                                    N_cov_available = nrow(current_data),
                                    note = "Insufficient sample size (<30)")
      next
    }
    
    tbl <- table(current_data$status, current_data$ethnicity_binary)
    if (any(rowSums(tbl) == 0) || any(colSums(tbl) == 0)) {
      results[[dsv_item]] <- tibble(dsv = dsv_item, p.value = NA,
                                    N_cov_available = nrow(current_data),
                                    note = "Invalid contingency table")
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
# 9. Logistic regression: all general metadata items as separate covariates
# ----------------------------------------------------------------------------
analyze_dsv_ethnicity_no_cov_filter <- function(target_eth, data,
                                                save_full = FALSE,
                                                full_out_dir = "results_ethnicity/full_models") {
  data <- data %>% mutate(ethnicity_binary = ifelse(ethnicity == target_eth, 1, 0))
  results <- list()
  full_results <- list()
  dsv_names <- unique(data$dsv)
  total_dsv <- length(dsv_names)
  
  cat(sprintf("\n  -> Multivariable logistic regression (individual general metadata items): %d candidate features\n",
              total_dsv))
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
             species_coverage  = cov_vec[data$dsv == !!dsv_item]) %>%
      drop_na(
        status, ethnicity_binary, Age, Sex, BMI, residency,
        diet, medication, urbanization,
        any_of(gen_meta_cols),
        species_abundance, species_coverage
      ) %>%
      droplevels()
    
    if (nrow(current_data) < 10) next
    
    status_table <- table(current_data$status)
    if (length(status_table) < 2 || any(status_table < 3)) next
    
    model_formula <- build_model_formula("status")
    
    model <- tryCatch({ glm(model_formula, data = current_data, family = binomial()) },
                      error = function(e) NULL)
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
        N_samples = nrow(current_data),
        N_positive = status_table["1"],
        N_negative = status_table["0"],
        or = exp(estimate)
      ) %>%
      select(dsv, species, N_samples, N_positive, N_negative,
             estimate, std.error, statistic, p.value,
             conf.low, conf.high, or)
    
    if (save_full) {
      full_results[[dsv_item]] <- tidy_res %>% mutate(dsv = dsv_item)
    }
  }
  close(pb)
  
  if (save_full && length(full_results) > 0) {
    if (!dir.exists(full_out_dir)) dir.create(full_out_dir, recursive = TRUE)
    openxlsx::write.xlsx(bind_rows(full_results),
                         file.path(full_out_dir,
                                   paste0("full_model_", target_eth, ".xlsx")),
                         rowNames = FALSE)
  }
  
  if (length(results) == 0) return(NULL)
  bind_rows(results) %>% mutate(ethnicity = target_eth) %>%
    select(ethnicity, dsv, species, N_samples, N_positive, N_negative, everything())
}

# ----------------------------------------------------------------------------
# 10. Run Version A Analysis
# ----------------------------------------------------------------------------
if (!dir.exists("results_ethnicity")) dir.create("results_ethnicity")

ethnicities <- as.character(unique(combined_metadata$ethnicity))
all_ethnicity_results <- list()
all_screening_results <- list()

for (eth in ethnicities) {
  cat(sprintf("\n==================================================\n"))
  cat(sprintf("[%s] Analyzing %s (1 vs Others) - Version A\n", Sys.time(), eth))
  
  screen_res <- screen_dsv_ethnicity_no_cov_filter(combined_data, eth) %>%
    mutate(target_ethnicity = eth)
  all_screening_results[[eth]] <- screen_res
  
  passed_dsv <- screen_res %>% filter(!is.na(p.value), p.value < 0.05) %>% pull(dsv)
  cat(sprintf("\n  -> Number of dSVs passing screening: %d\n", length(passed_dsv)))
  
  if (length(passed_dsv) == 0) next
  
  eth_data_filtered <- combined_data %>% filter(dsv %in% passed_dsv)
  current_result <- analyze_dsv_ethnicity_no_cov_filter(
    eth, eth_data_filtered,
    save_full = TRUE,
    full_out_dir = "results_ethnicity/full_models_verA"
  )
  
  if (!is.null(current_result) && nrow(current_result) > 0) {
    current_result <- current_result %>%
      mutate(p_adj = p.adjust(p.value, method = "BH"),
             fdr_note = ifelse(p_adj < 0.05,
                               "Significant after FDR correction",
                               "Not significant after FDR correction"))
    all_ethnicity_results[[eth]] <- current_result
    write.xlsx(current_result,
               paste0("results_ethnicity/", eth,
                      "_vs_Others_dsv_results_verA_coverage_as_covariate.xlsx"),
               rowNames = FALSE)
  }
}

final_results <- bind_rows(all_ethnicity_results)
if (nrow(final_results) > 0) {
  significant_results <- final_results %>%
    filter(!is.na(p_adj), p_adj < 0.05) %>% arrange(p_adj)
  write.xlsx(final_results,
             "results_ethnicity/all_ethnicities_dsv_results_verA_coverage_as_covariate.xlsx",
             rowNames = FALSE)
  write.xlsx(significant_results,
             "results_ethnicity/all_ethnicities_significant_verA_coverage_as_covariate.xlsx",
             rowNames = FALSE)
}

# ----------------------------------------------------------------------------
# 11. Version B: Sensitivity analysis restricted to sequencing coverage >= 5X
# ----------------------------------------------------------------------------
analyze_dsv_ethnicity_cov5x <- function(target_eth, data,
                                        save_full = FALSE,
                                        full_out_dir = "results_ethnicity/full_models_verB") {
  data <- data %>% mutate(ethnicity_binary = ifelse(ethnicity == target_eth, 1, 0))
  results <- list()
  full_results <- list()
  dsv_names <- unique(data$dsv)
  total_dsv <- length(dsv_names)
  
  cat(sprintf("\n  -> Multivariable logistic regression (coverage >= 5X): %d candidate features\n",
              total_dsv))
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
             species_coverage  = cov_vec[data$dsv == !!dsv_item]) %>%
      drop_na(
        status, ethnicity_binary, Age, Sex, BMI, residency,
        diet, medication, urbanization,
        any_of(gen_meta_cols),
        species_abundance, species_coverage
      ) %>%
      filter(species_coverage >= 5) %>%
      droplevels()
    
    if (nrow(current_data) < 10) next
    
    status_table <- table(current_data$status)
    if (length(status_table) < 2 || any(status_table < 3)) next
    
    model_formula <- build_model_formula("status")
    
    model <- tryCatch({ glm(model_formula, data = current_data, family = binomial()) },
                      error = function(e) NULL)
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
      select(dsv, species, N_samples_cov5x, N_positive, N_negative,
             estimate, std.error, statistic, p.value,
             conf.low, conf.high, or)
    
    if (save_full) {
      full_results[[dsv_item]] <- tidy_res %>% mutate(dsv = dsv_item)
    }
  }
  close(pb)
  
  if (save_full && length(full_results) > 0) {
    if (!dir.exists(full_out_dir)) dir.create(full_out_dir, recursive = TRUE)
    openxlsx::write.xlsx(bind_rows(full_results),
                         file.path(full_out_dir,
                                   paste0("full_model_", target_eth, "_cov5x.xlsx")),
                         rowNames = FALSE)
  }
  
  if (length(results) == 0) return(NULL)
  bind_rows(results) %>% mutate(ethnicity = target_eth) %>%
    select(ethnicity, dsv, species, N_samples_cov5x, N_positive, N_negative, everything())
}


# ----------------------------------------------------------------------------
# 12. 10-fold cross-validation for significant dSVs
#     + McFadden R2
#     + Marginal Delta R2
#     + Random Forest variable importance
#
# Reference logic:
#   - Only dSVs significant after Version A BH-FDR correction are validated.
#   - Outcome: dSV presence/absence (Class_0 / Class_1).
#   - Predictors: ethnicity + the same covariates used in the main model.
#   - 10-fold CV uses Random Forest and ROC as the optimization metric.
# ----------------------------------------------------------------------------

cat("\n========== 12. 10-fold CV + McFadden R2 + Marginal Delta R2 ==========\n")

if (!exists("significant_results") || nrow(significant_results) == 0) {
  cat("No Version A FDR-significant dSVs are available; 10-fold CV is skipped.\n")
} else {
  
  significant_dsv_list <- unique(significant_results$dsv)
  
  # dsv_filtered was created in Section 7 with sample IDs stored as row names.
  if (!"sample_id" %in% colnames(dsv_filtered)) {
    dsv_filtered_df <- dsv_filtered %>% rownames_to_column("sample_id")
  } else {
    dsv_filtered_df <- dsv_filtered
  }
  
  # Extract the species name using the same rule as the association analysis.
  extract_species <- function(dsv_name) {
    dsv_species <- str_extract(
      dsv_name,
      "^[A-Za-z0-9.]+\\.[a-zA-Z0-9_.]+|^[A-Za-z0-9]+_[a-zA-Z0-9_]+"
    )
    if (is.na(dsv_species)) {
      dsv_species <- str_split(dsv_name, "_")[[1]][1]
    }
    dsv_species
  }
  
  run_cv_for_dsv_classification <- function(dsv_name, data_meta, dsv_df) {
    
    if (!dsv_name %in% colnames(dsv_df)) return(NULL)
    
    dsv_status_df <- dsv_df %>%
      select(sample_id, status_val = all_of(dsv_name))
    
    dsv_species <- as.character(extract_species(dsv_name))[1]
    coverage_col <- if (!is.na(dsv_species)) {
      paste0(dsv_species, "_coverage")
    } else {
      NA_character_
    }
    
    df <- data_meta %>%
      inner_join(dsv_status_df, by = "sample_id") %>%
      mutate(raw_status = as.numeric(status_val)) %>%
      filter(raw_status %in% c(0, 1)) %>%
      mutate(
        status = factor(
          paste0("Class_", raw_status),
          levels = c("Class_0", "Class_1")
        )
      )
    
    if (!is.na(dsv_species) && dsv_species %in% colnames(df)) {
      df <- df %>% mutate(species_abundance = .data[[dsv_species]])
    } else {
      df <- df %>% mutate(species_abundance = 0)
    }
    
    if (!is.na(coverage_col) && coverage_col %in% colnames(df)) {
      df <- df %>% mutate(species_coverage = .data[[coverage_col]])
    } else {
      df <- df %>% mutate(species_coverage = 0)
    }
    
    df <- df %>%
      drop_na(
        status, ethnicity, residency, Age, Sex, BMI,
        diet, medication, urbanization,
        any_of(gen_meta_cols),
        species_abundance, species_coverage
      ) %>%
      droplevels()
    
    class_counts <- table(df$status)
    
    if (
      nrow(df) < 30 ||
      length(class_counts) < 2 ||
      any(class_counts < 5)
    ) {
      return(NULL)
    }
    
    # ------------------------------------------------------------------------
    # 12.1 McFadden R2 and Marginal Delta R2
    # ------------------------------------------------------------------------
    
    df_glm <- df %>%
      mutate(
        status_numeric = ifelse(status == "Class_1", 1, 0)
      )
    
    calc_r2 <- function(model) {
      if (is.null(model)) return(NA_real_)
      1 - (model$deviance / model$null.deviance)
    }
    
    base_vars <- c(
      "ethnicity",
      "residency",
      "Age",
      "Sex",
      "BMI",
      "diet",
      "medication",
      "urbanization",
      gen_meta_cols,
      "species_coverage",
      "species_abundance"
    )
    
    full_formula <- as.formula(
      paste(
        "status_numeric ~",
        paste(sprintf("`%s`", base_vars), collapse = " + ")
      )
    )
    
    full_glm <- tryCatch(
      glm(
        full_formula,
        data = df_glm,
        family = binomial()
      ),
      error = function(e) NULL
    )
    
    if (is.null(full_glm)) return(NULL)
    
    full_r2 <- calc_r2(full_glm)
    
    get_delta_r2 <- function(drop_var) {
      keep_vars <- setdiff(base_vars, drop_var)
      
      reduced_formula <- as.formula(
        paste(
          "status_numeric ~",
          paste(sprintf("`%s`", keep_vars), collapse = " + ")
        )
      )
      
      reduced_model <- tryCatch(
        glm(
          reduced_formula,
          data = df_glm,
          family = binomial()
        ),
        error = function(e) NULL
      )
      
      if (is.null(reduced_model)) return(NA_real_)
      
      full_r2 - calc_r2(reduced_model)
    }
    
    r2_ethnicity    <- get_delta_r2("ethnicity")
    r2_residency    <- get_delta_r2("residency")
    r2_urbanization <- get_delta_r2("urbanization")
    r2_diet         <- get_delta_r2("diet")
    r2_coverage     <- get_delta_r2("species_coverage")
    r2_species      <- get_delta_r2("species_abundance")
    r2_medication   <- get_delta_r2("medication")
    r2_age          <- get_delta_r2("Age")
    r2_sex          <- get_delta_r2("Sex")
    r2_bmi          <- get_delta_r2("BMI")
    
    r2_general <- setNames(
      lapply(gen_meta_cols, get_delta_r2),
      gen_meta_cols
    )
    
    # ------------------------------------------------------------------------
    # 12.2 10-fold cross-validation
    # ------------------------------------------------------------------------
    
    covariate_cols <- c(
      "ethnicity",
      "residency",
      "Age",
      "Sex",
      "BMI",
      "diet",
      "medication",
      "urbanization",
      gen_meta_cols,
      "species_coverage",
      "species_abundance"
    )
    
    model_data <- df %>%
      select(all_of(c("status", covariate_cols)))
    
    formula_str <- as.formula(
      paste(
        "status ~",
        paste(sprintf("`%s`", covariate_cols), collapse = " + ")
      )
    )
    
    # Keep the same global seed used by the analysis script.
    set.seed(12345)
    
    ctrl <- trainControl(
      method = "cv",
      number = 10,
      classProbs = TRUE,
      summaryFunction = twoClassSummary,
      savePredictions = "final"
    )
    
    model_rf <- tryCatch(
      {
        train(
          formula_str,
          data = model_data,
          method = "rf",
          metric = "ROC",
          trControl = ctrl,
          importance = TRUE
        )
      },
      error = function(e) {
        message("10-fold CV failed for ", dsv_name, ": ", conditionMessage(e))
        NULL
      }
    )
    
    if (is.null(model_rf)) return(NULL)
    
    cv_results <- model_rf$results %>%
      filter(mtry == model_rf$bestTune$mtry)
    
    if (nrow(cv_results) == 0) return(NULL)
    
    res_row <- data.frame(
      dsv                 = dsv_name,
      species             = dsv_species,
      N                   = nrow(model_data),
      N_Class_0           = unname(class_counts["Class_0"]),
      N_Class_1           = unname(class_counts["Class_1"]),
      
      Full_Model_R2       = round(as.numeric(full_r2), 4),
      R2_Ethnicity        = round(as.numeric(r2_ethnicity), 4),
      R2_Residency        = round(as.numeric(r2_residency), 4),
      R2_Urbanization     = round(as.numeric(r2_urbanization), 4),
      R2_Diet             = round(as.numeric(r2_diet), 4),
      R2_SpeciesCoverage  = round(as.numeric(r2_coverage), 4),
      R2_SpeciesAbundance = round(as.numeric(r2_species), 4),
      R2_Medication       = round(as.numeric(r2_medication), 4),
      R2_Age              = round(as.numeric(r2_age), 4),
      R2_Sex              = round(as.numeric(r2_sex), 4),
      R2_BMI              = round(as.numeric(r2_bmi), 4),
      
      CV_Mean_ROC         = round(cv_results$ROC[1], 4),
      CV_ROC_SD           = round(cv_results$ROCSD[1], 4),
      CV_Mean_Sens        = round(cv_results$Sens[1], 4),
      CV_Sens_SD          = round(cv_results$SensSD[1], 4),
      CV_Mean_Spec        = round(cv_results$Spec[1], 4),
      CV_Spec_SD          = round(cv_results$SpecSD[1], 4),
      
      stringsAsFactors = FALSE
    )
    
    for (g in gen_meta_cols) {
      res_row[[paste0("R2_", g)]] <- round(
        as.numeric(r2_general[[g]]),
        4
      )
    }
    
    # ------------------------------------------------------------------------
    # 12.3 Random Forest variable importance (0-100 scale)
    # ------------------------------------------------------------------------
    
    var_imp <- varImp(
      model_rf,
      scale = TRUE
    )$importance
    
    imp_vec <- if ("Overall" %in% colnames(var_imp)) {
      setNames(var_imp$Overall, rownames(var_imp))
    } else {
      setNames(rowMeans(var_imp), rownames(var_imp))
    }
    
    for (cov_name in covariate_cols) {
      
      matched_keys <- grep(
        paste0("^", cov_name),
        names(imp_vec),
        value = TRUE,
        ignore.case = TRUE
      )
      
      score <- if (length(matched_keys) > 0) {
        max(imp_vec[matched_keys], na.rm = TRUE)
      } else {
        0
      }
      
      clean_name <- case_when(
        tolower(cov_name) == "ethnicity" ~ "Imp_Ethnicity",
        tolower(cov_name) == "residency" ~ "Imp_Residency",
        tolower(cov_name) == "age" ~ "Imp_Age",
        tolower(cov_name) == "sex" ~ "Imp_Sex",
        tolower(cov_name) == "bmi" ~ "Imp_BMI",
        tolower(cov_name) == "species_coverage" ~ "Imp_SpeciesCoverage",
        tolower(cov_name) == "species_abundance" ~ "Imp_SpeciesAbundance",
        tolower(cov_name) == "diet" ~ "Imp_Diet",
        tolower(cov_name) == "medication" ~ "Imp_Medication",
        tolower(cov_name) == "urbanization" ~ "Imp_Urbanization",
        TRUE ~ paste0("Imp_", tools::toTitleCase(cov_name))
      )
      
      res_row[[clean_name]] <- round(
        as.numeric(score),
        4
      )
    }
    
    res_row
  }
  
  # --------------------------------------------------------------------------
  # 12.4 Run 10-fold CV for all Version A FDR-significant dSVs
  # --------------------------------------------------------------------------
  
  cv_results_list <- list()
  
  cat(
    sprintf(
      "Number of Version A FDR-significant dSVs for CV: %d\n",
      length(significant_dsv_list)
    )
  )
  
  for (i in seq_along(significant_dsv_list)) {
    
    dsv_name <- significant_dsv_list[i]
    
    cat(
      sprintf(
        "[%d/%d] 10-fold CV + R2: %s\n",
        i,
        length(significant_dsv_list),
        dsv_name
      )
    )
    
    res <- run_cv_for_dsv_classification(
      dsv_name,
      combined_metadata,
      dsv_filtered_df
    )
    
    if (!is.null(res)) {
      cv_results_list[[dsv_name]] <- res
    }
  }
  
  cv_results_df <- bind_rows(cv_results_list)
  
  cat(
    sprintf(
      "\n10-fold CV completed successfully for %d dSVs.\n",
      nrow(cv_results_df)
    )
  )
  
  # --------------------------------------------------------------------------
  # 12.5 Export Excel
  # --------------------------------------------------------------------------
  
  if (nrow(cv_results_df) > 0) {
    
    cv_output_file <- file.path(
      "results_ethnicity",
      "cv_10fold_results_for_significant_dsv_with_varImp.xlsx"
    )
    
    wb <- createWorkbook()
    
    addWorksheet(
      wb,
      "CV_and_Marginal_R2"
    )
    
    writeData(
      wb,
      "CV_and_Marginal_R2",
      cv_results_df
    )
    
    setColWidths(
      wb,
      "CV_and_Marginal_R2",
      cols = 1:ncol(cv_results_df),
      widths = "auto"
    )
    
    saveWorkbook(
      wb,
      cv_output_file,
      overwrite = TRUE
    )
    
    cat(
      "10-fold CV + McFadden R2 + Marginal Delta R2 results saved to: ",
      cv_output_file,
      "\n",
      sep = ""
    )
  }
}

cat("\n========== Analysis including 10-fold CV completed ==========\n")
