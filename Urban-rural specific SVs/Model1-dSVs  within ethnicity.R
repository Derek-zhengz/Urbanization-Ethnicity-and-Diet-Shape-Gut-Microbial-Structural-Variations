
# ============================================================================
# Ethnicity-stratified Rural vs Urban dSV analysis
# No dSV prevalence filtering is applied
# Species abundance and species coverage are included as covariates
# ============================================================================

set.seed(12345)

library(tidyverse)
library(broom)
library(openxlsx)
library(brglm2)
library(caret)
library(pROC)

data_dir <- "."
results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir)

# ============================================================================
# 1. Shared data preparation
# ============================================================================

cat(sprintf("[%s] Reading input files...\n", Sys.time()))

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

dsv_data <- read.csv(file.path(data_dir, "final_filtered_dsv.csv"), row.names = 1, check.names = FALSE)
colnames(dsv_data) <- make.unique(colnames(dsv_data), sep = "__dup")

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
  numeric_cols <- data %>% select(-sample_id) %>% select(where(is.numeric))
  if (ncol(numeric_cols) == 0) {
    return(data %>% select(sample_id) %>% mutate(!!paste0(covariate_name, "_total_score") := NA_real_))
  }
  data %>% select(sample_id) %>%
    mutate(!!paste0(covariate_name, "_total_score") := rowSums(numeric_cols, na.rm = TRUE))
}

reverse_ordinal <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  if (all(is.na(x_num))) return(x)
  rng <- range(x_num, na.rm = TRUE)
  rng[1] + rng[2] - x_num
}

animal_col <- names(urbanization_data)[tolower(trimws(names(urbanization_data))) == "contact with animals"]
if (length(animal_col) == 1) urbanization_data[[animal_col]] <- reverse_ordinal(urbanization_data[[animal_col]])

diet_total <- process_covariate(diet_data, "diet")
medication_total <- process_covariate(medication_data, "medication")
urbanization_total <- process_covariate(urbanization_data, "urbanization")

# ----------------------------------------------------------------------------
# General metadata: keep each variable as an individual covariate
# ----------------------------------------------------------------------------

general_metadata_cov <- general_metadata_data
general_cols_original <- setdiff(names(general_metadata_cov), "sample_id")
general_covars <- paste0("general__", make.names(general_cols_original, unique = TRUE))
names(general_metadata_cov)[match(general_cols_original, names(general_metadata_cov))] <- general_covars

general_metadata_cov <- general_metadata_cov %>%
  mutate(across(all_of(general_covars), function(x) {
    if (is.numeric(x)) x else as.factor(trimws(as.character(x)))
  }))

combined_metadata <- metadata %>%
  inner_join(diet_total, by = "sample_id") %>%
  inner_join(medication_total, by = "sample_id") %>%
  inner_join(urbanization_total, by = "sample_id") %>%
  inner_join(general_metadata_cov, by = "sample_id") %>%
  left_join(abundance_matrix, by = "sample_id") %>%
  left_join(coverage_raw, by = "sample_id") %>%
  mutate(Sex = as.factor(Sex), residency = as.factor(residency), ethnicity = as.factor(ethnicity))

cat("Samples after metadata merge:", nrow(combined_metadata), "\n")

# ============================================================================
# 2. dSV preprocessing
# No prevalence-based filtering is applied
# ============================================================================

dsv_long <- dsv_data %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(cols = -sample_id, names_to = "dsv", values_to = "status") %>%
  mutate(status = suppressWarnings(as.numeric(status)))

extract_species <- function(dsv_vec) {
  s <- str_extract(dsv_vec, "^[A-Za-z0-9_]+(?=\\.)")
  na_idx <- is.na(s)
  if (any(na_idx)) s[na_idx] <- str_extract(dsv_vec[na_idx], "^[A-Za-z0-9]+_[A-Za-z0-9]+")
  s
}

dsv_species_map <- data.frame(
  dsv = colnames(dsv_data),
  species = extract_species(colnames(dsv_data))
)

combined_data <- dsv_long %>%
  inner_join(combined_metadata, by = "sample_id") %>%
  left_join(dsv_species_map, by = "dsv")

cat("dSVs included:", ncol(dsv_data), "\n")

# ============================================================================
# 3. Initial Rural vs Urban screening within each ethnicity
# ============================================================================

screen_dsv_within_ethnicity <- function(data, target_ethnicity) {
  ethnicity_data <- data %>% filter(ethnicity == target_ethnicity)
  dsv_names <- unique(ethnicity_data$dsv)
  results <- list()
  
  cat(sprintf("[%s] Screening ethnicity %s: %d dSVs\n",
              Sys.time(), target_ethnicity, length(dsv_names)))
  pb <- txtProgressBar(min = 0, max = length(dsv_names), style = 3)
  
  for (i in seq_along(dsv_names)) {
    setTxtProgressBar(pb, i)
    dsv_name <- dsv_names[i]
    
    current_data <- ethnicity_data %>%
      filter(dsv == dsv_name) %>%
      drop_na(status, residency) %>%
      filter(status %in% c(0, 1))
    
    if (nrow(current_data) < 30 ||
        length(unique(current_data$residency)) < 2 ||
        length(unique(current_data$status)) < 2) next
    
    contingency_table <- table(dsv_status = current_data$status, residency = current_data$residency)
    if (nrow(contingency_table) < 2 || ncol(contingency_table) < 2 ||
        any(rowSums(contingency_table) == 0) || any(colSums(contingency_table) == 0)) next
    
    chi_test <- try(chisq.test(contingency_table), silent = TRUE)
    
    if (inherits(chi_test, "try-error") || any(chi_test$expected < 5)) {
      fisher_test <- try(fisher.test(contingency_table), silent = TRUE)
      if (!inherits(fisher_test, "try-error")) {
        results[[dsv_name]] <- tibble(
          ethnicity = target_ethnicity, dsv = dsv_name,
          p.value = fisher_test$p.value, test = "Fisher exact test"
        )
      }
    } else {
      results[[dsv_name]] <- tibble(
        ethnicity = target_ethnicity, dsv = dsv_name,
        p.value = chi_test$p.value, test = "Chi-squared test"
      )
    }
  }
  
  close(pb)
  bind_rows(results)
}

# ============================================================================
# 4. Adjusted logistic regression within each ethnicity
# Rural is the reference group; Urban is the comparison group
# ============================================================================

analyze_dsv_within_ethnicity <- function(data, target_ethnicity, selected_dsv) {
  ethnicity_data <- data %>%
    filter(ethnicity == target_ethnicity) %>%
    mutate(residency = factor(residency, levels = c("Rural", "Urban")))
  
  results <- list()
  
  cat(sprintf("[%s] Adjusted models for ethnicity %s: %d candidate dSVs\n",
              Sys.time(), target_ethnicity, length(selected_dsv)))
  pb <- txtProgressBar(min = 0, max = length(selected_dsv), style = 3)
  
  for (i in seq_along(selected_dsv)) {
    setTxtProgressBar(pb, i)
    dsv_name <- selected_dsv[i]
    
    current_data <- ethnicity_data %>%
      filter(dsv == dsv_name) %>%
      drop_na(status, residency, Age, Sex, BMI,
              diet_total_score, medication_total_score,
              urbanization_total_score, all_of(general_covars)) %>%
      droplevels()
    
    dsv_species <- as.character(extract_species(dsv_name))[1]
    coverage_col <- if (!is.na(dsv_species)) paste0(dsv_species, "_coverage") else NA_character_
    
    if (!is.na(dsv_species) && dsv_species %in% colnames(current_data)) {
      current_data <- current_data %>% mutate(species_abundance = .data[[dsv_species]])
      abundance_included <- TRUE
    } else {
      current_data <- current_data %>% mutate(species_abundance = NA_real_)
      abundance_included <- FALSE
    }
    
    if (!is.na(coverage_col) && coverage_col %in% colnames(current_data)) {
      current_data <- current_data %>% mutate(species_coverage = .data[[coverage_col]])
      coverage_included <- TRUE
    } else {
      current_data <- current_data %>% mutate(species_coverage = NA_real_)
      coverage_included <- FALSE
    }
    
    current_data <- current_data %>% drop_na(species_coverage) %>% droplevels()
    
    if (nrow(current_data) < 30) next
    if (length(unique(current_data$residency)) < 2) next
    if (length(unique(current_data$status)) < 2) next
    
    model_covars <- c(
      "residency", "Age", "Sex", "BMI",
      "diet_total_score", "medication_total_score",
      "urbanization_total_score", general_covars
    )
    
    if (coverage_included) model_covars <- c(model_covars, "species_coverage")
    if (abundance_included) model_covars <- c(model_covars, "species_abundance")
    
    valid_covars <- model_covars[sapply(model_covars, function(v) {
      if (!v %in% colnames(current_data)) return(FALSE)
      x <- current_data[[v]]
      if (is.factor(x) || is.character(x)) {
        length(unique(na.omit(as.character(x)))) >= 2
      } else {
        length(unique(na.omit(x))) >= 2
      }
    })]
    
    if (!"residency" %in% valid_covars) next
    
    model_formula <- as.formula(
      paste("status ~", paste(sprintf("`%s`", valid_covars), collapse = " + "))
    )
    
    model <- tryCatch({
      sep_test <- try(
        detect_separation(model_formula, data = current_data, family = binomial()),
        silent = TRUE
      )
      
      if (!inherits(sep_test, "try-error") && isTRUE(sep_test$separation)) {
        brglm(model_formula, data = current_data, family = binomial(), method = "firth")
      } else {
        glm(model_formula, data = current_data, family = binomial())
      }
    }, error = function(e) {
      message("\nModel failed for ", dsv_name, ": ", e$message)
      NULL
    })
    
    if (is.null(model)) next
    
    note <- if (inherits(model, "brglmFit")) {
      "Firth regression"
    } else {
      "Standard logistic regression"
    }
    
    model_result <- tryCatch({
      tidy(model) %>%
        filter(str_detect(term, "^residency")) %>%
        mutate(
          ethnicity = target_ethnicity,
          comparison = "Urban_vs_Rural",
          dsv = dsv_name,
          reference_group = "Rural",
          comparison_group = "Urban",
          OR = exp(estimate),
          note = note,
          species = dsv_species,
          abundance_included = abundance_included,
          coverage_included = coverage_included,
          N_samples = nrow(current_data)
        ) %>%
        select(
          ethnicity, comparison, dsv, reference_group, comparison_group,
          N_samples, estimate, std.error, statistic, p.value, OR,
          note, species, abundance_included, coverage_included
        )
    }, error = function(e) {
      message("\nResult extraction failed for ", dsv_name, ": ", e$message)
      NULL
    })
    
    if (!is.null(model_result) && nrow(model_result) > 0) {
      results[[dsv_name]] <- model_result
    }
  }
  
  close(pb)
  bind_rows(results)
}

# ============================================================================
# 5. Run Rural vs Urban analysis separately within each ethnicity
# ============================================================================

ethnicities <- levels(combined_metadata$ethnicity)
screening_results <- list()
all_results <- list()

for (eth in ethnicities) {
  eth_meta <- combined_metadata %>% filter(ethnicity == eth)
  
  if (nrow(eth_meta) == 0 ||
      length(unique(na.omit(as.character(eth_meta$residency)))) < 2) {
    cat(sprintf("Skipping %s: fewer than two residency groups available.\n", eth))
    next
  }
  
  cat("\n============================================================\n")
  cat(sprintf("Ethnicity: %s | Rural vs Urban\n", eth))
  cat("============================================================\n")
  
  screen_res <- screen_dsv_within_ethnicity(combined_data, eth)
  screening_results[[eth]] <- screen_res
  
  selected_dsv <- screen_res %>%
    filter(!is.na(p.value), p.value < 0.05) %>%
    pull(dsv)
  
  cat(sprintf("dSVs passing raw screening (P < 0.05): %d\n", length(selected_dsv)))
  if (length(selected_dsv) == 0) next
  
  current_result <- analyze_dsv_within_ethnicity(combined_data, eth, selected_dsv)
  
  if (!is.null(current_result) && nrow(current_result) > 0) {
    current_result <- current_result %>%
      mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
      arrange(p_adj)
    
    all_results[[eth]] <- current_result
    
    write.xlsx(
      current_result,
      file.path(results_dir, paste0(make.names(eth), "_Rural_vs_Urban_dsv_results.xlsx")),
      rowNames = FALSE
    )
  }
}

screening_results_df <- bind_rows(screening_results)
final_results <- bind_rows(all_results)

if (nrow(final_results) > 0) {
  significant_results <- final_results %>%
    filter(!is.na(p_adj), p_adj < 0.05) %>%
    arrange(ethnicity, p_adj)
} else {
  significant_results <- tibble()
}

write.xlsx(
  screening_results_df,
  file.path(results_dir, "all_ethnicities_Rural_vs_Urban_dsv_screening.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  final_results,
  file.path(results_dir, "all_ethnicities_Rural_vs_Urban_dsv_results.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  significant_results,
  file.path(results_dir, "all_ethnicities_significant_Rural_vs_Urban_dsv.xlsx"),
  rowNames = FALSE
)

# ============================================================================
# 6. Ten-fold CV, McFadden R2, and marginal Delta R2
# ============================================================================

run_cv_within_ethnicity <- function(target_ethnicity, dsv_name, data_meta, dsv_df) {
  dsv_status_df <- dsv_df %>% select(sample_id, status_val = all_of(dsv_name))
  dsv_species <- as.character(extract_species(dsv_name))[1]
  coverage_col <- if (!is.na(dsv_species)) paste0(dsv_species, "_coverage") else NA_character_
  
  df <- data_meta %>%
    filter(ethnicity == target_ethnicity) %>%
    inner_join(dsv_status_df, by = "sample_id") %>%
    mutate(
      raw_status = suppressWarnings(as.numeric(status_val)),
      residency = factor(residency, levels = c("Rural", "Urban"))
    ) %>%
    filter(raw_status %in% c(0, 1)) %>%
    mutate(status = factor(paste0("Class_", raw_status), levels = c("Class_0", "Class_1")))
  
  if (!is.na(dsv_species) && dsv_species %in% colnames(df)) {
    df$species_abundance <- df[[dsv_species]]
  } else {
    df$species_abundance <- NA_real_
  }
  
  if (!is.na(coverage_col) && coverage_col %in% colnames(df)) {
    df$species_coverage <- df[[coverage_col]]
  } else {
    df$species_coverage <- NA_real_
  }
  
  df <- df %>%
    drop_na(status, residency, Age, Sex, BMI,
            diet_total_score, medication_total_score,
            urbanization_total_score, all_of(general_covars),
            species_abundance, species_coverage) %>%
    droplevels()
  
  class_counts <- table(df$status)
  if (nrow(df) < 30 || length(class_counts) < 2 || any(class_counts < 5)) return(NULL)
  if (length(unique(df$residency)) < 2) return(NULL)
  
  df_glm <- df %>% mutate(status_numeric = ifelse(status == "Class_1", 1, 0))
  
  base_vars <- c(
    "residency", "Age", "Sex", "BMI",
    "diet_total_score", "medication_total_score",
    "urbanization_total_score", general_covars,
    "species_coverage", "species_abundance"
  )
  
  valid_vars <- base_vars[sapply(base_vars, function(v) {
    if (!v %in% colnames(df_glm)) return(FALSE)
    x <- df_glm[[v]]
    if (is.factor(x) || is.character(x)) {
      length(unique(na.omit(as.character(x)))) >= 2
    } else {
      length(unique(na.omit(x))) >= 2
    }
  })]
  
  if (!"residency" %in% valid_vars) return(NULL)
  
  calc_r2 <- function(m) {
    if (is.null(m)) return(NA_real_)
    1 - (m$deviance / m$null.deviance)
  }
  
  full_formula <- as.formula(
    paste("status_numeric ~", paste(sprintf("`%s`", valid_vars), collapse = " + "))
  )
  
  full_glm <- tryCatch(
    glm(full_formula, data = df_glm, family = binomial()),
    error = function(e) NULL
  )
  if (is.null(full_glm)) return(NULL)
  
  full_r2 <- calc_r2(full_glm)
  
  get_delta_r2 <- function(drop_var) {
    if (!drop_var %in% valid_vars) return(0)
    keep_vars <- setdiff(valid_vars, drop_var)
    if (length(keep_vars) == 0) return(full_r2)
    
    reduced_formula <- as.formula(
      paste("status_numeric ~", paste(sprintf("`%s`", keep_vars), collapse = " + "))
    )
    
    reduced_model <- tryCatch(
      glm(reduced_formula, data = df_glm, family = binomial()),
      error = function(e) NULL
    )
    
    if (is.null(reduced_model)) return(NA_real_)
    full_r2 - calc_r2(reduced_model)
  }
  
  r2_values <- setNames(lapply(valid_vars, get_delta_r2), valid_vars)
  
  model_data <- df %>% select(status, all_of(valid_vars))
  ctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)
  
  model_rf <- tryCatch(
    train(status ~ ., data = model_data, method = "rf",
          metric = "ROC", trControl = ctrl, importance = TRUE),
    error = function(e) NULL
  )
  if (is.null(model_rf)) return(NULL)
  
  cv_results <- model_rf$results %>% filter(mtry == model_rf$bestTune$mtry)
  
  result <- data.frame(
    ethnicity = target_ethnicity,
    comparison = "Rural_vs_Urban",
    dsv = dsv_name,
    species = dsv_species,
    N = nrow(model_data),
    Full_Model_R2 = round(full_r2, 4),
    CV_Mean_ROC = round(cv_results$ROC, 4),
    CV_ROC_SD = round(cv_results$ROCSD, 4),
    CV_Mean_Sens = round(cv_results$Sens, 4),
    CV_Mean_Spec = round(cv_results$Spec, 4),
    stringsAsFactors = FALSE
  )
  
  for (nm in names(r2_values)) {
    result[[paste0("R2_", make.names(nm))]] <- round(as.numeric(r2_values[[nm]]), 4)
  }
  
  result
}

if (nrow(significant_results) > 0) {
  dsv_df <- dsv_data %>% rownames_to_column("sample_id")
  cv_list <- vector("list", nrow(significant_results))
  
  cat("\nRunning 10-fold CV for significant ethnicity-specific Rural vs Urban dSV associations...\n")
  pb <- txtProgressBar(min = 0, max = nrow(significant_results), style = 3)
  
  for (i in seq_len(nrow(significant_results))) {
    cv_list[[i]] <- run_cv_within_ethnicity(
      significant_results$ethnicity[i],
      significant_results$dsv[i],
      combined_metadata,
      dsv_df
    )
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  cv_results <- bind_rows(cv_list)
  
  if (nrow(cv_results) > 0) {
    write.xlsx(
      cv_results,
      file.path(results_dir, "ethnicity_specific_Rural_vs_Urban_dsv_cv_marginal_r2.xlsx"),
      rowNames = FALSE
    )
  }
}

cat("\nAnalysis completed.\n")
cat("Ethnicities with valid Rural/Urban results:", length(all_results), "\n")
cat("Significant ethnicity-specific Rural vs Urban dSV associations:", nrow(significant_results), "\n")