# ============================================================================
# Residency-associated dSV analysis
# Rural vs Urban comparison across the full cohort
# Species abundance and species coverage are included as covariates
# ============================================================================

set.seed(12345)
library(tidyverse)
library(broom)
library(openxlsx)
library(brglm2)
library(caret)
library(pROC)


# ============================================================================
# 0. 

# ============================================================================

cat("========== 1.  ==========\n")

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

dsv_data <- read.csv("final_filtered_dsv.csv", row.names = 1, check.names = FALSE)
colnames(dsv_data) <- make.unique(colnames(dsv_data), sep = "__dup")
metadata <- read.csv("metadata.csv")
diet_data <- read.csv("final_filtered_diet.csv", check.names = FALSE)
urbanization_data <- read.csv("Urbanization.csv", check.names = FALSE)
medication_data <- read.csv("Medication.csv", check.names = FALSE)
general_metadata_data <- read.csv("General metadata.csv", check.names = FALSE)

required_meta_cols <- c("sample_id", "residency", "Age", "Sex", "BMI", "ethnicity")
missing_cols <- setdiff(required_meta_cols, names(metadata))
if (length(missing_cols) > 0) {
  stop(paste("metadata :", paste(missing_cols, collapse = ", ")))
}

process_covariate <- function(data, covariate_name) {
  numeric_cols <- data %>% select(-sample_id) %>% select(where(is.numeric))
  if (ncol(numeric_cols) == 0) {
    return(data %>% select(sample_id) %>% mutate(!!paste0(covariate_name, "_total_score") := NA_real_))
  }
  result <- data %>% select(sample_id) %>% mutate(!!paste0(covariate_name, "_total_score") := rowSums(numeric_cols, na.rm = TRUE))
  return(result)
}

diet_total <- process_covariate(diet_data, "diet")
medication_total <- process_covariate(medication_data, "medication")
urbanization_total <- process_covariate(urbanization_data, "urbanization")

# General metadata: keep each variable as an individual covariate
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
  mutate(
    Sex = as.factor(Sex),
    residency = as.factor(residency),
    ethnicity = as.factor(ethnicity)
  )

cat("Samples after metadata merge:", nrow(combined_metadata), "\n")

# ============================================================================
# 2. dSV  ( )

# ============================================================================

dsv_data_original <- dsv_data
dsv_data_prevalence <- dsv_data %>% mutate(across(everything(), ~replace(., is.na(.), 0)))

dsv_prevalence <- data.frame(
  dsv = colnames(dsv_data_prevalence),
  count_1 = colSums(dsv_data_prevalence == 1, na.rm = TRUE),
  count_0 = colSums(dsv_data_prevalence == 0, na.rm = TRUE),
  prevalence = colMeans(dsv_data_prevalence == 1, na.rm = TRUE)
)

valid_dsv <- dsv_prevalence %>%
  filter(prevalence >= 0.1) %>%
  filter(count_1 >= 10 & count_0 >= 10) %>%
  pull(dsv)

dsv_filtered <- dsv_data_original[, valid_dsv, drop = FALSE]

dsv_long <- dsv_filtered %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(cols = -sample_id, names_to = "dsv", values_to = "status") %>%
  mutate(status = as.numeric(status))

#  ：
extract_species <- function(dsv_vec) {
  s <- str_extract(dsv_vec, "^[A-Za-z0-9_]+(?=\\.)")
  na_idx <- is.na(s)
  if (any(na_idx)) {
    s[na_idx] <- str_extract(dsv_vec[na_idx], "^[A-Za-z0-9]+_[A-Za-z0-9]+")
  }
  return(s)
}

dsv_species_map <- data.frame(
  dsv = colnames(dsv_filtered),
  species = extract_species(colnames(dsv_filtered))
)

combined_data <- dsv_long %>%
  inner_join(combined_metadata, by = "sample_id") %>%
  left_join(dsv_species_map, by = "dsv")

cat(" dSV :", ncol(dsv_filtered), "\n")


# ============================================================================
# 3.  dSV (/Fisher)

# ============================================================================

screen_dsv_for_residency <- function(data) {
  results <- list()
  dsv_names <- unique(data$dsv)
  
  for (i in seq_along(dsv_names)) {
    dsv <- dsv_names[i]
    current_data <- filter(data, dsv == !!dsv) %>% drop_na(status, residency)
    
    if (nrow(current_data) < 30 || length(unique(current_data$residency)) < 2) next
    contingency_table <- table(dsv_status = current_data$status, residency = current_data$residency)
    if (any(rowSums(contingency_table) == 0) || any(colSums(contingency_table) == 0)) next
    
    chi_test <- try(chisq.test(contingency_table), silent = TRUE)
    if (inherits(chi_test, "try-error") || any(chi_test$expected < 5)) {
      fisher_test <- try(fisher.test(contingency_table), silent = TRUE)
      if (!inherits(fisher_test, "try-error")) {
        results[[dsv]] <- tibble(dsv = dsv, p.value = fisher_test$p.value)
      }
    } else {
      results[[dsv]] <- tibble(dsv = dsv, p.value = chi_test$p.value)
    }
  }
  
  all_res <- bind_rows(results)
  selected_dsv <- all_res %>% filter(!is.na(p.value), p.value < 0.05) %>% pull(dsv)
  cat(sprintf("/Fisher  p < 0.05  dSV : %d\n", length(selected_dsv)))
  return(selected_dsv)
}


# ============================================================================
# 4. Logistic 

# ============================================================================

analyze_residency_both_references <- function(data, selected_dsv) {
  all_results <- list()
  reference_levels <- list(
    Rural_ref = c("Rural", "Urban"),
    Urban_ref = c("Urban", "Rural")
  )
  
  for (ref_name in names(reference_levels)) {
    ref_levels <- reference_levels[[ref_name]]
    reference_group <- ref_levels[1]
    comparison_group <- ref_levels[2]
    
    data_ref <- data %>%
      mutate(residency = factor(residency, levels = ref_levels))
    
    results <- list()
    
    cat(sprintf("\nRunning %s (%s vs %s): %d dSVs\n",
                ref_name, comparison_group, reference_group, length(selected_dsv)))
    
    pb <- txtProgressBar(min = 0, max = length(selected_dsv), style = 3)
    
    for (i in seq_along(selected_dsv)) {
      setTxtProgressBar(pb, i)
      
      dsv <- selected_dsv[i]
      
      current_data <- filter(data_ref, dsv == !!dsv) %>%
        drop_na(
          status, residency, Age, Sex, BMI, ethnicity,
          diet_total_score, medication_total_score,
          urbanization_total_score, all_of(general_covars)
        ) %>%
        droplevels()
      
      dsv_species <- as.character(extract_species(dsv))[1]
      coverage_col <- if (!is.na(dsv_species)) paste0(dsv_species, "_coverage") else NA_character_
      
      if (!is.na(dsv_species) && dsv_species %in% colnames(current_data)) {
        current_data <- current_data %>%
          mutate(species_abundance = .data[[dsv_species]])
        abundance_included <- TRUE
      } else {
        current_data <- current_data %>%
          mutate(species_abundance = NA_real_)
        abundance_included <- FALSE
      }
      
      if (!is.na(coverage_col) && coverage_col %in% colnames(current_data)) {
        current_data <- current_data %>%
          mutate(species_coverage = .data[[coverage_col]])
        coverage_included <- TRUE
      } else {
        current_data <- current_data %>%
          mutate(species_coverage = NA_real_)
        coverage_included <- FALSE
      }
      
      current_data <- current_data %>%
        drop_na(species_coverage) %>%
        droplevels()
      
      if (nrow(current_data) < 30) next
      if (length(unique(current_data$residency)) < 2) next
      if (length(unique(current_data$status)) < 2) next
      
      model_covars <- c(
        "residency", "Age", "Sex", "BMI", "ethnicity",
        "diet_total_score", "medication_total_score",
        "urbanization_total_score",
        general_covars
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
          detect_separation(
            model_formula,
            data = current_data,
            family = binomial()
          ),
          silent = TRUE
        )
        
        if (!inherits(sep_test, "try-error") && isTRUE(sep_test$separation)) {
          brglm(
            model_formula,
            data = current_data,
            family = binomial(),
            method = "firth"
          )
        } else {
          glm(
            model_formula,
            data = current_data,
            family = binomial()
          )
        }
      }, error = function(e) {
        message("\nModel failed for ", dsv, ": ", e$message)
        NULL
      })
      
      if (is.null(model)) next
      
      note <- if (inherits(model, "brglmFit")) {
        paste0(
          "Firth regression",
          ifelse(abundance_included, " + abundance", " (without abundance)"),
          ifelse(coverage_included, " + coverage", " (without coverage)")
        )
      } else {
        paste0(
          "Standard logistic regression",
          ifelse(abundance_included, " + abundance", " (without abundance)"),
          ifelse(coverage_included, " + coverage", " (without coverage)")
        )
      }
      
      tryCatch({
        model_result <- tidy(model) %>%
          filter(str_detect(term, "^residency")) %>%
          mutate(
            dsv = dsv,
            reference_group = reference_group,
            comparison_group = comparison_group,
            or = exp(estimate),
            note = note,
            term = str_remove(term, "residency"),
            species = dsv_species,
            abundance_included = abundance_included,
            coverage_included = coverage_included
          ) %>%
          select(
            dsv, reference_group, comparison_group, term,
            estimate, std.error, statistic, p.value, or,
            note, species, abundance_included, coverage_included
          )
        
        if (nrow(model_result) > 0) {
          results[[dsv]] <- model_result
        }
      }, error = function(e) {
        message("\nResult extraction failed for ", dsv, ": ", e$message)
      })
    }
    
    close(pb)
    
    current_res <- bind_rows(results)
    
    if (nrow(current_res) > 0) {
      current_res <- current_res %>%
        mutate(p_adj = p.adjust(p.value, method = "BH"))
    }
    
    all_results[[ref_name]] <- current_res
    
    cat(sprintf(
      "\n%s completed: %d result rows generated.\n",
      ref_name, nrow(current_res)
    ))
  }
  
  bind_rows(all_results)
}

selected_dsv <- screen_dsv_for_residency(combined_data)

if (!dir.exists("results")) dir.create("results")

if (length(selected_dsv) > 0) {
  residency_results <- analyze_residency_both_references(combined_data, selected_dsv)
  write.xlsx(residency_results, "results/residency_dsv_logistic_results_with_abundance.xlsx", rowNames = FALSE)
  
  significant_results <- residency_results %>% filter(!is.na(p_adj), p_adj < 0.05)
  write.xlsx(significant_results, "results/significant_residency_dsv_results_with_abundance.xlsx", rowNames = FALSE)
  cat("Logistic ， dSV :", length(unique(significant_results$dsv)), "\n")
} else {
  stop(" dSV，。")
}


# ============================================================================
# 5. 10  CV 、 R2  Marginal Delta R2  ( )

# ============================================================================

cat("\n========== 5.  10  CV 、 R2  R2  ==========\n")

significant_dsv_list <- unique(significant_results$dsv)

if (!"sample_id" %in% colnames(dsv_filtered)) {
  dsv_filtered_df <- dsv_filtered %>% rownames_to_column("sample_id")
} else {
  dsv_filtered_df <- dsv_filtered
}

run_cv_for_dsv_classification <- function(dsv_name, data_meta, dsv_df) {
  
  dsv_status_df <- dsv_df %>% select(sample_id, status_val = !!sym(dsv_name))
  dsv_species <- as.character(extract_species(dsv_name))[1]
  coverage_col <- if (!is.na(dsv_species)) paste0(dsv_species, "_coverage") else NA_character_
  
  df <- data_meta %>%
    inner_join(dsv_status_df, by = "sample_id") %>%
    mutate(raw_status = as.numeric(status_val)) %>%
    filter(raw_status %in% c(0, 1)) %>%
    mutate(status = factor(paste0("Class_", raw_status), levels = c("Class_0", "Class_1")))
  
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
    drop_na(status, residency, Age, Sex, BMI, ethnicity,
            diet_total_score, medication_total_score,
            urbanization_total_score, all_of(general_covars),
            species_abundance, species_coverage)
  
  class_counts <- table(df$status)
  if (nrow(df) < 30 || length(class_counts) < 2 || any(class_counts < 5)) return(NULL)
  
  # --------------------------------------------------------------------------
  #  ： McFadden R2  Marginal Delta R2
  # --------------------------------------------------------------------------
  
  df_glm <- df %>% mutate(status_numeric = ifelse(status == "Class_1", 1, 0))
  
  calc_r2 <- function(m) {
    if (is.null(m)) return(NA_real_)
    1 - (m$deviance / m$null.deviance)
  }
  
  base_vars <- c(
    "ethnicity", "residency", "Age", "Sex", "BMI",
    "diet_total_score", "medication_total_score", "urbanization_total_score",
    general_covars, "species_coverage", "species_abundance"
  )
  
  full_formula <- as.formula(
    paste("status_numeric ~", paste(sprintf("`%s`", base_vars), collapse = " + "))
  )
  
  full_glm <- tryCatch(
    glm(full_formula, data = df_glm, family = binomial()),
    error = function(e) NULL
  )
  if (is.null(full_glm)) return(NULL)
  full_r2 <- calc_r2(full_glm)
  
  get_delta_r2 <- function(drop_var) {
    keep_vars <- setdiff(base_vars, drop_var)
    f <- as.formula(paste("status_numeric ~", paste(keep_vars, collapse = " + ")))
    m_red <- tryCatch(glm(f, data = df_glm, family = binomial()), error = function(e) NULL)
    if (is.null(m_red)) return(NA_real_)
    return(full_r2 - calc_r2(m_red))
  }
  
  #  R2 
  r2_ethnicity    <- get_delta_r2("ethnicity")
  r2_residency    <- get_delta_r2("residency")
  r2_urbanization <- get_delta_r2("urbanization_total_score")
  r2_diet         <- get_delta_r2("diet_total_score")
  r2_coverage     <- get_delta_r2("species_coverage")
  r2_species      <- get_delta_r2("species_abundance")
  r2_medication   <- get_delta_r2("medication_total_score")
  r2_age          <- get_delta_r2("Age")
  r2_bmi          <- get_delta_r2("BMI")
  r2_general <- setNames(lapply(general_covars, get_delta_r2), general_covars)
  
  # --------------------------------------------------------------------------
  # 10  CV 
  # --------------------------------------------------------------------------
  
  covariate_cols <- c(
    "ethnicity", "residency", "Age", "Sex", "BMI",
    "diet_total_score", "medication_total_score", "urbanization_total_score",
    general_covars, "species_coverage", "species_abundance"
  )
  
  model_data <- df %>% select(all_of(c("status", covariate_cols)))
  formula_str <- as.formula(paste("status ~", paste(covariate_cols, collapse = " + ")))
  
  ctrl <- trainControl(
    method = "cv",
    number = 10,
    classProbs = TRUE,
    summaryFunction = twoClassSummary
  )
  
  model_rf <- tryCatch({
    train(formula_str, data = model_data, method = "rf", metric = "ROC", trControl = ctrl, importance = TRUE)
  }, error = function(e) return(NULL))
  
  if (is.null(model_rf)) return(NULL)
  
  cv_results <- model_rf$results %>% filter(mtry == model_rf$bestTune$mtry)
  
  res_row <- data.frame(
    dsv                   = dsv_name,
    species               = dsv_species,
    N                     = nrow(model_data),
    
    #   R2 
    Full_Model_R2         = round(as.numeric(full_r2), 4),
    R2_Ethnicity          = round(as.numeric(r2_ethnicity), 4),
    R2_Residency          = round(as.numeric(r2_residency), 4),
    R2_Urbanization       = round(as.numeric(r2_urbanization), 4),
    R2_Diet               = round(as.numeric(r2_diet), 4),
    R2_SpeciesCoverage    = round(as.numeric(r2_coverage), 4),
    R2_SpeciesAbundance   = round(as.numeric(r2_species), 4),
    R2_Medication         = round(as.numeric(r2_medication), 4),
    R2_Age                = round(as.numeric(r2_age), 4),
    R2_BMI                = round(as.numeric(r2_bmi), 4),
    
    # CV 
    CV_Mean_ROC           = round(cv_results$ROC, 4),
    CV_ROC_SD             = round(cv_results$ROCSD, 4),
    CV_Mean_Sens          = round(cv_results$Sens, 4),
    CV_Mean_Spec          = round(cv_results$Spec, 4),
    stringsAsFactors      = FALSE
  )
  
  for (g in general_covars) {
    res_row[[paste0("R2_", g)]] <- round(as.numeric(r2_general[[g]]), 4)
  }
  
  #  (0-100)
  var_imp <- varImp(model_rf, scale = TRUE)$importance
  imp_vec <- if ("Overall" %in% colnames(var_imp)) setNames(var_imp$Overall, rownames(var_imp)) else setNames(rowMeans(var_imp), rownames(var_imp))
  
  for (cov_name in covariate_cols) {
    matched_keys <- grep(paste0("^", cov_name), names(imp_vec), value = TRUE, ignore.case = TRUE)
    score <- if (length(matched_keys) > 0) max(imp_vec[matched_keys], na.rm = TRUE) else 0
    clean_name <- case_when(
      tolower(cov_name) == "ethnicity" ~ "Imp_Ethnicity",
      tolower(cov_name) == "residency" ~ "Imp_Residency",
      tolower(cov_name) == "age" ~ "Imp_Age",
      tolower(cov_name) == "sex" ~ "Imp_Sex",
      tolower(cov_name) == "bmi" ~ "Imp_BMI",
      tolower(cov_name) == "species_coverage" ~ "Imp_SpeciesCoverage",
      tolower(cov_name) == "species_abundance" ~ "Imp_SpeciesAbundance",
      tolower(cov_name) == "diet_total_score" ~ "Imp_Diet",
      tolower(cov_name) == "medication_total_score" ~ "Imp_Medication",
      tolower(cov_name) == "urbanization_total_score" ~ "Imp_Urbanization",
      TRUE ~ paste0("Imp_", tools::toTitleCase(cov_name))
    )
    res_row[[clean_name]] <- round(score, 4)
  }
  
  return(res_row)
}

#  10  CV 
cv_results_list <- list()
cat(sprintf(" dSV : %d\n", length(significant_dsv_list)))

for (i in seq_along(significant_dsv_list)) {
  dsv <- significant_dsv_list[i]
  cat(sprintf("[%d/%d]  R2: %s\n", i, length(significant_dsv_list), dsv))
  
  res <- run_cv_for_dsv_classification(dsv, combined_metadata, dsv_filtered_df)
  if (!is.null(res)) cv_results_list[[dsv]] <- res
}

cv_results_df <- bind_rows(cv_results_list)

cat(sprintf("\n ， %d ！\n", nrow(cv_results_df)))

#   Excel 
if (nrow(cv_results_df) > 0) {
  wb <- createWorkbook()
  addWorksheet(wb, "CV_and_Marginal_R2")
  writeData(wb, "CV_and_Marginal_R2", cv_results_df)
  setColWidths(wb, "CV_and_Marginal_R2", cols = 1:ncol(cv_results_df), widths = "auto")
  
  saveWorkbook(wb, "results/cv_10fold_results_for_significant_dsv_with_varImp.xlsx", overwrite = TRUE)
  cat("  R2  Marginal Delta R2  Excel : results/cv_10fold_results_for_significant_dsv_with_varImp.xlsx\n")
}


# ============================================================================
# VERSION B
# Sensitivity analysis restricted to species coverage >= 5X
# Coverage remains included as a covariate
# ============================================================================

cat("\n================ VERSION B ================\n")

screen_dsv_for_residency_verB <- function(data) {
  results <- list()
  dsv_names <- unique(data$dsv)
  
  cat(sprintf("Version B screening: %d dSVs\n", length(dsv_names)))
  pb <- txtProgressBar(min = 0, max = length(dsv_names), style = 3)
  
  for (i in seq_along(dsv_names)) {
    setTxtProgressBar(pb, i)
    dsv_name <- dsv_names[i]
    dsv_species <- as.character(extract_species(dsv_name))[1]
    coverage_col <- if (!is.na(dsv_species)) paste0(dsv_species, "_coverage") else NA_character_
    
    current_data <- data %>% filter(dsv == dsv_name)
    
    if (!is.na(coverage_col) && coverage_col %in% colnames(current_data)) {
      current_data <- current_data %>%
        mutate(species_coverage = .data[[coverage_col]])
    } else {
      next
    }
    
    current_data <- current_data %>%
      drop_na(status, residency, species_coverage) %>%
      filter(status %in% c(0, 1), species_coverage >= 5)
    
    if (nrow(current_data) < 30 ||
        length(unique(current_data$residency)) < 2 ||
        length(unique(current_data$status)) < 2) next
    
    contingency_table <- table(
      dsv_status = current_data$status,
      residency = current_data$residency
    )
    
    if (nrow(contingency_table) < 2 || ncol(contingency_table) < 2 ||
        any(rowSums(contingency_table) == 0) ||
        any(colSums(contingency_table) == 0)) next
    
    chi_test <- try(chisq.test(contingency_table), silent = TRUE)
    
    if (inherits(chi_test, "try-error") || any(chi_test$expected < 5)) {
      fisher_test <- try(fisher.test(contingency_table), silent = TRUE)
      if (!inherits(fisher_test, "try-error")) {
        results[[dsv_name]] <- tibble(
          dsv = dsv_name,
          p.value = fisher_test$p.value,
          test = "Fisher exact test",
          N_cov5x = nrow(current_data)
        )
      }
    } else {
      results[[dsv_name]] <- tibble(
        dsv = dsv_name,
        p.value = chi_test$p.value,
        test = "Chi-squared test",
        N_cov5x = nrow(current_data)
      )
    }
  }
  
  close(pb)
  
  all_res <- bind_rows(results)
  selected_dsv <- all_res %>%
    filter(!is.na(p.value), p.value < 0.05) %>%
    pull(dsv)
  
  cat(sprintf("Version B screening P < 0.05: %d dSVs\n", length(selected_dsv)))
  selected_dsv
}


analyze_residency_both_references_verB <- function(data, selected_dsv) {
  all_results <- list()
  reference_levels <- list(
    Rural_ref = c("Rural", "Urban"),
    Urban_ref = c("Urban", "Rural")
  )
  
  for (ref_name in names(reference_levels)) {
    ref_levels <- reference_levels[[ref_name]]
    reference_group <- ref_levels[1]
    comparison_group <- ref_levels[2]
    
    data_ref <- data %>%
      mutate(residency = factor(residency, levels = ref_levels))
    
    results <- list()
    
    cat(sprintf("\nVersion B %s (%s vs %s): %d dSVs\n",
                ref_name, comparison_group, reference_group, length(selected_dsv)))
    pb <- txtProgressBar(min = 0, max = length(selected_dsv), style = 3)
    
    for (i in seq_along(selected_dsv)) {
      setTxtProgressBar(pb, i)
      dsv_name <- selected_dsv[i]
      
      current_data <- data_ref %>%
        filter(dsv == dsv_name) %>%
        drop_na(status, residency, Age, Sex, BMI, ethnicity,
                diet_total_score, medication_total_score,
                urbanization_total_score, all_of(general_covars)) %>%
        droplevels()
      
      dsv_species <- as.character(extract_species(dsv_name))[1]
      coverage_col <- if (!is.na(dsv_species)) paste0(dsv_species, "_coverage") else NA_character_
      
      if (!is.na(dsv_species) && dsv_species %in% colnames(current_data)) {
        current_data <- current_data %>%
          mutate(species_abundance = .data[[dsv_species]])
        abundance_included <- TRUE
      } else {
        current_data <- current_data %>%
          mutate(species_abundance = NA_real_)
        abundance_included <- FALSE
      }
      
      if (!is.na(coverage_col) && coverage_col %in% colnames(current_data)) {
        current_data <- current_data %>%
          mutate(species_coverage = .data[[coverage_col]])
      } else {
        next
      }
      
      current_data <- current_data %>%
        drop_na(species_coverage) %>%
        filter(species_coverage >= 5) %>%
        droplevels()
      
      if (nrow(current_data) < 30) next
      if (length(unique(current_data$residency)) < 2) next
      if (length(unique(current_data$status)) < 2) next
      
      model_covars <- c(
        "residency", "Age", "Sex", "BMI", "ethnicity",
        "diet_total_score", "medication_total_score",
        "urbanization_total_score", general_covars,
        "species_coverage"
      )
      
      if (abundance_included) {
        model_covars <- c(model_covars, "species_abundance")
      }
      
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
          brglm(
            model_formula,
            data = current_data,
            family = binomial(),
            method = "firth"
          )
        } else {
          glm(
            model_formula,
            data = current_data,
            family = binomial()
          )
        }
      }, error = function(e) {
        message("\nVersion B model failed for ", dsv_name, ": ", e$message)
        NULL
      })
      
      if (is.null(model)) next
      
      note <- if (inherits(model, "brglmFit")) {
        "Firth regression + coverage >= 5X"
      } else {
        "Standard logistic regression + coverage >= 5X"
      }
      
      model_result <- tryCatch({
        tidy(model) %>%
          filter(str_detect(term, "^residency")) %>%
          mutate(
            dsv = dsv_name,
            reference_group = reference_group,
            comparison_group = comparison_group,
            OR = exp(estimate),
            note = note,
            species = dsv_species,
            abundance_included = abundance_included,
            coverage_included = TRUE,
            N_cov5x = nrow(current_data)
          ) %>%
          select(
            dsv, reference_group, comparison_group, N_cov5x,
            estimate, std.error, statistic, p.value, OR,
            note, species, abundance_included, coverage_included
          )
      }, error = function(e) {
        message("\nVersion B result extraction failed for ", dsv_name, ": ", e$message)
        NULL
      })
      
      if (!is.null(model_result) && nrow(model_result) > 0) {
        results[[dsv_name]] <- model_result
      }
    }
    
    close(pb)
    current_res <- bind_rows(results)
    
    if (nrow(current_res) > 0) {
      current_res <- current_res %>%
        mutate(p_adj = p.adjust(p.value, method = "BH"))
    }
    
    all_results[[ref_name]] <- current_res
  }
  
  bind_rows(all_results)
}


selected_dsv_verB <- screen_dsv_for_residency_verB(combined_data)

if (length(selected_dsv_verB) > 0) {
  residency_results_verB <- analyze_residency_both_references_verB(
    combined_data,
    selected_dsv_verB
  )
  
  write.xlsx(
    residency_results_verB,
    file.path(results_dir, "residency_dsv_results_verB_cov5x.xlsx"),
    rowNames = FALSE
  )
  
  significant_results_verB <- residency_results_verB %>%
    filter(!is.na(p_adj), p_adj < 0.05)
  
  write.xlsx(
    significant_results_verB,
    file.path(results_dir, "significant_residency_dsv_results_verB_cov5x.xlsx"),
    rowNames = FALSE
  )
  
  cat("Version B significant dSVs after FDR:",
      length(unique(significant_results_verB$dsv)), "\n")
} else {
  residency_results_verB <- tibble()
  significant_results_verB <- tibble()
  cat("No dSVs passed Version B screening.\n")
}


# ----------------------------------------------------------------------------
# Version B: ten-fold CV, McFadden R2, and marginal Delta R2
# ----------------------------------------------------------------------------

run_cv_for_dsv_classification_verB <- function(dsv_name, data_meta, dsv_df) {
  dsv_status_df <- dsv_df %>%
    select(sample_id, status_val = all_of(dsv_name))
  
  dsv_species <- as.character(extract_species(dsv_name))[1]
  coverage_col <- if (!is.na(dsv_species)) paste0(dsv_species, "_coverage") else NA_character_
  
  df <- data_meta %>%
    inner_join(dsv_status_df, by = "sample_id") %>%
    mutate(raw_status = suppressWarnings(as.numeric(status_val))) %>%
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
    return(NULL)
  }
  
  df <- df %>%
    drop_na(status, residency, Age, Sex, BMI, ethnicity,
            diet_total_score, medication_total_score,
            urbanization_total_score, all_of(general_covars),
            species_abundance, species_coverage) %>%
    filter(species_coverage >= 5) %>%
    droplevels()
  
  class_counts <- table(df$status)
  if (nrow(df) < 30 || length(class_counts) < 2 || any(class_counts < 5)) return(NULL)
  
  base_vars <- c(
    "ethnicity", "residency", "Age", "Sex", "BMI",
    "diet_total_score", "medication_total_score",
    "urbanization_total_score", general_covars,
    "species_coverage", "species_abundance"
  )
  
  valid_vars <- base_vars[sapply(base_vars, function(v) {
    if (!v %in% colnames(df)) return(FALSE)
    x <- df[[v]]
    
    if (is.factor(x) || is.character(x)) {
      length(unique(na.omit(as.character(x)))) >= 2
    } else {
      length(unique(na.omit(x))) >= 2
    }
  })]
  
  df_glm <- df %>%
    mutate(status_numeric = ifelse(status == "Class_1", 1, 0))
  
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
  
  model_data <- df %>%
    select(status, all_of(valid_vars))
  
  ctrl <- trainControl(
    method = "cv",
    number = 10,
    classProbs = TRUE,
    summaryFunction = twoClassSummary
  )
  
  model_rf <- tryCatch(
    train(
      status ~ .,
      data = model_data,
      method = "rf",
      metric = "ROC",
      trControl = ctrl,
      importance = TRUE
    ),
    error = function(e) NULL
  )
  
  if (is.null(model_rf)) return(NULL)
  
  cv_results <- model_rf$results %>%
    filter(mtry == model_rf$bestTune$mtry)
  
  result <- data.frame(
    dsv = dsv_name,
    species = dsv_species,
    N_cov5x = nrow(model_data),
    Full_Model_R2 = round(full_r2, 4),
    CV_Mean_ROC = round(cv_results$ROC, 4),
    CV_ROC_SD = round(cv_results$ROCSD, 4),
    CV_Mean_Sens = round(cv_results$Sens, 4),
    CV_Mean_Spec = round(cv_results$Spec, 4),
    stringsAsFactors = FALSE
  )
  
  for (nm in names(r2_values)) {
    result[[paste0("R2_", make.names(nm))]] <- round(
      as.numeric(r2_values[[nm]]),
      4
    )
  }
  
  result
}


if (nrow(significant_results_verB) > 0) {
  dsv_verB_df <- dsv_filtered %>%
    rownames_to_column("sample_id")
  
  significant_dsv_verB <- unique(significant_results_verB$dsv)
  cv_verB_list <- vector("list", length(significant_dsv_verB))
  
  cat(sprintf("Version B CV: %d significant dSVs\n", length(significant_dsv_verB)))
  pb <- txtProgressBar(min = 0, max = length(significant_dsv_verB), style = 3)
  
  for (i in seq_along(significant_dsv_verB)) {
    cv_verB_list[[i]] <- run_cv_for_dsv_classification_verB(
      significant_dsv_verB[i],
      combined_metadata,
      dsv_verB_df
    )
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  cv_verB <- bind_rows(cv_verB_list)
  
  if (nrow(cv_verB) > 0) {
    write.xlsx(
      cv_verB,
      file.path(results_dir, "dsv_verB_cov5x_cv_marginal_r2.xlsx"),
      rowNames = FALSE
    )
  }
}

cat("\nVersion B analysis completed.\n")
