# ============================================================================
# Reproducible vSV ethnicity association models
#
# Version A:
#   sequencing coverage included as a covariate (no coverage cutoff)
#
# Validation:
#   10-fold cross-validation for EVERY FDR-significant ethnicity-vSV pair
#   (association-level validation; repeated vSVs across ethnicities are retained)
#
# Version B:
#   original sequencing coverage >= 5X sensitivity analysis retained below
# ============================================================================

# ============================================================================
# Reproducible vSV ethnicity association models
# Version A: sequencing coverage included as a covariate (no coverage cutoff)
# ============================================================================

library(tidyverse)
library(broom)
library(openxlsx)
library(lmtest)
library(caret)

set.seed(12345)

# ----------------------------------------------------------------------------
# 1. Paths and data input
# ----------------------------------------------------------------------------
data_dir <- "."

cat(sprintf("[%s] Reading input files...\n", Sys.time()))
vsv_data <- read.csv(file.path(data_dir, "final_filtered_vsv.csv"), row.names = 1)
metadata <- read.csv("metadata.csv")
diet_data <- read.csv("final_filtered_diet.csv", check.names = FALSE)
medication_data <- read.csv("Medication.csv", check.names = FALSE)
urbanization_data <- read.csv("Urbanization.csv", check.names = FALSE)
general_metadata_data <- read.csv("General metadata.csv", check.names = FALSE)

# 2. Check required metadata columns
required_cols <- c("sample_id", "residency", "Age", "Sex", "BMI", "ethnicity")
missing_cols <- setdiff(required_cols, names(metadata))
if (length(missing_cols) > 0) stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))

# 3. Construct covariates
process_covariate <- function(data, covariate_name) {
  vars_df <- data %>% select(-sample_id)
  numeric_df <- vars_df %>% mutate(across(everything(), function(col) {
    if (is.numeric(col)) return(col)
    col_char <- trimws(as.character(col)); col_lower <- tolower(col_char)
    if (all(col_lower %in% c("yes", "no", "y", "n", "true", "false", "1", "0", "", NA))) {
      return(case_when(col_lower %in% c("yes", "y", "true", "1") ~ 1,
                       col_lower %in% c("no", "n", "false", "0") ~ 0,
                       TRUE ~ NA_real_))
    }
    num_col <- suppressWarnings(as.numeric(col_char))
    if (all(is.na(num_col[!is.na(col_char) & col_char != ""]))) num_col <- as.numeric(as.factor(col_char))
    num_col
  }))
  data %>% select(sample_id) %>% mutate(!!covariate_name := rowSums(numeric_df, na.rm = TRUE))
}

process_individual_covariates <- function(data, prefix = "general") {
  vars_df <- data %>% select(-sample_id) %>% mutate(across(everything(), function(col) {
    if (is.numeric(col)) return(col)
    col_char <- trimws(as.character(col)); col_lower <- tolower(col_char)
    if (all(col_lower %in% c("yes", "no", "y", "n", "true", "false", "1", "0", "", NA))) {
      return(case_when(col_lower %in% c("yes", "y", "true", "1") ~ 1,
                       col_lower %in% c("no", "n", "false", "0") ~ 0,
                       TRUE ~ NA_real_))
    }
    num_col <- suppressWarnings(as.numeric(col_char))
    if (sum(!is.na(num_col)) == sum(!is.na(col_char) & col_char != "")) return(num_col)
    factor(col_char)
  }))
  names(vars_df) <- paste0(prefix, "__", make.names(names(vars_df), unique = TRUE))
  bind_cols(data %>% select(sample_id), vars_df)
}

reverse_ordinal <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(x)
  min(x, na.rm = TRUE) + max(x, na.rm = TRUE) - x
}

animal_col <- names(urbanization_data)[tolower(trimws(names(urbanization_data))) == "contact with animals"]
if (length(animal_col) != 1) stop("Expected exactly one 'Contact with Animals' column in Urbanization.csv.")
urbanization_scoring_data <- urbanization_data
urbanization_scoring_data[[animal_col]] <- reverse_ordinal(urbanization_scoring_data[[animal_col]])

cat(sprintf("[%s] Constructing covariates...
", Sys.time()))
diet_covariate <- process_covariate(diet_data, "diet")
medication_covariate <- process_covariate(medication_data, "medication")
urbanization_covariate <- process_covariate(urbanization_scoring_data, "urbanization")
general_metadata_covariates <- process_individual_covariates(general_metadata_data, "general")
general_covariates <- setdiff(names(general_metadata_covariates), "sample_id")

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
  inner_join(diet_covariate, by = "sample_id") %>% inner_join(medication_covariate, by = "sample_id") %>%
  inner_join(urbanization_covariate, by = "sample_id") %>% inner_join(general_metadata_covariates, by = "sample_id") %>%
  left_join(abundance_matrix, by = "sample_id") %>% left_join(coverage_raw, by = "sample_id") %>%
  mutate(Sex = as.factor(Sex), residency = as.factor(residency), ethnicity = as.factor(ethnicity))

cat("\n================ Cohort summary ================\n")
cat("Sample counts by ethnicity:\n")
print(table(combined_metadata$ethnicity))
cat("==============================================\n\n")

# 6. vSV preprocessing
cat(sprintf("[%s] Preprocessing vSV profiles...\n", Sys.time()))
vsv_na <- vsv_data
vsv_na[vsv_na == 0] <- NA
vsv_scaled <- as.data.frame(scale(vsv_na, center = TRUE, scale = TRUE))

combined_data <- vsv_scaled %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(cols = -sample_id, names_to = "vsv", values_to = "expression") %>%
  inner_join(combined_metadata, by = "sample_id")

# ----------------------------------------------------------------------------
# 7. Screening: no coverage cutoff, coverage data required
# ----------------------------------------------------------------------------
screen_vsv_ethnicity_no_cov_filter <- function(data, target_eth) {
  data <- data %>% mutate(ethnicity_binary = ifelse(ethnicity == target_eth, 1, 0))
  results <- list()
  vsv_names <- unique(data$vsv)
  total_vsv <- length(vsv_names)
  
  cat(sprintf("  -> Spearman screening (no coverage cutoff), %d features\n", total_vsv))
  pb <- txtProgressBar(min = 0, max = total_vsv, style = 3)
  
  for (i in seq_along(vsv_names)) {
    setTxtProgressBar(pb, i)
    vsv_item <- vsv_names[i]
    vsv_species <- str_extract(vsv_item, "^[A-Za-z0-9.]+\\.[a-zA-Z0-9_.]+|^[A-Za-z0-9]+_[a-zA-Z0-9_]+")
    if (is.na(vsv_species)) vsv_species <- str_split(vsv_item, "_")[[1]][1]
    
    cov_col_name <- paste0(vsv_species, "_coverage")
    cov_vec <- if (!is.na(cov_col_name) && cov_col_name %in% colnames(data)) data[[cov_col_name]] else 0
    
    current_data <- data %>%
      filter(vsv == !!vsv_item) %>%
      mutate(species_coverage = cov_vec[data$vsv == !!vsv_item]) %>%
      drop_na(expression, ethnicity_binary, species_coverage)  # no cutoff; require non-missing data
    
    if (nrow(current_data) < 30) {
      results[[vsv_item]] <- tibble(vsv = vsv_item, p.value = NA, rho = NA, N_cov_available = nrow(current_data))
      next
    }
    
    cor_test <- cor.test(current_data$expression, current_data$ethnicity_binary, method = "spearman", exact = FALSE)
    results[[vsv_item]] <- tibble(
      vsv = vsv_item, p.value = cor_test$p.value, rho = cor_test$estimate, N_cov_available = nrow(current_data)
    )
  }
  close(pb)
  bind_rows(results)
}

# ----------------------------------------------------------------------------
# 8. Multivariable linear regression with coverage as a covariate
# ----------------------------------------------------------------------------
analyze_ethnicity_no_cov_filter <- function(target_eth, data) {
  data <- data %>% mutate(ethnicity_binary = ifelse(ethnicity == target_eth, 1, 0))
  results <- list()
  vsv_names <- unique(data$vsv)
  total_vsv <- length(vsv_names)
  
  cat(sprintf("\n  -> Multivariable linear regression (coverage as covariate, no cutoff), %d candidate features\n", total_vsv))
  pb <- txtProgressBar(min = 0, max = total_vsv, style = 3)
  
  for (i in seq_along(vsv_names)) {
    setTxtProgressBar(pb, i)
    vsv_item <- vsv_names[i]
    vsv_species <- str_extract(vsv_item, "^[A-Za-z0-9.]+\\.[a-zA-Z0-9_.]+|^[A-Za-z0-9]+_[a-zA-Z0-9_]+")
    if (is.na(vsv_species)) vsv_species <- str_split(vsv_item, "_")[[1]][1]
    
    cov_col_name <- paste0(vsv_species, "_coverage")
    abun_vec <- if (!is.na(vsv_species) && vsv_species %in% colnames(data)) data[[vsv_species]] else 0
    cov_vec  <- if (!is.na(cov_col_name) && cov_col_name %in% colnames(data)) data[[cov_col_name]] else abun_vec
    
    current_data <- data %>%
      filter(vsv == !!vsv_item) %>%
      mutate(species_abundance = abun_vec[data$vsv == !!vsv_item], 
             species_coverage = cov_vec[data$vsv == !!vsv_item]) %>%
      drop_na(any_of(c("expression", "ethnicity_binary", "Age", "Sex", "BMI", "residency",
                       "diet", "medication", "urbanization", "species_abundance", "species_coverage",
                       general_covariates)))  # no coverage cutoff
    
    if (nrow(current_data) < 10) next
    
    # include coverage as a covariate without applying a cutoff
    model_terms <- c("ethnicity_binary", "species_coverage", "species_abundance", "Age", "Sex", "BMI",
                     "residency", "diet", "medication", "urbanization", general_covariates)
    model_formula <- reformulate(model_terms, response = "expression")
    
    model <- tryCatch({ lm(model_formula, data = current_data) }, error = function(e) NULL)
    if (is.null(model)) next
    
    results[[vsv_item]] <- tidy(model, conf.int = TRUE) %>%
      filter(term == "ethnicity_binary") %>%
      mutate(vsv = vsv_item, species = ifelse(is.na(vsv_species), "Unknown", vsv_species), 
             N_samples = nrow(current_data)) %>%
      select(vsv, species, N_samples, estimate, std.error, statistic, p.value, conf.low, conf.high)
  }
  close(pb)
  
  if(length(results) == 0) return(NULL)
  bind_rows(results) %>% mutate(ethnicity = target_eth) %>% 
    select(ethnicity, vsv, species, N_samples, everything())
}

# ----------------------------------------------------------------------------
# 9. Run analyses for each ethnicity
# ----------------------------------------------------------------------------
if (!dir.exists("results_ethnicity")) dir.create("results_ethnicity")

ethnicities <- as.character(unique(combined_metadata$ethnicity))
all_ethnicity_results <- list()
all_screening_results <- list()

for (eth in ethnicities) {
  cat(sprintf("\n==================================================\n"))
  cat(sprintf("[%s] Analyzing %s ethnicity (1 vs others) - Version A (coverage as covariate)\n", Sys.time(), eth))
  
  screen_res <- screen_vsv_ethnicity_no_cov_filter(combined_data, eth) %>% mutate(target_ethnicity = eth)
  all_screening_results[[eth]] <- screen_res
  
  passed_vsv <- screen_res %>% filter(!is.na(p.value), p.value < 0.05) %>% pull(vsv)
  cat(sprintf("\n  -> Number of vSVs passing screening: %d\n", length(passed_vsv)))
  
  if(length(passed_vsv) == 0) next
  
  eth_data_filtered <- combined_data %>% filter(vsv %in% passed_vsv)
  current_result <- analyze_ethnicity_no_cov_filter(eth, eth_data_filtered)
  
  if (!is.null(current_result) && nrow(current_result) > 0) {
    current_result <- current_result %>% 
      mutate(p_adj = p.adjust(p.value, method = "BH"),
             fdr_note = ifelse(p_adj < 0.05, "FDR-significant", "Not significant after FDR correction"))
    all_ethnicity_results[[eth]] <- current_result
    write.xlsx(current_result, paste0("results_ethnicity/", eth, "_vs_Others_vsv_results_verA_coverage_as_covariate.xlsx"), rowNames = FALSE)
  }
}

final_results <- bind_rows(all_ethnicity_results)
significant_results <- final_results %>% filter(!is.na(p_adj), p_adj < 0.05) %>% arrange(p_adj)

write.xlsx(final_results, "results_ethnicity/all_ethnicities_vsv_results_verA_coverage_as_covariate.xlsx", rowNames = FALSE)
write.xlsx(significant_results, "results_ethnicity/all_ethnicities_significant_verA_coverage_as_covariate.xlsx", rowNames = FALSE)

# ============================================================================
# 10-fold cross-validation for significant ethnicity-vSV associations
# ============================================================================

# 10. Cross-validation and variable importance function
run_cv_with_varimp_ethnicity_no_cov_filter <- function(target_ethnicity, vsv_name, data_meta, vsv_scaled_data) {
  
  if (!vsv_name %in% colnames(vsv_scaled_data)) return(NULL)
  
  current_vsv_df <- vsv_scaled_data %>%
    select(all_of(vsv_name)) %>%
    rownames_to_column("sample_id") %>%
    rename(expression = !!sym(vsv_name))
  
  vsv_species <- str_extract(vsv_name, "^[A-Za-z0-9.]+\\.[a-zA-Z0-9_.]+|^[A-Za-z0-9]+_[a-zA-Z0-9_]+")
  if (is.na(vsv_species)) {
    vsv_species <- str_split(vsv_name, "_")[[1]][1]
  }
  cov_col_name <- paste0(vsv_species, "_coverage")
  
  df <- data_meta %>%
    inner_join(current_vsv_df, by = "sample_id") %>%
    mutate(
      ethnicity_binary = ifelse(ethnicity == target_ethnicity, 1, 0),
      species_abundance = if (!is.na(vsv_species) && vsv_species %in% colnames(.)) .[[vsv_species]] else 0,
      species_coverage = if (!is.na(cov_col_name) && cov_col_name %in% colnames(.)) .[[cov_col_name]] else species_abundance
    ) %>%
    drop_na(any_of(c("expression", "ethnicity_binary", "Age", "Sex", "BMI", "residency",
                     "diet", "medication", "urbanization", "species_abundance", "species_coverage",
                     general_covariates)))
  
  if (nrow(df) < 20) return(NULL)
  
  covars_all <- c("ethnicity_binary", "species_coverage", "species_abundance", "residency",
                  "urbanization", "diet", "medication", "Age", "Sex", "BMI", general_covariates)
  
  valid_covars <- covars_all[sapply(df[covars_all], function(x) length(unique(x)) > 1)]
  
  if (!"ethnicity_binary" %in% valid_covars) return(NULL)
  
  full_formula <- as.formula(paste("expression ~", paste(valid_covars, collapse = " + ")))
  fit_full <- tryCatch({ lm(full_formula, data = df) }, error = function(e) NULL)
  
  if (is.null(fit_full)) return(NULL)
  r2_full <- summary(fit_full)$r.squared
  
  get_marginal_delta_r2 <- function(drop_var) {
    if (!drop_var %in% valid_covars) return(0)
    sub_covars <- setdiff(valid_covars, drop_var)
    if (length(sub_covars) == 0) return(r2_full)
    f_drop <- as.formula(paste("expression ~", paste(sub_covars, collapse = " + ")))
    fit_drop <- lm(f_drop, data = df)
    return(max(0, r2_full - summary(fit_drop)$r.squared))
  }
  
  marginal_r2_ethnicity   <- get_marginal_delta_r2("ethnicity_binary")
  marginal_r2_coverage    <- get_marginal_delta_r2("species_coverage")
  marginal_r2_residency   <- get_marginal_delta_r2("residency")
  marginal_r2_urban       <- get_marginal_delta_r2("urbanization")
  marginal_r2_diet        <- get_marginal_delta_r2("diet")
  marginal_r2_spec        <- get_marginal_delta_r2("species_abundance")
  
  ctrl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)
  
  model_cv <- tryCatch({
    suppressWarnings(train(full_formula, data = df, method = "lm", trControl = ctrl, metric = "Rsquared"))
  }, error = function(e) NULL)
  
  if (is.null(model_cv)) return(NULL)
  
  resample_stats <- model_cv$resample %>%
    summarise(
      RMSE_SD = sd(RMSE, na.rm = TRUE),
      Rsquared_SD = sd(Rsquared, na.rm = TRUE)
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
    vsv = vsv_name,
    species = ifelse(is.na(vsv_species), "Unknown", vsv_species),
    N_no_filter = nrow(df),
    Full_Model_R2 = round(r2_full, 4),
    CV_Mean_Rsquared = round(mean(model_cv$results$Rsquared, na.rm = TRUE), 4),
    CV_Mean_RMSE = round(mean(model_cv$results$RMSE, na.rm = TRUE), 4),
    CV_Rsquared_SD = round(resample_stats$Rsquared_SD, 4),
    CV_RMSE_SD = round(resample_stats$RMSE_SD, 4),
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


# ----------------------------------------------------------------------------
# 11. Extract significant ethnicity-vSV PAIRS (association-level, not unique vSV)
# ----------------------------------------------------------------------------

sig_ethnicity_pairs <- significant_results %>%
  filter(!is.na(p_adj), p_adj < 0.05) %>%
  transmute(
    ethnicity = as.character(ethnicity),
    vsv = as.character(vsv)
  ) %>%
  distinct()

cat("\n============================================================\n")
cat("Version A: 10-fold CV for significant ethnicity-vSV associations\n")
cat("============================================================\n")
cat("FDR-significant result rows: ", nrow(significant_results), "\n", sep = "")
cat("Unique ethnicity-vSV pairs: ", nrow(sig_ethnicity_pairs), "\n", sep = "")
cat("Unique vSVs represented: ", n_distinct(sig_ethnicity_pairs$vsv), "\n", sep = "")

# ----------------------------------------------------------------------------
# 12. Run 10-fold CV for every significant ethnicity-vSV association
# ----------------------------------------------------------------------------

cv_eth_list <- list()
cv_skipped_list <- list()

if (nrow(sig_ethnicity_pairs) > 0) {
  
  for (i in seq_len(nrow(sig_ethnicity_pairs))) {
    
    eth <- sig_ethnicity_pairs$ethnicity[i]
    vsv_name <- sig_ethnicity_pairs$vsv[i]
    
    cat(
      sprintf(
        "[%d/%d] 10-fold CV: %s | %s\n",
        i,
        nrow(sig_ethnicity_pairs),
        eth,
        vsv_name
      )
    )
    
    if (!vsv_name %in% colnames(vsv_scaled)) {
      cv_skipped_list[[i]] <- data.frame(
        Target_Ethnicity = eth,
        vsv = vsv_name,
        Reason = "vSV not found in Version A scaled vSV matrix",
        stringsAsFactors = FALSE
      )
      next
    }
    
    res_cv <- tryCatch(
      run_cv_with_varimp_ethnicity_no_cov_filter(
        eth,
        vsv_name,
        combined_metadata,
        vsv_scaled
      ),
      error = function(e) {
        cv_skipped_list[[i]] <<- data.frame(
          Target_Ethnicity = eth,
          vsv = vsv_name,
          Reason = paste0("CV error: ", conditionMessage(e)),
          stringsAsFactors = FALSE
        )
        NULL
      }
    )
    
    if (!is.null(res_cv)) {
      cv_eth_list[[i]] <- res_cv
    } else if (is.null(cv_skipped_list[[i]])) {
      cv_skipped_list[[i]] <- data.frame(
        Target_Ethnicity = eth,
        vsv = vsv_name,
        Reason = paste(
          "Returned NULL:",
          "N < 20, ethnicity_binary has one level,",
          "linear model failed, or caret CV failed"
        ),
        stringsAsFactors = FALSE
      )
    }
  }
}

cv_eth_df <- bind_rows(cv_eth_list)
cv_skipped_df <- bind_rows(cv_skipped_list)

if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

cv_output_file <- "results/cv_10fold_ethnicity_vsv_significant_verA.xlsx"

wb_cv <- createWorkbook()

addWorksheet(wb_cv, "CV_results")
writeData(wb_cv, "CV_results", cv_eth_df)
if (ncol(cv_eth_df) > 0) {
  setColWidths(
    wb_cv,
    "CV_results",
    cols = seq_len(ncol(cv_eth_df)),
    widths = "auto"
  )
}

addWorksheet(wb_cv, "Skipped_pairs")
writeData(wb_cv, "Skipped_pairs", cv_skipped_df)
if (ncol(cv_skipped_df) > 0) {
  setColWidths(
    wb_cv,
    "Skipped_pairs",
    cols = seq_len(ncol(cv_skipped_df)),
    widths = "auto"
  )
}

cv_summary <- data.frame(
  Metric = c(
    "FDR-significant result rows",
    "Unique ethnicity-vSV pairs requested",
    "Unique vSVs represented",
    "CV completed",
    "CV skipped"
  ),
  N = c(
    nrow(significant_results),
    nrow(sig_ethnicity_pairs),
    n_distinct(sig_ethnicity_pairs$vsv),
    nrow(cv_eth_df),
    nrow(cv_skipped_df)
  )
)

addWorksheet(wb_cv, "Run_summary")
writeData(wb_cv, "Run_summary", cv_summary)
setColWidths(wb_cv, "Run_summary", cols = 1:2, widths = "auto")

saveWorkbook(
  wb_cv,
  cv_output_file,
  overwrite = TRUE
)

cat("\nVersion A 10-fold CV completed.\n")
cat("CV completed: ", nrow(cv_eth_df), "\n", sep = "")
cat("CV skipped: ", nrow(cv_skipped_df), "\n", sep = "")
cat("Saved: ", cv_output_file, "\n", sep = "")


# ============================================================================
# vSV multi-ethnicity analysis - Version B: sequencing coverage >= 5X
# Sensitivity analysis with a strict sequencing coverage cutoff
# ============================================================================

library(tidyverse)
library(broom)
library(openxlsx)
library(lmtest)
library(caret)

set.seed(12345)

# ----------------------------------------------------------------------------
# 1. Paths and data input
# ----------------------------------------------------------------------------
data_dir <- "."

cat(sprintf("[%s] Reading input files...\n", Sys.time()))
vsv_data <- read.csv(file.path(data_dir, "final_filtered_vsv.csv"), row.names = 1)
metadata <- read.csv("metadata.csv")
diet_data <- read.csv("final_filtered_diet.csv", check.names = FALSE)
medication_data <- read.csv("Medication.csv", check.names = FALSE)
urbanization_data <- read.csv("Urbanization.csv", check.names = FALSE)
general_metadata_data <- read.csv("General metadata.csv", check.names = FALSE)

# 2. Check required metadata columns
required_cols <- c("sample_id", "residency", "Age", "Sex", "BMI", "ethnicity")
missing_cols <- setdiff(required_cols, names(metadata))
if (length(missing_cols) > 0) stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))

# 3. Construct covariates
process_covariate <- function(data, covariate_name) {
  vars_df <- data %>% select(-sample_id)
  numeric_df <- vars_df %>% mutate(across(everything(), function(col) {
    if (is.numeric(col)) return(col)
    col_char <- trimws(as.character(col)); col_lower <- tolower(col_char)
    if (all(col_lower %in% c("yes", "no", "y", "n", "true", "false", "1", "0", "", NA))) {
      return(case_when(col_lower %in% c("yes", "y", "true", "1") ~ 1,
                       col_lower %in% c("no", "n", "false", "0") ~ 0,
                       TRUE ~ NA_real_))
    }
    num_col <- suppressWarnings(as.numeric(col_char))
    if (all(is.na(num_col[!is.na(col_char) & col_char != ""]))) num_col <- as.numeric(as.factor(col_char))
    num_col
  }))
  data %>% select(sample_id) %>% mutate(!!covariate_name := rowSums(numeric_df, na.rm = TRUE))
}

process_individual_covariates <- function(data, prefix = "general") {
  vars_df <- data %>% select(-sample_id) %>% mutate(across(everything(), function(col) {
    if (is.numeric(col)) return(col)
    col_char <- trimws(as.character(col)); col_lower <- tolower(col_char)
    if (all(col_lower %in% c("yes", "no", "y", "n", "true", "false", "1", "0", "", NA))) {
      return(case_when(col_lower %in% c("yes", "y", "true", "1") ~ 1,
                       col_lower %in% c("no", "n", "false", "0") ~ 0,
                       TRUE ~ NA_real_))
    }
    num_col <- suppressWarnings(as.numeric(col_char))
    if (sum(!is.na(num_col)) == sum(!is.na(col_char) & col_char != "")) return(num_col)
    factor(col_char)
  }))
  names(vars_df) <- paste0(prefix, "__", make.names(names(vars_df), unique = TRUE))
  bind_cols(data %>% select(sample_id), vars_df)
}

reverse_ordinal <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(x)
  min(x, na.rm = TRUE) + max(x, na.rm = TRUE) - x
}

animal_col <- names(urbanization_data)[tolower(trimws(names(urbanization_data))) == "contact with animals"]
if (length(animal_col) != 1) stop("Expected exactly one 'Contact with Animals' column in Urbanization.csv.")
urbanization_scoring_data <- urbanization_data
urbanization_scoring_data[[animal_col]] <- reverse_ordinal(urbanization_scoring_data[[animal_col]])

cat(sprintf("[%s] Constructing covariates...
", Sys.time()))
diet_covariate <- process_covariate(diet_data, "diet")
medication_covariate <- process_covariate(medication_data, "medication")
urbanization_covariate <- process_covariate(urbanization_scoring_data, "urbanization")
general_metadata_covariates <- process_individual_covariates(general_metadata_data, "general")
general_covariates <- setdiff(names(general_metadata_covariates), "sample_id")

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
  inner_join(diet_covariate, by = "sample_id") %>% inner_join(medication_covariate, by = "sample_id") %>%
  inner_join(urbanization_covariate, by = "sample_id") %>% inner_join(general_metadata_covariates, by = "sample_id") %>%
  left_join(abundance_matrix, by = "sample_id") %>% left_join(coverage_raw, by = "sample_id") %>%
  mutate(Sex = as.factor(Sex), residency = as.factor(residency), ethnicity = as.factor(ethnicity))

cat("\n================ Cohort summary ================\n")
cat("Sample counts by ethnicity:\n")
print(table(combined_metadata$ethnicity))
cat("==============================================\n\n")

# 6. vSV preprocessing
cat(sprintf("[%s] Preprocessing vSV profiles...\n", Sys.time()))
vsv_prevalence <- colMeans(vsv_data != 0, na.rm = TRUE)
vsv_filtered <- vsv_data[, names(vsv_prevalence[vsv_prevalence >= 0.1]), drop = FALSE]
vsv_na <- vsv_filtered
vsv_na[vsv_na == 0] <- NA
vsv_scaled <- as.data.frame(scale(vsv_na, center = TRUE, scale = TRUE))

combined_data <- vsv_scaled %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(cols = -sample_id, names_to = "vsv", values_to = "expression") %>%
  inner_join(combined_metadata, by = "sample_id")

# ----------------------------------------------------------------------------
# 7. Screening: coverage >= 5X, coverage data required
# ----------------------------------------------------------------------------
screen_vsv_ethnicity_cov5x <- function(data, target_eth) {
  data <- data %>% mutate(ethnicity_binary = ifelse(ethnicity == target_eth, 1, 0))
  results <- list()
  vsv_names <- unique(data$vsv)
  total_vsv <- length(vsv_names)
  
  cat(sprintf("  -> Spearman screening (coverage >= 5X), %d features\n", total_vsv))
  pb <- txtProgressBar(min = 0, max = total_vsv, style = 3)
  
  for (i in seq_along(vsv_names)) {
    setTxtProgressBar(pb, i)
    vsv_item <- vsv_names[i]
    vsv_species <- str_extract(vsv_item, "^[A-Za-z0-9.]+\\.[a-zA-Z0-9_.]+|^[A-Za-z0-9]+_[a-zA-Z0-9_]+")
    if (is.na(vsv_species)) vsv_species <- str_split(vsv_item, "_")[[1]][1]
    
    cov_col_name <- paste0(vsv_species, "_coverage")
    cov_vec <- if (!is.na(cov_col_name) && cov_col_name %in% colnames(data)) data[[cov_col_name]] else 0
    
    current_data <- data %>%
      filter(vsv == !!vsv_item) %>%
      mutate(species_coverage = cov_vec[data$vsv == !!vsv_item]) %>%
      drop_na(expression, ethnicity_binary, species_coverage) %>%
      filter(species_coverage >= 5)
    
    if (nrow(current_data) < 30) {
      results[[vsv_item]] <- tibble(vsv = vsv_item, p.value = NA, rho = NA, N_cov5x = nrow(current_data))
      next
    }
    
    cor_test <- cor.test(current_data$expression, current_data$ethnicity_binary, method = "spearman", exact = FALSE)
    results[[vsv_item]] <- tibble(
      vsv = vsv_item, p.value = cor_test$p.value, rho = cor_test$estimate, N_cov5x = nrow(current_data)
    )
  }
  close(pb)
  bind_rows(results)
}

# ----------------------------------------------------------------------------
# 8. Multivariable linear regression with coverage as a covariate
# ----------------------------------------------------------------------------
analyze_ethnicity_cov5x <- function(target_eth, data) {
  data <- data %>% mutate(ethnicity_binary = ifelse(ethnicity == target_eth, 1, 0))
  results <- list()
  vsv_names <- unique(data$vsv)
  total_vsv <- length(vsv_names)
  
  cat(sprintf("\n  -> Multivariable linear regression (coverage as covariate, coverage >= 5X), %d candidate features\n", total_vsv))
  pb <- txtProgressBar(min = 0, max = total_vsv, style = 3)
  
  for (i in seq_along(vsv_names)) {
    setTxtProgressBar(pb, i)
    vsv_item <- vsv_names[i]
    vsv_species <- str_extract(vsv_item, "^[A-Za-z0-9.]+\\.[a-zA-Z0-9_.]+|^[A-Za-z0-9]+_[a-zA-Z0-9_]+")
    if (is.na(vsv_species)) vsv_species <- str_split(vsv_item, "_")[[1]][1]
    
    cov_col_name <- paste0(vsv_species, "_coverage")
    abun_vec <- if (!is.na(vsv_species) && vsv_species %in% colnames(data)) data[[vsv_species]] else 0
    cov_vec  <- if (!is.na(cov_col_name) && cov_col_name %in% colnames(data)) data[[cov_col_name]] else abun_vec
    
    current_data <- data %>%
      filter(vsv == !!vsv_item) %>%
      mutate(species_abundance = abun_vec[data$vsv == !!vsv_item], 
             species_coverage = cov_vec[data$vsv == !!vsv_item]) %>%
      drop_na(any_of(c("expression", "ethnicity_binary", "Age", "Sex", "BMI", "residency",
                       "diet", "medication", "urbanization", "species_abundance", "species_coverage",
                       general_covariates))) %>%
      filter(species_coverage >= 5)
    
    if (nrow(current_data) < 10) next
    
    # include coverage as a covariate without applying a cutoff
    model_terms <- c("ethnicity_binary", "species_coverage", "species_abundance", "Age", "Sex", "BMI",
                     "residency", "diet", "medication", "urbanization", general_covariates)
    model_formula <- reformulate(model_terms, response = "expression")
    
    model <- tryCatch({ lm(model_formula, data = current_data) }, error = function(e) NULL)
    if (is.null(model)) next
    
    results[[vsv_item]] <- tidy(model, conf.int = TRUE) %>%
      filter(term == "ethnicity_binary") %>%
      mutate(vsv = vsv_item, species = ifelse(is.na(vsv_species), "Unknown", vsv_species), 
             N_samples = nrow(current_data)) %>%
      select(vsv, species, N_samples, estimate, std.error, statistic, p.value, conf.low, conf.high)
  }
  close(pb)
  
  if(length(results) == 0) return(NULL)
  bind_rows(results) %>% mutate(ethnicity = target_eth) %>% 
    select(ethnicity, vsv, species, N_samples, everything())
}

# ----------------------------------------------------------------------------
# 9. Run analyses for each ethnicity
# ----------------------------------------------------------------------------
if (!dir.exists("results_ethnicity")) dir.create("results_ethnicity")

ethnicities <- as.character(unique(combined_metadata$ethnicity))
all_ethnicity_results <- list()
all_screening_results <- list()

for (eth in ethnicities) {
  cat(sprintf("\n==================================================\n"))
  cat(sprintf("[%s] Analyzing %s ethnicity (1 vs others) - Version B (coverage >= 5X)\n", Sys.time(), eth))
  
  screen_res <- screen_vsv_ethnicity_cov5x(combined_data, eth) %>% mutate(target_ethnicity = eth)
  all_screening_results[[eth]] <- screen_res
  
  passed_vsv <- screen_res %>% filter(!is.na(p.value), p.value < 0.05) %>% pull(vsv)
  cat(sprintf("\n  -> Number of vSVs passing screening: %d\n", length(passed_vsv)))
  
  if(length(passed_vsv) == 0) next
  
  eth_data_filtered <- combined_data %>% filter(vsv %in% passed_vsv)
  current_result <- analyze_ethnicity_cov5x(eth, eth_data_filtered)
  
  if (!is.null(current_result) && nrow(current_result) > 0) {
    current_result <- current_result %>% 
      mutate(p_adj = p.adjust(p.value, method = "BH"),
             fdr_note = ifelse(p_adj < 0.05, "FDR-significant", "Not significant after FDR correction"))
    all_ethnicity_results[[eth]] <- current_result
    write.xlsx(current_result, paste0("results_ethnicity/", eth, "_vs_Others_vsv_results_verB_cov5x.xlsx"), rowNames = FALSE)
  }
}

final_results <- bind_rows(all_ethnicity_results)
significant_results <- final_results %>% filter(!is.na(p_adj), p_adj < 0.05) %>% arrange(p_adj)

write.xlsx(final_results, "results_ethnicity/all_ethnicities_vsv_results_verB_cov5x.xlsx", rowNames = FALSE)
write.xlsx(significant_results, "results_ethnicity/all_ethnicities_significant_verB_cov5x.xlsx", rowNames = FALSE)

# ============================================================================
# 10-fold cross-validation for significant ethnicity-vSV associations
# ============================================================================

# 10. Cross-validation and variable importance function
run_cv_with_varimp_ethnicity_cov5x <- function(target_ethnicity, vsv_name, data_meta, vsv_scaled_data) {
  
  if (!vsv_name %in% colnames(vsv_scaled_data)) return(NULL)
  
  current_vsv_df <- vsv_scaled_data %>%
    select(all_of(vsv_name)) %>%
    rownames_to_column("sample_id") %>%
    rename(expression = !!sym(vsv_name))
  
  vsv_species <- str_extract(vsv_name, "^[A-Za-z0-9.]+\\.[a-zA-Z0-9_.]+|^[A-Za-z0-9]+_[a-zA-Z0-9_]+")
  if (is.na(vsv_species)) {
    vsv_species <- str_split(vsv_name, "_")[[1]][1]
  }
  cov_col_name <- paste0(vsv_species, "_coverage")
  
  df <- data_meta %>%
    inner_join(current_vsv_df, by = "sample_id") %>%
    mutate(
      ethnicity_binary = ifelse(ethnicity == target_ethnicity, 1, 0),
      species_abundance = if (!is.na(vsv_species) && vsv_species %in% colnames(.)) .[[vsv_species]] else 0,
      species_coverage = if (!is.na(cov_col_name) && cov_col_name %in% colnames(.)) .[[cov_col_name]] else species_abundance
    ) %>%
    drop_na(any_of(c("expression", "ethnicity_binary", "Age", "Sex", "BMI", "residency",
                     "diet", "medication", "urbanization", "species_abundance", "species_coverage",
                     general_covariates))) %>%
    filter(species_coverage >= 5)
  
  if (nrow(df) < 20) return(NULL)
  
  covars_all <- c("ethnicity_binary", "species_coverage", "species_abundance", "residency",
                  "urbanization", "diet", "medication", "Age", "Sex", "BMI", general_covariates)
  
  valid_covars <- covars_all[sapply(df[covars_all], function(x) length(unique(x)) > 1)]
  
  if (!"ethnicity_binary" %in% valid_covars) return(NULL)
  
  full_formula <- as.formula(paste("expression ~", paste(valid_covars, collapse = " + ")))
  fit_full <- tryCatch({ lm(full_formula, data = df) }, error = function(e) NULL)
  
  if (is.null(fit_full)) return(NULL)
  r2_full <- summary(fit_full)$r.squared
  
  get_marginal_delta_r2 <- function(drop_var) {
    if (!drop_var %in% valid_covars) return(0)
    sub_covars <- setdiff(valid_covars, drop_var)
    if (length(sub_covars) == 0) return(r2_full)
    f_drop <- as.formula(paste("expression ~", paste(sub_covars, collapse = " + ")))
    fit_drop <- lm(f_drop, data = df)
    return(max(0, r2_full - summary(fit_drop)$r.squared))
  }
  
  marginal_r2_ethnicity   <- get_marginal_delta_r2("ethnicity_binary")
  marginal_r2_coverage    <- get_marginal_delta_r2("species_coverage")
  marginal_r2_residency   <- get_marginal_delta_r2("residency")
  marginal_r2_urban       <- get_marginal_delta_r2("urbanization")
  marginal_r2_diet        <- get_marginal_delta_r2("diet")
  marginal_r2_spec        <- get_marginal_delta_r2("species_abundance")
  
  ctrl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)
  
  model_cv <- tryCatch({
    suppressWarnings(train(full_formula, data = df, method = "lm", trControl = ctrl, metric = "Rsquared"))
  }, error = function(e) NULL)
  
  if (is.null(model_cv)) return(NULL)
  
  resample_stats <- model_cv$resample %>%
    summarise(
      RMSE_SD = sd(RMSE, na.rm = TRUE),
      Rsquared_SD = sd(Rsquared, na.rm = TRUE)
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
    vsv = vsv_name,
    species = ifelse(is.na(vsv_species), "Unknown", vsv_species),
    N_cov5x = nrow(df),
    Full_Model_R2 = round(r2_full, 4),
    CV_Mean_Rsquared = round(mean(model_cv$results$Rsquared, na.rm = TRUE), 4),
    CV_Mean_RMSE = round(mean(model_cv$results$RMSE, na.rm = TRUE), 4),
    CV_Rsquared_SD = round(resample_stats$Rsquared_SD, 4),
    CV_RMSE_SD = round(resample_stats$RMSE_SD, 4),
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

# 11. Extract significant ethnicity-vSV pairs
sig_ethnicity_pairs <- if (exists("significant_results")) {
  significant_results %>% 
    filter(!is.na(p_adj), p_adj < 0.05) %>% 
    select(ethnicity, vsv)
} else if (exists("final_results")) {
  final_results %>% 
    filter(!is.na(p_adj), p_adj < 0.05) %>% 
    select(ethnicity, vsv)
} else {
  tibble(ethnicity = character(), vsv = character())
}

# 12. Run cross-validation
if (nrow(sig_ethnicity_pairs) > 0) {
  cat(sprintf("\nRunning 10-fold CV for %d significant ethnicity-vSV pairs (coverage >= 5X)...\n", nrow(sig_ethnicity_pairs)))
  cv_eth_list <- list()
  
  for (i in 1:nrow(sig_ethnicity_pairs)) {
    eth <- sig_ethnicity_pairs$ethnicity[i]
    vsv <- sig_ethnicity_pairs$vsv[i]
    
    res_cv <- run_cv_with_varimp_ethnicity_cov5x(eth, vsv, combined_metadata, vsv_scaled)
    if (!is.null(res_cv)) cv_eth_list[[i]] <- res_cv
  }
  
  cv_eth_df <- bind_rows(cv_eth_list)
  
  # Create output directory if needed
  if (!dir.exists("results")) {
    dir.create("results")
    cat("Created results directory\n")
  }
  
  # Save cross-validation results
  openxlsx::write.xlsx(cv_eth_df, "results/cv_10fold_ethnicity_vsv_sig_svs_cov5x.xlsx", rowNames = FALSE)
  cat("Saved: results/cv_10fold_ethnicity_vsv_sig_svs_cov5x.xlsx\n")
}
