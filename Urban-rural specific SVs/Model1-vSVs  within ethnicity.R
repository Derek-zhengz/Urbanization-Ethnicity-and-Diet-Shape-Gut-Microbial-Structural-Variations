# ============================================================================
# Ethnicity-stratified Rural vs Urban vSV analysis
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
vsv_filtered <- vsv_scaled

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
# Ethnicity-stratified Rural vs Urban vSV analysis
# Coverage is included as a covariate; no coverage cutoff is applied.
# ============================================================================

cat("\n================ Ethnicity-stratified Rural vs Urban analysis ================\n")

screen_vsv_within_ethnicity <- function(data, target_ethnicity) {
  ethnicity_data <- data %>%
    filter(ethnicity == target_ethnicity) %>%
    mutate(residency_binary = ifelse(residency == "Urban", 1, 0))
  
  results <- list()
  vsv_names <- unique(ethnicity_data$vsv)
  
  cat(sprintf("[%s] Screening ethnicity %s: %d vSVs\n",
              Sys.time(), target_ethnicity, length(vsv_names)))
  
  pb <- txtProgressBar(min = 0, max = length(vsv_names), style = 3)
  
  for (i in seq_along(vsv_names)) {
    setTxtProgressBar(pb, i)
    vsv_item <- vsv_names[i]
    
    vsv_species <- extract_species(vsv_item)
    cov_col_name <- paste0(vsv_species, "_coverage")
    
    current_data <- ethnicity_data %>%
      filter(vsv == !!vsv_item)
    
    if (!is.na(cov_col_name) && cov_col_name %in% colnames(current_data)) {
      current_data <- current_data %>%
        mutate(species_coverage = .data[[cov_col_name]])
    } else {
      current_data <- current_data %>%
        mutate(species_coverage = NA_real_)
    }
    
    current_data <- current_data %>%
      drop_na(expression, residency_binary, species_coverage)
    
    if (nrow(current_data) < 30 ||
        length(unique(current_data$residency_binary)) < 2) {
      results[[vsv_item]] <- tibble(
        vsv = vsv_item,
        p.value = NA_real_,
        rho = NA_real_,
        N = nrow(current_data)
      )
      next
    }
    
    cor_test <- tryCatch(
      cor.test(
        current_data$expression,
        current_data$residency_binary,
        method = "spearman",
        exact = FALSE
      ),
      error = function(e) NULL
    )
    
    if (is.null(cor_test)) next
    
    results[[vsv_item]] <- tibble(
      vsv = vsv_item,
      p.value = cor_test$p.value,
      rho = unname(cor_test$estimate),
      N = nrow(current_data)
    )
  }
  
  close(pb)
  bind_rows(results)
}


analyze_vsv_within_ethnicity <- function(data, target_ethnicity, selected_vsv) {
  ethnicity_data <- data %>%
    filter(ethnicity == target_ethnicity) %>%
    mutate(residency_binary = ifelse(residency == "Urban", 1, 0))
  
  results <- list()
  
  cat(sprintf("[%s] Adjusted models for ethnicity %s: %d candidate vSVs\n",
              Sys.time(), target_ethnicity, length(selected_vsv)))
  
  pb <- txtProgressBar(min = 0, max = length(selected_vsv), style = 3)
  
  for (i in seq_along(selected_vsv)) {
    setTxtProgressBar(pb, i)
    vsv_item <- selected_vsv[i]
    
    vsv_species <- extract_species(vsv_item)
    cov_col_name <- paste0(vsv_species, "_coverage")
    
    current_data <- ethnicity_data %>%
      filter(vsv == !!vsv_item)
    
    if (!is.na(vsv_species) && vsv_species %in% colnames(current_data)) {
      current_data <- current_data %>%
        mutate(species_abundance = .data[[vsv_species]])
    } else {
      current_data <- current_data %>%
        mutate(species_abundance = NA_real_)
    }
    
    if (!is.na(cov_col_name) && cov_col_name %in% colnames(current_data)) {
      current_data <- current_data %>%
        mutate(species_coverage = .data[[cov_col_name]])
    } else {
      current_data <- current_data %>%
        mutate(species_coverage = NA_real_)
    }
    
    current_data <- current_data %>%
      drop_na(
        expression, residency_binary, Age, Sex, BMI,
        diet, medication, urbanization,
        all_of(general_covars),
        species_abundance, species_coverage
      ) %>%
      droplevels()
    
    if (nrow(current_data) < 10 ||
        length(unique(current_data$residency_binary)) < 2 ||
        length(unique(current_data$expression)) < 2) {
      next
    }
    
    model_vars <- c(
      "residency_binary", "species_coverage", "species_abundance",
      "Age", "Sex", "BMI", "diet", "medication", "urbanization",
      general_covars
    )
    
    valid_vars <- model_vars[sapply(model_vars, function(v) {
      if (!v %in% colnames(current_data)) return(FALSE)
      x <- current_data[[v]]
      if (is.factor(x) || is.character(x)) {
        length(unique(na.omit(as.character(x)))) >= 2
      } else {
        length(unique(na.omit(x))) >= 2
      }
    })]
    
    if (!"residency_binary" %in% valid_vars) next
    
    model <- tryCatch(
      lm(make_formula("expression", valid_vars), data = current_data),
      error = function(e) NULL
    )
    
    if (is.null(model)) next
    
    model_result <- tidy(model, conf.int = TRUE) %>%
      filter(term == "residency_binary")
    
    if (nrow(model_result) == 0) next
    
    results[[vsv_item]] <- model_result %>%
      mutate(
        ethnicity = target_ethnicity,
        comparison = "Rural_vs_Urban",
        vsv = vsv_item,
        species = ifelse(is.na(vsv_species), "Unknown", vsv_species),
        N_samples = nrow(current_data)
      ) %>%
      select(
        ethnicity, comparison, vsv, species, N_samples,
        estimate, std.error, statistic, p.value, conf.low, conf.high
      )
  }
  
  close(pb)
  
  if (length(results) == 0) return(NULL)
  bind_rows(results)
}


# ----------------------------------------------------------------------------
# Run Rural vs Urban analysis separately within each ethnicity
# ----------------------------------------------------------------------------

ethnicities <- levels(combined_metadata$ethnicity)
all_results <- list()
screening_results <- list()

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
  
  screen_res <- screen_vsv_within_ethnicity(combined_data, eth) %>%
    mutate(ethnicity = eth)
  
  screening_results[[eth]] <- screen_res
  
  selected_vsv <- screen_res %>%
    filter(!is.na(p.value), p.value < 0.05) %>%
    pull(vsv)
  
  cat(sprintf("vSVs passing raw screening (P < 0.05): %d\n", length(selected_vsv)))
  
  if (length(selected_vsv) == 0) next
  
  current_result <- analyze_vsv_within_ethnicity(
    combined_data,
    eth,
    selected_vsv
  )
  
  if (!is.null(current_result) && nrow(current_result) > 0) {
    current_result <- current_result %>%
      mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
      arrange(p_adj)
    
    all_results[[eth]] <- current_result
    
    write.xlsx(
      current_result,
      file.path(
        results_dir,
        paste0(make.names(eth), "_Rural_vs_Urban_vsv_results.xlsx")
      ),
      rowNames = FALSE
    )
  }
}

screening_results_df <- bind_rows(screening_results)
final_results <- bind_rows(all_results)

significant_results <- final_results %>%
  filter(!is.na(p_adj), p_adj < 0.05) %>%
  arrange(ethnicity, p_adj)

write.xlsx(
  screening_results_df,
  file.path(results_dir, "all_ethnicities_Rural_vs_Urban_vsv_screening.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  final_results,
  file.path(results_dir, "all_ethnicities_Rural_vs_Urban_vsv_results.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  significant_results,
  file.path(results_dir, "all_ethnicities_significant_Rural_vs_Urban_vsv.xlsx"),
  rowNames = FALSE
)


# ----------------------------------------------------------------------------
# Ten-fold cross-validation and marginal R2 for significant associations
# ----------------------------------------------------------------------------

run_cv_within_ethnicity <- function(target_ethnicity, vsv_name, data_meta, vsv_df) {
  if (!vsv_name %in% colnames(vsv_df)) return(NULL)
  
  vsv_expr_df <- vsv_df %>%
    select(sample_id, expression = all_of(vsv_name))
  
  vsv_species <- extract_species(vsv_name)
  cov_col_name <- paste0(vsv_species, "_coverage")
  
  df <- data_meta %>%
    filter(ethnicity == target_ethnicity) %>%
    inner_join(vsv_expr_df, by = "sample_id") %>%
    mutate(
      expression = as.numeric(expression),
      residency_binary = ifelse(residency == "Urban", 1, 0)
    )
  
  if (!is.na(vsv_species) && vsv_species %in% colnames(df)) {
    df$species_abundance <- df[[vsv_species]]
  } else {
    df$species_abundance <- NA_real_
  }
  
  if (!is.na(cov_col_name) && cov_col_name %in% colnames(df)) {
    df$species_coverage <- df[[cov_col_name]]
  } else {
    df$species_coverage <- NA_real_
  }
  
  df <- df %>%
    drop_na(
      expression, residency_binary, Age, Sex, BMI,
      diet, medication, urbanization,
      all_of(general_covars),
      species_abundance, species_coverage
    ) %>%
    droplevels()
  
  if (nrow(df) < 20 ||
      length(unique(df$expression)) < 2 ||
      length(unique(df$residency_binary)) < 2) {
    return(NULL)
  }
  
  covars <- c(
    "residency_binary", "Age", "Sex", "BMI",
    "diet", "medication", "urbanization",
    general_covars, "species_abundance", "species_coverage"
  )
  
  valid_covars <- covars[sapply(covars, function(v) {
    if (!v %in% colnames(df)) return(FALSE)
    x <- df[[v]]
    if (is.factor(x) || is.character(x)) {
      length(unique(na.omit(as.character(x)))) >= 2
    } else {
      length(unique(na.omit(x))) >= 2
    }
  })]
  
  if (!"residency_binary" %in% valid_covars) return(NULL)
  
  full_lm <- tryCatch(
    lm(make_formula("expression", valid_covars), data = df),
    error = function(e) NULL
  )
  
  if (is.null(full_lm)) return(NULL)
  
  full_r2 <- summary(full_lm)$r.squared
  
  get_delta_r2 <- function(drop_var) {
    if (!drop_var %in% valid_covars) return(0)
    
    keep_vars <- setdiff(valid_covars, drop_var)
    if (length(keep_vars) == 0) return(full_r2)
    
    fit_red <- tryCatch(
      lm(make_formula("expression", keep_vars), data = df),
      error = function(e) NULL
    )
    
    if (is.null(fit_red)) return(NA_real_)
    max(0, full_r2 - summary(fit_red)$r.squared)
  }
  
  r2_vals <- setNames(lapply(valid_covars, get_delta_r2), valid_covars)
  
  model_data <- df %>%
    select(expression, all_of(valid_covars))
  
  ctrl <- trainControl(
    method = "cv",
    number = 10,
    savePredictions = TRUE
  )
  
  model_lm <- tryCatch(
    train(
      expression ~ .,
      data = model_data,
      method = "lm",
      trControl = ctrl,
      metric = "Rsquared"
    ),
    error = function(e) NULL
  )
  
  if (is.null(model_lm)) return(NULL)
  
  res_row <- data.frame(
    ethnicity = target_ethnicity,
    comparison = "Rural_vs_Urban",
    vsv = vsv_name,
    species = vsv_species,
    N = nrow(model_data),
    Full_Model_R2 = round(full_r2, 4),
    CV_Mean_Rsquared = round(model_lm$results$Rsquared, 4),
    CV_Rsquared_SD = round(sd(model_lm$resample$Rsquared, na.rm = TRUE), 4),
    CV_Mean_RMSE = round(model_lm$results$RMSE, 4),
    CV_RMSE_SD = round(sd(model_lm$resample$RMSE, na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  )
  
  for (nm in names(r2_vals)) {
    res_row[[paste0("R2_", make.names(nm))]] <- round(
      as.numeric(r2_vals[[nm]]),
      4
    )
  }
  
  res_row
}


if (nrow(significant_results) > 0) {
  vsv_scaled_df <- vsv_scaled %>%
    rownames_to_column("sample_id")
  
  cv_list <- vector("list", nrow(significant_results))
  
  cat("\nRunning 10-fold CV for significant ethnicity-specific Rural vs Urban associations...\n")
  pb <- txtProgressBar(min = 0, max = nrow(significant_results), style = 3)
  
  for (i in seq_len(nrow(significant_results))) {
    cv_list[[i]] <- run_cv_within_ethnicity(
      significant_results$ethnicity[i],
      significant_results$vsv[i],
      combined_metadata,
      vsv_scaled_df
    )
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  cv_results <- bind_rows(cv_list)
  
  if (nrow(cv_results) > 0) {
    write.xlsx(
      cv_results,
      file.path(results_dir, "ethnicity_specific_Rural_vs_Urban_vsv_cv_marginal_r2.xlsx"),
      rowNames = FALSE
    )
  }
}

cat("\nAnalysis completed.\n")
cat("Ethnicities with valid Rural/Urban analyses:", length(all_results), "\n")
cat("Significant ethnicity-specific Rural vs Urban vSV associations:",
    nrow(significant_results), "\n")
