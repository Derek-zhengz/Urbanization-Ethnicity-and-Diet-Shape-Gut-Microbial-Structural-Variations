# =============================================================================
# Ethnicity -> Diet -> SV mediation analysis
#  A / B / C / C' path
# Primary: Coverage + covariates
# Sensitivity: Primary + species abundance
# dSV: logistic regression
# vSV: linear regression
# =============================================================================

set.seed(12345)
library(tidyverse)
library(readxl)
library(openxlsx)
library(mediation)
library(broom)

if (!dir.exists("mediation_results")) dir.create("mediation_results")
N_SIMS <- 5000

# =============================================================================
# 1. 
# =============================================================================

clean_diet <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\\.+", " ")
  x <- stringr::str_replace_all(x, "_+", " ")
  x <- stringr::str_squish(x)
  stringr::str_to_lower(x)
}

clean_sv <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, ":", ".")
  x <- stringr::str_replace_all(x, ";", ".")
  x <- stringr::str_replace_all(x, "\\s+", "")
  x
}

extract_species <- function(x) {
  x <- as.character(x)
  x <- sub("__dup[0-9]+$", "", x)
  out <- stringr::str_extract(x, "^[A-Za-z0-9]+_[A-Za-z0-9]+(?=[\\.:])")
  ifelse(is.na(out), stringr::str_extract(x, "^[A-Za-z0-9]+_[A-Za-z0-9]+"), out)
}

process_covariate <- function(data, name) {
  data$sample_id <- as.character(data$sample_id)
  num <- data %>%
    dplyr::select(-dplyr::any_of(c("sample_id","ethnicity","residency","region"))) %>%
    dplyr::select(where(is.numeric))
  if (ncol(num) == 0) {
    out <- data.frame(sample_id = data$sample_id, score = NA_real_)
  } else {
    out <- data.frame(sample_id = data$sample_id, score = rowSums(num, na.rm = TRUE))
  }
  names(out)[2] <- paste0(name, "_total_score")
  out
}


reverse_ordinal <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  if (all(is.na(x_num))) return(x)
  rng <- range(x_num, na.rm = TRUE)
  rng[1] + rng[2] - x_num
}

prepare_general_metadata <- function(data) {
  data$sample_id <- trimws(as.character(data$sample_id))
  general_cols_original <- setdiff(names(data), "sample_id")
  general_covars <- paste0("general__", make.names(general_cols_original, unique = TRUE))
  names(data)[match(general_cols_original, names(data))] <- general_covars
  
  data <- data %>% dplyr::mutate(dplyr::across(dplyr::all_of(general_covars), function(col) {
    if (is.numeric(col)) return(col)
    col_char <- trimws(as.character(col))
    col_lower <- tolower(col_char)
    if (all(col_lower %in% c("yes","no","y","n","true","false","1","0","",NA))) {
      return(dplyr::case_when(
        col_lower %in% c("yes","y","true","1") ~ 1,
        col_lower %in% c("no","n","false","0") ~ 0,
        TRUE ~ NA_real_
      ))
    }
    num_col <- suppressWarnings(as.numeric(col_char))
    nonmissing_original <- !is.na(col_char) & col_char != ""
    if (sum(!is.na(num_col)) == sum(nonmissing_original)) return(num_col)
    factor(col_char)
  }))
  
  list(data = data, covars = general_covars)
}

# =============================================================================
# 2. 5significantresults
# =============================================================================

cat("\n========== 1. Read significant results ==========\n")

eth_diet <- readxl::read_excel("all significant results.xlsx") %>% as.data.frame()
eth_dsv  <- readxl::read_excel("all_ethnicities_significant_dSV.xlsx") %>% as.data.frame()
eth_vsv  <- readxl::read_excel("all_ethnicities_significant_vSV.xlsx") %>% as.data.frame()
diet_dsv <- readxl::read_excel("significant_diet_factors_dsv.xlsx") %>% as.data.frame()
diet_vsv <- readxl::read_excel("significant_diet_factors_vsv.xlsx") %>% as.data.frame()

cat("Ethnicity-Diet:", nrow(eth_diet), "\n")
cat("Ethnicity-dSV:", nrow(eth_dsv), "\n")
cat("Ethnicity-vSV:", nrow(eth_vsv), "\n")
cat("Diet-dSV:", nrow(diet_dsv), "\n")
cat("Diet-vSV:", nrow(diet_vsv), "\n")

# =============================================================================
# 3. results
# =============================================================================

eth_diet <- eth_diet %>%
  dplyr::mutate(Ethnicity = trimws(as.character(ethnicity)),
                Diet = as.character(diet_variable),
                Diet_key = clean_diet(diet_variable))

eth_dsv <- eth_dsv %>%
  dplyr::mutate(Ethnicity = trimws(as.character(ethnicity)),
                SV = as.character(dsv), SV_key = clean_sv(dsv),
                Species = if ("species" %in% colnames(.)) as.character(species) else extract_species(dsv))

eth_vsv <- eth_vsv %>%
  dplyr::mutate(Ethnicity = trimws(as.character(ethnicity)),
                SV = as.character(vsv), SV_key = clean_sv(vsv),
                Species = if ("species" %in% colnames(.)) as.character(species) else extract_species(vsv))

diet_dsv <- diet_dsv %>%
  dplyr::mutate(Diet = as.character(diet_factor), Diet_key = clean_diet(diet_factor),
                SV = as.character(dsv), SV_key = clean_sv(dsv),
                Species_diet = if ("species" %in% colnames(.)) as.character(species) else extract_species(dsv))

diet_vsv <- diet_vsv %>%
  dplyr::mutate(Diet = as.character(diet_factor), Diet_key = clean_diet(diet_factor),
                SV = as.character(vsv), SV_key = clean_sv(vsv),
                Species_diet = if ("species" %in% colnames(.)) as.character(species) else extract_species(vsv))

# =============================================================================
# 4. Candidate pathways
# =============================================================================

cat("\n========== 2. Build candidate pathways ==========\n")

candidate_dsv <- eth_diet %>%
  dplyr::select(Ethnicity, Diet, Diet_key) %>% dplyr::distinct() %>%
  dplyr::inner_join(
    eth_dsv %>% dplyr::select(Ethnicity, SV_eth = SV, SV_key, Species) %>% dplyr::distinct(),
    by = "Ethnicity"
  ) %>%
  dplyr::inner_join(
    diet_dsv %>% dplyr::select(Diet_key, SV_key, SV_diet = SV, Species_diet) %>% dplyr::distinct(),
    by = c("Diet_key","SV_key")
  ) %>%
  dplyr::mutate(Species = dplyr::coalesce(Species, Species_diet, extract_species(SV_diet))) %>%
  dplyr::transmute(Ethnicity, Diet, Diet_key, SV = SV_diet, SV_key, Species, SV_Type = "dSV") %>%
  dplyr::distinct()

candidate_vsv <- eth_diet %>%
  dplyr::select(Ethnicity, Diet, Diet_key) %>% dplyr::distinct() %>%
  dplyr::inner_join(
    eth_vsv %>% dplyr::select(Ethnicity, SV_eth = SV, SV_key, Species) %>% dplyr::distinct(),
    by = "Ethnicity"
  ) %>%
  dplyr::inner_join(
    diet_vsv %>% dplyr::select(Diet_key, SV_key, SV_diet = SV, Species_diet) %>% dplyr::distinct(),
    by = c("Diet_key","SV_key")
  ) %>%
  dplyr::mutate(Species = dplyr::coalesce(Species, Species_diet, extract_species(SV_diet))) %>%
  dplyr::transmute(Ethnicity, Diet, Diet_key, SV = SV_diet, SV_key, Species, SV_Type = "vSV") %>%
  dplyr::distinct()

mediation_candidates <- dplyr::bind_rows(candidate_dsv, candidate_vsv)

cat("dSV candidate:", nrow(candidate_dsv), "\n")
cat("vSV candidate:", nrow(candidate_vsv), "\n")
cat("Total:", nrow(mediation_candidates), "\n")

openxlsx::write.xlsx(mediation_candidates,
                     "mediation_results/mediation_candidate_triplets.xlsx",
                     rowNames = FALSE)

if (nrow(mediation_candidates) == 0) stop("No candidate mediation pathways were identified.")

# =============================================================================
# 5. samples
# =============================================================================

cat("\n========== 3. Read sample-level data ==========\n")

metadata <- read.csv("metadata.csv", check.names = FALSE, stringsAsFactors = FALSE)
diet_data <- read.csv("final_filtered_diet.csv", check.names = FALSE, stringsAsFactors = FALSE)

dsv_data <- read.csv(
  "final_filtered_dsv.csv",
  row.names = 1, check.names = FALSE
)

vsv_data <- read.csv(
  "final_filtered_vsv.csv",
  row.names = 1, check.names = FALSE
)

urbanization_data <- read.csv("Urbanisation.csv", check.names = FALSE)
medication_data <- read.csv("Medication.csv", check.names = FALSE)
general_metadata_data <- read.csv("General metadata.csv", check.names = FALSE)

metadata$sample_id <- trimws(as.character(metadata$sample_id))
diet_data$sample_id <- trimws(as.character(diet_data$sample_id))
rownames(dsv_data) <- trimws(as.character(rownames(dsv_data)))
rownames(vsv_data) <- trimws(as.character(rownames(vsv_data)))

# =============================================================================
# 6. Species abundance
# =============================================================================

cat("\n========== 4. Species abundance ==========\n")

species_map <- read.csv("sv_to_metaphlan_final_mapping.csv", check.names = FALSE)
metaphlan_raw <- readr::read_tsv(
  "yunnan_7_years_ago_combined_metaphlan_profile.tsv",
  comment = "", show_col_types = FALSE
)

if (grepl("#", colnames(metaphlan_raw)[1]))
  colnames(metaphlan_raw)[1] <- "clade_name"

abundance_matrix <- metaphlan_raw %>%
  dplyr::filter(clade_name %in% unique(species_map$metaphlan_clade)) %>%
  dplyr::left_join(species_map, by = c("clade_name" = "metaphlan_clade")) %>%
  dplyr::select(-dplyr::any_of(c("clade_name","metaphlan_species"))) %>%
  tidyr::pivot_longer(cols = -sv_species, names_to = "sample_id", values_to = "abundance") %>%
  dplyr::mutate(sample_id = trimws(as.character(sample_id)),
                abundance = suppressWarnings(as.numeric(abundance))) %>%
  dplyr::group_by(sample_id, sv_species) %>%
  dplyr::summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(id_cols = sample_id, names_from = sv_species,
                     values_from = abundance, values_fill = 0)

# =============================================================================
# 7. Species coverage
# =============================================================================

cat("\n========== 5. Species coverage ==========\n")

coverage0 <- read.csv("merged_species_mean_coverage_matrix.csv",
                      check.names = FALSE, stringsAsFactors = FALSE)
colnames(coverage0)[1] <- "sample_id"

coverage_matrix <- coverage0 %>%
  tidyr::pivot_longer(cols = -sample_id, names_to = "sv_species", values_to = "coverage") %>%
  dplyr::mutate(sample_id = trimws(as.character(sample_id)),
                coverage = suppressWarnings(as.numeric(coverage))) %>%
  dplyr::group_by(sample_id, sv_species) %>%
  dplyr::summarise(coverage = mean(coverage, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(id_cols = sample_id, names_from = sv_species,
                     values_from = coverage, values_fill = 0)

colnames(coverage_matrix)[-1] <- paste0(colnames(coverage_matrix)[-1], "_coverage")

# =============================================================================
# 8. Environmental covariates
# =============================================================================

cat("\n========== 6. Covariates ==========\n")

medication_total <- process_covariate(medication_data, "medication")

animal_col <- names(urbanization_data)[tolower(trimws(names(urbanization_data))) == "contact with animals"]
if (length(animal_col) == 1) {
  urbanization_data[[animal_col]] <- reverse_ordinal(urbanization_data[[animal_col]])
  cat("Contact with Animals reverse-coded before Urbanization composite.\n")
} else if (length(animal_col) == 0) {
  stop("Column 'Contact with Animals' was not found in Urbanisation.csv.")
} else {
  stop("Multiple columns matching 'Contact with Animals' were found in Urbanisation.csv.")
}

urbanization_total <- process_covariate(urbanization_data, "urbanization")
general_prepared <- prepare_general_metadata(general_metadata_data)
general_cov <- general_prepared$data
general_covars <- general_prepared$covars

# =============================================================================
# 9. Diet map + master data
# =============================================================================

diet_exclude <- c("sample_id","ethnicity","residency","region")

diet_column_map <- data.frame(
  Diet_col = setdiff(colnames(diet_data), diet_exclude),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(Diet_key = clean_diet(Diet_col)) %>%
  dplyr::distinct(Diet_key, .keep_all = TRUE)

metadata_base <- metadata %>%
  dplyr::select(-dplyr::any_of(intersect(diet_column_map$Diet_col, colnames(metadata))))

diet_join <- diet_data %>%
  dplyr::select(sample_id, dplyr::all_of(diet_column_map$Diet_col))

sample_data <- metadata_base %>%
  dplyr::inner_join(diet_join, by = "sample_id") %>%
  dplyr::inner_join(medication_total, by = "sample_id") %>%
  dplyr::inner_join(urbanization_total, by = "sample_id") %>%
  dplyr::inner_join(general_cov, by = "sample_id") %>%
  dplyr::left_join(abundance_matrix, by = "sample_id") %>%
  dplyr::left_join(coverage_matrix, by = "sample_id") %>%
  dplyr::mutate(
    ethnicity = trimws(as.character(ethnicity)),
    Sex = as.factor(Sex),
    residency = as.factor(residency)
  )

cat("sample_data N =", nrow(sample_data), "\n")

# =============================================================================
# 10. SV map
# =============================================================================

dsv_map <- data.frame(SV_col = colnames(dsv_data), stringsAsFactors = FALSE) %>%
  dplyr::mutate(SV_key = clean_sv(SV_col)) %>%
  dplyr::distinct(SV_key, .keep_all = TRUE)

vsv_map <- data.frame(SV_col = colnames(vsv_data), stringsAsFactors = FALSE) %>%
  dplyr::mutate(SV_key = clean_sv(SV_col)) %>%
  dplyr::distinct(SV_key, .keep_all = TRUE)

# =============================================================================
# 11. Single mediation analysis
# =============================================================================

run_one_mediation <- function(eth, diet_name, diet_key, sv_key, species,
                              sv_type, add_abundance = FALSE, sims = N_SIMS) {
  
  diet_col <- diet_column_map %>%
    dplyr::filter(Diet_key == .env$diet_key) %>%
    dplyr::pull(Diet_col)
  if (length(diet_col) == 0) return(NULL)
  diet_col <- diet_col[1]
  
  if (sv_type == "dSV") {
    sv_col <- dsv_map %>% dplyr::filter(SV_key == .env$sv_key) %>% dplyr::pull(SV_col)
    sv_matrix <- dsv_data
  } else {
    sv_col <- vsv_map %>% dplyr::filter(SV_key == .env$sv_key) %>% dplyr::pull(SV_col)
    sv_matrix <- vsv_data
  }
  
  if (length(sv_col) == 0) return(NULL)
  sv_col <- sv_col[1]
  
  abundance_candidates <- c(species, paste0("s__",species), gsub("_"," ",species))
  abundance_col <- abundance_candidates[abundance_candidates %in% colnames(sample_data)][1]
  
  coverage_candidates <- c(paste0(species,"_coverage"),
                           paste0("s__",species,"_coverage"),
                           paste0(gsub("_"," ",species),"_coverage"))
  coverage_col <- coverage_candidates[coverage_candidates %in% colnames(sample_data)][1]
  
  if (is.na(coverage_col)) return(NULL)
  if (add_abundance && is.na(abundance_col)) return(NULL)
  
  sv_df <- data.frame(
    sample_id = rownames(sv_matrix),
    outcome = suppressWarnings(as.numeric(sv_matrix[[sv_col]])),
    stringsAsFactors = FALSE
  )
  
  df <- sample_data %>%
    dplyr::inner_join(sv_df, by = "sample_id") %>%
    dplyr::mutate(
      ethnicity_binary = ifelse(ethnicity == eth, 1, 0),
      mediator = suppressWarnings(as.numeric(.data[[diet_col]])),
      species_coverage = suppressWarnings(as.numeric(.data[[coverage_col]]))
    )
  
  if (add_abundance)
    df$species_abundance <- suppressWarnings(as.numeric(df[[abundance_col]]))
  
  candidate_covars <- c("Age","Sex","BMI","residency",
                        "medication_total_score","urbanization_total_score",
                        general_covars)
  
  covars <- candidate_covars[candidate_covars %in% colnames(df)]
  
  needed <- c("outcome","ethnicity_binary","mediator","species_coverage",covars)
  if (add_abundance) needed <- c(needed,"species_abundance")
  
  df <- df %>% tidyr::drop_na(dplyr::all_of(needed))
  if (sv_type == "dSV") df <- df %>% dplyr::filter(outcome %in% c(0,1))
  
  if (nrow(df) < 30 || length(unique(df$ethnicity_binary)) < 2 ||
      length(unique(df$mediator)) < 2 || length(unique(df$outcome)) < 2) return(NULL)
  
  df <- droplevels(df)
  
  covars <- covars[sapply(df[covars], function(z) length(unique(z[!is.na(z)])) > 1)]
  covars <- covars[sapply(covars, function(v) {
    z <- df[[v]]
    if (is.factor(z) || is.character(z)) {
      tab <- table(z)
      length(tab) > 1 && min(tab) >= 5
    } else {
      length(unique(z)) > 1
    }
  })]
  
  eth_tab <- table(df$ethnicity_binary)
  if (length(eth_tab) < 2 || any(eth_tab < 5)) return(NULL)
  
  if (sv_type == "dSV") {
    ytab <- table(df$outcome)
    if (length(ytab) < 2 || any(ytab < 5)) return(NULL)
  }
  
  # ===========================================================================
  # A path: Ethnicity -> Diet
  # ===========================================================================
  
  mediator_formula <- stats::as.formula(
    paste("mediator ~", paste(c("ethnicity_binary",covars), collapse = " + "))
  )
  
  model_m <- tryCatch(
    stats::lm(mediator_formula, data = df),
    error = function(e) NULL
  )
  if (is.null(model_m)) return(NULL)
  
  a_result <- broom::tidy(model_m) %>%
    dplyr::filter(term == "ethnicity_binary")
  if (nrow(a_result) == 0) return(NULL)
  
  A_beta <- a_result$estimate[1]
  A_P <- a_result$p.value[1]
  
  # ===========================================================================
  # C path: Ethnicity -> SV，Diet
  # outcome model
  # ===========================================================================
  
  c_rhs <- c("ethnicity_binary","species_coverage",covars)
  if (add_abundance) c_rhs <- c(c_rhs,"species_abundance")
  
  c_formula <- stats::as.formula(
    paste("outcome ~", paste(c_rhs, collapse = " + "))
  )
  
  if (sv_type == "dSV") {
    model_c <- tryCatch(
      stats::glm(c_formula, data = df, family = stats::binomial()),
      error = function(e) NULL
    )
  } else {
    model_c <- tryCatch(
      stats::lm(c_formula, data = df),
      error = function(e) NULL
    )
  }
  
  if (is.null(model_c)) return(NULL)
  
  c_result <- broom::tidy(model_c) %>%
    dplyr::filter(term == "ethnicity_binary")
  if (nrow(c_result) == 0) return(NULL)
  
  C_beta <- c_result$estimate[1]
  C_SE <- c_result$std.error[1]
  C_P <- c_result$p.value[1]
  
  # ===========================================================================
  # B + C' model: Ethnicity + Diet -> SV
  # ===========================================================================
  
  outcome_rhs <- c("ethnicity_binary","mediator","species_coverage",covars)
  if (add_abundance) outcome_rhs <- c(outcome_rhs,"species_abundance")
  
  outcome_formula <- stats::as.formula(
    paste("outcome ~", paste(outcome_rhs, collapse = " + "))
  )
  
  if (sv_type == "dSV") {
    model_y <- tryCatch(
      stats::glm(outcome_formula, data = df, family = stats::binomial()),
      error = function(e) NULL
    )
  } else {
    model_y <- tryCatch(
      stats::lm(outcome_formula, data = df),
      error = function(e) NULL
    )
  }
  
  if (is.null(model_y)) return(NULL)
  
  b_result <- broom::tidy(model_y) %>%
    dplyr::filter(term == "mediator")
  
  cprime_result <- broom::tidy(model_y) %>%
    dplyr::filter(term == "ethnicity_binary")
  
  if (nrow(b_result) == 0 || nrow(cprime_result) == 0) return(NULL)
  
  B_beta <- b_result$estimate[1]
  B_SE <- b_result$std.error[1]
  B_P <- b_result$p.value[1]
  
  Cprime_beta <- cprime_result$estimate[1]
  Cprime_SE <- cprime_result$std.error[1]
  Cprime_P <- cprime_result$p.value[1]
  
  # dSV logisticOR；vSVOR
  B_OR <- ifelse(sv_type == "dSV", exp(B_beta), NA_real_)
  B_OR_low <- ifelse(sv_type == "dSV", exp(B_beta - 1.96*B_SE), NA_real_)
  B_OR_high <- ifelse(sv_type == "dSV", exp(B_beta + 1.96*B_SE), NA_real_)
  
  C_OR <- ifelse(sv_type == "dSV", exp(C_beta), NA_real_)
  C_OR_low <- ifelse(sv_type == "dSV", exp(C_beta - 1.96*C_SE), NA_real_)
  C_OR_high <- ifelse(sv_type == "dSV", exp(C_beta + 1.96*C_SE), NA_real_)
  
  Cprime_OR <- ifelse(sv_type == "dSV", exp(Cprime_beta), NA_real_)
  Cprime_OR_low <- ifelse(sv_type == "dSV", exp(Cprime_beta - 1.96*Cprime_SE), NA_real_)
  Cprime_OR_high <- ifelse(sv_type == "dSV", exp(Cprime_beta + 1.96*Cprime_SE), NA_real_)
  
  # ===========================================================================
  # Mediation
  # ===========================================================================
  
  med <- tryCatch(
    mediation::mediate(
      model.m = model_m,
      model.y = model_y,
      treat = "ethnicity_binary",
      mediator = "mediator",
      treat.value = 1,
      control.value = 0,
      boot = FALSE,
      sims = sims
    ),
    error = function(e) {
      cat("\nMEDIATION ERROR:", eth,"|",diet_name,"|",sv_type,"|",sv_col,"\n")
      cat("Error:", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (is.null(med)) return(NULL)
  
  ACME <- if (!is.null(med$d.avg)) med$d.avg else mean(c(med$d0,med$d1),na.rm=TRUE)
  ACME_P <- if (!is.null(med$d.avg.p)) med$d.avg.p else mean(c(med$d0.p,med$d1.p),na.rm=TRUE)
  ACME_CI <- if (!is.null(med$d.avg.ci)) med$d.avg.ci else
    c(mean(c(med$d0.ci[1],med$d1.ci[1]),na.rm=TRUE),
      mean(c(med$d0.ci[2],med$d1.ci[2]),na.rm=TRUE))
  
  ADE <- if (!is.null(med$z.avg)) med$z.avg else mean(c(med$z0,med$z1),na.rm=TRUE)
  ADE_P <- if (!is.null(med$z.avg.p)) med$z.avg.p else mean(c(med$z0.p,med$z1.p),na.rm=TRUE)
  ADE_CI <- if (!is.null(med$z.avg.ci)) med$z.avg.ci else
    c(mean(c(med$z0.ci[1],med$z1.ci[1]),na.rm=TRUE),
      mean(c(med$z0.ci[2],med$z1.ci[2]),na.rm=TRUE))
  
  PROP <- if (!is.null(med$n.avg)) med$n.avg else mean(c(med$n0,med$n1),na.rm=TRUE)
  PROP_P <- if (!is.null(med$n.avg.p)) med$n.avg.p else mean(c(med$n0.p,med$n1.p),na.rm=TRUE)
  PROP_CI <- if (!is.null(med$n.avg.ci)) med$n.avg.ci else c(NA_real_,NA_real_)
  
  data.frame(
    Ethnicity = eth,
    Comparison = paste0(eth,"_vs_Others"),
    Diet = diet_name,
    SV = sv_col,
    Species = species,
    SV_Type = sv_type,
    Model = ifelse(add_abundance,"Coverage + Species abundance","Coverage"),
    N = nrow(df),
    N_TargetEthnicity = sum(df$ethnicity_binary == 1),
    N_Others = sum(df$ethnicity_binary == 0),
    
    A_path_beta = as.numeric(A_beta),
    A_path_P = as.numeric(A_P),
    
    B_path_beta = as.numeric(B_beta),
    B_path_P = as.numeric(B_P),
    B_path_OR_dSV = as.numeric(B_OR),
    B_path_OR_low_dSV = as.numeric(B_OR_low),
    B_path_OR_high_dSV = as.numeric(B_OR_high),
    
    C_path_beta = as.numeric(C_beta),
    C_path_P = as.numeric(C_P),
    C_path_OR_dSV = as.numeric(C_OR),
    C_path_OR_low_dSV = as.numeric(C_OR_low),
    C_path_OR_high_dSV = as.numeric(C_OR_high),
    
    Cprime_path_beta = as.numeric(Cprime_beta),
    Cprime_path_P = as.numeric(Cprime_P),
    Cprime_path_OR_dSV = as.numeric(Cprime_OR),
    Cprime_path_OR_low_dSV = as.numeric(Cprime_OR_low),
    Cprime_path_OR_high_dSV = as.numeric(Cprime_OR_high),
    
    ACME = as.numeric(ACME),
    ACME_CI_low = as.numeric(ACME_CI[1]),
    ACME_CI_high = as.numeric(ACME_CI[2]),
    ACME_P = as.numeric(ACME_P),
    
    ADE = as.numeric(ADE),
    ADE_CI_low = as.numeric(ADE_CI[1]),
    ADE_CI_high = as.numeric(ADE_CI[2]),
    ADE_P = as.numeric(ADE_P),
    
    Total_Effect = as.numeric(med$tau.coef),
    Total_Effect_CI_low = as.numeric(med$tau.ci[1]),
    Total_Effect_CI_high = as.numeric(med$tau.ci[2]),
    Total_Effect_P = as.numeric(med$tau.p),
    
    Prop_Mediated = as.numeric(PROP),
    Prop_Mediated_CI_low = as.numeric(PROP_CI[1]),
    Prop_Mediated_CI_high = as.numeric(PROP_CI[2]),
    Prop_Mediated_P = as.numeric(PROP_P),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# 12. Primary
# =============================================================================

cat("\n========== 7. Primary mediation ==========\n")

primary_list <- vector("list", nrow(mediation_candidates))

for (i in seq_len(nrow(mediation_candidates))) {
  x <- mediation_candidates[i, ]
  cat(sprintf("[Primary %d/%d] %s -> %s -> %s | %s\n",
              i,nrow(mediation_candidates),x$Ethnicity,x$Diet,x$SV_Type,x$SV))
  flush.console()
  
  primary_list[[i]] <- run_one_mediation(
    eth = x$Ethnicity,
    diet_name = x$Diet,
    diet_key = x$Diet_key,
    sv_key = x$SV_key,
    species = x$Species,
    sv_type = x$SV_Type,
    add_abundance = FALSE,
    sims = N_SIMS
  )
}

mediation_primary <- dplyr::bind_rows(primary_list)

if (nrow(mediation_primary) > 0) {
  mediation_primary <- mediation_primary %>%
    dplyr::mutate(
      ACME_FDR = p.adjust(ACME_P, method = "BH"),
      Mediation_Significant = !is.na(ACME_FDR) & ACME_FDR < 0.05 &
        ((ACME_CI_low > 0 & ACME_CI_high > 0) |
           (ACME_CI_low < 0 & ACME_CI_high < 0))
    ) %>%
    dplyr::arrange(ACME_FDR)
}

openxlsx::write.xlsx(
  mediation_primary,
  "mediation_results/mediation_primary_coverage_adjusted_with_ABCCprime.xlsx",
  rowNames = FALSE
)

# =============================================================================
# 13. Sensitivity
# =============================================================================

cat("\n========== 8. Sensitivity: + species abundance ==========\n")

sensitivity_list <- vector("list", nrow(mediation_candidates))

for (i in seq_len(nrow(mediation_candidates))) {
  x <- mediation_candidates[i, ]
  cat(sprintf("[Sensitivity %d/%d] %s -> %s -> %s | %s\n",
              i,nrow(mediation_candidates),x$Ethnicity,x$Diet,x$SV_Type,x$SV))
  flush.console()
  
  sensitivity_list[[i]] <- run_one_mediation(
    eth = x$Ethnicity,
    diet_name = x$Diet,
    diet_key = x$Diet_key,
    sv_key = x$SV_key,
    species = x$Species,
    sv_type = x$SV_Type,
    add_abundance = TRUE,
    sims = N_SIMS
  )
}

mediation_sensitivity <- dplyr::bind_rows(sensitivity_list)

if (nrow(mediation_sensitivity) > 0) {
  mediation_sensitivity <- mediation_sensitivity %>%
    dplyr::mutate(
      ACME_FDR = p.adjust(ACME_P, method = "BH"),
      Mediation_Significant = !is.na(ACME_FDR) & ACME_FDR < 0.05 &
        ((ACME_CI_low > 0 & ACME_CI_high > 0) |
           (ACME_CI_low < 0 & ACME_CI_high < 0))
    ) %>%
    dplyr::arrange(ACME_FDR)
}

openxlsx::write.xlsx(
  mediation_sensitivity,
  "mediation_results/mediation_sensitivity_abundance_coverage_with_ABCCprime.xlsx",
  rowNames = FALSE
)

# =============================================================================
# 14. Significant mediation paths
# =============================================================================

sig_primary <- if (nrow(mediation_primary) > 0) {
  mediation_primary %>% dplyr::filter(Mediation_Significant)
} else data.frame()

sig_sensitivity <- if (nrow(mediation_sensitivity) > 0) {
  mediation_sensitivity %>% dplyr::filter(Mediation_Significant)
} else data.frame()

openxlsx::write.xlsx(
  sig_primary,
  "mediation_results/significant_mediation_primary_with_ABCCprime.xlsx",
  rowNames = FALSE
)

openxlsx::write.xlsx(
  sig_sensitivity,
  "mediation_results/significant_mediation_sensitivity_with_ABCCprime.xlsx",
  rowNames = FALSE
)

# =============================================================================
# 15. Separate dSV / vSV outputs
# =============================================================================

if (nrow(mediation_primary) > 0) {
  openxlsx::write.xlsx(
    mediation_primary %>% dplyr::filter(SV_Type == "dSV"),
    "mediation_results/mediation_primary_dSV_with_ABCCprime.xlsx",
    rowNames = FALSE
  )
  
  openxlsx::write.xlsx(
    mediation_primary %>% dplyr::filter(SV_Type == "vSV"),
    "mediation_results/mediation_primary_vSV_with_ABCCprime.xlsx",
    rowNames = FALSE
  )
}

# =============================================================================
# 16. Primary vs Sensitivity
# =============================================================================

if (nrow(mediation_primary) > 0 && nrow(mediation_sensitivity) > 0) {
  
  comparison <- mediation_primary %>%
    dplyr::select(
      Ethnicity,Diet,SV,Species,SV_Type,
      N_primary = N,
      A_primary = A_path_beta,
      B_primary = B_path_beta,
      C_primary = C_path_beta,
      Cprime_primary = Cprime_path_beta,
      ACME_primary = ACME,
      ACME_FDR_primary = ACME_FDR,
      ADE_primary = ADE,
      Total_Effect_primary = Total_Effect,
      Prop_Mediated_primary = Prop_Mediated
    ) %>%
    dplyr::full_join(
      mediation_sensitivity %>%
        dplyr::select(
          Ethnicity,Diet,SV,Species,SV_Type,
          N_sensitivity = N,
          A_sensitivity = A_path_beta,
          B_sensitivity = B_path_beta,
          C_sensitivity = C_path_beta,
          Cprime_sensitivity = Cprime_path_beta,
          ACME_sensitivity = ACME,
          ACME_FDR_sensitivity = ACME_FDR,
          ADE_sensitivity = ADE,
          Total_Effect_sensitivity = Total_Effect,
          Prop_Mediated_sensitivity = Prop_Mediated
        ),
      by = c("Ethnicity","Diet","SV","Species","SV_Type")
    )
  
  openxlsx::write.xlsx(
    comparison,
    "mediation_results/mediation_primary_vs_sensitivity_with_ABCCprime.xlsx",
    rowNames = FALSE
  )
}
# =============================================================================
# Nominal mediation signals: ACME P < 0.05 and 95% CI does not cross zero
# =============================================================================

nominal_primary <- mediation_primary %>%
  filter(
    !is.na(ACME_P),
    ACME_P < 0.05,
    !is.na(ACME_CI_low),
    !is.na(ACME_CI_high),
    (ACME_CI_low > 0 & ACME_CI_high > 0) |
      (ACME_CI_low < 0 & ACME_CI_high < 0)
  ) %>%
  arrange(ACME_FDR)

nominal_sensitivity <- mediation_sensitivity %>%
  filter(
    !is.na(ACME_P),
    ACME_P < 0.05,
    !is.na(ACME_CI_low),
    !is.na(ACME_CI_high),
    (ACME_CI_low > 0 & ACME_CI_high > 0) |
      (ACME_CI_low < 0 & ACME_CI_high < 0)
  ) %>%
  arrange(ACME_FDR)

cat("\nNominal mediation signals in Primary analysis:", nrow(nominal_primary), "\n")
cat("Nominal mediation signals in Sensitivity analysis:", nrow(nominal_sensitivity), "\n")

openxlsx::write.xlsx(
  nominal_primary,
  "mediation_results/nominal_mediation_primary_ACME_P_lt_0.05.xlsx",
  rowNames = FALSE
)

openxlsx::write.xlsx(
  nominal_sensitivity,
  "mediation_results/nominal_mediation_sensitivity_ACME_P_lt_0.05.xlsx",
  rowNames = FALSE
)

# =============================================================================
# 17. Summary output
# =============================================================================

cat("\n============================================================\n")
cat("Mediation analysis completed\n")
cat("============================================================\n")

cat("Candidate paths:", nrow(mediation_candidates), "\n")
cat("dSV candidates:", nrow(candidate_dsv), "\n")
cat("vSV candidates:", nrow(candidate_vsv), "\n")
cat("Primary models completed:", nrow(mediation_primary), "\n")
cat("Primary significant mediation paths:", nrow(sig_primary), "\n")
cat("Sensitivity models completed:", nrow(mediation_sensitivity), "\n")
cat("Sensitivity significant mediation paths:", nrow(sig_sensitivity), "\n")

if (nrow(mediation_primary) > 0) {
  cat("\nPrimary results preview:\n")
  
  print(
    mediation_primary %>%
      dplyr::select(
        Ethnicity,Diet,SV_Type,SV,N,
        A_path_beta,A_path_P,
        B_path_beta,B_path_P,
        C_path_beta,C_path_P,
        Cprime_path_beta,Cprime_path_P,
        ACME,ACME_CI_low,ACME_CI_high,ACME_P,ACME_FDR,
        ADE,Total_Effect,Prop_Mediated,Mediation_Significant
      ) %>%
      head(30)
  )
}

cat("\nResults directory: mediation_results/\n")
