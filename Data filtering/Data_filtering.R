library(tidyverse)

# --------------------------
# Part 1: Process VSV/DSV Data
# --------------------------

# 1. Read sample information without antibiotic filtering
urbanisation_data <- read_csv("matched_urbanisation.csv")
sample_ids <- urbanisation_data %>% select(sample_id)

# 2. Read other datasets and match samples
diet_data <- read_csv("diet.csv") %>%
  inner_join(sample_ids, by = "sample_id") %>%
  select(sample_id, ethnicity) %>%
  mutate(ethnicity = as.character(ethnicity))

vs_data <- read_csv("vsv.csv") %>%
  inner_join(sample_ids, by = "sample_id")

ds_data <- read_csv("dsv.csv") %>%
  inner_join(sample_ids, by = "sample_id")

# 3. Filter features by ethnicity-specific prevalence
filter_by_ethnic_prevalence <- function(feature_data, meta_data, threshold = 10) {
  merged_data <- feature_data %>% inner_join(meta_data, by = "sample_id")
  
  if (nrow(merged_data) == 0) stop("Data merging failed. Please check sample_id.")
  
  feature_cols <- setdiff(names(merged_data), c("sample_id", "ethnicity"))
  if (length(feature_cols) == 0) stop("No feature columns were found.")
  
  ethnic_prevalence <- merged_data %>%
    group_by(ethnicity) %>%
    summarise(
      across(
        all_of(feature_cols),
        ~ sum(.x != 0 & !is.na(.x)) / n() * 100
      ),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = -ethnicity,
      names_to = "feature",
      values_to = "prevalence"
    )
  
  features_to_keep <- ethnic_prevalence %>%
    group_by(feature) %>%
    summarise(
      max_prevalence = max(prevalence),
      .groups = "drop"
    ) %>%
    filter(max_prevalence >= threshold) %>%
    pull(feature)
  
  list(
    filtered_data = feature_data %>%
      select(sample_id, all_of(features_to_keep)),
    prevalence_stats = ethnic_prevalence,
    dropped_features = setdiff(feature_cols, features_to_keep)
  )
}

# 4. Apply filtering and save results
vs_result <- filter_by_ethnic_prevalence(vs_data, diet_data, 10)
ds_result <- filter_by_ethnic_prevalence(ds_data, diet_data, 10)

write_csv(vs_result$filtered_data, "final_filtered_vsv.csv")
write_csv(ds_result$filtered_data, "final_filtered_dsv.csv")

write_csv(
  vs_result$prevalence_stats,
  "vsv_ethnic_prevalence_details.csv"
)

write_csv(
  ds_result$prevalence_stats,
  "dsv_ethnic_prevalence_details.csv"
)

write_lines(
  vs_result$dropped_features,
  "dropped_vsv_features.txt"
)

write_lines(
  ds_result$dropped_features,
  "dropped_dsv_features.txt"
)

# 5. Report filtering results
cat("=== VSV filtering results ===\n")
cat(
  "Original features:", ncol(vs_data) - 1,
  "-> Retained:", ncol(vs_result$filtered_data) - 1,
  paste0(" (Removed ", length(vs_result$dropped_features), ")\n")
)

cat("=== DSV filtering results ===\n")
cat(
  "Original features:", ncol(ds_data) - 1,
  "-> Retained:", ncol(ds_result$filtered_data) - 1,
  paste0(" (Removed ", length(ds_result$dropped_features), ")\n\n")
)


# --------------------------
# Part 2: Process Dietary Data
# --------------------------

# 1. Read dietary data and match samples
diet_data_full <- read_csv("matched_dietary_habits.csv") %>%
  inner_join(sample_ids, by = "sample_id") %>%
  rename_with(~ make.names(.), .cols = everything())

# 2. Filter dietary variables by ethnicity-specific prevalence
clean_diet_by_ethnicity <- function(data) {
  ethnic_col <- names(data)[
    grep("ethnic", names(data), ignore.case = TRUE)
  ]
  
  if (length(ethnic_col) == 0) {
    stop(
      "Ethnicity column was not found. Available columns: ",
      paste(names(data), collapse = ", ")
    )
  }
  
  non_diet_cols <- intersect(
    c("sample_id", ethnic_col, "residency", "region"),
    names(data)
  )
  
  ethnic_prevalence <- data %>%
    group_by_at(ethnic_col) %>%
    summarise(
      across(
        -any_of(non_diet_cols),
        ~ sum(.x > 0 & !is.na(.x)) / n() * 100
      ),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = -all_of(ethnic_col),
      names_to = "diet_item",
      values_to = "prevalence"
    )
  
  diet_to_keep <- ethnic_prevalence %>%
    group_by(diet_item) %>%
    summarise(
      keep = any(prevalence >= 1),
      .groups = "drop"
    ) %>%
    filter(keep) %>%
    pull(diet_item)
  
  list(
    filtered_data = data %>%
      select(
        all_of(non_diet_cols),
        all_of(diet_to_keep)
      ),
    prevalence_stats = ethnic_prevalence
  )
}

# 3. Apply filtering and save results
diet_result <- clean_diet_by_ethnicity(diet_data_full)

write_csv(
  diet_result$filtered_data,
  "final_filtered_diet.csv"
)

write_csv(
  diet_result$prevalence_stats,
  "diet_prevalence_by_ethnicity.csv"
)

# 4. Report filtering results
original_diets <- setdiff(
  names(diet_data_full),
  c("sample_id", "ethnicity", "residency", "region")
)

filtered_diets <- setdiff(
  names(diet_result$filtered_data),
  c("sample_id", "ethnicity", "residency", "region")
)

cat("=== Dietary data filtering results ===\n")
cat("Original dietary variables:", length(original_diets), "\n")
cat("Retained dietary variables:", length(filtered_diets), "\n")
cat(
  "Removed dietary variables:",
  length(original_diets) - length(filtered_diets),
  "\n\n"
)

cat("Processing completed successfully.\n")