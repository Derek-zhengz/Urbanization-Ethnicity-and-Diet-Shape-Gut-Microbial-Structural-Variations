
library(readxl)
library(dplyr)
library(tidyr)
library(UpSetR)

#=============================
# 1. Ethnicity order and colors
#=============================

ethnicity_order <- c("Bai", "Dai", "Han", "Hani", "Miao", "Zang")

ethnicity_colors <- c(
  "Bai"  = "#E2B0BA",
  "Dai"  = "#87AFD5",
  "Han"  = "#6EB9C3",
  "Hani" = "#C98B88",
  "Miao" = "#93C89A",
  "Zang" = "#FFCC98"
)

#=============================
# 2. Read data
#=============================

vsv <- read_excel("all_ethnicities_significant_Rural_vs_Urban_vsv.xlsx")
dsv <- read_excel("all_ethnicities_significant_Rural_vs_Urban_dsv.xlsx")

#=============================
# 3. Prepare vSV data
#=============================

vsv2 <- vsv %>%
  transmute(
    ethnicity = ethnicity,
    SV_ID = vsv,
    SV_type = "vSV",
    species = species,
    p_adj_bh = p_adj
  )

#=============================
# 4. Prepare dSV data
#=============================

dsv2 <- dsv %>%
  transmute(
    ethnicity = ethnicity,
    SV_ID = dsv,
    SV_type = "dSV",
    species = species,
    p_adj_bh = p_adj
  )

#=============================
# 5. Combine datasets
#=============================

sv_all <- bind_rows(vsv2, dsv2) %>%
  filter(
    !is.na(SV_ID),
    !is.na(ethnicity),
    ethnicity %in% ethnicity_order
  ) %>%
  mutate(
    # Prevent accidental overlap between identical vSV and dSV names
    SV_unique_ID = paste(SV_type, SV_ID, sep = "__")
  ) %>%
  distinct(ethnicity, SV_unique_ID, .keep_all = TRUE)

#=============================
# 6. Build SV sets for each ethnicity
#=============================

sv_list <- lapply(
  ethnicity_order,
  function(x) {
    sv_all %>%
      filter(ethnicity == x) %>%
      pull(SV_unique_ID) %>%
      unique()
  }
)

names(sv_list) <- ethnicity_order

# Check the number of SVs in each ethnicity
sapply(sv_list, length)

#=============================
# 7. Convert to UpSetR format
#=============================

upset_data <- fromList(sv_list)

#=============================
# 8. Export to PDF
#=============================

pdf(
  "Urban_rural_SV_ethnicity_UpSet_colored.pdf",
  width = 9,
  height = 6
)

upset(
  upset_data,
  
  # Ethnicity order
  sets = ethnicity_order,
  keep.order = TRUE,
  
  # Assign a different color to each ethnicity
  sets.bar.color = unname(
    ethnicity_colors[ethnicity_order]
  ),
  
  # Intersection bar color
  main.bar.color = "#666666",
  
  # Matrix point color
  matrix.color = "#444444",
  
  # Text size
  text.scale = c(
    1.4,  # Intersection size title
    1.2,  # Intersection size tick labels
    1.3,  # Set size title
    1.1,  # Set size tick labels
    1.3,  # Set names
    1.1   # Numbers above bars
  ),
  
  # Maximum y-axis value for the main bar plot
  mainbar.y.max = 10,
  
  # Sort intersections by frequency
  order.by = "freq",
  
  # Keep all existing intersections
  empty.intersections = NULL,
  
  # Display numbers above bars
  show.numbers = "yes",
  
  # Matrix point and line sizes
  point.size = 3.5,
  line.size = 0.8,
  
  # Relative size of the main bar and set size panels
  mb.ratio = c(0.65, 0.35),
  
  # Set size x-axis label
  sets.x.label = "Number of urban/rural-associated SVs",
  
  # Intersection y-axis label
  mainbar.y.label = "Intersection size"
)

dev.off()
```
