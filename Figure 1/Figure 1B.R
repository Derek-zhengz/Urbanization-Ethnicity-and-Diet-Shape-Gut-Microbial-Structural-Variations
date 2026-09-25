
# ==============================================================================
# Publication-Quality Circos Plot
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(dplyr)
library(circlize)
library(readxl)

# 1. Read data and inspect column names
subject_info = readxl::read_xlsx("metadata.xlsx")
cat("Column names in the dataset:", colnames(subject_info), "\n")

# 2. Define colors
age_color <- '#E480A4'
bmi_color <- '#1473B2'
sex_color <- c("0" = '#FBD38D', "1" = '#A7F3D0')
ethnicity_color <- c("Hani" = '#E2B0BA', "Han" = '#87AFD5',
                     "Zang" = '#6EB9C3', "Bai" = '#C98B88',
                     "Dai" = '#93C89A', "Miao" = '#FFCC98')
residency_color <- c("Urban" = '#8ECAE6', "Rural" = '#F4ACB7')

# 3. Prepare data
df <-
  data.frame(
    factors = subject_info$sample_id,    # Sample ID
    x = 1, y = 1, subject_info, stringsAsFactors = TRUE
  ) %>%
  mutate(
    # Convert "Unknown" values in Age/BMI to NA
    Age = case_when(
      str_to_lower(trimws(as.character(Age))) == "unknown" ~ NA_real_,
      TRUE ~ as.numeric(Age)
    ),
    BMI = case_when(
      str_to_lower(trimws(as.character(BMI))) == "unknown" ~ NA_real_,
      TRUE ~ as.numeric(BMI)
    ),
    # Convert residency to character
    residency = as.character(residency)
  ) %>%
  filter(!is.na(Age) & !is.na(BMI)) %>%  # Remove samples with missing Age/BMI
  arrange(Age, BMI) %>%
  slice(c(1:50, (nrow(.)-49):nrow(.))) %>%  # Select the first and last 50 samples
  mutate(factors = factor(factors, levels = factors))  # Preserve sample order

# 4. Define plotting functions
plot_circos_track <- function(temp_value, color, ylab, ylim_mult = c(0.8, 1.1), track_height = 0.2) {
  value_range <- range(temp_value, na.rm = TRUE)
  value_max <- value_range[2] * ylim_mult[2]
  
  circos.track(
    factors = df$factors,
    y = temp_value,
    ylim = value_range * ylim_mult,
    bg.border = "black",
    track.height = track_height,
    panel.fun = function(x, y) {
      name = get.cell.meta.data("sector.index")
      i = get.cell.meta.data("sector.numeric.index")
      xlim = get.cell.meta.data("xlim")
      
      # 1. Draw y-axis ticks
      circos.yaxis(
        side = "left",
        at = c(ceiling(0.8*value_range[1]),
               round(mean(value_range), 0),
               round(value_range[2], 0)),
        sector.index = get.all.sector.index()[1],
        labels.cex = 0.8,
        labels.niceFacing = FALSE
      )
      
      # 2. Draw horizontal lines
      circos.lines(
        x = mean(xlim, na.rm = TRUE),
        y = temp_value[i],
        type = "h",
        col = color,
        lwd = 2
      )
      
      # 3. Add sample ID labels
      if (ylab == "Age") {
        circos.text(
          x = 1,
          y = value_max + (value_range[2] * 0.2),
          labels = name,
          facing = "clockwise",
          niceFacing = TRUE,
          cex = 1,
          col = "black"
        )
      }
      
      # 4. Draw data points
      circos.points(
        x = mean(xlim),
        y = temp_value[i],
        pch = 16,
        cex = 0.8,
        col = color
      )
    }
  )
}

# Function for categorical tracks
plot_categorical_track <- function(temp_var, colors, track_height = 0.1) {
  temp_var[is.na(temp_var) | !temp_var %in% names(colors)] <- "unmatched"
  colors_extended <- c(colors, "unmatched" = "#666666")
  mapped_colors <- colors_extended[temp_var]
  
  circos.track(
    factors = df$factors,
    y = df$y,
    ylim = c(0, 1),
    bg.border = "black",
    track.height = track_height,
    panel.fun = function(x, y) {
      i = get.cell.meta.data("sector.numeric.index")
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      
      circos.rect(
        xleft = xlim[1],
        ybottom = ylim[1],
        xright = xlim[2],
        ytop = ylim[2],
        col = mapped_colors[i],
        border = "black"
      )
    }
  )
}

# 5. Open the PDF device and initialize Circos
pdf("circos_plot.pdf", width = 12, height = 12, pointsize = 10)

circos.clear()
circos.par(
  "track.height" = 0.2,
  start.degree = 90,
  clock.wise = TRUE,
  gap.after = c(rep(0, nrow(df) - 1), 90),
  circle.margin = c(0.1, 0.1, 0.1, 0.1),
  cell.padding = c(0, 0, 0, 0)
)

# Use a fixed global x-axis range to avoid zero-width sector errors
circos.initialize(factors = df$factors, xlim = c(0.5, 1.5))

# 6. Draw all tracks
plot_circos_track(df$Age, age_color, "Age")                    # Age
plot_circos_track(df$BMI, bmi_color, "BMI")                    # BMI
plot_categorical_track(df$Sex, sex_color)                      # Sex
plot_categorical_track(df$ethnicity, ethnicity_color)          # Ethnicity
plot_categorical_track(df$residency, residency_color)          # Residency

# 7. Close the PDF device
dev.off()

# Print output file location
cat("Circos plot PDF exported successfully!\n")
cat("File path:", file.path(getwd(), "circos_plot.pdf"), "\n")
```
