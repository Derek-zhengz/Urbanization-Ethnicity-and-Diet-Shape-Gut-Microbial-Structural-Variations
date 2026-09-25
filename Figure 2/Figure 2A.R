
library(vegan)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(tibble)
library(grid)

## 1. Data loading and filtering ------------------------
# Read grouping information, including sample_id, ethnicity group, and residency
groups <- read.csv(
  'group.csv',
  header = TRUE,
  row.names = 1,  # Use the first column as row names (sample_id)
  sep = ",",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Ensure sample_id is retained as a column
if(!"sample_id" %in% colnames(groups)) {
  groups <- groups %>% rownames_to_column("sample_id")
}

# Inspect unique values in the residency column
unique_residency <- unique(groups$residency)
print("Unique values in the residency column:")
print(unique_residency)

# Retain Urban/Rural samples across all ethnic groups
filtered <- groups %>%
  mutate(
    residency = trimws(residency)  # Remove extra spaces
  ) %>%
  filter(
    residency %in% c("urban", "rural", "Urban", "Rural")
  ) %>%
  mutate(residency = tolower(residency))  # Standardize to lowercase


## 2. Read SV data and match samples ------------------------
# Read vSV and dSV data (rows = samples, columns = SV features)
vSV <- read.csv(
  'final_filtered_vsv.csv',
  header = TRUE, row.names = 1, sep = ",", comment.char = "", check.names = FALSE
)
dSV <- read.csv(
  'final_filtered_dsv.csv',
  header = TRUE, row.names = 1, sep = ",", comment.char = "", check.names = FALSE
)

# Identify samples present in vSV, dSV, and the filtered Urban/Rural sample set
common <- Reduce(intersect, list(
  rownames(vSV),
  rownames(dSV),
  filtered$sample_id
))

# Subset SV data using common samples and preserve sample order
vSV <- vSV[common, ]
dSV <- dSV[common, ]

# Match grouping data to SV data and preserve the same sample order
matched <- filtered %>%
  filter(sample_id %in% common) %>%
  arrange(match(sample_id, common))

# Verify sample matching
cat("\nSample matching check:\n")
cat("Number of vSV samples:", nrow(vSV), "\n")
cat("Number of dSV samples:", nrow(dSV), "\n")
cat("Number of matched Urban/Rural samples:", nrow(matched), "\n")


## 3. Data cleaning ------------------------
# vSV processing: replace NA with 0 without normalization
vSV_clean <- vSV
vSV_clean[is.na(vSV_clean)] <- 0
rownames(vSV_clean) <- rownames(vSV)

# dSV processing: replace NA with 0
dSV_clean <- dSV
dSV_clean[is.na(dSV_clean)] <- 0

# Remove all-zero samples
valid_rows <- (rowSums(vSV_clean, na.rm = TRUE) > 0) &
  (rowSums(dSV_clean, na.rm = TRUE) > 0)

vSV_clean <- vSV_clean[valid_rows, ]
dSV_clean <- dSV_clean[valid_rows, ]

# Synchronize grouping data after removing all-zero samples
matched <- matched %>%
  filter(sample_id %in% rownames(vSV_clean)) %>%
  arrange(match(sample_id, rownames(vSV_clean)))

# Verify cleaned data
cat("\nData cleaning check:\n")
cat("Number of vSV samples after cleaning:", nrow(vSV_clean), "\n")
cat("Number of dSV samples after cleaning:", nrow(dSV_clean), "\n")
cat("Number of matched grouping samples:", nrow(matched), "\n")


## 4. PCoA and distance calculation ------------------------
# Calculate distance matrices
d_vsv  <- vegdist(vSV_clean, method = "canberra")
d_dsv  <- vegdist(dSV_clean, method = "jaccard")
d_comb <- (d_vsv + d_dsv) / 2

# Perform PCoA using capscale
pcoa_cap <- capscale(d_comb ~ 1, add = TRUE)
pcoa_pts <- scores(pcoa_cap, display = "sites")[, 1:2]

# Calculate explained variance for the first two PCoA axes
eig <- pcoa_cap$CA$eig
pc1 <- round(100 * eig[1] / sum(eig), 2)
pc2 <- round(100 * eig[2] / sum(eig), 2)

cat("\nExplained variance of the first two PCoA axes:\n")
cat("PCoA1:", pc1, "%, PCoA2:", pc2, "%\n")


## 5. Build plotting data and perform statistical analysis ------------------------
# Combine PCoA coordinates with grouping information
points <- tibble(
  sample_id = rownames(pcoa_pts),
  dim1 = pcoa_pts[, 1],
  dim2 = pcoa_pts[, 2]
) %>%
  left_join(matched, by = "sample_id")

# PERMANOVA analysis for residency differences
set.seed(123)
adonis_result <- adonis2(
  d_comb ~ residency,
  data = matched,
  permutations = 999
)

# Extract PERMANOVA statistics
R2 <- adonis_result$R2[1]
pval <- adonis_result$`Pr(>F)`[1]
stat_text <- sprintf("PERMANOVA:\nR² = %.4f\np = %.4f", R2, pval)


## 6. Plotting settings ------------------------
# Avoid using built-in function names as variable names
n_samples <- length(unique(matched$sample_id))
n_residency <- length(unique(matched$residency))

# Grouping information used for plotting
plot_groups <- matched %>%
  select(sample_id, residency) %>%
  rename(group = residency)

# Color settings
mycol <- c("#E2B0BA","#87AFD5","#6EB9C3","#C98B88",
           "#93C89A","#FFCC98","#E1D1BA","#EDAFA9",
           "#0AB0C8","#9781BB","#E8BD65","#E39844","#ADD1E5")

times_col <- n_residency %/% length(mycol)
res_col <- n_residency %% length(mycol)
col <- c(rep(mycol, times_col), mycol[1:res_col])

# Shape settings
pich <- rep(16, n_residency)
names(pich) <- unique(plot_groups$group)

# Number of legend columns
legend_ncol <- ifelse(n_residency > 8, 2, 1)


## 7. Plotting data frame ------------------------
plotdata <- points %>%
  select(sample_id, dim1, dim2, residency) %>%
  rename(group = residency)

# Main PCoA plot
pcoa_plot <- ggplot(plotdata, aes(x = dim1, y = dim2)) +
  geom_point(aes(colour = group, shape = group), size = 2.5) +
  stat_ellipse(aes(color = group), level = 0.95, linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.5) +
  xlab(paste0("PCoA1 (", pc1, "%)")) +
  ylab(paste0("PCoA2 (", pc2, "%)")) +
  scale_colour_manual(values = col) +
  scale_shape_manual(values = pich) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0),
    axis.title = element_text(color = "black", size = 14),
    axis.text = element_text(colour = "black", size = 12,
                             margin = margin(0.6, 0.6, 0.6, 0.6, "lines")),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key = element_blank(),
    legend.position = c(0.92, 0.85),
    legend.background = element_rect(fill = "transparent", colour = NA)
  ) +
  guides(
    col = guide_legend(ncol = legend_ncol),
    shape = guide_legend(ncol = legend_ncol)
  )

## Top boxplot
box_top <- ggplot(plotdata, aes(x = group, y = dim1, fill = group)) +
  geom_boxplot(
    show.legend = FALSE,
    outlier.colour = "gray50",
    outlier.size = 0.8,
    outlier.alpha = 0.5,
    linewidth = 0.5
  ) +
  stat_boxplot(geom = "errorbar", width = 0.1, linewidth = 0.4) +
  geom_jitter(show.legend = FALSE, color = "gray50", alpha = 0.5, size = 0.8) +
  scale_fill_manual(values = col) +
  coord_flip() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 13)
  )

## Right-side boxplot
box_right <- ggplot(plotdata, aes(x = group, y = dim2, fill = group)) +
  geom_boxplot(
    show.legend = FALSE,
    outlier.colour = "gray50",
    outlier.size = 0.8,
    outlier.alpha = 0.5,
    linewidth = 0.5
  ) +
  stat_boxplot(geom = "errorbar", width = 0.1, linewidth = 0.4) +
  geom_jitter(show.legend = FALSE, color = "gray50", alpha = 0.5, size = 0.8, width = 0.1) +
  scale_fill_manual(values = col) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0),
    axis.title = element_blank(),
    axis.text.x = element_text(colour = "black", size = 13, angle = 45, hjust = 1),
    axis.text.y = element_blank()
  )

## Statistics text box
stat_box <- ggplot(plotdata, aes(dim1, dim2)) +
  annotate(
    "rect",
    xmin = mean(range(plotdata$dim1)) - 0.18 * diff(range(plotdata$dim1)),
    xmax = mean(range(plotdata$dim1)) + 0.18 * diff(range(plotdata$dim1)),
    ymin = mean(range(plotdata$dim2)) - 0.18 * diff(range(plotdata$dim2)),
    ymax = mean(range(plotdata$dim2)) + 0.18 * diff(range(plotdata$dim2)),
    fill = "white",
    color = "black",
    linewidth = 0.3,
    alpha = 0.8
  ) +
  annotate(
    "text",
    x = mean(range(plotdata$dim1)),
    y = mean(range(plotdata$dim2)),
    label = stat_text,
    size = 4.2,
    hjust = 0.5,
    vjust = 0.5
  ) +
  theme_void()

## Combine plots
final_plot <- box_top + stat_box +
  pcoa_plot + box_right +
  plot_layout(
    heights = c(1.2, 4),
    widths = c(4, 1.2),
    ncol = 2,
    nrow = 2
  ) +
  plot_annotation(title = NULL)

# Display the combined plot
print(final_plot)

# Save as PDF
pdf("PCoA_AllEthnic_residency_vSV_noNorm.pdf", width = 10, height = 8.5)
print(final_plot)
dev.off()
```
