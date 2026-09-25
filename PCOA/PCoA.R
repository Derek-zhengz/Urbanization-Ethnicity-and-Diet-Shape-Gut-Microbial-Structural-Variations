library(vegan)
library(ggplot2)
library(patchwork)
library(dplyr)

## Plotting settings ------------------------
groups <- read.csv(
  "group.csv",
  skip = 1,                 # skip the actual header (row 1)
  header = FALSE            # do not use row 2 as column names automatically
) %>% 
  setNames(c("sample", "group", paste0("V", 3:ncol(.)))) %>%  # name only first two columns
  select(sample = sample, group = group) %>%                  # select only needed columns
  filter(!is.na(sample), sample != "",
         !is.na(group),  group  != "")

## PCoA analysis with vegan ------------------------
vSV <- read.csv('vsv.csv', header = TRUE, 
                row.names = 1, sep = ",", comment.char = "", check.names = FALSE)

dSV <- read.csv('dsv.csv', header = TRUE, 
                row.names = 1, sep = ",", comment.char = "", check.names = FALSE)

# Intersect samples
common <- Reduce(intersect, list(rownames(vSV), rownames(dSV), groups$sample))
vSV <- vSV[common, ]
dSV <- dSV[common, ]
groups <- groups %>% slice(match(common, sample))

# VSV scaling (replace NA with 0)
vSV_scaled <- as.data.frame(apply(vSV, 2, function(x) {
  x_non_na <- x[!is.na(x)]
  if(length(x_non_na) > 0) {
    rng <- range(x_non_na)
    if(rng[2] != rng[1]) {
      x_scaled <- (x - rng[1]) / (rng[2] - rng[1])
    } else {
      x_scaled <- rep(0, length(x)) # constant columns set to 0
    }
    x_scaled[is.na(x_scaled)] <- 0 # replace NA with 0
    return(x_scaled)
  } else {
    return(rep(0, length(x))) # all-NA columns set to 0
  }
}))
rownames(vSV_scaled) <- rownames(vSV)

# DSV processing (replace NA with 0)
dSV[is.na(dSV)] <- 0

# Remove all-zero rows
valid_rows <- (rowSums(vSV_scaled, na.rm = TRUE) > 0) & 
  (rowSums(dSV, na.rm = TRUE) > 0)
vSV_scaled <- vSV_scaled[valid_rows, ]
dSV <- dSV[valid_rows, ]
groups <- groups %>% slice(match(rownames(vSV_scaled), sample))

## 5. Distance calculation ---------------------------------------------------------
d_vsv  <- vegdist(vSV_scaled, method = "canberra")
d_dsv  <- vegdist(dSV,        method = "jaccard")
d_comb <- (d_vsv + d_dsv) / 2                # combined distance

## 6. PCoA --------------------------------------------------------------
pcoa_cap <- capscale(d_comb ~ 1, add = TRUE)  # equivalent to PCoA
pcoa_pts <- scores(pcoa_cap, display = "sites")[, 1:2]

# Explained variance
eig <- pcoa_cap$CA$eig
pc1 <- round(100 * eig[1] / sum(eig), 2)
pc2 <- round(100 * eig[2] / sum(eig), 2)

## 7. Build plotting data frame ---------------------------------------------------
points <- tibble(
  sample = rownames(pcoa_pts),
  dim1   = pcoa_pts[, 1],
  dim2   = pcoa_pts[, 2]
) %>%
  left_join(groups, by = "sample")

plotdata <- points  # reuse directly

## 7  PERMANOVA -----------------------------------------
set.seed(123)
adonis_result_dis <- adonis2(d_comb ~ group, data = groups, permutations = 999)
R2    <- adonis_result_dis$R2[1]
pval  <- adonis_result_dis$`Pr(>F)`[1]
stat_text <- sprintf("PERMANOVA:\nR² = %.4f\np = %.4f", R2, pval)


length  <- length(unique(as.character(groups$sample)))
length1 <- length(unique(as.character(groups$group)))
times1  <- length%/%8
res1    <- length%%8
times2  <- length%/%5
res2    <- length%%5


groups <- data.frame(
  sample = rownames(groups),
  group = as.character(groups$group),
  stringsAsFactors = FALSE
)
length_group <- length(unique(groups$group))  # number of groups

# Color settings (dynamically match group count)
mycol <- c("#E2B0BA","#87AFD5","#6EB9C3","#C98B88",
           "#93C89A","#FFCC98","#E1D1BA","#EDAFA9",
           "#0AB0C8","#9781BB","#E8BD65","#E39844","#ADD1E5")
times_col <- length_group %/% length(mycol)
res_col <- length_group %% length(mycol)
col <- c(rep(mycol, times_col), mycol[1:res_col])  # generate color vector matching group count

# Shape settings (dynamically match group count)
length_group <- length(unique(groups$group))
shapes <- c(15:18, 19, 7:14, 0:6) # original shape library
pich <- rep(16, length_group)  # force all groups to use solid circle (16)
names(pich) <- unique(groups$group)

# Legend column settings
if (length_group > 8) {
  legend_ncol <- 2
} else {
  legend_ncol <- 1
}



## PCoA plot using ggplot2 ------------------------
plotdata <- data.frame(
  sample = points$sample,  # use existing sample column from points
  dim1 = points$dim1,      # correct dim1 values
  dim2 = points$dim2,      # correct dim2 values
  group = points$group     # matched group column after merge
)


## p value ------------------------
adonis_result_dis <- adonis2(d_comb ~ group, data = groups, permutations = 999)

R2 = adonis_result_dis$R2[1]
pvalue = adonis_result_dis$`Pr(>F)`[1]

adonis <- paste("PERMANOVA:\nR2 = ", round(R2,4), "\nP-value = ", pvalue)


## 2. Main PCoA plot (reference style) ----------------------------
## PCoA main plot (fully aligned with CCA code format) ----------------------------

pcoa_plot <- ggplot(plotdata, aes(x = dim1, y = dim2)) +
  geom_point(aes(colour = group, shape = group), size = 3) +
  stat_ellipse(aes(x = dim1, y = dim2, color = group), level = 0.95) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  xlab(paste0("PCoA1 (", pc1, "%)")) +
  ylab(paste0("PCoA2 (", pc2, "%)")) +
  scale_colour_manual(values = col) +
  scale_shape_manual(values = pich) +
  
  theme_classic(base_size = 14) +   # start with a clean classic theme
  theme(
    # 1. remove all four borders first
    axis.line = element_blank(),
    
    # 2. then draw only needed borders (bottom / left)
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left   = element_line(colour = "black"),
    
    # 3. remove top/right ticks and text
    axis.ticks.x.top   = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.x.top    = element_blank(),
    axis.text.y.right  = element_blank(),
    
    # 4. keep other original styles
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid       = element_blank(),
    axis.title       = element_text(color = "black", size = 18),
    axis.text        = element_text(colour = "black", size = 16,
                                    margin = margin(0.6, 0.6, 0.6, 0.6, "lines")),
    legend.title     = element_blank(),
    legend.text      = element_text(size = 10),
    legend.key       = element_blank(),
    legend.position  = c(0.95, 0.8),
    legend.background = element_rect(fill = "transparent", colour = NA)
  ) +
  guides(
    col   = guide_legend(ncol = legend_ncol),
    shape = guide_legend(ncol = legend_ncol)
  )

print(pcoa_plot)

## Boxplot drawing (top and right, strictly aligned with CCA format) ------------------------
# Top boxplot (PCoA1 dimension, adjust outlier appearance)
box_top <- ggplot(plotdata, aes(x = group, y = dim1, fill = group)) +
  geom_boxplot(
    show.legend = FALSE,
    # set outlier appearance: semi-transparent gray, smaller size
    outlier.colour = "gray50",    # outlier color
    outlier.size = 1,             # outlier size
    outlier.alpha = 0.5           # outlier alpha
  ) + 
  stat_boxplot(geom = "errorbar", width = 0.1, size = 0.1) +
  geom_jitter(
    show.legend = FALSE, 
    color = "gray50",  # jitter point color
    alpha = 0.5,       # jitter point alpha
    size = 1           # jitter point size
  ) + 
  scale_fill_manual(values = col) +
  coord_flip() + 
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 14)
  )

# Right boxplot (PCoA2 dimension, with similar outlier appearance)
box_right <- ggplot(plotdata, aes(x = group, y = dim2, fill = group)) + 
  geom_boxplot(
    show.legend = FALSE,
    outlier.colour = "gray50",
    outlier.size = 1,
    outlier.alpha = 0.5
  ) +
  stat_boxplot(geom = "errorbar", width = 0.1, size = 0.1) +
  geom_jitter(
    show.legend = FALSE, 
    color = "gray50",
    alpha = 0.5,
    size = 1,
    width = 0.1,  # reduce this value to make points more clustered (default 0.4)
    height = 0    # no vertical jitter
  ) + 
  scale_fill_manual(values = col) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(colour = "black", size = 14, angle = 45, hjust = 1),  # Group labels on x-axis
    axis.text.y = element_blank()  # hide y-axis labels
  )

# Text box (aligned with stat_box logic in CCA)
stat_text <- "PERMANOVA:\nR² = 0.0657\np = 0.0010"  # replace with actual results
stat_box <- ggplot(plotdata, aes(dim1, dim2)) +
  annotate(
    "text",
    x = mean(range(plotdata$dim1)),
    y = mean(range(plotdata$dim2)),
    label = stat_text,
    size = 4,
    hjust = 0.5,
    vjust = 0.5
  ) +
  theme_void()

# Keep layout unchanged (fully replicate CCA layout logic)
final_plot <- box_top + plot_spacer() + 
  pcoa_plot + box_right +
  plot_layout(
    ncol = 2, 
    nrow = 2, 
    widths = c(4, 1),  # main plot width : right boxplot width
    heights = c(1, 4)  # top boxplot height : main plot height
  )
# Statistics text box (aligned with cca_with_stat logic in CCA)
pcoa_with_stat <- pcoa_plot + stat_box + 
  plot_layout(widths = c(3, 1))  # width ratio between main plot and statistics text

## Combine plots (fully replicate CCA format) ------------------------
# Using patchwork layout: left boxplot + main plot + bottom boxplot + statistics text
final_plot <- box_top + 
  (pcoa_with_stat | box_right) +  # main plot + statistics text side by side with right boxplot
  plot_layout(
    ncol = 1,  # overall 1 column (top boxplot on top, bottom combination below)
    heights = c(1, 4)  # height ratio between top boxplot and bottom combination
  ) +
  patchwork::plot_annotation(title = "PCoA Analysis with Group Comparisons")

# Statistics text box (adjust position to top-right, correct variable reference, align with CCA format)
stat_box <- ggplot(plotdata, aes(dim1, dim2)) +
  annotate(  # replace geom_text with annotate
    "text",  # specify annotation type as text
    x = mean(range(plotdata$dim1)),  # correctly reference dim1 in plotdata
    y = mean(range(plotdata$dim2)),  # correctly reference dim2 in plotdata
    label = stat_text,      # text content
    size = 4, 
    hjust = 0.5, 
    vjust = 0.5
  ) +
  theme_bw() +
  xlab("") + ylab("") +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )

# Re-assemble plot (ensure layout logic consistency, fully replicate CCA format)
final_plot <- box_top + stat_box +
  pcoa_plot + box_right +
  plot_layout(
    heights = c(1, 4),  # height ratio between top boxplot and main plot area
    widths = c(4, 1),   # width ratio between main plot area and right boxplot
    ncol = 2,
    nrow = 2
  ) +
  patchwork::plot_annotation(title = NULL)

print(final_plot)

pdf("PCoA_result.pdf", width = 9, height = 8)
print(final_plot)
dev.off()