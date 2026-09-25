library(circlize)
library(tidyr)
library(dplyr)

# 1. 读取 CSV 文件与数据预处理
# check.names = FALSE 可保持原始列名（如 sig vsv）
data <- read.csv("data - 审稿修改.csv", check.names = FALSE)

# 转换数据格式并统一指标名称
data_long <- pivot_longer(
  data, 
  cols = -Ethnicity, 
  names_to = "indicator", 
  values_to = "value"
) %>%
  mutate(
    # 统一指标名称格式（将空格替换为点号，与颜色映射函数一致）
    indicator = gsub(" ", ".", indicator),
    value_cutoff = ifelse(Ethnicity == "Miao" & indicator == "vsv" & value > 500, 500, value),
    is_cutoff = ifelse(Ethnicity == "Miao" & indicator == "vsv" & value > 500, TRUE, FALSE)
  )

# 2. 颜色方案
ethnicity_colors <- c(
  "#F9A6B0",   # 浅粉加深
  "#A0CEF7",   # 浅蓝加深
  "#F9B08A",   # 浅橙加深
  "#F78A8C",   # 浅红加深
  "#A4C6B4",   # 淡绿加深
  "#B3ACA9"    # 浅灰加深
)

get_indicator_color <- function(indicator_name) {
  colors <- c(
    "vsv"     = "#FCC5CE",
    "dsv"     = "#A0CEF7",
    "sig.vsv" = "#F78A8C",
    "sig.dsv" = "#6BAED6"
  )
  return(colors[indicator_name])
}

# 开始导出 PDF
pdf("ethnicity_circos_plot_revised.pdf", width = 12, height = 10)

# 3. 初始化环形图参数
circos.clear()
circos.par(
  start.degree = 90,
  gap.degree = 6,
  track.margin = c(0.01, 0.01),
  cell.padding = c(0, 0, 0, 0),
  points.overflow.warning = FALSE
)

circos.initialize(factors = data$Ethnicity, xlim = c(0, 1))

# 4. 轨道 1：民族标签
circos.track(
  factors = data$Ethnicity,
  ylim = c(0, 1),
  track.height = 0.08,
  bg.col = ethnicity_colors,
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[2] - 0.18,
      CELL_META$sector.index,
      facing = "bending.inside",
      cex = 1.0,
      adj = c(0.5, 0),
      col = "black",
      font = 2
    )
  }
)

# 5. 轨道 2：刻度线和条形图
circos.track(
  factors = data$Ethnicity,
  ylim = c(0, 1),
  track.height = 0.65,
  bg.col = adjustcolor(ethnicity_colors, alpha.f = 0.2),
  panel.fun = function(x, y) {
    ethnicity <- CELL_META$sector.index
    sector_data <- data_long[data_long$Ethnicity == ethnicity, ]
    
    max_value <- max(sector_data$value_cutoff) * 1.15
    
    # 动态添加刻度线
    if (max_value < 100) {
      tick_interval <- ceiling(max_value / 4)
    } else if (max_value < 500) {
      tick_interval <- ceiling(max_value / 40) * 10
    } else {
      tick_interval <- 120
    }
    
    at <- seq(tick_interval, max_value - tick_interval/2, by = tick_interval)
    for (a in at) {
      circos.lines(
        CELL_META$xlim,
        c(a, a) / max_value,
        lty = 3,
        col = "black",
        lwd = 1
      )
      circos.text(
        CELL_META$cell.xlim[1] - mm_h(1.5),
        a / max_value,
        labels = round(a),
        facing = "clockwise",
        adj = c(0.5, 0),
        cex = 0.8,
        col = "black"
      )
      circos.lines(
        CELL_META$cell.xlim[1] + c(-0.015, 0),
        c(a, a) / max_value,
        lwd = 1.2,
        col = "black"
      )
    }
    
    # 绘制条形图
    bar_width <- 0.16
    gap_width <- 0.045
    n_bars <- nrow(sector_data)
    total_width <- n_bars * bar_width + (n_bars - 1) * gap_width
    start_x <- (1 - total_width) / 2
    
    for (i in 1:n_bars) {
      xleft <- start_x + (i - 1) * (bar_width + gap_width)
      xright <- xleft + bar_width
      ytop <- sector_data$value_cutoff[i] / max_value
      
      bar_color <- get_indicator_color(sector_data$indicator[i])
      
      circos.rect(
        xleft = xleft,
        ybottom = 0,
        xright = xright,
        ytop = ytop,
        col = bar_color,
        border = "black",
        lwd = 1.2
      )
      
      # 截断标记
      if (sector_data$is_cutoff[i]) {
        circos.polygon(
          x = c(xleft + bar_width/4, (xleft + xright)/2, xright - bar_width/4),
          y = c(ytop, ytop + 0.025, ytop),
          col = "white",
          border = "black",
          lwd = 1.2
        )
      }
      
      # 数值标签
      if (ytop > 0.06) {
        circos.text(
          x = (xleft + xright)/2,
          y = ytop + 0.015,
          labels = sector_data$value[i],
          cex = 1,
          facing = "inside",
          adj = c(0.5, 0),
          col = "black"
        )
      }
    }
  }
)

# 6. 添加图例
legend_labels <- c("vsv", "dsv", "sig.vsv", "sig.dsv")
legend_colors <- get_indicator_color(legend_labels)

legend(
  x = 1.15,
  y = 0.5,
  legend = legend_labels,
  fill = legend_colors,
  border = "black",
  bty = "n",
  title.adj = 0.5,
  cex = 1.2,
  pt.cex = 1.2,
  xpd = TRUE
)

dev.off()
circos.clear()