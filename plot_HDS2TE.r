# Ensure required packages are installed and loaded

# 加载所需的库
# ggpubr 用于创建美观的散点图并添加回归线和相关系数
library(ggpubr)
# optparse 用于处理命令行参数，使脚本更灵活
library(optparse)

# 定义命令行选项
option_list <- list(
  # 定义输入HDS密度文件的选项
  make_option(c("--inHDS"), type="character", default=NULL,
              help="Input file for HDS density (e.g., HDS.density)",
              metavar="character"),
  # 定义输入TE密度文件的选项
  make_option(c("--inTE"), type="character", default=NULL,
              help="Input file for TE density (e.g., TE.density)",
              metavar="character"),
  # 定义输出PDF文件名的选项
  make_option(c("--out"), type="character", default="HDS2TEs.pdf",
              help="Output PDF file name for the plot (default: HDS2TEs.pdf)",
              metavar="character")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 检查必要的输入文件是否提供
# 修正了 is.is.null 为 is.null
if (is.null(opt$inHDS) || is.null(opt$inTE)){
  print_help(opt_parser)
  stop("Both --inHDS and --inTE arguments must be supplied.", call.=FALSE)
}

# 读取HDS注释文件
# header=F 表示文件没有标题行
# stringsAsFactors=F 避免将字符串自动转换为因子
HDS_annotation <- read.table(opt$inHDS, header = F, stringsAsFactors = F)
# 设置HDS注释数据的列名
colnames(HDS_annotation) <- c("Chr", "Start", "End", "Value")

# 读取TE注释文件
# header=F 表示文件没有标题行
# stringsAsFactors=F 避免将字符串自动转换为因子
TE_annotation <- read.table(opt$inTE, header = F, stringsAsFactors = F)
# 设置TE注释数据的列名
colnames(TE_annotation) <- c("Chr", "Start", "End", "Value")

# 为HDS和TE数据创建副本，方便后续操作
aa <- HDS_annotation
bb <- TE_annotation

# 为HDS数据创建唯一ID，通过连接染色体、起始和结束位置
aa$ID <- paste(aa$Chr, aa$Start, aa$End, sep = "-")
# 为TE数据创建唯一ID，通过连接染色体、起始和结束位置
bb$ID <- paste(bb$Chr, bb$Start, bb$End, sep = "-")

# 创建HDS值的临时数据框
tmp1.df <- data.frame("ID" = aa$ID, "notalign.value" = aa$Value)
# 创建TE值的临时数据框
tmp2.df <- data.frame("ID" = bb$ID, "repeat.value" = bb$Value)

# 根据ID合并两个数据框
# by="ID" 表示以ID列作为合并键
cor.df <- merge(tmp1.df, tmp2.df, by = "ID")

# 使用ggscatter创建散点图
P1 <- ggscatter(cor.df, x = "notalign.value", y = "repeat.value", size = 1.5,
                add = "reg.line", # 添加回归线
                color = "#449945", # 设置点和回归线的颜色
                conf.int = T, # 显示置信区间
                linewidth = 0.8, # 设置回归线宽度
                cor.coef = TRUE, # 显示相关系数
                cor.coeff.args = list(method = "pearson")) + # 指定皮尔逊相关系数
  # 设置X轴和Y轴的标签
  labs(x = "The distribution of HDS",
       y = "The distribution of TEs") +
  # 设置X轴的连续刻度范围
  scale_x_continuous(limits = c(0, 1)) +
  # 设置Y轴的连续刻度范围
  scale_y_continuous(limits = c(0, 1)) +
  # 使用经典主题
  theme_classic() +
  # 自定义主题元素
  theme(axis.ticks = element_blank(), # 移除坐标轴刻度线
        axis.text.x = element_text(colour = "black", size = 12), # X轴文本颜色和大小
        axis.text.y = element_text(colour = "black", size = 12, face = "plain"), # Y轴文本颜色、大小和字体
        axis.title.y = element_text(colour = "black", size = 12, face = "plain"), # Y轴标题颜色、大小和字体
        axis.title.x = element_text(colour = "black", size = 12, face = "plain"), # X轴标题颜色、大小和字体
        panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid"), # 面板边框
        panel.grid.major = element_line(size = 0.5, color = "grey80")) # 主要网格线

# 保存图表为PDF文件
# file 指定输出文件名，从命令行参数获取
# width 和 height 设置输出PDF的宽度和高度
ggsave(P1, file = opt$out, width = 5, height = 4.8)

