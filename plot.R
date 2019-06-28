#Load required libraries
library(CB2)
library(tibble)
library(magrittr)
library(dplyr)
library(ggplot2)

#Read input arguments (counts table and sample map)
args = commandArgs(trailingOnly=TRUE)
file_1<-read.table(args[1],header=TRUE)
file_2<-read.table(args[2],header=TRUE)

#Define functions to get counts
counts_design <- function (sgcount, df_design)
{
  df_design$group <- factor(df_design$group,levels=unique(df_design$group))
  cols <- colnames(sgcount)
  sgcount %>% as.data.frame(stringsAsFactors = F) %>%
    tibble::rownames_to_column("sgRNA") %>%
    tidyr::gather_(key_col = "sample_name", value_col = "count", gather_cols = cols) %>%
    dplyr::left_join(df_design, by = "sample_name")
}

plot_cd <- function (sgcount, df_design)
{
  lvl<-levels(df_design$group)
  cdata <- counts_design(sgcount, df_design)
  cdata$count <- 1 + log2(cdata$count)
  cdata$sample_name <- factor(cdata$sample_name,levels=unique(cdata$sample_name))
  df_design$group <- factor(df_design$group,levels=unique(df_design$group))

  ggplot2::ggplot(data = cdata, ggplot2::aes_string(x = "count")) +
  ggplot2::scale_fill_brewer(palette="Accent") +
  ggplot2::geom_density(ggplot2::aes_string(fill = "group")) +
  ggplot2::facet_wrap(~sample_name, ncol = 1) +
  ggplot2::xlab("log2(1+count)") +
  ggplot2::ggtitle(paste("Count Density Plot for ",lvl[2]," vs ",lvl[1],sep="")) +
  ggplot2::theme(plot.title = element_text(hjust = 0.5))
}

plot_heatmap <- function (sgcount, df_design, cor_method = "spearman")
{
  df_design$group <- factor(df_design$group,levels=unique(df_design$group))

  ann_colors <- list(group = c("#7FC97F","#BEAED4"))
  names(ann_colors$group) = as.character(levels(df_design$group))
  sgcount %>% cor(method = cor_method) %>%
  pheatmap::pheatmap(display_numbers = T, number_format = "%.2f", number_color = "snow1",color = colorRampPalette(c("royalblue4", "firebrick1"))(100), annotation_colors = ann_colors, annotation_col = df_design %>% tibble::column_to_rownames("sample_name") %>% dplyr::select_("group"))
}


data=read.table(file_1,header=TRUE)
df_design=read.table(file_2,header=TRUE)

rownames(data)<-data[,1]
data<-data[,-c(1:2)]
data <- data[as.character(df_design$sample_name)]

#Filename for outputs
lvl<-levels(design_mat$group)
count_dist_file=paste("temp/plots/count_distribution_",lvl[2],"_vs_",lvl[1],".png",sep="")
corr_matrix_file=paste("temp/plots/corr_matrix_",lvl[2],"_vs_",lvl[1],".png",sep="")

#Store outputs as pngs
png(count_dist_file)
plot_cd(data,df_design)
dev.off()

png(corr_matrix_file)
plot_heatmap(data, df_design)
dev.off()
