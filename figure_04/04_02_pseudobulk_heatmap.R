##*************************##
##  Allelome.PRO scRNAseq  ##
##*************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 05.2023
## Creation: 05.2023
## Visualize result by plotting the AR 
library(tidyverse)
library(ggplot2)
library(data.table)
library(Seurat)
library(RColorBrewer)
library(gplots)
library(ComplexHeatmap)
library(writexl)


######------  set environment  ------###### 
input <- "./"
output <- "./"


######------  get data  ------###### 
read <- function(x) read.table(x, header=T)
dataFiles <- lapply(Sys.glob(paste0(input,"*_Xa_annotation_us.bed_1/locus_table.txt")), read)
# name each data frame
sample_names <- Sys.glob(paste0(input,"*_Xa_annotation_us.bed_1/locus_table.txt"))
sample_names <- gsub("_annotation_us.bed_1/locus_table.txt","",sample_names)
sample_names <- gsub(paste0(input,"2023_01_16_"),"",sample_names)
names(dataFiles) <- sample_names


######------  filter data  ------###### 
filter <- function(x) { 
  x <- x[x$chr == "chrX",]
  x <- x[x$total_reads >= 30,]
}
dataFiles_filter <- lapply(dataFiles, filter)
rm(dataFiles,sample_names)


######------  generate matrix ------###### 
res <- bind_rows(dataFiles_filter, .id = "sample")
res$start <- as.numeric(res$start)
res$end <- as.numeric(res$end)
res <- res%>%arrange(start)
res <- res%>%select(name,sample,allelic_ratio)
res <- res%>%pivot_wider(names_from=sample, values_from = allelic_ratio, values_fill = NA)
res <- res%>%drop_na()
write_xlsx(res, paste0(output,"AR_scRNA_pseudobulk.xlsx"), col_names = TRUE)


######------  plot ARs ------###### 
options(scipen = 9999999)
plot <- res%>%column_to_rownames("name")
plot <- plot%>%select(WT_cast_Xa,KO_cast_Xa,WT_bl6_Xa,KO_bl6_Xa)
pdf(file=paste0(output,"Pseudobulk_AP_genes_heatmap_new.pdf"),width=5, height=10)
colfunc <- circlize::colorRamp2(c(0, 0.5, 1), c("black","white","#B87333")) 
ht <- Heatmap(as.matrix(plot), 
              cluster_rows = F, 
              cluster_columns = F, 
              column_title = NULL, 
              show_row_names = T, 
              row_names_gp =  gpar(fontsize = 4),
              col = colfunc, 
              width = unit(10, "cm"), border = T,
              row_gap = unit(0.5, "mm"), 
              row_title_gp = gpar(fontsize = 6, col="black"),
              column_names_gp = gpar(fontsize = 10),
              heatmap_legend_param = list(at = c(0,0.25,0.5,0.75,1))) 
draw(ht, heatmap_legend_side = "left")
dev.off()


######------  plot AR deltas ------###### 
delta <- res
res$delta_Xi <- abs(res$WT_cast_Xa - res$KO_cast_Xa)
res$delta_Xa <- abs(res$WT_bl6_Xa - res$KO_bl6_Xa)
plot_delta <- res%>%select(name,delta_Xi,delta_Xa)
plot_delta <- plot_delta%>%column_to_rownames("name")
colfunc <- circlize::colorRamp2(c(0, 0.05,0.3), c("white","white","red")) 
pdf(file=paste0(output,"Pseudobulk_AP_delta_genes_heatmap_new.pdf"),width=5, height=10)
ht <- Heatmap(as.matrix(plot_delta), 
              cluster_rows = F, 
              cluster_columns = F, 
              column_title = NULL, 
              show_row_names = T, 
              col = colfunc,
              row_names_gp =  gpar(fontsize = 4),
              width = unit(10, "cm"), border = T,
              row_gap = unit(0.5, "mm"), 
              row_title_gp = gpar(fontsize = 6, col="black"),
              column_names_gp = gpar(fontsize = 10)) 
draw(ht, heatmap_legend_side = "left")
dev.off()




