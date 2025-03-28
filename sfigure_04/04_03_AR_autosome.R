##*****************************##
##  Allelome.PRO Violin plots  ##
##*****************************##
## Project: dissectX
## Tim Hasenbein 
## Last modification 08.2023
## Creation: 04.2022 
## Violin plot for AR of autosomal genes
library(tidyverse) 
library(ggplot2)
library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(gdata)
library(ggpubr)
library(ggrepel)
library(circlize)
library(ggbeeswarm)


######------  set environment  ------###### 
input <- "./02_AllelomePRO/"
output <- "./02_supplementary_figure_4/"
rds <- "./03_srat_clean_integrated_reg_MT_labelled_XX_KO_WT.rds"


######------  get data of AP run  ------###### 
read <- function(x) read.table(x, header=T)
dataFiles <- lapply(Sys.glob(paste0(input,"*_chr_annotation_mm10.bed_1/locus_table.txt")), read)
# name each data frame
sample_names <- Sys.glob(paste0(input,"*_chr_annotation_mm10.bed_1/locus_table.txt"))
sample_names <- gsub("_chr_annotation_mm10.bed_1/locus_table.txt","",sample_names)
sample_names <- gsub(paste0(input,"2023_07_24_"),"",sample_names)
sample_names <- gsub(paste0(input,"2023_07_25_"),"",sample_names)
names(dataFiles) <- sample_names
# filter AP results, genes with >= 10 reads
filter <- function(x) {
  x <- x[x$total_reads >= 10,]
}
dataFiles_filter <- lapply(dataFiles, filter)
dataFiles_filter <- dataFiles_filter[which(lapply(dataFiles_filter, nrow) != 0)]
rm(dataFiles,sample_names)
# select allelic ratio 
select <- function(x) {
  x <- x[,c("name","allelic_ratio","total_reads")]
}
dataFiles_filter_select <- lapply(dataFiles_filter, select)
res <- bind_rows(dataFiles_filter_select, .id = "sample")
res <- res%>%separate(sample, into=c("barcode","GT"),sep="_")
rm(dataFiles_filter,dataFiles_filter_select)


######------  filter cells that pass QC  ------###### 
srat <- readRDS(paste0(rds)) # 3694 cells
cells <- srat@meta.data$cell
ap_x <- res%>%dplyr::filter(barcode %in% cells) 


######------   violin plot combined ------###### 
plot <- ap_x%>%dplyr::filter(name!="chrX")
dodge <- position_dodge(width = 1)
pdf(file=paste0(output,"AR_autosomes_combined.pdf"),width=12, height=5)
ggplot(plot, aes(x=name, y=allelic_ratio,fill=GT)) + 
  geom_hline(yintercept=0.7, linetype="dashed", color = "black") + 
  geom_hline(yintercept=0.3, linetype="dashed", color = "black") +
 #geom_quasirandom(size = 0.1,color="darkgrey",position = dodge) +
  geom_violin(width = 0.5, size = 0.1, scale = "width", position = dodge,color="black") +
  labs(x="Chromosomes",y="Allelic ratio") +
  geom_boxplot(width = 0.3, size=0.2, alpha = 0.2, color = "white", outlier.shape=NA,position = dodge) +
  scale_x_discrete(limits=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                            "chr17","chr18","chr19"), guide = guide_axis(angle = 0), position = "bottom") +
  scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0,0.3,0.5,0.7,1)) +
  scale_fill_manual(values = c(alpha("darkgrey",alpha=0.6),alpha("#93CDB5",alpha=0.6))) +
  theme(
    axis.title.x=element_blank(),
    axis.text.x = element_text(angle = 45, hjust=0.1),
    axis.line = element_line(colour = 'darkgrey', size = 0.2),
    axis.ticks.y=element_line(colour = "grey40",size=0.2),
    axis.title.y=element_text(size=12),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    legend.position="none") 
dev.off()
