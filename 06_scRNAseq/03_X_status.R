 ##*****************************##
##  Allelome.PRO Violin plots  ##
##*****************************##
## Project: dissectX
## Tim Hasenbein 
## Last modification 08.2023
## Creation: 04.2022 
## Sort cells according to X status
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
input <- "/./"
output <- "./"
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
# filter AP results, X genes with >= 10 reads
filter <- function(x) {
  x <- x[x$chr == "chrX",]
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
srat <- readRDS(paste0(rds)) .
cells <- srat@meta.data$cell
ap_x <- res%>%dplyr::filter(barcode %in% cells) .


######------  filter bl6 Xi and cast Xi cells  ------###### 
bl6_ko <- ap_x%>%dplyr::filter(GT == "KO" & allelic_ratio <= 0.3) 
cast_ko <- ap_x%>%dplyr::filter(GT == "KO" & allelic_ratio >= 0.7) 
react_ko <- ap_x%>%dplyr::filter(GT == "KO" & allelic_ratio < 0.7 & allelic_ratio > 0.3) 
bl6_wt <- ap_x%>%dplyr::filter(GT == "WT" & allelic_ratio <= 0.3) 
cast_wt <- ap_x%>%dplyr::filter(GT == "WT" & allelic_ratio >= 0.7) 
react_wt <- ap_x%>%dplyr::filter(GT == "WT" & allelic_ratio < 0.7 & allelic_ratio > 0.3) 
write.table(bl6_ko,paste0(output,"bl6_+-_Xa.txt"),quote = F, row.names = F, sep="\t")
write.table(cast_ko,paste0(output,"cast_+-_Xa.txt"),quote = F, row.names = F, sep="\t")
write.table(react_ko,paste0(output,"react_+-.txt"),quote = F, row.names = F, sep="\t")
write.table(bl6_wt,paste0(output,"bl6_++_Xa.txt"),quote = F, row.names = F, sep="\t")
write.table(cast_wt,paste0(output,"cast_++_Xa.txt"),quote = F, row.names = F, sep="\t")
write.table(react_wt,paste0(output,"react_++.txt"),quote = F, row.names = F, sep="\t")


######------ add Xa status ------###### 
meta <- srat@meta.data
meta <- meta%>%mutate(Xa=ifelse(cell %in% bl6_ko$barcode, "bl6_Xa",
                                ifelse(cell %in% bl6_wt$barcode, "bl6_Xa_wt",
                                       ifelse(cell %in% cast_ko$barcode, "cast_Xa",
                                              ifelse(cell %in% cast_wt$barcode, "cast_Xa_wt",
                                                     ifelse(cell %in% react_wt$barcode, "wt_react",
                                                            ifelse(cell %in% react_ko$barcode, "ko_react",
                                                                   "not_assigned")))))))
srat@meta.data <- meta


######------ Split UMAP plot ------###### 
pdf(file = paste0(output,"03_umap_split.pdf"),width=15,height=5)
DimPlot(srat, reduction = "umap", label = TRUE, label.size = 2, split.by = "Xa")
dev.off()




