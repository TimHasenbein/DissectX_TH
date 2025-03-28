##***********************##
##  Correlation matrix   ##
##***********************##
## Project: dissectX
## Tim Hasenbein
## Last modification 01.2022
## Creation: 10.2021
## plot heatmap of unsupervised clustering of a pearson correlation matrix from expression data, to confirm the expected relationships
library(tidyverse)
library(ggplot2)
library(biomaRt)
library(data.table)
library(stats)
library(pheatmap)
library(RColorBrewer)
library(GenomicFeatures)
library(openxlsx)


######------ Set environment  ------###### 
input <- "./TPMs/"
output <- "./12_supplements/"


######------ Get data and make count matrix  ------######
adult <- fread(paste0(input,"TPM_table_adult_replicates.txt"), sep="\t")
TPM_table <- adult
TPM_table <- TPM_table%>%dplyr::select("chromosome_name","mgi_symbol","start_position","end_position","strand","name","transcript_length",
  
  "Br_DFxDF_++_1","Br_DFxDF_++_2","Br_DFxDF_++_3","Br_XFLxXFL_++_1",
  "Br_DFxDF_--_1","Br_DFxDF_--_2","Br_DFxDF_--_3","Br_DFxDF_--_4",
  "Br_LFDxLFD_--_1","Br_LFDxLFD_--_2","Br_LFDxLFD_--_3",
  
  "He_DFxDF_++_1","He_DFxDF_++_2","He_DFxDF_++_3","He_XFLxXFL_++_1",
  "He_DFxDF_--_1","He_DFxDF_--_2","He_DFxDF_--_3","He_DFxDF_--_4",
  "He_LFDxLFD_--_1","He_LFDxLFD_--_2","He_LFDxLFD_--_3",
                                         
  "Ki_DFxDF_++_1","Ki_DFxDF_++_2","Ki_DFxDF_++_3","Ki_XFLxXFL_++_1",
  "Ki_DFxDF_--_1","Ki_DFxDF_--_2","Ki_DFxDF_--_3","Ki_DFxDF_--_4",
  "Ki_LFDxLFD_--_1","Ki_LFDxLFD_--_2","Ki_LFDxLFD_--_3",
  
  "Li_DFxDF_++_1","Li_DFxDF_++_2","Li_DFxDF_++_3","Li_XFLxXFL_++_1",
  "Li_DFxDF_--_1","Li_DFxDF_--_2","Li_DFxDF_--_3","Li_DFxDF_--_4",                                  
  "Li_LFDxLFD_--_1","Li_LFDxLFD_--_2","Li_LFDxLFD_--_3",
                                         
  "Lu_DFxDF_++_1","Lu_DFxDF_++_2","Lu_DFxDF_++_3","Lu_XFLxXFL_++_1",
  "Lu_DFxDF_--_1","Lu_DFxDF_--_2","Lu_DFxDF_--_3","Lu_DFxDF_--_4",
  "Lu_LFDxLFD_--_1","Lu_LFDxLFD_--_2","Lu_LFDxLFD_--_3",
                                         
  "Sp_DFxDF_++_1","Sp_DFxDF_++_2","Sp_DFxDF_++_3","Sp_XFLxXFL_++_1",
  "Sp_DFxDF_--_1","Sp_DFxDF_--_2","Sp_DFxDF_--_3","Sp_DFxDF_--_4",
  "Sp_LFDxLFD_--_1","Sp_LFDxLFD_--_2","Sp_LFDxLFD_--_3",
  "Sp_LFxLF_--_1","Sp_LFxLF_--_2",
  "Sp_DxD_--_1","Sp_DxD_--_2",
  "Sp_FixFi_--_1","Sp_FixFi_--_2","Sp_FixFi_--_3",
  "Sp_XFLxXFL_--_1","Sp_XFLxXFL_--_2","Sp_XFLxXFL_--_3")


######------ Calculate pearson correlation matrix  ------######
# filter rows with only 0 TPMs
TPM_table <- as.data.frame(TPM_table[apply(TPM_table[,8:83],1,function(x) !all(x==0)),]) 
# filter for autosomes and sex chromosomes (remove MT and other patches)
chr <- c(1:19,"X","Y")
TPM_table <- as.data.frame(TPM_table)%>%filter(chromosome_name %in% chr)


######------ Plot correlation matrix  ------######
breaklist <- seq(0, 1, length.out=161)
col <- colorRampPalette(c("white","#B2182B"))(200)


pdf(file = paste0(output,"pearson_correlation_heatmap-paper-samples.pdf"), pointsize=2)
pheatmap(cor(TPM_table[8:83],method="pearson"),col=col,fontsize = 4, show_colnames = F, breaks=breaklist) 
dev.off()




