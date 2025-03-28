##***************************************##
##   Figure 5 Number of DEGs in spleen   ##
##***************************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 09.2022
## Creation: 08.2021
## Get numbers and direction of degs in spleen
## Cut-off degs: |log2FC| >= 1, FDR <= 0.01
library(ggplot2)
library(tidyverse)
library(data.table)
library(biomaRt)
FDR=0.01
LFC=1


######------ set environment  ------###### 
input <- "./02_RNAseq_DESeq2/"
output <- "./number_degs/"
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="http://aug2020.archive.ensembl.org")


######------ get degs  ------######
samples <- c("XFL","Fi","D","DF","LF","LFD")
for (sample in samples){
  tissue="spleen"
  GT=sample
  # import 
  data <- read.table(paste0(Sys.glob(paste0(input,tissue,"/",GT,"/",GT,"_shrunk_diffexpr-results*"))),
                     header = T)
  data <- data%>%dplyr::filter(abs(log2FoldChange) >= LFC & padj <= FDR)
  # write degs for GT
  write.table(data, paste0(output,"DEG_",GT,"_",tissue,"_FDR",FDR,"_LFC",LFC,".txt"),
              row.names = F,col.names = T,quote = F,sep="\t")
  # get numbers of up and downregulated degs and write out
  data <- data%>%dplyr::select("log2FoldChange")
  up <- sum(data>0); down <- sum(data<0); up_down <- data.frame(up=up,down=down)
  write.table(up_down, paste0(output,"DEG_",GT,"_",tissue,"_FDR",FDR,"_LFC",LFC,"_direction_info.txt"),
              row.names = F,col.names = T,quote = F,sep="\t")  
  assign(paste0(GT), data)
}
