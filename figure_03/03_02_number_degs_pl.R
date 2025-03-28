##**************************************##
##   Figure Number of DEGs placenta   ##
##**************************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 12.2021
## Creation: 12.2021
## Get number of significant DEGs in placenta
## Cut-off degs: |log2FC| >= 1, FDR <= 0.01
library(ggplot2)
library(tidyverse)
library(data.table)
library(biomaRt)
FDR=0.01
LFC=1


######------ Set environment  ------###### 
input <- paste0("./02_RNAseq_DESeq2/placenta/")
output <- paste0("./number_degs/")
mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="http://aug2020.archive.ensembl.org")


######------ get degs  ------######
samples <- c("XFLxC","FixC","DFxC","LFxC","LFDxC", 
             "CxXFL","CxFi","CxD","CxLF","CxLFD") # "DxC" and  "CxDF" have 0
sample <- "DxC"
# run loop
for (sample in samples){
  GT=sample
  # import 
  data <- read.table(paste0(Sys.glob(paste0(input,GT,"/",GT,"_shrunk_diffexpr-results*"))))
  data <- data%>%dplyr::filter(abs(log2FoldChange) >= LFC & padj <= FDR)
  # change name and add position
  data <- setDT(data, keep.rownames = TRUE)[]
  data$rn <- sapply(strsplit(data$rn, ".", fixed=T), function(x) x[1])
  gm35612 <- data%>%dplyr::filter(rn=="Gm35612")
  IDs <- getBM(attributes = c("ensembl_gene_id","mgi_symbol"),
               filters = "ensembl_gene_id", values = data[,1],
               mart = mart)
  data <- data%>%inner_join(IDs, by=c("rn"="ensembl_gene_id"))
  data <- data%>%dplyr::filter(mgi_symbol != "")
  data$rn <- data$mgi_symbol
  data <- data%>%dplyr::select(-"mgi_symbol")
  data <- data%>%rbind(gm35612)
  # write degs for GT
  write.table(data, paste0(output,"DEG_",GT,"_FDR",FDR,"_LFC",LFC,".txt"),
             row.names = T,col.names = T,quote = F,sep="\t")
  # get numbers of up and downregulated degs and write out
  data <- data%>%dplyr::select("log2FoldChange")
  up <- sum(data>0); down <- sum(data<0); up_down <- data.frame(up=up,down=down)
  write.table(up_down, paste0(output,"DEG_",GT,"_FDR",FDR,"_LFC",LFC,"_direction.txt"),
              row.names = T,col.names = T,quote = F,sep="\t")  
  assign(paste0(GT), data)
}

