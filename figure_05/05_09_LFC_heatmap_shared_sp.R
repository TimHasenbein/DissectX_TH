##******************************************************##
##   Figure 5 Heatmap for LFD/LF females degs overlap   ##
##******************************************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 09.2022
## Creation: 12.2021
## Heatmap for overlapping degs for LFD/LF as in venn diagram
## Cut-off degs: |log2FC| >= 1, FDR <= 0.01
library(ggplot2)
library(tidyverse)
library(data.table)
library(biomaRt)
library(pheatmap)
LFC=1
FDR=0.01


######------ Set environment  ------###### 
input <- "./02_RNAseq_DESeq2/"
output <- "/LFC_heatmap/"
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="http://aug2020.archive.ensembl.org")


######------ Get data  ------######
samples <- c("XFL","Fi","D","DF","LF","LFD")
for (sample in samples){
  tissue="spleen"
  GT=sample
  data <- read.table(paste0(Sys.glob(paste0(input,tissue,"/",GT,"/",GT,"_shrunk_diffexpr-results*"))),
                     header=T)
  data <- data %>% dplyr::filter(abs(log2FoldChange) >= LFC & padj <= FDR)%>%dplyr::select("name","padj")
  colnames(data) <- c("name",paste0(GT)); data[,2] <- 1
  assign(paste0(GT), data)
}


######------ get gene name of shared degs between LF/LFD females ------######
degs <- LFD%>%merge(LF, by=c("name"), all = T)%>%merge(XFL,by=c("name"), all = T)%>%
  merge(DF, by=c("name"), all = T)%>%merge(D,by=c("name"), all = T)%>%merge(Fi,by=c("name"), all = T)
degs_shared <- degs%>%filter(LFD==1 & LF==1)
degs_shared <- degs_shared$name


######------ get LFC data for all shared samples ------###### 
for (sample in samples){
  GT=sample
  data <- read.table(paste0(Sys.glob(paste0(input,tissue,"/",GT,"/",GT,"_shrunk_diffexpr-results*"))),
                     header=T)
  data <- data%>%dplyr::select("name","log2FoldChange")
  colnames(data) <- c("name",paste0(GT))
  assign(paste0(GT), data)
}
# merge
LFC <- XFL%>%merge(Fi, by=c("name"), all = T)%>%merge(D,by=c("name"), all = T)%>%
  merge(DF, by=c("name"), all = T)%>%merge(LF,by=c("name"), all = T)%>%merge(LFD,by=c("name"), all = T)


######------ Get values for shared degs (filter)  ------######
LFC_shared <- LFC%>%filter(name %in% degs_shared) 


######------ Filter LF male/female and LFD columns  ------######
LFC_shared <- LFC_shared[,c("name","LFD","LF")] # 74 genes
write.table(LFC_shared, paste0(output,"LFC_heatmap_lfc_values.txt"), col.names = T, row.names = T, sep = "\t", quote = F)


######------ Generate plot matrix  ------######
plot <- LFC_shared %>% column_to_rownames("name") %>% as.matrix()


######------ Adjust LFC scale  ------######
for(r in 1:nrow(plot)){
  for(c in 1:ncol(plot)){
    ifelse(plot[r,c] >= 2, 
           plot[r,c] <- 2,
           ifelse(plot[r,c] <= -2, 
                  plot[r,c] <- -2,
                  plot[r,c] <- plot[r,c]))
  }
}


######------ heatmap ------######
cols = colfunc_heat <- colorRampPalette(c("black", "white","#C04000"))(100)
pdf(file=paste0(output,"LFC_heatmap_sharedLFDLF_females_shrunk.pdf"), paper="a4",width=5, height=20) 
pheatmap(plot,fontsize_row = 4,color=cols, cluster_cols = F,legend_breaks=c(-2,-1,0,1,2), show_rownames=T) 
dev.off()





