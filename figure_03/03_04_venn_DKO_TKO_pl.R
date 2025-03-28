##**************************************************##
##   Venn diagrams for Xa/Xi DKO and TKO overlaps   ##
##**************************************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 09.2022
## Creation: 12.2021
## Venn diagram for overlaps in DKO and TKO 
## Cut-off degs: |log2FC| >= 1, FDR <= 0.01
library(tidyverse)
library(data.table)
library(biomaRt)
library(eulerr)
library(gridExtra)
LFC=1
FDR=0.01


######------ Set environment  ------###### 
input <- "./placenta/"
output <- "./venn_diagram/"
mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="http://aug2020.archive.ensembl.org")


######------ Get data  ------######
samples <- c("DFxC","LFxC","LFDxC",
             "CxLF","CxLFD") # "CxDF" has 0 degs
for (sample in samples){
  GT=sample
  # import 
  data <- read.table(paste0(Sys.glob(paste0(input,"/",GT,"/",GT,"_shrunk_diffexpr-results.txt"))))
  data <- data %>% dplyr::filter(abs(log2FoldChange) >= LFC & padj <= FDR)
  data <- setDT(as.data.frame(data), keep.rownames = TRUE)[]
  # get gene names
  data$rn <- sapply(strsplit(data$rn, ".", fixed=T), function(x) x[1])
  gm35612 <- data%>%filter(rn=="Gm35612") 
  IDs <- getBM(attributes = c("ensembl_gene_id","mgi_symbol"),
               filters = "ensembl_gene_id", values = data[,1],
               mart = mart)
  data <- data%>%inner_join(IDs, by=c("rn"="ensembl_gene_id"))
  data$rn <- data$mgi_symbol
  data <- data%>%dplyr::select(-"mgi_symbol")
  data <- data%>%filter(rn!="")
  data <- data%>%rbind(gm35612)
  assign(paste0(GT), data)
}


######------ Venn plot fo Xa ------######
list_venn <- list("LFDxC"=LFDxC$rn, "LFxC"=LFxC$rn, "DFxC"=DFxC$rn)
set.seed(6)
venn <- plot(euler(list_venn, shape = c("ellipse")),input = "disjoint", quantities = TRUE, fills=c("#84C9B0","#156B37","#1F75B9"),edges = "black", alpha=0.4, lwd = 0.5)
pdf(file=paste(output,"venn_pl_Xa_FDR0.01_LFC1_shrunk.pdf",sep=""), paper="a4r",width=8, height=8)
grid.arrange(venn, ncol = 1)
dev.off()
DFxC$DF <- 1 
LFxC$LF <- 1 
LFDxC$LFD <- 1 
overlap <- DFxC%>%merge(LFDxC,by=c("rn"), all=T)%>%merge(LFxC,by=c("rn"), all=T)
overlap_xa = overlap%>%dplyr::select("rn","DF","LF","LFD")
write.table(overlap_xa, paste0(output,"overlap_xa.txt"), row.names = F, col.names = T, sep = "\t", quote = F)


######------ Venn plot fo Xi ------######
list_venn <- list("CxLFD"=CxLFD$rn, "CxLF"=CxLF$rn)
set.seed(2)
venn <- plot(euler(list_venn, shape = c("ellipse")),input = "disjoint", quantities = TRUE, fills=c("#84C9B0","#156B37"),edges = "black", alpha=0.4, lwd = 0.5)
pdf(file=paste(output,"venn_pl_Xi_FDR0.01_LFC1_shrunk.pdf",sep=""), paper="a4r",width=8, height=8)
grid.arrange(venn, ncol = 1)
dev.off()




