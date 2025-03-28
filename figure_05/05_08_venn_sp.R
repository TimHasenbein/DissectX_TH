##**********************************##
##   Figure 5 Venn diagram spleen   ##
##**********************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 09.2022
## Creation: 08.2021
## Cut-off degs: |log2FC| >= 1, FDR <= 0.01
library(UpSetR)
library(ggplot2)
library(tidyverse)
library(data.table)
library(eulerr)
library(gridExtra)
FDR=0.01
LFC=1


######------ set environment  ------###### 
input <- "./02_RNAseq_DESeq2/"
output <- "./venn_plot/"


######------ get degs  ------######
samples <- c("XFL","Fi","D","DF","LF","LFD")
for (sample in samples){
  tissue="spleen"
  GT=sample
  # import 
  data <- read.table(paste0(Sys.glob(paste0(input,tissue,"/",GT,"/",GT,"_shrunk_diffexpr-results*"))),
                     header=T)
  data <- data%>%dplyr::filter(abs(log2FoldChange) >= LFC & padj <= FDR)
  data <- data%>%dplyr::select("name","log2FoldChange")
  colnames(data) <- c("name",paste0(GT)); data[,2] <- 1
  assign(paste0(GT), data)
}


######------ make plot dataframe  ------######
plot <- XFL%>%merge(Fi, by=c("name"), all = T)%>%merge(D,by=c("name"), all = T)%>%
  merge(DF, by=c("name"), all = T)%>%merge(LF,by=c("name"), all = T)%>%merge(LFD,by=c("name"), all = T)
plot[is.na(plot)] <- 0
sum(duplicated(plot$name))


######------ export results  ------######
write.table(plot, paste0(output,"overlapping_degs.txt"), col.names = T, row.names = F, sep = "\t", quote = F)


######------ Venn plot ------######
# select intersections
intersections <- list(list("LFD"),list("LF"),list("DF"),list("D"),list("XFL"),list("Fi"),
                       list("Fi","D"),list("DF","XFL"),list("DF","LF"),list("Fi","LFD"),
                       list("D","LFD"),list("DF","LFD"),list("XFL","LFD"),list("LF","LFD"),
                       list("DF", "XFL","LF"), list("LFD","DF","D"), # Dxz4
                       list("Fi","LFD","LF","DF")) # Firre
# plot
pdf(file=paste0(output,"venn_spleen_FDR",FDR,"_lfc",LFC,"_shrunk.pdf"), paper="a4",width=6, height=4)
  upset(plot, point.size = 1.5, line.size = 0.4,  
        main.bar.color="black" ,order.by = c("degree"),decreasing = F,
        keep.order = T, mainbar.y.label="DE genes shared",matrix.dot.alpha=0,
        sets=c("LFD","LF","XFL","DF","D","Fi"),
        sets.bar.color= c("#93CDB9","#12672F","#3C73B9","#808181","#B81B1A","#DB5D14"), 
        nintersects = NA, intersections = intersections,
        matrix.color="black",set_size.show=T) 
dev.off()