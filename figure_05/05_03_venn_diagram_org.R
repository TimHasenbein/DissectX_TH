##**********************************************************##
##   Figure 5 Venn diagram for shared genes across tissues  ##
##**********************************************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 09.2022
## Creation: 12.2021
## Show DEG that are shared across different tissues in upset plot
## Cut-off degs: |log2FC| >= 1, FDR <= 0.01
library(UpSetR)
library(ggplot2)
library(tidyverse)
library(data.table)
library(biomaRt)
FDR=0.01
LFC=1


######------ set environment  ------###### 
input <- "./"
output <- "./venn_plot/"


######------ get degs  ------######
samples <- c("spleen_LFD","heart_LFD","brain_LFD","lung_LFD","liver_LFD","kidney_LFD")
for (sample in samples){
  split_sample=unlist(strsplit(sample,"_"))
  tissue=split_sample[1]
  if (length(split_sample)==3){
    GT=paste(split_sample[2],split_sample[3],sep="_")
  } else{
    GT=split_sample[2]
  }
  # import 
  data <- read.table(paste0(Sys.glob(paste0(input,tissue,"/",GT,"/",GT,"_shrunk_diffexpr-results*"))),
                     header=T)
  data <- data%>%dplyr::filter(abs(log2FoldChange) >= LFC & padj <= FDR)
  data <- data%>%dplyr::select("name","log2FoldChange")
  colnames(data) <- c("name",paste0(tissue)); data[,2] <- 1
  assign(paste0(tissue), data)
}


######------ make plot dataframe  ------######
plot <- spleen%>%merge(kidney, by=c("name"), all = T)%>%merge(lung,by=c("name"), all = T)%>%
  merge(heart, by=c("name"), all = T)%>%merge(liver,by=c("name"), all = T)%>%merge(brain,by=c("name"), all = T)
plot[is.na(plot)] <- 0
sum(duplicated(plot$name))


######------ Venn plot ------######
# select intersections 
intersections <- list(list("lung","spleen"),
                      list("lung","kidney"),
                      list("lung","heart"),
                      list("spleen","kidney"),
                      list("heart","kidney"),
                      list("spleen","heart"),
                      list("liver","kidney"),
                      list("lung","brain"),
                      list("liver","heart"),
                      list("liver","spleen"),
                      list("lung","kidney","spleen"),
                      list("lung","heart","spleen"),
                      list("heart","kidney","spleen"),
                      list("liver","heart","kidney"),
                      list("liver","lung","kidney"),
                      list("lung","heart","kidney"),
                      list("brain","spleen","kidney"),
                      list("heart","kidney","spleen","lung"),
                      list("heart","kidney","spleen","brain"),
                      list("lung","kidney","spleen","brain"),
                      list("lung","kidney","brain","liver"),
                      list("kidney","spleen","liver","heart"),
                      list("lung","kidney","spleen","liver"), 
                      list("lung","kidney","spleen","brain","heart"),
                      list("lung","kidney","spleen","liver","heart"),
                      list("lung","kidney","spleen","brain","liver","heart"))        
pdf(file=paste0(output,"venn_spleen_FDR",FDR,"_lfc",LFC,"_shrunk.pdf"), paper="a4",width=25, height=5)
upset(plot, point.size = 1.5, line.size = 0.4,
      main.bar.color="black" ,order.by = c("freq","degree"),
      keep.order = T, mainbar.y.label="DE genes shared",matrix.dot.alpha=0,
      nintersects = NA,
      matrix.color="black",intersections = intersections,set_size.show = T,
      sets=c("spleen","kidney","lung","heart","liver","brain")) 
dev.off()


######------ Save shared genes ------######
shared_genes <- plot[rowSums(plot[,2:ncol(plot)])>=2,]
saveRDS(shared_genes,paste0(output,"shared_genes.rds"))
