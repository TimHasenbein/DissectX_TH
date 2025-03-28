##**********************************************##
##   Figure 5 Heatmap for LFD/LF degs overlap   ##
##**********************************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 12.2022
## Creation: 12.2022
## Heatmap for overlapping degs for LFD/LF as in venn diagram
## Cut-off degs: |log2FC| >= 1, FDR <= 0.01
library(ggplot2)
library(tidyverse)
library(data.table)
library(ComplexHeatmap)


######------ Set environment  ------###### 
input <- "./02_RNAseq_DESeq2/"
output <- "./heatmap_LFC/"


######------ get LFC data for heatmap ------###### 
samples <- c("spleen","heart","brain","lung","liver","kidney")
for (sample in samples){
  tissue=sample
  GT="LFD"
  # import 
  data <- read.table(paste0(Sys.glob(paste0(input,tissue,"/",GT,"/",GT,"_shrunk_diffexpr-results*"))),
                     header=T)
  data <- data%>%filter(padj <= 0.01 & abs(log2FoldChange) >= 1)
  data <- data%>%dplyr::select("name","log2FoldChange")
  colnames(data) <- c("name",paste0(tissue))
  assign(paste0(tissue), data)
}


######------ Merge ------######
LFC <- spleen%>%merge(kidney, by="name", all = T)%>%merge(lung,by=c("name"), all = T)%>%
  merge(heart, by=c("name"), all = T)%>%merge(liver,by=c("name"), all = T)%>%merge(brain,by=c("name"), all = T)


######------ Get index matrix ------###### 
LFC[is.na(LFC)] <- 0
LFC_index <- sign(LFC[,2:7])
LFC_index$name <- LFC$name
LFC_index$rowsum <- rowSums(abs(LFC_index[,1:6]))
LFC_index <- LFC_index%>%filter(rowsum >= 2 )
LFC_index$direction <- ifelse(rowSums(LFC_index[,1:6] == 1) >= rowSums(LFC_index[,1:6] == -1),
                              "up","down")
# filter for number shared
LFC_2_up <- LFC_index%>%filter(rowsum==2 & direction=="up")
LFC_2_down <- LFC_index%>%filter(rowsum==2 & direction=="down")
LFC_3_up <- LFC_index%>%filter(rowsum==3 & direction=="up")
LFC_3_down <- LFC_index%>%filter(rowsum==3 & direction=="down")
LFC_4_up <- LFC_index%>%filter(rowsum==4 & direction=="up")
LFC_4_down <- LFC_index%>%filter(rowsum==4 & direction=="down")
LFC_5_up <- LFC_index%>%filter(rowsum==5 & direction=="up")
LFC_5_down <- LFC_index%>%filter(rowsum==5 & direction=="down")
LFC_6_up <- LFC_index%>%filter(rowsum==6 & direction=="up")
LFC_6_down <- LFC_index%>%filter(rowsum==6 & direction=="down")
# order 
LFC_2_up_sort <- LFC_2_up[order(-LFC_2_up$spleen, -LFC_2_up$kidney,-LFC_2_up$lung,-LFC_2_up$heart,-LFC_2_up$liver,-LFC_2_up$brain),]
LFC_2_down_sort <- LFC_2_down[order(LFC_2_down$spleen, LFC_2_down$kidney,LFC_2_down$lung,LFC_2_down$heart,LFC_2_down$liver,LFC_2_down$brain),]
# 3
LFC_3_up_sort <- LFC_3_up[order(-LFC_3_up$spleen, -LFC_3_up$kidney,-LFC_3_up$lung,-LFC_3_up$heart,-LFC_3_up$liver,-LFC_3_up$brain),]
LFC_3_down_sort <- LFC_3_down[order(-LFC_3_down$spleen, -LFC_3_down$kidney,-LFC_3_down$lung,-LFC_3_down$heart,-LFC_3_down$liver,-LFC_3_down$brain),]
# 4
LFC_4_up_sort <- LFC_4_up[order(-LFC_4_up$spleen, -LFC_4_up$kidney,-LFC_4_up$lung,-LFC_4_up$heart,-LFC_4_up$liver,-LFC_4_up$brain),]
LFC_4_down_sort <- LFC_4_down[order(-LFC_4_down$spleen, -LFC_4_down$kidney,-LFC_4_down$lung,-LFC_4_down$heart,-LFC_4_down$liver,-LFC_4_down$brain),]
# 5
LFC_5_up_sort <- LFC_5_up[order(-LFC_5_up$spleen, -LFC_5_up$kidney,-LFC_5_up$lung,-LFC_5_up$heart,-LFC_5_up$liver,-LFC_5_up$brain),]
LFC_5_down_sort <- LFC_5_down[order(-LFC_5_down$spleen, -LFC_5_down$kidney,-LFC_5_down$lung,-LFC_5_down$heart,-LFC_5_down$liver,-LFC_5_down$brain),]
# 6
LFC_6_up_sort <- LFC_6_up[order(-LFC_6_up$spleen, -LFC_6_up$kidney,-LFC_6_up$lung,-LFC_6_up$heart,-LFC_6_up$liver,-LFC_6_up$brain),]
LFC_6_down_sort <- LFC_6_down[order(-LFC_6_down$spleen, -LFC_6_down$kidney,-LFC_6_down$lung,-LFC_6_down$heart,-LFC_6_down$liver,-LFC_6_down$brain),]

all_index <- LFC_6_up_sort%>%rbind(LFC_6_down_sort)%>%rbind(LFC_5_up_sort)%>%rbind(LFC_5_down_sort)%>%rbind(LFC_4_up_sort)%>%rbind(LFC_4_down_sort)%>%rbind(LFC_3_up_sort)%>%rbind(LFC_3_down_sort)%>%rbind(LFC_2_up_sort)%>%rbind(LFC_2_down_sort)


######------ get LFC data for heatmap ------###### 
samples <- c("spleen","heart","brain","lung","liver","kidney")
for (sample in samples){
  sample="spleen"
  tissue=sample
  GT="LFD"
  # import 
  data <- read.table(paste0(Sys.glob(paste0(input,tissue,"/",GT,"/",GT,"_shrunk_diffexpr-results*"))),
                     header=T)
  data <- data%>%filter(padj <= 0.01)
  # modify table for downstream analysis
  data <- data%>%dplyr::select("name","log2FoldChange")
  colnames(data) <- c("name",paste0(tissue))
  assign(paste0(tissue), data)
}
LFC <- spleen%>%merge(kidney, by=c("name"), all = T)%>%merge(lung,by=c("name"), all = T)%>%
  merge(heart, by=c("name"), all = T)%>%merge(liver,by=c("name"), all = T)%>%merge(brain,by=c("name"), all = T)


######------ Add info about up/down regulated  ------######
LFC_mat <- LFC %>% dplyr::select(-name) %>% as.matrix()
rownames(LFC_mat) <- LFC$name
LFC_mat[is.na(LFC_mat)] <- 0
plot <- LFC_mat[all_index$name,]


######------ Order plot again  ------######
plot <- as.data.frame(plot)


######------ Adjust LFC scale  ------######
for(r in 1:nrow(plot)){
  for(c in 1:ncol(plot)){
    ifelse(plot[r,c] >= 2, 
           plot[r,c] <- 2,
           ifelse(plot[r,c] <= -2, 
                  plot[r,c] <- -2,plot[r,c] <- plot[r,c]))
  }
}


######------ heatmap ------######
colfunc_heat <- colorRampPalette(c("black", "white","#C04000"))(100)
cols <- c("#000000", "#484848","white", "white","#D27648","#C04000")
pdf(file=paste0(output,"LFC_heatmap_shared2.pdf"), paper="a4",width=10, height=20) 
Heatmap(as.matrix(plot), name = "logpadj", cluster_rows = F, cluster_columns = F,
        column_title = NULL, na_col = "white", show_row_names = T, row_names_gp =  gpar(fontsize = 4),
        col = cols, 
        width = unit(10, "cm"), border = T,
        row_gap = unit(0.5, "mm"), 
        row_title_gp = gpar(fontsize = 6, col="black"),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(at = c(-2,-1,-0.5,0,0.5,1,2)))
dev.off()



                      
                      
                      
                      
                      