##*************************************##
##  Figure 5 GSEA heatmap plot spleen  ##
##*************************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 09.2022
## Creation: 09.2021
## Visualize GSEA results for the different KO models in spleen
## Heatmap for top 100 GS with FDR>0.1 and min10 max 500 genes
library(tidyverse)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(writexl)
tissue <- "spleen"


######------ Set environment  ------###### 
input <- "./04_RNAseq_GSEA/"
output <- "./GSEA_heatmap/"


######------ Gene set file  ------###### 
GO_gmt <- paste0(input,"c5.go.v7.4.symbols.gmt")
GO_gmt <- read.table(GO_gmt, sep="\n", header=F, stringsAsFactors = F)
GO_gmt <- lapply(GO_gmt[,1], function(x) {strsplit(x, split="\t")[[1]]})
names(GO_gmt) <- sapply(GO_gmt, function(x) {x[1]})
GO_gmt <- lapply(GO_gmt, function(x) {x[-c(1,2)]})


######------ get all results of GSEA for spleen GTs ------######
GTs = c("XFL","Fi","D","LF","LFD") 
for(GT in GTs){
  data <- read.table(paste0(input,tissue,"/",GT,"/",tissue,"_",GT,"_GSEA.txt"), sep="\t")
  data <- setDT(as.data.frame(data), keep.rownames = TRUE)[]
  data <- data%>%dplyr::select("rn", "log_padj")
  colnames(data) <- c("gene_set",paste0(GT))
  assign(paste0(GT), data)
} 
GS_res <- LFD%>%merge(LF, by=c("gene_set"), all = T)%>%merge(D,by=c("gene_set"), all = T)%>%merge(Fi,by=c("gene_set"), all = T)%>%merge(XFL,by=c("gene_set"), all = T) 
rm(XFL,Fi,D,LF,LFD)


######------ get organ GS top 100 for plotting ------###### 
gs <- read.table("./top_100_genesets.txt")["rn"]


######------ filter gene-sets from data ------###### 
GS_res_filter <- GS_res%>%filter(gene_set %in% gs$rn) # 100
# supplementary output
suppl <- setDT(GS_res_filter, keep.rownames = TRUE)[]
write_xlsx(GS_res_filter, paste0(output,"GSEA_spleen.xlsx"), col_names = TRUE)
GS_res_filter[is.na(GS_res_filter)] <- 0
GS_res_filter <- GS_res_filter%>%column_to_rownames("gene_set")


######------ scale ------###### 
logpadj = 10
for(r in 1:nrow(GS_res_filter)){
  for(c in 1:5){
    ifelse(GS_res_filter[r,c] >= logpadj, 
           GS_res_filter[r,c] <- logpadj,
           ifelse(GS_res_filter[r,c] <= -logpadj, 
                  GS_res_filter[r,c] <- -logpadj,
                  GS_res_filter[r,c] <- GS_res_filter[r,c]))
  }
}


######------ specify clustering ------###### 
lf_lfd <- GS_res_filter[1:2]
dist_matrix <- dist(lf_lfd, method = "euclidean")
cluster_result <- hclust(dist_matrix, method="average")
reordered_matrix <- GS_res_filter[rev(cluster_result$order), ]


######------ plot heatmap ------######
colfunc_heat <- circlize::colorRamp2(c(-10, 0, 10), c("black","white","#C04000")) 
pdf(file=paste0(output,"/heatmap_tissue_pheno_embo.pdf"), paper="a4",width=10, height=15)
ht <- Heatmap(as.matrix(reordered_matrix), name = "logpadj", 
              cluster_rows = F,
              cluster_columns = F, clustering_distance_columns = "euclidean",
              clustering_method_columns = "complete",
              column_title = NULL, na_col = "white", show_row_names = T, 
              row_names_gp =  gpar(fontsize = 4),
              col = colfunc_heat, 
              width = unit(10, "cm"), border = T,
              row_gap = unit(0.5, "mm"), 
              row_title_gp = gpar(fontsize = 6, col="black"),
              column_names_gp = gpar(fontsize = 10),
              heatmap_legend_param = list(at = c(-10,-5,0,5,10))) 
draw(ht, heatmap_legend_side = "left")
dev.off()
