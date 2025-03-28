##******************************##
##  Figure 5 GSEA heatmap plot  ##
##******************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 09.2022
## Creation: 09.2021
## Visualize GSEA results per tissue in a heatmap to identify common features
library(tidyverse)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(igraph)
library(simplifyEnrichment)
gt = "LFD"
top_rank = 100


######------ Set environment  ------###### 
input <- "./04_RNAseq_GSEA/"
output <- "./03_dissectX_RNA/"


######------ Gene set file  ------###### 
GO_gmt <- paste0(input,"c5.go.v7.4.symbols.gmt")
GO_gmt <- read.table(GO_gmt, sep="\n", header=F, stringsAsFactors = F)
GO_gmt <- lapply(GO_gmt[,1], function(x) {strsplit(x, split="\t")[[1]]})
names(GO_gmt) <- sapply(GO_gmt, function(x) {x[1]})
GO_gmt <- lapply(GO_gmt, function(x) {x[-c(1,2)]})


######------ get sig GS for organs ------###### 
samples <- c("spleen","heart","kidney","lung","liver","brain")
for (sample in samples){
  data <- read.table(paste0(input,sample,"/LFD/",sample,"_",gt,"_GSEA_FDR0.1_min10max500.txt"))
  data <- data%>%filter(padj<=0.1)
  data <- data%>%dplyr::select("rn", "log_padj")
  colnames(data) <- c("gene_set",paste0(sample))
  assign(paste(sample),data)
} 
sigGS <- spleen%>%merge(heart, by=c("gene_set"), all = T)%>%merge(kidney, by=c("gene_set"), all = T)%>%
  merge(lung,by=c("gene_set"), all = T)%>%merge(liver,by=c("gene_set"), all = T)%>%merge(brain,by=c("gene_set"), all = T) 
rm(spleen,heart,kidney,lung,liver,brain)


######------ get all results of GSEA for spleen GTs ------###### 
for(sample in samples){
  data <- read.table(paste0(input,sample,"/LFD/",sample,"_",gt,"_GSEA.txt"))
  data <- setDT(as.data.frame(data), keep.rownames = TRUE)[]
  data <- data%>%dplyr::select("rn", "log_padj")
  colnames(data) <- c("gene_set",paste0(sample))
  assign(paste0(sample), data)
} 
GS_res <- spleen%>%merge(heart, by=c("gene_set"), all = T)%>%merge(kidney, by=c("gene_set"), all = T)%>%
  merge(lung,by=c("gene_set"), all = T)%>%merge(liver,by=c("gene_set"), all = T)%>%merge(brain,by=c("gene_set"), all = T)


######------ filter results sig in at least one GT  ------###### 
GS_res_filter <- GS_res%>%filter(gene_set %in% sigGS$gene_set) # 1187 GS


######------ filter top 100 GS based on log_padj  ------###### 
GS_res_filter[is.na(GS_res_filter)] <- 0
GS_res_filter$max <- apply(abs(GS_res_filter[,2:7]), 1, max)
GS_res_filter <- GS_res_filter%>%arrange(desc(max))
GS_res_filter <- GS_res_filter%>%column_to_rownames("gene_set") 
GS_res_top <- GS_res_filter[1:100,1:6]


######------ Calculate similarity matrix for topGS  ------###### 
sim_matrix <- term_similarity(GO_gmt[rownames(GS_res_top)])
# filter for edges (edges is drawn when similarity (# of genes shared btw 2 GOs) is >0)
sim_matrix[sim_matrix<0.2] <- 0


######------ Calculate cluster for topGS  ------###### 
cluster <- cluster_terms(sim_matrix, method = "walktrap")
# get a table to annotate the clusters
topGS_cluster <- data.frame(rownames(GS_res_top), cluster)
topGS_cluster <- topGS_cluster[order(cluster),] 


######------ assign clusters to heatmap matrix  ------###### 
GS_res_top <- setDT(as.data.frame(GS_res_top), keep.rownames = TRUE)[]
GS_res_top <- GS_res_top%>%left_join(topGS_cluster, by = c("rn"="rownames.GS_res_top."))
write.table(GS_res_top,paste0(output,"top_",top_rank,"_genesets.txt"),sep="\t",quote=F)


######------ row annotation of clusters ------###### 
# assign annotation color to ra
annotation <- GS_res_top%>%dplyr::select("cluster")
anno_col <- unique(annotation)
col <- colorRampPalette(brewer.pal(9, "Spectral"))(length(anno_col$cluster)) 
names(col) <- anno_col$cluster
col <- list(group = col)
# save colors
col_df <- setDT(as.data.frame(col), keep.rownames = TRUE)[]
saveRDS(col_df,paste0(output,"cluster_colors.rds"))
# get annotation bar
row_annotation = rowAnnotation(group = annotation$cluster, 
                               col = col, border = T,show_legend=T)


######------ plot heatmap ------######
colfunc_heat <- circlize::colorRamp2(c(-10, 0, 10), c("black","white","#C04000")) 
plot <- GS_res_top%>%column_to_rownames("rn")
# scale
logpadj = 10
for(r in 1:nrow(plot)){
  for(c in 1:6){
    ifelse(plot[r,c] >= logpadj, 
           plot[r,c] <- logpadj,
           ifelse(plot[r,c] <= -logpadj, 
                  plot[r,c] <- -logpadj,
                  plot[r,c] <- plot[r,c]))
  }
}
pdf(file=paste0(output,"/",gt,"_heatmap_top",top_rank,"_GS_with_nonsig.pdf"), 
    paper="a4",width=10, height=15)
ht <- Heatmap(as.matrix(plot[,1:6]), name = "logpadj", 
              cluster_rows = T, clustering_distance_rows = "euclidean",
              clustering_method_rows = "complete",
              cluster_columns = T, clustering_distance_columns = "euclidean",
              clustering_method_columns = "complete",
              column_title = NULL, na_col = "white", show_row_names = T, 
              left_annotation = row_annotation, 
              row_names_gp =  gpar(fontsize = 4),
              col = colfunc_heat, 
              width = unit(10, "cm"), border = T,
              row_gap = unit(0.5, "mm"), 
              row_title_gp = gpar(fontsize = 6, col="black"),
              column_names_gp = gpar(fontsize = 10),
              heatmap_legend_param = list(at = c(-10,-5,0,5,10))) 
draw(ht, heatmap_legend_side = "left")
dev.off()
write.table(GS_res_top,paste0(output,"top_",top_rank,"_genesets_cluster.txt"),sep="\t",quote=F,row.names = F)