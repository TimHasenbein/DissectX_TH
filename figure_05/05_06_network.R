##******************************##
##  Figure 5 GSEA network plot  ##
##******************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 01.2023
## Creation: 07.2021
## Make GSEA network plots for gene sets of top GS performed for all tissues LFD samples
library(tidyverse)
library(ggplot2)
library(data.table)
library(igraph)
library(simplifyEnrichment)
library(RColorBrewer)
sample = "spleen"


######------ Set environment  ------###### 
input <- "./heatmap_GS/"
output <- "./network_plot/"
topGS <- read.table(paste0(input,"top_100_genesets_cluster.txt"),sep="\t", header=T)
colors <- readRDS(paste0(input,"cluster_colors.rds"))
colors$rn <- as.integer(colors$rn) 


######------ Gene set file  ------###### 
GO_gmt <- "./c5.go.v7.4.symbols.gmt"
GO_gmt <- read.table(GO_gmt, sep="\n", header=F, stringsAsFactors = F)
GO_gmt <- lapply(GO_gmt[,1], function(x) {strsplit(x, split="\t")[[1]]})
names(GO_gmt) <- sapply(GO_gmt, function(x) {x[1]})
GO_gmt <- lapply(GO_gmt, function(x) {x[-c(1,2)]})


######------ Calculate similarity matrix for topGS  ------###### 
topGS_sample <- topGS[c("rn",paste0(sample),"cluster")]
sim_matrix <- term_similarity(GO_gmt[topGS_sample$rn])
sim_matrix[sim_matrix<0.2] <- 0


######------ Calculate network plots ------###### 
# generate igraph object from similarity matrix
network <- graph_from_adjacency_matrix(sim_matrix, weighted=T, mode="undirected", diag=FALSE)
# simplify annotation names
V(network)$name <- tolower(V(network)$name)
V(network)$name <- gsub("gomf_|gocc_|gobp_","",V(network)$name)
V(network)$name <- gsub("_"," ",V(network)$name)
# add pval to make circel size according to plot 
V(network)$pval <- log(abs(topGS_sample[,2])+1)*2


######------ Add colors to network object matching by cluster id  ------###### 
# add cluster col to network
V(network)$cluster <- topGS_sample$cluster
# match colors to clusters
cluster_df <- data.frame(cluster=as.integer(V(network)$cluster))
network_color <- left_join(cluster_df, colors, by = c("cluster"="rn")) # keep same order
# add color col to network
V(network)$color <- network_color$group

 
######------ plot network ------###### 
pdf(file=paste0(output,sample,"_network_plot.pdf"), 
    paper="a4r",width=8, height=15)
set.seed(1)
plot(network, vertex.size=V(network)$pval, vertex.label.cex = 0.4, 
     vertex.label.color="black", vertex.label.dist=0, edge.color="grey",
     vertex.color=V(network)$color,vertex.label = "",edge.lty=1, xlim=c(-1,1),ylim=c(-1,1),axes=F,edge.width=0.05)
# labeled
set.seed(1)
plot(network, vertex.size=V(network)$pval, vertex.label.cex = 0.4, 
     vertex.label.color="black", vertex.label.dist=0, edge.color="grey",
     vertex.color=V(network)$col,edge.lty=1, xlim=c(-1,1),ylim=c(-1,1),axes=F,edge.width=0.05)
dev.off()





