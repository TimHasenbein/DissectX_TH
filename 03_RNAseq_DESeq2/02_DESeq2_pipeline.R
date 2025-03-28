#**********************************************##
## DESeq2 Differential Gene Expression Analysis ##
##**********************************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 09.2022
## Creation: 06.2021
## Performs differential gene-expression analysis using DESeq2 per tissue individually
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(biomaRt)
library(data.table)
library(fdrtool)
library(EnhancedVolcano)


######------ Sample info  ------###### 
sample <- "LFD"
tissue <- "spleen"
abbr <- "Sp"


######------ Set environment  ------###### 
output <-  paste0("./",tissue,"/",sample,"/")
count_dir <- "./"
mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="http://aug2020.archive.ensembl.org")

######-------------------------------------------------------------------######
##                        Differential gene expression                       ##
######-------------------------------------------------------------------######

######------  Get input data  ------###### 
samples <- c(paste0(abbr,"_DFxDF_\\++"),paste0(abbr,"_",sample,"x",sample,"_\\--"),paste0(abbr,"_XFLxXFL_\\++_1"))
sampleFiles <- grep(paste(samples, collapse = "|"),list.files(count_dir),value=TRUE)
sampleNames <- sub("*_unstranded.count","\\1",sampleFiles)
sampleCondition <- ifelse(str_detect(sampleFiles, "DFxDF_\\++|XFLxXFL_\\++_1"), "WT",paste0(sample))
sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleCondition)
sampleTable$condition <- factor(sampleTable$condition)
write.table(sampleTable,paste0(output,"sample_table.txt"),sep = "\t",row.names =F, quote = F)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = count_dir,
                                       design= ~ condition)


######------  Differential expression analysis  ------###### 
# pre-filtering counts to reduce memory: remove rows with reads <= 1
keep <- rowSums(counts(ddsHTSeq)) >= 1 
dds <- ddsHTSeq[keep,]
# specifying factors so DESeq2 knows which levels to compare
dds$condition <- factor(dds$condition, levels = c("WT", paste0(sample)))
# perform differential gene expression
dds <- DESeq(dds)


######-------------------------------------------------------------------######
##                              Export results                               ##
######-------------------------------------------------------------------######


######------  Raw results  ------###### 
res_raw <- results(dds, format="DataFrame") 
res_raw <- res_raw[order(res_raw$padj),]
# remove genes filtered out by independent filtering and the dispersion outlier
res_raw <- res_raw[ !is.na(res_raw$padj), ]
res_raw <- setDT(as.data.frame(res_raw), keep.rownames = TRUE)[]; res_raw <- res_raw%>%dplyr::rename("name"="rn")
res_raw$name <- sapply(strsplit(res_raw$name, ".", fixed=T), function(x) x[1])
gm35612_raw <- res_raw%>%filter(name=="Gm35612")
IDs_raw <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
             filters = "ensembl_gene_id", values = res_raw[,1],
             mart = mart)
res_raw <- res_raw%>%inner_join(IDs_raw, by=c("name"="ensembl_gene_id"))
res_raw$name <- res_raw$mgi_symbol
res_raw <- res_raw[,1:7]
res_raw <- res_raw%>%filter(name!="")
res_raw <- res_raw%>%rbind(gm35612_raw)
write.table(res_raw, file=paste0(output,sample,"_diffexpr-results.txt"),sep = "\t",row.names =F, quote = F)


######------  Shrunk results  ------###### 
res_shrunk <- lfcShrink(dds, coef=paste0("condition_",sample,"_vs_WT"), type="apeglm")
res_shrunk <- res_shrunk[order(res_shrunk$pvalue),]
res_shrunk <- res_shrunk[ !is.na(res_shrunk$padj), ]
res_shrunk <- setDT(as.data.frame(res_shrunk), keep.rownames = TRUE)[]; res_shrunk <- res_shrunk%>%dplyr::rename("name"="rn")
res_shrunk$name <- sapply(strsplit(res_shrunk$name, ".", fixed=T), function(x) x[1])
gm35612_shrunk <- res_shrunk%>%filter(name=="Gm35612")
IDs <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
             filters = "ensembl_gene_id", values = res_shrunk[,1],
             mart = mart)
res_shrunk <- res_shrunk%>%inner_join(IDs, by=c("name"="ensembl_gene_id"))
res_shrunk$name <- res_shrunk$mgi_symbol
res_shrunk <- res_shrunk[,1:6]
res_shrunk <- res_shrunk%>%filter(name!="")
res_shrunk <- res_shrunk%>%rbind(gm35612_shrunk)
write.table(res_shrunk, file=paste0(output,sample,"_shrunk_diffexpr-results.txt"),sep = "\t",row.names =F, quote = F)


######------  Count table  ------###### 
write.table(as.data.frame(counts(dds, normalized=TRUE)), paste0(output,sample,"_count_table_normalizedTRUE.txt"))


######-------------------------------------------------------------------######
##                          Quality assessment                               ##
######-------------------------------------------------------------------######


######------  Plot dispersion  ------###### 
pdf(file = paste0(output,"qc-dispersions.pdf"), pointsize=2)
plotDispEsts(dds, main="Dispersion plot")
dev.off()


######------  Data quality by clustering  ------######
# Regularized log transformation for clustering/heatmaps; for variance stabilizing effects
rld <- rlogTransformation(dds)
head(assay(rld))
# heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf(file = paste0(output,"qc-heatmap-samples.pdf"), pointsize=2)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, main = "sample-to-sample distances", treeheight_row = 12, treeheight_col = 12, fontsize = 6)
dev.off()
# pca
pcaData <- plotPCA(rld, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf(file = paste0(output,"qc-pca-samples.pdf") , pointsize=2)
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + 
  ggtitle("Principal Component Analysis") +
  geom_text(aes(label=row.names(pcaData), vjust=2, hjust=0.5), colour = "black", size = 2) +
  theme_minimal()
dev.off()


######------  Check pvalue distribution of raw results  ------######
# use z–scores returned by DESeq2 as input to fdrtool
res_corrected <- fdrtool(res_raw$stat, statistic= "normal", plot = F)
sigma <- round(as.numeric(res_corrected$param[1, "sd"]), digits = 2)
# plot variance
pdf(file=paste0(output,"pvalue_hist_after_correction.pdf"), paper="a4r")
par(mfrow=c(1,2),mar=c(4,2,4,4))
hist(res_raw$pvalue, main=paste0("data without NA FDRtool sigma: ",sigma), xlab="p-values",cex.main=1,breaks = 50)
hist(res_corrected$pval, main="data p-value corrected", xlab="p-values",cex.main=1,breaks = 50)
title(main="Fixing overestimation of the variance", outer=T, line=-1.4, cex.main=1)
dev.off()


######------------------------------------------------------------------######
##                             Visualize results                            ##
######------------------------------------------------------------------######
plot <- res_shrunk


######------  MA plot  ------###### 
maplot <- function (res, sigthresh=0.1, lfcthresh=1, labelsig=TRUE, textcx=0.8,legendpos="bottomleft", ...) {
  with(subset(res, padj>sigthresh), plot(baseMean, log2FoldChange, pch=20,col="lightgrey", cex=1, log="x", ...))
  with(subset(res, padj<=sigthresh & abs(log2FoldChange)<1), points(baseMean, log2FoldChange, col="grey30", pch=20, cex=1.7))
  with(subset(res, padj<=sigthresh & abs(log2FoldChange)>=1 & abs(log2FoldChange) < 2), points(baseMean, log2FoldChange, col="orange", pch=20, cex=1.7))
  with(subset(res, padj<=sigthresh & abs(log2FoldChange)>=lfcthresh), points(baseMean,    log2FoldChange, col="darkorange", pch=20, cex=1.7))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh), textxy(baseMean, log2FoldChange, labs=name, cex=textcx, col="black"))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh," |LogFC|<1",sep=""), paste("1>=|LogFC|<",lfcthresh,sep=""), paste("|LogFC|>=",lfcthresh,sep="")), pch=16, col=c("grey30","orange","darkorange"))
}
pdf(file = paste(output,"diffexpr-maplott.pdf",sep=""), pointsize=2)
y_range <- round(max(abs(plot[!is.na(plot$padj<=0.1),"log2FoldChange"])),digits = 0)+1
maplot(plot,lfcthresh=2,sigthresh=0.1, ylim=c(-y_range,y_range), main=paste("Nr. DEG:",sum(plot$padj<=0.1)))
dev.off()


######------  Volcano plot  ------###### 
# replace 0 p value of KO genes with 1e-10
plot[!is.na(plot$padj) & plot$padj<1e-10,"padj"]<-1e-10
plot <- plot[abs(plot$log2FoldChange)>=0, ]
# volcano plot
pdf(file = paste0(output,"diffexpr-volcanoplot.pdf"), pointsize=2,paper="a4r",width=6, height=8)
options(ggrepel.max.overlaps = Inf)
x_range <- round(max(abs(plot[plot$padj<=0.1,"log2FoldChange"])),digits = 0)+1
y_range <- round(max(-log10(plot$padj)),digits = 0)+1
EnhancedVolcano(plot,
                lab = plot$name,
                x = "log2FoldChange",
                y = "padj",
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                axisLabSize = 8,
                pCutoff = 0.01,     
                FCcutoff = 1,    
                xlim=c(-x_range, x_range),
                ylim=c(0, y_range),
                title = "",
                legendPosition = "top",
                legendLabSize = 10,
                labSize = 2,
                col = c("lightgrey", "grey30", "grey30", "#84C9B0"),
                colAlpha=0.7,
                pointSize = 4,
                widthConnectors=0.2,
                gridlines.major=F,
                gridlines.minor=F,
                borderWidth=1,
                cutoffLineWidth=0.2,
                shape=16,
                drawConnectors=T,
                lengthConnectors = unit(0.001, "npc"),typeConnectors = "open")
dev.off()

