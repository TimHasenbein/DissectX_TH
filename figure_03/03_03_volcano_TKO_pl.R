##***************************************##
## Volcano plot for TKO results placenta ##
##***************************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 09.2022
## Creation: 12.2021
## Volcano plot for the TKO results in placenta
library(tidyverse)
library(ggplot2)
library(biomaRt)
library(data.table)
library(EnhancedVolcano)


######------ Set environment  ------###### 
input <- "./02_RNAseq_DESeq2/placenta/"
output <- "./"
mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="http://aug2020.archive.ensembl.org")


######------ Get data  ------######
# Xa
GT <- "LFDxC"
# Xi
GT <- "CxLFD"
# import 
data <- read.table(paste0(Sys.glob(paste0(input,"/",GT,"/",GT,"_shrunk_diffexpr-results*"))))


######------  Convert ids to symbol  ------###### 
data <- setDT(as.data.frame(data), keep.rownames = TRUE)[]
data$rn <- sapply(strsplit(data$rn, ".", fixed=T), function(x) x[1])
gm35612 <- data%>%filter(rn=="Gm35612")
IDs <- getBM(attributes = c("ensembl_gene_id","mgi_symbol"),
             filters = "ensembl_gene_id", values = data[,1],
             mart = mart)
data <- data%>%inner_join(IDs, by=c("rn"="ensembl_gene_id"))
data <- data%>%filter(mgi_symbol!="")
#data <- data%>%rbind(gm35612)


######------  Volcano plot  ------###### 
# replace 0 p value of KO genes with 1e-10
data[!is.na(data$padj) & data$padj<1e-10,"padj"]<-1e-10
data <- data[abs(data$log2FoldChange)>=0, ]
# volcano plot
pdf(file = paste0(output,"Xa_volcanoplot.pdf"), pointsize=2,paper="a4r",width=6, height=8)
options(ggrepel.max.overlaps = Inf)
x_range <- round(max(abs(data[data$padj<=0.1,"log2FoldChange"])),digits = 0)+1
y_range <- round(max(-log10(data$padj)),digits = 0)+1
EnhancedVolcano(data,
                lab = data$mgi_symbol,
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

