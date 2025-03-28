##*************************************##
##  Calculate TPMs for adult tissues   ##
##*************************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 01.2022
## Creation: 10.2021
## output TPM table for adult tissues
library(tidyverse)
library(ggplot2)
library(biomaRt)
library(data.table)
library(stats)
library(GenomicFeatures)
library(openxlsx)


######------ Set environment  ------###### 
input <- "./"
output <- "./TPMs/"
annotation <- "./annotation.gtf"
mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="http://aug2020.archive.ensembl.org")


######------ Get data and make count matrix  ------######
count_table <- data.frame(name=character())
samples <- c((list.files(input)))
for (sample in samples){
  try(data <- read.table(paste0(input,sample)))
  try(colnames(data) <- c("name",gsub("_unstranded.count","",sample)))
  try(count_table <- count_table %>% merge(data, by="name",all=T))
}
count_table <- count_table[-c(1:5),]


######------ Get gene lengths  ------######
# makes TxDb object
txdb <- makeTxDbFromGFF(annotation,format="gtf") 
# collect exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")
# for each gene, reduce all exons to a set of non overlapping exons, calculate their lengths and sum then
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))
# merge count_table with extracted sizes
exonic.gene.sizes <- data.frame(exonic.gene.sizes)
exonic.gene.sizes$refseq_id <- row.names(exonic.gene.sizes)
colnames(exonic.gene.sizes) <- c("size","name")
count_table <- count_table%>%inner_join(exonic.gene.sizes,by="name")


######------ Calculate TPMs  ------######
tpm <- function(counts, lengths) {
  rate <- counts / lengths 
  rate / sum(rate) * 1e6
}
genes <- data.frame(
  Gene = count_table$name,
  Length = count_table$size)
TPM_table <- as.data.frame(apply(count_table[2:95], 2, function(x) tpm(x, genes$Length)))  # adapt for new data (length -1 which is size)
TPM_table <- round(TPM_table,digits = 3)
TPM_table <- cbind(count_table$name,count_table$size,TPM_table)
colnames(TPM_table)[1:2] <- c("name","transcript_length")


######------ Add gene symbol  ------######
TPM_table$name <- sapply(strsplit(TPM_table$name, ".", fixed=T), function(x) x[1])
gm35612 <- TPM_table%>%filter(name=="Gm35612")
IDs <- getBM(attributes = c("ensembl_gene_id","mgi_symbol","chromosome_name","strand","start_position","end_position"),
             filters = "ensembl_gene_id", values = TPM_table[,1],
             mart = mart)
sum(duplicated(IDs$mgi_symbol))
TPM_table <- TPM_table%>%inner_join(IDs, by=c("name"="ensembl_gene_id"))
TPM_table <- TPM_table%>%dplyr::select("chromosome_name","mgi_symbol","start_position","end_position","strand", everything())
# add Gm35612
gm35612 <- gm35612%>%mutate("chromosome_name"="X", "mgi_symbol"="Gm35612","start_position"=50537389, "end_position"=50562598, "strand"=1) 
gm35612 <- gm35612%>%dplyr::select("chromosome_name","mgi_symbol","start_position","end_position","strand","name","transcript_length",everything())
TPM_table <- TPM_table%>%rbind(gm35612)
# rm genes with no mgi symbol
TPM_table <- TPM_table%>%filter(mgi_symbol != "")
# substitute strand -1 and 1 to + and -
TPM_table$strand <- gsub("\\-1","-",TPM_table$strand)
TPM_table$strand <- gsub("1","+",TPM_table$strand)
# sort table
TPM_table <- TPM_table%>%arrange(chromosome_name,start_position)


######------ Write output  ------######
write.table(TPM_table, paste0(output,"TPM_table_adult_replicates.txt"), quote=FALSE, sep = "\t", col.names = TRUE, row.names =FALSE)


######------ Calculate TPM means  ------######
TPM_table <- fread(paste0(output,"TPM_table_adult_replicates.txt"),header = T)
TPM_groups <- pivot_longer(TPM_table,cols=8:101,names_to="sample",values_to="TPMs")
TPM_groups$sample <- gsub("_1","",TPM_groups$sample)
TPM_groups$sample <- gsub("_2","",TPM_groups$sample)
TPM_groups$sample <- gsub("_3","",TPM_groups$sample)
TPM_groups$sample <- gsub("_4","",TPM_groups$sample)
TPM_groups$sample <- gsub("_XFLxXFL_\\+\\+","_DFxDF_\\+\\+",TPM_groups$sample)
TPM_groups$sample <- gsub("_LFxLF_\\+\\y","_DFxDF_\\+\\y",TPM_groups$sample)
unique(TPM_groups$sample)
# group by sample and gene and calculate TPM means
TPM_groups <- TPM_groups%>%dplyr::select("mgi_symbol","sample","TPMs")
mean <- TPM_groups %>% group_by(sample, mgi_symbol) %>% summarize(mean = mean(TPMs))
mean <- mean%>%pivot_wider(names_from = sample,values_from = mean)
write.table(mean, paste0(output,"TPM_table_adult_means.txt"), quote=FALSE, sep = "\t", col.names = TRUE, row.names =FALSE)











