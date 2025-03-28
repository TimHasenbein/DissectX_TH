##*********************************************##
## scRNA-seq: Split bam file into per cell bam ##
##**********'**********************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 01.2023
## Creation: 04.2022
## Perform on raw data, so we can select clean cells later, regardless of the QC parameters and 
## do not have to run AP all over again
library(Seurat)
library(tidyverse)


######------ Set environment  ------###### 
input <- "./00_cellranger_data/"
output <- "./AllelomePRO_2023/"


######------ Read data and make seurat object ------###### 
samples <- c("LFDxCast_WT","LFDxCast_KO")
for (sample in samples){
  seurat_data <- Read10X(data.dir = paste0(input,sample))
  seurat_obj <- CreateSeuratObject(counts = seurat_data$`Gene Expression`, 
                                   project = sample) 
  assign(sample, seurat_obj)}


######------  Merge data ------###### 
merged_seurat <- merge(x = LFDxCast_WT, y = LFDxCast_KO, project = "dissectX")
merged_seurat$cell <- rownames(merged_seurat@meta.data)
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'orig.ident', into = c('cross', 'condition'), sep = '_')
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)


######------ create "cell" file for KO/WT bam females ------###### 
meta <- merged_seurat@meta.data 
# KO bam cell file
KO <- meta%>%dplyr::filter(condition=="KO")%>%dplyr::select("cell")
KO$ID <- paste0(KO$cell,"_KO") # add iD column
write.table(KO, paste0(output,"cellfile_female_KO.txt"), quote=FALSE, sep = "\t", col.names = F, row.names =F)
# WT bam cell file
WT <- meta%>%dplyr::filter(condition=="WT")%>%dplyr::select("cell")
WT$ID <- paste0(WT$cell,"_WT")
write.table(WT, paste0(output,"cellfile_female_WT.txt"), quote=FALSE, sep = "\t", col.names = F, row.names =F)




