##******************************************************##
## scRNA-seq: Split bam file according to Xa Xi status ##
##**********'*******************************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 01.2023
## Creation: 04.2022
## split bam file according to Xa Xi status
library(Seurat)
library(tidyverse)


######------ set environment ------###### 
input <- "./"
output <- "./"


######------ Read data ------###### 
bl6_Xa_ko <- read.table(paste0(input,"bl6_+-_Xa.txt"), header=T)
bl6_Xa_wt <- read.table(paste0(input,"bl6_++_Xa.txt"), header=T)
cast_Xa_ko <- read.table(paste0(input,"cast_+-_Xa.txt"), header=T)
cast_Xa_wt <- read.table(paste0(input,"cast_++_Xa.txt"), header=T)


######------ create "cell" file for KO/WT bam females ------###### 
bl6_Xa_ko <- bl6_Xa_ko %>% select("barcode", "GT")
bl6_Xa_ko$GT <- paste0(bl6_Xa_ko$GT,"_bl6_Xa")
write.table(bl6_Xa_ko, paste0(output,"cellfile_bl6_Xa_ko.txt"), quote=FALSE, sep = "\t", col.names = F, row.names =F)
  
bl6_Xa_wt <- bl6_Xa_wt %>% select("barcode", "GT")
bl6_Xa_wt$GT <- paste0(bl6_Xa_wt$GT,"_bl6_Xa")
write.table(bl6_Xa_wt, paste0(output,"cellfile_bl6_Xa_wt.txt"), quote=FALSE, sep = "\t", col.names = F, row.names =F)

cast_Xa_ko <- cast_Xa_ko %>% select("barcode", "GT")
cast_Xa_ko$GT <- paste0(cast_Xa_ko$GT,"_cast_Xa")
write.table(cast_Xa_ko, paste0(output,"cellfile_cast_Xa_ko.txt"), quote=FALSE, sep = "\t", col.names = F, row.names =F)

cast_Xa_wt <- cast_Xa_wt %>% select("barcode", "GT")
cast_Xa_wt$GT <- paste0(cast_Xa_wt$GT,"_cast_Xa")
write.table(cast_Xa_wt, paste0(output,"cellfile_cast_Xa_wt.txt"), quote=FALSE, sep = "\t", col.names = F, row.names =F)


