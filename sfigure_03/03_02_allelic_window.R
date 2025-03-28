##*****************************##
##  Allelome.PRO run placenta  ##
##*****************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 09.2022
## Creation: 08.2021
library(tidyverse)
library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(gdata)
library(ggpubr)
library(ggrepel)


######------ Set environment  ------###### 
input <- "./placenta/"
output <- "./allelic_window/"


######----------------------------------------------------------------######
##                             Data preparation                           ##
######----------------------------------------------------------------######


######------ Get sample results  ------###### 
# read in locus_table.txt for each sample
read <- function(x) read.table(x, header=T)
dataFiles <- lapply(Sys.glob(paste0(input,"*/locus_table.txt")), read)
# name each data frame
sample_names <- Sys.glob(paste0(input,"*/locus_table.txt"))
sample_names <- gsub(paste0(input,"2021_08_11_"),"",sample_names)
sample_names <- gsub(paste0(input,"2021_08_10_"),"",sample_names)
sample_names <- gsub(paste0(input,"2021_08_09_"),"",sample_names)
sample_names <- gsub("_annotation_us.bed_1/locus_table.txt","",sample_names)
names(dataFiles) <- sample_names
# check for duplicates 
sample_names[duplicated(sample_names)]
# filter for chrX genes and genes including >= 50 reads
filter <- function(x) {
  x <- x[x$chr == "chrX",]
  x <- x[x$total_reads >= 50,]
}
dataFiles_filter <- lapply(dataFiles, filter)


######------ Generate matrix for plots  ------###### 
# get gene position start 
gene_start <- dataFiles_filter$`Pl_CxDF_++_1`%>%dplyr::select("name","start")
# select allelic ratio 
select <- function(x) {
  x <- x[,c("name","allelic_ratio")]
}
dataFiles_filter_select <- lapply(dataFiles_filter, select)
# tidy data frame
res <- bind_rows(dataFiles_filter_select, .id = "sample")
res <- res %>% pivot_wider(names_from = name, values_from = allelic_ratio)
res <- res %>% remove_rownames %>% column_to_rownames(var="sample")
res_tidy <- t(res)
res_tidy <- as.data.frame(res_tidy)
rm(dataFiles_filter,dataFiles_filter_select,sample_names)


######------ Calculate medians and error bars  ------###### 
# filter for rows with no NA present
res_tidy <- na.omit(res_tidy)
# calculate median for each gene from the pooled subgroups
# group replicates for same GT
group <- setDT(res_tidy, keep.rownames = TRUE)[]
group <- group%>%pivot_longer(!rn, names_to="GT",values_to="allelic_ratio" )
group$GT <- gsub("_1","",group$GT)
group$GT <- gsub("_2","",group$GT)
group$GT <- gsub("_3","",group$GT)
group$GT <- gsub("\\+\\+","WT",group$GT)
# same name for all controls from forward/reverse
for(r in 1:nrow(group)){
  group[r,"GT"] <- ifelse(str_detect(group[r,"GT"],"Cx")&&str_detect(group[r,"GT"],"WT"), "CxGT_WT", ifelse(str_detect(group[r,"GT"],"xC")&&str_detect(group[r,"GT"],"WT"), "GTxC_WT", group[r,"GT"]))
}
# calculate medians of replicates
median <- group %>% dplyr::group_by(GT,rn) %>%dplyr::summarise(AR_median=median(allelic_ratio), AR_sd=sd(allelic_ratio))
# reorder data frame
median <- median %>% inner_join(gene_start, by = c("rn"="name"))
median$start <- as.numeric(median$start)


######------------------------------------------------------######
##                             GTxC                             ##
######------------------------------------------------------###### 


######------   select samples with deletion on Xa   ------###### 
samples <- c("Pl_LFDxC_-+","GTxC_WT")
Xa <- median%>%dplyr::filter(GT %in% samples)
Xa <- Xa[order(Xa$start),] # sort according to position


######------   Filter 3mb window around Firre   ------###### 
Xa_window_firre <- Xa%>%dplyr::filter(start >= 47563119 & start <= 53563119) 
Xa_window_firre$AR_median <- 1-Xa_window_firre$AR_median # match maternal color

######------   Filter 3mb window around Dxz4  ------###### 
Xa_window_dxz4 <- Xa%>%dplyr::filter(start >= 73525458 & start <= 77925458) 
Xa_window_dxz4$AR_median <- 1-Xa_window_dxz4$AR_median # match maternal color


######------   Plot   ------###### 
pdf(file = paste0(output,"Xa_crossfirreFirre_AR_3mb_window_min50reads_median.pdf"), pointsize=4,paper="a4",width=5, height=4)
ggplot(Xa_window_firre,aes(x=start,y=AR_median,label=rn, colour=GT)) +
  geom_errorbar(aes(ymin = AR_median - AR_sd, 
                    ymax = AR_median + AR_sd), size = 0.1, color="black", alpha=1) +
  geom_point(aes(fill=GT),size=4) + 
  geom_text(size=1,color="black") +
  labs(x="position on X",y="Allelic Ratio") +
  theme(axis.text.x=element_text(size=9,angle=45,vjust=1,hjust=1,colour = 'black'),
        axis.line = element_line(colour = 'black', size = 0.2),
        axis.ticks=element_line(colour = "black",size=0.2),
        axis.title=element_text(size=9),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.spacing = unit(1.0, 'cm')) +
  scale_color_manual(values=c("black", "#93CDB5")) + 
  scale_x_continuous(limits = c(47563119,53563119), breaks = c(47563119,48563119, 49563119,  50563119, 51563119, 52563119,53563119)) +
  scale_y_continuous(limits = c(0,1), breaks=c(0,0.25,0.5,0.75,1)) +
  geom_vline(xintercept = 50562787,size = 0.2,linetype="dashed") 
dev.off()

pdf(file = paste0(output,"Xa_dxz4_AR_2_2mb_window_min50reads__median.pdf"), pointsize=4,paper="a4",width=6, height=4)
ggplot(Xa_window_dxz4,aes(x=start,y=AR_median,label=rn, colour=GT)) +
  geom_errorbar(aes(ymin = AR_median - AR_sd, 
                    ymax = AR_median + AR_sd), size = 0.1, color="black", alpha=1) +
  geom_point(aes(fill=GT),size=4) + 
  geom_text(size=1,color="black") +
  labs(x="position on X",y="Allelic Ratio") +
  theme(axis.text.x=element_text(size=9,angle=45,vjust=1,hjust=1,colour = 'black'),
        axis.line = element_line(colour = 'black', size = 0.2),
        axis.ticks=element_line(colour = "black",size=0.2),
        axis.title=element_text(size=9),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.spacing = unit(1.0, 'cm')) +
  scale_color_manual(values=c("black", "#93CDB5")) + 
  scale_x_continuous(limits = c(73525458,77925458), breaks = c(73725458,74725458,75725458,76725458,77725458)) +
  scale_y_continuous(limits = c(0,1), breaks=c(0,0.25,0.5,0.75,1)) +
  geom_vline(xintercept = 75725458,size = 0.2,linetype="dashed") 
dev.off()

######------------------------------------------------------######
##                             CxGT                             ##
######------------------------------------------------------###### 


######------   select samples with deletion on Xi   ------###### 
samples <- c("Pl_CxLFD_+-","CxGT_WT")
Xi <- median%>%dplyr::filter(GT %in% samples)
Xi <- Xi[order(Xi$start),] # sort according to position


######------   Filter 3mb window around Firre   ------###### 
Xi_window_firre <- Xi%>%dplyr::filter(start >= 47563119 & start <= 53563119) 
Xi_window_firre$AR_median <- 1-Xi_window_firre$AR_median # match maternal color


######------   Filter 3mb window around Dxz4  ------###### 
Xi_window_dxz4 <- Xi%>%dplyr::filter(start >= 73525458 & start <= 77925458) 
Xi_window_dxz4$AR_median <- 1-Xi_window_dxz4$AR_median # match maternal color


######------   Plot   ------###### 
pdf(file = paste0(output,"Xi_crossfirreFirre_AR_3mb_window_min50reads_median.pdf"), pointsize=4,paper="a4",width=5, height=4)
ggplot(Xi_window_firre,aes(x=start,y=1-AR_median,label=rn, colour=GT)) +
  geom_errorbar(aes(ymin = (1-AR_median) - AR_sd, 
                    ymax = (1-AR_median) + AR_sd),size = 0.1, color="black", alpha=1) +
  geom_point(aes(fill=GT),size=4) + 
  geom_text(size=1,color="black") +
  labs(x="position on X",y="Allelic Ratio") +
  theme(axis.text.x=element_text(size=9,angle=45,vjust=1,hjust=1,colour = 'black'),
        axis.line = element_line(colour = 'black', size = 0.2),
        axis.ticks=element_line(colour = "black",size=0.2),
        axis.title=element_text(size=9),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.spacing = unit(1.0, 'cm')) +
  scale_color_manual(values=c("black", "#93CDB5")) + 
  scale_x_continuous(limits = c(47563119,53563119), breaks = c(47563119,48563119, 49563119,  50563119, 51563119, 52563119,53563119)) +
  scale_y_continuous(limits = c(0,1), breaks=c(0,0.25,0.5,0.75,1)) +
  geom_vline(xintercept = 50562787,size = 0.2,linetype="dashed") 
dev.off()

pdf(file = paste0(output,"Xi_dxz4_AR_2_2mb_window_min50reads_median.pdf"), pointsize=4,paper="a4",width=6, height=4)
ggplot(Xi_window_dxz4,aes(x=start,y=1-AR_median,label=rn, colour=GT)) +
  geom_errorbar(aes(ymin = (1-AR_median) - AR_sd, 
                    ymax = (1-AR_median) + AR_sd), size = 0.1, color="black", alpha=1) +
  geom_point(aes(fill=GT),size=4) + 
  geom_text(size=1,color="black") +
  labs(x="position on X",y="Allelic Ratio") +
  theme(axis.text.x=element_text(size=9,angle=45,vjust=1,hjust=1,colour = 'black'),
        axis.line = element_line(colour = 'black', size = 0.2),
        axis.ticks=element_line(colour = "black",size=0.2),
        axis.title=element_text(size=9),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.spacing = unit(1.0, 'cm')) +
  scale_color_manual(values=c("black", "#93CDB5")) + 
  scale_x_continuous(limits = c(73525458,77925458), breaks = c(73725458,  74725458,75725458,76725458,77725458)) +
  scale_y_continuous(limits = c(0,1), breaks=c(0,0.25,0.5,0.75,1)) +
  geom_vline(xintercept = 75725458,size = 0.2,linetype="dashed") 
dev.off()


