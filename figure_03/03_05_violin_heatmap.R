##*****************************##
##  Allelome.PRO Violin plots  ##
##*****************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 09.2022
## Creation: 08.2021
## Visualize results from the Allelome.PRO runs of the placenta samples
## Filter results for min reads overlapping a gene >30; omit rows with NA
library(tidyverse)
library(ggplot2)
library(data.table)
library(gdata)
library(ComplexHeatmap)


######------  set environment  ------###### 
input <- "./03_RNAseq_AllelomePRO/placenta/"
output <- "./violin_plots/"


######------  get AP results  ------###### 
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
# filter for chrX genes and genes including >= 30 reads
filter <- function(x) {
  x <- x[x$chr == "chrX",]
  x <- x[x$total_reads >= 30,]
}
dataFiles_filter <- lapply(dataFiles, filter)


######------ Generate matrix for plots  ------###### 
# get gene position start 
gene_start <- dataFiles$`Pl_CxDF_++_1`%>%dplyr::select("name","start")
# select allelic ratio 
select <- function(x) {
  x <- x[,c("name","allelic_ratio")]
}
dataFiles_filter_select <- lapply(dataFiles_filter, select)
# tidy data frame
res <- bind_rows(dataFiles_filter_select, .id = "sample")
res <- res %>% pivot_wider(names_from = name, values_from = allelic_ratio)
res <- res %>% remove_rownames %>% column_to_rownames(var="sample")
res_tidy <- as.data.frame(t(res))
res_tidy <- res_tidy%>%dplyr::select(-"Pl_CxXFp_++_1",-"Pl_CxXFp_++_2",-"Pl_CxXFp_++_3",
                                     -"Pl_CxXFp_+-_1",-"Pl_CxXFp_+-_2",-"Pl_CxXFp_+-_3",
                                     -"Pl_XFpxC_++_1",-"Pl_XFpxC_++_2",-"Pl_XFpxC_++_3",
                                     -"Pl_XFpxC_-+_1",-"Pl_XFpxC_-+_2",-"Pl_XFpxC_-+_3")
# filter rows with no NA 
res_tidy <- na.omit(res_tidy) # 315 out of 551 genes
saveRDS(res_tidy, file = paste0(output,"replicates_APrun_pl.rds"))
rm(dataFiles_filter,dataFiles_filter_select,sample_names)


######------  get median of values  ------###### 
# calculate median for each gene from the pooled subgroups
# group replicates with same GT
group <- setDT(res_tidy, keep.rownames = TRUE)[]
group <- group%>%pivot_longer(!rn, names_to="GT",values_to="allelic_ratio" )
group$GT <- gsub("_1","",group$GT)
group$GT <- gsub("_2","",group$GT)
group$GT <- gsub("_3","",group$GT)
group$GT <- gsub("\\+\\+","WT",group$GT)
unique(group$GT)
# same name for all controls from forward/reverse cross
for(r in 1:nrow(group)){
  group[r,"GT"] <- ifelse(str_detect(group[r,"GT"],"Cx")&&str_detect(group[r,"GT"],"WT"), "CxGT_WT", ifelse(str_detect(group[r,"GT"],"xC")&&str_detect(group[r,"GT"],"WT"), "GTxC_WT", group[r,"GT"]))
}
unique(group$GT)
# calculate medians of replicates
median <- group %>% dplyr::group_by(GT,rn) %>%dplyr::summarise(median_ratio=median(allelic_ratio))
# reorder data frame and sort according to position
median <- median %>% pivot_wider(names_from = GT, values_from = median_ratio)
median <- median %>% inner_join(gene_start, by = c("rn"="name")) 
median$start <- as.numeric(median$start)
median <- median %>% remove_rownames %>% column_to_rownames(var="rn")
median <- median[order(median$start),] 


######------    Violin plot: Deletion on Xa (GTxC)    ------###### 
# generate matrix 
xa <- median[ ,str_detect(colnames(median),"xC|GTxC_WT|start")] 
xa <- xa%>%dplyr::select(-"start") 
xa <- 1-xa # match maternal color
xa <- as.data.frame(t(xa)) 
xa <- setDT(xa, keep.rownames = TRUE)[]
order <- c("GTxC_WT","Pl_FixC_-+","Pl_DxC_-+","Pl_DFxC_-+","Pl_XFLxC_-+","Pl_LFxC_-+","Pl_LFDxC_-+")
# specify row order
xa$rn <- reorder.factor(xa$rn, new.order=order) 
xa <- xa %>% arrange(rn)
xa$rn <- gsub("Pl_","",xa$rn)
xa$rn <- gsub("xC_\\-\\+","",xa$rn)
xa$rn <- gsub("GTxC_","",xa$rn)
xa <- xa %>% dplyr::select(-"Tsix") # remove Tsix
saveRDS(xa, file = paste0(output,"median_Xa_APrun_pl.rds"))
# Violin plot   
violin_xa <- xa %>% pivot_longer(!rn,names_to = "gene", values_to = "Allelic_ratio")
violin_xa$rn <- factor(violin_xa$rn, levels = c("WT","Fi","D","DF","XFL","LF","LFD")) 
pdf(file=paste0(output,"violin_xa.pdf"), paper="a4",width=10, height=4)
ggplot(violin_xa, aes(x=rn, y=Allelic_ratio,fill=rn)) + 
  geom_violin() +
  labs(x="",y="Allelic ratio") +
  scale_fill_manual(values = c("black","#DB5D14","#808181",
                               "#3C73B9","#B81B1A","#12672F","#93CDB5")) +
  scale_x_discrete(limits=c("WT","Fi","D","DF","XFL","LF","LFD"),position = "top") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0,0.5,1)) +
  geom_boxplot(width=0.15,fill="white") + 
  geom_point(data=violin_xa%>%dplyr::filter(gene == "Xist"), 
             aes(x=rn, y=Allelic_ratio,fill=rn), 
             color="#053061",
             size=3) +
  theme(axis.ticks.x=element_blank(),
    axis.title.x=element_blank(),
    axis.text.x = element_text(angle = 45, hjust=0.1),
    axis.line = element_line(colour = 'black', size = 0.2),
    axis.ticks.y=element_line(colour = "black",size=0.2),
    axis.title.y=element_text(size=12),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    legend.position="none") 
dev.off()


######------    Violin plot: Deletion on Xi (CxGT)    ------###### 
# generate matrix 
xi <- median[ ,str_detect(colnames(median),"Cx|CxGT_WT|start")] 
xi <- xi%>%dplyr::select(-"start") 
xi <- as.data.frame(t(xi)) 
xi <- setDT(xi, keep.rownames = TRUE)[]
order <- c("CxGT_WT","Pl_CxFi_+-","Pl_CxD_+-","Pl_CxDF_+-","Pl_CxXFL_+-","Pl_CxLF_+-","Pl_CxLFD_+-")
# specify row order
xi$rn <- reorder.factor(xi$rn, new.order=order) 
xi <- xi %>% arrange(rn)
xi$rn <- gsub("Pl_","",xi$rn)
xi$rn <- gsub("Cx","",xi$rn)
xi$rn <- gsub("_\\+\\-","",xi$rn)
xi$rn <- gsub("GT_","",xi$rn)
xi <- xi %>% dplyr::select(-"Tsix") # remove Tsix
saveRDS(xi, file = paste0(output,"median_Xi_APrun_pl.rds"))
# Violin plot   
violin_xi <- xi %>% pivot_longer(!rn,names_to = "gene", values_to = "Allelic_ratio")
violin_xi$rn <- factor(violin_xi$rn, levels = c("WT","Fi","D","DF","XFL","LF","LFD")) 
pdf(file=paste0(output,"violin_xi.pdf"), paper="a4",width=10, height=4)
ggplot(violin_xi, aes(x=rn, y=Allelic_ratio,fill=rn)) + 
  geom_violin() +
  labs(x="",y="Allelic ratio") +
  scale_fill_manual(values = c("black","#DB5D14","#808181",
                               "#3C73B9","#B81B1A","#12672F","#93CDB5")) +
  scale_x_discrete(limits=c("WT","Fi","D","DF","XFL","LF","LFD"),position = "top") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0,0.5,1)) +
  geom_boxplot(width=0.15,fill="white") + 
  geom_point(data=violin_xi%>%dplyr::filter(gene == "Xist"), 
             aes(x=rn, y=Allelic_ratio,fill=rn), 
             color="#053061",
             size=3) +
  theme(axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45, hjust=0.1),
        axis.line = element_line(colour = 'black', size = 0.2),
        axis.ticks.y=element_line(colour = "black",size=0.2),
        axis.title.y=element_text(size=12),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position="none") 
dev.off()



######------  Plot AR heatmap  ------###### 
# generate matrix Xa
heatmap <- median
heatmap <- heatmap%>%dplyr::select(-"start") 
# specify row order
heatmap <- heatmap%>%dplyr::select(c("GTxC_WT","Pl_FixC_-+","Pl_DxC_-+","Pl_DFxC_-+","Pl_XFLxC_-+","Pl_LFxC_-+","Pl_LFDxC_-+",
                                     "CxGT_WT","Pl_CxFi_+-","Pl_CxD_+-","Pl_CxDF_+-","Pl_CxXFL_+-","Pl_CxLF_+-","Pl_CxLFD_+-"))
pdf(file=paste0(output,"AR_heatmap.pdf"),width=5, height=10)
colfunc <- circlize::colorRamp2(c(0, 0.5, 1), c("black","white","#B87333")) 
ht <- Heatmap(as.matrix(heatmap), 
              cluster_rows = F, 
              cluster_columns = F, 
              column_title = NULL, 
              show_row_names = T, 
              row_names_gp =  gpar(fontsize = 4),
              col = colfunc, 
              width = unit(10, "cm"), border = T,
              row_gap = unit(0.5, "mm"), 
              row_title_gp = gpar(fontsize = 6, col="black"),
              column_names_gp = gpar(fontsize = 10),
              heatmap_legend_param = list(at = c(0,0.25,0.5,0.75,1))) 
draw(ht, heatmap_legend_side = "left")
dev.off()
