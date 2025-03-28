##*****************************##
##  Allelome.PRO Violin plots  ##
##*****************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 09.2022
## Creation: 08.2021
library(tidyverse)
library(ggplot2)
library(data.table)
library(gdata)
library(ComplexHeatmap)


######------  set environment  ------###### 
input <- "./placenta/"
output <- "./"


######------  get AP results  ------###### 
read <- function(x) read.table(x, header=T)
dataFiles <- lapply(Sys.glob(paste0(input,"*/locus_table.txt")), read)
sample_names <- Sys.glob(paste0(input,"*/locus_table.txt"))
sample_names <- gsub(paste0(input,"2021_08_11_"),"",sample_names)
sample_names <- gsub(paste0(input,"2021_08_10_"),"",sample_names)
sample_names <- gsub(paste0(input,"2021_08_09_"),"",sample_names)
sample_names <- gsub("_annotation_us.bed_1/locus_table.txt","",sample_names)
names(dataFiles) <- sample_names
sample_names[duplicated(sample_names)]
filter <- function(x) {
  x <- x[x$chr == "chrX",]
  x <- x[x$total_reads >= 30,]
}
dataFiles_filter <- lapply(dataFiles, filter)


######------ Generate matrix for plots  ------###### 
gene_start <- dataFiles$`Pl_CxDF_++_1`%>%dplyr::select("name","start")
select <- function(x) {
  x <- x[,c("name","allelic_ratio")]
}
dataFiles_filter_select <- lapply(dataFiles_filter, select)
res <- bind_rows(dataFiles_filter_select, .id = "sample")
res <- res %>% pivot_wider(names_from = name, values_from = allelic_ratio)
res <- res %>% remove_rownames %>% column_to_rownames(var="sample")
res_tidy <- as.data.frame(t(res))
# filter rows with no NA 
res_tidy <- na.omit(res_tidy) # 214 genes
rm(dataFiles_filter,dataFiles_filter_select,sample_names)


######------  select WT samples only   ------###### 
group <- setDT(res_tidy, keep.rownames = TRUE)[]
group <- group%>%dplyr::select("rn","Pl_CxDF_++_1","Pl_CxDF_++_2","Pl_CxDF_++_3",,"Pl_CxD_++_1",
                               "Pl_CxD_++_2","Pl_CxD_++_3","Pl_CxFi_++_1","Pl_CxFi_++_2","Pl_CxFi_++_3",
                               "Pl_DFxC_++_1","Pl_DFxC_++_2","Pl_DFxC_++_3","Pl_DxC_++_1","Pl_DxC_++_2",
                               "Pl_FixC_++_1","Pl_FixC_++_2","Pl_FixC_++_3")
group <- group%>%pivot_longer(!rn, names_to="GT",values_to="allelic_ratio" )
group$GT <- gsub("\\+\\+","WT",group$GT)
unique(group$GT)


######------  plot Xa samples only   ------###### 
xa <- group%>%dplyr::filter(GT %in% c("Pl_DFxC_WT_1","Pl_DFxC_WT_2","Pl_DFxC_WT_3","Pl_DxC_WT_1","Pl_DxC_WT_2",
                                      "Pl_FixC_WT_1","Pl_FixC_WT_2","Pl_FixC_WT_3"))
xa$GT <- factor(xa$GT, levels = c("Pl_DFxC_WT_1","Pl_DFxC_WT_2","Pl_DFxC_WT_3","Pl_DxC_WT_1","Pl_DxC_WT_2",
                                  "Pl_FixC_WT_1","Pl_FixC_WT_2","Pl_FixC_WT_3"))
xa$allelic_ratio <- 1 - xa$allelic_ratio
# violin plot
pdf(file=paste0(output,"violin_plot_wt_xa_GTxC.pdf"), width=8, height=5)
ggplot(xa, aes(x=GT, y=allelic_ratio,fill=GT)) + 
  geom_boxplot(width = 0, size=0, alpha = 0.8, fill = NA, outlier.size = 1.2) +
  geom_violin(width = 0.5, size = 0.1, scale = "width", fill = alpha("black",alpha=0.4)) +
  labs(x="",y="Allelic ratio") +
  geom_point(data=xa%>%dplyr::filter(rn == "Xist"), 
             aes(x=GT, y=allelic_ratio,fill=GT), 
             color="darkred",
             size=1.2) +
  scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0,0.5,1)) +
  theme(axis.line = element_line(colour = 'black', size = 0.5),
        axis.text.y=element_text(size=10, colour = 'black'),
        axis.text.x=element_text(size=10, angle = 45, colour = 'black', vjust = 1, hjust=1),
        axis.title.y=element_text(size=12, angle = 90, colour = 'black'),
        axis.title.x=element_text(size=12, colour = 'black'),
        panel.background = element_blank(),
        panel.grid.minor =  element_line(colour = "grey",size=0.1),
        panel.grid.major.x = element_line(colour = "grey",size=0.1),
        panel.grid.major =  element_line(colour = "grey",size=0.1),
        legend.position="none") 
dev.off()


######------  plot Xi samples only   ------###### 
xi <- group%>%dplyr::filter(GT %in% c("Pl_CxDF_WT_1","Pl_CxDF_WT_2","Pl_CxDF_WT_3","Pl_CxD_WT_1",
                                      "Pl_CxD_WT_2","Pl_CxD_WT_3","Pl_CxFi_WT_1","Pl_CxFi_WT_2","Pl_CxFi_WT_3"))
xi$GT <- factor(xi$GT, levels =c("Pl_CxDF_WT_2","Pl_CxDF_WT_3","Pl_CxD_WT_1",
                                 "Pl_CxD_WT_2","Pl_CxD_WT_3","Pl_CxFi_WT_1","Pl_CxFi_WT_2","Pl_CxDF_WT_1","Pl_CxFi_WT_3"))
xi$allelic_ratio <- 1 - xi$allelic_ratio
pdf(file=paste0(output,"violin_plot_wt_xi_CxGT.pdf"), width=8, height=5)
# violin plot
ggplot(xi, aes(x=GT, y=allelic_ratio,fill=GT)) + 
  geom_violin(width = 0.5, size = 0.1, scale = "width", fill = alpha("#C6874F",0.5)) +
  geom_boxplot(width = 0, size=0., alpha = 0.8, fill = NA, outlier.size = 1.2) +
  labs(x="",y="Allelic ratio") +
  geom_point(data=xi%>%dplyr::filter(rn == "Xist"), 
             aes(x=GT, y=allelic_ratio,fill=GT), 
             color="darkred",
             size=1.2) +
  scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0,0.5,1)) +
  theme(axis.line = element_line(colour = 'black', size = 0.5),
        axis.text.y=element_text(size=10, colour = 'black'),
        axis.text.x=element_text(size=10, angle = 45, colour = 'black', vjust = 1, hjust=1),
        axis.title.y=element_text(size=12, angle = 90, colour = 'black'),
        axis.title.x=element_text(size=12, colour = 'black'),
        panel.background = element_blank(),
        panel.grid.minor =  element_line(colour = "grey",size=0.1),
        panel.grid.major.x = element_line(colour = "grey",size=0.1),
        panel.grid.major =  element_line(colour = "grey",size=0.1),
        legend.position="none") 
dev.off()



######------  plot heatmap  ------###### 
heat <- group %>% pivot_wider(names_from = GT, values_from = allelic_ratio)
heat <- heat%>%dplyr::select("rn","Pl_CxDF_WT_2","Pl_CxDF_WT_3","Pl_CxD_WT_1",
                             "Pl_CxD_WT_2","Pl_CxD_WT_3","Pl_CxFi_WT_1","Pl_CxFi_WT_2","Pl_CxDF_WT_1","Pl_CxFi_WT_3",
                             "Pl_DFxC_WT_1","Pl_DFxC_WT_2","Pl_DFxC_WT_3","Pl_DxC_WT_1","Pl_DxC_WT_2",
                             "Pl_FixC_WT_1","Pl_FixC_WT_2","Pl_FixC_WT_3")
heat <- heat %>% column_to_rownames("rn")
# filter rows for AR > 0.1 & < 0.9
#heat_fil <- heat[apply(heat, 1, function(row) any(row > 0.2 & row < 0.8)),]
pdf(file=paste0(output,"AR_heatmap_wt_pl.pdf"),width=5, height=10)
colfunc <- circlize::colorRamp2(c(0, 0.5, 1), c("black","white","#B87333")) 
ht <- Heatmap(as.matrix(heat), 
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


######------  calculate median and plot AR delta 0.1  ------###### 
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
median$GTxC_WT <- 1- median$GTxC_WT
median$delta <- median$CxGT_WT - median$GTxC_WT
# add label if abs(delta) >= 0.1
median$label <- rownames(median)
median <- median%>%mutate(label=ifelse(abs(delta)>=0.1,label,""))
# plot delta
colfunc <- circlize::colorRamp2(c(-0.5,0,0.5), c("black","white","#B87333")) 
pdf(file=paste0(output,"AR_delta_genes_heatmap_wt.pdf"),width=5, height=10)
ht <- Heatmap(as.matrix(median%>%dplyr::select("delta")), 
              cluster_rows = F, 
              cluster_columns = F, 
              column_title = NULL, 
              show_row_names = T, 
              col = colfunc,
              row_names_gp =  gpar(fontsize = 4),
              width = unit(10, "cm"), border = T,
              row_gap = unit(0.5, "mm"), 
              row_title_gp = gpar(fontsize = 6, col="black"),
              column_names_gp = gpar(fontsize = 10)) 
draw(ht, heatmap_legend_side = "left")
dev.off()
