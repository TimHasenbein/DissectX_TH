##*****************************##
##  Allelome.PRO Violin plots  ##
##*****************************##
## Project: dissectX
## Tim Hasenbein 
## Last modification 04.2024
## Creation: 04.2022 
## AR violin plot for X-linked genes
library(tidyverse)
library(ggplot2)
library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(gdata)
library(ggpubr)
library(ggrepel)
library(circlize)
library(xlsx)


######------  set environment  ------###### 
input <- "./08_AR_pseudobulk_TKO_WT/"
output <- "./08_AR_pseudobulk_TKO_WT/"

 
######------  get AP results and filter  ------######
samples <- c("WT","KO")
for (sample in samples){
  locus <- fread(paste0(input,"2024_04_24_",sample,"_annotation_us.bed_1/locus_table.txt"), header = T)
  locus <- locus%>%dplyr::filter(chr == "chrX" & total_reads >= 30)
  locus <- locus%>%dplyr::select("name","allelic_ratio")
  colnames(locus) <- c("name",paste0(sample))
  assign(paste0(sample),locus)
}


######------ Generate matrix for plots  ------###### 
res <- WT%>%merge(KO,by="name",all=T)
# filter rows with no NA 
res <- na.omit(res) # 302 genes
res <- res%>%dplyr::filter(!name %in% c("Firre","Gm35612","4933407K13Rik")) # 300 genes
write.xlsx(res, paste0(output,"pseudobulk_AR.xlsx"), 
           col.names = TRUE, row.names = FALSE, append = FALSE)
rm(locus,sample,samples,WT,KO)


#####------  Violin plot X-linked genes  ------###### 
plot <- res %>% pivot_longer(!name,names_to = "sample", values_to = "Allelic_ratio")
plot$sample <- factor(plot$sample, levels = c("WT","KO")) 
pdf(file=paste0(output,"violin_plot_xgenes_totalreads.pdf"), width=6, height=5)
ggplot(plot, aes(x=sample, y=Allelic_ratio,fill=sample)) + 
  geom_hline(yintercept=0.7, linetype="dashed", color = "grey") + 
  geom_hline(yintercept=0.3, linetype="dashed", color = "grey") +
  geom_violin(width = 0.5, size = 0.1, scale = "width") +
  labs(x="Genotype",y="Allelic ratio") +
  geom_boxplot(width = 0.3, size=0.2, alpha = 0.8, outlier.shape=NA) +
  labs(x="",y="Allelic ratio") +
  scale_fill_manual(values = c(alpha("grey",alpha=0.6),alpha("#93CDB5",alpha=0.6))) +
  scale_x_discrete(limits=c("WT","KO")) +
  geom_point(data=plot%>%dplyr::filter(name == "Xist"), 
             aes(x=sample, y=Allelic_ratio,fill=sample), 
             color="black",
             size=1.5) +
  scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0,0.3,0.5,0.7,1)) +
  theme(axis.line = element_line(colour = 'black', size = 0.5),
        axis.text.y=element_text(size=10, colour = 'black'),
        axis.text.x=element_text(size=10, colour = 'black'),
        axis.title.y=element_text(size=12, angle = 90, colour = 'black'),
        axis.title.x=element_text(size=12, colour = 'black'),
        panel.background = element_blank(),
        panel.grid.minor =  element_line(colour = "grey",size=0.1),
        panel.grid.major.x = element_line(colour = "grey",size=0.1),
        panel.grid.major =  element_line(colour = "grey",size=0.1),
        legend.position="none") 
dev.off()

