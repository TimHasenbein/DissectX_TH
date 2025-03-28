##*****************************##
##  Allelome.PRO Violin plots  ##
##*****************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 04.2024
## Creation: 04.2024
library(tidyverse)
library(ggplot2)
library(data.table)
library(gdata)
library(ggpubr) 
library(xlsx)


######------  set environment  ------###### 
input <- "./"
output <- "./"


######------  get AP results and filter  ------######
samples <- samples <- c("Sp_6w_RNA_LFDxC_XX_-+_1","Sp_6w_RNA_LFDxC_XX_-+_2","Sp_6w_RNA_LFDxC_XX_-+_3","Sp_6w_RNA_LFDxC_XX_++_1","Sp_6w_RNA_LFDxC_XX_++_2","Sp_6w_RNA_LFDxC_XX_++_3")
for (sample in samples){
  locus <- fread(paste0(input,"2022_12_22_",sample,"_annotation_us.bed_1/locus_table.txt"), header = T)
  locus <- locus%>%dplyr::filter(chr == "chrX" & total_reads >= 30)
  locus <- locus%>%dplyr::select("name","allelic_ratio")
  colnames(locus) <- c("name",paste0(sample))
  assign(paste0(sample),locus)
}


######------ Generate matrix for plots  ------###### 
res <- `Sp_6w_RNA_LFDxC_XX_-+_1`%>%merge(`Sp_6w_RNA_LFDxC_XX_-+_2`,by="name",all=T)%>%merge(`Sp_6w_RNA_LFDxC_XX_-+_3`,by="name",all=T)%>%merge(`Sp_6w_RNA_LFDxC_XX_++_1`,by="name",all=T)%>%merge(`Sp_6w_RNA_LFDxC_XX_++_2`,by="name",all=T)%>%merge(`Sp_6w_RNA_LFDxC_XX_++_3`,by="name",all=T)
# filter rows with no NA 
res <- na.omit(res) # 288 genes
res <- res%>%dplyr::filter(!name %in% c("Firre","Gm35612","4933407K13Rik"))
write.xlsx(res, paste0(output,"replicates_APrun_sp_het.xlsx"), 
           col.names = TRUE, row.names = FALSE, append = FALSE)
saveRDS(res, file = paste0(output,"replicates_APrun_sp.rds"))
keep(res,input,output,sure=T)


#####------  Violin plot X-linked genes  ------###### 
plot <- res %>% pivot_longer(!name,names_to = "sample", values_to = "Allelic_ratio")
plot$sample <- factor(plot$sample, levels = c("Sp_6w_RNA_LFDxC_XX_++_1","Sp_6w_RNA_LFDxC_XX_++_2","Sp_6w_RNA_LFDxC_XX_++_3","Sp_6w_RNA_LFDxC_XX_-+_1","Sp_6w_RNA_LFDxC_XX_-+_2","Sp_6w_RNA_LFDxC_XX_-+_3")) 
pdf(file=paste0(output,"violin_plot_xgenes_totalreads_30_new.pdf"), width=8, height=5)
ggplot(plot, aes(x=sample, y=Allelic_ratio,fill=sample)) + 
  geom_violin(width = 0.5, size = 0.1, scale = "width") +
  labs(x="Genotype",y="Allelic ratio") +
  #geom_boxplot(width = 0.3, size=0.2, alpha = 0.8, outlier.shape=NA) +
  labs(x="",y="Allelic ratio") +
  scale_fill_manual(values = c(alpha("grey",alpha=0.6),alpha("grey",alpha=0.6),alpha("grey",alpha=0.6),
                               alpha("#93CDB5",alpha=0.6),alpha("#93CDB5",alpha=0.6),alpha("#93CDB5",alpha=0.6))) +
  scale_x_discrete(limits=c("Sp_6w_RNA_LFDxC_XX_++_1","Sp_6w_RNA_LFDxC_XX_++_2","Sp_6w_RNA_LFDxC_XX_++_3","Sp_6w_RNA_LFDxC_XX_-+_1","Sp_6w_RNA_LFDxC_XX_-+_2","Sp_6w_RNA_LFDxC_XX_-+_3"),labels=c("WT_1","WT_2","WT_3","LFD_1","LFD_2","LFD_3")) +
  geom_point(data=plot%>%dplyr::filter(name == "Xist"), 
             aes(x=sample, y=Allelic_ratio,fill=sample), 
             color="black",
             size=1.5) +
  geom_hline(yintercept=0.7, linetype="dashed", color = "grey") + 
  geom_hline(yintercept=0.3, linetype="dashed", color = "grey") +
  scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0,0.3,0.5,0.7,1)) +
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


######------  Boxplot Xist  ------###### 
plot_xist <- plot%>%dplyr::filter(name == "Xist")
plot_xist$sample <- gsub("_1","",plot_xist$sample)
plot_xist$sample <- gsub("_2","",plot_xist$sample)
plot_xist$sample <- gsub("_3","",plot_xist$sample)
pdf(file=paste0(output,"boxplot_xist_new.pdf"), width=5, height=4)
ggplot(plot_xist, aes(x=sample, y=Allelic_ratio,color=sample)) + 
  #geom_hline(yintercept=0.7, linetype="dashed", color = "grey") + 
  #geom_hline(yintercept=0.3, linetype="dashed", color = "grey") +
  stat_boxplot(geom ='errorbar',width= 0.2, color= "black") +
  labs(x="",y="Allelic ratio") +
  geom_point(size=2.5) +
  scale_color_manual(values = c("#93CDB5","black")) +
  scale_x_discrete(limits=c("Sp_6w_RNA_LFDxC_XX_++","Sp_6w_RNA_LFDxC_XX_-+"),labels=c("WT","LFD")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0,0.25,0.5,0.75,1),expand=c(0,0)) +
  theme(axis.line = element_line(colour = 'black', size = 0.5),
        axis.text.y=element_text(size=10, colour = 'black'),
        axis.text.x=element_text(size=10, colour = 'black'),
        axis.title.y=element_text(size=12, angle = 90, colour = 'black'),
        axis.title.x=element_text(size=12, colour = 'black'),
        panel.background = element_blank(),
        panel.grid.minor =  element_line(colour = "grey",size=0.1),
        panel.grid.major.x = element_line(colour = "grey",size=0.1),
        panel.grid.major =  element_line(colour = "grey",size=0.1),
        legend.position="none") +
  stat_compare_means(comparisons = list(c("Sp_6w_RNA_LFDxC_XX_++","Sp_6w_RNA_LFDxC_XX_-+")), method = "t.test",label.y = 0.9) 
dev.off()
