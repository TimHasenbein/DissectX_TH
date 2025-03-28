##*********************************##
##  scRNA-seq celltype composition ##
##*********************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 04.2024
## Creation: 04.2024
## Count the cells per celltype
library(tidyverse)
library(ggplot2)
library(data.table)
library(gdata)
library(Seurat)
library(ComplexHeatmap)


######------  set environment  ------###### 
output <- "./"
srat <- readRDS(paste0("./03_srat_clean_integrated_reg_MT_labelled_XX_KO_WT_Xa_status.rds")) # 3694 cells
DimPlot(srat, reduction = "umap", label = TRUE, label.size = 6, split.by = "Xa")

######------  plot celltype count ------###### 
meta <- srat@meta.data
cells <- as.data.frame(Idents(srat))
cells <- rownames_to_column(cells,"cells")
colnames(cells) <- c("cell","celltype")
meta <-  meta%>%merge(cells, by="cell")
# remove unassigned cells to match the numbers
meta <- meta %>% filter(Xa != "not_assigned") # 3685
count <- meta%>%dplyr::select(Xa,celltype)
count <- data.frame(table(count$Xa,count$celltype))
# get percentages
count <- count%>%pivot_wider(names_from = "Var2",values_from = "Freq")
count <- count%>%dplyr::filter(Var1!="not_assigned")
count$total <- rowSums(count[,2:8])
counts_per <- count
counts_per$`B cells` <- (counts_per$`B cells`/count$total)*100
counts_per$`T cells CD4` <- (counts_per$`T cells CD4`/count$total)*100
counts_per$`T cells CD8` <- (counts_per$`T cells CD8`/count$total)*100
counts_per$Monocytes <- (counts_per$Monocytes/count$total)*100
counts_per$`T cells` <- (counts_per$`T cells`/count$total)*100
counts_per$`Nk cells` <- (counts_per$`Nk cells`/count$total)*100
counts_per$Granulocytes <- (counts_per$Granulocytes/count$total)*100
plot <- counts_per%>%pivot_longer(cols = 2:8, names_to = "celltype")
plot <- plot%>%dplyr::filter(Var1 != c("wt_react"))
plot <- plot%>%dplyr::filter(Var1 != c("ko_react")) 
plot$celltype <- factor(plot$celltype, levels=c("Nk cells","Granulocytes","T cells","Monocytes","T cells CD8","T cells CD4","B cells"))  # define order 
plot$Var1 <- factor(plot$Var1, levels=c("cast_Xa_wt","cast_Xa","bl6_Xa_wt","bl6_Xa"))  
plot$Xa_status <- c(rep("BL6 Xa",14),rep("CAST Xa",14))
pdf(file=paste0(output,"celltype_composition_Xa_Xipdf"),width=4, height=7)
ggplot(plot, aes(x=Var1, y=value,fill=celltype,color="black")) + 
  geom_bar(stat="identity") +
  labs(x="",y="Percentage of Celltype") +
  scale_color_manual(values="black") +
  scale_fill_manual(values=rev(c("#F37A71","#C770AD","#00B9ED","#8784BF","#31B995","#BD9627","#56B947"))) +
  scale_y_continuous(limits = c(-0, 100.5), expand = c(0,0)) +
  theme_bw() +
  facet_wrap(~Xa_status, scales="free")
dev.off()
"#F37A71" # b cells 
"#BD9627" # cd4 
"#56B947" # cd8 
"#31B995" # t cells 
"#C770AD" # granu 
"#00B9ED" # mono 
"#8784BF" # nk 


count[1:2,]
######------  fisher test per cell type BL6 Xa------###### 
# B cells
fisher_mat_bcells_BL6_Xa <- data.frame(row.names = c(paste0("B cells BL6 Xa"),paste0("no B cells BL6 Xa")), 
                                KO=as.numeric(c(117,126)), # 243 - 117
                                WT = as.numeric(c(363,277))) # 640 - 363
fisher_res_bcells_BL6_Xa <- fisher.test(fisher_mat_bcells_BL6_Xa)
fisher_res_bcells_BL6_Xa <- data.frame(fisher_pvalue = as.numeric(unlist(fisher_res_bcells_BL6_Xa["p.value"])),
                                fisher_odds = as.numeric(unlist(fisher_res_bcells_BL6_Xa["estimate"])))
fisher_res_bcells_BL6_Xa$sample <- "B cells BL6 Xa"
# CD4
fisher_mat_cd4_BL6_Xa <- data.frame(row.names = c(paste0("CD4 BL6 Xa"),paste0("no CD4 BL6 Xa")), 
                             KO=as.numeric(c(49,194)), # 243 - 49
                             WT = as.numeric(c(115,525))) # 640 - 115
fisher_res_cd4_BL6_Xa <- fisher.test(fisher_mat_cd4_BL6_Xa)
fisher_res_cd4_BL6_Xa <- data.frame(fisher_pvalue = as.numeric(unlist(fisher_res_cd4_BL6_Xa["p.value"])),
                             fisher_odds = as.numeric(unlist(fisher_res_cd4_BL6_Xa["estimate"])))
fisher_res_cd4_BL6_Xa$sample <- "T cells CD4 BL6 Xa"
# CD8
fisher_mat_cd8_BL6_Xa <- data.frame(row.names = c(paste0("CD8 BL6 Xa"),paste0("no CD8 BL6 Xa")), 
                             KO=as.numeric(c(36,207)), # 243 - 36
                             WT = as.numeric(c(92,548))) # 640 - 92
fisher_res_cd8_BL6_Xa <- fisher.test(fisher_mat_cd8_BL6_Xa)
fisher_res_cd8_BL6_Xa <- data.frame(fisher_pvalue = as.numeric(unlist(fisher_res_cd8_BL6_Xa["p.value"])),
                             fisher_odds = as.numeric(unlist(fisher_res_cd8_BL6_Xa["estimate"])))
fisher_res_cd8_BL6_Xa$sample <- "T cells CD8 BL6 Xa"
# NK
fisher_mat_nk_BL6_Xa <- data.frame(row.names = c(paste0("Nk cells BL6 Xa"),paste0("no Nk cells BL6 Xa")), 
                            KO = as.numeric(c(8,235)), # 243 - 8
                            WT = as.numeric(c(13,627))) # 640 - 13
fisher_res_nk_BL6_Xa <- fisher.test(fisher_mat_nk_BL6_Xa)
fisher_res_nk_BL6_Xa <- data.frame(fisher_pvalue = as.numeric(unlist(fisher_res_nk_BL6_Xa["p.value"])),
                            fisher_odds = as.numeric(unlist(fisher_res_nk_BL6_Xa["estimate"])))
fisher_res_nk_BL6_Xa$sample <- "Nk cells BL6 Xa"
# Monocytes
fisher_mat_mono_BL6_Xa <- data.frame(row.names = c(paste0("Monocytes BL6 Xa"),paste0("no Monocytes BL6 Xa")), 
                            KO = as.numeric(c(11,232)), # 243 - 11
                            WT = as.numeric(c(25,615))) # 640 - 25
fisher_res_mono_BL6_Xa <- fisher.test(fisher_mat_mono_BL6_Xa)
fisher_res_mono_BL6_Xa <- data.frame(fisher_pvalue = as.numeric(unlist(fisher_res_mono_BL6_Xa["p.value"])),
                              fisher_odds = as.numeric(unlist(fisher_res_mono_BL6_Xa["estimate"])))
fisher_res_mono_BL6_Xa$sample <- "Monocytes BL6 Xa"
# Granulocytes
fisher_mat_granu_BL6_Xa <- data.frame(row.names = c(paste0("Granulocytes BL6 Xa"),paste0("no Granulocytes BL6 Xa")), 
                               KO=as.numeric(c(5,238)), # 243 - 5
                               WT = as.numeric(c(10,630))) # 640 - 10
fisher_res_granu_BL6_Xa <- fisher.test(fisher_mat_granu_BL6_Xa)
fisher_res_granu_BL6_Xa <- data.frame(fisher_pvalue = as.numeric(unlist(fisher_res_granu_BL6_Xa["p.value"])),
                               fisher_odds = as.numeric(unlist(fisher_res_granu_BL6_Xa["estimate"])))
fisher_res_granu_BL6_Xa$sample <- "Granulocytes BL6 Xa"
# T cells
fisher_mat_t_BL6_Xa <- data.frame(row.names = c(paste0("T cells BL6 Xa"),paste0("no T cells BL6 Xa")), 
                           KO=as.numeric(c(17,226)), # 243 - 17
                           WT = as.numeric(c(22,618))) # 640 - 22
fisher_res_t_BL6_Xa <- fisher.test(fisher_mat_t_BL6_Xa)
fisher_res_t_BL6_Xa <- data.frame(fisher_pvalue = as.numeric(unlist(fisher_res_t_BL6_Xa["p.value"])),
                           fisher_odds = as.numeric(unlist(fisher_res_t_BL6_Xa["estimate"])))
fisher_res_t_BL6_Xa$sample <- "T cells BL6 Xa"

count[3:4,]
######------  fisher test per cell type CAST Xa------###### 
# B cells
fisher_mat_bcells_CAST_Xa <- data.frame(row.names = c(paste0("B cells CAST Xa"),paste0("no B cells CAST Xa")), 
                                       KO=as.numeric(c(911,448)), # 1359 - 911
                                       WT = as.numeric(c(885,457))) # 1342 - 885
fisher_res_bcells_CAST_Xa <- fisher.test(fisher_mat_bcells_CAST_Xa)
fisher_res_bcells_CAST_Xa <- data.frame(fisher_pvalue = as.numeric(unlist(fisher_res_bcells_CAST_Xa["p.value"])),
                                       fisher_odds = as.numeric(unlist(fisher_res_bcells_CAST_Xa["estimate"])))
fisher_res_bcells_CAST_Xa$sample <- "B cells CAST Xa"
# CD4
fisher_mat_cd4_CAST_Xa <- data.frame(row.names = c(paste0("CD4 CAST Xa"),paste0("no CD4 CAST Xa")), 
                                    KO=as.numeric(c(143,1216)), # 1359 - 143
                                    WT = as.numeric(c(191,1151))) # 1342 - 191
fisher_res_cd4_CAST_Xa <- fisher.test(fisher_mat_cd4_CAST_Xa)
fisher_res_cd4_CAST_Xa <- data.frame(fisher_pvalue = as.numeric(unlist(fisher_res_cd4_CAST_Xa["p.value"])),
                                    fisher_odds = as.numeric(unlist(fisher_res_cd4_CAST_Xa["estimate"])))
fisher_res_cd4_CAST_Xa$sample <- "T cells CD4 CAST Xa"
# CD8
fisher_mat_cd8_CAST_Xa <- data.frame(row.names = c(paste0("CD8 CAST Xa"),paste0("no CD8 CAST Xa")), 
                                    KO=as.numeric(c(133,1226)), # 1359 - 133
                                    WT = as.numeric(c(121,1221))) # 1342 - 121
fisher_res_cd8_CAST_Xa <- fisher.test(fisher_mat_cd8_CAST_Xa)
fisher_res_cd8_CAST_Xa <- data.frame(fisher_pvalue = as.numeric(unlist(fisher_res_cd8_CAST_Xa["p.value"])),
                                    fisher_odds = as.numeric(unlist(fisher_res_cd8_CAST_Xa["estimate"])))
fisher_res_cd8_CAST_Xa$sample <- "T cells CD8 CAST Xa"
# NK
fisher_mat_nk_CAST_Xa <- data.frame(row.names = c(paste0("Nk cells CAST Xa"),paste0("no Nk cells CAST Xa")), 
                                   KO = as.numeric(c(39,1320)), # 1359 - 39
                                   WT = as.numeric(c(33,1309))) # 1342 - 33
fisher_res_nk_CAST_Xa <- fisher.test(fisher_mat_nk_CAST_Xa)
fisher_res_nk_CAST_Xa <- data.frame(fisher_pvalue = as.numeric(unlist(fisher_res_nk_CAST_Xa["p.value"])),
                                   fisher_odds = as.numeric(unlist(fisher_res_nk_CAST_Xa["estimate"])))
fisher_res_nk_CAST_Xa$sample <- "Nk cells CAST Xa"
# Monocytes
fisher_mat_mono_CAST_Xa <- data.frame(row.names = c(paste0("Monocytes CAST Xa"),paste0("no Monocytes CAST Xa")), 
                                     KO = as.numeric(c(60,1299)), # 1359 - 60
                                     WT = as.numeric(c(45,1297))) # 1342 - 45
fisher_res_mono_CAST_Xa <- fisher.test(fisher_mat_mono_CAST_Xa)
fisher_res_mono_CAST_Xa <- data.frame(fisher_pvalue = as.numeric(unlist(fisher_res_mono_CAST_Xa["p.value"])),
                                     fisher_odds = as.numeric(unlist(fisher_res_mono_CAST_Xa["estimate"])))
fisher_res_mono_CAST_Xa$sample <- "Monocytes CAST Xa"
# Granulocytes
fisher_mat_granu_CAST_Xa <- data.frame(row.names = c(paste0("Granulocytes CAST Xa"),paste0("no Granulocytes CAST Xa")), 
                                      KO=as.numeric(c(33,1326)), # 1359 - 33
                                      WT = as.numeric(c(32,1310))) # 1342 - 32
fisher_res_granu_CAST_Xa <- fisher.test(fisher_mat_granu_CAST_Xa)
fisher_res_granu_CAST_Xa <- data.frame(fisher_pvalue = as.numeric(unlist(fisher_res_granu_CAST_Xa["p.value"])),
                                      fisher_odds = as.numeric(unlist(fisher_res_granu_CAST_Xa["estimate"])))
fisher_res_granu_CAST_Xa$sample <- "Granulocytes CAST Xa"
# T cells
fisher_mat_t_CAST_Xa <- data.frame(row.names = c(paste0("T cells CAST Xa"),paste0("no T cells CAST Xa")), 
                                  KO=as.numeric(c(40,1319)), # 1359 - 40
                                  WT = as.numeric(c(35,1307))) # 1342 - 35
fisher_res_t_CAST_Xa <- fisher.test(fisher_mat_t_CAST_Xa)
fisher_res_t_CAST_Xa <- data.frame(fisher_pvalue = as.numeric(unlist(fisher_res_t_CAST_Xa["p.value"])),
                                  fisher_odds = as.numeric(unlist(fisher_res_t_CAST_Xa["estimate"])))
fisher_res_t_CAST_Xa$sample <- "T cells CAST Xa"


######------  plot fisher results BL6 Xa ------###### 
plot <- fisher_res_bcells_BL6_Xa%>%rbind(fisher_res_cd4_BL6_Xa,fisher_res_cd8_BL6_Xa,fisher_res_t_BL6_Xa,fisher_res_nk_BL6_Xa,fisher_res_mono_BL6_Xa,fisher_res_granu_BL6_Xa)
plot$log_p <- -log10(plot$fisher_pvalue)
# filter for total cells
plot <- plot %>% dplyr::filter(sample %in% c("B cells BL6 Xa","T cells CD4 BL6 Xa","T cells CD8 BL6 Xa"))
pdf(file=paste0(output,"/fisher_plot_BL6_Xa.pdf"),width=5, height=6)
ggplot(plot, aes(x = log_p, y = sample, size = 2,color=sample)) +
  geom_point() +
  scale_color_manual(values=rev(c("#F37A71","#C770AD","#00B9ED","#8784BF","#31B995","#BD9627","#56B947","#F37A71","#C770AD","#00B9ED","#8784BF","#31B995","#BD9627","#56B947"))) +
  geom_vline(xintercept=1.30103, linetype="dashed", color = "darkred") + 
  theme(axis.line = element_line(colour = 'black', size = 0.5),
        axis.text.y=element_text(size=10, colour = 'black'),
        axis.text.x=element_text(size=10, colour = 'black'),
        axis.title.y=element_text(size=12, angle = 90, colour = 'black'),
        axis.title.x=element_text(size=12, colour = 'black'),
        panel.background = element_blank(),
        panel.grid.minor =  element_line(colour = "grey",size=0.1),
        panel.grid.major.x = element_line(colour = "grey",size=0.1),
        panel.grid.major =  element_line(colour = "grey",size=0.1),
        legend.position = "none") +
  scale_x_continuous(limits=c(0,2.5),expand = c(0,0)) 
dev.off()


######------  plot fisher results CAST Xa ------###### 
plot <- fisher_res_bcells_CAST_Xa%>%rbind(fisher_res_cd4_CAST_Xa,fisher_res_cd8_CAST_Xa,fisher_res_t_CAST_Xa,fisher_res_nk_CAST_Xa,fisher_res_mono_CAST_Xa,fisher_res_granu_CAST_Xa)
plot$log_p <- -log10(plot$fisher_pvalue)
pdf(file=paste0(output,"/fisher_plot_CAST_Xa.pdf"),width=5, height=6)
ggplot(plot, aes(x = log_p, y = sample, size = 2,color=sample)) +
  geom_point() +
  scale_color_manual(values=rev(c("#F37A71","#C770AD","#00B9ED","#8784BF","#31B995","#BD9627","#56B947","#F37A71","#C770AD","#00B9ED","#8784BF","#31B995","#BD9627","#56B947"))) +
  geom_vline(xintercept=1.30103, linetype="dashed", color = "darkred") + 
  theme(axis.line = element_line(colour = 'black', size = 0.5),
        axis.text.y=element_text(size=10, colour = 'black'),
        axis.text.x=element_text(size=10, colour = 'black'),
        axis.title.y=element_text(size=12, angle = 90, colour = 'black'),
        axis.title.x=element_text(size=12, colour = 'black'),
        panel.background = element_blank(),
        panel.grid.minor =  element_line(colour = "grey",size=0.1),
        panel.grid.major.x = element_line(colour = "grey",size=0.1),
        panel.grid.major =  element_line(colour = "grey",size=0.1),
        legend.position = "none") +
  scale_x_continuous(limits=c(0,2.5),expand = c(0,0)) 
dev.off()
