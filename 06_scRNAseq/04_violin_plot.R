##*****************************##
##  Allelome.PRO Violin plots  ##
##*****************************##
## Project: dissectX
## Tim Hasenbein 
## Last modification 04.2022
## Creation: 04.2022 
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


######------  set environment  ------###### 
input <- "./"
output <- "./"


######------  get sorted cells cells  ------###### 
bl6_ko <- read.table(paste0(input,"bl6_+-_Xa.txt"),header=T) 
cast_ko <- read.table(paste0(input,"cast_+-_Xa.txt"),header=T)
react_ko <- read.table(paste0(input,"react_+-.txt"),header=T) 
bl6_wt <- read.table(paste0(input,"bl6_++_Xa.txt"),header=T) 
cast_wt <- read.table(paste0(input,"cast_++_Xa.txt"),header=T) 
react_wt <- read.table(paste0(input,"react_++.txt"),header=T) 
ap_x <- bl6_ko%>%rbind(cast_ko)%>%rbind(react_ko)%>%rbind(bl6_wt)%>%rbind(cast_wt)%>%rbind(react_wt)


######------   violin plot   ------###### 
# prepare matrix
ap_x$GT <- factor(ap_x$GT, levels = c("WT", "KO")) 
# define dodge in order that barplot and violin plot overlap
dodge <- position_dodge(width = 1)
pdf(file=paste0(output,"04_Xi_cast_vs_bl6.pdf"),width=5, height=5)
ggplot(ap_x, aes(x=GT, y=allelic_ratio,fill=GT)) + 
  geom_violin(width = 0.5, size = 0.1, scale = "width", position = dodge) +
  labs(x="Genotype",y="Allelic ratio") +
  geom_boxplot(width = 0.3, size=0.2, alpha = 0.2, color = "white", outlier.shape=NA,position = dodge) +
  scale_fill_manual(values = c("#232b2b","#90CAB2")) +
  scale_x_discrete(limits=c("WT", "KO"), guide = guide_axis(angle = 0), position = "bottom") +
  scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0,0.3,0.5,0.7,1)) +
  theme(
    axis.ticks.x=element_blank(),
    axis.title.x=element_blank(),
    axis.text.x = element_text(angle = 45, hjust=0.1),
    axis.line = element_line(colour = 'darkgrey', size = 0.2),
    axis.ticks.y=element_line(colour = "grey40",size=0.2),
    axis.title.y=element_text(size=12),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.border = element_rect(colour = "darkgrey", fill=NA, size=0.5),
    legend.position="none") +
  geom_hline(yintercept=0.7, linetype="dashed", color = "grey") + 
  geom_hline(yintercept=0.3, linetype="dashed", color = "grey") +
  annotate("text", x="WT", y=-0.02, label= paste0("n=",nrow(bl6_wt)," (",round((nrow(bl6_wt)/(nrow(bl6_wt)+nrow(cast_wt)+nrow(react_wt)))*100,digits = 1),"%)")) +
  annotate("text", x="WT", y=1.02, label= paste0("n=",nrow(cast_wt)," (",round((nrow(cast_wt)/(nrow(bl6_wt)+nrow(cast_wt)+nrow(react_wt)))*100,digits = 1),"%)")) +
  annotate("text", x="WT", y=0.5, label= paste0("n=",nrow(react_wt)," (",round((nrow(react_wt)/(nrow(bl6_wt)+nrow(cast_wt)+nrow(react_wt)))*100,digits = 1),"%)")) +
  annotate("text", x="KO", y=-0.02, label= paste0("n=",nrow(bl6_ko)," (",round((nrow(bl6_ko)/(nrow(bl6_ko)+nrow(cast_ko)+nrow(react_ko)))*100,digits = 1),"%)")) +
  annotate("text", x="KO", y=1.02, label= paste0("n=",nrow(cast_ko)," (",round((nrow(cast_ko)/(nrow(bl6_ko)+nrow(cast_ko)+nrow(react_ko)))*100,digits = 1),"%)")) +
  annotate("text", x="KO", y=0.5, label= paste0("n=",nrow(react_ko)," (",round((nrow(react_ko)/(nrow(bl6_ko)+nrow(cast_ko)+nrow(react_ko)))*100,digits = 1),"%)")) 
dev.off()







