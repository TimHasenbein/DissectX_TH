##*************************************##
##  Coexpression plot Firre/CrossFirre ##
##*************************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 09.2022
## Creation: 10.2021
## Create coexpression plot for Crossfirre, Firre and Dxz4 for different mouse tissues
library(tidyverse)
library(ggplot2)
library(data.table)


######------ Set environment  ------###### 
output <- "./"


######------ Get data and make count matrix  ------######
TPMs <- read.table("./TPM_table_adult_means.txt",header=T,check.names = FALSE)


######------ Filter genes/tissues of interest ------######
TPMs <- TPMs%>%filter(mgi_symbol=="Firre"|mgi_symbol=="Gm35612"|mgi_symbol=="4933407K13Rik")
TPMs <- TPMs%>%dplyr::select("mgi_symbol","Br_DFxDF_++","He_DFxDF_++","Ki_DFxDF_++","Li_DFxDF_++","Lu_DFxDF_++","Sp_DFxDF_++")


######------ plot scatter plot ------######
plot <- pivot_longer(TPMs, cols = 2:7, names_to = "sample",values_to = "TPMs")
plot$tissue <- rep(c("Brain","Heart","Kidney","Liver","Lung","Spleen"),3)
plot$log_TPMs <- log10(plot$TPMs+1)
pdf(file=paste0(output,"coexpression_plot.pdf"),width=7, height=4)
ggplot(plot,aes(x=tissue,y=log_TPMs,color=mgi_symbol,group=mgi_symbol)) +
  geom_line(color="grey",size=0.2) +
  geom_point(size=4) +
  labs(x="",y="TPM") +
  theme(axis.text.x=element_text(size=9,angle = 90),
        axis.ticks.x=element_blank(),
        axis.line = element_line(colour = 'black', size = 0.5),
        axis.ticks.y=element_line(colour = "black",size=0.5),
        axis.title.y=element_text(size=9),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.spacing = unit(1.0, 'cm')) +
  scale_y_continuous(limits = c(0,2), breaks = c(0,1,2),expand = c(0,0)) +
  scale_color_manual(values = c("#6B6766","#DB5D12","#B81918"), labels = c("Dxz4","Firre","crossFirre"))
dev.off()


