##************************************##
##  WT/LFD KO bar plot for TKO spleen  ##
##************************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 09.2022
## Creation: 11.2021
## Plot TPMs for samples for the KO genes in WT and mutant
library(tidyverse)
library(ggplot2)
library(biomaRt)
library(data.table)
library(stats)
library(pheatmap)
library(RColorBrewer)
library(GenomicFeatures)
library(openxlsx)


######------ Set environment and get data  ------###### 
input <- "./TPMs/"
output <- "./"


######------ Get TPM data  ------###### 
TPMs <- fread(paste0(input,"TPM_table_adult_replicates.txt"),header = T)


######------ Helper function  ------######
# for calculating the SD
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


######------------------------------------------------######
######        Spleen: XFL/WT + LF/WT + LFD/WT         ######
######------------------------------------------------######


######------ Filter spleen and group replicates  ------######
TPMs_sp <- TPMs%>%dplyr::select("mgi_symbol",starts_with("Sp_"))
TPMs_sp <- pivot_longer(TPMs_sp,cols=2:30,names_to="sample",values_to="TPMs")
TPMs_sp$sample <- gsub("Sp_","",TPMs_sp$sample)
TPMs_sp$sample <- gsub("_1","",TPMs_sp$sample)
TPMs_sp$sample <- gsub("_2","",TPMs_sp$sample)
TPMs_sp$sample <- gsub("_3","",TPMs_sp$sample)
TPMs_sp$sample <- gsub("_4","",TPMs_sp$sample)
TPMs_sp$sample <- gsub("XFLxXFL_\\+\\+","DFxDF_++",TPMs_sp$sample)


######------ Calculate mean TPMs  ------######
TPMs_sp <- TPMs_sp%>%filter(mgi_symbol=="Firre"|mgi_symbol=="Gm35612"|mgi_symbol=="4933407K13Rik")

# WT samples
sp_WT <- TPMs_sp%>%filter(sample=="DFxDF_++")
mean_sp_WT <- data_summary(sp_WT, varname="TPMs",groupnames=c("mgi_symbol"))
mean_sp_WT$mgi_symbol<- factor(mean_sp_WT$mgi_symbol, levels=c("Gm35612", "Firre", "4933407K13Rik"))
mean_sp_WT$GT <- "WT"

# TKO samples
sp_LFD <- TPMs_sp%>%filter(sample=="LFDxLFD_--")
mean_sp_LFD <- data_summary(sp_LFD, varname="TPMs", groupnames=c("mgi_symbol"))
mean_sp_LFD$mgi_symbol<- factor(mean_sp_LFD$mgi_symbol, levels=c("Gm35612", "Firre", "4933407K13Rik"))
mean_sp_LFD$GT <- "LFD"

# XFL samples
sp_XFL <- TPMs_sp%>%filter(sample=="XFLxXFL_--")
mean_sp_XFL <- data_summary(sp_XFL, varname="TPMs", groupnames=c("mgi_symbol"))
mean_sp_XFL$mgi_symbol<- factor(mean_sp_XFL$mgi_symbol, levels=c("Gm35612", "Firre", "4933407K13Rik"))
mean_sp_XFL$GT <- "XFL"

# LF samples
sp_LF <- TPMs_sp%>%filter(sample=="LFxLF_--")
mean_sp_LF <- data_summary(sp_LF, varname="TPMs", groupnames=c("mgi_symbol"))
mean_sp_LF$mgi_symbol<- factor(mean_sp_LF$mgi_symbol, levels=c("Gm35612", "Firre", "4933407K13Rik"))
mean_sp_LF$GT <- "LF"


######------ Plot TPMs for LFD, LF, XFL and WT  ------######
WT_LFD_sp <- mean_sp_WT%>%rbind(mean_sp_LFD)
WT_XFL_sp <- mean_sp_WT%>%rbind(mean_sp_XFL)
WT_LF_sp <- mean_sp_WT%>%rbind(mean_sp_LF)
pdf(file=paste0(output,"spleen_WT_LFD_TPMs.pdf"),paper="a4",width=25, height=6)
# TKO
ggplot(WT_LFD_sp, mapping = aes(x=mgi_symbol,y=TPMs, fill=mgi_symbol, ymin=TPMs-sd, ymax=TPMs+sd)) + 
  geom_bar(position = position_dodge2(reverse = T), stat = "identity", width = .9, color="black", size=0.2) +
  geom_errorbar(position=position_dodge2(reverse = T, width=.45), stat = "identity", colour="black") + 
  labs(x="",y="TPMs") +
  scale_y_continuous(expand = c(0, 0.1),limits = c(0., 20), breaks = c(0,5,10,15,10,20)) +
  scale_fill_manual(values=c("#B72227","#DA6627","#616262")) +
  theme(axis.text.x=element_text(size=9),
        axis.ticks.x=element_blank(),
        axis.line = element_line(colour = 'black', size = 0.2),
        axis.ticks.y=element_line(colour = "black",size=0.2),
        axis.title.y=element_text(size=9),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.spacing = unit(1.0, 'cm'))
dev.off()
# LF
pdf(file=paste0(output,"spleen_WT_LF_TPMs.pdf"),paper="a4",width=25, height=6)
ggplot(WT_LF_sp, mapping = aes(x=mgi_symbol,y=TPMs, fill=mgi_symbol, ymin=TPMs-sd, ymax=TPMs+sd)) + 
  geom_bar(position = position_dodge2(reverse = T), stat = "identity", width = .9, color="black", size=0.2) +
  geom_errorbar(position=position_dodge2(reverse = T, width=.45), stat = "identity", colour="black") + 
  labs(x="",y="TPMs") +
  scale_y_continuous(expand = c(0, 0.1),limits = c(0., 20), breaks = c(0,5,10,15,10,20)) +
  scale_fill_manual(values=c("#B72227","#DA6627","#616262")) +
  theme(axis.text.x=element_text(size=9),
        axis.ticks.x=element_blank(),
        axis.line = element_line(colour = 'black', size = 0.2),
        axis.ticks.y=element_line(colour = "black",size=0.2),
        axis.title.y=element_text(size=9),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.spacing = unit(1.0, 'cm'))
dev.off()

pdf(file=paste0(output,"spleen_WT_LF_TPMs.pdf"),paper="a4",width=25, height=6)
ggplot(WT_XFL_sp, mapping = aes(x=mgi_symbol,y=TPMs, fill=mgi_symbol, ymin=TPMs-sd, ymax=TPMs+sd)) + 
  geom_bar(position = position_dodge2(reverse = T), stat = "identity", width = .9, color="black", size=0.2) +
  geom_errorbar(position=position_dodge2(reverse = T, width=.45), stat = "identity", colour="black") + 
  labs(x="",y="TPMs") +
  scale_y_continuous(expand = c(0, 0.1),limits = c(0., 20), breaks = c(0,5,10,15,10,20)) +
  scale_fill_manual(values=c("#B72227","#DA6627","#616262")) +
  theme(axis.text.x=element_text(size=9),
        axis.ticks.x=element_blank(),
        axis.line = element_line(colour = 'black', size = 0.2),
        axis.ticks.y=element_line(colour = "black",size=0.2),
        axis.title.y=element_text(size=9),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.spacing = unit(1.0, 'cm'))
dev.off()

