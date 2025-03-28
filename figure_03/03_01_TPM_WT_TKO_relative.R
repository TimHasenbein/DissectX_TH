##*************************************##
##  WT/LFD KO bar plot for TKO spleen  ##
##*************************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 09.2022
## Creation: 11.2021
## Spleen: TKO
library(tidyverse)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(stats)


######------ Set environment and get data  ------###### 
input <- "./TPMs/"
output <- "./"


######------ Get TPM data  ------###### 
TPMs <- fread(paste0(input,"TPM_table_placenta_replicates.txt"),header = T)


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


######------ Filter Xa/Xi and group replicates  ------######
# Xi
TPMs_xi <- TPMs%>%dplyr::select("mgi_symbol","Pl_CxD_++_1","Pl_CxD_++_2","Pl_CxD_++_3","Pl_CxDF_++_1","Pl_CxDF_++_2","Pl_CxDF_++_3","Pl_CxFi_++_1","Pl_CxFi_++_2","Pl_CxFi_++_3",starts_with("Pl_CxLFD"))
TPMs_xi <- pivot_longer(TPMs_xi,cols=2:13,names_to="sample",values_to="TPMs")
TPMs_xi$sample <- gsub("Pl_","",TPMs_xi$sample)
TPMs_xi$sample <- gsub("_1","",TPMs_xi$sample)
TPMs_xi$sample <- gsub("_2","",TPMs_xi$sample)
TPMs_xi$sample <- gsub("_3","",TPMs_xi$sample)
TPMs_xi$sample <- gsub("CxD_\\++","WT",TPMs_xi$sample)
TPMs_xi$sample <- gsub("CxDF_\\++","WT",TPMs_xi$sample)
TPMs_xi$sample <- gsub("CxFi_\\++","WT",TPMs_xi$sample)
# Xa
TPMs_xa <- TPMs%>%dplyr::select("mgi_symbol","Pl_DxC_++_1","Pl_DxC_++_2","Pl_DFxC_++_1","Pl_DFxC_++_2","Pl_DFxC_++_3","Pl_FixC_++_1","Pl_FixC_++_2","Pl_FixC_++_3",starts_with("Pl_LFDxC"))
TPMs_xa <- pivot_longer(TPMs_xa,cols=2:12,names_to="sample",values_to="TPMs")
TPMs_xa$sample <- gsub("Pl_","",TPMs_xa$sample)
TPMs_xa$sample <- gsub("_1","",TPMs_xa$sample)
TPMs_xa$sample <- gsub("_2","",TPMs_xa$sample)
TPMs_xa$sample <- gsub("_3","",TPMs_xa$sample)
TPMs_xa$sample <- gsub("DxC_\\++","WT",TPMs_xa$sample)
TPMs_xa$sample <- gsub("DFxC_\\++","WT",TPMs_xa$sample)
TPMs_xa$sample <- gsub("FixC_\\++","WT",TPMs_xa$sample)


######------ Calculate mean TPMs  ------######
TPMs_xi <- TPMs_xi%>%filter(mgi_symbol=="Firre"|mgi_symbol=="Gm35612"|mgi_symbol=="4933407K13Rik")
# WT samples
xi_WT <- TPMs_xi%>%filter(sample=="WT")
mean_xi_WT <- data_summary(xi_WT, varname="TPMs",groupnames=c("mgi_symbol"))
mean_xi_WT$mgi_symbol<- factor(mean_xi_WT$mgi_symbol, levels=c("Gm35612", "Firre", "4933407K13Rik"))
mean_xi_WT$GT <- "WT"
# TKO samples
xi_LFD <- TPMs_xi%>%filter(sample=="CxLFD_+-")
mean_xi_LFD <- data_summary(xi_LFD, varname="TPMs", groupnames=c("mgi_symbol"))
mean_xi_LFD$mgi_symbol<- factor(mean_xi_LFD$mgi_symbol, levels=c("Gm35612", "Firre", "4933407K13Rik"))
mean_xi_LFD$GT <- "LFD"

TPMs_xa <- TPMs_xa%>%filter(mgi_symbol=="Firre"|mgi_symbol=="Gm35612"|mgi_symbol=="4933407K13Rik")
# WT samples
xa_WT <- TPMs_xa%>%filter(sample=="WT")
mean_xa_WT <- data_summary(xa_WT, varname="TPMs",groupnames=c("mgi_symbol"))
mean_xa_WT$mgi_symbol<- factor(mean_xa_WT$mgi_symbol, levels=c("Gm35612", "Firre", "4933407K13Rik"))
mean_xa_WT$GT <- "WT"
# TKO samples
xa_LFD <- TPMs_xa%>%filter(sample=="LFDxC_-+")
mean_xa_LFD <- data_summary(xa_LFD, varname="TPMs", groupnames=c("mgi_symbol"))
mean_xa_LFD$mgi_symbol<- factor(mean_xa_LFD$mgi_symbol, levels=c("Gm35612", "Firre", "4933407K13Rik"))
mean_xa_LFD$GT <- "LFD"


######------ Calculate relative expression  ------######
# Xa
mean_xa_WT$rel_exp <- mean_xa_WT$TPMs/mean_xa_WT$TPMs
mean_xa_LFD$rel_exp <- mean_xa_LFD$TPMs/mean_xa_WT$TPMs
mean_xa_WT$rel_sd <- mean_xa_WT$sd/mean_xa_WT$TPMs
mean_xa_LFD$rel_sd <- mean_xa_LFD$sd/mean_xa_WT$TPMs
WT_LFD_xa <- mean_xa_WT%>%rbind(mean_xa_LFD)
# Xi
mean_xi_WT$rel_exp <- mean_xi_WT$TPMs/mean_xi_WT$TPMs
mean_xi_LFD$rel_exp <- mean_xi_LFD$TPMs/mean_xi_WT$TPMs
mean_xi_WT$rel_sd <- mean_xi_WT$sd/mean_xi_WT$TPMs
mean_xi_LFD$rel_sd <- mean_xi_LFD$sd/mean_xi_WT$TPMs
WT_LFD_xi <- mean_xi_WT%>%rbind(mean_xi_LFD)


######------ Plot TPMs for LFD and WT  ------######
pdf(file=paste0(output,"/TPMs_WT_LFD_xa_relativeexpr.pdf"),paper="a4",width=25, height=6)
ggplot(WT_LFD_xa, mapping = aes(x=mgi_symbol,y=rel_exp, fill=GT, 
                                ymin=rel_exp-rel_sd, ymax=rel_exp+rel_sd)) + 
  geom_bar(position = position_dodge2(reverse = T), stat = "identity", width = .9, 
           color="black", size=0.2) +
  geom_errorbar(position=position_dodge2(reverse = T, width=.45), stat = "identity", colour="black") + 
  labs(x="",y="Relative expression (WT set to 1)") +
  scale_y_continuous(expand = c(0, 0.1),limits = c(0., 2), breaks = c(0,0.5,1,1.5,2)) +
  scale_fill_manual(values=c("#9FD0B8","black")) +
  theme(axis.text.x=element_text(size=9),
        axis.ticks.x=element_blank(),
        axis.line = element_line(colour = 'black', size = 0.5),
        axis.ticks.y=element_line(colour = "black",size=0.5),
        axis.title.y=element_text(size=9),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.spacing = unit(1.0, 'cm'))
dev.off()

# Xi
pdf(file=paste0(output,"/TPMs_WT_LFD_xi_relativeexpr.pdf"),paper="a4",width=25, height=6)
ggplot(WT_LFD_xi, mapping = aes(x=mgi_symbol,y=rel_exp, fill=GT, 
                                ymin=rel_exp-rel_sd, ymax=rel_exp+rel_sd)) + 
  geom_bar(position = position_dodge2(reverse = T), stat = "identity", width = .9, 
           color="black", size=0.2) +
  geom_errorbar(position=position_dodge2(reverse = T, width=.45), stat = "identity", colour="black") + 
  labs(x="",y="Relative expression (WT set to 1)") +
  scale_y_continuous(expand = c(0, 0.1),limits = c(0., 2.5), breaks = c(0,0.5,1,1.5,2,2.5)) +
  scale_fill_manual(values=c("#9FD0B8","black")) +
  theme(axis.text.x=element_text(size=9),
        axis.ticks.x=element_blank(),
        axis.line = element_line(colour = 'black', size = 0.5),
        axis.ticks.y=element_line(colour = "black",size=0.5),
        axis.title.y=element_text(size=9),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.spacing = unit(1.0, 'cm'))
dev.off()

