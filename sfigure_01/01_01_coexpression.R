##****************##
##  Coexpression  ##
##****************##
## Project: dissectX
## Tim Hasenbein
## Last modification 04.2024
## Creation: 04.2024
library(ggplot2)
library(data.table)
library(biomaRt)
library(GenomicFeatures)
library(openxlsx)
library(RColorBrewer)
library(pheatmap)


######------ Set environment  ------###### 
output <- "./"


######------ Get data and make count matrix  ------######
TPMs <- read.table("./TPM_table_adult_means.txt",header=T,check.names = FALSE)
TPMs <- read.table("./TPM_table_adult_replicates.txt",header=T,check.names = FALSE)

######------ Filter genes/tissues of interest ------######
TPMs <- TPMs%>%filter(mgi_symbol=="Firre"|mgi_symbol=="Gm35612"|mgi_symbol=="4933407K13Rik")
TPMs <- TPMs%>%dplyr::select("mgi_symbol","Br_DFxDF_++","He_DFxDF_++","Ki_DFxDF_++","Li_DFxDF_++","Lu_DFxDF_++","Sp_DFxDF_++")
colnames(TPMs) <- c("gene","brain","heart","kidney","liver","lung","spleen")
TPMs <- TPMs%>%pivot_longer(!gene, names_to="sample",values_to="TPM")
TPMs$log_TPMs <- log10(TPMs$TPM+1)


######------ Check distribution ------####
cor <- TPMs%>%pivot_wider(names_from="gene",values_from = c("TPM","log_TPMs"))
shapiro.test(cor$TPM_Firre) 
shapiro.test(cor$TPM_Gm35612) 
shapiro.test(cor$TPM_4933407K13Rik) 
shapiro.test(cor$log_TPMs_Firre) 
shapiro.test(cor$log_TPMs_4933407K13Rik) 
shapiro.test(cor$log_TPMs_Gm35612) 


######------ Correlation plot Crossfirre-Firre ------######
cor <- TPMs%>%filter(gene != "4933407K13Rik")
cor <- cor%>%pivot_wider(names_from="gene",values_from = c("TPM","log_TPMs"))
# raw TPMs
cor_test <- cor.test(cor$TPM_Firre, cor$TPM_Gm35612, method = "pearson")
fit <- lm(TPM_Firre ~ TPM_Gm35612, data = cor)
pdf(file=paste0(output,"/cor_TPM_CF.pdf"),width=6, height=5)
ggplot(cor, aes(x = TPM_Gm35612, y = TPM_Firre, color = sample)) +
  geom_point(size = 3) +
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], color = "black") +
  scale_color_manual(values = c("#FAA51C","black","#6B9E9C","#8B191B","#ABABAB","#3389B0")) +
  theme(axis.line = element_line(colour = 'black', size = 0.5),
        axis.text.y=element_text(size=10, colour = 'black'),
        axis.text.x=element_text(size=10, colour = 'black'),
        axis.title.y=element_text(size=12, angle = 90, colour = 'black'),
        axis.title.x=element_text(size=12, colour = 'black'),
        panel.background = element_blank(),
        panel.grid.minor =  element_line(colour = "grey",size=0.1),
        panel.grid.major.x = element_line(colour = "grey",size=0.1),
        panel.grid.major =  element_line(colour = "grey",size=0.1)) +
  annotate("text",x=10,y=70,label=paste0("R = ",round(cor_test$estimate,2),paste0(" p = 0.4893"))) +
  scale_y_continuous(limits=c(0,80),expand = c(0,0),breaks = c(20,40,60,80)) +
  scale_x_continuous(limits=c(0,12.5),expand = c(0,0)) +
  ggtitle("Pearson correlation Crossfirre-Firre raw TPMs")
dev.off()


######------ Correlation plot Crossfirre-Dxz4 ------######
cor <- TPMs%>%filter(gene != "Firre")
cor <- cor%>%pivot_wider(names_from="gene",values_from = c("TPM","log_TPMs"))
# raw TPMs
cor_test <- cor.test(cor$TPM_4933407K13Rik, cor$TPM_Gm35612, method = "pearson")
fit <- lm(TPM_4933407K13Rik ~ TPM_Gm35612, data = cor)
pdf(file=paste0(output,"/cor_TPM_CD.pdf"),width=6, height=5)
ggplot(cor, aes(x = TPM_Gm35612, y = TPM_4933407K13Rik, color = sample)) +
  geom_point(size = 3) +
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], color = "black") +
  scale_color_manual(values = c("#FAA51C","black","#6B9E9C","#8B191B","#ABABAB","#3389B0")) +
  theme(axis.line = element_line(colour = 'black', size = 0.5),
        axis.text.y=element_text(size=10, colour = 'black'),
        axis.text.x=element_text(size=10, colour = 'black'),
        axis.title.y=element_text(size=12, angle = 90, colour = 'black'),
        axis.title.x=element_text(size=12, colour = 'black'),
        panel.background = element_blank(),
        panel.grid.minor =  element_line(colour = "grey",size=0.1),
        panel.grid.major.x = element_line(colour = "grey",size=0.1),
        panel.grid.major =  element_line(colour = "grey",size=0.1)) +
  annotate("text",x=7.5,y=4,label=paste0("R = ",round(cor_test$estimate,2),paste0(" p = 0.5273"))) +
  scale_y_continuous(limits=c(0,4.5),expand = c(0,0),) +
  scale_x_continuous(limits=c(0,12.5),expand = c(0,0)) +
  ggtitle("Pearson correlation Crossfirre-Dxz4 raw TPMs")
dev.off()


######------ Correlation plot Firre-Dxz4 ------###### 
cor <- TPMs%>%filter(gene != "Gm35612")
cor <- cor%>%pivot_wider(names_from="gene",values_from = c("TPM","log_TPMs"))
# raw TPMs
cor_test <- cor.test(cor$TPM_4933407K13Rik, cor$TPM_Firre, method = "pearson")
fit <- lm(TPM_4933407K13Rik ~ TPM_Firre, data = cor)
pdf(file=paste0(output,"/cor_TPM_FD.pdf"),width=6, height=5)
ggplot(cor, aes(x = TPM_Firre, y = TPM_4933407K13Rik, color = sample)) +
  geom_point(size = 3) +
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], color = "black") +
  scale_color_manual(values = c("#FAA51C","black","#6B9E9C","#8B191B","#ABABAB","#3389B0")) +
  theme(axis.line = element_line(colour = 'black', size = 0.5),
        axis.text.y=element_text(size=10, colour = 'black'),
        axis.text.x=element_text(size=10, colour = 'black'),
        axis.title.y=element_text(size=12, angle = 90, colour = 'black'),
        axis.title.x=element_text(size=12, colour = 'black'),
        panel.background = element_blank(),
        panel.grid.minor =  element_line(colour = "grey",size=0.1),
        panel.grid.major.x = element_line(colour = "grey",size=0.1),
        panel.grid.major =  element_line(colour = "grey",size=0.1)) +
  annotate("text",x=20,y=4,label=paste0("R = ",round(cor_test$estimate,2),paste0(" p = 0.2419"))) +
  scale_x_continuous(limits=c(0,80),expand = c(0,0),breaks = c(20,40,60,80)) +
  scale_y_continuous(limits=c(0,4.5),expand = c(0,0)) +
  ggtitle("Pearson correlation Firre-Dxz4 raw TPMs")
dev.off()
