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


######------ Get degs  ------######
input <- "./02_RNAseq_DESeq2/"
all_degs <- data.frame(DEG=as.numeric(),sample=as.character())
FDR=0.01
LFC=1
samples <- c("spleen_LFD","heart_LFD","brain_LFD","lung_LFD","liver_LFD","kidney_LFD")
for (sample in samples){
  split_sample=unlist(strsplit(sample,"_"))
  tissue=split_sample[1]
  if (length(split_sample)==3){
    GT=paste(split_sample[2],split_sample[3],sep="_")
  } else{
    GT=split_sample[2]
  }
  # import 
  data <- read.table(paste0(Sys.glob(paste0(input,tissue,"/",GT,"/",GT,"_shrunk_diffexpr-results*"))),
                     header=T)
  data <- data%>%dplyr::filter(abs(log2FoldChange) >= LFC & padj <= FDR)
  data <- data.frame(DEG=nrow(data),sample=paste0(tissue))
  all_degs <- all_degs%>%rbind(data)
}


######------ plot scatter plot ------######
plot <- TPMs%>%merge(all_degs, by = "sample")
#plot$tissue <- rep(c("Brain","Heart","Kidney","Liver","Lung","Spleen"),3)
pdf(file=paste0(output,"coexpression_plot_deg_tpms.pdf"),width=6, height=5)
ggplot(plot,aes(x=DEG,y=log_TPMs,color=gene,group=gene)) +
  geom_line(color="black",size=0.2) +
  geom_point(size=4) +
  #scale_shape_manual(values = c(15, 16)) + 
  labs(x="No. of DEGs",y="log10(TPM+1)") +
  theme(axis.line = element_line(colour = 'black', size = 0.5),
        axis.text.y=element_text(size=10, colour = 'black'),
        axis.text.x=element_text(size=10, colour = 'black'),
        axis.title.y=element_text(size=12, angle = 90, colour = 'black'),
        axis.title.x=element_text(size=12, colour = 'black'),
        panel.background = element_blank(),
        panel.grid.minor =  element_line(colour = "grey",size=0.1),
        panel.grid.major.x = element_line(colour = "grey",size=0.1),
        panel.grid.major =  element_line(colour = "grey",size=0.1)) + 
  scale_y_continuous(limits = c(0,2), breaks = c(0,0.5,1,1.5,2),expand = c(0,0)) +
  scale_color_manual(values = c("#6C6867","#DB5E16","#B81D1B")) +
  geom_abline(intercept = coef(fit_firre)[1], slope = coef(fit_firre)[2], color = "black") +
  geom_abline(intercept = coef(fit_crossfirre)[1], slope = coef(fit_crossfirre)[2], color = "black") +
  geom_abline(intercept = coef(fit_dxz4)[1], slope = coef(fit_dxz4)[2], color = "black") 
dev.off()


######------ Correlation btw DEGs and TPMs ------######
plot <- plot%>%dplyr::select(-TPM)
cor <- plot%>%pivot_wider(names_from="gene",values_from = c("log_TPMs","DEG"))
shapiro.test(cor$DEG_Firre) 
# Firre
cor_firre <- cor.test(cor$log_TPMs_Firre, cor$DEG_Firre, method = "pearson")
fit_firre <- lm(log_TPMs_Firre ~ DEG_Firre, data = cor)
# Crossfirre
cor_crossfirre <- cor.test(cor$log_TPMs_Gm35612, cor$DEG_Gm35612, method = "pearson")
fit_crossfirre <-  lm(log_TPMs_Gm35612 ~ DEG_Gm35612, data = cor)
# Dxz4
cor_dxz4 <- cor.test(cor$log_TPMs_4933407K13Rik, cor$DEG_4933407K13Rik, method = "pearson")
fit_dxz4 <- lm(log_TPMs_4933407K13Rik ~ DEG_4933407K13Rik, data = cor)


