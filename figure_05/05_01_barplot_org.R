##***********************************************##
##   Figure 5 Number of degs for across organs   ##
##***********************************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 09.2022
## Creation: 09.2021
## Show number of DEG FDR <= 0.01, |log2FC| <= 1 
library(ggplot2)
library(tidyverse)
library(data.table)
FDR=0.01
LFC=1


######------ Set environment  ------######
input <- "./02_RNAseq_DESeq2/"
output <- "./number_degs/"


######------ get data for LFD and DF for each tissue  ------######
samples <- c("spleen_LFD","heart_LFD","brain_LFD","lung_LFD","liver_LFD","kidney_LFD",
             "spleen_DF","heart_DF","brain_DF","lung_DF","liver_DF","kidney_DF")
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
  # write degs for GT
  write.table(data, paste0(output,GT,"_",tissue,"_FDR",FDR,"_LFC",LFC,".txt"),
              row.names = F,col.names = T,quote = F,sep="\t")
  # get numbers of up and downregulated degs and write out
  data <- data%>%dplyr::select("log2FoldChange")
  up <- sum(data>0); down <- sum(data<0); up_down <- data.frame(tissue_name=tissue,up=up,down=down)
  assign(paste0(sample), up_down)
  write.table(up_down, paste0(output,GT,"_",tissue,"_FDR",FDR,"_LFC",LFC,"_direction.txt"),
              row.names = F,col.names = T,quote = F,sep="\t")  
}


######------ make matrix to count ------######
# LFD
LFD <- heart_LFD%>%rbind(lung_LFD)%>%rbind(liver_LFD)%>%
  rbind(spleen_LFD)%>%rbind(kidney_LFD)%>%rbind(brain_LFD)
LFD$tissue <- as.factor(c(1,2,3,4,5,6))
LFD$tissue <- factor(LFD$tissue, levels=c("4", "5","2", "1", "3","6"))         
LFD$GT <- "LFD"
LFD$total <- LFD$up + LFD$down

# DF
DF <- heart_DF%>%rbind(lung_DF)%>%rbind(liver_DF)%>%
  rbind(spleen_DF)%>%rbind(kidney_DF)%>%rbind(brain_DF)
DF$tissue <- as.factor(c(1,2,3,4,5,6))
DF$tissue <- factor(DF$tissue, levels=c("4", "5","2", "1", "3","6"))  
DF$GT <- "DF"
DF$total <- DF$up + DF$down


######------ plot stacked chart bar for DF and LFD ------######
pdf(file=paste0(output,"number_degs.pdf"),paper="a4",width=6, height=5)
ggplot() + 
   geom_bar(LFD, mapping = aes(x=as.numeric(tissue),y=total, fill=GT),
           position = "stack",stat = "identity", width = 0.35,linetype=1, color="black", size=0.2) +
   geom_bar(DF, mapping = aes(x=as.numeric(tissue)+ 0.35 + 0.02,y=total, fill=GT),
           position = "stack",stat = "identity", width = 0.35,linetype=1, color="black", size=0.2) +
   scale_x_discrete(limits = c("Spleen","Kidney","Lung","Heart","Liver","Brain")) +
  labs(x="",y="Number of DE genes") +
  theme(axis.text.x=element_text(size=9),
        axis.ticks.x=element_blank(),
        axis.line = element_line(colour = 'black', size = 0.5),
        axis.ticks.y=element_line(colour = "black",size=0.5),
        axis.title.y=element_text(size=9),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.spacing = unit(1.0, 'cm')) +
  scale_fill_manual(values = c("#0080ff","#92CDB5"),
                    labels = c("DF","LFD")) +
  scale_y_continuous(expand = c(0, 0),position = "left") +
  labs(fill = "Direction")
dev.off()
write.table(LFD, paste0(output,"overview_LFD.txt"),
                      row.names = F,col.names = T,quote = F,sep="\t")  
write.table(DF, paste0(output,"overview_DF.txt"),
            row.names = F,col.names = T,quote = F,sep="\t")  


