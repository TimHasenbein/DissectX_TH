##******************************************************##
##   Figure 5 Dysregulated genes from autosomes/X-chr   ##
##******************************************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 09.2022
## Creation: 01.2022
library(ggplot2)
library(data.table)
library(tidyverse)
library(xlsx)
library(biomaRt)
FDR=0.01
LFC=1
GT="LFD"


######------  set environment  ------######
input <- "./02_RNAseq_DESeq2/"
output <- "./pie_charts/"
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="http://aug2020.archive.ensembl.org")


######------  get TKO data, fisher test and plot  ------######
samples <- c("spleen","heart","brain","lung","liver","kidney")
for(sample in samples) {
  # import data and get DEGs
  data <- read.table(Sys.glob(paste0(input,sample,"/",GT,"/",GT,"_shrunk_diffexpr-results*")),header=T)
  data <- data%>%dplyr::select("name","log2FoldChange","padj")
  # remove KO genes
  data <- data%>%filter(name != "Firre" & name != "4933407K13Rik" & name != "Gm35612")
  # get chr info -> does not consider Gm35612 but is not in the data either way!!
  IDs <- getBM(attributes = c("mgi_symbol","chromosome_name"),
               filters = "mgi_symbol", values = data[,1],
               mart = mart)
  data  <- data%>%inner_join(IDs, by=c("name"="mgi_symbol"))
  data <- data[,c("name","log2FoldChange","padj","chromosome_name")]
  # remove genes from patches
  data <- data%>%filter(!grepl('CHR_|JH|GL|MT', chromosome_name))%>%distinct
  # assign chromosome status
  data$chr_status <- NA
  for(i in 1:nrow(data)) {
    ifelse(data$chromosome_name[i] != "X",
           data$chr_status[i] <- "autosome",
           data$chr_status[i] <- "chrX")}
  # get number of degs
  degs <- data %>% dplyr::filter(abs(log2FoldChange) >= LFC & padj <= FDR) 
  # get number of non degs
  non_degs <- data %>% filter(!name %in% degs$rn)
  # make fisher matrix and fisher test for degs/non_degs for autosome/X-chr
  fisher_mat <- data.frame(row.names = c(paste0(sample,"_autosome"),paste0(sample,"_chrX")), 
                     degs=as.numeric(c(sum(degs$chr_status=="autosome"),
                                   sum(degs$chr_status=="chrX"))),
                     non_degs = as.numeric(c(sum(non_degs$chr_status=="autosome"),
                                  sum(non_degs$chr_status=="chrX"))))
  fisher_res <- fisher.test(fisher_mat,alternative = "less")
  fisher_res <- data.frame(row.names = paste0(sample),
                           fisher_pvalue = as.numeric(unlist(fisher_res["p.value"])),
                           fisher_odds = as.numeric(unlist(fisher_res["estimate"])))
  # capture fisher result output
  assign(paste0("fisher_mat_",sample), fisher_mat)
  assign(paste0("fisher_res_",sample), fisher_res)
  # generate plot matrix
  plot <- fisher_mat
  plot$total <- sum(plot$degs)
  plot$percentage <- plot$degs/plot$total*100
  plot <- setDT(plot, keep.rownames = TRUE)[]
  assign(paste0("plot_",sample), plot)
}  


######------  make plot matrix  ------###### 
plot <- plot_brain%>%rbind(plot_heart)%>%rbind(plot_liver)%>%rbind(plot_lung)%>%rbind(plot_kidney)%>%rbind(plot_spleen)
plot <- plot%>%separate(rn,c("organ", "chr"), "_") 
plot$side1 <- c("A","A","A","A","A","A","B","B","B","B","B","B")
plot$side2 <- c("A","A","B","B","C","C","A","A","B","B","C","C")
write.table(plot, paste0(output,"plot_overview.txt"), sep = "\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
# get mean % across tissues
auto <- plot%>%filter(chr=="autosome"); auto_mean <- mean(auto$percentage) 
X <- plot%>%filter(chr=="chrX"); X_mean <- mean(X$percentage) 


######------  plot pie chart  ------###### 
pdf(file=paste0(output,"all_piechart_DGE_autoX_",GT,"_KOremoved_shrunk.pdf"), paper="a4r",width=10, height=10)
print(ggplot(plot, aes(x = total*10, y = degs, fill = chr, width=total*20)) + 
          geom_bar(position ="fill",stat = "identity", linetype=1, color = "black") + 
          labs(x="",y="% DE genes") +
          facet_grid(side1 ~ side2) +
          theme(axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.line = element_blank(),
                axis.line.x = element_blank(),
                axis.ticks.y=element_blank(),
                axis.title.y=element_blank(), 
                axis.text.y=element_blank(),
                panel.background = element_rect(fill = 'white', colour = 'white'),
                panel.grid  = element_blank()) +
          scale_fill_manual(NULL,values = c(paste0("#92CDB5"),"black"), labels = c("autosomes","X-chromosome")) +
          coord_polar("y")) 
dev.off()


######------  output fisher test results  ------###### 
fisher_res_all <- fisher_res_brain%>%rbind(fisher_res_heart)%>%rbind(fisher_res_liver)%>%rbind(fisher_res_kidney)%>%rbind(fisher_res_lung)%>%rbind(fisher_res_spleen)
fisher_mat_all <- fisher_mat_brain%>%rbind(fisher_mat_heart)%>%rbind(fisher_mat_liver)%>%rbind(fisher_mat_kidney)%>%rbind(fisher_mat_lung)%>%rbind(fisher_mat_spleen)
write.xlsx(fisher_mat_all, file=paste0(output,"fisher_test_TKO_less.xlsx"),sheetName="Fisher's exact test matrix")
write.xlsx(fisher_res_all, file=paste0(output,"fisher_test_TKO_less.xlsx"),sheetName="Fisher's exact test results",append=T)

