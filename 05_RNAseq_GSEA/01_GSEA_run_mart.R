##********************************##
##    GSE Anylsis using camera    ##
##********************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 09.2022
## Creation: 07.2021
## Perform GSEA for tissues and genotypes
library(tidyverse)
library(ggplot2)
library(data.table)
library(limma)
library(simplifyEnrichment)


######------ Set environment  ------###### 
input_folder <- "./"
output <- "./"
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="http://aug2020.archive.ensembl.org") 


######------ Gene set file  ------###### 
GO_gmt <- paste0(output,"c5.go.v7.4.symbols.gmt")
GO_gmt <- read.table(GO_gmt, sep="\n", header=F, stringsAsFactors = F)
GO_gmt <- lapply(GO_gmt[,1], function(x) {strsplit(x, split="\t")[[1]]})
names(GO_gmt) <- sapply(GO_gmt, function(x) {x[1]})
GO_gmt <- lapply(GO_gmt, function(x) {x[-c(1,2)]})
saveRDS(GO_gmt,paste0(output,"go_gmt_c5v7.4.rds"))


######------ Function for single gene set  ------###### 
runCamera <- function(geneset) {
  out <- cameraPR(input$stat, toupper(input$name) %in% geneset)
  colnames(out) <- c("numGenes", "direction", "pval")
  return(out)
}


######-------------------------------------------------------------######
##                             Perform GSEA                            ##
######-------------------------------------------------------------######
samples <- c("spleen_LFD","spleen_LF")


######------ Loop through all samples ------###### 
for (sample in samples){
  FDR=0.1
  split_sample=unlist(strsplit(sample,"_"))
  tissue=split_sample[1]
  if (length(split_sample)==3){
    GT=paste(split_sample[2],split_sample[3],sep="_")
  } else{
    GT=split_sample[2]
  }
  
  # import and change EnsembleID to gene symbol
  input <- read.table(paste0(Sys.glob(paste0(input_folder,tissue,"/",GT,"/",GT,"_diffexpr-results*"))),header = T)
  input <- input[,c("name","stat")]
  
  # run GSEA for every gene set
  res <- lapply(GO_gmt, runCamera)
  res <- do.call(rbind, res)
  res <- res[complete.cases(res), ]
  
  # Multiple test correction
  res <- res[order(res$pval), ]
  res$padj <- p.adjust(res$pval, method="BH")
  
  # Give log10 pval according to direction
  res$log_padj <- apply(res,1,function(x){
    if(x["direction"]=="Up") return(-log10(as.numeric(x["padj"]))) else return (log10(as.numeric(x["padj"])))
  })
  
  # export result
  write.table(res, paste0(output,tissue,"/",GT,"/",sample,"_GSEA.txt"), 
              sep = "\t",row.names = TRUE,col.names = TRUE,quote=FALSE)
  
  ######------ Export results  ------###### 
  
  # get significant results
  res_sig <- res[abs(res$padj)<=FDR,]
  res_sig <- res_sig[order(res_sig$padj),]
  res_sig <- res_sig%>%filter(numGenes >= 10 & numGenes <= 500 ) 
  res_sig <- setDT(res_sig, keep.rownames = TRUE)[]
  res_sig$anno <- tolower(res_sig$rn)
  res_sig$anno <- gsub("gomf_|gocc_|gobp_","",res_sig$anno)
  write.table(res_sig, paste0(output,tissue,"/",GT,"/",sample,"_GSEA_FDR",FDR,"_min10max500.txt"),
              sep = "\t",row.names = TRUE,col.names = TRUE,quote=FALSE)
  
  # visualize bar plot
  pdf(file=paste0(output,tissue,"/",GT,"/",sample,"_GSEA_barplot_FDR",FDR,"_min10max500_top50.pdf"), 
      paper="a4r",width=8, height=15)
  print(ggplot(res_sig[1:50], aes(log_padj, fct_reorder(anno, log_padj), fill = log_padj)) + 
          geom_bar(stat='identity') + 
          scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                               limits = c(min(res_sig$log_padj),max(res_sig$log_padj))) +
          theme_minimal() +
          ggtitle(paste0(tissue," ",GT)) +
          xlab("Log padj") + 
          ylab(NULL) + theme(axis.text.y = element_text(size=5)))
  dev.off() 
  
  # save RDS file
  try(assign(paste0("GSEA_",FDR,"_",sample), res))
  saveRDS(res, paste0(output,tissue,"/",GT,"/GSEA_",sample))
}

