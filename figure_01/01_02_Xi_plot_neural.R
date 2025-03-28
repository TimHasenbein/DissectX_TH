##*************************##
##  Allelome.PRO analysis  ##
##*************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 09.2022
## Creation: 02.2022
library(tidyverse)
library(ggplot2)
library(data.table)


######------  set environment  ------###### 
input <-"./02_allelomePRO/"
output <-"./Xi_specific_peaks/"


######------  get AP results  ------###### 
# read in locus_table.txt for each sample
read <- function(x) read.table(x, header=T)
dataFiles <- lapply(Sys.glob(paste0(input,"*/locus_table.txt")), read)
# name each data frame
sample_names <- Sys.glob(paste0(input,"*/locus_table.txt"))
sample_names <- gsub(paste0(input,"2022_03_22_"),"",sample_names)
sample_names <- gsub(paste0(input,"2022_04_19_"),"",sample_names)
sample_names <- gsub("_pass.rmChrM.rmAF.rmBl.rmDup_50kb_swindow.bed_1/locus_table.txt","",sample_names)
names(dataFiles) <- sample_names
# filter for chrX genes and genes including >= 50 reads
filter <- function(x) {
  x <- x[x$chr == "chrX",]
  x <- x[x$total_reads >= 50,]
}
dataFiles_filter <- lapply(dataFiles, filter)


######------ Generate matrix for plots  ------###### 
# select allelic ratio 
select <- function(x) {
  x <- x[,c("name","allelic_ratio","start","end")]
}
dataFiles_filter_select <- lapply(dataFiles_filter, select)
# bind results
res <- bind_rows(dataFiles_filter_select, .id = "sample")


######------ Select region of interest  ------###### 
res$locus <- NA
for(r in 1:nrow(res)){
  ifelse(res[r,"end"]>=50530000 & res[r,"start"]<=50640000, 
         res[r,"locus"] <- "Crossfirre/Firre",
         ifelse(res[r,"end"]>=75723000 & res[r,"start"]<=75766000, 
                res[r,"locus"] <- "Dxz4", 
                res[r,"locus"] <- res[r,"locus"]))
}
res$locus <- factor(res$locus, levels = c("Crossfirre/Firre","Dxz4")) 


######------  Violin plot with min value only  ------###### 
plot_1 <- res
# match direction
plot_1$allelic_ratio <- 1 - plot_1$allelic_ratio
# filter for min allelic ratio per sample and gene locus
min_cf_s1 <- plot_1%>%dplyr::filter(sample == "SRR3933589" & locus == "Crossfirre/Firre")
min_cf_s1 <- min_cf_s1%>%dplyr::filter(allelic_ratio == min(min_cf_s1$allelic_ratio))
min_cf_s2 <- plot_1%>%dplyr::filter(sample == "SRR3933595" & locus == "Crossfirre/Firre")
min_cf_s2 <- min_cf_s2%>%dplyr::filter(allelic_ratio == min(min_cf_s2$allelic_ratio))
min_d_s1 <- plot_1%>%dplyr::filter(sample == "SRR3933589" & locus == "Dxz4")
min_d_s1 <- min_d_s1%>%dplyr::filter(allelic_ratio == min(min_d_s1$allelic_ratio))
min_d_s2 <- plot_1%>%dplyr::filter(sample == "SRR3933595" & locus == "Dxz4")
min_d_s2 <- min_d_s2%>%dplyr::filter(allelic_ratio == min(min_d_s2$allelic_ratio))
label <- min_cf_s1%>%rbind(min_cf_s2)%>%rbind(min_d_s1)%>%rbind(min_d_s2)
# violin plot
plot_1$sample <- factor(plot_1$sample, levels = c("SRR3933589","SRR3933595")) 
pdf(file=paste0(output,"/Xi_plot_top_annotation.pdf"),width=10, height=8)
ggplot(plot_1, aes(x=sample, y=allelic_ratio, fill=sample)) + 
  geom_violin(width= 0.25) +
  geom_boxplot(width=0.2,fill="white") + 
  labs(x="",y="Allelic ratio") +
  scale_fill_manual(values = c("grey","grey")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0,0.5,1), minor_breaks = seq(0,1,0.25)) +
  geom_point(data=label,aes(x=sample, y=allelic_ratio),color="red", size=3) +
  geom_text(data=label, aes(x=sample, y=allelic_ratio, label=locus), size=3) +
  theme(axis.text.x=element_text(size=9),
        axis.ticks.x=element_blank(),
        axis.line = element_line(colour = 'black', size = 0.5),
        axis.ticks.y=element_line(colour = "black",size=0.5),
        axis.title.y=element_text(size=9),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        panel.grid.minor = element_line(colour = "darkgrey", size = 0.4),
        panel.grid.major.y = element_line(colour = "darkgrey", size = 0.4),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.spacing = unit(1.0, 'cm'))
dev.off()


