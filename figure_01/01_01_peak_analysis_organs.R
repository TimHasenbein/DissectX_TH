##*************************************##
##  Top female-specific loci analysis  ##
##*************************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 09.2022
## Creation: 03.2022
## Get differences btw males and females in terms of ATAC accessibility
library(plyr)
library(tidyverse)
library(gtools)
library(karyoploteR)
library(biomaRt)


######------  set environment  ------###### 
bedtag <- c("chr","start","end","name")
input <- "./02_intersect_peaks/broad/"
output <- "./output_broad_Liu_2019/"


######------  function to calculate allelic score ------###### 
calc.allelic.score <- function(x) {
   allelic.score <- log10(pbinom(min(matrix(x,ncol=2)), round(x[1]+x[2]), 0.5, lower.tail=TRUE, log = FALSE))
  ifelse(x[1] > x[2], -allelic.score,allelic.score)
}


######------  get data ------###### 
samples <- c("ce","he","ki","li","lu","sp")
for (sample in samples){
  # female
  data_female <- read.table(paste0(input,sample,"_female_count_broad.bed"), header=FALSE,
                            check.names = FALSE,stringsAsFactors = FALSE,comment.char = "",fill = TRUE)[,c(1:4,7)]
  colnames(data_female)=c(bedtag,paste0(sample,"_female_peak_counts_broad"))
  assign(paste0(sample,"_female"), data_female)
  # male
  data_male <-read.table(paste0(input,sample,"_male_count_broad.bed"), header=FALSE,
                         check.names = FALSE,stringsAsFactors = FALSE,comment.char = "",fill = TRUE)[,c(1:4,7)]
  colnames(data_male)=c(bedtag,paste0(sample,"_male_peak_counts_broad"))
  assign(paste0(sample,"_male"), data_male)
}
# merge female/male peak count tables
peak_counts_female <- ce_female%>%full_join(he_female,by=c("chr","start","end","name"))%>%full_join(li_female,by=c("chr","start","end","name"))%>%full_join(lu_female,by=c("chr","start","end","name"))%>%full_join(ki_female,by=c("chr","start","end","name"))%>%full_join(sp_female,by=c("chr","start","end","name"))
peak_counts_male <- ce_male%>%full_join(he_male,by=c("chr","start","end","name"))%>%full_join(li_male,by=c("chr","start","end","name"))%>%full_join(lu_male,by=c("chr","start","end","name"))%>%full_join(ki_male,by=c("chr","start","end","name"))%>%full_join(sp_male,by=c("chr","start","end","name"))


######------  get median values ------###### 
peak_counts_female$median_female <- round(apply(peak_counts_female[,5:ncol(peak_counts_female)], 1, median)) 
peak_counts_male$median_male <- round(apply(peak_counts_male[,5:ncol(peak_counts_male)], 1, median))
female_median <- peak_counts_female[,c(1:4,ncol(peak_counts_female))]
male_median <- peak_counts_male[,c(1:4,ncol(peak_counts_male))]
write.table(peak_counts_female, paste0(output,"peak_counts_female.txt"),sep="\t",quote=F,col.names=T,row.names=F)
write.table(peak_counts_male, paste0(output,"peak_counts_male.txt"),sep="\t",quote=F,col.names=T,row.names=F)


######------  join data ------###### 
peak_counts <- join(female_median,male_median)
peak_counts$score <- round(apply(peak_counts[,c("median_female","median_male")],1,calc.allelic.score),3)
# mark regions of interest 
peak_counts$label <- ""
peak_counts$color <- "darkgrey"
# firre
peak_counts[peak_counts$name %in% c("chrX:50550000-50650000","chrX:50500000-50600000","chrX:50600000-50700000"),"label"] <- "Firre"
peak_counts[peak_counts$name %in% c("chrX:50550000-50650000","chrX:50500000-50600000","chrX:50600000-50700000"),"color"] <- "darkred"
# dxz4
peak_counts[peak_counts$name %in% c("chrX:75700000-75800000","chrX:75650000-75750000"),"label"] <- "Dxz4"
peak_counts[peak_counts$name %in% c("chrX:75700000-75800000","chrX:75650000-75750000"),"color"] <- "black"
# xist
peak_counts[peak_counts$name %in% c("chrX:103400000-103500000","chrX:103450000-103550000"),"label"] <- "Xist"
peak_counts[peak_counts$name %in% c("chrX:103400000-103500000","chrX:103450000-103550000"),"color"] <- "darkblue"
# order by score
peak_counts <- peak_counts[order(peak_counts$score,decreasing = F),]


######------  plot data analysis results ------###### 
pdf(file=paste(output,"rank_plot.pdf",sep=""), paper="a4r",width=8, height=7)
par(mfrow=c(1,1),mar = c(5,10,5,5),cex.main=1,cex.axis=0.5,cex=1,pty ="s")
plot(peak_counts$score,ylim=c(-2,5),xlim=c(1,nrow(peak_counts)),pch = 17 , las=2,type = "n",
     col=peak_counts$color,cex=1+(abs(peak_counts$score)*2),
     ylab="Log10(p)",xlab="Rank")
segments(x0 = 1:nrow(peak_counts) - abs(peak_counts$score)*200, y0 = peak_counts$score,
         x1 = 1:nrow(peak_counts) + abs(peak_counts$score)*200, y1 = peak_counts$score, col = peak_counts$color)
x_pos=which(peak_counts$label %in% c("Firre","Xist","Dxz4"))
text(x = 1:nrow(peak_counts)-3000,y = peak_counts$score,labels = peak_counts$label)
abline(h = 0)
dev.off()
# karyotype plot
pdf(file=paste(output,"karyotype.pdf",sep=""), paper="a4r",width=10, height=12)
par(cex.main=1)
peak_counts_X <- peak_counts[peak_counts$chr=="chrX",]
peak_counts_X <- peak_counts_X[order(peak_counts_X$start),]
kp <- NULL
kp <- plotKaryotype(genome="mm10",chromosomes="chrX",plot.type = 4)
kpAddBaseNumbers(kp)
kpAxis(kp,r0= 0,r1=0.5,ymin=0, ymax=1,tick.pos = c(seq(0,1,0.2)), labels = c(seq(0,5,1)))
kpArea(kp, chr="chrX", x=peak_counts_X$start, y=(peak_counts_X$score/5),col = "black", border = "black", r0 = 0, r1 = 0.5)
kpAbline(kp, chr = "chrX",v = c(50550000,103458373,75700000), r0=1, r1=1.1,col = "red")
kpText(kp, chr="chrX", x=c(50550000,75700000,103458373), y=1.12, col="black", labels=c("Firre", "Dxz4", "Xist"), cex=0.7)
dev.off()


######------  save results ------###### 
peak_counts <- peak_counts[order(peak_counts$score,decreasing = TRUE),]
write.table(peak_counts, paste0(output,"peak_counts.txt"),sep="\t",quote=F,col.names=T,row.names=F)

