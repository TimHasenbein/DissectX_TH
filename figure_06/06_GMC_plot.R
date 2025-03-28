##************************##
##  Figure 6 GMC results  ## 
##************************##
## Project: dissectX
## Tim Hasenbein
## Last modification 05.2024
## Creation: 07.2023
## Plot the results of the GMC screening 
library(tidyverse)
library(ggplot2)
library(readxl)
library(ggpubr)
library(ComplexHeatmap)
library(xlsx)


######------ set environment -----###### 
input <- "./07_GMC/"
output <- "./07_GMC/"
abbr <- read.xlsx(paste0(input,"Abbreviations_parameter_names_GMC_2.xlsx"),sheetIndex = 1)[,2:4]
colnames(abbr) <- c("test","pname","long_name")


######------ read data in and add immuno screen -----###### 
# sex-independent
both <- read.xlsx(paste0(input,"effsize_hedge_all._08052024_1.xlsx"),sheetIndex = 1)[,2:8]
both_immu <- read.csv(paste0(input,"effsize_hedge_immuno_all.csv"))
both <- both%>%rbind(both_immu)
both$sample <- "sex-independent"
# female-specific
female <- read.xlsx(paste0(input,"effsize_hedge_f._08052024_1.xlsx"),sheetIndex = 1)[,2:8]
female_immu <- read.csv(paste0(input,"effsize_hedge_immuno_f.csv"))
female <- female%>%rbind(female_immu)
female$sample <- "female-specific"
# male-specific
male <- read.xlsx(paste0(input,"effsize_hedge_m._08052024_1.xlsx"),sheetIndex = 1)[,2:8]
male_immu <- read.csv(paste0(input,"effsize_hedge_immuno_m.csv"))
male <- male%>%rbind(male_immu)
male$sample <- "male-specific"


######------ filter screens and merge -----###### 
both <- both%>%filter(!pname %in% c("conc_IgE","conc_IL-6","TIBC","AUC_total","PLT","H","ThresholdAt18kHz","ClickABR","Heart_rate_Echo"))
male <- male%>%filter(!pname %in% c("conc_IgE","conc_IL-6","TIBC","AUC_total","PLT","H","ThresholdAt18kHz","ClickABR","Heart_rate_Echo"))
female <- female%>%filter(!pname %in% c("conc_IgE","conc_IL-6","TIBC","AUC_total","PLT","H","ThresholdAt18kHz","ClickABR","Heart_rate_Echo"))
merge <- both%>%rbind(female)%>%rbind(male)
rm(both_immu,female_immu,male_immu)


######------ plot data -----###### 
merge$pname <- factor(merge$pname, levels=c(both$pname))
merge$sample <- factor(merge$sample, levels=c("sex-independent","female-specific","male-specific"))
merge$p <- as.numeric(merge$p)
merge <- merge %>% mutate(new_p = ifelse(p<0.05,p,1))
merge$abs_d <- abs(as.numeric(merge$d))
merge$new_p <- as.numeric(merge$new_p)
merge <- merge %>% mutate(shape = ifelse(d>0,"up","down"))
pdf(file=paste0(output,"/overview_new.pdf"),width=15, height=12)
ggballoonplot(
  merge, x = "sample", y = "pname", shape = "shape",
  size = "abs_d", fill = "new_p",
    ggtheme = theme_bw()
) + facet_wrap("screen_name", scales = "free") +
  scale_fill_gradient(
    low = "#E60001",
    high = alpha("#FF8B01",0.5),
    breaks = c(0,0.01,0.02,0.03,0.04,0.05),
    labels = c(0,0.01,0.02,0.03,0.04,0.05),
    name = "Custom Gradient",
    limits = c(0, 0.05)
  ) + facet_wrap("screen_name", scales = "free") +
  scale_size_continuous(
    range = c(0,7),  
    breaks = c(0,0.5,1,1.5,2,2.5,3),
    name = "Cohen's d effect size"     
  ) +
  scale_shape_manual(values = c(25,24))
dev.off()


######------ plot number of sig. phenotypes -----###### 
plot <- merge%>%filter(p<0.05)
plot$heat <- 1
plot <- plot%>%dplyr::select(pname,sample,heat)
plot <- as.data.frame(plot)
plot <- plot%>%pivot_wider(names_from=sample, values_from = heat, values_fill = 0)
options(scipen = 9999999)
plot <- plot%>%column_to_rownames("pname")
pdf(file=paste0(output,"GMC_heatmap_new.pdf"),width=5, height=10)
ht <- Heatmap(as.matrix(plot), 
              cluster_rows = T, 
              cluster_columns = F, 
              column_title = NULL, 
              show_row_names = T, 
              row_names_gp =  gpar(fontsize = 4),
              col = c("grey","darkred"),
              width = unit(10, "cm"), border = T,
              row_gap = unit(0.5, "mm"), 
              row_title_gp = gpar(fontsize = 6, col="black"),
              column_names_gp = gpar(fontsize = 10),
              heatmap_legend_param = list(at = c(0,1))) 
draw(ht, heatmap_legend_side = "left")
dev.off()

