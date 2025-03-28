##****************************##
## scRNA-seq: Quality control ##
##**********'*****************##
## Project: dissectX
## Tim Hasenbein
## Last modification 05.2022
## Creation: 04.2022
library(Seurat)
library(tidyverse)
library(biomaRt)


######------ Set environment ------###### 
input <- "./00_cellranger_data/"
output <- "./01_preprocession/"
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="http://aug2020.archive.ensembl.org")


######------ Read data and make seurat object ------###### 
samples <- c("LFDxCast_KO","LFDxCast_WT")
for (sample in samples){
  seurat_data <- Read10X(data.dir = paste0(input,sample))
  seurat_obj <- CreateSeuratObject(counts = seurat_data$`Gene Expression`, 
                                   min.features = 100, 
                                   project = sample)
  assign(sample, seurat_obj)}


######------ Merge data ------###### 
merged_seurat <- merge(x = LFDxCast_KO, y = LFDxCast_WT, project = "dissectX")
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)
rm(seurat_data,seurat_obj,LFDxCast_KO,LFDxCast_WT)


######------ Add QC metrics of data to metadata ------###### 
# create a sample column
merged_seurat$cell <- rownames(merged_seurat@meta.data)
# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'orig.ident', into = c('cross', 'condition'), sep = '_')
# add number of genes per UMI for each cell to metadata
merged_seurat[["GenesPerUMI"]] <- merged_seurat$nFeature_RNA / merged_seurat$nCount_RNA
# add mt percentage per cell
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^mt-") 
# add ribosomal percentage (impression)
merged_seurat[["percent.rb"]] <- PercentageFeatureSet(merged_seurat, pattern = "^Rps|Rpl")
head(merged_seurat@meta.data)


######------ Check QC metrics of raw data ------###### 
pdf(file = paste0(output,"QC_metrics_raw.pdf"), pointsize=2)
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
merged_seurat@meta.data %>% 
  ggplot(aes(color=condition, x=nCount_RNA, fill= condition)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 2000) + 
  geom_vline(xintercept = 20000)
merged_seurat@meta.data %>% 
  ggplot(aes(color=condition, x=nFeature_RNA, fill= condition)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 500) + 
  geom_vline(xintercept = 5000)
merged_seurat@meta.data %>% 
  ggplot(aes(color=condition, x=percent.mt, fill=condition)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 10)
dev.off()
summary(merged_seurat@meta.data$nCount_RNA)
sd(merged_seurat@meta.data$nCount_RNA) 
summary(merged_seurat@meta.data$nFeature_RNA)
sd(merged_seurat@meta.data$nFeature_RNA) 


######------ Filter low quality data ------###### 
# cell-level filtering: base on min genes, transcripts and mt %
srat_clean <- subset(merged_seurat, subset = nCount_RNA > 2000 & nCount_RNA < 20000 & nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
# gene-level filtering: keep only genes which are expressed in 10 or more cells
counts <- GetAssayData(object = srat_clean, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
srat_clean <- CreateSeuratObject(filtered_counts, meta.data = srat_clean@meta.data)


######------ Check QC metrics of clean data ------###### 
Idents(srat_clean) <- srat_clean$condition
pdf(file = paste0(output,"QC_metrics_clean.pdf"), pointsize=2)
VlnPlot(srat_clean, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(srat_clean, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
srat_clean@meta.data %>% 
  ggplot(aes(color=condition, x=nCount_RNA, fill= condition)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() 
srat_clean@meta.data %>% 
  ggplot(aes(color=condition, x=nFeature_RNA, fill= condition)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() 
srat_clean@meta.data %>% 
  ggplot(aes(color=condition, x=percent.mt, fill=condition)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() 
dev.off()


######------ Save filtered seurat object ------###### 
saveRDS(srat_clean, file=paste0(output,"01_QC_srat_clean_XX_KO_WT.rds"))
nrow(srat_clean@meta.data %>% filter(condition == "KO")) 
nrow(srat_clean@meta.data %>% filter(condition == "WT")) 


######------ Normalization with SCTransform ------###### 
# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(srat_clean, split.by = "condition")
# use SCTransform
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], verbose = FALSE, vars.to.regress = "percent.mt")
} # default assay is now SCT


######------ Data integration ------###### 
pre_int <- FindVariableFeatures(object = srat_clean)
pre_int <- ScaleData(object = pre_int)
pre_int <- RunPCA(object = pre_int)
pre_int <- RunUMAP(pre_int, 
                   dims = 1:40,
                   reduction = "pca")
pdf(file = paste0(output,"02_pre_integration.pdf"),width=20,height=10)
DimPlot(pre_int, split.by = "condition")|DimPlot(pre_int)     
dev.off()
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000) 
split_seurat <- PrepSCTIntegration(object.list = split_seurat,anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
# Visualizate data integration
seurat_integrated <- RunPCA(object = seurat_integrated)
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")
pdf(file = paste0(output,"02_post_integration_data.pdf"),width=20,height=10)
DimPlot(seurat_integrated, split.by = "condition")|DimPlot(seurat_integrated)     
dev.off()


######------ Cluster identification ------######
srat_int <- seurat_integrated
srat_int <- FindNeighbors(srat_int, dims = 1:40, verbose = FALSE)
srat_int <- FindClusters(srat_int, resolution = c(0.4)) 
Idents(srat_int); head(srat_int@meta.data)
pdf(file = paste0(output,"03_cluster_umap.pdf"),width=15,height=10)
DimPlot(srat_int, reduction = "umap", label = TRUE, label.size = 6, split.by = "condition")
dev.off()


#---------------------------------------#
######------ Cell annotation ------###### 
#---------------------------------------#
DefaultAssay(srat_int) <- "RNA"


######------ Manual assignment ------######
DimPlot(srat_int,label=T)
DotPlot(srat_int,features = c("Itgam","Ccr2","S100a8","Cd3g","Cd8a","Cd4","Cxcr3","Ccl5","Cd19"),cols = c("grey", "red")) + coord_flip()
FeaturePlot(srat_int, feature="Itgam", split.by = "condition",max.cutoff = 3, cols = c("grey", "red")) # macrophages
FeaturePlot(srat_int, feature="Csf1r", split.by = "condition",max.cutoff = 3, cols = c("grey", "red")) # monocytes
FeaturePlot(srat_int, feature="S100a8", split.by = "condition",max.cutoff = 3, cols = c("grey", "red")) # granuloccytes
FeaturePlot(srat_int, feature="Cd3g", split.by = "condition",max.cutoff = 3, cols = c("grey", "red")) # T cells
FeaturePlot(srat_int, feature="Cd8a", split.by = "condition",max.cutoff = 3, cols = c("grey", "red")) # T cells CD8
FeaturePlot(srat_int, feature="Cd4", split.by = "condition",max.cutoff = 3, cols = c("grey", "red")) # T cells CD4
FeaturePlot(srat_int, feature="Cxcr3", split.by = "condition",max.cutoff = 3, cols = c("grey", "red")) # Thelper1
FeaturePlot(srat_int, feature="Ccl5", split.by = "condition",max.cutoff = 3, cols = c("grey", "red")) # Nk cells
FeaturePlot(srat_int, feature="Cd79a", split.by = "condition",max.cutoff = 3, cols = c("grey", "red")) # B cells


######------ Rename identities females ------###### 
srat_lab <- RenameIdents(object = srat_int, 
                         "0" = "B cells",
                         "1" = "B cells",
                         "2" = "T cells CD4",
                         "3" = "T cells CD8",
                         "4" = "B cells",
                         "5" = "B cells",
                         "6" = "Monocytes",
                         "7" = "T cells",
                         "8" = "Nk cells",
                         "9" = "Granulocytes",
                         "10"= "B cells")
pdf(file=paste0(output,"03_annotation_manual.pdf"),width=20,height=10)
DimPlot(srat_lab,label=T,split.by = "condition") 
dev.off()


######------ Save final R object ------###### 
saveRDS(srat_lab, file=paste0(output,"03_srat_clean_integrated_reg_MT_labelled_XX_KO_WT.rds"))
srat <- readRDS(paste0(output,"03_srat_clean_integrated_reg_MT_labelled_XX_KO_WT.rds"))


######------ UMAP plot ------######
pdf(file = paste0(output,"03_cluster_umap_labelled.pdf"),width=15,height=10)
DimPlot(srat, reduction = "umap", label = TRUE, label.size = 6, split.by = "condition")
dev.off()
