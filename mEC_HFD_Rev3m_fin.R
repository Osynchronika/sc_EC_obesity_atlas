library(Seurat)
library(dplyr)
library(Matrix)
require(gplots)
require(ggplot2)
library(RColorBrewer)
library(cowplot)
require(scales)
library(DoubletFinder)

# Load data
setwd("~/Desktop/HI-MAG OLGA/SeqData/Cellranger/2021_06_09_M_HFD_Rev_3m/analysis")

data.dir = "~/Desktop/HI-MAG OLGA/SeqData/Cellranger/2021_06_09_M_HFD_Rev_3m/data"
data_dirs = list.files(data.dir)
data_dirs
data = list()
tissue = list()

# A. Pre-processing and QC, filtering
for (i in seq(1,length(data_dirs),1)) {
  print(data_dirs[i])
  data[[i]] <- Read10X(data.dir = paste0(data.dir,"/",data_dirs[i],"/filtered_feature_bc_matrix"))
  print(ncol(data[[i]])) # filtered cells 
  
  print("create Seurat object")
  tissue[[i]] = CreateSeuratObject(counts = data[[i]], project = as.character(data_dirs[i]))
  
  print("clear memory")
  data[[i]] = NULL
  
  print("fraction mitochondrial genes")
  tissue[[i]][["percent.mt"]] <- PercentageFeatureSet(object = tissue[[i]], pattern = "^mt-")
  
  print("Visualize QC metrics as a violin plot")
  pdf(paste0("stats_",data_dirs[i],".pdf"),width=7,height=4,paper='special')
  # use print() for the loop to export figure
  print(VlnPlot(object = tissue[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)) 
  dev.off()
  
  print("filter droplets with nuclei")
  tissue[[i]] <- subset(tissue[[i]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA < 25000 &percent.mt < 20) #nFeature_RNA > 800 & nFeature_RNA < 5000 & nCount_RNA > 800 & nCount_RNA < 15000 & 
  
  print("Visualize QC metrics as a violin plot")
  pdf(paste0("stats_filt_",data_dirs[i],".pdf"),width=7,height=4,paper='special')
  # use print() for the loop to export figure
  print(VlnPlot(object = tissue[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)) 
  dev.off()
  
  print(paste0("Number of cells ",ncol(tissue[[i]])))
  print(paste0("Median genes ",median(tissue[[i]]@meta.data$nFeature_RNA)))
  
  print("Mark organs")
  tissue[[i]]$organ = sub("\\_.*", "", data_dirs[i])
  
  print("Run the standard workflow for visualization and clustering")
  tissue[[i]] <- NormalizeData(object = tissue[[i]])
  tissue[[i]] <- FindVariableFeatures(object = tissue[[i]], selection.method = "vst", nfeatures = 3000)
  tissue[[i]] <- ScaleData(object = tissue[[i]],  verbose = FALSE) #vars.to.regress = "stim",
  tissue[[i]] <- RunPCA(object = tissue[[i]], npcs = 30, verbose = FALSE)
}

data_dirs

print("Mark diet")
for (i in c(1,4,7,10,13,16,19)) {
  tissue[[i]]$diet = "chow"
}
for (i in c(2,5,8,11,14,17,20)) {
  tissue[[i]]$diet = "Western"
}
for (i in c(3,6,9,12,15,18,21)) {
  tissue[[i]]$diet = "Rev"
}

print("Mark time")
for (i in c(1,2,4,5,7,8,10,11,13,14,16,17,19,20)) {
  tissue[[i]]$time = "6m"
}
for (i in c(3,6,9,12,15,18,21)) {
  tissue[[i]]$time = "3m"
}


cell.num = NULL
for (i in seq(1,length(data_dirs),1)) {
  print("sample")
  print(data_dirs[i])
  print(paste0("Number of cells ",ncol(tissue[[i]])))
  cell.num[i] = ncol(tissue[[i]])
  print(paste0("Median genes ",median(tissue[[i]]@meta.data$nFeature_RNA)))
}
names(cell.num) = data_dirs
cell.num

ElbowPlot(object = tissue[[10]], ndims = 30)  # n=15

# do clustering
for (i in seq(1,length(data_dirs),1)) {
  print(data_dirs[i])
  print("UMAP and Clustering")
  tissue[[i]] <- RunUMAP(object = tissue[[i]], reduction = "pca", dims = 1:30)
  tissue[[i]] <- FindNeighbors(object = tissue[[i]], reduction = "pca", dims = 1:30)
  tissue[[i]] <- FindClusters(tissue[[i]], resolution = 0.5)
  DefaultAssay(object = tissue[[i]]) <- "RNA"
  print("Plot UMAP")
  pdf(paste0("UMAP_",data_dirs[i],".pdf"),width=12,height=6,paper='special')
  p1=DimPlot(object = tissue[[i]], reduction = "umap", label=T)
  p2=FeaturePlot(tissue[[i]], features = 'nCount_RNA', order = T)
  # use print() for the loop to export figure
  print(plot_grid(p1, p2))
  dev.off()
}

plot1 <- FeatureScatter(tissue[[11]], feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(tissue[[11]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

FeaturePlot(tissue[[2]], features = 'Pecam1', order = T)
rm(data)

# B. Removing doublets with Doublet Finder
sweep.res.list = list()
sweep.stats = list()
bcmvn = list()
doublet.rate = read.table("~/Desktop/HI-MAG OLGA/10x coding/10x_multiplet.txt", header = T)
doublet.rate
# equation = 8E-06*n.cells + 0.0005
rate = 8E-06*cell.num + 0.0005
rate
for (i in seq(1,length(data_dirs),1)) {
  print(data_dirs[i])
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.list[[i]] <- paramSweep_v3(tissue[[i]], PCs = 1:30, sct = FALSE)
  sweep.stats[[i]] <- summarizeSweep(sweep.res.list[[i]], GT = FALSE)
  bcmvn[[i]] <- find.pK(sweep.stats[[i]])
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  nExp_poi <- round(rate[i]*length(rownames(tissue[[i]]@meta.data)))  ## Assuming 4% doublet formation rate for 5k cells in 10X
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  tissue[[i]] <- doubletFinder_v3(tissue[[i]], PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = F, sct = F)
  
  # result of DoubletFinder
  print("Save DotPlot with Doublets")
  # nExp_poi.adj:  DF.classifications_0.25_0.09_356
  pdf(paste0("Doublets_",data_dirs[i],".pdf"),width=7,height=7,paper='special')
  # use print() for the loop to export figure
  print(DimPlot(tissue[[i]], reduction = "umap", 
                group.by = colnames(tissue[[i]]@meta.data)[grep("DF.classifications", colnames(tissue[[i]]@meta.data))]))
  dev.off()
  print("Write Doublet Meta column")
  tissue[[i]]@meta.data$doublet = tissue[[i]]@meta.data[,grep("DF.classifications", colnames(tissue[[i]]@meta.data))]
}

DimPlot(tissue[[1]], reduction = "umap", group.by = "doublet")

for (i in seq(1,length(data_dirs),1)) {
  print("Statistics")
  write.table(table(Idents(tissue[[i]]), tissue[[i]]$doublet), file = paste0("stat_",data_dirs[i],"_doublet.txt"))
  
  print("Save RDS object")
  saveRDS(tissue[[i]],file = paste0(data_dirs[i],"_doublet.rds"))
}

rm(bcmvn, data, sweep.res.list, sweep.stats)


# combine samples
# 2. Just combine
tissue = list()
for (i in seq(1,length(data_dirs),1)) {
  tissue[[i]] = readRDS(file = paste0("QC-preproc/rmDup/",data_dirs[i],"_doublet.rds"))
}

combined = merge(tissue[[1]], y = unlist(tissue)[2:21], add.cell.ids = data_dirs, project = "mouse_HFD_Rev1m")

rm(tissue)

#combined = readRDS("combined/combined_rev3m_merged_organs.rds")
# remove genes not expressed in any cell
keep_feature <- rownames(combined)[rowSums(GetAssayData(combined,slot = "counts")) > 0]
length(keep_feature) # 25355 genes
# subset only expressed genes
combined <- combined[keep_feature, ]


# Run the standard workflow for visualization and clustering
combined <- ScaleData(object = combined, verbose = T) #, vars.to.regress = "nCount_RNA"
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 3000)
combined <- RunPCA(object = combined, npcs = 40, verbose = T)
ElbowPlot(object = combined,  ndims = 40) 

# UMAP and Clustering
combined <- RunUMAP(object = combined, reduction = "pca", dims = 1:40)
combined <- FindNeighbors(object = combined, reduction = "pca", dims = 1:40)
combined <- FindClusters(object = combined, resolution = 0.3)

# Visualization
DimPlot(object = combined, reduction = "umap", label = TRUE)
DimPlot(object = combined, reduction = "umap", group.by = "organ") # "stim"
DimPlot(object = combined, reduction = "umap", group.by = "diet")
DimPlot(object = combined, reduction = "umap", group.by = "doublet")
DimPlot(object = combined, reduction = "umap", split.by = "diet", group.by = "organ")

# Endothelial subtypes Good Markers
FeaturePlot(combined, features = c("Vwf", "Vcam1", #  large vessels
                                   "Fbln5", "Cytl1", # large artery
                                   "Hey1", "Gkn3",# artery 
                                   "Tgfb2", "Glul", # CapA
                                   "Rgcc", "Mfsd2a", # Cap
                                   "Car4", "Tfrc", # CapV
                                   "Lcn2", "Slc38a5", # large vein
                                   "Hba-a1", "Hbb-bs", 
                                   "Plvap", "Plpp3",
                                   "Flt4", "Ccl21a"),
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 4) #

# Pericyte amd SMC markers
FeaturePlot(combined, features = c("Flt1", "Pecam1","Cdh5",# EC
                                   "Pdgfrb", "Cspg4","Anpep",  #  Pericyte
                                   "Myh11", "Acta2", "Des",   # SMC
                                   "Dcn", "Col1a1", "Pdgfra"), #FB
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10") #, min.cutoff = "q10"

FeaturePlot(combined, features = c("Prox1", "Pdpn","Lyve1"), # lymph
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10") #, min.cutoff = "q10"

# 3. Find marker genes
Idents(combined) = "organ"
combined.markers = FindAllMarkers(combined, only.pos = T)
write.table(combined.markers, "combined_rev3m_organ_markers.txt", sep='\t')

top5.mark <- as.data.frame(combined.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "combined_rev3m_top5_organ_markers.txt", sep='\t')
DotPlot(combined, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

# number of cells per cluster per sample
table(combined$orig.ident)
table(Idents(combined), combined$orig.ident)
table(Idents(combined), combined$organ)
table(Idents(combined), combined$diet)

# rename clusters
Idents(combined) = "seurat_clusters"
combined <- RenameIdents(object = combined,`10` = "Lymph-EC", `18` = "FB",
                         `6` = "SMC-Mesench", `19` = "Pericyte-lung",`23` = "Pericyte-brain") 

## Relevel object@ident
combined@active.ident <- factor(x = combined@active.ident, levels = sort(levels(Idents(object = combined))))
combined$celltype <- Idents(object = combined)
DimPlot(object = combined, reduction = "umap", label = TRUE)

saveRDS(combined, "combined/combined_rev3m_merged_organs.rds")

# remove doublets
DimPlot(combined, reduction = "umap", group.by = "doublet",label = F, raster = F)  # 10.5 x 10
Idents(combined) = combined$doublet
combined = subset(combined, idents = "Singlet")
# remove unused metadata:
combined@meta.data[,grep("DF.", colnames(combined@meta.data))] = NULL
combined@meta.data[,grep("pANN", colnames(combined@meta.data))] = NULL
saveRDS(combined, "combined/combined_rev3m_merged_organs_singlet.rds")

combined = readRDS("combined/combined_rev3m_merged_organs_singlet.rds")

Idents(combined) = "celltype"
plot <- DimPlot(combined, reduction = "umap")
select.cells <- CellSelector(plot = plot) # 
Idents(combined, cells = select.cells) <- "FB"

plot <- DimPlot(combined, reduction = "umap")
select.cells <- CellSelector(plot = plot) # 
Idents(combined, cells = select.cells) <- "outlier"
write.table(WhichCells(combined, idents = "outlier"), "outliers.txt")

combined$celltype <- Idents(object = combined)

# select ECs 1
combined$ecs = "Neg"
ec.cells = WhichCells(combined, expression = Pecam1>0 | Cdh5>0)
combined$ecs[ec.cells] = "Pecam1-or-Cdh5-pos"
non.ecs = WhichCells(combined, idents = c("FB","SMC-Mesench","Pericyte-lung","Pericyte-brain", "outlier"))
combined$ecs[non.ecs] = "Neg"
DimPlot(object = combined, reduction = "umap", group.by = "ecs", cols = c("darkgrey","blue"))

#saveRDS(combined, "combined_merged_organs_singlet.rds")

# subset only ECs
Idents(combined) = combined$ecs
combined = subset(combined, idents = "Pecam1-or-Cdh5-pos")
saveRDS(combined, "combined/combined_rev3m_merged_organs_ECs.rds")

Idents(combined) = combined$celltype
DimPlot(combined, reduction = "umap", group.by = "organ",label = TRUE)  # 10.5 x 10
DimPlot(combined, reduction = "umap", group.by = "diet",label = F, cols = c("#165bc9","#14e39e","#f79525"), 
        raster = F)
DimPlot(combined, reduction = "umap", split.by = "diet",label = F, group.by = "organ",
        raster = F)

combined$diet_time = paste0(combined$diet, combined$time)

DimPlot(combined, reduction = "umap", split.by = "diet_time",label = F, group.by = "organ",
        raster = F)

# 4.  split by organ, save
Idents(combined) = combined$organ
for (i in unique(combined$organ)) {
  saveRDS(subset(combined,idents = i), paste0(i,'_subset.rds'))
}

# save lymphatic EC separate
Idents(combined) = combined$celltype
saveRDS(subset(combined,idents = "Lymph-EC"), 'lymphatic_subset.rds')

