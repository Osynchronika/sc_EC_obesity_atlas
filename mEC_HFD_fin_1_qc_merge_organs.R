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
setwd("~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis")

data.dir = "~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/data"
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
for (i in c(1,3,4,7,9,11,13,15)) {
  tissue[[i]]$diet = "chow"
}
for (i in c(2,5,6,8,10,12,14,16)) {
  tissue[[i]]$diet = "Western"
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

FeaturePlot(tissue[[11]], features = 'Pecam1', order = T)
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
  saveRDS(tissue[[i]], file = paste0(data_dirs[i],"_doublet.rds"))
}

rm(bcmvn, data, sweep.res.list, sweep.stats)

# combine samples
# 2. Just combine
combined = merge(tissue[[1]], y = unlist(tissue)[2:16], add.cell.ids = data_dirs, project = "mouse_EC_HFD")

rm(tissue)
rm(p1,p2,plot1,plot2)

combined = readRDS("combined/combined_merged_organs.rds")
# remove genes not expressed in any cell
keep_feature <- rownames(combined)[rowSums(GetAssayData(combined,slot = "counts")) > 0]
length(keep_feature) # 24333 genes
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
                                    "Myh11", "Acta2", "Tagln",   # SMC
                                    "Dcn", "Col1a1", "Pdgfra"), #FB
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10", pt.size = 1) #, min.cutoff = "q10"

FeaturePlot(combined, features = c("Prox1", "Pdpn","Lyve1",
                                   "Flt4", "Ccl12","Fgl2"), # lymph
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10", pt.size = 1) #, min.cutoff = "q10"

FeaturePlot(combined, features = c("Ptprc", "Mrc1","Cd74","Igkc", "Cd52"), # immune
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10", pt.size = 1) 

FeaturePlot(combined, features = c("Cd19", "Bcl11a","Nkg7","Pax5"), # immune
            order = T, cols=c("grey", "red"), ncol=2, min.cutoff = "q10") 



# 3. Find marker genes

combined.markers = FindAllMarkers(combined, only.pos = T)
write.table(combined.markers, "combined_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(combined.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "combined_top5_cluster_markers.txt", sep='\t')
DotPlot(combined, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

# number of cells per cluster per sample
table(Idents(combined), combined$orig.ident)
table(Idents(combined), combined$organ)
table(Idents(combined), combined$diet)

# rename clusters
combined <- RenameIdents(object = combined,`7` = "LECs", `11` = "SMC", 
                               `17` = "Pericyte-brain", `20` = "FB",
                                `21` = "Pericyte-lung",
                               `16` = "FB-heart",
                               `24` = "Immune") 

## Relevel object@ident
combined@active.ident <- factor(x = combined@active.ident, levels = sort(levels(Idents(object = combined))))
combined$celltype <- Idents(object = combined)
DimPlot(object = combined, reduction = "umap", label = TRUE)

# remove doublets
DimPlot(combined, reduction = "umap", group.by = "doublet",label = F, raster = F)  # 10.5 x 10
Idents(combined) = combined$doublet
combined = subset(combined, idents = "Singlet")
# remove unused metadata:
combined@meta.data[,grep("DF.", colnames(combined@meta.data))] = NULL
combined@meta.data[,grep("pANN", colnames(combined@meta.data))] = NULL

saveRDS(combined, "combined/combined_merged_organs.rds")

combined = readRDS("combined/combined_merged_organs.rds")


# select ECs 
combined$ecs = "Neg"
ec.cells = WhichCells(combined, expression = Pecam1>0 | Cdh5>0)
combined$ecs[ec.cells] = "Pecam1-or-Cdh5-pos"
non.ecs = WhichCells(combined, idents = c("SMC","Pericyte-brain", "Pericyte-lung","FB","FB-heart","Immune"))
combined$ecs[non.ecs] = "Neg"
DimPlot(object = combined, reduction = "umap", group.by = "ecs", cols = c("darkgrey","blue"))

# select cluster by hand (first time)
plot <- DimPlot(combined, reduction = "umap")
select.cells <- CellSelector(plot = plot) # 
#outliers = select.cells
outliers = c(outliers, select.cells)
outliers = as.vector(outliers)
write(outliers, "combined/outlier_cells_combined.txt")

excl.cells = read.table("combined/outlier_cells_combined.txt")
DimPlot(combined, cells.highlight = c(excl.cells$V1))
combined$ecs[excl.cells$V1] = "Neg"
DimPlot(object = combined, reduction = "umap", group.by = "ecs", cols = c("darkgrey","blue"))


# subset only ECs
Idents(combined) = combined$ecs
combined = subset(combined, idents = "Pecam1-or-Cdh5-pos")
DimPlot(object = combined, reduction = "umap", group.by = "organ")

saveRDS(combined, "combined/combined_merged_organs_ECs_singlet.rds")

Idents(combined) = combined$celltype
DimPlot(combined, reduction = "umap", group.by = "organ",label = TRUE)  # 10.5 x 10
DimPlot(combined, reduction = "umap", group.by = "diet",label = F, cols = c("#165bc9","#f79525"))


Idents(combined) = combined$organ
combined.markers = FindAllMarkers(combined, only.pos = T)
write.table(combined.markers, "combined_organ_markers.txt", sep='\t')

top5.mark <- as.data.frame(combined.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
write.table(top5.mark, "combined_top5_organ_markers.txt", sep='\t')

# lipid
FeaturePlot(combined, features = c("Cd36", "Fabp1","Fabp4",
                                   "Fabp5", "Alb","Abca1",
                                   "Apoe","Apoa1","Apoa2",
                                   "Apoc1","Mgll","Lpl"), 
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10", pt.size = 1) #, min.cutoff = "q10"

FeaturePlot(combined, features = c("Dbi","Lipe","Pnpla2","Nr1h4"), 
            order = T, cols=c("grey", "red"), ncol=2, min.cutoff = "q10", pt.size = 1) #, min.cutoff = "q10"

FeaturePlot(combined, features = c("Slc16a1"), 
            order = T, cols=c("grey", "red"), min.cutoff = "q10", pt.size = 1, raster = T) #, min.cutoff = "q10"

FeaturePlot(combined, features = c("Npr3"), 
            order = T, cols=c("grey", "red"), min.cutoff = "q10", pt.size = 1, raster = T) #, min.cutoff = "q10"

FeaturePlot(combined, features = c("Mt1","Mt2","Sgms1","Degs2"), 
            order = T, cols=c("grey", "red"),ncol=2, min.cutoff = "q10", pt.size = 1, raster = T) #, min.cutoff = "q10"

FeaturePlot(combined, features = c("Slc34a1",
                                   "Slc4a4",
                                   "Slc6a19",
                                   "Slc5a2",
                                   "Slc22a8",
                                   "Slc13a1",
                                   "Slc22a18",
                                   "Slc5a12",
                                   "Slc2a2"), 
            order = T, cols=c("grey", "red"),ncol=3, min.cutoff = "q10", pt.size = 1, raster = T) #, min.cutoff = "q10"

FeaturePlot(brain, features = c("Slc16a1"), 
            order = T, cols=c("grey", "red"), min.cutoff = "q10", pt.size = 1, raster = T) #, min.cutoff = "q10"


VlnPlot(combined, features = c("Cd36", "Fabp1","Fabp4",
                                   "Fabp5", "Alb","Abca1",
                                   "Apoe","Apoa1","Apoa2",
                                   "Apoc1","Mgll","Lpl"), 
           group.by = "organ", split.by = "diet", ncol=3, pt.size = 0) #, min.cutoff = "q10"

FeaturePlot(combined, features = c("Klf2","Klf4","Klf6","Nfkbia","Meox2","Tcf15"), 
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10", pt.size = 1) #, min.cutoff = "q10"
VlnPlot(combined, features = c("Meox2","Tcf15"), group.by = "organ", split.by = "diet", pt.size = 0, ncol = 1)



# 4.  split by organ, save
Idents(combined) = combined$organ
for (i in unique(combined$organ)) {
  saveRDS(subset(combined,idents = i), paste0(i,'_subset.rds'))
}

# save lymphatic EC separate
Idents(combined) = combined$celltype
saveRDS(subset(combined,idents = "Lymph-EC"), 'lymphatic_subset.rds')


########################################
# analyze organs separately
########################################
# 1. Liver
############
setwd("liver")
liver = readRDS("../liver_subset.rds")

Idents(liver) = liver$celltype
DimPlot(liver, reduction = "umap", label = TRUE)
DimPlot(liver, reduction = "umap", group.by = "ecs")

# subset only ECs
Idents(liver) = liver$ecs
liver = subset(liver, idents = "Pecam1-or-Cdh5-pos")
saveRDS(liver, "liver_ECs.rds")

Idents(liver) = liver$celltype
DimPlot(liver, reduction = "umap", group.by = "diet")

# Endothelial subtypes Good Markers
FeaturePlot(liver, features = c("Vwf", "Vcam1", #  large vessels
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
FeaturePlot(liver, features = c("Flt1", "Pecam1","Cdh5",# EC
                                   "Pdgfrb", "Cspg4","Anpep",  #  Pericyte
                                   "Myh11", "Acta2", "Des",   # SMC
                                   "Dcn", "Col1a1", "Pdgfra"), #FB
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10") #, min.cutoff = "q10"

FeaturePlot(liver, features = c("Prox1", "Pdpn","Lyve1"), # lymph
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10") #, min.cutoff = "q10"


liver <- FindVariableFeatures(liver, selection.method = "vst", nfeatures = 3000)
liver <- RunPCA(object = liver, npcs = 30, verbose = T)
ElbowPlot(object = liver,  ndims = 30) 

# UMAP and Clustering
liver <- RunUMAP(object = liver, reduction = "pca", dims = 1:12)
liver <- FindNeighbors(object = liver, reduction = "pca", dims = 1:12)
liver <- FindClusters(object = liver, resolution = 0.3)

# Visualization
DimPlot(object = liver, reduction = "umap", label = TRUE)
DimPlot(object = liver, reduction = "umap", group.by = "diet")
DimPlot(object = liver, reduction = "umap", group.by = "celltype", label = T)

# Find marker genes
Idents(liver) = liver$seurat_clusters
liver.markers = FindAllMarkers(liver, only.pos = T)
write.table(liver.markers, "liver_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(liver.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "liver_top5_cluster_markers.txt", sep='\t')
DotPlot(liver, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

# number of cells per cluster per sample
table(Idents(liver), liver$diet)

# DE genes
liver$clust.diet <- paste(liver$seurat_clusters, liver$diet, sep = "_")
Idents(object = liver) <- "clust.diet"
liver@active.ident <- factor(x = liver@active.ident, levels = sort(levels(Idents(object = liver))))
levels(Idents(object = liver))

lev1 <- levels(liver@active.ident)[seq(2,(length(levels(Idents(object = liver)))),2)] # condition 1=HFD
lev2 <- levels(liver@active.ident)[seq(1,(length(levels(Idents(object = liver)))),2)] # condition 2=chow

# for all the clusters
for (i in seq(1,length(lev1),1)) { # 16 for not-filtered
  print(lev1[i])
  print(lev2[i])
  # TG 19m vs WT 19m
  response2 <- FindMarkers(object = liver, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = FALSE) #, min.pct = 0.25
  write.table(response2, file = paste(lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
}

# find DEG for all ECs (like bulk), all cells are capillary here
Idents(object = liver) <- "diet"
response2 <- FindMarkers(object = liver, ident.1 = "Western", 
                         ident.2 = "chow", 
                         logfc.threshold = 0,
                         verbose = FALSE) #, min.pct = 0.25
write.table(response2, "all_EC_liver_WD_vs_chow.txt", sep = "\t")


# markers for liver from Murine Atlas
FeaturePlot(liver, features = c("Vwf", "Edn1","Bgn", #  vein
                                "Clu", "Plac8","Lrg1", # large vessels
                                 "Efnb1", "Glul", "Gja4","Sox17",# CapA
                                "Stab2", "Kdr", # Cap
                                 "Thbd"), # CapV
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 4) #


saveRDS(liver, "liver_EC.rds")


########################################
# 2. Fat
############
sc = readRDS("sc_subset.rds")
vis = readRDS("vis_subset.rds")

##############
# A. sc separate
#############
setwd("sc")

Idents(sc) = sc$celltype
DimPlot(sc, reduction = "umap", label = TRUE)
DimPlot(sc, reduction = "umap", group.by = "ecs", split.by = "diet")

# subset only ECs
Idents(sc) = sc$ecs
sc = subset(sc, idents = "Pecam1-or-Cdh5-pos")
saveRDS(sc, "sc_ECs.rds")

Idents(sc) = sc$celltype
DimPlot(sc, reduction = "umap", group.by = "diet")

# Endothelial subtypes Good Markers
FeaturePlot(sc, features = c("Vwf", "Vcam1", #  large vessels
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
FeaturePlot(sc, features = c("Flt1", "Pecam1","Cdh5",# EC
                                "Pdgfrb", "Cspg4","Anpep",  #  Pericyte
                                "Myh11", "Acta2", "Des",   # SMC
                                "Dcn", "Col1a1", "Pdgfra"), #FB
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10") #, min.cutoff = "q10"

FeaturePlot(sc, features = c("Prox1", "Pdpn","Lyve1", "Mki67"), # lymph
            order = T, cols=c("grey", "red"), ncol=2, min.cutoff = "q10") #, min.cutoff = "q10"



sc <- ScaleData(object = sc, verbose = T)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)
sc <- RunPCA(object = sc, npcs = 30, verbose = T)
ElbowPlot(object = sc,  ndims = 30) 

# UMAP and Clustering
sc <- RunUMAP(object = sc, reduction = "pca", dims = 1:28)
sc <- FindNeighbors(object = sc, reduction = "pca", dims = 1:28)
sc <- FindClusters(object = sc, resolution = 0.4)

# Visualization
DimPlot(object = sc, reduction = "umap", label = TRUE)
DimPlot(object = sc, reduction = "umap", group.by = "diet")
#DimPlot(object = sc, reduction = "umap", group.by = "organ", cols = c("purple", "orange"))
DimPlot(object = sc, reduction = "umap", group.by = "celltype", label = T)

# rename clusters
sc <- RenameIdents(object = sc,`8` = "Lymph-EC") 

# select cluster by hand
plot <- DimPlot(sc, reduction = "umap")
select.cells <- CellSelector(plot = plot) # 
Idents(sc, cells = select.cells) <- "Prolif"

## Relevel object@ident
sc@active.ident <- factor(x = sc@active.ident, levels = sort(levels(Idents(sc))))
sc$sc_celltype <- Idents(sc)
DimPlot(sc, reduction = "umap", label = TRUE)

# Find marker genes
Idents(sc) = sc$sc_celltype
sc.markers = FindAllMarkers(sc, only.pos = T)
write.table(sc.markers, "sc_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(sc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "sc_top5_cluster_markers.txt", sep='\t')
DotPlot(sc, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

# number of cells per cluster per sample
table(Idents(sc), sc$diet)

DimPlot(sc, reduction = "umap", label = TRUE, split.by = "diet", ncol = 2) + NoLegend()

# DE genes
sc$clust.diet <- paste(sc$seurat_clusters, sc$diet, sep = "_")
Idents(object = sc) <- "clust.diet"
sc@active.ident <- factor(x = sc@active.ident, levels = sort(levels(Idents(object = sc))))
levels(Idents(object = sc))

lev1 <- levels(sc@active.ident)[seq(2,(length(levels(Idents(object = sc)))),2)] # condition 1=HFD
lev2 <- levels(sc@active.ident)[seq(1,(length(levels(Idents(object = sc)))),2)] # condition 2=chow

# for all the clusters
for (i in seq(1,length(lev1),1)) { # 16 for not-filtered
  print(lev1[i])
  print(lev2[i])
  # TG 19m vs WT 19m
  response2 <- FindMarkers(object = sc, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = FALSE) #, min.pct = 0.25
  write.table(response2, file = paste(lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
}

# Endothelial subtypes Good Markers
FeaturePlot(sc, features = c("Vwf", "Vcam1", #  large vessels
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


saveRDS(sc, "sc_EC.rds")

###############
# B. vis separate
##############
setwd("vis")

Idents(vis) = vis$celltype
DimPlot(vis, reduction = "umap", label = TRUE)
DimPlot(vis, reduction = "umap", group.by = "ecs", split.by = "diet")

# subset only ECs
Idents(vis) = vis$ecs
vis = subset(vis, idents = "Pecam1-or-Cdh5-pos")
saveRDS(vis, "vis_ECs.rds")

Idents(vis) = vis$celltype
DimPlot(vis, reduction = "umap", group.by = "diet")

# Endothelial subtypes Good Markers
FeaturePlot(vis, features = c("Vwf", "Vcam1", #  large vessels
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
FeaturePlot(vis, features = c("Flt1", "Pecam1","Cdh5",# EC
                                "Pdgfrb", "Cspg4","Anpep",  #  Pericyte
                                "Myh11", "Acta2", "Des",   # SMC
                                "Dcn", "Col1a1", "Pdgfra"), #FB
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10") #, min.cutoff = "q10"

FeaturePlot(vis, features = c("Prox1", "Pdpn","Lyve1", "Mki67"), # lymph
            order = T, cols=c("grey", "red"), ncol=2, min.cutoff = "q10") #, min.cutoff = "q10"



vis <- ScaleData(object = vis, verbose = T)
vis <- FindVariableFeatures(vis, selection.method = "vst", nfeatures = 3000)
vis <- RunPCA(object = vis, npcs = 30, verbose = T)
ElbowPlot(object = vis,  ndims = 30) 

# UMAP and Clustering
vis <- RunUMAP(object = vis, reduction = "pca", dims = 1:25)
vis <- FindNeighbors(object = vis, reduction = "pca", dims = 1:25)
vis <- FindClusters(object = vis, resolution = 0.4)

# Visualization
DimPlot(object = vis, reduction = "umap", label = TRUE)
DimPlot(object = vis, reduction = "umap", group.by = "diet")
DimPlot(object = vis, reduction = "umap", group.by = "celltype", label = T)

# select cluster by hand
plot <- DimPlot(vis, reduction = "umap")
select.cells <- CellSelector(plot = plot) # 
Idents(vis, cells = select.cells) <- "Prolif"

vis <- RenameIdents(vis, `1` = "lymph-EC1",`4` = "lymph-EC2",`8` = "lymph-EC3",`9` = "lymph-EC4") 
## Relevel object@ident
vis@active.ident <- factor(x = vis@active.ident, levels = sort(levels(Idents(vis))))
vis$vis_celltype <- Idents(vis)
DimPlot(vis, reduction = "umap", label = TRUE)

# Find marker genes
Idents(vis) = vis$vis_celltype
vis.markers = FindAllMarkers(vis, only.pos = T)
write.table(vis.markers, "vis_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(vis.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "vis_top5_cluster_markers.txt", sep='\t')
DotPlot(vis, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

# number of cells per cluster per sample
table(Idents(vis), vis$diet)

DimPlot(vis, reduction = "umap", label = TRUE, split.by = "diet", ncol = 2) + NoLegend()

# DE genes
vis$clust.diet <- paste(vis$vis_celltype, vis$diet, sep = "_")
Idents(object = vis) <- "clust.diet"
vis@active.ident <- factor(x = vis@active.ident, levels = sort(levels(Idents(object = vis))))
levels(Idents(object = vis))

lev1 <- levels(vis@active.ident)[seq(2,(length(levels(Idents(object = vis)))),2)] # condition 1=HFD
lev2 <- levels(vis@active.ident)[seq(1,(length(levels(Idents(object = vis)))),2)] # condition 2=chow

# for all the clusters
for (i in seq(1,length(lev1),1)) { # 16 for not-filtered
  print(lev1[i])
  print(lev2[i])
  # TG 19m vs WT 19m
  response2 <- FindMarkers(object = vis, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = FALSE) #, min.pct = 0.25
  write.table(response2, file = paste(lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
}

# Endothelial subtypes Good Markers
FeaturePlot(vis, features = c("Vwf", "Vcam1", #  large vessels
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

#  tip cell-enriched, angiogenic
FeaturePlot(vis, features = c("Apln", "Col4a2", "Trp53i11"),
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 4)

Idents(vis) = vis$vis_celltype
DimPlot(object = vis, reduction = "umap", label = TRUE) + NoLegend()
vis <- RenameIdents(vis, `3` = "EC-ven",`0` = "EC-cap1",`2` = "EC-cap2",`5` = "EC-art",
                    `7` = "EC-ang") 
## Relevel object@ident
vis@active.ident <- factor(x = vis@active.ident, levels = sort(levels(Idents(vis))))
vis$vis_celltype <- Idents(vis)
DimPlot(vis, reduction = "umap", label = TRUE)

saveRDS(vis, "vis_EC.rds")

#########################
# 3. brain
#########################
brain = readRDS("../brain_subset.rds")
setwd("../brain")

Idents(brain) = brain$celltype
DimPlot(brain, reduction = "umap", label = TRUE)
DimPlot(brain, reduction = "umap", group.by = "ecs", split.by = "diet")

# subset only ECs
Idents(brain) = brain$ecs
brain = subset(brain, idents = "Pecam1-or-Cdh5-pos")
saveRDS(brain, "brain_ECs.rds")


# Pericyte amd SMC markers
FeaturePlot(brain, features = c("Flt1", "Pecam1","Cdh5",# EC
                              "Pdgfrb", "Cspg4","Anpep",  #  Pericyte
                              "Myh11", "Acta2", "Des",   # SMC
                              "Dcn", "Col1a1", "Pdgfra"), #FB
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10") #, min.cutoff = "q10"

FeaturePlot(brain, features = c("Prox1", "Pdpn","Lyve1", "Mki67"), # lymph
            order = T, cols=c("grey", "red"), ncol=2, min.cutoff = "q10") #, min.cutoff = "q10"


brain <- ScaleData(object = brain, verbose = T)
brain <- FindVariableFeatures(brain, selection.method = "vst", nfeatures = 3000)
brain <- RunPCA(object = brain, npcs = 30, verbose = T)
ElbowPlot(object = brain,  ndims = 30) 

# UMAP and Clustering
brain <- RunUMAP(object = brain, reduction = "pca", dims = 1:25)
brain <- FindNeighbors(object = brain, reduction = "pca", dims = 1:25)
brain <- FindClusters(object = brain, resolution = 0.4)

# visualization
DimPlot(object = brain, reduction = "umap", label = TRUE)
DimPlot(object = brain, reduction = "umap", group.by = "diet")
DimPlot(object = brain, reduction = "umap", group.by = "celltype", label = T)

# select cluster by hand
#plot <- DimPlot(brain, reduction = "umap")
#select.cells <- CellSelector(plot = plot) # 
#Idents(brain, cells = select.cells) <- "Prolif"

# Find marker genes
brain.markers = FindAllMarkers(brain, only.pos = T)
write.table(brain.markers, "brain_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(brain.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "brain_top5_cluster_markers.txt", sep='\t')
DotPlot(brain, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

# Endothelial subtypes Good Markers
FeaturePlot(brain, features = c("Vwf", "Vcam1", #  large vessels
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

#  tip cell-enriched, angiogenic
FeaturePlot(brain, features = c("Apln", "Col4a2", "Trp53i11"),
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 3)

DimPlot(object = brain, reduction = "umap", label = TRUE) + NoLegend()
brain <- RenameIdents(brain, `2` = "EC-ven",`0` = "EC-cap1",`2` = "EC-cap2",`9` = "EC-art",
                      `1` = "EC-capA1",`3` = "EC-capA2",
                      `10` = "EC-fenestr",`7` = "EC-Hb",`5` = "EC-AP1",`4` = "EC-IFN") 
## Relevel object@ident
brain@active.ident <- factor(x = brain@active.ident, levels = sort(levels(Idents(brain))))
brain$brain_celltype <- Idents(brain)
DimPlot(brain, reduction = "umap", label = TRUE)

# number of cells per cluster per sample
table(Idents(brain), brain$diet)

DimPlot(brain, reduction = "umap", label = TRUE, split.by = "diet", ncol = 2) + NoLegend()

# DE genes
brain$clust.diet <- paste(brain$brain_celltype, brain$diet, sep = "_")
Idents(object = brain) <- "clust.diet"
brain@active.ident <- factor(x = brain@active.ident, levels = sort(levels(Idents(object = brain))))
levels(Idents(object = brain))

lev1 <- levels(brain@active.ident)[seq(2,(length(levels(Idents(object = brain)))),2)] # condition 1=HFD
lev2 <- levels(brain@active.ident)[seq(1,(length(levels(Idents(object = brain)))),2)] # condition 2=chow

# for all the clusters
for (i in seq(1,length(lev1),1)) { # 16 for not-filtered
  print(lev1[i])
  print(lev2[i])
  # TG 19m vs WT 19m
  response2 <- FindMarkers(object = brain, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = FALSE) #, min.pct = 0.25
  write.table(response2, file = paste(lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
}


saveRDS(brain, "brain_ECs.rds")

#########################
# 4. heart
#########################
heart = readRDS("../heart_subset.rds")
setwd("../heart")

Idents(heart) = heart$celltype
DimPlot(heart, reduction = "umap", label = TRUE)
DimPlot(heart, reduction = "umap", group.by = "ecs", split.by = "diet")

# subset only ECs
Idents(heart) = heart$ecs
heart = subset(heart, idents = "Pecam1-or-Cdh5-pos")
saveRDS(heart, "heart_ECs.rds")


# Pericyte amd SMC markers
FeaturePlot(heart, features = c("Flt1", "Pecam1","Cdh5",# EC
                                "Pdgfrb", "Cspg4","Anpep",  #  Pericyte
                                "Myh11", "Acta2", "Des",   # SMC
                                "Dcn", "Col1a1", "Pdgfra"), #FB
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10") #, min.cutoff = "q10"

FeaturePlot(heart, features = c("Prox1", "Pdpn","Lyve1", "Mki67"), # lymph
            order = T, cols=c("grey", "red"), ncol=2, min.cutoff = "q10") #, min.cutoff = "q10"


heart <- ScaleData(object = heart, verbose = T)
heart <- FindVariableFeatures(heart, selection.method = "vst", nfeatures = 3000)
heart <- RunPCA(object = heart, npcs = 30, verbose = T)
ElbowPlot(object = heart,  ndims = 30) 

# UMAP and Clustering
heart <- RunUMAP(object = heart, reduction = "pca", dims = 1:25)
heart <- FindNeighbors(object = heart, reduction = "pca", dims = 1:25)
heart <- FindClusters(object = heart, resolution = 0.4)

# visualization
DimPlot(object = heart, reduction = "umap", label = TRUE)
DimPlot(object = heart, reduction = "umap", group.by = "diet")
DimPlot(object = heart, reduction = "umap", group.by = "celltype", label = T)
DimPlot(object = heart, reduction = "umap", group.by = "heart_celltype", 
        split.by = "orig.ident",label = T, ncol = 2)
table(Idents(heart), heart$orig.ident)

# select cluster by hand
plot <- DimPlot(heart, reduction = "umap")
select.cells <- CellSelector(plot = plot) # 
Idents(heart, cells = select.cells) <- "11"

plot <- DimPlot(heart, reduction = "umap")
select.cells <- CellSelector(plot = plot) # 
Idents(heart, cells = select.cells) <- "12"
#heart = subset(heart, idents = "outlier", invert = T)
heart$seurat_clusters = Idents(heart)

# Find marker genes
heart.markers = FindAllMarkers(heart, only.pos = T)
write.table(heart.markers, "heart_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(heart.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "heart_top5_cluster_markers.txt", sep='\t')
DotPlot(heart, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

# Endothelial subtypes Good Markers
FeaturePlot(heart, features = c("Vwf", "Vcam1", #  large vessels
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

#  tip cell-enriched, angiogenic
FeaturePlot(heart, features = c("Apln", "Col4a2", "Trp53i11"),
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 3)

DimPlot(object = heart, reduction = "umap", label = TRUE) + NoLegend()
heart <- RenameIdents(heart, `4` = "EC-ven",`0` = "EC-cap1",`1` = "EC-cap2",`3` = "EC-art",
                      `9` = "EC-ecm", `8` = "EC-Hb",
                      `7` = "EC-lymph", `10` = "Prolif",`2` = "EC-AP1",`5` = "EC-IFN") 
## Relevel object@ident
heart@active.ident <- factor(x = heart@active.ident, levels = sort(levels(Idents(heart))))
heart$heart_celltype <- Idents(heart)
DimPlot(heart, reduction = "umap", label = TRUE) + NoLegend()

# number of cells per cluster per sample
table(Idents(heart), heart$diet)

DimPlot(heart, reduction = "umap", label = TRUE, split.by = "diet", ncol = 2) + NoLegend()

# DE genes
heart$clust.diet <- paste(heart$heart_celltype, heart$diet, sep = "_")
Idents(object = heart) <- "clust.diet"
heart@active.ident <- factor(x = heart@active.ident, levels = sort(levels(Idents(object = heart))))
levels(Idents(object = heart))

lev1 <- levels(heart@active.ident)[seq(2,(length(levels(Idents(object = heart)))),2)] # condition 1=HFD
lev2 <- levels(heart@active.ident)[seq(1,(length(levels(Idents(object = heart)))),2)] # condition 2=chow

# for all the clusters
for (i in seq(1,length(lev1),1)) { # 16 for not-filtered
  print(lev1[i])
  print(lev2[i])
  # TG 19m vs WT 19m
  response2 <- FindMarkers(object = heart, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = FALSE) #, min.pct = 0.25
  write.table(response2, file = paste(lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
}


saveRDS(heart, "heart_ECs.rds")

#########################
# 5. kidney
#########################
kidney = readRDS("kidney_subset.rds")
setwd("kidney")

Idents(kidney) = kidney$celltype
DimPlot(kidney, reduction = "umap", label = TRUE)
DimPlot(kidney, reduction = "umap", group.by = "ecs", split.by = "diet")

# subset only ECs
Idents(kidney) = kidney$ecs
kidney = subset(kidney, idents = "Pecam1-or-Cdh5-pos")
saveRDS(kidney, "kidney_ECs.rds")


# Pericyte amd SMC markers
FeaturePlot(kidney, features = c("Flt1", "Pecam1","Cdh5",# EC
                                "Pdgfrb", "Cspg4","Anpep",  #  Pericyte
                                "Myh11", "Acta2", "Des",   # SMC
                                "Dcn", "Col1a1", "Pdgfra"), #FB
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10") #, min.cutoff = "q10"

FeaturePlot(kidney, features = c("Prox1", "Pdpn","Lyve1", "Mki67"), # lymph
            order = T, cols=c("grey", "red"), ncol=2, min.cutoff = "q10") #, min.cutoff = "q10"


kidney <- ScaleData(object = kidney, verbose = T)
kidney <- FindVariableFeatures(kidney, selection.method = "vst", nfeatures = 3000)
kidney <- RunPCA(object = kidney, npcs = 30, verbose = T)
ElbowPlot(object = kidney,  ndims = 30) 

# UMAP and Clustering
kidney <- RunUMAP(object = kidney, reduction = "pca", dims = 1:29)
kidney <- FindNeighbors(object = kidney, reduction = "pca", dims = 1:29)
kidney <- FindClusters(object = kidney, resolution = 0.4)

# visualization
DimPlot(object = kidney, reduction = "umap", label = TRUE)
DimPlot(object = kidney, reduction = "umap", group.by = "diet")
DimPlot(object = kidney, reduction = "umap", group.by = "celltype", label = T)

# select cluster by hand
plot <- DimPlot(kidney, reduction = "umap")
select.cells <- CellSelector(plot = plot) # 
Idents(kidney, cells = select.cells) <- "12"
#kidney = subset(kidney, idents = "outlier", invert = T)
kidney$RNA_snn_res.0.4 = Idents(kidney)

# Find marker genes
kidney.markers = FindAllMarkers(kidney, only.pos = T)
write.table(kidney.markers, "kidney_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(kidney.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "kidney_top5_cluster_markers.txt", sep='\t')
DotPlot(kidney, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

# Endothelial subtypes Good Markers
FeaturePlot(kidney, features = c("Vwf", "Vcam1", #  large vessels
                                "Fbln5", "Cytl1", # large artery
                                "Hey1", "Gkn3",# artery 
                                "Tgfb2", "Glul", # CapA
                                "Rgcc", "Mfsd2a", # Cap
                                "Car4", "Tfrc", # CapV
                                "Lcn2", #"Slc38a5", # large vein
                                "Hba-a1", "Hbb-bs", 
                                "Plvap", "Plpp3",
                                "Flt4", "Ccl21a"),
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 4) #

#  tip cell-enriched, angiogenic
FeaturePlot(kidney, features = c("Apln", "Col4a2", "Trp53i11"),
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 3)

DimPlot(object = kidney, reduction = "umap", label = TRUE) + NoLegend()
kidney <- RenameIdents(kidney, `6` = "EC-ven",`0` = "EC-cap1",`1` = "EC-cap2",`8` = "EC-art",
                      `9` = "EC-AP1", `7` = "mEC1", `10` = "mEC2",  # m = medular
                      `11` = "EC-glomerular", `12` = "EC-lymph", `3` = "EC-Aqp1") 
# select cluster by hand
plot <- DimPlot(kidney, reduction = "umap", cells = WhichCells(kidney, idents = "4"))
select.cells <- CellSelector(plot = plot) # 
Idents(kidney, cells = select.cells) <- "EC-ang"

## Relevel object@ident
kidney@active.ident <- factor(x = kidney@active.ident, levels = sort(levels(Idents(kidney))))
kidney$kidney_celltype <- Idents(kidney)
DimPlot(kidney, reduction = "umap", label = TRUE) + NoLegend()

# number of cells per cluster per sample
table(Idents(kidney), kidney$diet)

DimPlot(kidney, reduction = "umap", label = TRUE, split.by = "diet", ncol = 2) + NoLegend()

# DE genes
kidney$clust.diet <- paste(kidney$kidney_celltype, kidney$diet, sep = "_")
Idents(object = kidney) <- "clust.diet"
kidney@active.ident <- factor(x = kidney@active.ident, levels = sort(levels(Idents(object = kidney))))
levels(Idents(object = kidney))

lev1 <- levels(kidney@active.ident)[seq(2,(length(levels(Idents(object = kidney)))),2)] # condition 1=HFD
lev2 <- levels(kidney@active.ident)[seq(1,(length(levels(Idents(object = kidney)))),2)] # condition 2=chow

# for all the clusters
for (i in seq(1,length(lev1),1)) { # 16 for not-filtered
  print(lev1[i])
  print(lev2[i])
  # TG 19m vs WT 19m
  response2 <- FindMarkers(object = kidney, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = FALSE) #, min.pct = 0.25
  write.table(response2, file = paste(lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
}


saveRDS(kidney, "kidney_ECs.rds")

#########################
# 6. lung
#########################
lung = readRDS("../lung_subset.rds")
setwd("../lung")

Idents(lung) = lung$celltype
DimPlot(lung, reduction = "umap", label = TRUE)
DimPlot(lung, reduction = "umap", group.by = "ecs", split.by = "diet")

# subset only ECs
Idents(lung) = lung$ecs
lung = subset(lung, idents = "Pecam1-or-Cdh5-pos")
saveRDS(lung, "lung_ECs.rds")


# Pericyte amd SMC markers
FeaturePlot(lung, features = c("Flt1", "Pecam1","Cdh5",# EC
                                "Pdgfrb", "Cspg4","Anpep",  #  Pericyte
                                "Myh11", "Acta2", "Des",   # SMC
                                "Dcn", "Col1a1", "Pdgfra"), #FB
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10") #, min.cutoff = "q10"

FeaturePlot(lung, features = c("Prox1", "Pdpn","Lyve1", "Mki67"), # lymph
            order = T, cols=c("grey", "red"), ncol=2, min.cutoff = "q10") #, min.cutoff = "q10"


lung <- ScaleData(object = lung, verbose = T)
lung <- FindVariableFeatures(lung, selection.method = "vst", nfeatures = 3000)
lung <- RunPCA(object = lung, npcs = 30, verbose = T)
ElbowPlot(object = lung,  ndims = 30) 

# UMAP and Clustering
lung <- RunUMAP(object = lung, reduction = "pca", dims = 1:25)
lung <- FindNeighbors(object = lung, reduction = "pca", dims = 1:25)
lung <- FindClusters(object = lung, resolution = 0.4)

# visualization
DimPlot(object = lung, reduction = "umap", label = TRUE)
DimPlot(object = lung, reduction = "umap", group.by = "diet")
DimPlot(object = lung, reduction = "umap", group.by = "celltype", label = T)

# select cluster by hand
#plot <- DimPlot(lung, reduction = "umap")
#select.cells <- CellSelector(plot = plot) # 
#Idents(lung, cells = select.cells) <- "11"

#plot <- DimPlot(lung, reduction = "umap")
#select.cells <- CellSelector(plot = plot) # 
#Idents(lung, cells = select.cells) <- "12"
#lung = subset(lung, idents = "outlier", invert = T)
#lung$seurat_clusters = Idents(lung)

# Find marker genes
lung.markers = FindAllMarkers(lung, only.pos = T)
write.table(lung.markers, "lung_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(lung.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "lung_top5_cluster_markers.txt", sep='\t')
DotPlot(lung, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

# Endothelial subtypes Good Markers
FeaturePlot(lung, features = c("Vwf", "Vcam1", #  large vessels
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

#  tip cell-enriched, angiogenic
FeaturePlot(lung, features = c("Apln", "Col4a2", "Trp53i11"),
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 3)

DimPlot(object = lung, reduction = "umap", label = TRUE) + NoLegend()
lung <- RenameIdents(lung, `1` = "EC-ven",`0` = "EC-cap1",`2` = "EC-cap2",`6` = "EC-art",
                     `8` = "aEC1",`3` = "aEC2",
                     `4` = "pulmEC-a",`10` = "pulmEC-b",`14` = "EC-platelet",
                     `5` = "EC-Aqp5a",`11` = "EC-Aqp5b" ,`9` = "EC-pneumocyte",
                      `7` = "EC-lymph", `16` = "Prolif") 

## Relevel object@ident
lung@active.ident <- factor(x = lung@active.ident, levels = sort(levels(Idents(lung))))
lung$lung_celltype <- Idents(lung)
DimPlot(lung, reduction = "umap", label = TRUE) + NoLegend()

# number of cells per cluster per sample
table(Idents(lung), lung$diet)

DimPlot(lung, reduction = "umap", label = TRUE, split.by = "diet", ncol = 2) + NoLegend()

# DE genes
lung$clust.diet <- paste(lung$lung_celltype, lung$diet, sep = "_")
Idents(object = lung) <- "clust.diet"
lung@active.ident <- factor(x = lung@active.ident, levels = sort(levels(Idents(object = lung))))
levels(Idents(object = lung))

lev1 <- levels(lung@active.ident)[seq(2,(length(levels(Idents(object = lung)))),2)] # condition 1=HFD
lev2 <- levels(lung@active.ident)[seq(1,(length(levels(Idents(object = lung)))),2)] # condition 2=chow

# for all the clusters
for (i in seq(1,length(lev1),1)) { # 16 for not-filtered
  print(lev1[i])
  print(lev2[i])
  # TG 19m vs WT 19m
  response2 <- FindMarkers(object = lung, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = FALSE) #, min.pct = 0.25
  write.table(response2, file = paste(lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
}


saveRDS(lung, "lung_ECs.rds")


#########################
# 7. lymphatics
#########################
lymph = readRDS("lymphatic_subset.rds")
setwd("lymph")

Idents(lymph) = lymph$celltype
DimPlot(lymph, reduction = "umap", label = TRUE)
DimPlot(lymph, reduction = "umap", group.by = "ecs", split.by = "diet")

# subset only ECs
Idents(lymph) = lymph$ecs
lymph = subset(lymph, idents = "Pecam1-or-Cdh5-pos")
saveRDS(lymph, "lymph_ECs.rds")


# Pericyte amd SMC markers
FeaturePlot(lymph, features = c("Flt1", "Pecam1","Cdh5",# EC
                               "Pdgfrb", "Cspg4","Anpep",  #  Pericyte
                               "Myh11", "Acta2", "Des",   # SMC
                               "Dcn", "Col1a1", "Pdgfra"), #FB
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10") #, min.cutoff = "q10"

FeaturePlot(lymph, features = c("Prox1", "Pdpn","Lyve1", "Mki67"), # lymph
            order = T, cols=c("grey", "red"), ncol=2, min.cutoff = "q10") #, min.cutoff = "q10"


lymph <- ScaleData(object = lymph, verbose = T)
lymph <- FindVariableFeatures(lymph, selection.method = "vst", nfeatures = 3000)
lymph <- RunPCA(object = lymph, npcs = 30, verbose = T)
ElbowPlot(object = lymph,  ndims = 30) 

# UMAP and Clustering
lymph <- RunUMAP(object = lymph, reduction = "pca", dims = 1:25)
lymph <- FindNeighbors(object = lymph, reduction = "pca", dims = 1:25)
lymph <- FindClusters(object = lymph, resolution = 0.4)

# visualization
DimPlot(object = lymph, reduction = "umap", label = TRUE)
DimPlot(object = lymph, reduction = "umap", group.by = "diet")
DimPlot(object = lymph, reduction = "umap", group.by = "organ", label = T)

# Find marker genes
lymph.markers = FindAllMarkers(lymph, only.pos = T)
write.table(lymph.markers, "lymph_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(lymph.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "lymph_top5_cluster_markers.txt", sep='\t')
DotPlot(lymph, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

#  tip cell-enriched, angiogenic
FeaturePlot(lymph, features = c("Apln", "Col4a2", "Trp53i11"),
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 3)

DimPlot(object = lymph, reduction = "umap", label = TRUE) + NoLegend()
# rename
lymph <- RenameIdents(lymph,  `5` = "EC-FA-transp",`6` = "EC-Pltp", 
                     `7` = "EC-antigen", `8` = "EC-Nfat5") 

## Relevel object@ident
lymph@active.ident <- factor(x = lymph@active.ident, levels = sort(levels(Idents(lymph))))
lymph$lymph_celltype <- Idents(lymph)
DimPlot(lymph, reduction = "umap", label = TRUE) + NoLegend()

# number of cells per cluster per sample
table(Idents(lymph), lymph$diet)
table(Idents(lymph), lymph$organ)

DimPlot(lymph, reduction = "umap", label = TRUE, split.by = "diet", ncol = 2) + NoLegend()

# DE genes
lymph$clust.diet <- paste(lymph$lymph_celltype, lymph$diet, sep = "_")
Idents(object = lymph) <- "clust.diet"
lymph@active.ident <- factor(x = lymph@active.ident, levels = sort(levels(Idents(object = lymph))))
levels(Idents(object = lymph))

lev1 <- levels(lymph@active.ident)[seq(2,(length(levels(Idents(object = lymph)))),2)] # condition 1=HFD
lev2 <- levels(lymph@active.ident)[seq(1,(length(levels(Idents(object = lymph)))),2)] # condition 2=chow

# for all the clusters
for (i in seq(1,length(lev1),1)) { # 16 for not-filtered
  print(lev1[i])
  print(lev2[i])
  # TG 19m vs WT 19m
  response2 <- FindMarkers(object = lymph, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = FALSE) #, min.pct = 0.25
  write.table(response2, file = paste(lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
}

# pseudo-bulk
Idents(object = lymph) <- "diet"
response2 <- FindMarkers(object = lymph, ident.1 = "Western", 
                         ident.2 = "chow", 
                         logfc.threshold = 0,
                         verbose = FALSE) #, min.pct = 0.25
write.table(response2, file ="lymph_all_WD_vs_chow.txt", sep = "\t")

saveRDS(lymph, "lymph_ECs.rds")

#####################################
# 8. Compare all organs by Venous / Capillary / Arterial DEGs
#####################################

setwd('comparisons')

brain = readRDS("../brain/brain_ECs.rds")
lung = readRDS("../lung/lung_ECs.rds")
heart = readRDS("../heart/heart_ECs.rds")
liver = readRDS("../liver/liver_ECs.rds")
kidney = readRDS("../kidney/kidney_ECs.rds")
vis = readRDS("../vis/vis_ECs.rds")
sc = readRDS("../sc/sc_ECs.rds")

col.gen = c("#e30026","#fcacac","#4287f5","#e02d8a","#46d130", "#f09000", "#bed647","#47d6c1",
            "#3037c7", "#006887","#87b5ff","#ff5e5e")
names(col.gen) = c("Art","Cap","Ven","Prolif","Lymph","Angiog","Sp2","Sp3","Sp4","Sp5","Venule","Arteriole")
col.gen  
  
##############
# 1. liver
Idents(liver) = liver$seurat_clusters
DimPlot(object = liver, reduction = "umap", label = TRUE) + NoLegend()
# Endothelial subtypes Good Markers
FeaturePlot(liver, features = c("Vwf", "Vcam1", #  large vessels
                                   "Fbln5", "Cytl1", # large artery
                                   "Hey1", "Gkn3",# artery 
                                   "Tgfb2", "Glul", # CapA
                                   "Rgcc", "Mfsd2a", # Cap
                                   "Car4", "Tfrc", # CapV
                                   "Lcn2", "Slc38a5", # large vein
                                   "Hba-a1", "Hbb-bs", 
                                   "Plvap", "Plpp3", "Mki67",
                                    "Ccl21a","Prox1","Lyve1",  # lymph
                                   "Apln", "Col4a2", "Trp53i11"), # angiogen
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 5) #

# select cluster by hand
plot <- DimPlot(liver, reduction = "umap")
select.cells <- CellSelector(plot = plot) # 
Idents(liver, cells = select.cells) <- "EC-art"

# rename
liver <- RenameIdents(liver,  `4` = "EC-ven",`0` = "EC-cap", `1` = "EC-cap", `2` = "EC-cap", 
                      `6` = "Prolif", `3` = "EC-liver1", `5` = "EC-liver2") 

## Relevel object@ident
liver@active.ident <- factor(x = liver@active.ident, levels = sort(levels(Idents(liver))))
liver$gen_celltype <- Idents(liver)
col.gen
DimPlot(liver, reduction = "umap", label = TRUE, cols = as.vector(col.gen[c(1,2,6,7,3,4)]), pt.size = 1)

# number of cells per cluster per sample
table(Idents(liver), liver$diet)

# Find marker genes
Idents(liver) = liver$gen_celltype
liver.markers = FindAllMarkers(liver, only.pos = T)
write.table(liver.markers, "liver_gen_markers.txt", sep='\t')

top5.mark <- as.data.frame(liver.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "liver_top5_gen_markers.txt", sep='\t')
DotPlot(liver, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

# Liver-specific Capillaries : Stab2 gene (liver sinusoidal capillaries)
FeaturePlot(liver, features = c("Vwf", "Vcam1", #  large vessels
                                 "Rspo3","Lhx6", # artery??
                                "Fbln5","Hey1", # artery 
                                "Stab2", # liver sinusoidal capillaries
                                "Adgrg6", "Ednrb"), # venous? 
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 3) #

# number of cells per cluster per sample
table(Idents(liver), liver$diet)

# DE genes
liver$gen.diet <- paste(liver$gen_celltype, liver$diet, sep = "_")
Idents(object = liver) <- "gen.diet"
liver@active.ident <- factor(x = liver@active.ident, levels = sort(levels(Idents(object = liver))))
levels(Idents(object = liver))

lev1 <- levels(liver@active.ident)[seq(2,(length(levels(Idents(object = liver)))),2)] # condition 1=HFD
lev2 <- levels(liver@active.ident)[seq(1,(length(levels(Idents(object = liver)))),2)] # condition 2=chow

# for all the clusters
for (i in seq(1,length(lev1),1)) { # 16 for not-filtered
  print(lev1[i])
  print(lev2[i])
  # TG 19m vs WT 19m
  response2 <- FindMarkers(object = liver, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = FALSE) #, min.pct = 0.25
  write.table(response2, file = paste("liver_",lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
}

# subclusters
DimPlot(liver, reduction = "umap", label = TRUE, cols = as.vector(col.gen[c(1,2,6,7,3,4)]), pt.size = 1)
DimPlot(liver, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1)

# rename
ac = WhichCells(liver, idents = "EC-art")
Idents(liver) = liver$seurat_clusters
liver <- RenameIdents(liver,  `4` = "EC-ven",`0` = "EC-cap1", `1` = "EC-cap2", `2` = "EC-cap3", 
                      `6` = "Prolif", `3` = "EC-liver1", `5` = "EC-liver2") 
Idents(liver, cells = ac) <- "EC-art"
levels(Idents(liver))
## Relevel object@ident
liver@active.ident <- factor(x = liver@active.ident, levels = sort(levels(Idents(liver))))
liver$liver_celltype <- Idents(liver)
levels(Idents(liver))
col.gen2 = c("#e30026","#fcacac","#ff7d61","#ffb5db","#cf91cc","#cf91b1",
             "#4287f5","#e02d8a","#46d130", 
             "#f09000", "#bed647","#47d6c1",
            "#3037c7", "#006887","#87b5ff","#ff5e5e")
names(col.gen2) = c("Art","Cap","Cap","Cap","Cap","Cap",
                    "Ven","Prolif","Lymph",
                    "Angiog","Sp2","Sp3",
                    "Sp4","Sp5","Venule","Arteriole")
col.gen2  

DimPlot(liver, reduction = "umap", label = TRUE, cols = as.vector(col.gen2[c(1,2,3,4,10,11,7,8)]), pt.size = 1)

# number of cells per cluster per sample
table(Idents(liver), liver$diet)

saveRDS(liver, "liver_ECs.rds")

##############
# 2. vis fat 
Idents(vis) = vis$vis_celltype
DimPlot(object = vis, reduction = "umap", label = TRUE) + NoLegend()
# Endothelial subtypes Good Markers
FeaturePlot(vis, features = c("Vwf", "Vcam1", #  large vessels
                                "Fbln5", "Cytl1", # large artery
                                "Hey1", "Gkn3",# artery 
                                "Tgfb2", "Glul", # CapA
                                "Rgcc", "Mfsd2a", # Cap
                                "Car4", "Tfrc", # CapV
                                "Lcn2", "Slc38a5", # large vein
                                "Hba-a1", "Hbb-bs", 
                                "Plvap", "Plpp3", "Mki67",
                                "Ccl21a","Prox1","Lyve1",  # lymph
                                "Apln", "Col4a2", "Trp53i11"), # angiogen
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 5) #

# rename
vis <- RenameIdents(vis,  `6` = "EC-venule",`10` = "EC-vis", 
                    `EC-cap1` = "EC-cap", `EC-cap2` = "EC-cap",
                    `lymph-EC1` = "EC-lymph", `lymph-EC2` = "EC-lymph") 

## Relevel object@ident
vis@active.ident <- factor(x = vis@active.ident, levels = sort(levels(Idents(vis))))
vis$gen_celltype <- Idents(vis)
col.gen
levels(Idents(vis))
DimPlot(vis, reduction = "umap", label = TRUE, cols = as.vector(col.gen[c(6,7,1,2,5,3,11,8,4)]), pt.size = 1)

# number of cells per cluster per sample
table(Idents(vis), vis$diet)

# Find marker genes
Idents(vis) = vis$gen_celltype
vis.markers = FindAllMarkers(vis, only.pos = T)
write.table(vis.markers, "vis_gen_markers.txt", sep='\t')

top5.mark <- as.data.frame(vis.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "vis_top5_gen_markers.txt", sep='\t')
DotPlot(vis, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

# vis-specific Capillaries : Rgcc, Car4, Cd36
FeaturePlot(vis, features = c("Vwf", "Cfh", "Cpe",#  vein
                                "Cytl1","Gkn3","Eln", # artery 
                                "Rgcc","Car4","Cd36", # vis capillaries
                                "Flt4", "Ccl21a"), # venous? 
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 3) #

# number of cells per cluster per sample
table(Idents(vis), vis$diet)

# DE genes
vis$gen.diet <- paste(vis$gen_celltype, vis$diet, sep = "_")
Idents(object = vis) <- "gen.diet"
vis@active.ident <- factor(x = vis@active.ident, levels = sort(levels(Idents(object = vis))))
levels(Idents(object = vis))

lev1 <- levels(vis@active.ident)[seq(2,(length(levels(Idents(object = vis)))),2)] # condition 1=HFD
lev2 <- levels(vis@active.ident)[seq(1,(length(levels(Idents(object = vis)))),2)] # condition 2=chow

# for all the clusters
for (i in seq(1,length(lev1),1)) { # 16 for not-filtered
  print(lev1[i])
  print(lev2[i])
  # TG 19m vs WT 19m
  response2 <- FindMarkers(object = vis, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = FALSE) #, min.pct = 0.25
  write.table(response2, file = paste("vis_",lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
}

# subclusters
Idents(vis) = vis$vis_celltype
DimPlot(object = vis, reduction = "umap", label = TRUE) + NoLegend()
DimPlot(object = vis, reduction = "umap",group.by = "gen_celltype", label = TRUE) + NoLegend()

# rename
vis <- RenameIdents(vis,  `6` = "EC-venule",`10` = "EC-vis") 

## Relevel object@ident
vis@active.ident <- factor(x = vis@active.ident, levels = sort(levels(Idents(vis))))
vis$vis_celltype <- Idents(vis)

levels(Idents(vis))
 
DimPlot(vis, reduction = "umap", label = TRUE, cols = as.vector(col.gen2[c(10,11,1,2,3,7,15,12,9,17,18,19,8)]), pt.size = 1)
# 7.6 x 7
table(Idents(vis), vis$diet)

saveRDS(vis, "vis_ECs.rds")

###############
# 3. sc fat

Idents(sc) = sc$sc_celltype
DimPlot(object = sc, reduction = "umap", label = TRUE) + NoLegend()
# Endothelial subtypes Good Markers
FeaturePlot(sc, features = c("Vwf", "Vcam1", #  large vessels
                              "Fbln5", "Cytl1", # large artery
                              "Hey1", "Gkn3",# artery 
                              "Tgfb2", "Glul", # CapA
                              "Rgcc", "Mfsd2a", # Cap
                              "Car4", "Tfrc", # CapV
                              "Lcn2", "Slc38a5", # large vein
                              "Hba-a1", "Hbb-bs", 
                              "Plvap", "Plpp3", "Mki67",
                              "Ccl21a","Prox1","Lyve1",  # lymph
                              "Apln", "Col4a2", "Trp53i11"), # angiogen
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 5) #

# rename
sc <- RenameIdents(sc, `4` = "EC-art", `9` = "EC-ven",`5` = "EC-venule",
                   `0` = "EC-cap",`1` = "EC-cap",`2` = "EC-cap",`3` = "EC-cap",`7` = "EC-cap",`6` = "EC-cap",
                   `10` = "EC-Hb", `11` = "EC-antigen",
                    `Lymph-EC` = "EC-lymph") 

## Relevel object@ident
sc@active.ident <- factor(x = sc@active.ident, levels = sort(levels(Idents(sc))))
sc$gen_celltype <- Idents(sc)
col.gen
levels(Idents(sc))
DimPlot(sc, reduction = "umap", label = TRUE, cols = as.vector(col.gen[c(7,1,2,6,5,3,11,4)]), pt.size = 1)

# number of cells per cluster per sample
table(Idents(sc), sc$diet)

# Find marker genes
Idents(sc) = sc$gen_celltype
sc.markers = FindAllMarkers(sc, only.pos = T)
write.table(sc.markers, "sc_gen_markers.txt", sep='\t')

top5.mark <- as.data.frame(sc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "sc_top5_gen_markers.txt", sep='\t')
DotPlot(sc, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

# sc-specific Capillaries : Rgcc, Car4, Cd36
FeaturePlot(sc, features = c("Vwf", "Cfh", "Cpe",#  vein
                              "Clu","Gkn3","Fn1", # artery 
                              "Rgcc","Car4","Cd36", # sc capillaries
                              "Flt4", "Ccl21a"), # venous? 
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 3) #

# number of cells per cluster per sample
table(Idents(sc), sc$diet)

# DE genes
sc$gen.diet <- paste(sc$gen_celltype, sc$diet, sep = "_")
Idents(object = sc) <- "gen.diet"
sc@active.ident <- factor(x = sc@active.ident, levels = sort(levels(Idents(object = sc))))
levels(Idents(object = sc))

lev1 <- levels(sc@active.ident)[seq(2,(length(levels(Idents(object = sc)))),2)] # condition 1=HFD
lev2 <- levels(sc@active.ident)[seq(1,(length(levels(Idents(object = sc)))),2)] # condition 2=chow

# for all the clusters
for (i in seq(1,length(lev1),1)) { # 16 for not-filtered
  print(lev1[i])
  print(lev2[i])
  # TG 19m vs WT 19m
  response2 <- FindMarkers(object = sc, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = FALSE) #, min.pct = 0.25
  write.table(response2, file = paste("sc_",lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
}

# subclusters
Idents(sc) = sc$sc_celltype
DimPlot(object = sc, reduction = "umap", label = TRUE) + NoLegend()
DimPlot(object = sc, reduction = "umap",group.by = "gen_celltype", label = TRUE) + NoLegend()

# rename
sc <- RenameIdents(sc, `4` = "EC-art", `9` = "EC-ven",`5` = "EC-venule",
                   `0` = "EC-cap1",`1` = "EC-cap2",`2` = "EC-cap3",`3` = "EC-cap4",`7` = "EC-cap5",`6` = "EC-cap6",
                   `10` = "EC-Hb", `Lymph-EC` = "EC-lymph") 

## Relevel object@ident
sc@active.ident <- factor(x = sc@active.ident, levels = sort(levels(Idents(sc))))
sc$sc_celltype <- Idents(sc)

levels(Idents(sc))
DimPlot(sc, reduction = "umap", label = TRUE, cols = as.vector(col.gen2[c(11,1,2,16,3,4,5,6,19,9,7,15,8)]), pt.size = 1)

table(Idents(sc), sc$diet)

saveRDS(sc, "sc_ECs.rds")

###############
# 4. heart

Idents(heart) = heart$heart_celltype
DimPlot(object = heart, reduction = "umap", label = TRUE) + NoLegend()
# Endothelial subtypes Good Markers
FeaturePlot(heart, features = c("Vwf", "Vcam1", #  large vessels
                             "Fbln5", "Cytl1", # large artery
                             "Hey1", "Gkn3",# artery 
                             "Tgfb2", "Glul", # CapA
                             "Rgcc", "Mfsd2a", # Cap
                             "Car4", "Tfrc", # CapV
                             "Lcn2", "Slc38a5", # large vein
                             "Hba-a1", "Hbb-bs", 
                             "Plvap", "Plpp3", "Mki67",
                             "Ccl21a","Prox1","Lyve1",  # lymph
                             "Apln", "Col4a2", "Trp53i11"), # angiogen
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 5) #

# rename
heart <- RenameIdents(heart, #`4` = "EC-art", `9` = "EC-ven",`5` = "EC-venule",
                   `EC-cap1` = "EC-cap",`EC-cap2` = "EC-cap",`EC-IFN` = "EC-cap",
                   `EC-AP1` = "EC-cap",`6` = "EC-ang") 

## Relevel object@ident
heart@active.ident <- factor(x = heart@active.ident, levels = sort(levels(Idents(heart))))
heart$gen_celltype <- Idents(heart)
col.gen
levels(Idents(heart))
DimPlot(heart, reduction = "umap", label = TRUE, cols = as.vector(col.gen[c(6,7,1,2,8,9,10,5,3,4)]), pt.size = 1)

# number of cells per cluster per sample
table(Idents(heart), heart$diet)

# Find marker genes
Idents(heart) = heart$gen_celltype
heart.markers = FindAllMarkers(heart, only.pos = T)
write.table(heart.markers, "heart_gen_markers.txt", sep='\t')

top5.mark <- as.data.frame(heart.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "heart_top5_gen_markers.txt", sep='\t')
DotPlot(heart, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

# heart-specific Capillaries : Rgcc, Car4
FeaturePlot(heart, features = c("Vwf", "Vcam1", "Nr2f2",#  vein
                             "Hey1","Fbln5","Tgfb2", # artery 
                             "Rgcc","Car4","Cd300lg", # heart capillaries
                             "Flt4", "Ccl21a", "Prox1"), # lymph
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 3) #

# DE genes
heart$gen.diet <- paste(heart$gen_celltype, heart$diet, sep = "_")
Idents(object = heart) <- "gen.diet"
heart@active.ident <- factor(x = heart@active.ident, levels = sort(levels(Idents(object = heart))))
levels(Idents(object = heart))
x = levels(Idents(object = heart))[-3]

lev1 <- x[seq(2,(length(x)),2)] # condition 1=HFD
lev2 <- x[seq(1,(length(x)),2)] # condition 2=chow

# for all the clusters
for (i in seq(1,length(lev1),1)) { # 16 for not-filtered
  print(lev1[i])
  print(lev2[i])
  # TG 19m vs WT 19m
  response2 <- FindMarkers(object = heart, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = FALSE) #, min.pct = 0.25
  write.table(response2, file = paste("heart_",lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
}


saveRDS(heart, "heart_ECs.rds")

###############
# 5. lung

Idents(lung) = lung$lung_celltype
DimPlot(object = lung, reduction = "umap", label = TRUE) + NoLegend()
# Endothelial subtypes Good Markers
FeaturePlot(lung, features = c("Vwf", "Vcam1", #  large vessels
                                "Fbln5", "Cytl1", # large artery
                                "Hey1", "Gkn3",# artery 
                                "Tgfb2", "Glul", # CapA
                                "Rgcc", "Mfsd2a", # Cap
                                "Car4", "Tfrc", # CapV
                                "Lcn2", "Slc38a5", # large vein
                                "Hba-a1", "Hbb-bs", 
                                "Plvap", "Plpp3", "Mki67",
                                "Ccl21a","Prox1","Lyve1",  # lymph
                                "Apln", "Col4a2", "Trp53i11"), # angiogen
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 5) #


# rename
lung <- RenameIdents(lung, #`4` = "EC-art", `9` = "EC-ven",`5` = "EC-venule",
                      `EC-cap1` = "EC-cap",`EC-cap2` = "EC-cap",
                     `8` = "aEC",`3` = "aEC",
                     `4` = "pulmEC",`10` = "pulmEC",`14` = "EC-platelet",
                     `5` = "EC-Aqp5",`11` = "EC-Aqp5") 

## Relevel object@ident
lung@active.ident <- factor(x = lung@active.ident, levels = sort(levels(Idents(lung))))
lung$gen_celltype <- Idents(lung)
col.gen
levels(Idents(lung))
DimPlot(lung, reduction = "umap", label = TRUE, cols = as.vector(col.gen[c(12,7,8,1,2,5,6,9,11,3,4)]), pt.size = 1)

# number of cells per cluster per sample
table(Idents(lung), lung$diet)

# Find marker genes
Idents(lung) = lung$gen_celltype
lung.markers = FindAllMarkers(lung, only.pos = T)
write.table(lung.markers, "lung_gen_markers.txt", sep='\t')

top5.mark <- as.data.frame(lung.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "lung_top5_gen_markers.txt", sep='\t')
DotPlot(lung, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

# lung-specific Capillaries : Rgcc, Car4
FeaturePlot(lung, features = c("Vwf", "Vcam1", "Selp",#  vein
                                "Hey1","Atp13a3","Tgfb2", # artery 
                                "Rgcc","Car4","Kit", # lung capillaries
                                "Aqp5","Ednrb","Hpgd", # cap-aerocyte
                               "Slc2a3","Spp1","Postn", # bronchial vessel
                                "Flt4", "Ccl21a", "Prox1"), # lymph
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 3) #


# DE genes
lung$gen.diet <- paste(lung$gen_celltype, lung$diet, sep = "_")
Idents(object = lung) <- "gen.diet"
lung@active.ident <- factor(x = lung@active.ident, levels = sort(levels(Idents(object = lung))))
levels(Idents(object = lung))

lev1 <- levels(lung@active.ident)[seq(2,(length(levels(Idents(object = lung)))),2)] # condition 1=HFD
lev2 <- levels(lung@active.ident)[seq(1,(length(levels(Idents(object = lung)))),2)] # condition 2=chow

# for all the clusters
for (i in seq(1,length(lev1),1)) { # 16 for not-filtered
  print(lev1[i])
  print(lev2[i])
  # TG 19m vs WT 19m
  response2 <- FindMarkers(object = lung, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = FALSE) #, min.pct = 0.25
  write.table(response2, file = paste("lung_",lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
}


# subclusters
Idents(lung) = lung$lung_celltype
DimPlot(object = lung, reduction = "umap", label = TRUE) + NoLegend()
DimPlot(object = lung, reduction = "umap",group.by = "gen_celltype", label = TRUE) + NoLegend()

# rename
lung <- RenameIdents(lung, 
                     `8` = "aEC1",`3` = "aEC2",
                     `4` = "pulmEC-a",`10` = "pulmEC-b",`14` = "EC-platelet",
                     `5` = "EC-Aqp5a",`11` = "EC-Aqp5b") 

## Relevel object@ident
lung@active.ident <- factor(x = lung@active.ident, levels = sort(levels(Idents(lung))))
lung$lung_celltype <- Idents(lung)

levels(Idents(lung))

DimPlot(lung, reduction = "umap", label = TRUE, cols = as.vector(col.gen2[c(10,11,4,5,
                                                                            1,2,3,9,13,
                                                                            12,19,18,17,
                                                                            20,16,7,8)]), pt.size = 1)
# 7.7 x 7
table(Idents(lung), lung$diet)

saveRDS(lung, "lung_ECs.rds")

###############
# 6. brain

Idents(brain) = brain$brain_celltype
DimPlot(object = brain, reduction = "umap", label = TRUE) + NoLegend()
# Endothelial subtypes Good Markers
FeaturePlot(brain, features = c("Vwf", "Vcam1", #  large vessels
                                "Fbln5", "Cytl1", # large artery
                                "Hey1", "Gkn3",# artery 
                                "Tgfb2", "Glul", # CapA
                                "Rgcc", "Mfsd2a", # Cap
                                "Car4", "Tfrc", # CapV
                                "Lcn2", "Slc38a5", # large vein
                                "Hba-a1", "Hbb-bs", 
                                "Plvap", "Plpp3", "Mki67",
                                "Ccl21a","Lyve1",  # lymph
                                "Apln", "Col4a2", "Trp53i11"), # angiogen
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 5) #

# rename
brain <- RenameIdents(brain, #`4` = "EC-art", `9` = "EC-ven",`5` = "EC-venule",
                      `EC-cap1` = "EC-cap",`EC-IFN` = "EC-cap",
                      `EC-capA1` = "EC-cap",`EC-capA2` = "EC-cap",`6` = "EC-cap",
                      `EC-AP1` = "EC-cap",`8` = "EC-brain") 

## Relevel object@ident
brain@active.ident <- factor(x = brain@active.ident, levels = sort(levels(Idents(brain))))
brain$gen_celltype <- Idents(brain)
col.gen
levels(Idents(brain))
DimPlot(brain, reduction = "umap", label = TRUE, cols = as.vector(col.gen[c(1,8,2,9,6,3)]), pt.size = 1)

# number of cells per cluster per sample
table(Idents(brain), brain$diet)

# Find marker genes
Idents(brain) = brain$gen_celltype
brain.markers = FindAllMarkers(brain, only.pos = T)
write.table(brain.markers, "brain_gen_markers.txt", sep='\t')

top5.mark <- as.data.frame(brain.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "brain_top5_gen_markers.txt", sep='\t')
DotPlot(brain, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

# brain-specific Capillaries : "Rgcc", Glul,"Mfsd2a"
FeaturePlot(brain, features = c("Vwf", "Cfh", "Lcn2",#  vein
                                "Gkn3","Fbln5","Cytl1", # artery 
                                "Rgcc","Car4","Mfsd2a"), # brain capillaries
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 3) #

# DE genes
brain$gen.diet <- paste(brain$gen_celltype, brain$diet, sep = "_")
Idents(object = brain) <- "gen.diet"
brain@active.ident <- factor(x = brain@active.ident, levels = sort(levels(Idents(object = brain))))
levels(Idents(object = brain))

lev1 <- levels(brain@active.ident)[seq(2,(length(levels(Idents(object = brain)))),2)] # condition 1=HFD
lev2 <- levels(brain@active.ident)[seq(1,(length(levels(Idents(object = brain)))),2)] # condition 2=chow

# for all the clusters
for (i in seq(1,length(lev1),1)) { # 16 for not-filtered
  print(lev1[i])
  print(lev2[i])
  # TG 19m vs WT 19m
  response2 <- FindMarkers(object = brain, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = FALSE) #, min.pct = 0.25
  write.table(response2, file = paste("brain_",lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
}

# subclusters
Idents(object = brain) <- "brain_celltype"
DimPlot(object = brain, reduction = "umap", label = TRUE) + NoLegend()
DimPlot(object = brain, reduction = "umap",group.by = "gen_celltype", label = TRUE) + NoLegend()

# rename
brain <- RenameIdents(brain, `6` = "EC-cap2",`8` = "EC-brain") 

## Relevel object@ident
brain@active.ident <- factor(x = brain@active.ident, levels = sort(levels(Idents(brain))))
brain$brain_celltype <- Idents(brain)

levels(Idents(brain))
DimPlot(brain, reduction = "umap", label = TRUE, cols = as.vector(col.gen2[c(16,1,10,
                                                                             4,5,2,3,
                                                                              12,18,20,7)]), pt.size = 1)
table(Idents(brain), brain$diet)

saveRDS(brain, "brain_ECs.rds")

###############
# 7. kidney

Idents(kidney) = kidney$kidney_celltype
DimPlot(object = kidney, reduction = "umap", label = TRUE) + NoLegend()
# Endothelial subtypes Good Markers
FeaturePlot(kidney, features = c("Vwf", "Vcam1", #  large vessels
                                "Fbln5", "Cytl1", # large artery
                                "Hey1", "Gkn3",# artery 
                                "Tgfb2", "Glul", # CapA
                                "Rgcc", "Mfsd2a", # Cap
                                "Car4", "Tfrc", # CapV
                                "Lcn2", "Slc38a5", # large vein
                                "Hba-a1", "Hbb-bs", 
                                "Plvap", "Plpp3", "Mki67",
                                "Ccl21a","Prox1","Lyve1",  # lymph
                                "Apln", "Col4a2", "Trp53i11"), # angiogen
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 5) #

# rename
kidney <- RenameIdents(kidney, `2` = "EC-venule",
                      `EC-cap1` = "EC-cap",`EC-cap2` = "EC-cap",`4` = "EC-cap",
                      `EC-AP1` = "EC-cap",`EC-Aqp1` = "EC-Aqp1-arteriole",
                      `EC-enzym` = "EC-glomeruli",`5` = "mEC1",`7` = "mEC2") 

## Relevel object@ident
kidney@active.ident <- factor(x = kidney@active.ident, levels = sort(levels(Idents(kidney))))
kidney$gen_celltype <- Idents(kidney)
col.gen
levels(Idents(kidney))
DimPlot(kidney, reduction = "umap", label = TRUE, cols = as.vector(col.gen[c(9,7,12,1,2,8,6,10,5,3,11)]), pt.size = 1)

# number of cells per cluster per sample
table(Idents(kidney), kidney$diet)

# Find marker genes
Idents(kidney) = kidney$gen_celltype
kidney.markers = FindAllMarkers(kidney, only.pos = T)
write.table(kidney.markers, "kidney_gen_markers.txt", sep='\t')

top5.mark <- as.data.frame(kidney.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "kidney_top5_gen_markers.txt", sep='\t')
DotPlot(kidney, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

# kidney-specific Capillaries : Rgcc, "Igfbp5","Kdr"
FeaturePlot(kidney, features = c("Vwf", "Cfh", "Nr2f2","Igf1","Bgn",#  vein
                                "Hey1","Fbln5","Tgfb2", "Eln","Gja5",# artery 
                                "Rgcc","Igfbp5","Kdr", # kidney capillaries
                                "Aqp1", "Pi16", "Fxyd2",# kidney-specif
                                "Flt4", "Ccl21a", "Prox1"), # lymph
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 4) #

# DE genes
kidney$gen.diet <- paste(kidney$gen_celltype, kidney$diet, sep = "_")
Idents(object = kidney) <- "gen.diet"
kidney@active.ident <- factor(x = kidney@active.ident, levels = sort(levels(Idents(object = kidney))))
levels(Idents(object = kidney))

lev1 <- levels(kidney@active.ident)[seq(2,(length(levels(Idents(object = kidney)))),2)] # condition 1=HFD
lev2 <- levels(kidney@active.ident)[seq(1,(length(levels(Idents(object = kidney)))),2)] # condition 2=chow

# for all the clusters
for (i in seq(1,length(lev1),1)) { # 16 for not-filtered
  print(lev1[i])
  print(lev2[i])
  # TG 19m vs WT 19m
  response2 <- FindMarkers(object = kidney, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = FALSE) #, min.pct = 0.25
  write.table(response2, file = paste("kidney_",lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
}

# subclusters
Idents(kidney) = kidney$kidney_celltype
DimPlot(object = kidney, reduction = "umap", label = TRUE) + NoLegend()
DimPlot(object = kidney, reduction = "umap",group.by = "gen_celltype", label = TRUE) + NoLegend()

# rename
kidney <- RenameIdents(kidney, `2` = "EC-venule",`4` = "EC-cap3",
                       `EC-Aqp1` = "EC-Aqp1-arteriole",
                       `EC-enzym` = "EC-glomeruli",`5` = "mEC1",`7` = "mEC2") 

## Relevel object@ident
kidney@active.ident <- factor(x = kidney@active.ident, levels = sort(levels(Idents(kidney))))
kidney$kidney_celltype <- Idents(kidney)

levels(Idents(kidney))

DimPlot(kidney, reduction = "umap", label = TRUE, cols = as.vector(col.gen2[c(12,10,20,16,
                                                                            1,2,3,4,
                                                                            17,13,14,
                                                                            9,7,15)]), pt.size = 1)
# 7.9 x 7
table(Idents(kidney), kidney$diet)

saveRDS(kidney, "kidney_ECs.rds")


