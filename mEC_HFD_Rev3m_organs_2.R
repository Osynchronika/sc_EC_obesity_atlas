library(Seurat)
library(dplyr)
library(Matrix)
require(gplots)
require(ggplot2)
library(RColorBrewer)
library(cowplot)
require(scales)
library(DoubletFinder)
library(Rfast)
library(corrplot)

# Load data
setwd("~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_06_09_M_HFD_Rev_3m/analysis/NEW_SEQ/organs")
file.dir = "~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_06_09_M_HFD_Rev_3m/analysis/organs"
data.dir = "~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_06_09_M_HFD_Rev_3m/data"
data_dirs = list.files(data.dir)
data_dirs

data.dir.3m = "~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/organs-3rd-processing"
data.dir.4m = "~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_04_15_M_HFD_Rev_1m/analysis/combined/"
data.dir.6m = "~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_06_09_M_HFD_Rev_3m/analysis/NEW_SEQ/combined/"

########################################
# analyze organs separately
########################################
#########################
# 3. lymph
#########################
setwd("lymph")

lymph.3m = readRDS(paste0("~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis//combined/lymphatic_subset.rds"))
lymph.4m = readRDS(paste0(data.dir.4m,"lymphatic_subset.rds"))
lymph.6m = readRDS(paste0(data.dir.6m,"lymphatic_subset.rds"))

lymph.3m$time = "3m"

DimPlot(lymph.3m, reduction = "umap", label = TRUE)

# merge objects
lymph = merge(lymph.3m, list(lymph.4m, lymph.6m))
rm(lymph.3m, lymph.4m, lymph.6m)

lymph <- ScaleData(object = lymph, verbose = T)
lymph <- FindVariableFeatures(lymph, selection.method = "vst", nfeatures = 3000)
lymph <- RunPCA(object = lymph, npcs = 30, verbose = T)
ElbowPlot(object = lymph,  ndims = 30) 

lymph$diet_time = paste0(lymph$diet, lymph$time)

# UMAP and Clustering
lymph <- RunUMAP(object = lymph, reduction = "pca", dims = 1:21)
lymph <- FindNeighbors(object = lymph, reduction = "pca", dims = 1:21, nn.method = "rann")
lymph <- FindClusters(object = lymph, resolution = 0.2)


# visualization
DimPlot(object = lymph, reduction = "umap", label = TRUE)
DimPlot(object = lymph, reduction = "umap", group.by = "orig.ident")
DimPlot(object = lymph, reduction = "umap", split.by = "organ", label = T, ncol = 3) + NoLegend()
DimPlot(object = lymph, reduction = "umap", split.by = "diet_time", ncol = 3, pt.size = 1.1)
DimPlot(object = lymph, reduction = "umap", group.by = "organ", label = T,
        split.by = "diet_time", ncol = 3)


# Find marker genes
lymph.markers = FindAllMarkers(lymph, only.pos = T)
write.table(lymph.markers, "lymph_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(lymph.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
write.table(top5.mark, "lymph_top5_cluster_markers.txt", sep='\t')
DotPlot(lymph, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

# Endothelial subtypes Good Markers
FeaturePlot(lymph, features = c("Vwf", "Vcam1", #  large vessels
                                "Fbln5", "Cytl1", # large artery
                                "Hey1", "Gkn3",# artery 
                                "Tgfb2", "Glul", # CapA
                                "Rgcc",  # Cap
                                "Car4", "Tfrc", # CapV
                                "Lcn2",  # large vein
                                "Hba-a1", "Hbb-bs", 
                                "Plvap", "Plpp3",
                                "Flt4", "Ccl21a", "Prox1"),
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 4) #

FeaturePlot(lymph, features = c("Prox1", "Pdpn","Lyve1", "Flt4", "Ccl21a", # lymph
                                "Mrc1", "Ptprc"), 
            order = T, cols=c("grey", "red"), ncol=2, min.cutoff = "q10") #, min.cutoff = "q10"

#  tip cell-enriched, angiogenic
FeaturePlot(lymph, features = c("Apln", "Col4a2", "Trp53i11"),
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 3)

lymph <- RenameIdents(lymph, `2` = "2-lung",`6` = "6-lung",`7` = "7-immune?") 
## Relevel object@ident
lymph@active.ident <- factor(x = lymph@active.ident, levels = sort(levels(Idents(lymph))))
lymph$lymph_celltype2 <- Idents(lymph)
DimPlot(lymph, reduction = "umap", label = TRUE)

lymph$clust_organ = paste0(lymph$lymph_celltype2, "_",lymph$organ)

# number of cells per cluster per sample
table(Idents(lymph), lymph$diet_time)
table(lymph$clust_organ, lymph$diet_time)
table(Idents(lymph), lymph$organ)
table(lymph$diet_time, lymph$organ)

DimPlot(lymph, reduction = "umap", label = TRUE, split.by = "diet", ncol = 2) + NoLegend()
DimPlot(object = lymph, reduction = "umap", group.by = "lymph_celltype2", label = T,
        split.by = "diet_time", ncol = 3)+ NoLegend()

saveRDS(lymph,"lymph_3m_4m_6m.rds")

lymph = readRDS("lymph_3m_4m_6m.rds")

# DE genes clusters
setwd("deg-clust")
lymph$clust.diet_time <- paste(lymph$lymph_celltype2, lymph$diet_time, sep = "_")
Idents(object = lymph) <- "clust.diet_time"
lymph@active.ident <- factor(x = lymph@active.ident, levels = sort(levels(Idents(object = lymph))))
levels(Idents(object = lymph))

lev1 <- levels(lymph@active.ident)[seq(7,(length(levels(Idents(lymph)))),8)] # condition 1=HFD_4m
lev2 <- levels(lymph@active.ident)[seq(2,(length(levels(Idents(lymph)))),8)] # condition 2=chow_4m
lev3 <- levels(lymph@active.ident)[seq(4,(length(levels(Idents(lymph)))),8)] # condition 3=Rev_1m
lev4 <- levels(lymph@active.ident)[seq(8,(length(levels(Idents(lymph)))),8)] # condition 4=HFD_6m
lev5 <- levels(lymph@active.ident)[seq(3,(length(levels(Idents(lymph)))),8)] # condition 5=chow_6m
lev6 <- levels(lymph@active.ident)[seq(5,(length(levels(Idents(lymph)))),8)] # condition 6=Rev_3m
lev7 <- levels(lymph@active.ident)[seq(6,(length(levels(Idents(lymph)))),8)] # condition 7=HFD_3m
lev8 <- levels(lymph@active.ident)[seq(1,(length(levels(Idents(lymph)))),8)] # condition 8=chow_3m

# for all the clusters
for (i in seq(1,length(lev1),1)) {
  # 4 months
  # Western
  print(lev1[i])
  print(lev2[i])
  response2 <- FindMarkers(object = lymph, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # Reverse
  print(lev3[i])
  print(lev2[i])
  response2 <- FindMarkers(object = lymph, ident.1 = lev3[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev3[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # 6 months
  # Western
  print(lev4[i])
  print(lev5[i])
  response2 <- FindMarkers(object = lymph, ident.1 = lev4[i], 
                           ident.2 = lev5[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev4[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # Reverse
  print(lev6[i])
  print(lev5[i])
  response2 <- FindMarkers(object = lymph, ident.1 = lev6[i], 
                           ident.2 = lev5[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev6[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # 3 months
  # Western
  print(lev7[i])
  print(lev8[i])
  response2 <- FindMarkers(object = brain, ident.1 = lev7[i], 
                           ident.2 = lev8[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev7[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
}



########################################
# 2. Fat
############


sc.3m = readRDS(paste0(data.dir.3m,"/sc/sc_ECs_singlet.rds"))
sc.4m = readRDS(paste0(data.dir.4m,"sc_subset.rds"))
sc.6m = readRDS(paste0(data.dir.6m,"sc_subset.rds"))
sc.3m$time = "3m"

##############
# A. sc separate
#############
setwd("sc")

# merge objects
sc = merge(sc.3m, list(sc.4m, sc.6m))
rm(sc.3m, sc.4m, sc.6m)
sc$diet_time = paste0(sc$diet, sc$time)

sc <- ScaleData(object = sc, verbose = T)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)
sc <- RunPCA(object = sc, npcs = 30, verbose = T)
ElbowPlot(object = sc,  ndims = 30) 

# UMAP and Clustering
sc <- RunUMAP(object = sc, reduction = "pca", dims = 1:30)

sc = readRDS("sc_all_3m_4m_6m.rds")
#nn.method from "annoy" to "rann"
sc <- FindNeighbors(object = sc, reduction = "pca", dims = 1:30, nn.method = "rann") 
sc <- FindClusters(object = sc, resolution = 0.3)

# Visualization
#Idents(sc) = sc$RNA_snn_res.0.3
DimPlot(object = sc, reduction = "umap", label = TRUE)
DimPlot(object = sc, reduction = "umap", group.by = "diet_time")
#DimPlot(object = sc, reduction = "umap", group.by = "organ", cols = c("purple", "orange"))
DimPlot(object = sc, reduction = "umap", group.by = "sc_celltype", label = T)

# select cluster by hand
plot <- DimPlot(sc, reduction = "umap")
select.cells <- CellSelector(plot = plot) # 
Idents(sc, cells = select.cells) <- "16"

## Relevel object@ident
sc@active.ident <- factor(x = sc@active.ident, levels = sort(levels(Idents(sc))))
sc$seurat_clusters <- Idents(sc)
DimPlot(sc, reduction = "umap", label = TRUE)

# Find marker genes
Idents(sc) = sc$seurat_clusters
sc.markers = FindAllMarkers(sc, only.pos = T)
write.table(sc.markers, "sc_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(sc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "sc_top5_cluster_markers.txt", sep='\t')
DotPlot(sc, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

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
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 4, raster = T) #

# Pericyte amd SMC markers
FeaturePlot(sc, features = c("Flt1", "Pecam1","Cdh5",# EC
                             "Pdgfrb", "Cspg4","Anpep",  #  Pericyte
                             "Myh11", "Acta2", "Des", # SMC
                             "Prox1", "Pdpn","Lyve1", # lymph
                             "Dcn", "Col1a1", "Pdgfra"), #FB
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10", raster = T) #, min.cutoff = "q10"


#  tip cell-enriched, angiogenic
FeaturePlot(sc, features = c("Apln", "Col4a2", "Trp53i11", 
                             "Mki67", "Cenpf"), # prolif
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 3)



# rename clusters
sc <- RenameIdents(object = sc, `13` = "Immune1",`12` = "Immune2",
                   `8` = "Mesench",`9` = "EC-lymph",
                   `11` = "Prolif",`10` = "EC-Hb",
                   `5` = "EC-ang",`1` = "EC-art") 
## Relevel object@ident
sc@active.ident <- factor(x = sc@active.ident, levels = sort(levels(Idents(sc))))
sc$sc_celltype2 <- Idents(sc)
DimPlot(sc, reduction = "umap", label = TRUE)
# number of cells per cluster per sample
table(Idents(sc), sc$diet_time)

saveRDS(sc, "sc_all_3m_4m_6m.rds")

####
sc = readRDS("sc_all_3m_4m_6m.rds")
# remove immune cells and mesenchymal cells?
sc = subset(sc, idents = c("Immune1","Immune2","Mesench"), invert = T)
DimPlot(sc, reduction = "umap")

# select cluster by hand
plot <- DimPlot(sc, reduction = "umap")
select.cells <- CellSelector(plot = plot) # 
Idents(sc, cells = select.cells) <- "outlier"

sc = subset(sc, idents = c("outlier"), invert = T)

# re-cluster
sc <- ScaleData(object = sc, verbose = T)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)
sc <- RunPCA(object = sc, npcs = 30, verbose = T)
ElbowPlot(object = sc,  ndims = 30) 

# UMAP and Clustering
sc <- RunUMAP(object = sc, reduction = "pca", dims = 1:28)
sc <- FindNeighbors(object = sc, reduction = "pca", dims = 1:28, nn.method = "rann")
sc <- FindClusters(object = sc, resolution = 0.3)

# Visualization
DimPlot(object = sc, reduction = "umap", label = TRUE)
DimPlot(object = sc, reduction = "umap", group.by = "diet_time")
#DimPlot(object = sc, reduction = "umap", group.by = "organ", cols = c("purple", "orange"))
DimPlot(object = sc, reduction = "umap", group.by = "sc_celltype", label = T)
DimPlot(object = sc, reduction = "umap", group.by = "sc_celltype2", label = T)

# select cluster by hand
plot <- DimPlot(sc, reduction = "umap")
select.cells <- CellSelector(plot = plot) # 
Idents(sc, cells = select.cells) <- "lymph2"
Idents(sc, cells = select.cells) <- "EC-ven"

# rename clusters
sc <- RenameIdents(object = sc, `8` = "EC-lymph",
                   `11` = "Prolif",`10` = "EC-Hb",
                   `6` = "EC-ang",`9` = "EC-art", `1` = "EC-arteriole",
               `7` = "EC-fenestr",`2` = "EC-venule",
             `0` = "EC-cap1", `4` = "EC-cap2",`3` = "EC-cap3",`5` = "EC-cap4") 

## Relevel object@ident
sc@active.ident <- factor(x = sc@active.ident, levels = sort(levels(Idents(sc))))
sc$sc_celltype2 <- Idents(sc)

DimPlot(sc, reduction = "umap", label = TRUE)

# number of cells per cluster per sample
table(Idents(sc), sc$diet_time) 

DimPlot(sc, reduction = "umap", label = TRUE, split.by = "diet_time", ncol = 3) + NoLegend()

saveRDS(sc, "sc_subset_3m_4m_6m.rds")

###########################
# access Cell Cycle

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase - ONLY FOR HUMAN
# Load from file - for human and mouse
cell.cycle <- read.table("~/Desktop/HI-MAG_OLGA/R scripts/cell_cycle_Seurat.txt", header = T, sep = '\t')
s.genes <- cell.cycle[cell.cycle$type=="s.genes","mouse_corr"]
g2m.genes <- cell.cycle[cell.cycle$type=="g2m.genes","mouse_corr"]
sc <- CellCycleScoring(sc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
sc$sc_celltype2.diet <- paste(sc$sc_celltype2, sc$diet, sep = "_")
Idents(sc) <- "sc_celltype2.diet"

# view cell cycle scores and phase assignments
sc$sc_celltype2.diet_time <- paste(sc$sc_celltype2, sc$diet_time, sep = "_")
Idents(sc) <- "sc_celltype2.diet_time"

# get the fraction of the phases per cluster
prop.table(table(sc@meta.data[,c("Phase","sc_celltype2.diet_time")]), margin = 2)

prop.table(table(sc@meta.data[,c("Phase","sc_celltype2")]), margin = 2)

prop.table(table(sc@meta.data[,c("Phase","diet_time")]), margin = 2)

FeaturePlot(sc, features = c("Mki67", "Pcna"), 
            cols = c("grey", "red"), reduction = "umap", split.by = "diet",label = F, order = T)

FeaturePlot(sc, features = c("S.Score", "G2M.Score"), 
            cols = c("grey", "red"), reduction = "umap", split.by = "diet",label = F, order = T, min.cutoff = "q5")

FeaturePlot(sc, features = c("S.Score", "G2M.Score"), 
            reduction = "umap", label = F, order = T, blend = T)

table(sc@meta.data[,c("diet_time","sc_celltype2")])

##################

sc = readRDS("sc_subset_3m_4m_6m.rds")

# find marker genes
Idents(sc) <- "sc_celltype2"
sc.markers = FindAllMarkers(sc, only.pos = T)
write.table(sc.markers, "sc_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(sc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
write.table(top5.mark, "sc_top5_cluster_markers.txt", sep='\t')
DotPlot(sc, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()


# DE genes
sc$clust.diet_time <- paste(sc$sc_celltype2, sc$diet_time, sep = "_")
Idents(object = sc) <- "clust.diet_time"
sc@active.ident <- factor(x = sc@active.ident, levels = sort(levels(Idents(object = sc))))
levels(Idents(object = sc))

lev1 <- levels(sc@active.ident)[seq(7,(length(levels(Idents(sc)))),8)] # condition 1=HFD_4m
lev2 <- levels(sc@active.ident)[seq(2,(length(levels(Idents(sc)))),8)] # condition 2=chow_4m
lev3 <- levels(sc@active.ident)[seq(4,(length(levels(Idents(sc)))),8)] # condition 3=Rev_1m
lev4 <- levels(sc@active.ident)[seq(8,(length(levels(Idents(sc)))),8)] # condition 4=HFD_6m
lev5 <- levels(sc@active.ident)[seq(3,(length(levels(Idents(sc)))),8)] # condition 5=chow_6m
lev6 <- levels(sc@active.ident)[seq(5,(length(levels(Idents(sc)))),8)] # condition 6=Rev_3m
lev7 <- levels(sc@active.ident)[seq(6,(length(levels(Idents(sc)))),8)] # condition 7=HFD_3m
lev8 <- levels(sc@active.ident)[seq(1,(length(levels(Idents(sc)))),8)] # condition 8=chow_3m

# for all the clusters
for (i in seq(1,length(lev1),1)) {
  # 4 months
  # Western
  print(lev1[i])
  print(lev2[i])
  response2 <- FindMarkers(object = sc, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # Reverse
  print(lev3[i])
  print(lev2[i])
  response2 <- FindMarkers(object = sc, ident.1 = lev3[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev3[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # 6 months
  # Western
  print(lev4[i])
  print(lev5[i])
  response2 <- FindMarkers(object = sc, ident.1 = lev4[i], 
                           ident.2 = lev5[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev4[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # Reverse
  print(lev6[i])
  print(lev5[i])
  response2 <- FindMarkers(object = sc, ident.1 = lev6[i], 
                           ident.2 = lev5[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev6[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # 3 months
  # Western
  print(lev7[i])
  print(lev8[i])
  response2 <- FindMarkers(object = sc, ident.1 = lev7[i], 
                           ident.2 = lev8[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev7[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
}


# DEG for capillaries combined
# 4m
response2 <- FindMarkers(object = sc, 
                         ident.1 = c("EC-cap1_Western4m","EC-cap2_Western4m",
                                     "EC-cap3_Western4m","EC-cap4_Western4m"), 
                         ident.2 = c("EC-cap1_chow4m","EC-cap2_chow4m",
                                     "EC-cap3_chow4m","EC-cap4_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "sc_Cap_WD_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = sc, 
                         ident.1 = c("EC-cap1_Rev1m","EC-cap2_Rev1m",
                                     "EC-cap3_Rev1m","EC-cap4_Rev1m"), 
                         ident.2 = c("EC-cap1_chow4m","EC-cap2_chow4m",
                                     "EC-cap3_chow4m","EC-cap4_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "sc_Cap_Rev1m_vs_chow_4m.txt", sep = "\t")

# 6m
response2 <- FindMarkers(object = sc, 
                         ident.1 = c("EC-cap1_Western6m","EC-cap2_Western6m",
                                     "EC-cap3_Western6m","EC-cap4_Western6m"), 
                         ident.2 = c("EC-cap1_chow6m","EC-cap2_chow6m",
                                     "EC-cap3_chow6m","EC-cap4_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "sc_Cap_WD_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = sc, 
                         ident.1 = c("EC-cap1_Rev3m","EC-cap2_Rev3m",
                                     "EC-cap3_Rev3m","EC-cap4_Rev3m"), 
                         ident.2 = c("EC-cap1_chow6m","EC-cap2_chow6m",
                                     "EC-cap3_chow6m","EC-cap4_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "sc_Cap_Rev3m_vs_chow_6m.txt", sep = "\t")

# DEG for vein and artery combined
# 4m
response2 <- FindMarkers(object = sc, 
                         ident.1 = c("EC-art_Western4m","EC-arteriole_Western4m"), 
                         ident.2 = c("EC-art_chow4m","EC-arteriole_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "sc_Art_WD_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = sc, 
                         ident.1 = c("EC-ven_Western4m","EC-venule_Western4m"), 
                         ident.2 = c("EC-ven_chow4m","EC-venule_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "sc_Ven_WD_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = sc, 
                         ident.1 = c("EC-art_Rev1m","EC-arteriole_Rev1m"), 
                         ident.2 = c("EC-art_chow4m","EC-arteriole_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "sc_Art_Rev1m_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = sc, 
                         ident.1 = c("EC-ven_Rev1m","EC-venule_Rev1m"), 
                         ident.2 = c("EC-ven_chow4m","EC-venule_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "sc_Ven_Rev1m_vs_chow_4m.txt", sep = "\t")

# 6m
response2 <- FindMarkers(object = sc, 
                         ident.1 = c("EC-art_Western6m","EC-arteriole_Western6m"), 
                         ident.2 = c("EC-art_chow6m","EC-arteriole_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "sc_Art_WD_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = sc, 
                         ident.1 = c("EC-ven_Western6m","EC-venule_Western6m"), 
                         ident.2 = c("EC-ven_chow6m","EC-venule_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "sc_Ven_WD_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = sc, 
                         ident.1 = c("EC-art_Rev3m","EC-arteriole_Rev3m"), 
                         ident.2 = c("EC-art_chow6m","EC-arteriole_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "sc_Art_Rev3m_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = sc, 
                         ident.1 = c("EC-ven_Rev3m","EC-venule_Rev3m"), 
                         ident.2 = c("EC-ven_chow6m","EC-venule_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "sc_Ven_Rev3m_vs_chow_6m.txt", sep = "\t")


#### HERE
####
# combine Cap/ Art/ Ven in 1 table
setwd("deg-gen")
# 1. Capillaries
tabs = list.files()[grep("sc_Cap_", list.files())]
deg.capillaries = list()
deg.i = NULL
for (i in seq(1,4,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  # change the base
  deg.i$avg_logFC = deg.i$avg_log2FC/log(exp(1),2)
  print(head(deg.i))
  deg.capillaries[[i]] = deg.i
}

deg.capillaries[[5]] = read.table("/Users/olga.bondareva/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new/deg_gen/sc_EC-cap_Western_vs_chow.txt", sep = '\t', header = T)

head(deg.capillaries[[5]])

# combine a table of logFC
tab.deg.cap = merge(deg.capillaries[[5]][,c(2,5)], deg.capillaries[[3]][,c(6,5)], by="row.names", all=TRUE)
colnames(tab.deg.cap) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[4]][,c(6,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[1]][,c(6,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[2]][,c(6,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.cap)

rownames(tab.deg.cap) = tab.deg.cap$gene
tab.deg.cap = tab.deg.cap[,-1]
nrow(tab.deg.cap) # 7674 genes
write.table(tab.deg.cap,"table_deg_sc_Cap_all_TP.txt", sep='\t')



# plot correlations

pdf(paste0("corr-cap-sc.pdf"),width=10,height=10,paper='special')
p1 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,WD_4m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_4m") + theme_linedraw() + geom_smooth(method = "lm")
p2 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,WD_6m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_6m") + theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,Rev_1m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_1m") + theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,Rev_3m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_31m") + theme_linedraw() + geom_smooth(method = "lm")

print((p1 | p2 ) /
        ( p3 |p4)) # 10 x 10
dev.off()

# 2. Artery
tabs = list.files()[grep("sc_Art_", list.files())]
deg.artery = list()
deg.i = NULL
for (i in seq(1,4,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  # change the base
  deg.i$avg_logFC = deg.i$avg_log2FC/log(exp(1),2)
  print(head(deg.i))
  deg.artery[[i]] = deg.i
}

deg.artery[[5]] = read.table("/Users/olga.bondareva/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new/deg_gen/sc_EC-art_Western_vs_chow.txt", sep = '\t', header = T)

head(deg.artery[[1]])

# combine a table of logFC
tab.deg.art = merge(deg.artery[[5]][,c(2,5)], deg.artery[[3]][,c(6,5)], by="row.names", all=TRUE)
colnames(tab.deg.art) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.art = merge(tab.deg.art, deg.artery[[4]][,c(6,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.art = merge(tab.deg.art, deg.artery[[1]][,c(6,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.art = merge(tab.deg.art, deg.artery[[2]][,c(6,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.art)

rownames(tab.deg.art) = tab.deg.art$gene
tab.deg.art = tab.deg.art[,-1]
nrow(tab.deg.art) # 8171 genes
write.table(tab.deg.art,"table_deg_sc_art_all_TP.txt", sep='\t')

# 3. Vein
tabs = list.files()[grep("sc_Ven_", list.files())]
deg.vein = list()
deg.i = NULL
for (i in seq(1,4,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  # change the base
  deg.i$avg_logFC = deg.i$avg_log2FC/log(exp(1),2)
  print(head(deg.i))
  deg.vein[[i]] = deg.i
}

deg.vein[[5]] = read.table("/Users/olga.bondareva/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new/deg_gen/sc_EC-ven_Western_vs_chow.txt", sep = '\t', header = T)

head(deg.vein[[3]])

# combine a table of logFC
tab.deg.ven = merge(deg.vein[[4]][,c(2,5)], deg.vein[[3]][,c(6,5)], by="row.names", all=TRUE)
colnames(tab.deg.ven) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.ven = merge(tab.deg.ven, deg.vein[[4]][,c(6,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.ven = merge(tab.deg.ven, deg.vein[[1]][,c(6,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.ven = merge(tab.deg.ven, deg.vein[[2]][,c(6,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.ven)

rownames(tab.deg.ven) = tab.deg.ven$gene
tab.deg.ven = tab.deg.ven[,-1]
nrow(tab.deg.ven) # 7752 genes
write.table(tab.deg.ven,"table_deg_sc_ven_all_TP.txt", sep='\t')

##########
# Find marker genes
Idents(sc) = sc$sc_celltype2
sc.markers = FindAllMarkers(sc, only.pos = T)
write.table(sc.markers, "sc_celltype2_markers.txt", sep='\t')

top5.mark <- as.data.frame(sc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "sc_top5_celltype2_markers.txt", sep='\t')
DotPlot(sc, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

saveRDS(sc, "sc_subset_3m_4m_6m.rds")

###############
# B. vis separate
##############
vis.3m = readRDS(paste0(data.dir.3m,"/vis/vis_ECs_singlet.rds"))
vis.4m = readRDS(paste0(data.dir.4m,"vis_subset.rds"))
vis.6m = readRDS(paste0(data.dir.6m,"vis_subset.rds"))
vis.3m$time = "3m"

setwd("../vis")

# merge objects
vis = merge(vis.3m, list(vis.4m,vis.6m))
rm(vis.3m, vis.4m, vis.6m)
vis$diet_time = paste0(vis$diet, vis$time)

vis <- ScaleData(object = vis, verbose = T)
vis <- FindVariableFeatures(vis, selection.method = "vst", nfeatures = 3000)
vis <- RunPCA(object = vis, npcs = 30, verbose = T)
ElbowPlot(object = vis,  ndims = 30) 

# UMAP and Clustering
vis <- RunUMAP(object = vis, reduction = "pca", dims = 1:25)
vis <- FindNeighbors(object = vis, reduction = "pca", dims = 1:25, nn.method = "rann")
vis <- FindClusters(object = vis, resolution = 0.4)

# Visualization
DimPlot(object = vis, reduction = "umap", label = TRUE)
DimPlot(object = vis, reduction = "umap", group.by = "diet_time")
DimPlot(object = vis, reduction = "umap", group.by = "vis_celltype", label = T)

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
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 4, raster = T) #

# Pericyte amd SMC markers
FeaturePlot(vis, features = c("Flt1", "Pecam1","Cdh5",# EC
                              "Pdgfrb", "Cspg4","Anpep",  #  Pericyte
                              "Myh11", "Acta2", "Des",   # SMC
                              "Prox1", "Pdpn","Lyve1",
                              "Ptprc", "Mrc1", "Cd52",
                              "Hba-a1", "Hbb-bs", "Hba-a2",
                              "Dcn", "Col1a1", "Pdgfra"), #FB
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10", raster = T) #, min.cutoff = "q10"

#  tip cell-enriched, angiogenic
FeaturePlot(vis, features = c("Apln", "Col4a2", "Trp53i11", 
                             "Mki67", "Cenpf"), # prolif
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 3)

# select cluster by hand
plot <- DimPlot(vis, reduction = "umap")
select.cells <- CellSelector(plot = plot) # 
Idents(vis, cells = select.cells) <- "Prolif"


# Find marker genes
#Idents(vis) = vis$seurat_clusters
vis.markers = FindAllMarkers(vis, only.pos = T)
write.table(vis.markers, "vis_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(vis.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "vis_top5_cluster_markers.txt", sep='\t')
DotPlot(vis, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()


vis <- RenameIdents(vis, `2` = "EC-lymph",`12` = "Immune",
                    `9` = "Mesench",`11` = "Mesench",`13` = "EC-vis",
                    `5` = "EC-ang",`4` = "EC-art",`10` = "EC-Hb",
                    `7` = "EC-ven",`6` = "EC-venule", `0` = "EC-cap1",`1` = "EC-cap2",
                    `3` = "EC-cap3", `8` = "EC-cap4") 
## Relevel object@ident
vis@active.ident <- factor(x = vis@active.ident, levels = sort(levels(Idents(vis))))
vis$vis_celltype2 <- Idents(vis)
DimPlot(vis, reduction = "umap", label = TRUE)

# number of cells per cluster per sample
table(vis$vis_celltype2, vis$diet_time)

saveRDS(vis, "vis_all_3m_4m_6m.rds")

vis = readRDS("vis_all_3m_4m_6m.rds")
###
# remove immune and mesenchymal cells
vis = subset(vis, idents = c("Immune", "Mesench"), invert = T)
vis <- ScaleData(object = vis, verbose = T)
vis <- FindVariableFeatures(vis, selection.method = "vst", nfeatures = 3000)
vis <- RunPCA(object = vis, npcs = 30, verbose = T)
ElbowPlot(object = vis,  ndims = 30) 

# UMAP and Clustering
vis <- RunUMAP(object = vis, reduction = "pca", dims = 1:25)
vis <- FindNeighbors(object = vis, reduction = "pca", dims = 1:25, nn.method = "rann")
vis <- FindClusters(object = vis, resolution = 0.4)

# Visualization
DimPlot(object = vis, reduction = "umap", label = TRUE, group.by = "vis_celltype2")
DimPlot(object = vis, reduction = "umap", label = TRUE)

vis <- RenameIdents(vis, `2` = "EC-lymph1",`11` = "EC-lymph2",
                    `6` = "EC-ang",`9` = "EC-art",`5` = "EC-arteriole",
                    `10` = "EC-Hb",`12` = "Prolif",
                    `8` = "EC-ven",`7` = "EC-venule", `0` = "EC-cap1",`1` = "EC-cap2",
                    `3` = "EC-cap3", `4` = "EC-cap4") 
## Relevel object@ident
vis@active.ident <- factor(x = vis@active.ident, levels = sort(levels(Idents(vis))))
vis$vis_celltype3 <- Idents(vis)
DimPlot(vis, reduction = "umap", label = TRUE)

# number of cells per cluster per sample
table(vis$vis_celltype3, vis$diet_time)

DimPlot(vis, reduction = "umap", label = TRUE, split.by = "diet_time", 
        ncol = 3, group.by = "vis_celltype3", pt.size = 1.3) + NoLegend()

###########################
# access Cell Cycle

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase - ONLY FOR HUMAN
# Load from file - for human and mouse
cell.cycle <- read.table("~/Desktop/HI-MAG_OLGA/R scripts/cell_cycle_Seurat.txt", header = T, sep = '\t')
s.genes <- cell.cycle[cell.cycle$type=="s.genes","mouse_corr"]
g2m.genes <- cell.cycle[cell.cycle$type=="g2m.genes","mouse_corr"]
vis <- CellCycleScoring(vis, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
vis$vis_celltype3.diet <- paste(vis$vis_celltype3, vis$diet, sep = "_")
Idents(vis) <- "vis_celltype3.diet"

# view cell cycle scores and phase assignments
vis$vis_celltype3.diet_time <- paste(vis$vis_celltype3, vis$diet_time, sep = "_")
Idents(vis) <- "vis_celltype3.diet_time"

# get the fraction of the phases per cluster
prop.table(table(vis@meta.data[,c("Phase","vis_celltype3.diet_time")]), margin = 2)

prop.table(table(vis@meta.data[,c("Phase","vis_celltype3")]), margin = 2)

prop.table(table(vis@meta.data[,c("Phase","diet_time")]), margin = 2)

saveRDS(vis, "vis_subset_3m_4m_6m.rds")

vis = readRDS("vis_subset_3m_4m_6m.rds")

#########
#  tip cell-enriched, angiogenic
FeaturePlot(vis, features = c("Apln", "Col4a2", "Trp53i11", 
                              "Mki67", "Cenpf", # prolif
                              "S100a4","Col1a1", "Col3a1","Acta2","Vim","Tgfb1","Snai1","Snai2","Twist1"), #EndoMT
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 3, raster = T)

FeaturePlot(vis, features = c("Col1a1","Acta2", "Tagln","Apln"), 
                              order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 2, raster = T)
Idents(vis) <- "vis_celltype3"
VlnPlot(vis, features = c("Apln", "Col4a2", "Trp53i11", 
                          #"Mki67", "Cenpf", # prolif
                          "Cdh5","Pecam1","Tek","Vwf", # EC
                          "S100a4","Col1a1", "Col3a1","Acta2","Vim","Tgfb1","Snai1","Snai2","Twist1"), #EndoMT
        idents = c("EC-ang"), split.by = "diet_time", stack = T)

DotPlot(vis, features = c("Apln", "Col4a2", "Trp53i11", 
                          #"Mki67", "Cenpf", # prolif
                          "Cdh5","Pecam1","Tek","Vwf", # EC
                          "S100a4","Col1a1", "Col3a1","Acta2","Vim","Tgfb1","Snai1","Snai2","Twist1"), #EndoMT
        idents = c("EC-ang"), split.by = "diet_time",
        cols = c("red","red","red","red","red","red","red","red"))+RotatedAxis()

ec.migration = read.table("/Users/olga.bondareva/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/GO_0043542_ec_migration.txt", sep = "\t")

ec.migration = unique(ec.migration$V2)
ec.migration = ec.migration[order(ec.migration)]
ec.migration = ec.migration[-c(9,67,52)] # no genes

vis = AddModuleScore(vis, features = list(c("Apln", "Col4a2", "Trp53i11",ec.migration)), name = "migration", ctrl = 5)
vis = AddModuleScore(vis, features = list(c("S100a4","Col1a1", "Col3a1","Acta2","Vim","Tgfb1","Snai1","Snai2","Twist1")), 
                     name = "EndoMT")

Idents(vis) <- "vis_celltype3"
VlnPlot(vis, features = c("migration1", "EndoMT1"),idents = c("EC-ang"),split.by = "diet_time")
DotPlot(vis, features = c("migration1", "EndoMT1"),
        idents = c("EC-ang"), split.by = "diet_time",
        cols = c("red","red","red","red","red","red","red","red"))+RotatedAxis()
FeaturePlot(vis, features = "Slc2a1", order = T, cols = c("grey","red"))

# find marker genes
vis.markers = FindAllMarkers(vis, only.pos = T)
write.table(vis.markers, "vis_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(vis.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
write.table(top5.mark, "vis_top5_cluster_markers.txt", sep='\t')
DotPlot(vis, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

##########

# DE genes
vis$clust.diet_time <- paste(vis$vis_celltype3, vis$diet_time, sep = "_")
Idents(object = vis) <- "clust.diet_time"
vis@active.ident <- factor(x = vis@active.ident, levels = sort(levels(Idents(object = vis))))
levels(Idents(object = vis))

lev1 <- levels(vis@active.ident)[seq(7,(length(levels(Idents(vis)))),8)] # condition 1=HFD_4m
lev2 <- levels(vis@active.ident)[seq(2,(length(levels(Idents(vis)))),8)] # condition 2=chow_4m
lev3 <- levels(vis@active.ident)[seq(4,(length(levels(Idents(vis)))),8)] # condition 3=Rev_1m
lev4 <- levels(vis@active.ident)[seq(8,(length(levels(Idents(vis)))),8)] # condition 4=HFD_6m
lev5 <- levels(vis@active.ident)[seq(3,(length(levels(Idents(vis)))),8)] # condition 5=chow_6m
lev6 <- levels(vis@active.ident)[seq(5,(length(levels(Idents(vis)))),8)] # condition 6=Rev_3m
lev7 <- levels(vis@active.ident)[seq(6,(length(levels(Idents(vis)))),8)] # condition 7=HFD_3m
lev8 <- levels(vis@active.ident)[seq(1,(length(levels(Idents(vis)))),8)] # condition 8=chow_3m

# for all the clusters
for (i in seq(1,length(lev1),1)) {
  # 4 months
  # Western
  print(lev1[i])
  print(lev2[i])
  response2 <- FindMarkers(object = vis, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # Reverse
  print(lev3[i])
  print(lev2[i])
  response2 <- FindMarkers(object = vis, ident.1 = lev3[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev3[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # 6 months
  # Western
  print(lev4[i])
  print(lev5[i])
  response2 <- FindMarkers(object = vis, ident.1 = lev4[i], 
                           ident.2 = lev5[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev4[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # Reverse
  print(lev6[i])
  print(lev5[i])
  response2 <- FindMarkers(object = vis, ident.1 = lev6[i], 
                           ident.2 = lev5[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev6[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # 3 months
  # Western
  print(lev7[i])
  print(lev8[i])
  response2 <- FindMarkers(object = vis, ident.1 = lev7[i], 
                           ident.2 = lev8[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev7[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
}


# DEG for capillaries combined
# 4m
response2 <- FindMarkers(object = vis, 
                         ident.1 = c("EC-cap1_Western4m","EC-cap2_Western4m",
                                     "EC-cap3_Western4m","EC-cap4_Western4m"), 
                         ident.2 = c("EC-cap1_chow4m","EC-cap2_chow4m",
                                     "EC-cap3_chow4m","EC-cap4_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "vis_Cap_WD_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = vis, 
                         ident.1 = c("EC-cap1_Rev1m","EC-cap2_Rev1m",
                                     "EC-cap3_Rev1m","EC-cap4_Rev1m"), 
                         ident.2 = c("EC-cap1_chow4m","EC-cap2_chow4m",
                                     "EC-cap3_chow4m","EC-cap4_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "vis_Cap_Rev1m_vs_chow_4m.txt", sep = "\t")

# 6m
response2 <- FindMarkers(object = vis, 
                         ident.1 = c("EC-cap1_Western6m","EC-cap2_Western6m",
                                     "EC-cap3_Western6m","EC-cap4_Western6m"), 
                         ident.2 = c("EC-cap1_chow6m","EC-cap2_chow6m",
                                     "EC-cap3_chow6m","EC-cap4_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "vis_Cap_WB_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = vis, 
                         ident.1 = c("EC-cap1_Rev3m","EC-cap2_Rev3m",
                                     "EC-cap3_Rev3m","EC-cap4_Rev3m"), 
                         ident.2 = c("EC-cap1_chow6m","EC-cap2_chow6m",
                                     "EC-cap3_chow6m","EC-cap4_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "vis_Cap_Rev3m_vs_chow_6m.txt", sep = "\t")

# DEG for vein and artery combined
# 4m
response2 <- FindMarkers(object = vis, 
                         ident.1 = c("EC-art_Western4m","EC-arteriole_Western4m"), 
                         ident.2 = c("EC-art_chow4m","EC-arteriole_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "vis_Art_WB_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = vis, 
                         ident.1 = c("EC-ven_Western4m","EC-venule_Western4m"), 
                         ident.2 = c("EC-ven_chow4m","EC-venule_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "vis_Ven_WB_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = vis, 
                         ident.1 = c("EC-art_Rev1m","EC-arteriole_Rev1m"), 
                         ident.2 = c("EC-art_chow4m","EC-arteriole_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "vis_Art_Rev1m_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = vis, 
                         ident.1 = c("EC-ven_Rev1m","EC-venule_Rev1m"), 
                         ident.2 = c("EC-ven_chow4m","EC-venule_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "vis_Ven_Rev1m_vs_chow_4m.txt", sep = "\t")

# 6m
response2 <- FindMarkers(object = vis, 
                         ident.1 = c("EC-art_Western6m","EC-arteriole_Western6m"), 
                         ident.2 = c("EC-art_chow6m","EC-arteriole_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "vis_Art_WB_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = vis, 
                         ident.1 = c("EC-ven_Western6m","EC-venule_Western6m"), 
                         ident.2 = c("EC-ven_chow6m","EC-venule_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "vis_Ven_WB_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = vis, 
                         ident.1 = c("EC-art_Rev3m","EC-arteriole_Rev3m"), 
                         ident.2 = c("EC-art_chow6m","EC-arteriole_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "vis_Art_Rev3m_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = vis, 
                         ident.1 = c("EC-ven_Rev3m","EC-venule_Rev3m"), 
                         ident.2 = c("EC-ven_chow6m","EC-venule_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "vis_Ven_Rev3m_vs_chow_6m.txt", sep = "\t")


####
# combine Cap/ Art/ Ven in 1 table

# 1. Capillaries
tabs = list.files()[grep("vis_Cap_", list.files())]
deg.capillaries = list()
deg.i = NULL
for (i in seq(1,4,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  # change the base
  deg.i$avg_logFC = deg.i$avg_log2FC/log(exp(1),2)
  print(head(deg.i))
  deg.capillaries[[i]] = deg.i
}

deg.capillaries[[5]] = read.table("/Users/olga.bondareva/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new/deg_gen/vis_EC-cap_Western_vs_chow.txt", sep = '\t', header = T)

head(deg.capillaries[[1]])

# combine a table of logFC
tab.deg.cap = merge(deg.capillaries[[5]][,c(2,5)], deg.capillaries[[3]][,c(6,5)], by="row.names", all=TRUE)
colnames(tab.deg.cap) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[4]][,c(6,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[1]][,c(6,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[2]][,c(6,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.cap)

rownames(tab.deg.cap) = tab.deg.cap$gene
tab.deg.cap = tab.deg.cap[,-1]
nrow(tab.deg.cap) # 7683 genes
write.table(tab.deg.cap,"table_deg_vis_Cap_all_TP.txt", sep='\t')


############
setwd(paste0(file.dir,"/vis/deg-gen"))
# 2. Artery
tabs = list.files()[grep("vis_Art_", list.files())]
deg.artery = list()
deg.i = NULL
for (i in seq(1,4,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  # change the base
  deg.i$avg_logFC = deg.i$avg_log2FC/log(exp(1),2)
  print(head(deg.i))
  deg.artery[[i]] = deg.i
}

deg.artery[[5]] = read.table("/Users/olga.bondareva/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new/deg_gen/vis_EC-art_Western_vs_chow.txt", sep = '\t', header = T)

head(deg.artery[[1]])

# combine a table of logFC
tab.deg.art = merge(deg.artery[[5]][,c(2,5)], deg.artery[[3]][,c(6,5)], by="row.names", all=TRUE)
colnames(tab.deg.art) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.art = merge(tab.deg.art, deg.artery[[4]][,c(6,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.art = merge(tab.deg.art, deg.artery[[1]][,c(6,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.art = merge(tab.deg.art, deg.artery[[2]][,c(6,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.art)

rownames(tab.deg.art) = tab.deg.art$gene
tab.deg.art = tab.deg.art[,-1]
nrow(tab.deg.art) # 6848 genes
write.table(tab.deg.art,"table_deg_vis_art_all_TP.txt", sep='\t')



############
# 3. Vein
tabs = list.files()[grep("vis_Ven_", list.files())]
deg.vein = list()
deg.i = NULL
for (i in seq(1,4,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  # change the base
  deg.i$avg_logFC = deg.i$avg_log2FC/log(exp(1),2)
  print(head(deg.i))
  deg.vein[[i]] = deg.i
}

deg.vein[[5]] = read.table("/Users/olga.bondareva/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new/deg_gen/vis_EC-ven_Western_vs_chow.txt", sep = '\t', header = T)

head(deg.vein[[1]])

# combine a table of logFC
tab.deg.ven = merge(deg.vein[[5]][,c(2,5)], deg.vein[[3]][,c(6,5)], by="row.names", all=TRUE)
colnames(tab.deg.ven) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.ven = merge(tab.deg.ven, deg.vein[[4]][,c(6,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.ven = merge(tab.deg.ven, deg.vein[[1]][,c(6,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.ven = merge(tab.deg.ven, deg.vein[[2]][,c(6,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.ven)

rownames(tab.deg.ven) = tab.deg.ven$gene
tab.deg.ven = tab.deg.ven[,-1]
nrow(tab.deg.ven) # 6848 genes
write.table(tab.deg.ven,"table_deg_vis_ven_all_TP.txt", sep='\t')




###########
VlnPlot(vis, features = c("Hspa1a", "Hspa1b","Hsp90aa1","Hspa5","Hspb1","Jun","Junb","Fos", "Egr1"), group.by = "diet_time", pt.size = 0)

saveRDS(vis, "vis_subset_3m_4m_6m.rds")




col.gen2 = c("#e30026","#fcacac","#ff7d61","#ffb5db","#cf91cc","#cf91b1",
             "#4287f5","#e02d8a","#46d130", 
             "#f09000", "#bed647","#47d6c1",
             "#3037c7", "#006887","#87b5ff","#ff5e5e",
             "#4bbd8f","#84bf86","#1f8500","#bf57f2")
names(col.gen2) = c("Art","Cap","Cap","Cap","Cap","Cap",
                    "Ven","Prolif","Lymph",
                    "Angiog","Antigen","Sp3",
                    "Sp4","Sp5","Venule","Arteriole",
                    "Lymph","Lymph","Lymph","Sp6")
col.gen2  
DimPlot(lung, reduction = "umap", label = TRUE, cols = as.vector(col.gen2[c(10,11,4,5,
                                                                            1,2,3,9,13,
                                                                            12,19,18,17,
                                                                            20,16,7,8)]), pt.size = 1)



################
###########
###
# Quantify the % of Hbb-positive cells per organ
# 1. In NormCounts
Hbb.cells = WhichCells(brain, slot = "data",expression = `Hbb-bs`>1 | `Hba-a1`>1 | `Hba-a2`>1 | `Hbb-bt`>1)
100*length(Hbb.cells)/ncol(brain)  # % of Hbb-bs cells
# Any Hb 19.76261 %
DotPlot(brain, group.by = "diet_time", features = c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt"), 
        cols = c("blue","red"))

Hbb.cells = WhichCells(heart, slot = "data",expression = `Hbb-bs`>1 | `Hba-a1`>1 | `Hba-a2`>1 | `Hbb-bt`>1)
100*length(Hbb.cells)/ncol(heart)  # % of Hbb-cells
# heart 12.5066 %

Hbb.cells = WhichCells(lung, slot = "data",expression = `Hbb-bs`>1 | `Hba-a1`>1 | `Hba-a2`>1 | `Hbb-bt`>1)
100*length(Hbb.cells)/ncol(lung)  # % of Hbb-cells
# lung 2.592805 %

Hbb.cells = WhichCells(liver, slot = "data",expression = `Hbb-bs`>1 | `Hba-a1`>1 | `Hba-a2`>1 | `Hbb-bt`>1)
100*length(Hbb.cells)/ncol(liver)  # % of Hbb-cells
# liver 1.233211 %

Hbb.cells = WhichCells(kidney, slot = "data",expression = `Hbb-bs`>1 | `Hba-a1`>1 | `Hba-a2`>1 | `Hbb-bt`>1)
100*length(Hbb.cells)/ncol(kidney)  # % of Hbb-cells
# kidney 0.5688889 %

Hbb.cells = WhichCells(vis, slot = "data",expression = `Hbb-bs`>1 | `Hba-a1`>1 | `Hba-a2`>1 | `Hbb-bt`>1)
100*length(Hbb.cells)/ncol(vis)  # % of Hbb-cells
# vis 2.540071 %

Hbb.cells = WhichCells(sc, slot = "data",expression = `Hbb-bs`>1 | `Hba-a1`>1 | `Hba-a2`>1 | `Hbb-bt`>1)
100*length(Hbb.cells)/ncol(sc)  # % of Hbb-cells
# sc 4.847016 %
DotPlot(sc, group.by = "diet", features = c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt"), 
        cols = c("blue","red"))

# 1. In RawCounts
Hbb.cells = WhichCells(brain, slot = "counts",expression = `Hbb-bs`>1 | `Hba-a1`>1 | `Hba-a2`>1 | `Hbb-bt`>1)
100*length(Hbb.cells)/ncol(brain)  # % of Hbb-bs cells
# Any Hb 5.749035 %
table(brain$diet_time[Hbb.cells])
table(brain$diet_time)

Hbb.cells = WhichCells(heart, slot = "counts",expression = `Hbb-bs`>1 | `Hba-a1`>1 | `Hba-a2`>1 | `Hbb-bt`>1)
100*length(Hbb.cells)/ncol(heart)  # % of Hbb-cells
# heart 1.75068 %
table(heart$diet_time[Hbb.cells])
table(heart$diet_time)

Hbb.cells = WhichCells(lung, slot = "counts",expression = `Hbb-bs`>1 | `Hba-a1`>1 | `Hba-a2`>1 | `Hbb-bt`>1)
100*length(Hbb.cells)/ncol(lung)  # % of Hbb-cells
# lung 0.5414693 %
table(lung$diet_time[Hbb.cells])
table(lung$diet_time)

Hbb.cells = WhichCells(liver, slot = "counts",expression = `Hbb-bs`>1 | `Hba-a1`>1 | `Hba-a2`>1 | `Hbb-bt`>1)
100*length(Hbb.cells)/ncol(liver)  # % of Hbb-cells
# liver 0.5616606 %
table(liver$diet_time[Hbb.cells])
table(liver$diet_time)

Hbb.cells = WhichCells(kidney, slot = "counts",expression = `Hbb-bs`>1 | `Hba-a1`>1 | `Hba-a2`>1 | `Hbb-bt`>1)
100*length(Hbb.cells)/ncol(kidney)  # % of Hbb-cells
# kidney 0.04740741 %
table(kidney$diet_time[Hbb.cells])
table(kidney$diet_time)

Hbb.cells = WhichCells(vis, slot = "counts",expression = `Hbb-bs`>1 | `Hba-a1`>1 | `Hba-a2`>1 | `Hbb-bt`>1)
100*length(Hbb.cells)/ncol(vis)  # % of Hbb-cells
# vis 0.9372453 %
table(vis$diet_time[Hbb.cells])
table(vis$diet)

Hbb.cells = WhichCells(sc, slot = "counts",expression = `Hbb-bs`>1 | `Hba-a1`>1 | `Hba-a2`>1 | `Hbb-bt`>1)
100*length(Hbb.cells)/ncol(sc)  # % of Hbb-cells
# sc 1.312734 %
table(sc$diet[Hbb.cells])
table(sc$diet_time)

