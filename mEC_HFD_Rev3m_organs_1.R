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
setwd("~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_06_09_M_HFD_Rev_3m/analysis/organs")
file.dir = "~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_06_09_M_HFD_Rev_3m/analysis/organs"
data.dir = "~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_06_09_M_HFD_Rev_3m/data"
data_dirs = list.files(data.dir)
data_dirs

data.dir.3m = "~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/organs-3rd-processing"
data.dir.4m = "~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_04_15_M_HFD_Rev_1m/analysis/combined/"
data.dir.6m = "~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_06_09_M_HFD_Rev_3m/analysis/combined/"

########################################
# analyze organs separately
########################################
#########################
# 3. brain
#########################
setwd("brain")

brain.3m = readRDS(paste0(data.dir.3m,"/brain/brain_ECs_singlet.rds"))
brain.4m = readRDS(paste0(data.dir.4m,"brain_subset.rds"))
brain.6m = readRDS(paste0(data.dir.6m,"brain_subset.rds"))

brain.3m$time = "3m"

Idents(brain.3m) = brain.3m$brain_celltype
DimPlot(brain.3m, reduction = "umap", label = TRUE)
#DimPlot(brain.4m, reduction = "umap", label = TRUE)
#DimPlot(brain.6m, reduction = "umap", label = TRUE)

# merge objects
brain = merge(brain.3m, list(brain.4m, brain.6m))
rm(brain.3m, brain.4m, brain.6m)

brain <- ScaleData(object = brain, verbose = T)
brain <- FindVariableFeatures(brain, selection.method = "vst", nfeatures = 3000)
brain <- RunPCA(object = brain, npcs = 30, verbose = T)
ElbowPlot(object = brain,  ndims = 30) 

# UMAP and Clustering
brain <- RunUMAP(object = brain, reduction = "pca", dims = 1:26)
brain <- FindNeighbors(object = brain, reduction = "pca", dims = 1:26)
brain <- FindClusters(object = brain, resolution = 0.4)

brain$diet_time = paste0(brain$diet, brain$time)

# visualization
DimPlot(object = brain, reduction = "umap", label = TRUE)
DimPlot(object = brain, reduction = "umap", group.by = "orig.ident")
DimPlot(object = brain, reduction = "umap", group.by = "brain_celltype", label = T) + NoLegend()
DimPlot(object = brain, reduction = "umap", split.by = "diet_time", ncol = 3, pt.size = 1.3)
DimPlot(object = brain, reduction = "umap", group.by = "brain_celltype", label = T,
        split.by = "diet_time", ncol = 3)


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

# Pericyte amd SMC markers
FeaturePlot(brain, features = c("Flt1", "Pecam1","Cdh5",# EC
                                "Pdgfrb", "Cspg4","Anpep",  #  Pericyte
                                "Myh11", "Acta2", "Des",   # SMC
                                "Dcn", "Col1a1", "Pdgfra"), #FB
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10") #, min.cutoff = "q10"

FeaturePlot(brain, features = c("Prox1", "Pdpn","Lyve1", "Mki67", "Cenpf"), # lymph
            order = T, cols=c("grey", "red"), ncol=2, min.cutoff = "q10") #, min.cutoff = "q10"

#  tip cell-enriched, angiogenic
FeaturePlot(brain, features = c("Apln", "Col4a2", "Trp53i11"),
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 3)

# select atrerial cells
plot <- DimPlot(brain, reduction = "umap", cells = WhichCells(brain, idents = "6"))
select.cells <- CellSelector(plot = plot) # 
Idents(brain, cells = select.cells) <- "EC-art"

brain <- RenameIdents(brain, `5` = "EC-ven",`0` = "EC-cap1",`1` = "EC-cap2",#`6` = "EC-art",
                      `3` = "EC-capA1",`2` = "EC-capA2",`4` = "EC-capA3",`9` = "EC-platelet",
                      `10` = "EC-fenestr",`7` = "EC-Hb",`6` = "EC-AP1",
                      `8` = "EC-Myl9") 
## Relevel object@ident
brain@active.ident <- factor(x = brain@active.ident, levels = sort(levels(Idents(brain))))
brain$brain_celltype2 <- Idents(brain)
DimPlot(brain, reduction = "umap", label = TRUE)

# number of cells per cluster per sample
table(Idents(brain), brain$diet_time)

DimPlot(brain, reduction = "umap", label = TRUE, split.by = "diet", ncol = 2) + NoLegend()
DimPlot(object = brain, reduction = "umap", group.by = "brain_celltype2", label = T,
        split.by = "diet_time", ncol = 3)+ NoLegend()



brain = readRDS("brain_3m_4m_6m.rds")

# DE genes clusters
setwd("deg-clust")
brain$clust.diet_time <- paste(brain$brain_celltype2, brain$diet_time, sep = "_")
Idents(object = brain) <- "clust.diet_time"
brain@active.ident <- factor(x = brain@active.ident, levels = sort(levels(Idents(object = brain))))
levels(Idents(object = brain))

lev1 <- levels(brain@active.ident)[seq(7,(length(levels(Idents(brain)))),8)] # condition 1=HFD_4m
lev2 <- levels(brain@active.ident)[seq(2,(length(levels(Idents(brain)))),8)] # condition 2=chow_4m
lev3 <- levels(brain@active.ident)[seq(4,(length(levels(Idents(brain)))),8)] # condition 3=Rev_1m
lev4 <- levels(brain@active.ident)[seq(8,(length(levels(Idents(brain)))),8)] # condition 4=HFD_6m
lev5 <- levels(brain@active.ident)[seq(3,(length(levels(Idents(brain)))),8)] # condition 5=chow_6m
lev6 <- levels(brain@active.ident)[seq(5,(length(levels(Idents(brain)))),8)] # condition 6=Rev_3m
lev7 <- levels(brain@active.ident)[seq(6,(length(levels(Idents(brain)))),8)] # condition 7=HFD_3m
lev8 <- levels(brain@active.ident)[seq(1,(length(levels(Idents(brain)))),8)] # condition 8=chow_3m

# for all the clusters
for (i in seq(1,length(lev1),1)) {
  # 4 months
  # Western
  print(lev1[i])
  print(lev2[i])
  response2 <- FindMarkers(object = brain, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # Reverse
  print(lev3[i])
  print(lev2[i])
  response2 <- FindMarkers(object = brain, ident.1 = lev3[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev3[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # 6 months
  # Western
  print(lev4[i])
  print(lev5[i])
  response2 <- FindMarkers(object = brain, ident.1 = lev4[i], 
                           ident.2 = lev5[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev4[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # Reverse
  print(lev6[i])
  print(lev5[i])
  response2 <- FindMarkers(object = brain, ident.1 = lev6[i], 
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

####

# DEG for capillaries combined
# 4m
response2 <- FindMarkers(object = brain, 
                         ident.1 = c("EC-AP1_Western4m","EC-cap1_Western4m","EC-cap2_Western4m",
                                     "EC-capA1_Western4m","EC-capA2_Western4m","EC-capA3_Western4m"), 
                         ident.2 = c("EC-AP1_chow4m","EC-cap1_chow4m","EC-cap2_chow4m",
                                     "EC-capA1_chow4m","EC-capA2_chow4m","EC-capA3_chow4m"), 
                         logfc.threshold = 0,
                         verbose = FALSE) #, min.pct = 0.25
write.table(response2, file = "brain_Cap_WB_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = brain, 
                         ident.1 = c("EC-AP1_Rev1m","EC-cap1_Rev1m","EC-cap2_Rev1m",
                                     "EC-capA1_Rev1m","EC-capA2_Rev1m","EC-capA3_Rev1m"), 
                         ident.2 = c("EC-AP1_chow4m","EC-cap1_chow4m","EC-cap2_chow4m",
                                     "EC-capA1_chow4m","EC-capA2_chow4m","EC-capA3_chow4m"), 
                         logfc.threshold = 0,
                         verbose = FALSE) #, min.pct = 0.25
write.table(response2, file = "brain_Cap_Rev1m_vs_chow_4m.txt", sep = "\t")

# 6m
response2 <- FindMarkers(object = brain, 
                         ident.1 = c("EC-AP1_Western6m","EC-cap1_Western6m","EC-cap2_Western6m",
                                     "EC-capA1_Western6m","EC-capA2_Western6m","EC-capA3_Western6m"), 
                         ident.2 = c("EC-AP1_chow6m","EC-cap1_chow6m","EC-cap2_chow6m",
                                     "EC-capA1_chow6m","EC-capA2_chow6m","EC-capA3_chow6m"), 
                         logfc.threshold = 0,
                         verbose = FALSE) #, min.pct = 0.25
write.table(response2, file = "brain_Cap_WB_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = brain, 
                         ident.1 = c("EC-AP1_Rev3m","EC-cap1_Rev3m","EC-cap2_Rev3m",
                                     "EC-capA1_Rev3m","EC-capA2_Rev3m","EC-capA3_Rev3m"), 
                         ident.2 = c("EC-AP1_chow6m","EC-cap1_chow6m","EC-cap2_chow6m",
                                     "EC-capA1_chow6m","EC-capA2_chow6m","EC-capA3_chow6m"), 
                         logfc.threshold = 0,
                         verbose = FALSE) #, min.pct = 0.25
write.table(response2, file = "brain_Cap_Rev3m_vs_chow_6m.txt", sep = "\t")

# DEG for vein and artery combined
# 4m
response2 <- FindMarkers(object = brain, 
                         ident.1 = c("EC-art_Western4m"), 
                         ident.2 = c("EC-art_chow4m"), 
                         logfc.threshold = 0,
                         verbose = FALSE) #, min.pct = 0.25
write.table(response2, file = "brain_Art_WB_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = brain, 
                         ident.1 = c("EC-ven_Western4m"), 
                         ident.2 = c("EC-ven_chow4m"), 
                         logfc.threshold = 0,
                         verbose = FALSE) #, min.pct = 0.25
write.table(response2, file = "brain_Ven_WB_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = brain, 
                         ident.1 = c("EC-art_Rev1m"), 
                         ident.2 = c("EC-art_chow4m"), 
                         logfc.threshold = 0,
                         verbose = FALSE) #, min.pct = 0.25
write.table(response2, file = "brain_Art_Rev1m_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = brain, 
                         ident.1 = c("EC-ven_Rev1m"), 
                         ident.2 = c("EC-ven_chow4m"), 
                         logfc.threshold = 0,
                         verbose = FALSE) #, min.pct = 0.25
write.table(response2, file = "brain_Ven_Rev1m_vs_chow_4m.txt", sep = "\t")

# 6m
response2 <- FindMarkers(object = brain, 
                         ident.1 = c("EC-art_Western6m"), 
                         ident.2 = c("EC-art_chow6m"), 
                         logfc.threshold = 0,
                         verbose = FALSE) #, min.pct = 0.25
write.table(response2, file = "brain_Art_WB_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = brain, 
                         ident.1 = c("EC-ven_Western6m"), 
                         ident.2 = c("EC-ven_chow6m"), 
                         logfc.threshold = 0,
                         verbose = FALSE) #, min.pct = 0.25
write.table(response2, file = "brain_Ven_WB_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = brain, 
                         ident.1 = c("EC-art_Rev3m"), 
                         ident.2 = c("EC-art_chow6m"), 
                         logfc.threshold = 0,
                         verbose = FALSE) #, min.pct = 0.25
write.table(response2, file = "brain_Art_Rev3m_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = brain, 
                         ident.1 = c("EC-ven_Rev3m"), 
                         ident.2 = c("EC-ven_chow6m"), 
                         logfc.threshold = 0,
                         verbose = FALSE) #, min.pct = 0.25
write.table(response2, file = "brain_Ven_Rev3m_vs_chow_6m.txt", sep = "\t")


saveRDS(brain, "brain_3m_4m_6m.rds")

####
# combine Cap/ Art/ Ven in 1 table

# 1. Capillaries
tabs = list.files()[grep("brain_Cap_", list.files())]
deg.capillaries = list()
deg.i = NULL
for (i in seq(1,4,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.capillaries[[i]] = deg.i
}

deg.capillaries[[5]] = read.table("/Users/olga.bondareva/Desktop/HI-MAG OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new/deg_gen/brain_EC-cap_Western_vs_chow.txt", sep = '\t', header = T)

head(deg.capillaries[[1]])

# combine a table of logFC
tab.deg.cap = merge(deg.capillaries[[5]][,c(2,5)], deg.capillaries[[3]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.cap) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[4]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[1]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.cap)

rownames(tab.deg.cap) = tab.deg.cap$gene
tab.deg.cap = tab.deg.cap[,-1]
nrow(tab.deg.cap) # 7842 genes
write.table(tab.deg.cap,"table_deg_brain_Cap_all_TP.txt", sep='\t')

# signif DEGs
tab.deg.cap.sign = tab.deg.cap[(tab.deg.cap$WD_3m_Q.val<0.05)|(tab.deg.cap$WD_4m_Q.val<0.05)|(tab.deg.cap$WD_6m_Q.val<0.05)|(tab.deg.cap$Rev_1m_Q.val<0.05)|(tab.deg.cap$Rev_3m_Q.val<0.05),]
# rm NAs
tab.deg.cap.sign = tab.deg.cap.sign[complete.cases(tab.deg.cap.sign),]
nrow(tab.deg.cap.sign) # 1657 genes with signif Q-val
head(tab.deg.cap.sign)
write.table(tab.deg.cap.sign, "table_deg_brain_Cap_all_TP_sign.txt", sep = '\t')

# signif only at 3m TP
tab.deg.cap.sign = tab.deg.cap[(tab.deg.cap$WD_3m_Q.val<0.05),]
# rm NAs
tab.deg.cap.sign = tab.deg.cap.sign[complete.cases(tab.deg.cap.sign),]
nrow(tab.deg.cap.sign) # 274 genes with signif Q-val
head(tab.deg.cap.sign)

# keep only logFC
tab.deg.cap.logFC = tab.deg.cap.sign[,seq(1,10,2)] # tab.deg.cap[,seq(1,14,2)]
head(tab.deg.cap.logFC)

# replace NA with 0s
tab.deg.cap.logFC[is.na(tab.deg.cap.logFC)] = 0
head(tab.deg.cap.logFC)

# remove rows with all 0s
tab.deg.cap.logFC = tab.deg.cap.logFC[!rowSums(tab.deg.cap.logFC)==0,]
head(tab.deg.cap.logFC)
nrow(tab.deg.cap.logFC) # 4237 genes with logFC
colnames(tab.deg.cap.logFC) = c("WD_3m","WD_4m","WD_6m","Rev_1m","Rev_3m")
write.table(tab.deg.cap.logFC, "table_deg_brain_Cap_all_TP_sign_LOG.txt", sep = '\t')

b <- tab.deg.cap.logFC
# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]
# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5
mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma
ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
          main = "Brain Cap DEG, logFC(x/chow)", 
          Rowv= T,Colv= T, 
          ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="both", cexRow=0.8, cexCol=0.8) 

# plot correlations

pdf(paste0("corr-cap-brain.pdf"),width=10,height=10,paper='special')
p1 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,WD_4m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_4m") + theme_linedraw() + geom_smooth(method = "lm")
p2 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,WD_6m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_6m") + theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,Rev_1m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_1m") + theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,Rev_3m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_31m") + theme_linedraw() + geom_smooth(method = "lm")

print((p1 | p2 ) /
        ( p3 |p4)) # 10 x 10
dev.off()

# all clusters
library(corrplot)
m <- cor(tab.deg.cap.logFC)
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black",
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))

############
setwd(paste0(file.dir,"/brain/deg-gen"))
# 2. Artery
tabs = list.files()[grep("brain_Art_", list.files())]
deg.artery = list()
deg.i = NULL
for (i in seq(1,4,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.artery[[i]] = deg.i
}

deg.artery[[5]] = read.table("/Users/olga.bondareva/Desktop/HI-MAG OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new/deg_gen/brain_EC-art_Western_vs_chow.txt", sep = '\t', header = T)

head(deg.artery[[1]])

# combine a table of logFC
tab.deg.art = merge(deg.artery[[5]][,c(2,5)], deg.artery[[3]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.art) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.art = merge(tab.deg.art, deg.artery[[4]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.art = merge(tab.deg.art, deg.artery[[1]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.art = merge(tab.deg.art, deg.artery[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.art)

rownames(tab.deg.art) = tab.deg.art$gene
tab.deg.art = tab.deg.art[,-1]
nrow(tab.deg.art) # 6848 genes
write.table(tab.deg.art,"table_deg_brain_art_all_TP.txt", sep='\t')

# signif DEGs
tab.deg.art.sign = tab.deg.art[(tab.deg.art$WD_3m_Q.val<0.05)|(tab.deg.art$WD_4m_Q.val<0.05)|(tab.deg.art$WD_6m_Q.val<0.05)|(tab.deg.art$Rev_1m_Q.val<0.05)|(tab.deg.art$Rev_3m_Q.val<0.05),]
# rm NAs
tab.deg.art.sign = tab.deg.art.sign[complete.cases(tab.deg.art.sign),]
nrow(tab.deg.art.sign) # 1671 genes with signif Q-val
head(tab.deg.art.sign)
write.table(tab.deg.art.sign, "table_deg_brain_art_all_TP_sign.txt", sep = '\t')

# keep only logFC
tab.deg.art.logFC = tab.deg.art.sign[,seq(1,10,2)] # tab.deg.art[,seq(1,14,2)]
head(tab.deg.art.logFC)

# replace NA with 0s
tab.deg.art.logFC[is.na(tab.deg.art.logFC)] = 0
head(tab.deg.art.logFC)

# remove rows with all 0s
tab.deg.art.logFC = tab.deg.art.logFC[!rowSums(tab.deg.art.logFC)==0,]
head(tab.deg.art.logFC)
nrow(tab.deg.art.logFC) # 591 genes with logFC
colnames(tab.deg.art.logFC) = c("WD_3m","WD_4m","WD_6m","Rev_1m","Rev_3m")
write.table(tab.deg.art.logFC, "table_deg_brain_art_all_TP_sign_LOG.txt", sep = '\t')

b <- tab.deg.art.logFC
# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]
# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5
mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma
ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))

pdf(paste0("heatmap-deg-art-sign-brain.pdf"),width=10,height=15,paper='special')
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
          main = "brain art DEG, logFC(x/chow)", 
          Rowv= T,Colv= T, 
          ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="both", cexRow=0.8, cexCol=0.8) 
dev.off()
# plot correlations

pdf(paste0("corr-art-brain.pdf"),width=10,height=10,paper='special')
p1 <- ggplot(tab.deg.art.logFC, aes(WD_3m,WD_4m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_4m") + theme_linedraw() + geom_smooth(method = "lm")
p2 <- ggplot(tab.deg.art.logFC, aes(WD_3m,WD_6m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_6m") + theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.art.logFC, aes(WD_3m,Rev_1m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_1m") + theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.art.logFC, aes(WD_3m,Rev_3m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_31m") + theme_linedraw() + geom_smooth(method = "lm")

print((p1 | p2 ) /
        ( p3 |p4)) # 10 x 10
dev.off()

# all clusters
m <- cor(tab.deg.art.logFC)
pdf(paste0("corr-art-tab-brain.pdf"),width=7,height=7,paper='special')
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black",
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))
dev.off()

############
# 3. Vein
tabs = list.files()[grep("brain_Ven_", list.files())]
deg.vein = list()
deg.i = NULL
for (i in seq(1,4,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.vein[[i]] = deg.i
}

deg.vein[[5]] = read.table("/Users/olga.bondareva/Desktop/HI-MAG OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new/deg_gen/brain_EC-ven_Western_vs_chow.txt", sep = '\t', header = T)

head(deg.vein[[1]])

# combine a table of logFC
tab.deg.ven = merge(deg.vein[[5]][,c(2,5)], deg.vein[[3]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.ven) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.ven = merge(tab.deg.ven, deg.vein[[4]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.ven = merge(tab.deg.ven, deg.vein[[1]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.ven = merge(tab.deg.ven, deg.vein[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.ven)

rownames(tab.deg.ven) = tab.deg.ven$gene
tab.deg.ven = tab.deg.ven[,-1]
nrow(tab.deg.ven) # 6848 genes
write.table(tab.deg.ven,"table_deg_brain_ven_all_TP.txt", sep='\t')

# signif DEGs
tab.deg.ven.sign = tab.deg.ven[(tab.deg.ven$WD_3m_Q.val<0.05)|(tab.deg.ven$WD_4m_Q.val<0.05)|(tab.deg.ven$WD_6m_Q.val<0.05)|(tab.deg.ven$Rev_1m_Q.val<0.05)|(tab.deg.ven$Rev_3m_Q.val<0.05),]
# rm NAs
tab.deg.ven.sign = tab.deg.ven.sign[complete.cases(tab.deg.ven.sign),]
nrow(tab.deg.ven.sign) # 1671 genes with signif Q-val
head(tab.deg.ven.sign)
write.table(tab.deg.ven.sign, "table_deg_brain_ven_all_TP_sign.txt", sep = '\t')

# keep only logFC
tab.deg.ven.logFC = tab.deg.ven.sign[,seq(1,10,2)] # tab.deg.ven[,seq(1,14,2)]
head(tab.deg.ven.logFC)

# replace NA with 0s
tab.deg.ven.logFC[is.na(tab.deg.ven.logFC)] = 0
head(tab.deg.ven.logFC)

# remove rows with all 0s
tab.deg.ven.logFC = tab.deg.ven.logFC[!rowSums(tab.deg.ven.logFC)==0,]
head(tab.deg.ven.logFC)
nrow(tab.deg.ven.logFC) # 591 genes with logFC
colnames(tab.deg.ven.logFC) = c("WD_3m","WD_4m","WD_6m","Rev_1m","Rev_3m")
write.table(tab.deg.ven.logFC, "table_deg_brain_ven_all_TP_sign_LOG.txt", sep = '\t')

b <- tab.deg.ven.logFC
# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]
# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5
mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma
ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))

pdf(paste0("heatmap-deg-ven-sign-brain.pdf"),width=10,height=15,paper='special')
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
          main = "brain ven DEG, logFC(x/chow)", 
          Rowv= T,Colv= T, 
          ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="both", cexRow=0.8, cexCol=0.8) 
dev.off()
# plot correlations

pdf(paste0("corr-ven-brain.pdf"),width=10,height=10,paper='special')
p1 <- ggplot(tab.deg.ven.logFC, aes(WD_3m,WD_4m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_4m") + theme_linedraw() + geom_smooth(method = "lm")
p2 <- ggplot(tab.deg.ven.logFC, aes(WD_3m,WD_6m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_6m") + theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.ven.logFC, aes(WD_3m,Rev_1m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_1m") + theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.ven.logFC, aes(WD_3m,Rev_3m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_31m") + theme_linedraw() + geom_smooth(method = "lm")

print((p1 | p2 ) /
        ( p3 |p4)) # 10 x 10
dev.off()

# all clusters
m <- cor(tab.deg.ven.logFC)
pdf(paste0("corr-ven-tab-brain.pdf"),width=7,height=7,paper='special')
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black",
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))
dev.off()

### fenestrated EC
setwd("brain/deg-clust")
tabs = list.files()[grep("EC-fenestr_", list.files())]
deg.capillaries = list()
deg.i = NULL
for (i in seq(1,5,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.capillaries[[i]] = deg.i
}

head(deg.capillaries[[1]])

# combine a table of logFC
tab.deg.cap = merge(deg.capillaries[[3]][,c(2,5)], deg.capillaries[[4]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.cap) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[5]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[1]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.cap)

rownames(tab.deg.cap) = tab.deg.cap$gene
tab.deg.cap = tab.deg.cap[,-1]
nrow(tab.deg.cap) # 8950 genes
write.table(tab.deg.cap,"table_deg_brain_fenestr_all_TP.txt", sep='\t')




############
# 1. Liver
############
setwd("liver")
liver.3m = readRDS(paste0(data.dir.3m,"/liver/liver_ECs_singlet.rds"))
liver.4m = readRDS(paste0(data.dir.4m,"liver_subset.rds"))
liver.6m = readRDS(paste0(data.dir.6m,"liver_subset.rds"))

liver.3m$time = "3m"

# merge objects
liver = merge(liver.3m, list(liver.4m, liver.6m))
rm(liver.3m, liver.4m, liver.6m)

liver$diet_time = paste0(liver$diet, liver$time)

liver <- ScaleData(object = liver, verbose = T)
liver <- FindVariableFeatures(liver, selection.method = "vst", nfeatures = 3000)
liver <- RunPCA(object = liver, npcs = 30, verbose = T)
ElbowPlot(object = liver,  ndims = 30) 

# UMAP and Clustering
liver <- RunUMAP(object = liver, reduction = "pca", dims = 1:28)
liver <- FindNeighbors(object = liver, reduction = "pca", dims = 1:28)
liver <- FindClusters(object = liver, resolution = 0.3)

# Visualization
DimPlot(object = liver, reduction = "umap", label = TRUE)
DimPlot(object = liver, reduction = "umap", split.by = "diet_time", pt.size = 1.3)
DimPlot(object = liver, reduction = "umap", group.by = "liver_celltype", label = T)

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

FeaturePlot(liver, features = c("Prox1", "Pdpn","Lyve1", # lymph
                                "Cenpf", "Mki67"), # prolif
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10") #, min.cutoff = "q10"

#  tip cell-enriched, angiogenic
FeaturePlot(liver, features = c("Apln", "Col4a2", "Trp53i11"),
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 3)

# Find marker genes
Idents(liver) = liver$seurat_clusters
liver.markers = FindAllMarkers(liver, only.pos = T)
write.table(liver.markers, "liver_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(liver.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "liver_top5_cluster_markers.txt", sep='\t')
DotPlot(liver, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()

# Rename
liver <- RenameIdents(liver, `0` = "EC-cap1",`1` = "EC-cap2",`2` = "EC-cap3",`11` = "EC-AP1",
                      `7` = "EC-art",`3` = "EC-ven",
                      `10` = "EC-Hb",`8` = "EC-platelet",
                      `12` = "Immune1",`13` = "Immune2",`9` = "Immune3",
                      `4` = "EC-liver1",`5` = "EC-liver2",`6` = "Prolif") 

## Relevel object@ident
liver@active.ident <- factor(x = liver@active.ident, levels = sort(levels(Idents(liver))))
liver$liver_celltype2 <- Idents(liver)
DimPlot(liver, reduction = "umap", group.by = "liver_celltype2",label = TRUE)



# number of cells per cluster per sample
table(Idents(liver), liver$diet_time)

# remove immune cells
liver = subset(liver, idents = c("Immune1", "Immune2","Immune3"), invert = T)

# select cluster by hand
plot <- DimPlot(liver, reduction = "umap")
select.cells <- CellSelector(plot = plot) # 
Idents(liver, cells = select.cells) <- "outlier"
liver = subset(liver, idents = c("outlier"), invert = T)

DimPlot(object = liver, reduction = "umap", group.by = "liver_celltype2", split.by = "diet_time", label = T, ncol =3, pt.size = 1.3) + NoLegend()

# 15 x 11

liver = readRDS("liver_3m_4m_6m.rds")
# DE genes
liver$clust.diet_time <- paste(liver$liver_celltype2, liver$diet_time, sep = "_")
Idents(liver) <- "clust.diet_time"
liver@active.ident <- factor(x = liver@active.ident, levels = sort(levels(Idents(liver))))
l = levels(Idents(liver))[-seq(1,7,1)]

lev1 <- l[seq(7,(length(l)),8)] # condition 1=HFD_4m
lev2 <- l[seq(2,(length(l)),8)] # condition 2=chow_4m
lev3 <- l[seq(4,(length(l)),8)] # condition 3=Rev_1m
lev4 <- l[seq(8,(length(l)),8)] # condition 4=HFD_6m
lev5 <- l[seq(3,(length(l)),8)] # condition 5=chow_6m
lev6 <- l[seq(5,(length(l)),8)] # condition 6=Rev_3m
lev7 <- l[seq(1,(length(l)),8)] # condition 7=HFD_3m
lev8 <- l[seq(6,(length(l)),8)] # condition 8=chow_3m

# for all the clusters
for (i in rev(seq(1,length(lev1),1))) {
  # 4 months
  # Western
  print(lev1[i])
  print(lev2[i])
  response2 <- FindMarkers(object = liver, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # Reverse
  print(lev3[i])
  print(lev2[i])
  response2 <- FindMarkers(object = liver, ident.1 = lev3[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev3[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # 6 months
  # Western
  print(lev4[i])
  print(lev5[i])
  response2 <- FindMarkers(object = liver, ident.1 = lev4[i], 
                           ident.2 = lev5[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev4[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # Reverse
  print(lev6[i])
  print(lev5[i])
  response2 <- FindMarkers(object = liver, ident.1 = lev6[i], 
                           ident.2 = lev5[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev6[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # 3 months
  # Western
  print(lev7[i])
  print(lev8[i])
  response2 <- FindMarkers(object = liver, ident.1 = lev7[i], 
                           ident.2 = lev8[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev7[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
}


# DEG for capillaries combined
# 4m
response2 <- FindMarkers(object = liver, 
                         ident.1 = c("EC-AP1_Western4m","EC-cap1_Western4m","EC-cap2_Western4m",
                                     "EC-cap3_Western4m"), 
                         ident.2 = c("EC-AP1_chow4m","EC-cap1_chow4m","EC-cap2_chow4m",
                                     "EC-cap3_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "liver_Cap_WB_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = liver, 
                         ident.1 = c("EC-AP1_Rev1m","EC-cap1_Rev1m","EC-cap2_Rev1m",
                                     "EC-cap3_Rev1m"), 
                         ident.2 = c("EC-AP1_chow4m","EC-cap1_chow4m","EC-cap2_chow4m",
                                     "EC-cap3_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "liver_Cap_Rev1m_vs_chow_4m.txt", sep = "\t")

# 6m
response2 <- FindMarkers(object = liver, 
                         ident.1 = c("EC-AP1_Western6m","EC-cap1_Western6m","EC-cap2_Western6m",
                                     "EC-cap3_Western6m"), 
                         ident.2 = c("EC-AP1_chow6m","EC-cap1_chow6m","EC-cap2_chow6m",
                                     "EC-cap3_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "liver_Cap_WB_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = liver, 
                         ident.1 = c("EC-cap1_Rev3m","EC-cap2_Rev3m",
                                     "EC-cap3_Rev3m"), 
                         ident.2 = c("EC-AP1_chow6m","EC-cap1_chow6m","EC-cap2_chow6m",
                                     "EC-cap3_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "liver_Cap_Rev3m_vs_chow_6m.txt", sep = "\t")

# DEG for vein and artery combined
# 4m
response2 <- FindMarkers(object = liver, 
                         ident.1 = c("EC-art_Western4m"), 
                         ident.2 = c("EC-art_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "liver_Art_WB_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = liver, 
                         ident.1 = c("EC-ven_Western4m"), 
                         ident.2 = c("EC-ven_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "liver_Ven_WB_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = liver, 
                         ident.1 = c("EC-art_Rev1m"), 
                         ident.2 = c("EC-art_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "liver_Art_Rev1m_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = liver, 
                         ident.1 = c("EC-ven_Rev1m"), 
                         ident.2 = c("EC-ven_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "liver_Ven_Rev1m_vs_chow_4m.txt", sep = "\t")

# 6m
response2 <- FindMarkers(object = liver, 
                         ident.1 = c("EC-art_Western6m"), 
                         ident.2 = c("EC-art_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "liver_Art_WB_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = liver, 
                         ident.1 = c("EC-ven_Western6m"), 
                         ident.2 = c("EC-ven_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "liver_Ven_WB_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = liver, 
                         ident.1 = c("EC-art_Rev3m"), 
                         ident.2 = c("EC-art_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "liver_Art_Rev3m_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = liver, 
                         ident.1 = c("EC-ven_Rev3m"), 
                         ident.2 = c("EC-ven_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "liver_Ven_Rev3m_vs_chow_6m.txt", sep = "\t")


####
# combine Cap/ Art/ Ven in 1 table

# 1. Capillaries
tabs = list.files()[grep("liver_Cap_", list.files())]
deg.capillaries = list()
deg.i = NULL
for (i in seq(1,4,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.capillaries[[i]] = deg.i
}

deg.capillaries[[5]] = read.table("/Users/olga.bondareva/Desktop/HI-MAG OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new/deg_gen/liver_EC-cap_Western_vs_chow.txt", sep = '\t', header = T)

head(deg.capillaries[[1]])

# combine a table of logFC
tab.deg.cap = merge(deg.capillaries[[5]][,c(2,5)], deg.capillaries[[3]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.cap) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[4]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[1]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.cap)

rownames(tab.deg.cap) = tab.deg.cap$gene
tab.deg.cap = tab.deg.cap[,-1]
nrow(tab.deg.cap) # 7842 genes
write.table(tab.deg.cap,"table_deg_liver_Cap_all_TP.txt", sep='\t')

# signif DEGs
tab.deg.cap.sign = tab.deg.cap[(tab.deg.cap$WD_3m_Q.val<0.05)|(tab.deg.cap$WD_4m_Q.val<0.05)|(tab.deg.cap$WD_6m_Q.val<0.05)|(tab.deg.cap$Rev_1m_Q.val<0.05)|(tab.deg.cap$Rev_3m_Q.val<0.05),]
# rm NAs
tab.deg.cap.sign = tab.deg.cap.sign[complete.cases(tab.deg.cap.sign),]
nrow(tab.deg.cap.sign) # 1657 genes with signif Q-val
head(tab.deg.cap.sign)
write.table(tab.deg.cap.sign, "table_deg_liver_Cap_all_TP_sign.txt", sep = '\t')

#####
# signif only at 3m TP
tab.deg.cap.sign = tab.deg.cap[(tab.deg.cap$WD_3m_Q.val<0.05),]
# rm NAs
tab.deg.cap.sign = tab.deg.cap.sign[complete.cases(tab.deg.cap.sign),]
nrow(tab.deg.cap.sign) # 274 genes with signif Q-val
head(tab.deg.cap.sign)
####

# keep only logFC
tab.deg.cap.logFC = tab.deg.cap.sign[,seq(1,10,2)] # tab.deg.cap[,seq(1,14,2)]
head(tab.deg.cap.logFC)

# replace NA with 0s
tab.deg.cap.logFC[is.na(tab.deg.cap.logFC)] = 0
head(tab.deg.cap.logFC)

# remove rows with all 0s
tab.deg.cap.logFC = tab.deg.cap.logFC[!rowSums(tab.deg.cap.logFC)==0,]
head(tab.deg.cap.logFC)
nrow(tab.deg.cap.logFC) # 4237 genes with logFC
colnames(tab.deg.cap.logFC) = c("WD_3m","WD_4m","WD_6m","Rev_1m","Rev_3m")
write.table(tab.deg.cap.logFC, "table_deg_liver_Cap_all_TP_sign_LOG.txt", sep = '\t')

b <- tab.deg.cap.logFC
# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]
# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5
mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma
ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
          main = "liver Cap DEG, logFC(x/chow)", 
          Rowv= T,Colv= T, 
          ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="both", cexRow=0.8, cexCol=0.8) 

# plot correlations

pdf(paste0("corr-cap-liver.pdf"),width=10,height=10,paper='special')
p1 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,WD_4m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_4m") + theme_linedraw() + geom_smooth(method = "lm")
p2 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,WD_6m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_6m") + theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,Rev_1m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_1m") + theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,Rev_3m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_31m") + theme_linedraw() + geom_smooth(method = "lm")

print((p1 | p2 ) /
        ( p3 |p4)) # 10 x 10
dev.off()


# all clusters
m <- cor(tab.deg.cap.logFC)
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black",
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))

############
setwd(paste0(file.dir,"/liver/deg-gen"))
# 2. Artery
tabs = list.files()[grep("liver_Art_", list.files())]
deg.artery = list()
deg.i = NULL
for (i in seq(1,4,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.artery[[i]] = deg.i
}

deg.artery[[5]] = read.table("/Users/olga.bondareva/Desktop/HI-MAG OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new/deg_gen/liver_EC-art_Western_vs_chow.txt", sep = '\t', header = T)

head(deg.artery[[1]])

# combine a table of logFC
tab.deg.art = merge(deg.artery[[5]][,c(2,5)], deg.artery[[3]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.art) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.art = merge(tab.deg.art, deg.artery[[4]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.art = merge(tab.deg.art, deg.artery[[1]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.art = merge(tab.deg.art, deg.artery[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.art)

rownames(tab.deg.art) = tab.deg.art$gene
tab.deg.art = tab.deg.art[,-1]
nrow(tab.deg.art) # 6848 genes
write.table(tab.deg.art,"table_deg_liver_art_all_TP.txt", sep='\t')

# signif DEGs
tab.deg.art.sign = tab.deg.art[(tab.deg.art$WD_3m_Q.val<0.05)|(tab.deg.art$WD_4m_Q.val<0.05)|(tab.deg.art$WD_6m_Q.val<0.05)|(tab.deg.art$Rev_1m_Q.val<0.05)|(tab.deg.art$Rev_3m_Q.val<0.05),]
# rm NAs
tab.deg.art.sign = tab.deg.art.sign[complete.cases(tab.deg.art.sign),]
nrow(tab.deg.art.sign) # 1671 genes with signif Q-val
head(tab.deg.art.sign)
write.table(tab.deg.art.sign, "table_deg_liver_art_all_TP_sign.txt", sep = '\t')

# keep only logFC
tab.deg.art.logFC = tab.deg.art.sign[,seq(1,10,2)] # tab.deg.art[,seq(1,14,2)]
head(tab.deg.art.logFC)

# replace NA with 0s
tab.deg.art.logFC[is.na(tab.deg.art.logFC)] = 0
head(tab.deg.art.logFC)

# remove rows with all 0s
tab.deg.art.logFC = tab.deg.art.logFC[!rowSums(tab.deg.art.logFC)==0,]
head(tab.deg.art.logFC)
nrow(tab.deg.art.logFC) # 591 genes with logFC
colnames(tab.deg.art.logFC) = c("WD_3m","WD_4m","WD_6m","Rev_1m","Rev_3m")
write.table(tab.deg.art.logFC, "table_deg_liver_art_all_TP_sign_LOG.txt", sep = '\t')

b <- tab.deg.art.logFC
# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]
# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5
mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma
ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))

pdf(paste0("heatmap-deg-art-sign-liver.pdf"),width=10,height=15,paper='special')
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
          main = "liver art DEG, logFC(x/chow)", 
          Rowv= T,Colv= T, 
          ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="both", cexRow=0.8, cexCol=0.8) 
dev.off()
# plot correlations

pdf(paste0("corr-art-liver.pdf"),width=10,height=10,paper='special')
p1 <- ggplot(tab.deg.art.logFC, aes(WD_3m,WD_4m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_4m") + theme_linedraw() + geom_smooth(method = "lm")
p2 <- ggplot(tab.deg.art.logFC, aes(WD_3m,WD_6m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_6m") + theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.art.logFC, aes(WD_3m,Rev_1m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_1m") + theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.art.logFC, aes(WD_3m,Rev_3m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_31m") + theme_linedraw() + geom_smooth(method = "lm")

print((p1 | p2 ) /
        ( p3 |p4)) # 10 x 10
dev.off()

# all clusters
m <- cor(tab.deg.art.logFC)
pdf(paste0("corr-art-tab-liver.pdf"),width=7,height=7,paper='special')
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black",
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))
dev.off()

############
# 3. Vein
tabs = list.files()[grep("liver_Ven_", list.files())]
deg.vein = list()
deg.i = NULL
for (i in seq(1,4,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.vein[[i]] = deg.i
}

deg.vein[[5]] = read.table("/Users/olga.bondareva/Desktop/HI-MAG OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new/deg_gen/liver_EC-ven_Western_vs_chow.txt", sep = '\t', header = T)

head(deg.vein[[1]])

# combine a table of logFC
tab.deg.ven = merge(deg.vein[[5]][,c(2,5)], deg.vein[[3]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.ven) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.ven = merge(tab.deg.ven, deg.vein[[4]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.ven = merge(tab.deg.ven, deg.vein[[1]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.ven = merge(tab.deg.ven, deg.vein[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.ven)

rownames(tab.deg.ven) = tab.deg.ven$gene
tab.deg.ven = tab.deg.ven[,-1]
nrow(tab.deg.ven) # 6848 genes
write.table(tab.deg.ven,"table_deg_liver_ven_all_TP.txt", sep='\t')

# WRONG calculation
# signif DEGs
tab.deg.ven.sign = tab.deg.ven[(tab.deg.ven$WD_3m_Q.val<0.05)|(tab.deg.ven$WD_4m_Q.val<0.05)|(tab.deg.ven$WD_6m_Q.val<0.05)|(tab.deg.ven$Rev_1m_Q.val<0.05)|(tab.deg.ven$Rev_3m_Q.val<0.05),]
# rm NAs
tab.deg.ven.sign = tab.deg.ven.sign[complete.cases(tab.deg.ven.sign),]
nrow(tab.deg.ven.sign) # 1671 genes with signif Q-val
head(tab.deg.ven.sign)
write.table(tab.deg.ven.sign, "table_deg_liver_ven_all_TP_sign.txt", sep = '\t')

# keep only logFC
tab.deg.ven.logFC = tab.deg.ven.sign[,seq(1,10,2)] # tab.deg.ven[,seq(1,14,2)]
head(tab.deg.ven.logFC)

# replace NA with 0s
tab.deg.ven.logFC[is.na(tab.deg.ven.logFC)] = 0
head(tab.deg.ven.logFC)

# remove rows with all 0s
tab.deg.ven.logFC = tab.deg.ven.logFC[!rowSums(tab.deg.ven.logFC)==0,]
head(tab.deg.ven.logFC)
nrow(tab.deg.ven.logFC) # 591 genes with logFC
colnames(tab.deg.ven.logFC) = c("WD_3m","WD_4m","WD_6m","Rev_1m","Rev_3m")
write.table(tab.deg.ven.logFC, "table_deg_liver_ven_all_TP_sign_LOG.txt", sep = '\t')

b <- tab.deg.ven.logFC
# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]
# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5
mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma
ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))

pdf(paste0("heatmap-deg-ven-sign-liver.pdf"),width=10,height=15,paper='special')
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
          main = "liver ven DEG, logFC(x/chow)", 
          Rowv= T,Colv= T, 
          ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="both", cexRow=0.8, cexCol=0.8) 
dev.off()
# plot correlations

pdf(paste0("corr-ven-liver.pdf"),width=10,height=10,paper='special')
p1 <- ggplot(tab.deg.ven.logFC, aes(WD_3m,WD_4m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_4m") + theme_linedraw() + geom_smooth(method = "lm")
p2 <- ggplot(tab.deg.ven.logFC, aes(WD_3m,WD_6m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_6m") + theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.ven.logFC, aes(WD_3m,Rev_1m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_1m") + theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.ven.logFC, aes(WD_3m,Rev_3m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_31m") + theme_linedraw() + geom_smooth(method = "lm")

print((p1 | p2 ) /
        ( p3 |p4)) # 10 x 10
dev.off()

# all clusters
m <- cor(tab.deg.ven.logFC)
pdf(paste0("corr-ven-tab-liver.pdf"),width=7,height=7,paper='special')
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black",
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))
dev.off()



saveRDS(liver, "liver_3m_4m_6m.rds")

liver = readRDS("liver_3m_4m_6m.rds")

# Artery markers diff organs
FeaturePlot(liver, features = c("Fbln5", "Cytl1", "Gkn3", "Vegfc",# brain
                                "Stmn2", "Mgp", "Cldn5", "Eln", # kidney
                                "Scin", "Igfbp7", "Slc14a1","Hpgd", # kidney 2
                                "Atp13a3", "Adgrg6", "Plac8", "Ltbp4", # lung
                               "Clu", "Plac8", "Lrg1", "Efnb1", "Gja4", "Sox17", "Glul", # liver, Kalucka paper
                               "Ednrb", "Igfbp5", "Atp13a3", "Adgrg6", "Ly6a"), # liver artery??
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 4) #

# Vein markers diff organs
FeaturePlot(liver, features = c("Cfh","Scl38a5", "Nr2f2","Tmsb10",# brain
                                "Lbp", "Car8", "Serpinb6b",# kidney , Vcam1
                                "Igf1", "Cryab", "Fxyd6","Cd36","Bgn", # kidney 2
                                "Selp", "Vcam1", "Slc6a2", "Prss23", # lung , Vwf
                                "Thbd", "Edn1", "Vwf", "Bgn", # liver, Kalucka paper
                                "Rpos3", "Bmp4", "Fabp4", "Fbln2"), # liver vein?? , Vwf
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 4) #



#########################
# 4. heart
#########################

setwd("../heart")
heart.3m = readRDS(paste0(data.dir.3m,"/heart/heart_ECs_singlet.rds"))
heart.4m = readRDS(paste0(data.dir.4m,"heart_subset.rds"))
heart.6m = readRDS(paste0(data.dir.6m,"heart_subset.rds"))
heart.3m$time = "3m"

# merge objects
heart = merge(heart.3m, list(heart.4m, heart.6m))
rm(heart.3m, heart.4m,heart.6m)

heart$diet_time = paste0(heart$diet, heart$time)

# Pericyte amd SMC markers
FeaturePlot(heart, features = c("Flt1", "Pecam1","Cdh5",# EC
                                "Pdgfrb", "Cspg4","Anpep",  #  Pericyte
                                "Myh11", "Acta2", "Des",   # SMC
                                "Dcn", "Col1a1", "Pdgfra"), #FB
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10") #, min.cutoff = "q10"

FeaturePlot(heart, features = c("Prox1", "Pdpn","Lyve1", "Mki67", "Ptprc"), # lymph, prolif, immune
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
DimPlot(object = heart, reduction = "umap", group.by = "diet_time")
DimPlot(object = heart, reduction = "umap", group.by = "heart_celltype",label = TRUE)
table(Idents(heart), heart$orig.ident)

# select cluster by hand
#plot <- DimPlot(heart, reduction = "umap")
#select.cells <- CellSelector(plot = plot) # 
#Idents(heart, cells = select.cells) <- "11"


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

FeaturePlot(heart, features = c("Vwf", "Vcam1", #  large vessels
                                "Fbln5", "Cytl1", # large artery
                                "Hey1", "Gkn3"# artery 
),
order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 3,raster = T) #

#  tip cell-enriched, angiogenic
FeaturePlot(heart, features = c("Apln", "Col4a2", "Trp53i11"),
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 3)

DimPlot(object = heart, reduction = "umap", label = TRUE) + NoLegend()
heart <- RenameIdents(heart, `8` = "EC-ven",`0` = "EC-cap1",`2` = "EC-cap2",`3` = "EC-cap3",`6` = "EC-cap4",
                      `4` = "EC-art",`5` = "EC-venule", `10` = "EC-Hb", `12` = "Immune",
                      `9` = "EC-lymph", `11` = "Prolif",`1` = "EC-AP1",`7` = "EC-IFN") 
## Relevel object@ident
heart@active.ident <- factor(x = heart@active.ident, levels = sort(levels(Idents(heart))))
heart$heart_celltype <- Idents(heart)
DimPlot(heart, reduction = "umap", label = TRUE) + NoLegend()

# number of cells per cluster per sample
table(Idents(heart), heart$diet_time)

heart = subset(heart, idents = "Immune", invert = T)

DimPlot(heart, reduction = "umap", label = TRUE, split.by = "diet_time", ncol = 3, pt.size = 1.3) + NoLegend()


heart = readRDS("heart_3m_4m_6m.rds")

#####################
# subselect arterial cells tip
# select cluster by hand
plot <- DimPlot(heart, reduction = "umap")
select.cells <- CellSelector(plot = plot) # 
Idents(heart, cells = select.cells) <- "art-big"

## Relevel object@ident
heart@active.ident <- factor(x = heart@active.ident, levels = sort(levels(Idents(heart))))
heart$heart_celltype3 <- Idents(heart)
DimPlot(heart, reduction = "umap", label = TRUE)
table(heart$diet_time,heart$heart_celltype3)


# DE genes
heart$clust.diet_time <- paste(heart$heart_celltype3, heart$diet_time, sep = "_")
Idents(heart) <- "clust.diet_time"
heart@active.ident <- factor(x = heart@active.ident, levels = sort(levels(Idents(object = heart))))
levels(Idents(object = heart))

lev1 <- levels(heart@active.ident)[seq(7,(length(levels(Idents(heart)))),8)] # condition 1=HFD_4m
lev2 <- levels(heart@active.ident)[seq(2,(length(levels(Idents(heart)))),8)] # condition 2=chow_4m
lev3 <- levels(heart@active.ident)[seq(4,(length(levels(Idents(heart)))),8)] # condition 3=Rev_1m
lev4 <- levels(heart@active.ident)[seq(8,(length(levels(Idents(heart)))),8)] # condition 4=HFD_6m
lev5 <- levels(heart@active.ident)[seq(3,(length(levels(Idents(heart)))),8)] # condition 5=chow_6m
lev6 <- levels(heart@active.ident)[seq(5,(length(levels(Idents(heart)))),8)] # condition 6=Rev_3m
lev7 <- levels(heart@active.ident)[seq(6,(length(levels(Idents(heart)))),8)] # condition 7=HFD_3m
lev8 <- levels(heart@active.ident)[seq(1,(length(levels(Idents(heart)))),8)] # condition 8=chow_3m

# for all the clusters
for (i in seq(1,length(lev1),1)) {
  # 4 months
  # Western
  print(lev1[i])
  print(lev2[i])
  response2 <- FindMarkers(object = heart, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # Reverse
  print(lev3[i])
  print(lev2[i])
  response2 <- FindMarkers(object = heart, ident.1 = lev3[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev3[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # 6 months
  # Western
  print(lev4[i])
  print(lev5[i])
  response2 <- FindMarkers(object = heart, ident.1 = lev4[i], 
                           ident.2 = lev5[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev4[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # Reverse
  print(lev6[i])
  print(lev5[i])
  response2 <- FindMarkers(object = heart, ident.1 = lev6[i], 
                           ident.2 = lev5[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev6[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # 3 months
  # Western
  print(lev7[i])
  print(lev8[i])
  response2 <- FindMarkers(object = heart, ident.1 = lev7[i], 
                           ident.2 = lev8[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev7[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
}

### big artery
setwd("deg-clust")
tabs = list.files()[grep("art-big_", list.files())]
deg.capillaries = list()
deg.i = NULL
for (i in seq(1,5,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.capillaries[[i]] = deg.i
}

head(deg.capillaries[[1]])

# combine a table of logFC
tab.deg.cap = merge(deg.capillaries[[3]][,c(2,5)], deg.capillaries[[4]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.cap) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[5]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[1]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.cap)

rownames(tab.deg.cap) = tab.deg.cap$gene
tab.deg.cap = tab.deg.cap[,-1]
nrow(tab.deg.cap) # 8950 genes
write.table(tab.deg.cap,"table_deg_heart_art-big_all_TP.txt", sep='\t')


# DEG for capillaries combined
# 4m
response2 <- FindMarkers(object = heart, 
                         ident.1 = c("EC-AP1_Western4m","EC-cap1_Western4m","EC-cap2_Western4m",
                                     "EC-cap3_Western4m","EC-cap4_Western4m","EC-IFN_Western4m"), 
                         ident.2 = c("EC-AP1_chow4m","EC-cap1_chow4m","EC-cap2_chow4m",
                                     "EC-cap3_chow4m","EC-cap4_chow4m","EC-IFN_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "heart_Cap_WB_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = heart, 
                         ident.1 = c("EC-AP1_Rev1m","EC-cap1_Rev1m","EC-cap2_Rev1m",
                                     "EC-cap3_Rev1m","EC-cap4_Rev1m","EC-IFN_Rev1m"), 
                         ident.2 = c("EC-AP1_chow4m","EC-cap1_chow4m","EC-cap2_chow4m",
                                     "EC-cap3_chow4m","EC-cap4_chow4m","EC-IFN_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "heart_Cap_Rev1m_vs_chow_4m.txt", sep = "\t")

# 6m
response2 <- FindMarkers(object = heart, 
                         ident.1 = c("EC-AP1_Western6m","EC-cap1_Western6m","EC-cap2_Western6m",
                                     "EC-cap3_Western6m","EC-cap4_Western6m","EC-IFN_Western6m"), 
                         ident.2 = c("EC-AP1_chow6m","EC-cap1_chow6m","EC-cap2_chow6m",
                                     "EC-cap3_chow6m","EC-cap4_chow6m","EC-IFN_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "heart_Cap_WB_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = heart, 
                         ident.1 = c("EC-AP1_Rev3m","EC-cap1_Rev3m","EC-cap2_Rev3m",
                                     "EC-cap3_Rev3m","EC-cap4_Rev3m","EC-IFN_Rev3m"), 
                         ident.2 = c("EC-AP1_chow6m","EC-cap1_chow6m","EC-cap2_chow6m",
                                     "EC-cap3_chow6m","EC-cap4_chow6m","EC-IFN_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "heart_Cap_Rev3m_vs_chow_6m.txt", sep = "\t")

# DEG for vein and artery combined
# 4m
response2 <- FindMarkers(object = heart, 
                         ident.1 = c("EC-art_Western4m"), 
                         ident.2 = c("EC-art_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "heart_Art_WB_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = heart, 
                         ident.1 = c("EC-ven_Western4m","EC-venule_Western4m"), 
                         ident.2 = c("EC-ven_chow4m","EC-venule_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "heart_Ven_WB_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = heart, 
                         ident.1 = c("EC-art_Rev1m"), 
                         ident.2 = c("EC-art_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "heart_Art_Rev1m_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = heart, 
                         ident.1 = c("EC-ven_Rev1m","EC-venule_Rev1m"), 
                         ident.2 = c("EC-ven_chow4m","EC-venule_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "heart_Ven_Rev1m_vs_chow_4m.txt", sep = "\t")

# 6m
response2 <- FindMarkers(object = heart, 
                         ident.1 = c("EC-art_Western6m"), 
                         ident.2 = c("EC-art_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "heart_Art_WB_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = heart, 
                         ident.1 = c("EC-ven_Western6m","EC-venule_Western6m"), 
                         ident.2 = c("EC-ven_chow6m","EC-venule_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "heart_Ven_WB_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = heart, 
                         ident.1 = c("EC-art_Rev3m"), 
                         ident.2 = c("EC-art_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "heart_Art_Rev3m_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = heart, 
                         ident.1 = c("EC-ven_Rev3m","EC-venule_Rev3m"), 
                         ident.2 = c("EC-ven_chow6m","EC-venule_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "heart_Ven_Rev3m_vs_chow_6m.txt", sep = "\t")

####
# combine Cap/ Art/ Ven in 1 table

# 1. Capillaries
tabs = list.files()[grep("heart_Cap_", list.files())]
deg.capillaries = list()
deg.i = NULL
for (i in seq(1,4,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.capillaries[[i]] = deg.i
}

deg.capillaries[[5]] = read.table("/Users/olga.bondareva/Desktop/HI-MAG OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new/deg_gen/heart_EC-cap_Western_vs_chow.txt", sep = '\t', header = T)

head(deg.capillaries[[1]])

# combine a table of logFC
tab.deg.cap = merge(deg.capillaries[[5]][,c(2,5)], deg.capillaries[[3]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.cap) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[4]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[1]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.cap)

rownames(tab.deg.cap) = tab.deg.cap$gene
tab.deg.cap = tab.deg.cap[,-1]
nrow(tab.deg.cap) # 7842 genes
write.table(tab.deg.cap,"table_deg_heart_Cap_all_TP.txt", sep='\t')

# signif DEGs
tab.deg.cap.sign = tab.deg.cap[(tab.deg.cap$WD_3m_Q.val<0.05)|(tab.deg.cap$WD_4m_Q.val<0.05)|(tab.deg.cap$WD_6m_Q.val<0.05)|(tab.deg.cap$Rev_1m_Q.val<0.05)|(tab.deg.cap$Rev_3m_Q.val<0.05),]
# rm NAs
tab.deg.cap.sign = tab.deg.cap.sign[complete.cases(tab.deg.cap.sign),]
nrow(tab.deg.cap.sign) # 1657 genes with signif Q-val
head(tab.deg.cap.sign)
write.table(tab.deg.cap.sign, "table_deg_heart_Cap_all_TP_sign.txt", sep = '\t')

######
# signif only at 3m TP
tab.deg.cap.sign = tab.deg.cap[(tab.deg.cap$WD_3m_Q.val<0.05),]
# rm NAs
tab.deg.cap.sign = tab.deg.cap.sign[complete.cases(tab.deg.cap.sign),]
nrow(tab.deg.cap.sign) # 274 genes with signif Q-val
head(tab.deg.cap.sign)
######


# keep only logFC
tab.deg.cap.logFC = tab.deg.cap.sign[,seq(1,10,2)] # tab.deg.cap[,seq(1,14,2)]
head(tab.deg.cap.logFC)

# replace NA with 0s
tab.deg.cap.logFC[is.na(tab.deg.cap.logFC)] = 0
head(tab.deg.cap.logFC)

# remove rows with all 0s
tab.deg.cap.logFC = tab.deg.cap.logFC[!rowSums(tab.deg.cap.logFC)==0,]
head(tab.deg.cap.logFC)
nrow(tab.deg.cap.logFC) # 2090 genes with logFC
colnames(tab.deg.cap.logFC) = c("WD_3m","WD_4m","WD_6m","Rev_1m","Rev_3m")
write.table(tab.deg.cap.logFC, "table_deg_heart_Cap_all_TP_sign_LOG.txt", sep = '\t')

b <- tab.deg.cap.logFC
# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]
# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5
mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma
ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
          main = "heart Cap DEG, logFC(x/chow)", 
          Rowv= T,Colv= T, 
          ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="both", cexRow=0.8, cexCol=0.8) 

# plot correlations

pdf(paste0("corr-cap-heart.pdf"),width=10,height=10,paper='special')
p1 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,WD_4m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_4m") + theme_linedraw() + geom_smooth(method = "lm")
p2 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,WD_6m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_6m") + theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,Rev_1m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_1m") + theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,Rev_3m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_31m") + theme_linedraw() + geom_smooth(method = "lm")

print((p1 | p2 ) /
        ( p3 |p4)) # 10 x 10
dev.off()

m <- cor(tab.deg.cap.logFC)
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black",
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))

############
setwd(paste0(file.dir,"/heart/deg-gen"))
# 2. Artery
tabs = list.files()[grep("heart_Art_", list.files())]
deg.artery = list()
deg.i = NULL
for (i in seq(1,4,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.artery[[i]] = deg.i
}

deg.artery[[5]] = read.table("/Users/olga.bondareva/Desktop/HI-MAG OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new/deg_gen/heart_EC-art_Western_vs_chow.txt", sep = '\t', header = T)

head(deg.artery[[1]])

# combine a table of logFC
tab.deg.art = merge(deg.artery[[5]][,c(2,5)], deg.artery[[3]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.art) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.art = merge(tab.deg.art, deg.artery[[4]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.art = merge(tab.deg.art, deg.artery[[1]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.art = merge(tab.deg.art, deg.artery[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.art)

rownames(tab.deg.art) = tab.deg.art$gene
tab.deg.art = tab.deg.art[,-1]
nrow(tab.deg.art) # 6848 genes
write.table(tab.deg.art,"table_deg_heart_art_all_TP.txt", sep='\t')

# signif DEGs
tab.deg.art.sign = tab.deg.art[(tab.deg.art$WD_3m_Q.val<0.05)|(tab.deg.art$WD_4m_Q.val<0.05)|(tab.deg.art$WD_6m_Q.val<0.05)|(tab.deg.art$Rev_1m_Q.val<0.05)|(tab.deg.art$Rev_3m_Q.val<0.05),]
# rm NAs
tab.deg.art.sign = tab.deg.art.sign[complete.cases(tab.deg.art.sign),]
nrow(tab.deg.art.sign) # 1671 genes with signif Q-val
head(tab.deg.art.sign)
write.table(tab.deg.art.sign, "table_deg_heart_art_all_TP_sign.txt", sep = '\t')

# keep only logFC
tab.deg.art.logFC = tab.deg.art.sign[,seq(1,10,2)] # tab.deg.art[,seq(1,14,2)]
head(tab.deg.art.logFC)

# replace NA with 0s
tab.deg.art.logFC[is.na(tab.deg.art.logFC)] = 0
head(tab.deg.art.logFC)

# remove rows with all 0s
tab.deg.art.logFC = tab.deg.art.logFC[!rowSums(tab.deg.art.logFC)==0,]
head(tab.deg.art.logFC)
nrow(tab.deg.art.logFC) # 591 genes with logFC
colnames(tab.deg.art.logFC) = c("WD_3m","WD_4m","WD_6m","Rev_1m","Rev_3m")
write.table(tab.deg.art.logFC, "table_deg_heart_art_all_TP_sign_LOG.txt", sep = '\t')

b <- tab.deg.art.logFC
# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]
# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5
mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma
ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))

pdf(paste0("heatmap-deg-art-sign-heart.pdf"),width=10,height=15,paper='special')
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
          main = "heart art DEG, logFC(x/chow)", 
          Rowv= T,Colv= T, 
          ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="both", cexRow=0.8, cexCol=0.8) 
dev.off()
# plot correlations

pdf(paste0("corr-art-heart.pdf"),width=10,height=10,paper='special')
p1 <- ggplot(tab.deg.art.logFC, aes(WD_3m,WD_4m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_4m") + theme_linedraw() + geom_smooth(method = "lm")
p2 <- ggplot(tab.deg.art.logFC, aes(WD_3m,WD_6m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_6m") + theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.art.logFC, aes(WD_3m,Rev_1m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_1m") + theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.art.logFC, aes(WD_3m,Rev_3m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_31m") + theme_linedraw() + geom_smooth(method = "lm")

print((p1 | p2 ) /
        ( p3 |p4)) # 10 x 10
dev.off()

# all clusters
m <- cor(tab.deg.art.logFC)
pdf(paste0("corr-art-tab-heart.pdf"),width=7,height=7,paper='special')
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black",
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))
dev.off()

############
# 3. Vein
tabs = list.files()[grep("heart_Ven_", list.files())]
deg.vein = list()
deg.i = NULL
for (i in seq(1,4,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.vein[[i]] = deg.i
}

deg.vein[[5]] = read.table("/Users/olga.bondareva/Desktop/HI-MAG OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new/deg_gen/heart_EC-ven_Western_vs_chow.txt", sep = '\t', header = T)

head(deg.vein[[1]])

# combine a table of logFC
tab.deg.ven = merge(deg.vein[[5]][,c(2,5)], deg.vein[[3]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.ven) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.ven = merge(tab.deg.ven, deg.vein[[4]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.ven = merge(tab.deg.ven, deg.vein[[1]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.ven = merge(tab.deg.ven, deg.vein[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.ven)

rownames(tab.deg.ven) = tab.deg.ven$gene
tab.deg.ven = tab.deg.ven[,-1]
nrow(tab.deg.ven) # 6848 genes
write.table(tab.deg.ven,"table_deg_heart_ven_all_TP.txt", sep='\t')

# signif DEGs
tab.deg.ven.sign = tab.deg.ven[(tab.deg.ven$WD_3m_Q.val<0.05)|(tab.deg.ven$WD_4m_Q.val<0.05)|(tab.deg.ven$WD_6m_Q.val<0.05)|(tab.deg.ven$Rev_1m_Q.val<0.05)|(tab.deg.ven$Rev_3m_Q.val<0.05),]
# rm NAs
tab.deg.ven.sign = tab.deg.ven.sign[complete.cases(tab.deg.ven.sign),]
nrow(tab.deg.ven.sign) # 1671 genes with signif Q-val
head(tab.deg.ven.sign)
write.table(tab.deg.ven.sign, "table_deg_heart_ven_all_TP_sign.txt", sep = '\t')

# keep only logFC
tab.deg.ven.logFC = tab.deg.ven.sign[,seq(1,10,2)] # tab.deg.ven[,seq(1,14,2)]
head(tab.deg.ven.logFC)

# replace NA with 0s
tab.deg.ven.logFC[is.na(tab.deg.ven.logFC)] = 0
head(tab.deg.ven.logFC)

# remove rows with all 0s
tab.deg.ven.logFC = tab.deg.ven.logFC[!rowSums(tab.deg.ven.logFC)==0,]
head(tab.deg.ven.logFC)
nrow(tab.deg.ven.logFC) # 591 genes with logFC
colnames(tab.deg.ven.logFC) = c("WD_3m","WD_4m","WD_6m","Rev_1m","Rev_3m")
write.table(tab.deg.ven.logFC, "table_deg_heart_ven_all_TP_sign_LOG.txt", sep = '\t')

b <- tab.deg.ven.logFC
# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]
# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5
mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma
ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))

pdf(paste0("heatmap-deg-ven-sign-heart.pdf"),width=10,height=15,paper='special')
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
          main = "heart ven DEG, logFC(x/chow)", 
          Rowv= T,Colv= T, 
          ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="both", cexRow=0.8, cexCol=0.8) 
dev.off()
# plot correlations

pdf(paste0("corr-ven-heart.pdf"),width=10,height=10,paper='special')
p1 <- ggplot(tab.deg.ven.logFC, aes(WD_3m,WD_4m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_4m") + theme_linedraw() + geom_smooth(method = "lm")
p2 <- ggplot(tab.deg.ven.logFC, aes(WD_3m,WD_6m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_6m") + theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.ven.logFC, aes(WD_3m,Rev_1m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_1m") + theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.ven.logFC, aes(WD_3m,Rev_3m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_31m") + theme_linedraw() + geom_smooth(method = "lm")

print((p1 | p2 ) /
        ( p3 |p4)) # 10 x 10
dev.off()

# all clusters
m <- cor(tab.deg.ven.logFC)
pdf(paste0("corr-ven-tab-heart.pdf"),width=7,height=7,paper='special')
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black",
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))
dev.off()



saveRDS(heart, "heart_3m_4m_6m.rds")

### arterial only (w/o arterioles) EC
setwd("../heart/deg-clust")
tabs = list.files()[grep("EC-art_", list.files())]
deg.capillaries = list()
deg.i = NULL
for (i in seq(1,5,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.capillaries[[i]] = deg.i
}

head(deg.capillaries[[1]])

# combine a table of logFC
tab.deg.cap = merge(deg.capillaries[[3]][,c(2,5)], deg.capillaries[[4]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.cap) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[5]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[1]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.cap)

rownames(tab.deg.cap) = tab.deg.cap$gene
tab.deg.cap = tab.deg.cap[,-1]
nrow(tab.deg.cap) # 8950 genes
write.table(tab.deg.cap,"table_deg_brain_fenestr_all_TP.txt", sep='\t')

#########################
# 5. kidney
#########################

setwd("kidney")
kidney.3m = readRDS(paste0(data.dir.3m,"/kidney/kidney_ECs_singlet.rds"))
kidney.4m = readRDS(paste0(data.dir.4m,"kidney_subset.rds"))
kidney.6m = readRDS(paste0(data.dir.6m,"kidney_subset.rds"))

kidney.3m$time = "3m"

# merge objects
kidney = merge(kidney.3m, list(kidney.4m,kidney.6m))
rm(kidney.3m, kidney.4m,kidney.6m)

kidney$diet_time = paste0(kidney$diet, kidney$time)

kidney <- ScaleData(object = kidney, verbose = T)
kidney <- FindVariableFeatures(kidney, selection.method = "vst", nfeatures = 3000)
kidney <- RunPCA(object = kidney, npcs = 30, verbose = T)
ElbowPlot(object = kidney,  ndims = 30) 

# UMAP and Clustering
kidney <- RunUMAP(object = kidney, reduction = "pca", dims = 1:30)
kidney <- FindNeighbors(object = kidney, reduction = "pca", dims = 1:30)
kidney <- FindClusters(object = kidney, resolution = 0.4)

# visualization
DimPlot(object = kidney, reduction = "umap", label = TRUE)
DimPlot(object = kidney, reduction = "umap", group.by = "diet_time")
DimPlot(object = kidney, reduction = "umap", group.by = "kidney_celltype", label = T)


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
FeaturePlot(kidney, features = c("Apln", "Col4a2", "Trp53i11", "Mki67", "Cenpf"),
            order = T, cols=c("grey", "red"), min.cutoff = "q10", ncol = 3)



# Pericyte amd SMC markers
FeaturePlot(kidney, features = c("Flt1", "Pecam1","Cdh5",# EC
                                 "Pdgfrb", "Cspg4","Anpep",  #  Pericyte
                                 "Myh11", "Acta2", "Des",   # SMC
                                 "Dcn", "Col1a1", "Pdgfra"), #FB
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10") #, min.cutoff = "q10"

FeaturePlot(kidney, features = c("Prox1", "Pdpn","Lyve1"), # lymph
            order = T, cols=c("grey", "red"), ncol=2, min.cutoff = "q10") #, min.cutoff = "q10"

# select cluster by hand
Idents(kidney) = kidney$seurat_clusters
plot <- DimPlot(kidney, reduction = "umap", cells = WhichCells(kidney, idents = "4"))
select.cells <- CellSelector(plot = plot) # 
Idents(kidney, cells = select.cells) <- "EC-ang"


DimPlot(object = kidney, reduction = "umap", label = TRUE) + NoLegend()

kidney <- RenameIdents(kidney, `0` = "EC-cap1",`3` = "EC-cap2",`4` = "EC-cap3",
                       `5` = "EC-cap4",`10` = "EC-cap3",`12` = "EC-cap5",
                       `8` = "EC-art",`6` = "EC-ven",`14` = "proximal tubule",
                      `9` = "EC-AP1",`7` = "EC-kidney1", `13` = "EC-kidney2", `4` = "EC-ang", 
                      `11` = "EC-glomeruli", `1` = "EC-Aqp1-arteriole", `15` = "EC-lymph", 
                      `2` = "EC-venule",`16` = "Prolif") 


## Relevel object@ident
kidney@active.ident <- factor(x = kidney@active.ident, levels = sort(levels(Idents(kidney))))
kidney$kidney_celltype2 <- Idents(kidney)
DimPlot(kidney, reduction = "umap", label = TRUE) + NoLegend()

kidney = subset(kidney, idents = "proximal tubule", invert = T)

# UMAP and Clustering
kidney <- RunUMAP(object = kidney, reduction = "pca", dims = 1:30)
kidney <- FindNeighbors(object = kidney, reduction = "pca", dims = 1:30)
kidney <- FindClusters(object = kidney, resolution = 0.4)

DimPlot(object = kidney, reduction = "umap", label = T)
DimPlot(object = kidney, reduction = "umap", group.by = "kidney_celltype", label = T)

# select cluster by hand
Idents(kidney) = kidney$seurat_clusters

plot <- DimPlot(kidney, reduction = "umap", cells = WhichCells(kidney, idents = "4"))
select.cells <- CellSelector(plot = plot) # 
Idents(kidney, cells = select.cells) <- "EC-ang"

plot <- DimPlot(kidney, reduction = "umap")
select.cells <- CellSelector(plot = plot) # 
Idents(kidney, cells = select.cells) <- "Prolif"


kidney <- RenameIdents(kidney, `0` = "EC-cap1",`3` = "EC-cap2",`4` = "EC-cap3",
                       `11` = "EC-cap4",`8` = "EC-cap3",
                       `7` = "EC-art",`5` = "EC-ven",
                       `9` = "EC-AP1",`6` = "EC-kidney1", `12` = "EC-kidney2", # `4` = "EC-ang", 
                       `10` = "EC-glomeruli", `1` = "EC-Aqp1-arteriole", `13` = "EC-lymph", 
                       `2` = "EC-venule") 

kidney <- RenameIdents(kidney, `EC-kidney1` = "mEC1", `EC-kidney2` = "mEC2", # `4` = "EC-ang", 
                       `EC-glomeruli` = "glomerular-EC") 

## Relevel object@ident
kidney@active.ident <- factor(x = kidney@active.ident, levels = sort(levels(Idents(kidney))))
kidney$kidney_celltype2 <- Idents(kidney)
DimPlot(kidney, reduction = "umap", label = TRUE) + NoLegend()

# number of cells per cluster per sample
table(Idents(kidney), kidney$diet_time)

DimPlot(object = kidney, reduction = "umap", group.by = "kidney_celltype2", label = T)
DimPlot(kidney, reduction = "umap", label = TRUE, split.by = "diet_time", ncol = 3, pt.size = 1.3) + NoLegend()


kidney = readRDS("kidney_3m_4m_6m.rds")

# DE genes
kidney$clust.diet_time <- paste(kidney$kidney_celltype2, kidney$diet_time, sep = "_")
Idents(kidney) <- "clust.diet_time"
kidney@active.ident <- factor(x = kidney@active.ident, levels = sort(levels(Idents(kidney))))
levels(Idents( kidney))

lev1 <- levels(kidney@active.ident)[seq(7,(length(levels(Idents(kidney)))),8)] # condition 1=HFD_4m
lev2 <- levels(kidney@active.ident)[seq(2,(length(levels(Idents(kidney)))),8)] # condition 2=chow_4m
lev3 <- levels(kidney@active.ident)[seq(4,(length(levels(Idents(kidney)))),8)] # condition 3=Rev_1m
lev4 <- levels(kidney@active.ident)[seq(8,(length(levels(Idents(kidney)))),8)] # condition 4=HFD_6m
lev5 <- levels(kidney@active.ident)[seq(3,(length(levels(Idents(kidney)))),8)] # condition 5=chow_6m
lev6 <- levels(kidney@active.ident)[seq(5,(length(levels(Idents(kidney)))),8)] # condition 6=Rev_3m
lev7 <- levels(kidney@active.ident)[seq(6,(length(levels(Idents(kidney)))),8)] # condition 7=HFD_3m
lev8 <- levels(kidney@active.ident)[seq(1,(length(levels(Idents(kidney)))),8)] # condition 8=chow_3m

# for all the clusters
for (i in seq(1,length(lev1),1)) {
  # 4 months
  # Western
  print(lev1[i])
  print(lev2[i])
  response2 <- FindMarkers(object = kidney, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # Reverse
  print(lev3[i])
  print(lev2[i])
  response2 <- FindMarkers(object = kidney, ident.1 = lev3[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev3[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # 6 months
  # Western
  print(lev4[i])
  print(lev5[i])
  response2 <- FindMarkers(object = kidney, ident.1 = lev4[i], 
                           ident.2 = lev5[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev4[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # Reverse
  print(lev6[i])
  print(lev5[i])
  response2 <- FindMarkers(object = kidney, ident.1 = lev6[i], 
                           ident.2 = lev5[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev6[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # 3 months
  # Western
  print(lev7[i])
  print(lev8[i])
  response2 <- FindMarkers(object = kidney, ident.1 = lev7[i], 
                           ident.2 = lev8[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev7[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
}



VlnPlot(kidney, group.by = "kidney_celltype2", split.by = "diet_time", 
        features = c("Jun","Junb","Jund","Fos","Fosb","Egr1"), ncol = 1, pt.size = 0)

DotPlot(kidney, group.by = "kidney_celltype2", split.by = "diet_time", 
        features = c("Jun","Junb","Jund","Fos","Fosb","Egr1"), cols = c("red","red","red","red","red","red","red","red"))

# DEG for capillaries combined
# 4m
response2 <- FindMarkers(object = kidney, 
                         ident.1 = c("EC-AP1_Western4m","EC-cap1_Western4m","EC-cap2_Western4m",
                                     "EC-cap3_Western4m","EC-cap4_Western4m"), 
                         ident.2 = c("EC-AP1_chow4m","EC-cap1_chow4m","EC-cap2_chow4m",
                                     "EC-cap3_chow4m","EC-cap4_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "kidney_Cap_WB_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = kidney, 
                         ident.1 = c("EC-AP1_Rev1m","EC-cap1_Rev1m","EC-cap2_Rev1m",
                                     "EC-cap3_Rev1m","EC-cap4_Rev1m"), 
                         ident.2 = c("EC-AP1_chow4m","EC-cap1_chow4m","EC-cap2_chow4m",
                                     "EC-cap3_chow4m","EC-cap4_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "kidney_Cap_Rev1m_vs_chow_4m.txt", sep = "\t")

# 6m
response2 <- FindMarkers(object = kidney, 
                         ident.1 = c("EC-AP1_Western6m","EC-cap1_Western6m","EC-cap2_Western6m",
                                     "EC-cap3_Western6m","EC-cap4_Western6m"), 
                         ident.2 = c("EC-AP1_chow6m","EC-cap1_chow6m","EC-cap2_chow6m",
                                     "EC-cap3_chow6m","EC-cap4_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "kidney_Cap_WB_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = kidney, 
                         ident.1 = c("EC-AP1_Rev3m","EC-cap1_Rev3m","EC-cap2_Rev3m",
                                     "EC-cap3_Rev3m","EC-cap4_Rev3m"), 
                         ident.2 = c("EC-AP1_chow6m","EC-cap1_chow6m","EC-cap2_chow6m",
                                     "EC-cap3_chow6m","EC-cap4_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "kidney_Cap_Rev3m_vs_chow_6m.txt", sep = "\t")

# DEG for vein and artery combined
# 4m
response2 <- FindMarkers(object = kidney, 
                         ident.1 = c("EC-art_Western4m","EC-Aqp1-arteriole_Western4m"), 
                         ident.2 = c("EC-art_chow4m","EC-Aqp1-arteriole_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "kidney_Art_WB_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = kidney, 
                         ident.1 = c("EC-ven_Western4m","EC-venule_Western4m"), 
                         ident.2 = c("EC-ven_chow4m","EC-venule_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "kidney_Ven_WB_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = kidney, 
                         ident.1 = c("EC-art_Rev1m","EC-Aqp1-arteriole_Rev1m"), 
                         ident.2 = c("EC-art_chow4m","EC-Aqp1-arteriole_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "kidney_Art_Rev1m_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = kidney, 
                         ident.1 = c("EC-ven_Rev1m","EC-venule_Rev1m"), 
                         ident.2 = c("EC-ven_chow4m","EC-venule_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "kidney_Ven_Rev1m_vs_chow_4m.txt", sep = "\t")

# 6m
response2 <- FindMarkers(object = kidney, 
                         ident.1 = c("EC-art_Western6m","EC-Aqp1-arteriole_Western6m"), 
                         ident.2 = c("EC-art_chow6m","EC-Aqp1-arteriole_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "kidney_Art_WB_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = kidney, 
                         ident.1 = c("EC-ven_Western6m","EC-venule_Western6m"), 
                         ident.2 = c("EC-ven_chow6m","EC-venule_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "kidney_Ven_WB_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = kidney, 
                         ident.1 = c("EC-art_Rev3m","EC-Aqp1-arteriole_Rev3m"), 
                         ident.2 = c("EC-art_chow6m","EC-Aqp1-arteriole_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "kidney_Art_Rev3m_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = kidney, 
                         ident.1 = c("EC-ven_Rev3m","EC-venule_Rev3m"), 
                         ident.2 = c("EC-ven_chow6m","EC-venule_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "kidney_Ven_Rev3m_vs_chow_6m.txt", sep = "\t")

# 1. Kidney1 Capillaries
setwd("kidney/deg-clust")
tabs = list.files()[grep("EC-kidney1_", list.files())]
deg.capillaries = list()
deg.i = NULL
for (i in seq(1,5,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.capillaries[[i]] = deg.i
}

head(deg.capillaries[[1]])

# combine a table of logFC
tab.deg.cap = merge(deg.capillaries[[3]][,c(2,5)], deg.capillaries[[4]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.cap) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[5]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[1]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.cap)

rownames(tab.deg.cap) = tab.deg.cap$gene
tab.deg.cap = tab.deg.cap[,-1]
nrow(tab.deg.cap) # 8950 genes
write.table(tab.deg.cap,"table_deg_kidney_kidney1_all_TP.txt", sep='\t')

# 2. Kidney2 Capillaries
setwd("../organs/kidney/deg-clust")
tabs = list.files()[grep("EC-kidney2_", list.files())]
deg.capillaries = list()
deg.i = NULL
for (i in seq(1,5,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.capillaries[[i]] = deg.i
}

head(deg.capillaries[[1]])

# combine a table of logFC
tab.deg.cap = merge(deg.capillaries[[3]][,c(2,5)], deg.capillaries[[4]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.cap) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[5]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[1]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.cap)

rownames(tab.deg.cap) = tab.deg.cap$gene
tab.deg.cap = tab.deg.cap[,-1]
nrow(tab.deg.cap) # 9524 genes
write.table(tab.deg.cap,"table_deg_kidney_kidney2_all_TP.txt", sep='\t')

# 3. Glomeruli Capillaries
setwd("../organs/kidney/deg-clust")
tabs = list.files()[grep("EC-glomeruli_", list.files())]
deg.capillaries = list()
deg.i = NULL
for (i in seq(1,5,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.capillaries[[i]] = deg.i
}

head(deg.capillaries[[1]])

# combine a table of logFC
tab.deg.cap = merge(deg.capillaries[[3]][,c(2,5)], deg.capillaries[[4]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.cap) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[5]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[1]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.cap)

rownames(tab.deg.cap) = tab.deg.cap$gene
tab.deg.cap = tab.deg.cap[,-1]
nrow(tab.deg.cap) # 9524 genes
write.table(tab.deg.cap,"table_deg_kidney_glomeruli_all_TP.txt", sep='\t')

# bind all clusters, WD 6m
setwd("../organs/kidney/deg-clust")
tabs = list.files()[grep("Western6m_vs_chow", list.files())]
deg.capillaries = list()
deg.i = NULL
for (i in seq(1,length(tabs),1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.capillaries[[i]] = deg.i
}

head(deg.capillaries[[1]])

# combine a table of logFC
tab.deg.cap = merge(deg.capillaries[[1]][,c(2,5)], deg.capillaries[[2]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.cap) = c("gene",tabs[1],"Q.val1",tabs[2],"Q.val2")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[3]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[6:7] = c(tabs[3],"Q.val3")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[4]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[8:9] = c(tabs[4],"Q.val4")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[5]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[10:11] = c(tabs[5],"Q.val5")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[6]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[12:13] = c(tabs[6],"Q.val6")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[7]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[14:15] = c(tabs[7],"Q.val7")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[8]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[16:17] = c(tabs[8],"Q.val8")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[9]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[18:19] = c(tabs[9],"Q.val9")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[10]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[20:21] = c(tabs[10],"Q.val10")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[11]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[22:23] = c(tabs[11],"Q.val11")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[12]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[24:25] = c(tabs[12],"Q.val12")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[13]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[26:27] = c(tabs[13],"Q.val13")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[14]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[28:29] = c(tabs[14],"Q.val14")
head(tab.deg.cap)

rownames(tab.deg.cap) = tab.deg.cap$gene
tab.deg.cap = tab.deg.cap[,-1]
nrow(tab.deg.cap) # 8950 genes
head(tab.deg.cap)

# remove string from colnames 
library(stringr)
colnames(tab.deg.cap) = str_remove(colnames(tab.deg.cap),"_Western6m_vs_chow.txt")

write.table(tab.deg.cap,"table_deg_kidney_clust_all_WD6m.txt", sep='\t')



####
# combine Cap/ Art/ Ven in 1 table

# 1. Capillaries
tabs = list.files()[grep("kidney_Cap_", list.files())]
deg.capillaries = list()
deg.i = NULL
for (i in seq(1,4,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.capillaries[[i]] = deg.i
}

deg.capillaries[[5]] = read.table("/Users/olga.bondareva/Desktop/HI-MAG OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new/deg_gen/kidney_EC-cap_Western_vs_chow.txt", sep = '\t', header = T)

head(deg.capillaries[[1]])

# combine a table of logFC
tab.deg.cap = merge(deg.capillaries[[5]][,c(2,5)], deg.capillaries[[3]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.cap) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[4]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[1]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.cap)

rownames(tab.deg.cap) = tab.deg.cap$gene
tab.deg.cap = tab.deg.cap[,-1]
nrow(tab.deg.cap) # 7842 genes
write.table(tab.deg.cap,"table_deg_kidney_Cap_all_TP.txt", sep='\t')

# signif DEGs
tab.deg.cap.sign = tab.deg.cap[(tab.deg.cap$WD_3m_Q.val<0.05)|(tab.deg.cap$WD_4m_Q.val<0.05)|(tab.deg.cap$WD_6m_Q.val<0.05)|(tab.deg.cap$Rev_1m_Q.val<0.05)|(tab.deg.cap$Rev_3m_Q.val<0.05),]
# rm NAs
tab.deg.cap.sign = tab.deg.cap.sign[complete.cases(tab.deg.cap.sign),]
nrow(tab.deg.cap.sign) # 3175 genes with signif Q-val
head(tab.deg.cap.sign)
write.table(tab.deg.cap.sign, "table_deg_kidney_Cap_all_TP_sign.txt", sep = '\t')

# signif only at 3m TP
tab.deg.cap.sign = tab.deg.cap[(tab.deg.cap$WD_3m_Q.val<0.05),]
# rm NAs
tab.deg.cap.sign = tab.deg.cap.sign[complete.cases(tab.deg.cap.sign),]
nrow(tab.deg.cap.sign) # 274 genes with signif Q-val
head(tab.deg.cap.sign)

# keep only logFC
tab.deg.cap.logFC = tab.deg.cap.sign[,seq(1,10,2)] # tab.deg.cap[,seq(1,14,2)]
head(tab.deg.cap.logFC)

# replace NA with 0s
tab.deg.cap.logFC[is.na(tab.deg.cap.logFC)] = 0
head(tab.deg.cap.logFC)

# remove rows with all 0s
tab.deg.cap.logFC = tab.deg.cap.logFC[!rowSums(tab.deg.cap.logFC)==0,]
head(tab.deg.cap.logFC)
nrow(tab.deg.cap.logFC) # 4237 genes with logFC
colnames(tab.deg.cap.logFC) = c("WD_3m","WD_4m","WD_6m","Rev_1m","Rev_3m")
write.table(tab.deg.cap.logFC, "table_deg_kidney_Cap_all_TP_sign_LOG.txt", sep = '\t')

b <- tab.deg.cap.logFC
# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]
# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5
mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma
ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
          main = "kidney Cap DEG, logFC(x/chow)", 
          Rowv= T,Colv= T, 
          ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="both", cexRow=0.8, cexCol=0.8) 

# plot correlations

pdf(paste0("corr-cap-kidney.pdf"),width=10,height=10,paper='special')
p1 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,WD_4m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_4m") + theme_linedraw() + geom_smooth(method = "lm")
p2 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,WD_6m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_6m") + theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,Rev_1m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_1m") + theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,Rev_3m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_31m") + theme_linedraw() + geom_smooth(method = "lm")

print((p1 | p2 ) /
        ( p3 |p4)) # 10 x 10
dev.off()

m <- cor(tab.deg.cap.logFC)
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black",
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))


############
setwd(paste0(file.dir,"/kidney/deg-gen"))
# 2. Artery
tabs = list.files()[grep("kidney_Art_", list.files())]
deg.artery = list()
deg.i = NULL
for (i in seq(1,4,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.artery[[i]] = deg.i
}

deg.artery[[5]] = read.table("/Users/olga.bondareva/Desktop/HI-MAG OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new/deg_gen/kidney_EC-art_Western_vs_chow.txt", sep = '\t', header = T)

head(deg.artery[[1]])

# combine a table of logFC
tab.deg.art = merge(deg.artery[[5]][,c(2,5)], deg.artery[[3]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.art) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.art = merge(tab.deg.art, deg.artery[[4]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.art = merge(tab.deg.art, deg.artery[[1]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.art = merge(tab.deg.art, deg.artery[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.art)

rownames(tab.deg.art) = tab.deg.art$gene
tab.deg.art = tab.deg.art[,-1]
nrow(tab.deg.art) # 6848 genes
write.table(tab.deg.art,"table_deg_kidney_art_all_TP.txt", sep='\t')

# signif DEGs
tab.deg.art.sign = tab.deg.art[(tab.deg.art$WD_3m_Q.val<0.05)|(tab.deg.art$WD_4m_Q.val<0.05)|(tab.deg.art$WD_6m_Q.val<0.05)|(tab.deg.art$Rev_1m_Q.val<0.05)|(tab.deg.art$Rev_3m_Q.val<0.05),]
# rm NAs
tab.deg.art.sign = tab.deg.art.sign[complete.cases(tab.deg.art.sign),]
nrow(tab.deg.art.sign) # 1671 genes with signif Q-val
head(tab.deg.art.sign)
write.table(tab.deg.art.sign, "table_deg_kidney_art_all_TP_sign.txt", sep = '\t')

# keep only logFC
tab.deg.art.logFC = tab.deg.art.sign[,seq(1,10,2)] # tab.deg.art[,seq(1,14,2)]
head(tab.deg.art.logFC)

# replace NA with 0s
tab.deg.art.logFC[is.na(tab.deg.art.logFC)] = 0
head(tab.deg.art.logFC)

# remove rows with all 0s
tab.deg.art.logFC = tab.deg.art.logFC[!rowSums(tab.deg.art.logFC)==0,]
head(tab.deg.art.logFC)
nrow(tab.deg.art.logFC) # 591 genes with logFC
colnames(tab.deg.art.logFC) = c("WD_3m","WD_4m","WD_6m","Rev_1m","Rev_3m")
write.table(tab.deg.art.logFC, "table_deg_kidney_art_all_TP_sign_LOG.txt", sep = '\t')

b <- tab.deg.art.logFC
# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]
# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5
mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma
ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))

pdf(paste0("heatmap-deg-art-sign-kid.pdf"),width=10,height=15,paper='special')
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
          main = "kidney art DEG, logFC(x/chow)", 
          Rowv= T,Colv= T, 
          ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="both", cexRow=0.8, cexCol=0.8) 
dev.off()
# plot correlations

pdf(paste0("corr-art-kid.pdf"),width=10,height=10,paper='special')
p1 <- ggplot(tab.deg.art.logFC, aes(WD_3m,WD_4m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_4m") + theme_linedraw() + geom_smooth(method = "lm")
p2 <- ggplot(tab.deg.art.logFC, aes(WD_3m,WD_6m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_6m") + theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.art.logFC, aes(WD_3m,Rev_1m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_1m") + theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.art.logFC, aes(WD_3m,Rev_3m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_31m") + theme_linedraw() + geom_smooth(method = "lm")

print((p1 | p2 ) /
        ( p3 |p4)) # 10 x 10
dev.off()

# all clusters
m <- cor(tab.deg.art.logFC)
pdf(paste0("corr-art-tab-kid.pdf"),width=7,height=7,paper='special')
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black",
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))
dev.off()

############
# 3. Vein
tabs = list.files()[grep("kidney_Ven_", list.files())]
deg.vein = list()
deg.i = NULL
for (i in seq(1,4,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.vein[[i]] = deg.i
}

deg.vein[[5]] = read.table("/Users/olga.bondareva/Desktop/HI-MAG OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new/deg_gen/kidney_EC-ven_Western_vs_chow.txt", sep = '\t', header = T)

head(deg.vein[[1]])

# combine a table of logFC
tab.deg.ven = merge(deg.vein[[5]][,c(2,5)], deg.vein[[3]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.ven) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.ven = merge(tab.deg.ven, deg.vein[[4]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.ven = merge(tab.deg.ven, deg.vein[[1]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.ven = merge(tab.deg.ven, deg.vein[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.ven)

rownames(tab.deg.ven) = tab.deg.ven$gene
tab.deg.ven = tab.deg.ven[,-1]
nrow(tab.deg.ven) # 6848 genes
write.table(tab.deg.ven,"table_deg_kidney_ven_all_TP.txt", sep='\t')

# signif DEGs
tab.deg.ven.sign = tab.deg.ven[(tab.deg.ven$WD_3m_Q.val<0.05)|(tab.deg.ven$WD_4m_Q.val<0.05)|(tab.deg.ven$WD_6m_Q.val<0.05)|(tab.deg.ven$Rev_1m_Q.val<0.05)|(tab.deg.ven$Rev_3m_Q.val<0.05),]
# rm NAs
tab.deg.ven.sign = tab.deg.ven.sign[complete.cases(tab.deg.ven.sign),]
nrow(tab.deg.ven.sign) # 1671 genes with signif Q-val
head(tab.deg.ven.sign)
write.table(tab.deg.ven.sign, "table_deg_kidney_ven_all_TP_sign.txt", sep = '\t')

# keep only logFC
tab.deg.ven.logFC = tab.deg.ven.sign[,seq(1,10,2)] # tab.deg.ven[,seq(1,14,2)]
head(tab.deg.ven.logFC)

# replace NA with 0s
tab.deg.ven.logFC[is.na(tab.deg.ven.logFC)] = 0
head(tab.deg.ven.logFC)

# remove rows with all 0s
tab.deg.ven.logFC = tab.deg.ven.logFC[!rowSums(tab.deg.ven.logFC)==0,]
head(tab.deg.ven.logFC)
nrow(tab.deg.ven.logFC) # 591 genes with logFC
colnames(tab.deg.ven.logFC) = c("WD_3m","WD_4m","WD_6m","Rev_1m","Rev_3m")
write.table(tab.deg.ven.logFC, "table_deg_kidney_ven_all_TP_sign_LOG.txt", sep = '\t')

b <- tab.deg.ven.logFC
# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]
# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5
mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma
ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))

pdf(paste0("heatmap-deg-ven-sign-kid.pdf"),width=10,height=15,paper='special')
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
          main = "kidney ven DEG, logFC(x/chow)", 
          Rowv= T,Colv= T, 
          ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="both", cexRow=0.8, cexCol=0.8) 
dev.off()
# plot correlations

pdf(paste0("corr-ven-kid.pdf"),width=10,height=10,paper='special')
p1 <- ggplot(tab.deg.ven.logFC, aes(WD_3m,WD_4m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_4m") + theme_linedraw() + geom_smooth(method = "lm")
p2 <- ggplot(tab.deg.ven.logFC, aes(WD_3m,WD_6m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_6m") + theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.ven.logFC, aes(WD_3m,Rev_1m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_1m") + theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.ven.logFC, aes(WD_3m,Rev_3m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_31m") + theme_linedraw() + geom_smooth(method = "lm")

print((p1 | p2 ) /
        ( p3 |p4)) # 10 x 10
dev.off()

# all clusters
m <- cor(tab.deg.ven.logFC)
pdf(paste0("corr-ven-tab-kid.pdf"),width=7,height=7,paper='special')
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black",
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))
dev.off()

saveRDS(kidney, "kidney_3m_4m_6m.rds")

#########################
# 6. lung
#########################

setwd("lung")
lung.3m = readRDS(paste0(data.dir.3m,"/lung/lung_ECs_singlet.rds"))
lung.4m = readRDS(paste0(data.dir.4m,"lung_subset.rds"))
lung.6m = readRDS(paste0(data.dir.6m,"lung_subset.rds"))

lung.3m$time = "3m"

# merge objects
lung = merge(lung.3m, list(lung.4m, lung.6m))
rm(lung.3m, lung.4m, lung.6m)

lung$diet_time = paste0(lung$diet, lung$time)

lung <- ScaleData(object = lung, verbose = T)
lung <- FindVariableFeatures(lung, selection.method = "vst", nfeatures = 3000)
lung <- RunPCA(object = lung, npcs = 30, verbose = T)
ElbowPlot(object = lung,  ndims = 30) 

# UMAP and Clustering
lung <- RunUMAP(object = lung, reduction = "pca", dims = 1:25)
lung <- FindNeighbors(object = lung, reduction = "pca", dims = 1:25)
lung <- FindClusters(object = lung, resolution = 0.3)

# visualization
DimPlot(object = lung, reduction = "umap", label = TRUE)
DimPlot(object = lung, reduction = "umap", group.by = "diet_time")
DimPlot(object = lung, reduction = "umap", group.by = "lung_celltype", label = T)


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

# Pericyte amd SMC markers
FeaturePlot(lung, features = c("Flt1", "Pecam1","Cdh5",# EC
                               "Pdgfrb", "Cspg4","Anpep",  #  Pericyte
                               "Myh11", "Acta2", "Des",   # SMC
                               "Dcn", "Col1a1", "Pdgfra"), # FB
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10") #, min.cutoff = "q10"

FeaturePlot(lung, features = c("Flt1", "Pecam1","Cdh5",# EC
                               "Ptprc","Csf1r","Cd19",
                               "Nkg7","Itgam","Mrc1",
                               "Itgax","Itgae","H2-Eb1"), # immune
            order = T, cols=c("grey", "red"), ncol=3, min.cutoff = "q10") #, min.cutoff = "q10"


FeaturePlot(lung, features = c("Prox1", "Pdpn","Lyve1", "Mki67","Cenpf"), # lymph, prolif
            order = T, cols=c("grey", "red"), ncol=2, min.cutoff = "q10") #, min.cutoff = "q10"

lung = readRDS("lung_3m_4m_6m.rds")
Idents(lung) = lung$lung_celltype2 
DimPlot(object = lung, reduction = "umap", label = TRUE) + NoLegend()

FeaturePlot(lung, features = c("Aqp5", "Spp1","Pdgfrb","Anpep"), # 
            order = T, cols=c("grey", "red"), ncol=2, min.cutoff = "q10") #, min.cutoff = "q10"

lung <- RenameIdents(lung, `4` = "EC-ven",`0` = "EC-cap1",`2` = "EC-cap2",`5` = "EC-cap3",
                     `8` = "EC-art",`12` = "EC-platelet", `15` = "Immune",
                     `14` = "Prolif", `10` = "EC-lymph",`11` = "Pneumocyte Type II",
                      `9` = "EC-Aqp5",`7` = "Alveolar FB", `13` = "Pericyte",
                     `1` = "EC-pulm1", `3` = "EC-pulm2a",`6` = "EC-pulm2b") 

lung <- RenameIdents(lung, `Pneumocyte Type II` = "EC-pneumocyte",
                     `EC-pulm1` = "aEC", `EC-pulm2a` = "pulmEC_a",`EC-pulm2b` = "pulmEC_b") 

## Relevel object@ident
lung@active.ident <- factor(x = lung@active.ident, levels = sort(levels(Idents(lung))))
lung$lung_celltype2 <- Idents(lung)
DimPlot(lung, reduction = "umap", label = TRUE) + NoLegend()

# number of cells per cluster per sample
table(Idents(lung), lung$diet_time)

# remove immune cells
lung = subset(lung, idents = c("Immune","Alveolar FB","Pericyte"), invert = T)

lung <- ScaleData(object = lung, verbose = T)
lung <- FindVariableFeatures(lung, selection.method = "vst", nfeatures = 3000)
lung <- RunPCA(object = lung, npcs = 30, verbose = T)
ElbowPlot(object = lung,  ndims = 30) 

# UMAP and Clustering
lung <- RunUMAP(object = lung, reduction = "pca", dims = 1:25)
lung <- FindNeighbors(object = lung, reduction = "pca", dims = 1:25)
lung <- FindClusters(object = lung, resolution = 0.3)

# visualization
DimPlot(object = lung, reduction = "umap", label = TRUE)

Idents(lung) = lung$lung_celltype2
DimPlot(object = lung, reduction = "umap", label = TRUE)
DimPlot(lung, reduction = "umap", label = TRUE, split.by = "diet_time", ncol = 3, pt.size = 1.3) + NoLegend()


lung = readRDS("lung_3m_4m_6m.rds")
rm# DE genes
lung$clust.diet_time <- paste(lung$lung_celltype2, lung$diet_time, sep = "_")
Idents(object = lung) <- "clust.diet_time"
lung@active.ident <- factor(x = lung@active.ident, levels = sort(levels(Idents(object = lung))))
levels(Idents(object = lung))


lev1 <- levels(lung@active.ident)[seq(7,(length(levels(Idents(lung)))),8)] # condition 1=HFD_4m
lev2 <- levels(lung@active.ident)[seq(2,(length(levels(Idents(lung)))),8)] # condition 2=chow_4m
lev3 <- levels(lung@active.ident)[seq(4,(length(levels(Idents(lung)))),8)] # condition 3=Rev_1m
lev4 <- levels(lung@active.ident)[seq(8,(length(levels(Idents(lung)))),8)] # condition 4=HFD_6m
lev5 <- levels(lung@active.ident)[seq(3,(length(levels(Idents(lung)))),8)] # condition 5=chow_6m
lev6 <- levels(lung@active.ident)[seq(5,(length(levels(Idents(lung)))),8)] # condition 6=Rev_3m
lev7 <- levels(lung@active.ident)[seq(1,(length(levels(Idents(lung)))),8)] # condition 7=HFD_3m
lev8 <- levels(lung@active.ident)[seq(6,(length(levels(Idents(lung)))),8)] # condition 8=chow_3m

# for all the clusters
for (i in seq(1,length(lev1),1)) {
  # 4 months
  # Western
  print(lev1[i])
  print(lev2[i])
  response2 <- FindMarkers(object = lung, ident.1 = lev1[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev1[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # Reverse
  print(lev3[i])
  print(lev2[i])
  response2 <- FindMarkers(object = lung, ident.1 = lev3[i], 
                           ident.2 = lev2[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev3[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # 6 months
  # Western
  print(lev4[i])
  print(lev5[i])
  response2 <- FindMarkers(object = lung, ident.1 = lev4[i], 
                           ident.2 = lev5[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev4[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # Reverse
  print(lev6[i])
  print(lev5[i])
  response2 <- FindMarkers(object = lung, ident.1 = lev6[i], 
                           ident.2 = lev5[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev6[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
  # 3 months
  # Western
  print(lev7[i])
  print(lev8[i])
  response2 <- FindMarkers(object = lung, ident.1 = lev7[i], 
                           ident.2 = lev8[i], 
                           logfc.threshold = 0,
                           verbose = T) #, min.pct = 0.25
  write.table(response2, file = paste(lev7[i],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  
}


# DEG for capillaries combined
# 4m
response2 <- FindMarkers(object = lung, 
                         ident.1 = c("EC-cap1_Western4m","EC-cap2_Western4m"),  # Cap3 is different
                         ident.2 = c("EC-cap1_chow4m","EC-cap2_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "lung_Cap_WB_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = lung, 
                         ident.1 = c("EC-cap1_Rev1m","EC-cap2_Rev1m"), 
                         ident.2 = c("EC-cap1_chow4m","EC-cap2_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "lung_Cap_Rev1m_vs_chow_4m.txt", sep = "\t")

# 6m
response2 <- FindMarkers(object = lung, 
                         ident.1 = c("EC-cap1_Western6m","EC-cap2_Western6m"), 
                         ident.2 = c("EC-cap1_chow6m","EC-cap2_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "lung_Cap_WB_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = lung, 
                         ident.1 = c("EC-cap1_Rev3m","EC-cap2_Rev3m"), 
                         ident.2 = c("EC-cap1_chow6m","EC-cap2_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "lung_Cap_Rev3m_vs_chow_6m.txt", sep = "\t")

# DEG for vein and artery combined
# 4m
response2 <- FindMarkers(object = lung, 
                         ident.1 = c("EC-art_Western4m"), 
                         ident.2 = c("EC-art_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "lung_Art_WB_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = lung, 
                         ident.1 = c("EC-ven_Western4m"), 
                         ident.2 = c("EC-ven_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "lung_Ven_WB_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = lung, 
                         ident.1 = c("EC-art_Rev1m"), 
                         ident.2 = c("EC-art_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "lung_Art_Rev1m_vs_chow_4m.txt", sep = "\t")

response2 <- FindMarkers(object = lung, 
                         ident.1 = c("EC-ven_Rev1m"), 
                         ident.2 = c("EC-ven_chow4m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "lung_Ven_Rev1m_vs_chow_4m.txt", sep = "\t")

# 6m
response2 <- FindMarkers(object = lung, 
                         ident.1 = c("EC-art_Western6m"), 
                         ident.2 = c("EC-art_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "lung_Art_WB_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = lung, 
                         ident.1 = c("EC-ven_Western6m"), 
                         ident.2 = c("EC-ven_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "lung_Ven_WB_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = lung, 
                         ident.1 = c("EC-art_Rev3m"), 
                         ident.2 = c("EC-art_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "lung_Art_Rev3m_vs_chow_6m.txt", sep = "\t")

response2 <- FindMarkers(object = lung, 
                         ident.1 = c("EC-ven_Rev3m"), 
                         ident.2 = c("EC-ven_chow6m"), 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
write.table(response2, file = "lung_Ven_Rev3m_vs_chow_6m.txt", sep = "\t")


saveRDS(lung, "lung_3m_4m_6m.rds")


# 1. Pulm1 Capillaries
setwd("deg-clust")
tabs = list.files()[grep("EC-pulm1_", list.files())]
deg.capillaries = list()
deg.i = NULL
for (i in seq(2,5,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.capillaries[[i]] = deg.i
}
# 3m - change - to +, deg chow vs. WD
deg.capillaries[[1]] = read.table(tabs[1], sep = '\t', header = T)
deg.capillaries[[1]]$avg_logFC = -deg.capillaries[[1]]$avg_logFC

head(deg.capillaries[[1]])

# combine a table of logFC
tab.deg.cap = merge(deg.capillaries[[1]][,c(2,5)], deg.capillaries[[4]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.cap) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[5]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[3]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.cap)

rownames(tab.deg.cap) = tab.deg.cap$gene
tab.deg.cap = tab.deg.cap[,-1]
nrow(tab.deg.cap) # 8557 genes
write.table(tab.deg.cap,"table_deg_lung_pulm1_all_TP.txt", sep='\t')


####
# combine Cap/ Art/ Ven in 1 table

# 1. Capillaries
tabs = list.files()[grep("lung_Cap_", list.files())]
deg.capillaries = list()
deg.i = NULL
for (i in seq(1,4,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.capillaries[[i]] = deg.i
}

deg.capillaries[[5]] = read.table("/Users/olga.bondareva/Desktop/HI-MAG OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new/deg_gen/lung_EC-cap_Western_vs_chow.txt", sep = '\t', header = T)

head(deg.capillaries[[1]])

# combine a table of logFC
tab.deg.cap = merge(deg.capillaries[[5]][,c(2,5)], deg.capillaries[[3]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.cap) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[4]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[1]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.cap)

rownames(tab.deg.cap) = tab.deg.cap$gene
tab.deg.cap = tab.deg.cap[,-1]
nrow(tab.deg.cap) # 6848 genes
write.table(tab.deg.cap,"table_deg_lung_Cap_all_TP.txt", sep='\t')

# signif DEGs
tab.deg.cap.sign = tab.deg.cap[(tab.deg.cap$WD_3m_Q.val<0.05)|(tab.deg.cap$WD_4m_Q.val<0.05)|(tab.deg.cap$WD_6m_Q.val<0.05)|(tab.deg.cap$Rev_1m_Q.val<0.05)|(tab.deg.cap$Rev_3m_Q.val<0.05),]
# rm NAs
tab.deg.cap.sign = tab.deg.cap.sign[complete.cases(tab.deg.cap.sign),]
nrow(tab.deg.cap.sign) # 1671 genes with signif Q-val
head(tab.deg.cap.sign)
write.table(tab.deg.cap.sign, "table_deg_lung_Cap_all_TP_sign.txt", sep = '\t')

####
# signif only at 3m TP
tab.deg.cap.sign = tab.deg.cap[(tab.deg.cap$WD_3m_Q.val<0.05),]
# rm NAs
tab.deg.cap.sign = tab.deg.cap.sign[complete.cases(tab.deg.cap.sign),]
nrow(tab.deg.cap.sign) # 274 genes with signif Q-val
head(tab.deg.cap.sign)
#####

# keep only logFC
tab.deg.cap.logFC = tab.deg.cap.sign[,seq(1,10,2)] # tab.deg.cap[,seq(1,14,2)]
head(tab.deg.cap.logFC)

# replace NA with 0s
tab.deg.cap.logFC[is.na(tab.deg.cap.logFC)] = 0
head(tab.deg.cap.logFC)

# remove rows with all 0s
tab.deg.cap.logFC = tab.deg.cap.logFC[!rowSums(tab.deg.cap.logFC)==0,]
head(tab.deg.cap.logFC)
nrow(tab.deg.cap.logFC) # 591 genes with logFC
colnames(tab.deg.cap.logFC) = c("WD_3m","WD_4m","WD_6m","Rev_1m","Rev_3m")
write.table(tab.deg.cap.logFC, "table_deg_lung_Cap_all_TP_sign_LOG.txt", sep = '\t')

b <- tab.deg.cap.logFC
# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]
# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5
mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma
ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
          main = "lung Cap DEG, logFC(x/chow)", 
          Rowv= T,Colv= T, 
          ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="both", cexRow=0.8, cexCol=0.8) 

# plot correlations

pdf(paste0("corr-cap-lung.pdf"),width=10,height=10,paper='special')
p1 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,WD_4m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_4m") + theme_linedraw() + geom_smooth(method = "lm")
p2 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,WD_6m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_6m") + theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,Rev_1m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_1m") + theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.cap.logFC, aes(WD_3m,Rev_3m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_31m") + theme_linedraw() + geom_smooth(method = "lm")

print((p1 | p2 ) /
        ( p3 |p4)) # 10 x 10
dev.off()

# all clusters
library(corrplot)
m <- cor(tab.deg.cap.logFC)
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black",
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))

############
setwd("../lung/deg-gen")
# 2. Artery
tabs = list.files()[grep("lung_Art_", list.files())]
deg.artery = list()
deg.i = NULL
for (i in seq(1,4,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.artery[[i]] = deg.i
}

deg.artery[[5]] = read.table("/Users/olga.bondareva/Desktop/HI-MAG OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new/deg_gen/lung_EC-art_Western_vs_chow.txt", sep = '\t', header = T)

head(deg.artery[[1]])

# combine a table of logFC
tab.deg.art = merge(deg.artery[[5]][,c(2,5)], deg.artery[[3]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.art) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.art = merge(tab.deg.art, deg.artery[[4]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.art = merge(tab.deg.art, deg.artery[[1]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.art = merge(tab.deg.art, deg.artery[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.art)

rownames(tab.deg.art) = tab.deg.art$gene
tab.deg.art = tab.deg.art[,-1]
nrow(tab.deg.art) # 6848 genes
write.table(tab.deg.art,"table_deg_lung_art_all_TP.txt", sep='\t')

# signif DEGs
tab.deg.art.sign = tab.deg.art[(tab.deg.art$WD_3m_Q.val<0.05)|(tab.deg.art$WD_4m_Q.val<0.05)|(tab.deg.art$WD_6m_Q.val<0.05)|(tab.deg.art$Rev_1m_Q.val<0.05)|(tab.deg.art$Rev_3m_Q.val<0.05),]
# rm NAs
tab.deg.art.sign = tab.deg.art.sign[complete.cases(tab.deg.art.sign),]
nrow(tab.deg.art.sign) # 1671 genes with signif Q-val
head(tab.deg.art.sign)
write.table(tab.deg.art.sign, "table_deg_lung_art_all_TP_sign.txt", sep = '\t')

# keep only logFC
tab.deg.art.logFC = tab.deg.art.sign[,seq(1,10,2)] # tab.deg.art[,seq(1,14,2)]
head(tab.deg.art.logFC)

# replace NA with 0s
tab.deg.art.logFC[is.na(tab.deg.art.logFC)] = 0
head(tab.deg.art.logFC)

# remove rows with all 0s
tab.deg.art.logFC = tab.deg.art.logFC[!rowSums(tab.deg.art.logFC)==0,]
head(tab.deg.art.logFC)
nrow(tab.deg.art.logFC) # 591 genes with logFC
colnames(tab.deg.art.logFC) = c("WD_3m","WD_4m","WD_6m","Rev_1m","Rev_3m")
write.table(tab.deg.art.logFC, "table_deg_lung_art_all_TP_sign_LOG.txt", sep = '\t')

b <- tab.deg.art.logFC
# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]
# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5
mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma
ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
          main = "lung art DEG, logFC(x/chow)", 
          Rowv= T,Colv= T, 
          ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="both", cexRow=0.8, cexCol=0.8) 

# plot correlations

pdf(paste0("corr-art-lung.pdf"),width=10,height=10,paper='special')
p1 <- ggplot(tab.deg.art.logFC, aes(WD_3m,WD_4m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_4m") + theme_linedraw() + geom_smooth(method = "lm")
p2 <- ggplot(tab.deg.art.logFC, aes(WD_3m,WD_6m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_6m") + theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.art.logFC, aes(WD_3m,Rev_1m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_1m") + theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.art.logFC, aes(WD_3m,Rev_3m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_31m") + theme_linedraw() + geom_smooth(method = "lm")

print((p1 | p2 ) /
        ( p3 |p4)) # 10 x 10
dev.off()

# all clusters
m <- cor(tab.deg.art.logFC)
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black",
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))

############
# 3. Vein
tabs = list.files()[grep("lung_Ven_", list.files())]
deg.vein = list()
deg.i = NULL
for (i in seq(1,4,1)) {
  print(i)
  deg.i = read.table(tabs[i], sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.vein[[i]] = deg.i
}

deg.vein[[5]] = read.table("/Users/olga.bondareva/Desktop/HI-MAG OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new/deg_gen/lung_EC-ven_Western_vs_chow.txt", sep = '\t', header = T)

head(deg.vein[[1]])

# combine a table of logFC
tab.deg.ven = merge(deg.vein[[5]][,c(2,5)], deg.vein[[3]][,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.ven) = c("gene","WD_3m_logFC","WD_3m_Q.val","WD_4m_logFC","WD_4m_Q.val")
tab.deg.ven = merge(tab.deg.ven, deg.vein[[4]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[6:7] = c("WD_6m_logFC","WD_6m_Q.val")
tab.deg.ven = merge(tab.deg.ven, deg.vein[[1]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[8:9] = c("Rev_1m_logFC","Rev_1m_Q.val")
tab.deg.ven = merge(tab.deg.ven, deg.vein[[2]][,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[10:11] = c("Rev_3m_logFC","Rev_3m_Q.val")

head(tab.deg.ven)

rownames(tab.deg.ven) = tab.deg.ven$gene
tab.deg.ven = tab.deg.ven[,-1]
nrow(tab.deg.ven) # 6848 genes
write.table(tab.deg.ven,"table_deg_lung_ven_all_TP.txt", sep='\t')

# signif DEGs
tab.deg.ven.sign = tab.deg.ven[(tab.deg.ven$WD_3m_Q.val<0.05)|(tab.deg.ven$WD_4m_Q.val<0.05)|(tab.deg.ven$WD_6m_Q.val<0.05)|(tab.deg.ven$Rev_1m_Q.val<0.05)|(tab.deg.ven$Rev_3m_Q.val<0.05),]
# rm NAs
tab.deg.ven.sign = tab.deg.ven.sign[complete.cases(tab.deg.ven.sign),]
nrow(tab.deg.ven.sign) # 1671 genes with signif Q-val
head(tab.deg.ven.sign)
write.table(tab.deg.ven.sign, "table_deg_lung_ven_all_TP_sign.txt", sep = '\t')

# keep only logFC
tab.deg.ven.logFC = tab.deg.ven.sign[,seq(1,10,2)] # tab.deg.ven[,seq(1,14,2)]
head(tab.deg.ven.logFC)

# replace NA with 0s
tab.deg.ven.logFC[is.na(tab.deg.ven.logFC)] = 0
head(tab.deg.ven.logFC)

# remove rows with all 0s
tab.deg.ven.logFC = tab.deg.ven.logFC[!rowSums(tab.deg.ven.logFC)==0,]
head(tab.deg.ven.logFC)
nrow(tab.deg.ven.logFC) # 591 genes with logFC
colnames(tab.deg.ven.logFC) = c("WD_3m","WD_4m","WD_6m","Rev_1m","Rev_3m")
write.table(tab.deg.ven.logFC, "table_deg_lung_ven_all_TP_sign_LOG.txt", sep = '\t')

b <- tab.deg.ven.logFC
# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]
# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5
mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma
ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
          main = "lung ven DEG, logFC(x/chow)", 
          Rowv= T,Colv= T, 
          ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="both", cexRow=0.8, cexCol=0.8) 

# plot correlations

pdf(paste0("corr-ven-lung.pdf"),width=10,height=10,paper='special')
p1 <- ggplot(tab.deg.ven.logFC, aes(WD_3m,WD_4m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_4m") + theme_linedraw() + geom_smooth(method = "lm")
p2 <- ggplot(tab.deg.ven.logFC, aes(WD_3m,WD_6m)) + geom_point() + ggtitle("Correlation WD_3m vs WD_6m") + theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.ven.logFC, aes(WD_3m,Rev_1m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_1m") + theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.ven.logFC, aes(WD_3m,Rev_3m)) + geom_point() + ggtitle("Correlation WD_3m vs Rev_31m") + theme_linedraw() + geom_smooth(method = "lm")

print((p1 | p2 ) /
        ( p3 |p4)) # 10 x 10
dev.off()

# all clusters
m <- cor(tab.deg.ven.logFC)
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black",
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))


saveRDS(lung, "lung_3m_4m_6m.rds")
