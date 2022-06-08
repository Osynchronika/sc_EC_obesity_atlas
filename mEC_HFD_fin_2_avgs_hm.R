library(Seurat)
library(dplyr)
library(Matrix)
require(gplots)
require(ggplot2)
library(RColorBrewer)
library(cowplot)
require(scales)
library(monocle)
library(Rfast)

setwd("~/Desktop/HI-MAG OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/organs-2nd-processing")

# Load pre-processed data
brain = readRDS("../comparisons/brain_ECs.rds")
heart = readRDS("../comparisons/heart_ECs.rds")
lung = readRDS("../comparisons/lung_ECs.rds")
liver = readRDS("../comparisons/liver_ECs.rds")
kidney = readRDS("../comparisons/kidney_ECs.rds")
vis = readRDS("../comparisons/vis_ECs.rds")
sc = readRDS("../comparisons/sc_ECs.rds")


#1. Find marker genes for all clusters for all organs
# brain
Idents(brain) = brain$brain_celltype
brain.markers = FindAllMarkers(brain, only.pos = T)
write.table(brain.markers, "brain_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(brain.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "brain_top5_cluster_markers.txt", sep='\t')
pdf("dp_top5_brain.pdf",width=20,height=8,paper='special')
DotPlot(brain, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()
dev.off()

# heart
Idents(heart) = heart$heart_celltype
heart.markers = FindAllMarkers(heart, only.pos = T)
write.table(heart.markers, "heart_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(heart.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "heart_top5_cluster_markers.txt", sep='\t')
pdf("dp_top5_heart.pdf",width=20,height=8,paper='special')
DotPlot(heart, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()
dev.off()

# lung
Idents(lung) = lung$lung_celltype
lung.markers = FindAllMarkers(lung, only.pos = T)
write.table(lung.markers, "lung_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(lung.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "lung_top5_cluster_markers.txt", sep='\t')
pdf("dp_top5_lung.pdf",width=20,height=8,paper='special')
DotPlot(lung, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()
dev.off()

# liver
Idents(liver) = liver$liver_celltype
liver.markers = FindAllMarkers(liver, only.pos = T)
write.table(liver.markers, "liver_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(liver.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "liver_top5_cluster_markers.txt", sep='\t')
pdf("dp_top5_liver.pdf",width=20,height=8,paper='special')
DotPlot(liver, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()
dev.off()

# kidney
Idents(kidney) = kidney$kidney_celltype
kidney.markers = FindAllMarkers(kidney, only.pos = T)
write.table(kidney.markers, "kidney_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(kidney.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "kidney_top5_cluster_markers.txt", sep='\t')
pdf("dp_top5_kidney.pdf",width=20,height=8,paper='special')
DotPlot(kidney, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()
dev.off()

# vis
Idents(vis) = vis$vis_celltype
vis.markers = FindAllMarkers(vis, only.pos = T)
write.table(vis.markers, "vis_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(vis.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "vis_top5_cluster_markers.txt", sep='\t')
pdf("dp_top5_vis.pdf",width=20,height=8,paper='special')
DotPlot(vis, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()
dev.off()

# sc
Idents(sc) = sc$sc_celltype
sc.markers = FindAllMarkers(sc, only.pos = T)
write.table(sc.markers, "sc_cluster_markers.txt", sep='\t')

top5.mark <- as.data.frame(sc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC))
write.table(top5.mark, "sc_top5_cluster_markers.txt", sep='\t')
pdf("dp_top5_sc.pdf",width=20,height=8,paper='special')
DotPlot(sc, features = unique(top5.mark$gene), cols = c("grey", "red"))+RotatedAxis()
dev.off()

# 2. Find DEGs between main changing cluster for each organ

# liver
response <- FindMarkers(object = liver, ident.1 = "EC-cap1", 
                         ident.2 = "EC-cap2", 
                         logfc.threshold = 0,
                         verbose = FALSE) #, min.pct = 0.25
write.table(response, file = "DEG_cluster_shifts/liver_Cap1_vs_Cap2.txt", sep = "\t")

# vis AT
response <- FindMarkers(object = vis, ident.1 = "EC-cap1", 
                        ident.2 = "EC-cap2", 
                        logfc.threshold = 0,
                        verbose = FALSE) #, min.pct = 0.25
write.table(response, file = "DEG_cluster_shifts/vis_Cap1_vs_Cap2.txt", sep = "\t")


####
# calculate avg expression per cluster
names = c("brain","heart","lung","liver","kidney","vis","sc")
k = 0
for (i in c(brain,heart,lung,liver,kidney,vis,sc)) {
  k=k+1
  print(names[k])
  cluster.averages3 <- AverageExpression(object = i)
  #order the columns
  cluster.averages3[["RNA"]] <- cluster.averages3[["RNA"]][,order(colnames(cluster.averages3[["RNA"]]))]
  head(x = cluster.averages3[["RNA"]][, 1:5])  # all genes
  write.table(cluster.averages3[["RNA"]], paste0("cluster_avg_expr_",names[k],".txt"))
  
  cluster.averages4 <- AverageExpression(object = i, add.ident = "diet")
  #order the columns
  cluster.averages4[["RNA"]] <- cluster.averages4[["RNA"]][,order(colnames(cluster.averages4[["RNA"]]))]
  head(x = cluster.averages4[["RNA"]][, 1:5])  # all genes
  write.table(cluster.averages4[["RNA"]], paste0("cluster_avg_expr_diet_",names[k],".txt"))
  
  Idents(i) = "gen_celltype"
  cluster.averages1 <- AverageExpression(object = i)
  #order the columns
  cluster.averages1[["RNA"]] <- cluster.averages1[["RNA"]][,order(colnames(cluster.averages1[["RNA"]]))]
  head(x = cluster.averages1[["RNA"]][, 1:5])  # all genes
  write.table(cluster.averages1[["RNA"]], paste0("cluster_general_avg_expr_",names[k],".txt"))
  
  cluster.averages2 <- AverageExpression(object = i, add.ident = "diet")
  #order the columns
  cluster.averages2[["RNA"]] <- cluster.averages2[["RNA"]][,order(colnames(cluster.averages2[["RNA"]]))]
  head(x = cluster.averages2[["RNA"]][, 1:5])  # all genes
  write.table(cluster.averages2[["RNA"]], paste0("cluster_general_avg_expr_diet_",names[k],".txt"))
}

#####
# heatmaps of expression of GOs

# Look for inflammation signature
nfkb = read.table("../NFkB.txt", sep='\t', header = T)
nfkb = nfkb$Symbol
nfkb = nfkb[-67]
library(stringr)
nfkb = str_to_sentence(nfkb, locale = "en")

inflamm = read.table("../GO_0006954_inflammatory_response.txt", sep='\t')
inflamm = inflamm$V1

glycolysis = c("Foxk1", "Foxk2", "Hk1", "Gpi1", "Pfkl", "Pfkm", "Pfkp")
FA.oxid = read.table("../GO_0019395_FA_oxidation.txt", sep='\t')
FA.oxid = unique(FA.oxid$V1)
FA.trans = read.table("../GO_0015908_FA_transport.txt", sep='\t')
FA.trans = unique(FA.trans$V1)

ecm = read.table("../GO_0031012_ECM.txt", sep='\t')
ecm = unique(ecm$V1)
platelet = read.table("../GO_0070527_platelet_aggregation.txt", sep='\t')
platelet = unique(platelet$V1)
platelet.act = read.table("../GO_0030168_platelet_activation.txt", sep='\t')
platelet.act = unique(platelet.act$V1)
hemostasis = read.table("../GO_0007599_hemostasis.txt", sep='\t')
hemostasis = unique(hemostasis$V1)

### plot heatmaps for GOs in organs
term = c("NFkB targets","Inflammation","Glycolysis","FA oxidation","FA transport","ECM","Platelet aggregation","Platelet activation","Hemostasis")
go.list = list(nfkb,inflamm,glycolysis,FA.oxid,FA.trans,ecm,platelet,platelet.act,hemostasis)
ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'

clust.avg1 = read.table("clust_avgs/cluster_general_avg_expr_sc.txt",sep=' ',header=T)
head(clust.avg1)
clust.avg3 = read.table("clust_avgs/cluster_avg_expr_brain.txt",sep=' ',header=T)
head(clust.avg3)

# heatmap of the avg expr
for (i in seq(1,length(term),1)) {
  # 1= gen_celltype
  b <- clust.avg1[rownames(clust.avg1) %in% unlist(go.list[i]),] 
  #FA.trans platelet hemostasis platelet.act
  # exclude expr < 0.25 genes
  b = b[rowMaxs(as.matrix(b),value = T)>=(0.25),]
  # z-score
  z <- t(scale(t(b), center=TRUE, scale=TRUE))
  mi <- min(z,na.rm=TRUE) #l
  mi
  ma <- max(z,na.rm=TRUE) #l
  ma
  
  ColorLevels <- seq(mi, ma, length=length(ColorRamp))
  
  pdf(paste0("heatmaps_GO/hm-avg-expr-gen-",term[i],"-sc.pdf"),width=12,height=18,paper='special')
  heatmap.2(as.matrix(z), na.rm = TRUE, col= ColorRamp, 
            main = paste0(term[i],", z-score of natural log expr"), #FA transport Platelet aggr Hemostasis
            Rowv= T,Colv= F, 
            #ColSideColors= my_color_palette, # rep(cell_col.both, each=2)
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(5,10),     # widens margins around plot
            dendrogram="row", cexRow=0.8, cexCol=0.8) 
  dev.off()
  
  # 3= organ_celltype
  b <- clust.avg3[rownames(clust.avg3) %in% unlist(go.list[i]),] 
  #FA.trans platelet hemostasis platelet.act
  # exclude expr < 0.25 genes
  b = b[rowMaxs(as.matrix(b),value = T)>=(0.25),]
  # z-score
  z <- t(scale(t(b), center=TRUE, scale=TRUE))
  mi <- min(z,na.rm=TRUE) #l
  mi
  ma <- max(z,na.rm=TRUE) #l
  ma
  
  ColorLevels <- seq(mi, ma, length=length(ColorRamp))
  
  pdf(paste0("heatmaps_GO/hm-avg-expr-",term[i],"-sc.pdf"),width=15,height=18,paper='special')
  heatmap.2(as.matrix(z), na.rm = TRUE, col= ColorRamp, 
            main = paste0(term[i],", z-score of natural log expr"), #FA transport Platelet aggr Hemostasis
            Rowv= T,Colv= F, 
            #ColSideColors= my_color_palette, # rep(cell_col.both, each=2)
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(5,10),     # widens margins around plot
            dendrogram="row", cexRow=0.8, cexCol=0.8) 
  dev.off()
}


###########
###
# Quantify the % of Hbb-positive cells per organ
# 1. In NormCounts
Hbb.cells = WhichCells(brain, slot = "data",expression = `Hbb-bs`>1 | `Hba-a1`>1 | `Hba-a2`>1 | `Hbb-bt`>1)
100*length(Hbb.cells)/ncol(brain)  # % of Hbb-bs cells
# Any Hb 19.76261 %
DotPlot(brain, group.by = "diet", features = c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt"), 
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
table(brain$diet[Hbb.cells])
table(brain$diet)

Hbb.cells = WhichCells(heart, slot = "counts",expression = `Hbb-bs`>1 | `Hba-a1`>1 | `Hba-a2`>1 | `Hbb-bt`>1)
100*length(Hbb.cells)/ncol(heart)  # % of Hbb-cells
# heart 1.75068 %
table(heart$diet[Hbb.cells])
table(heart$diet)

Hbb.cells = WhichCells(lung, slot = "counts",expression = `Hbb-bs`>1 | `Hba-a1`>1 | `Hba-a2`>1 | `Hbb-bt`>1)
100*length(Hbb.cells)/ncol(lung)  # % of Hbb-cells
# lung 0.5414693 %
table(lung$diet[Hbb.cells])
table(lung$diet)

Hbb.cells = WhichCells(liver, slot = "counts",expression = `Hbb-bs`>1 | `Hba-a1`>1 | `Hba-a2`>1 | `Hbb-bt`>1)
100*length(Hbb.cells)/ncol(liver)  # % of Hbb-cells
# liver 0.5616606 %
table(liver$diet[Hbb.cells])
table(liver$diet)

Hbb.cells = WhichCells(kidney, slot = "counts",expression = `Hbb-bs`>1 | `Hba-a1`>1 | `Hba-a2`>1 | `Hbb-bt`>1)
100*length(Hbb.cells)/ncol(kidney)  # % of Hbb-cells
# kidney 0.04740741 %
table(kidney$diet[Hbb.cells])
table(kidney$diet)

Hbb.cells = WhichCells(vis, slot = "counts",expression = `Hbb-bs`>1 | `Hba-a1`>1 | `Hba-a2`>1 | `Hbb-bt`>1)
100*length(Hbb.cells)/ncol(vis)  # % of Hbb-cells
# vis 0.9372453 %
table(vis$diet[Hbb.cells])
table(vis$diet)

Hbb.cells = WhichCells(sc, slot = "counts",expression = `Hbb-bs`>1 | `Hba-a1`>1 | `Hba-a2`>1 | `Hbb-bt`>1)
100*length(Hbb.cells)/ncol(sc)  # % of Hbb-cells
# sc 1.312734 %
table(sc$diet[Hbb.cells])
table(sc$diet)

#############
# determine EC-platelet
ec.platelet = c("Pf4", "Ppbp", "Nrgn")

FeaturePlot(brain, features = ec.platelet)
FeaturePlot(heart, features = ec.platelet)
FeaturePlot(lung, features = ec.platelet)
FeaturePlot(liver, features = ec.platelet)
FeaturePlot(kidney, features = ec.platelet)
FeaturePlot(vis, features = ec.platelet)
FeaturePlot(sc, features = ec.platelet)

# coagulation cascade
coagul = c("F3","F7","F8","F9","F10","F11","F12","Vwf","F2","F2r","F2rl2","F2rl3",
           "C3","C5","A2m","Serpina1","Serpina5","Serpinf2",
           "Serpine1","Serpinb2","Fga","Fgb","Fgg", # pro-coag
           "Proc","Serpina5","Serpinc1","Serpind1","Tfpi","Plg","Plau","Plat") # anti-coag

fatty.acid.trans = c("Fabp1","Fabp2","Fabp3","Fabp4","Fabp5","Fabp6","Fabp7","Cd36","Lpl","Slc27a1","Slc27a4","Dbi",
                     "Acsbg1","Acsbg2","Acsl1","Acsl2","Acsl3","Acsl4","Acsl5","Acsl6") # enzymes

##############
# calculate avg gene expression for all ECs (pseudo-bulk) combines per organ, split by diet
# calculate avg expression per cluster
names = c("brain","heart","lung","liver","kidney","vis","sc")
k = 0
bulk.avg.exp.diet = as.data.frame(matrix(nrow = 31053))
bulk.avg.exp.diet$gene = rownames(liver)
for (i in c(brain,heart,lung,liver,kidney,vis,sc)) {
  k=k+1
  print(names[k])
  Idents(i) = "organ"
  cluster.averages5 <- AverageExpression(object = i, add.ident = "diet")
  #order the columns
  cluster.averages5[["RNA"]] <- cluster.averages5[["RNA"]][,order(colnames(cluster.averages5[["RNA"]]))]
  head(x = cluster.averages5[["RNA"]])  # all genes
  cluster.averages5[["RNA"]] = cluster.averages5[["RNA"]][order(rownames(cluster.averages5[["RNA"]]),decreasing = F),]
  write.table(cluster.averages5[["RNA"]], paste0("bulk_avg_expr_",names[k],"_diet.txt"))
  cluster.averages5[["RNA"]]$gene = rownames(cluster.averages5[["RNA"]])
  bulk.avg.exp.diet = merge(bulk.avg.exp.diet,cluster.averages5[["RNA"]], by = "gene", all = T)
}

head(bulk.avg.exp.diet)
rownames(bulk.avg.exp.diet) = bulk.avg.exp.diet$gene
bulk.avg.exp.diet = bulk.avg.exp.diet[,-1]
bulk.avg.exp.diet = bulk.avg.exp.diet[,-1]
write.table(bulk.avg.exp.diet, "bulk_avg_expr_all_organs_diet.txt", sep = '\t')
