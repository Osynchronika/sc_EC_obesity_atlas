library(Seurat)
library(dplyr)
library(Matrix)
require(gplots)
require(ggplot2)
library(RColorBrewer)
library(cowplot)
require(scales)
require(scales)
library(Rfast)
library(distinct)
library(scuttle)

setwd("~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/organs-3rd-processing")

###########################
# Load pre-processed data
brain = readRDS("../comparisons/brain_ECs.rds")
heart = readRDS("../comparisons/heart_ECs.rds")
lung = readRDS("../comparisons/lung_ECs.rds")
liver = readRDS("../comparisons/liver_ECs.rds")
kidney = readRDS("../comparisons/kidney_ECs.rds")
vis = readRDS("../comparisons/vis_ECs.rds")
sc = readRDS("../comparisons/sc_ECs.rds")

#############################

# downsampling for DEG per artery/ cap/ ven

table(brain$gen_celltype,brain$diet)
table(heart$gen_celltype,heart$diet)
table(lung$gen_celltype,lung$diet)
table(kidney$gen_celltype,kidney$diet)
table(liver$gen_celltype,liver$diet)
table(vis$gen_celltype,vis$diet)
table(sc$gen_celltype,sc$diet)

# Min Artery	73
# Min Vein	192
# Min Cap	1057

set.seed(1234)

# ARTERY

Idents(brain) <- "gen.diet"
Idents(heart) <- "gen.diet"
Idents(lung) <- "gen.diet"
Idents(liver) <- "gen.diet"
Idents(kidney) <- "gen.diet"
Idents(vis) <- "gen.diet"
Idents(sc) <- "gen.diet"

kidney <- RenameIdents(kidney, `EC-venule_Western` = "EC-ven_Western",`EC-venule_chow` = "EC-ven_chow",
                       `EC-Aqp1-arteriole_Western` = "EC-art_Western",`EC-Aqp1-arteriole_chow` = "EC-art_chow")
vis <- RenameIdents(vis, `EC-venule_Western` = "EC-ven_Western",`EC-venule_chow` = "EC-ven_chow")
sc <- RenameIdents(sc, `EC-venule_Western` = "EC-ven_Western",`EC-venule_chow` = "EC-ven_chow")                  
                     
# Downsample the clusters to a maximum of 73 cells each for ARTERY
brain.art <- subset(brain, downsample = 73)
heart.art <- subset(heart, downsample = 73)
lung.art <- subset(lung, downsample = 73)
liver.art <- subset(liver, downsample = 73)
kidney.art <- subset(kidney, downsample = 73)
vis.art <- subset(vis, downsample = 73)
sc.art <- subset(sc, downsample = 73)

names = c("brain","heart","lung","liver","kidney","vis","sc")
k=0
# find DEG for "EC-art"
for (i in c(brain.art,heart.art,lung.art,liver.art,kidney.art,vis.art,sc.art)) {
  k=k+1
  print(names[k])
  #Idents(i) <- "gen.diet"
  resp <- FindMarkers(i, ident.1 = "EC-art_Western", 
                         ident.2 = "EC-art_chow", 
                         logfc.threshold = 0,
                         verbose = T) #, min.pct = 0.25
  write.table(resp, file = paste("DS_",names[k],"_EC-art_WD_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  print(nrow(resp[resp$p_val_adj<0.05,]))
}

# cap
k=0
# find DEG for "EC-cap"
for (i in c(brain.art,heart.art,lung.art,liver.art,kidney.art,vis.art,sc.art)) {
  k=k+1
  print(names[k])
  #Idents(i) <- "gen.diet"
  resp <- FindMarkers(i, ident.1 = "EC-cap_Western", 
                      ident.2 = "EC-cap_chow", 
                      logfc.threshold = 0,
                      verbose = T) #, min.pct = 0.25
  write.table(resp, file = paste("DS_73",names[k],"_EC-cap_WD_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  print(nrow(resp[resp$p_val_adj<0.05,]))
}

# ven
k=0
# find DEG for "EC-ven"
for (i in c(brain.art,heart.art,lung.art,liver.art,kidney.art,vis.art,sc.art)) {
  k=k+1
  print(names[k])
  #Idents(i) <- "gen.diet"
  resp <- FindMarkers(i, ident.1 = "EC-ven_Western", 
                      ident.2 = "EC-ven_chow", 
                      logfc.threshold = 0,
                      verbose = T) #, min.pct = 0.25
  write.table(resp, file = paste("DS_73",names[k],"_EC-ven_WD_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  print(nrow(resp[resp$p_val_adj<0.05,]))
}

# CAPILLARY

# Downsample the clusters to a maximum of 1057 cells each for cap
brain.cap <- subset(brain, downsample = 1057)
heart.cap <- subset(heart, downsample = 1057)
lung.cap <- subset(lung, downsample = 1057)
liver.cap <- subset(liver, downsample = 1057)
kidney.cap <- subset(kidney, downsample = 1057)
vis.cap <- subset(vis, downsample = 1057)
sc.cap <- subset(sc, downsample = 1057)

names = c("brain","heart","lung","liver","kidney","vis","sc")
k=0
# find DEG for "EC-cap"
for (i in c(brain.cap,heart.cap,lung.cap,liver.cap,kidney.cap,vis.cap,sc.cap)) {
  k=k+1
  print(names[k])
  Idents(i) <- "gen.diet"
  resp <- FindMarkers(i, ident.1 = "EC-cap_Western", 
                      ident.2 = "EC-cap_chow", 
                      logfc.threshold = 0,
                      verbose = T) #, min.pct = 0.25
  write.table(resp, file = paste("DS_",names[k],"_EC-cap_WD_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  print(nrow(resp[resp$p_val_adj<0.05,]))
}

# VEIN

# Downsample the clusters to a maximum of 192 cells each for vein
brain.ven <- subset(brain, downsample = 192)
heart.ven <- subset(heart, downsample = 192)
lung.ven <- subset(lung, downsample = 192)
liver.ven <- subset(liver, downsample = 192)
kidney.ven <- subset(kidney, downsample = 192)
vis.ven <- subset(vis, downsample = 192)
sc.ven <- subset(sc, downsample = 192)

names = c("brain","heart","lung","liver","kidney","vis","sc")
k=0
# find DEG for "EC-ven"
for (i in c(brain.ven,heart.ven,lung.ven,liver.ven,kidney.ven,vis.ven,sc.ven)) {
  k=k+1
  print(names[k])
  #Idents(i) <- "gen.diet"
  resp <- FindMarkers(i, ident.1 = "EC-ven_Western", 
                      ident.2 = "EC-ven_chow", 
                      logfc.threshold = 0,
                      verbose = T) #, min.pct = 0.25
  write.table(resp, file = paste("DS_",names[k],"_EC-ven_WD_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
  print(nrow(resp[resp$p_val_adj<0.05,]))
}

####################

# DOWNSAMPLE clusters per organ for DEG

set.seed(1234)

Idents(brain) <- "clust.diet"
Idents(lung) <- "clust.diet"
Idents(kidney) <- "clust.diet"
Idents(sc) <- "clust.diet"

heart$clust.diet <- paste(heart$heart_celltype, heart$diet, sep = "_")
Idents(heart) <- "clust.diet"
liver$clust.diet <- paste(liver$liver_celltype, liver$diet, sep = "_")
Idents(liver) <- "clust.diet"
vis$clust.diet <- paste(vis$vis_celltype, vis$diet, sep = "_")
Idents(vis) <- "clust.diet"

# Downsample the clusters to a maximum of 100 cells each for cluster
brain.art <- subset(brain, downsample = 100)
heart.art <- subset(heart, downsample = 100)
lung.art <- subset(lung, downsample = 100)
liver.art <- subset(liver, downsample = 100)
kidney.art <- subset(kidney, downsample = 100)
vis.art <- subset(vis, downsample = 100)
sc.art <- subset(sc, downsample = 100)

dir.path = "/Users/olga.bondareva/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/organs-3rd-processing/downsample/organs_clusters/"
names = c("brain","heart","lung","liver","kidney","vis","sc")
names = c("heart","liver","vis")
k=0
# find DEG for "EC-art"
for (i in c(heart.art,liver.art,vis.art)) { # brain.art,heart.art,lung.art,liver.art,kidney.art,vis.art,sc.art
  k=k+1
  print(names[k])
  
  i@active.ident <- factor(x = i@active.ident, levels = sort(levels(Idents(i))))
  levels(Idents(i))
  
  lev1 <- levels(i@active.ident)[seq(2,(length(levels(Idents(i)))),2)] # condition 1=HFD
  lev2 <- levels(i@active.ident)[seq(1,(length(levels(Idents(i)))),2)] # condition 2=chow
  
  # for all the clusters
  for (j in seq(1,length(lev1),1)) { # 16 for not-filtered
    response2 = NULL
    print(lev1[j])
    print(lev2[j])
    # TG 19m vs WT 19m
    response2 <- FindMarkers(object = i, ident.1 = lev1[j], 
                             ident.2 = lev2[j], 
                             logfc.threshold = 0,
                             verbose = T) #, min.pct = 0.25
    write.table(response2, 
    file = paste(dir.path,"DS_",names[k],"_",lev1[j],"_vs_chow.txt", sep = "", collapse = NULL), sep = "\t")
    print(nrow(response2[response2$p_val<0.05,]))
    }
}


