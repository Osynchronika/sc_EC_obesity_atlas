#library(Seurat)
library(dplyr)
library(Matrix)
#library(igraph)
#library(iTALK)
require(gplots)
require(ggplot2)
library(RColorBrewer)
library(cowplot)
#library(biomaRt)
require(scales)
#library(monocle)
library(Rfast)


# load data
setwd("~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_06_09_M_HFD_Rev_3m/analysis/recovery")
setwd("~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_06_09_M_HFD_Rev_3m/analysis/NEW_SEQ/recovery-fat")
file.dir = "~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_06_09_M_HFD_Rev_3m/analysis/organs/"
file.dir = "~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_06_09_M_HFD_Rev_3m/analysis/NEW_SEQ/organs/"
# read in the DEG for ven-art-capillaries from all organs
organ = c("brain","lung", "heart", "liver", "kidney", "vis", "sc")

# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(organ))

tab.deg.cap = list()
tab.deg.art = list()
tab.deg.ven = list()

###
# all organs
for (i in seq(1,length(organ),1)) {
  tab.deg.cap[[i]] = read.table(paste0(file.dir,organ[i],"/deg-gen/table_deg_",organ[i],"_Cap_all_TP.txt"), sep='\t')
  tab.deg.art[[i]] = read.table(paste0(file.dir,organ[i],"/deg-gen/table_deg_",organ[i],"_art_all_TP.txt"), sep='\t')
  tab.deg.ven[[i]] = read.table(paste0(file.dir,organ[i],"/deg-gen/table_deg_",organ[i],"_ven_all_TP.txt"), sep='\t')
}
# fat new seq
for (i in seq(1,length(organ),1)) {
  tab.deg.cap[[i]] = read.table(paste0(file.dir,organ[i],"/deg-gen/table_deg_",organ[i],"_Cap_all_TP.txt"), sep='\t')
  tab.deg.art[[i]] = read.table(paste0(file.dir,organ[i],"/deg-gen/table_deg_",organ[i],"_art_all_TP.txt"), sep='\t')
  tab.deg.ven[[i]] = read.table(paste0(file.dir,organ[i],"/deg-gen/table_deg_",organ[i],"_ven_all_TP.txt"), sep='\t')
}



######
# 1. organs, sort by 6m time point
####
# 1. Capillaries

tab.deg.cap.up = list()
tab.deg.cap.down = list()

for (i in seq(1,length(organ),1)) {
  n = NULL
  cutoff = NULL
  fc = NULL
  # sort only 6m signif genes
  tab.deg.cap[[i]] = tab.deg.cap[[i]][tab.deg.cap[[i]]$WD_6m_Q.val<0.05,]
  # exclude at 6m  -0.1<logFC<0.1 genes
  n = rownames(tab.deg.cap[[i]][(tab.deg.cap[[i]]$WD_6m_logFC>=0.1)|(tab.deg.cap[[i]]$WD_6m_logFC<=(-0.1)),])
  # rm NAs
  n = n[-grep("NA",n)]

  tab.deg.cap[[i]] = tab.deg.cap[[i]][n,]
  # sort descending logFC
  tab.deg.cap[[i]] = tab.deg.cap[[i]][order(tab.deg.cap[[i]]$WD_6m_logFC, decreasing = T),]
  print(organ[i])
  print(nrow(tab.deg.cap[[i]]))
  
  # find recovered genes (50%)
  # UP
  tab.deg.cap.up[[i]] = tab.deg.cap[[i]][(tab.deg.cap[[i]]$WD_6m_logFC>0),]
  cutoff =tab.deg.cap.up[[i]]$WD_6m_logFC/2
  fc = (tab.deg.cap.up[[i]]$Rev_3m_logFC<=cutoff)
  
  # recovered 
  recov = tab.deg.cap.up[[i]][fc,]
  # rm NAs
  if (length(-grep("NA",rownames(recov)))!=0) {
    recov = recov[-grep("NA",rownames(recov)),]
  }
  
  write.table(recov, paste0(organ[i],"_cap_6m_UP_recovered.txt"), sep = '\t')
  
  # non-recovered
  nonrecov = tab.deg.cap.up[[i]][!fc,]
  # rm NAs
  if (length(-grep("NA",rownames(nonrecov)))!=0) {
    nonrecov = nonrecov[-grep("NA",rownames(nonrecov)),]
  }
  write.table(nonrecov, paste0(organ[i],"_cap_6m_UP_non-recovered.txt"), sep = '\t')
  
  # heatmaps
  if (nrow(recov)>1) {
  # recovered
  b <- recov[,c(1,3,5,7,9)]
  # make logFC max 0.5
  b[b>0.5] = 0.5
  b[b<(-0.5)] = -0.5
  # NA = 0
  b[is.na(b)] = 0
  
  mi <- min(b,na.rm=TRUE) #l
  mi
  ma <- max(b,na.rm=TRUE) #l
  ma
  ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
  ColorLevels <- seq(mi, ma, length=length(ColorRamp))
  
  pdf(paste0("heatmap-",organ[i],"-cap-UP-recov.pdf"),width=10,height=15,paper='special')
  heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
            main = paste0(organ[i]," cap UP recov DEG, logFC(x/chow)"), 
            Rowv= T,Colv= F, 
            ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(5,10),     # widens margins around plot
            keysize=0.75,
            dendrogram="row", cexRow=0.8, cexCol=0.8) 
  dev.off()
  }
  
  if (nrow(nonrecov)>1) {
  # non-recovered
  b <- nonrecov[,c(1,3,5,7,9)]
  # make logFC max 0.5
  b[b>0.5] = 0.5
  b[b<(-0.5)] = -0.5
  # NA = 0
  b[is.na(b)] = 0
  
  mi <- min(b,na.rm=TRUE) #l
  mi
  ma <- max(b,na.rm=TRUE) #l
  ma
  ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
  ColorLevels <- seq(mi, ma, length=length(ColorRamp))
  
  pdf(paste0("heatmap-",organ[i],"-cap-UP-non-recov.pdf"),width=10,height=15,paper='special')
  heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
            main = paste0(organ[i]," cap UP non-recov DEG, logFC(x/chow)"), 
            Rowv= T,Colv= F, 
            ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(5,10),     # widens margins around plot
            keysize=0.75,
            dendrogram="row", cexRow=0.8, cexCol=0.8) 
  dev.off()
  }
  
  recov = NULL
  
  # DOWN
  tab.deg.cap.down[[i]] = tab.deg.cap[[i]][(tab.deg.cap[[i]]$WD_6m_logFC<0),]
  cutoff =tab.deg.cap.down[[i]]$WD_6m_logFC/2
  fc = (tab.deg.cap.down[[i]]$Rev_3m_logFC>=cutoff)
  
  # recovered 
  recov = tab.deg.cap.down[[i]][fc,]
  # rm NAs
  if (length(-grep("NA",rownames(recov)))!=0) {
    recov = recov[-grep("NA",rownames(recov)),]
  }
  
  write.table(recov, paste0(organ[i],"_cap_6m_DOWN_recovered.txt"), sep = '\t')
  
  # non-recovered
  nonrecov = tab.deg.cap.down[[i]][!fc,]
  # rm NAs
  if (length(-grep("NA",rownames(nonrecov)))!=0) {
    nonrecov = nonrecov[-grep("NA",rownames(nonrecov)),]
  }
  write.table(nonrecov, paste0(organ[i],"_cap_6m_DOWN_non-recovered.txt"), sep = '\t')
  
  # heatmaps
  if (nrow(recov)>1) {# recovered
  b <- recov[,c(1,3,5,7,9)]
  # make logFC max 0.5
  b[b>0.5] = 0.5
  b[b<(-0.5)] = -0.5
  # NA = 0
  b[is.na(b)] = 0
  
  mi <- min(b,na.rm=TRUE) #l
  mi
  ma <- max(b,na.rm=TRUE) #l
  ma
  ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
  ColorLevels <- seq(mi, ma, length=length(ColorRamp))
  
  pdf(paste0("heatmap-",organ[i],"-cap-DOWN-recov.pdf"),width=10,height=15,paper='special')
  heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
            main = paste0(organ[i]," cap DOWN recov DEG, logFC(x/chow)"), 
            Rowv= T,Colv= F, 
            ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(5,10),     # widens margins around plot
            keysize=0.75,
            dendrogram="row", cexRow=0.8, cexCol=0.8) 
  dev.off()
  }
  if (nrow(nonrecov)>1) {
  # non-recovered
  b <- nonrecov[,c(1,3,5,7,9)]
  # make logFC max 0.5
  b[b>0.5] = 0.5
  b[b<(-0.5)] = -0.5
  # NA = 0
  b[is.na(b)] = 0
  
  mi <- min(b,na.rm=TRUE) #l
  mi
  ma <- max(b,na.rm=TRUE) #l
  ma
  ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
  ColorLevels <- seq(mi, ma, length=length(ColorRamp))
  
  pdf(paste0("heatmap-",organ[i],"-cap-DOWN-non-recov.pdf"),width=10,height=15,paper='special')
  heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
            main = paste0(organ[i]," cap DOWN non-recov DEG, logFC(x/chow)"), 
            Rowv= T,Colv= F, 
            ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(5,10),     # widens margins around plot
            keysize=0.75,
            dendrogram="row", cexRow=0.8, cexCol=0.8) 
  dev.off()
  }
}

###########
# liver
fatty.acid.trans = c("Fabp1","Fabp2","Fabp3","Fabp4","Fabp5","Fabp6","Fabp7","Cd36","Lpl","Slc27a1","Slc27a4","Dbi",
                     "Acsbg1","Acsbg2","Acsl1","Acsl2","Acsl3","Acsl4","Acsl5","Acsl6") # enzymes

b <- tab.deg.cap[[4]][c("Cd36","Fabp4"),c(1,3,5,7,9)]
b = b[complete.cases(b),]
# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5
# NA = 0
b[is.na(b)] = 0

mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma
ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))

heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
          main = paste0("liver cap recov DEG, logFC(x/chow)"), 
          Rowv= T,Colv= F, 
          ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="row", cexRow=0.8, cexCol=0.8) 


####
# 2. Artery

tab.deg.art.up = list()
tab.deg.art.down = list()

for (i in rev(seq(1,length(organ),1))) {
  # sort only 6m signif genes
  tab.deg.art[[i]] = tab.deg.art[[i]][tab.deg.art[[i]]$WD_6m_Q.val<0.05,]
  # exclude at 6m  -0.1<logFC<0.1 genes
  n = rownames(tab.deg.art[[i]][(tab.deg.art[[i]]$WD_6m_logFC>=0.1)|(tab.deg.art[[i]]$WD_6m_logFC<=(-0.1)),])
  # rm NAs
  n = n[-grep("NA",n)]
  
  tab.deg.art[[i]] = tab.deg.art[[i]][n,]
  # sort descending logFC
  tab.deg.art[[i]] = tab.deg.art[[i]][order(tab.deg.art[[i]]$WD_6m_logFC, decreasing = T),]
  print(organ[i])
  print(nrow(tab.deg.art[[i]]))
  # find recovered genes (50%)
  # UP
  tab.deg.art.up[[i]] = tab.deg.art[[i]][(tab.deg.art[[i]]$WD_6m_logFC>0),]
  cutoff =tab.deg.art.up[[i]]$WD_6m_logFC/2
  fc = (tab.deg.art.up[[i]]$Rev_3m_logFC<=cutoff)
  print(nrow(tab.deg.art.up[[i]]))
  # recovered 
  recov = tab.deg.art.up[[i]][fc,]
  # rm NAs
  if (length(-grep("NA",rownames(recov)))!=0) {
    recov = recov[-grep("NA",rownames(recov)),]
  }
  
  write.table(recov, paste0(organ[i],"_art_6m_UP_recovered.txt"), sep = '\t')
  
  # non-recovered
  nonrecov = tab.deg.art.up[[i]][!fc,]
  # rm NAs
  if (length(-grep("NA",rownames(nonrecov)))!=0) {
    nonrecov = nonrecov[-grep("NA",rownames(nonrecov)),]
  }
  
  write.table(nonrecov, paste0(organ[i],"_art_6m_UP_non-recovered.txt"), sep = '\t')
  
  # heatmaps
  if (nrow(recov)>1) {
    # recovered
    b <- recov[,c(1,3,5,7,9)]
    # make logFC max 0.5
    b[b>0.5] = 0.5
    b[b<(-0.5)] = -0.5
    # NA = 0
    b[is.na(b)] = 0
    
    mi <- min(b,na.rm=TRUE) #l
    mi
    ma <- max(b,na.rm=TRUE) #l
    ma
    ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    
    pdf(paste0("heatmap-",organ[i],"-art-UP-recov.pdf"),width=10,height=15,paper='special')
    heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
              main = paste0(organ[i]," art UP recov DEG, logFC(x/chow)"), 
              Rowv= T,Colv= F, 
              ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(5,10),     # widens margins around plot
              keysize=0.75,
              dendrogram="row", cexRow=0.8, cexCol=0.8) 
    dev.off()
  }
  
  if (nrow(nonrecov)>1) {
    # non-recovered
    b <- nonrecov[,c(1,3,5,7,9)]
    # make logFC max 0.5
    b[b>0.5] = 0.5
    b[b<(-0.5)] = -0.5
    # NA = 0
    b[is.na(b)] = 0
    
    mi <- min(b,na.rm=TRUE) #l
    mi
    ma <- max(b,na.rm=TRUE) #l
    ma
    ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    
    pdf(paste0("heatmap-",organ[i],"-art-UP-non-recov.pdf"),width=10,height=15,paper='special')
    heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
              main = paste0(organ[i]," art UP non-recov DEG, logFC(x/chow)"), 
              Rowv= T,Colv= F, 
              ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(5,10),     # widens margins around plot
              keysize=0.75,
              dendrogram="row", cexRow=0.8, cexCol=0.8) 
    dev.off()
  }
  
  # DOWN
  tab.deg.art.down[[i]] = tab.deg.art[[i]][(tab.deg.art[[i]]$WD_6m_logFC<0),]
  cutoff =tab.deg.art.down[[i]]$WD_6m_logFC/2
  fc = (tab.deg.art.down[[i]]$Rev_3m_logFC>=cutoff)
  
  # recovered 
  recov = tab.deg.art.down[[i]][fc,]
  # rm NAs
  if (length(-grep("NA",rownames(recov)))!=0) {
    recov = recov[-grep("NA",rownames(recov)),]
  }
  
  write.table(recov, paste0(organ[i],"_art_6m_DOWN_recovered.txt"), sep = '\t')
  
  # non-recovered
  nonrecov = tab.deg.art.down[[i]][!fc,]
  # rm NAs
  if (length(-grep("NA",rownames(nonrecov)))!=0) {
  nonrecov = nonrecov[-grep("NA",rownames(nonrecov)),]
  }
  write.table(nonrecov, paste0(organ[i],"_art_6m_DOWN_non-recovered.txt"), sep = '\t')
  
  # heatmaps
  if (nrow(recov)>1) {# recovered
    b <- recov[,c(1,3,5,7,9)]
    # make logFC max 0.5
    b[b>0.5] = 0.5
    b[b<(-0.5)] = -0.5
    # NA = 0
    b[is.na(b)] = 0
    
    mi <- min(b,na.rm=TRUE) #l
    mi
    ma <- max(b,na.rm=TRUE) #l
    ma
    ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    
    pdf(paste0("heatmap-",organ[i],"-art-DOWN-recov.pdf"),width=10,height=15,paper='special')
    heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
              main = paste0(organ[i]," art DOWN recov DEG, logFC(x/chow)"), 
              Rowv= T,Colv= F, 
              ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(5,10),     # widens margins around plot
              keysize=0.75,
              dendrogram="row", cexRow=0.8, cexCol=0.8) 
    dev.off()
  }
  if (nrow(nonrecov)>1) {
    # non-recovered
    b <- nonrecov[,c(1,3,5,7,9)]
    # make logFC max 0.5
    b[b>0.5] = 0.5
    b[b<(-0.5)] = -0.5
    # NA = 0
    b[is.na(b)] = 0
    
    mi <- min(b,na.rm=TRUE) #l
    mi
    ma <- max(b,na.rm=TRUE) #l
    ma
    ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    
    pdf(paste0("heatmap-",organ[i],"-art-DOWN-non-recov.pdf"),width=10,height=15,paper='special')
    heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
              main = paste0(organ[i]," art DOWN non-recov DEG, logFC(x/chow)"), 
              Rowv= T,Colv= F, 
              ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(5,10),     # widens margins around plot
              keysize=0.75,
              dendrogram="row", cexRow=0.8, cexCol=0.8) 
    dev.off()
  }
}

####
# 3. Vein

tab.deg.ven.up = list()
tab.deg.ven.down = list()

for (i in rev(seq(1,length(organ),1))) {
  # sort only 6m signif genes
  tab.deg.ven[[i]] = tab.deg.ven[[i]][tab.deg.ven[[i]]$WD_6m_Q.val<0.05,]
  # exclude at 6m  -0.1<logFC<0.1 genes
  n = rownames(tab.deg.ven[[i]][(tab.deg.ven[[i]]$WD_6m_logFC>=0.1)|(tab.deg.ven[[i]]$WD_6m_logFC<=(-0.1)),])
  # rm NAs
  n = n[-grep("NA",n)]
  
  tab.deg.ven[[i]] = tab.deg.ven[[i]][n,]
  # sort descending logFC
  tab.deg.ven[[i]] = tab.deg.ven[[i]][order(tab.deg.ven[[i]]$WD_6m_logFC, decreasing = T),]
  print(organ[i])
  print(nrow(tab.deg.ven[[i]]))
  
  # find recovered genes (50%)
  # UP
  tab.deg.ven.up[[i]] = tab.deg.ven[[i]][(tab.deg.ven[[i]]$WD_6m_logFC>0),]
  cutoff =tab.deg.ven.up[[i]]$WD_6m_logFC/2
  fc = (tab.deg.ven.up[[i]]$Rev_3m_logFC<=cutoff)
  
  # recovered 
  recov = tab.deg.ven.up[[i]][fc,]
  # rm NAs
  if (length(-grep("NA",rownames(recov)))!=0) {
    recov = recov[-grep("NA",rownames(recov)),]
  }
  
  write.table(recov, paste0(organ[i],"_ven_6m_UP_recovered.txt"), sep = '\t')
  
  # non-recovered
  nonrecov = tab.deg.ven.up[[i]][!fc,]
  # rm NAs
  if (length(-grep("NA",rownames(nonrecov)))!=0) {
    nonrecov = nonrecov[-grep("NA",rownames(nonrecov)),]
  }
  write.table(nonrecov, paste0(organ[i],"_ven_6m_UP_non-recovered.txt"), sep = '\t')
  
  # heatmaps
  if (nrow(recov)>1) {
    # recovered
    b <- recov[,c(1,3,5,7,9)]
    # make logFC max 0.5
    b[b>0.5] = 0.5
    b[b<(-0.5)] = -0.5
    # NA = 0
    b[is.na(b)] = 0
    
    mi <- min(b,na.rm=TRUE) #l
    mi
    ma <- max(b,na.rm=TRUE) #l
    ma
    ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    
    pdf(paste0("heatmap-",organ[i],"-ven-UP-recov.pdf"),width=10,height=15,paper='special')
    heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
              main = paste0(organ[i]," ven UP recov DEG, logFC(x/chow)"), 
              Rowv= T,Colv= F, 
              ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(5,10),     # widens margins around plot
              keysize=0.75,
              dendrogram="row", cexRow=0.8, cexCol=0.8) 
    dev.off()
  }
  
  if (nrow(nonrecov)>1) {
    # non-recovered
    b <- nonrecov[,c(1,3,5,7,9)]
    # make logFC max 0.5
    b[b>0.5] = 0.5
    b[b<(-0.5)] = -0.5
    # NA = 0
    b[is.na(b)] = 0
    
    mi <- min(b,na.rm=TRUE) #l
    mi
    ma <- max(b,na.rm=TRUE) #l
    ma
    ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    
    pdf(paste0("heatmap-",organ[i],"-ven-UP-non-recov.pdf"),width=10,height=15,paper='special')
    heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
              main = paste0(organ[i]," ven UP non-recov DEG, logFC(x/chow)"), 
              Rowv= T,Colv= F, 
              ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(5,10),     # widens margins around plot
              keysize=0.75,
              dendrogram="row", cexRow=0.8, cexCol=0.8) 
    dev.off()
  }
  
  # DOWN
  tab.deg.ven.down[[i]] = tab.deg.ven[[i]][(tab.deg.ven[[i]]$WD_6m_logFC<0),]
  cutoff =tab.deg.ven.down[[i]]$WD_6m_logFC/2
  fc = (tab.deg.ven.down[[i]]$Rev_3m_logFC>=cutoff)
  
  # recovered 
  recov = tab.deg.ven.down[[i]][fc,]
  # rm NAs
  if (length(-grep("NA",rownames(recov)))!=0) {
    recov = recov[-grep("NA",rownames(recov)),]
  }
  
  write.table(recov, paste0(organ[i],"_ven_6m_DOWN_recovered.txt"), sep = '\t')
  
  # non-recovered
  nonrecov = tab.deg.ven.down[[i]][!fc,]
  # rm NAs
  if (length(-grep("NA",rownames(nonrecov)))!=0) {
    nonrecov = nonrecov[-grep("NA",rownames(nonrecov)),]
  }
  write.table(nonrecov, paste0(organ[i],"_ven_6m_DOWN_non-recovered.txt"), sep = '\t')
  
  # heatmaps
  if (nrow(recov)>1) {# recovered
    b <- recov[,c(1,3,5,7,9)]
    # make logFC max 0.5
    b[b>0.5] = 0.5
    b[b<(-0.5)] = -0.5
    # NA = 0
    b[is.na(b)] = 0
    
    mi <- min(b,na.rm=TRUE) #l
    mi
    ma <- max(b,na.rm=TRUE) #l
    ma
    ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    
    pdf(paste0("heatmap-",organ[i],"-ven-DOWN-recov.pdf"),width=10,height=15,paper='special')
    heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
              main = paste0(organ[i]," ven DOWN recov DEG, logFC(x/chow)"), 
              Rowv= T,Colv= F, 
              ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(5,10),     # widens margins around plot
              keysize=0.75,
              dendrogram="row", cexRow=0.8, cexCol=0.8) 
    dev.off()
  }
  if (nrow(nonrecov)>1) {
    # non-recovered
    b <- nonrecov[,c(1,3,5,7,9)]
    # make logFC max 0.5
    b[b>0.5] = 0.5
    b[b<(-0.5)] = -0.5
    # NA = 0
    b[is.na(b)] = 0
    
    mi <- min(b,na.rm=TRUE) #l
    mi
    ma <- max(b,na.rm=TRUE) #l
    ma
    ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    
    pdf(paste0("heatmap-",organ[i],"-ven-DOWN-non-recov.pdf"),width=10,height=15,paper='special')
    heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
              main = paste0(organ[i]," ven DOWN non-recov DEG, logFC(x/chow)"), 
              Rowv= T,Colv= F, 
              ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), # rep(cell_col.both, each=2)
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(5,10),     # widens margins around plot
              keysize=0.75,
              dendrogram="row", cexRow=0.8, cexCol=0.8) 
    dev.off()
  }
}




for (i in seq(1,length(organ),1)) {
b <- tab.deg.cap[[i]][,c(1,3,5)] # 
b = b[rownames(b)%in%c("Klf2","Klf4","Klf6","Nfkbia","Meox2","Tcf15"),]
# make logFC max 0.5
#b[b>0.5] = 0.5
#b[b<(-0.5)] = -0.5
# NA = 0
b[is.na(b)] = 0

mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma
ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))

pdf(paste0("heatmap-",organ[i],"-cap-KLFs.pdf"),width=10,height=10,paper='special')
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
          main = paste0(organ[i]," cap KLFs, logFC(x/chow)"), 
          Rowv= F,Colv= F, 
          ColSideColors= c("#ffe0a1","#fc8b00","#f0e02e"), #,"#fc5000","#a2f02e"
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(10,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="none", cexRow=0.8, cexCol=0.8) 
dev.off()
}


####
electron = read.table("recovery/electron_transp.txt")

for (i in seq(1,length(organ),1)) {
  b <- tab.deg.cap[[i]][,c(1,3,5,7,9)] # 
  b = b[rownames(b)%in%electron$V1,]
  # make logFC max 0.5
  #b[b>0.5] = 0.5
  #b[b<(-0.5)] = -0.5
  # NA = 0
  b[is.na(b)] = 0
  
  mi <- min(b,na.rm=TRUE) #l
  mi
  ma <- max(b,na.rm=TRUE) #l
  ma
  ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
  ColorLevels <- seq(mi, ma, length=length(ColorRamp))
  
  pdf(paste0("heatmap-",organ[i],"-cap-ETC.pdf"),width=10,height=15,paper='special')
  heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
            main = paste0(organ[i]," cap electron transp chain, logFC(x/chow)"), 
            Rowv= T,Colv= F, 
            ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), #,"#fc5000","#a2f02e"
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(10,10),     # widens margins around plot
            keysize=0.75,
            dendrogram="row", cexRow=0.8, cexCol=0.8) 
  dev.off()
}

mt.genes = c("mt-Nd1", "mt-Nd2", "mt-Co2","mt-Atp6", "mt-Co3", 
             "mt-Co1","mt-Cytb","mt-Nd4","mt-Nd3","mt-Nd4l","mt-Atp8","mt-Nd6")

resp.chain = read.table("resp_chain.txt", sep='\t')
compl.col = c("cyan2","coral","deeppink1","darkturquoise","red")
names(compl.col) = unique(resp.chain$V2)
resp.chain$V3 = 1
resp.chain$V3 = compl.col[match(resp.chain$V2, names(compl.col))]

###########
organ = c("brain","lung", "heart", "liver", "kidney")
tab.deg.cap = list()
for (i in seq(1,length(organ),1)) {
  tab.deg.cap[[i]] = read.table(paste0(file.dir,organ[i],"/deg-gen/table_deg_",organ[i],"_Cap_all_TP.txt"), sep='\t')
  b <- tab.deg.cap[[i]][,c(1,3,5,7,9)] # 
  b = b[rownames(b)%in%klfs,] #resp.chain$V1
  b <- b[order(match(rownames(b),kk)),]
  
  # replase NA with 0
  b[is.na(b)] = 0
  # exclude -0.1<logFC<0.1 genes
  #b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]
  
  # make logFC max 0.5
  b[b>0.5] = 0.5
  b[b<(-0.5)] = -0.5
  
  mi <- min(b,na.rm=TRUE) #l
  ma <- max(b,na.rm=TRUE) #l
  
  ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
  ColorLevels <- seq(mi, ma, length=length(ColorRamp))
  pdf(paste0("heatmap-",organ[i],"-cap-kk.pdf"),width=10,height=15,paper='special')
  heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = paste0(organ[i],", kk in cap, logFC(WD/chow)"), #PPAR signaling
            Rowv= F,Colv= F, 
            ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"),
            #RowSideColors= resp.chain$V3,
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(5,10),     # widens margins around plot
            keysize=0.75,
            dendrogram="none", cexRow=0.8, cexCol=0.8) 
  dev.off()
}

tab.deg.art = list()
for (i in seq(1,length(organ),1)) {
  tab.deg.art[[i]] = read.table(paste0(file.dir,organ[i],"/deg-gen/table_deg_",organ[i],"_art_all_TP.txt"), sep='\t')
  b <- tab.deg.art[[i]][,c(1,3,5,7,9)] # 
  b = b[rownames(b)%in%hsp,] #resp.chain$V1
  b <- b[order(match(rownames(b),hsp)),]
  
  # replase NA with 0
  b[is.na(b)] = 0
  # exclude -0.1<logFC<0.1 genes
  #b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]
  
  # make logFC max 0.5
  b[b>0.5] = 0.5
  b[b<(-0.5)] = -0.5
  
  mi <- min(b,na.rm=TRUE) #l
  ma <- max(b,na.rm=TRUE) #l
  
  ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
  ColorLevels <- seq(mi, ma, length=length(ColorRamp))
  pdf(paste0("heatmap-",organ[i],"-art-hsp.pdf"),width=10,height=15,paper='special')
  heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = paste0(organ[i],", hsp in art, logFC(WD/chow)"), #PPAR signaling
            Rowv= F,Colv= F, 
            ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"),
            #RowSideColors= resp.chain$V3,
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(5,10),     # widens margins around plot
            keysize=0.75,
            dendrogram="none", cexRow=0.8, cexCol=0.8) 
  dev.off()
}

tab.deg.ven = list()
for (i in seq(1,length(organ),1)) {
  tab.deg.ven[[i]] = read.table(paste0(file.dir,organ[i],"/deg-gen/table_deg_",organ[i],"_ven_all_TP.txt"), sep='\t')
  b <- tab.deg.ven[[i]][,c(1,3,5,7,9)] # 
  b = b[rownames(b)%in%hsp,] #resp.chain$V1
  b <- b[order(match(rownames(b),hsp)),]
  
  # replase NA with 0
  b[is.na(b)] = 0
  # exclude -0.1<logFC<0.1 genes
  #b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]
  
  # make logFC max 0.5
  b[b>0.5] = 0.5
  b[b<(-0.5)] = -0.5
  
  mi <- min(b,na.rm=TRUE) #l
  ma <- max(b,na.rm=TRUE) #l
  
  ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
  ColorLevels <- seq(mi, ma, length=length(ColorRamp))
  pdf(paste0("heatmap-",organ[i],"-ven-hsp.pdf"),width=10,height=15,paper='special')
  heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = paste0(organ[i],", hsp in ven, logFC(WD/chow)"), #PPAR signaling
            Rowv= F,Colv= F, 
            ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"),
            #RowSideColors= resp.chain$V3,
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(5,10),     # widens margins around plot
            keysize=0.75,
            dendrogram="none", cexRow=0.8, cexCol=0.8) 
  dev.off()
}

### lung
hsp = c("Hsp90aa1","Hsp90ab1","Hspa1a","Hspa1b","Hspb1","Hspd1","Hspe1","Hsph1","Hspa8")

tab.deg.cap = read.table(paste0(file.dir,"lung/deg-gen/table_deg_lung_Cap_all_TP.txt"), sep='\t')
tab.deg.art = read.table(paste0(file.dir,"lung/deg-gen/table_deg_lung_art_all_TP.txt"), sep='\t')
tab.deg.ven = read.table(paste0(file.dir,"lung/deg-gen/table_deg_lung_ven_all_TP.txt"), sep='\t')


tab.deg.cap = read.table(paste0(file.dir,"lung/deg-clust/table_deg_lung_pulm1_all_TP.txt"), sep='\t')

a <- tab.deg.art[,c(1,3,5,7,9)] # mt.genes
a = a[rownames(a)%in%hsp,]  

b <- tab.deg.cap[,c(1,3,5,7,9)] # hsp
b = b[rownames(b)%in%resp.chain$V1,] 

d <- tab.deg.ven[,c(1,3,5,7,9)] # 
d = d[rownames(d)%in%hsp,] 

b = cbind(a,b,d)
# replase NA with 0
b[is.na(b)] = 0

b <- b[order(match(rownames(b),resp.chain$V1)),]

# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5

mi <- min(b,na.rm=TRUE) #l
ma <- max(b,na.rm=TRUE) #l

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
pdf(paste0("heatmap-lung-Pulm1-respchain-genes.pdf"),width=15,height=20,paper='special')
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = paste0("lung, Pulm1-resp chain, logFC(WD/chow)"), #PPAR signaling
          Rowv= F,Colv= F, 
          ColSideColors= rep(c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), times =1),
          #RowSideColors= resp.chain$V3,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="none", cexRow=0.5, cexCol=0.8) 
dev.off()

# Look for inflammation signature
nfkb = read.table("../NFkB.txt", sep='\t', header = T)
nfkb = nfkb$Symbol
nfkb = nfkb[-67]
library(stringr)
nfkb = str_to_sentence(nfkb, locale = "en")

inflamm = read.table("../GO_0006954_inflammatory_response.txt", sep='\t')
inflamm = inflamm$V1

FA.oxid = read.table("../GO_0019395_FA_oxidation.txt", sep='\t')
FA.oxid = unique(FA.oxid$V1)
FA.trans = read.table("../GO_0015908_FA_transport.txt", sep='\t')
FA.trans = unique(FA.trans$V1)

ecm = read.table("/Users/olga.bondareva/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/GO_0031012_ECM.txt", sep='\t')
ecm = unique(ecm$V1)
platelet = read.table("../GO_0070527_platelet_aggregation.txt", sep='\t')
platelet = unique(platelet$V1)
platelet.act = read.table("../GO_0030168_platelet_activation.txt", sep='\t')
platelet.act = unique(platelet.act$V1)
hemostasis = read.table("../GO_0007599_hemostasis.txt", sep='\t')
hemostasis = unique(hemostasis$V1)
calcium = read.table("../GO_0019722_calcium_mediated_signal.txt", sep="\t")
calcium = unique(calcium$V1)
ca2.ec = c("Anxa2","Mylk","Camk2d","F2r","Trpc6","Icam1","Vcam1","Stim1")
# coagulation cascade
coagul = c("F3","F7","F8","F9","F10","F11","F12","Vwf","F2","F2r","F2rl2","F2rl3",
           "C3","C5","A2m","Serpina1","Serpina5","Serpinf2",
           "Serpine1","Serpinb2","Fga","Fgb","Fgg", # pro-coag
           "Proc","Serpina5","Serpinc1","Serpind1","Tfpi","Plg","Plau","Plat") # anti-coag

fatty.acid.trans = c("Fabp1","Fabp2","Fabp3","Fabp4","Fabp5","Fabp6","Fabp7","Cd36","Lpl","Slc27a1","Slc27a4","Dbi",
                     "Acsbg1","Acsbg2","Acsl1","Acsl2","Acsl3","Acsl4","Acsl5","Acsl6") # enzymes
# coagulation cascade
coagul = c("F3","F7","F8","F9","F10","F11","F12","Vwf","F2","F2r","F2rl2","F2rl3",
           "C3","C5","A2m","Serpina1","Serpina5","Serpinf2",
           "Serpine1","Serpinb2","Fga","Fgb","Fgg", # pro-coag
           "Proc","Serpina5","Serpinc1","Serpind1","Tfpi","Plg","Plau","Plat") # anti-coag

leuk.recruit = c("Icam1","Icam2","Pecam1","Vcam1", "Sele", "Selp","Cd99","Cd47","Alcam","Hspg","Ptafr",
                 "Il1r","Tnfr1","Trap","Traf2","Traf6","Tirap","Tradd","Myd88","Irak1", "Irak4",
                 "Nfkb1","Nfbk2","Rela","Relb","Rel","Cox1","Cox2","Cxcl8","Cxcl10","Ccl26")

junctions = c("Cldn1","Cldn3","Cldn5","Cldn12","Tjp1","Tjp2","Tjp3","Ocln","Esam","F11r","Jam2","Jam3",
              "Nectin1","Nectin4","Afdn", # tight junctions
              "Gja1","Gja4","Gja5", # gap junctions
              "Jup","Ctnnb1","Ctnnd1","Cdh5","Pcdh12","Pecam1","Cd99", # adherent junctions
              "Icam1","Icam2","Pecam1","Vcam1", "Sele", "Selp","Cd99","Cd47","Alcam") # leukocyte adhesion
# fat Integrin signaling
integrin = c("Itgb1","Git2","Arpc1b","Pxn","Rhog","Ilk","Rhoc","Rhoa","Myl12b","Rhob","Arhgap10","Capns1","Rras",
             "Tspan6","Tspan7","Rapgef1","Fyn","Kras","Capn1","Rac1","Hras","Mapk3","Rap1a","Arpc3","Rock2",
             "Cav1","Myl12a")
# focal adhesion
focal = c("Itgb1","Rap1a","Rock2","Cav1","Itga1","Tnc","Rapgef1","Thbs1","Myl12a","Rhoa","Myl12b","Lama5","Cav2",
          "Lama4","Pxn","Pdgfb","Ilk","Parvb","Ccnd3","Col4a2","Col4a1","Itga6","Fyn","Rac1","Hras","Ppp1r12c","Mapk3")

# PPAR signaling
ppar = read.table("../KEGG_PPAR_signaling.txt", sep = '\t')
ppar = ppar$V2
ppar = c(ppar, "Abca1","Mgll","Acaa1a","Acaa1b","Apoa4")

#liver lipid digestion and absorbtion
lipid = c("Acaa1b","Abca1","Fabp1","Mgll", "Apoa2", "Apoa4","Cd36","HmgcS2", "Fabp4","Fabp5","Lpl","Dbi")

lipid = c("Fabp4","Cd36","Fabp5","Dbi","Lpl","Fabp1","Abca1")
# lipid mobilization
lip.mob = c("Saa1","Saa2","Saa4","Pon1","Pon3","Lcat","Ttr","Apoe", "Apoa1","Apoa2", "Apoa4","Apoc1","Apoc2",
            "Apob","Apod","Serpina1","Alb","Abca1")

# GO full focal adhesion
go.focal.adh = read.table("../GO_0005925_focal_adhesion.txt", sep='\t')
go.focal.adh = unique(go.focal.adh$V1)

# Kalucka 2018 Cell Metabolism
glycolysis = c("Pgm1","Pgm2","Gpi","Pfkl","Pfkm","Pfkp","Aldoa","Gapdh","Eno1","Pkm",
               "Pdha1","Pdha2","Pdhb","Dlat","Dld","Ldha","Ldhb")  #c("Foxk1", "Foxk2", "Hk1", "Gpi1", "Pfkl", "Pfkm", "Pfkp")
beta.oxid = c("Fabp4", "Fabp5","Acsl1","Acsl3","Acsl4","Acsl5","Cpt1a","Cpt1b","Cpt1c","Cpt2",
              "Acads","Acadm","Acadvl","Hadha","Hadhb","Acaa1", "Acaa2") #c("Cpt1", "Cpt2","Acadvl","Hadha","Hadhb","Acadm","Ehhadh", "Acat1", "Acaa1", "Acaa2")
oxid.PPP = c("G6phd","Pgls","H6pd","6pgdh")
non.oxid.PPP = c("Rpia","Rpe","Tkt")
oneC.metab = c("Dhfr","Mthfd1","Mthfd2")
serine = c("Phgdh","Psat1","Psph","Shmt1","Shmt2")
tca = c("Cs","Aco1","Aco2","Ogdh","Dlst","Dld","Sucla1","Suclg1","Suclg2","Sdha","Sdhb","Sdhc","Sdhd",
        "Fh","Mdh1","Mdh2","Acly","Acss2")
metabolism = c(beta.oxid,oxid.PPP, non.oxid.PPP,oneC.metab,serine,glycolysis,tca)



### for lungs
h2 = c("H2-Aa","H2-Ab1","H2-D1","H2-DMa","H2-DMb1","H2-Eb1","H2-K1","H2-Ke6","H2-M3","H2-Q4","H2-Q6","H2-Q7","H2-T22","H2-T23")
rps = c("Rps9","Rpl3","Rps8","Rpl11","Rpl7","Rpl7a","Rps15a","Rpl18a","Rps19","Rpl24","Rps20","Rpl26","Rpl28",
        "Rps21","Rps13")

klfs = c("Klf2","Klf3","Klf4","Klf6","Klf7","Klf9","Klf10","Klf12","Klf13","Klf16")
ap1 = c("Jun","Junb","Jund","Fos","Fosb","Egr1")

mt = rownames(tab.deg.cap)[grep("mt-",rownames(tab.deg.cap))]

slc.brain = c("Slc16a1","Slc22a8","Slc38a5","Slc2a1","Slc40a1","Slc1a1","Slc6a6","Slc39a10",
              "Slc7a1","Slc9a8","Slc25a32","Slc12a6","Slc39a7","Slc27a3","Slc7a5","Slc3a2",
              "Slc25a37","Slc25a30","Slc24a5","Slc43a2","Slc29a1","Slc52a2",
              "Slc25a24","Slc9a1","Slc6a8","Slc30a1","Slc4a7","Slc35f5")

g.br = c("Mt1","Mt2","Sgms1","Degs2","Id1","Id2","Jam2","Slc16a1")

cams = c("Vcam1","Pecam1","Alcam","Icam1")

resp.chain = read.table("../resp_chain.txt", sep='\t')
compl.col = c("cyan2","coral","deeppink1","darkturquoise","red")
names(compl.col) = unique(resp.chain$V2)
resp.chain$V3 = 1
resp.chain$V3 = compl.col[match(resp.chain$V2, names(compl.col))]



### brain
tab.deg.cap = read.table(paste0(file.dir,"brain/deg-gen/table_deg_brain_Cap_all_TP.txt"), sep='\t')
tab.deg.art = read.table(paste0(file.dir,"brain/deg-gen/table_deg_brain_art_all_TP.txt"), sep='\t')
tab.deg.ven = read.table(paste0(file.dir,"brain/deg-gen/table_deg_brain_ven_all_TP.txt"), sep='\t')

a <- tab.deg.art[,c(1,3,5,7,9)] # mt.genes
a = a[rownames(a)%in%cams,]  

b <- tab.deg.cap[,c(1,3,5,7,9)] # 
b = b[rownames(b)%in%g.br,] 

d <- tab.deg.ven[,c(1,3,5,7,9)] # 
d = d[rownames(d)%in%g.br,] 

b = cbind(a,b,d)
# replase NA with 0
b[is.na(b)] = 0

# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5

#b <- b[order(match(rownames(b),g.br)),]

mi <- min(b,na.rm=TRUE) #l
ma <- max(b,na.rm=TRUE) #l

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
pdf(paste0("heatmap-brain-art-cams.pdf"),width=15,height=10,paper='special')
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = paste0("brain, CAMs art, logFC(WD/chow)"), #PPAR signaling
          Rowv= T,Colv= F, 
          ColSideColors= rep(c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), times =3),
          #RowSideColors= resp.chain$V3,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="row", cexRow=0.8, cexCol=0.8) 
dev.off()

## fenestr
slc2 = c("Slc9a1","Slc29a1","Slc39a7","Slc12a6","Slc25a30","Slc6a8")
tjp = c("Gja1","Jam3","Pcdh12","Afdn","Tjp1","Jup","Tjp2")
tab.deg.cap = read.table(paste0(file.dir,"brain/deg-clust/table_deg_brain_fenestr_all_TP.txt"), sep='\t')

b <- tab.deg.cap[,c(1,3,5,7,9)] # 
b = b[rownames(b)%in%slc2,] 

# replase NA with 0
b[is.na(b)] = 0

# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5

#b <- b[order(match(rownames(b),g.br)),]
mi <- min(b,na.rm=TRUE) #l
ma <- max(b,na.rm=TRUE) #l

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
pdf(paste0("heatmap-brain-fenestr-slc.pdf"),width=15,height=10,paper='special')
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = paste0("brain, fenestr slcs, logFC(WD/chow)"), #PPAR signaling
          Rowv= T,Colv= F, 
          ColSideColors= rep(c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), times =1),
          #RowSideColors= resp.chain$V3,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="row", cexRow=0.8, cexCol=0.8) 
dev.off()

# kidney
tab.deg.cap = read.table(paste0(file.dir,"kidney/deg-gen/table_deg_kidney_Cap_all_TP.txt"), sep='\t')
tab.deg.art = read.table(paste0(file.dir,"kidney/deg-gen/table_deg_kidney_art_all_TP.txt"), sep='\t')
tab.deg.ven = read.table(paste0(file.dir,"kidney/deg-gen/table_deg_kidney_ven_all_TP.txt"), sep='\t')

tab.deg.cap = read.table(paste0(file.dir,"kidney/deg-clust/table_deg_kidney_glomeruli_all_TP.txt"), sep='\t')

a <- tab.deg.art[,c(1,3,5,7,9)] # mt.genes
a = a[rownames(a)%in%ap1,]  

b <- tab.deg.cap[,c(1,3,5,7,9)] # 
b = b[rownames(b)%in%ap1,] 

d <- tab.deg.ven[,c(1,3,5,7,9)] # 
d = d[rownames(d)%in%ap1,] 

b = cbind(a,b,d)
# replase NA with 0
b[is.na(b)] = 0

# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5

b <- b[order(match(rownames(b),c("Jund", "Junb", "Jun", "Fosb", "Fos", "Egr1"))),]

mi <- min(b,na.rm=TRUE) #l
ma <- max(b,na.rm=TRUE) #l

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
pdf(paste0("heatmap-kidney-glomer-Ap1.pdf"),width=15,height=10,paper='special')
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = paste0("kidney, AP1 glomer, logFC(WD/chow)"), #PPAR signaling
          Rowv= F,Colv= F, 
          ColSideColors= rep(c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), times =1),
          #RowSideColors= resp.chain$V3,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="none", cexRow=0.8, cexCol=0.8) 
dev.off()

# heart
tab.deg.cap = read.table(paste0(file.dir,"heart/deg-gen/table_deg_heart_Cap_all_TP.txt"), sep='\t')
tab.deg.art = read.table(paste0(file.dir,"heart/deg-gen/table_deg_heart_art_all_TP.txt"), sep='\t')
tab.deg.ven = read.table(paste0(file.dir,"heart/deg-gen/table_deg_heart_ven_all_TP.txt"), sep='\t')

a <- tab.deg.art[,c(1,3,5,7,9)] # mt.genes
a = a[rownames(a)%in%ap1,]  

b <- tab.deg.cap[,c(1,3,5,7,9)] # 
b = b[rownames(b)%in%ap1,] #klfs

d <- tab.deg.ven[,c(1,3,5,7,9)] # 
d = d[rownames(d)%in%klfs,] 

b = cbind(a,b,d)
# replase NA with 0
b[is.na(b)] = 0

# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5

#b <- b[order(match(rownames(b),g.br)),]

mi <- min(b,na.rm=TRUE) #l
ma <- max(b,na.rm=TRUE) #l

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
pdf(paste0("heatmap-heart-art-arteriole-AP1.pdf"),width=15,height=10,paper='special')
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = paste0("heart, AP1 art+arteriole, logFC(WD/chow)"), #PPAR signaling
          Rowv= T,Colv= F, 
          ColSideColors= rep(c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), times =1),
          #RowSideColors= resp.chain$V3,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="row", cexRow=0.8, cexCol=0.8) 
dev.off()

kk = c("Uba52","Rpl27","Gm10076","Ccdc85b","Plat","Trim37","Rhob","Grcc10","Gstp1",
       "Ccdc152","Dbp","Lars2","Rpl35","Rpl13a","Bloc1s1","Rpl36al","mt-Nd4l","Shank3","Tsc22d3",
       "Nr1d1","Ndufa3","Neurl3","H3f3b","Btg2","mt-Nd1","mt-Nd2","mt-Co2","mt-Atp6","mt-Co3")

# for fats
organ = c("vis","sc")
tab.deg.cap = list()
for (i in rev(seq(1,length(organ),1))) {
  tab.deg.cap[[i]] = read.table(paste0(file.dir,organ[i],"/deg-gen/table_deg_",organ[i],"_Cap_all_TP.txt"), sep='\t')
 # b <- tab.deg.cap[[i]][,c(1,3,7,9)] #only 4m time point
  b <- tab.deg.cap[[i]][,c(1,3,5,7,9)] #only 4m time point
  b = b[rownames(b)%in%focal,] #integrin focal ecm kk
  b <- b[order(match(rownames(b),focal)),]
  
  # replase NA with 0
  b[is.na(b)] = 0
  # exclude -0.1<logFC<0.1 genes
  #b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]
  
  # make logFC max 0.5
  b[b>0.5] = 0.5
  b[b<(-0.5)] = -0.5
  
  mi <- min(b,na.rm=TRUE) #l
  ma <- max(b,na.rm=TRUE) #l
  
  ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
  ColorLevels <- seq(mi, ma, length=length(ColorRamp))
  pdf(paste0("hm-",organ[i],"-cap-focal.pdf"),width=10,height=15,paper='special')
  heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = paste0(organ[i],", focal adh in cap, logFC(WD/chow)"), #PPAR signaling
            Rowv= T,Colv= F, 
            ColSideColors= c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), 
            #RowSideColors= resp.chain$V3,
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(5,10),     # widens margins around plot
            keysize=0.75,
            dendrogram="row", cexRow=0.8, cexCol=0.8) 
  dev.off()
  
  #### ACV
  head(tab.deg.art[[i]])
  a <- tab.deg.art[[i]][,c(1,3,5,7,9)] # 
  a = a[rownames(a)%in%lipid,] 
  
  head(tab.deg.cap[[i]])
  b <- tab.deg.cap[[i]][,c(1,3,5,7,9)] # 
  b = b[rownames(b)%in%lipid,] 
  
  head(tab.deg.ven[[i]])
  d <- tab.deg.ven[[i]][,c(1,3,5,7,9)] # 
  d = d[rownames(d)%in%lipid,] 
  
  b = cbind(a,b,d)
  b <- b[order(match(rownames(b),lipid)),]
  # replase NA with 0
  b[is.na(b)] = 0
  
  # make logFC max 0.5
  b[b>0.5] = 0.5
  b[b<(-0.5)] = -0.5
  
  mi <- min(b,na.rm=TRUE) #l
  ma <- max(b,na.rm=TRUE) #l
  
  ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
  ColorLevels <- seq(mi, ma, length=length(ColorRamp))
  pdf(paste0("hm-",organ[i],"-ACV-lipid.pdf"),width=15,height=7,paper='special')
  heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = paste0(organ[i],", ACV lipid, logFC(WD/chow)"), #PPAR signaling
            Rowv= F,Colv= F, 
            ColSideColors= rep(c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), times =3),
            #RowSideColors= resp.chain$V3,
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(5,10),     # widens margins around plot
            keysize=0.75,
            dendrogram="none", cexRow=0.8, cexCol=0.8) 
  dev.off()
  
}

### vis/ sc/ liver
tab.deg.cap = read.table(paste0(file.dir,"liver/deg-gen/table_deg_liver_Cap_all_TP.txt"), sep='\t')
tab.deg.art = read.table(paste0(file.dir,"liver/deg-gen/table_deg_liver_art_all_TP.txt"), sep='\t')
tab.deg.ven = read.table(paste0(file.dir,"liver/deg-gen/table_deg_liver_ven_all_TP.txt"), sep='\t')

tab.deg.cap = read.table(paste0(file.dir,"heart/deg-gen/table_deg_heart_Cap_all_TP.txt"), sep='\t')
g = c("Col15a1","Col4a2","Mmp15","Col4a1")

head(tab.deg.art)
a <- tab.deg.art[,c(1,3,5,7,9)] # 
a = a[rownames(a)%in%lipid,] 

head(tab.deg.cap)
b <- tab.deg.cap[,c(1,3,5,7,9)] # 
b = b[rownames(b)%in%g,] 

head(tab.deg.ven)
d <- tab.deg.ven[,c(1,3,5,7,9)] # 
d = d[rownames(d)%in%lipid,] 

b = cbind(a,b,d)
b <- b[order(match(rownames(b),lipid)),]
# replase NA with 0
b[is.na(b)] = 0

# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5

mi <- min(b,na.rm=TRUE) #l
ma <- max(b,na.rm=TRUE) #l

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
pdf(paste0("heatmap-liver-ACV-lipid.pdf"),width=15,height=10,paper='special')
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = paste0("vis ACV lipid, logFC(WD/chow)"), #PPAR signaling
          Rowv= T,Colv= F, 
          ColSideColors= rep(c("#ffe0a1","#fc8b00","#fc5000","#f0e02e","#a2f02e"), times =1),
          #RowSideColors= resp.chain$V3,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="row", cexRow=0.8, cexCol=0.8) 
dev.off()



###
# kidney cap 
kid.integr = c("Itgb1","Igta3","Lama5","Lamb2","Nid1","Tgm2")
kid.slc = c("Slc34a1","Slc4a4","Slc6a19","Slc5a2","Slc22a8","Slc13a1","Slc22a18","Slc5a12")

tab.deg.cap = read.table(paste0(file.dir,"kidney/deg-gen/table_deg_kidney_kidney1_all_TP.txt"), sep='\t')
tab.deg.cap = read.table(paste0(file.dir,"kidney/deg-gen/table_deg_kidney_kidney2_all_TP.txt"), sep='\t')
tab.deg.cap = read.table(paste0(file.dir,"kidney/deg-gen/table_deg_kidney_clust_all_WD6m.txt"), sep='\t')

head(tab.deg.cap)
b <- tab.deg.cap[,seq(1,ncol(tab.deg.cap),2)] # 
b = b[rownames(b)%in%resp.chain$V1,]  #glycolysis kid.slc

#b <- b[order(match(rownames(b),lipid)),]
# replase NA with 0
b[is.na(b)] = 0

# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5

mi <- min(b,na.rm=TRUE) #l
ma <- max(b,na.rm=TRUE) #l

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
pdf(paste0("heatmap-kidney-Kid2-respchain.pdf"),width=15,height=30,paper='special')
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = paste0("kidney mEC2 resp chain, logFC(WD/chow)"), #PPAR signaling
          Rowv= T,Colv= F, 
          ColSideColors= hue_pal()(ncol(b)),
          #RowSideColors= resp.chain$V3,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="row", cexRow=0.5, cexCol=0.8) 
dev.off()


######### combine ACV per organ per time point

organ = c("brain","lung", "heart", "liver", "kidney")
tab.deg.acv = list()
for (i in seq(1,length(organ),1)) {
  tab.deg.art[[i]]$gene = rownames(tab.deg.art[[i]])
  tab.deg.cap[[i]]$gene = rownames(tab.deg.cap[[i]])
  tab.deg.ven[[i]]$gene = rownames(tab.deg.ven[[i]])
  tab.deg.acv[[i]] = merge(tab.deg.art[[i]],tab.deg.cap[[i]], by = "gene")
  tab.deg.acv[[i]] = merge(tab.deg.acv[[i]],tab.deg.ven[[i]], by = "gene")
  
  write.table(tab.deg.acv[[i]], paste0("deg_ACV_3-4-6m_",organ[i],".txt"))
} 

head(tab.deg.acv[[1]]) 


organ = c("brain","lung", "heart", "liver", "kidney", "vis", "sc")
tab.deg.acv = list()

for (i in seq(1,length(organ),1)) {
  tab.deg.art[[i]]$gene = rownames(tab.deg.art[[i]])
  tab.deg.cap[[i]]$gene = rownames(tab.deg.cap[[i]])
  tab.deg.ven[[i]]$gene = rownames(tab.deg.ven[[i]])
}

# combine arteries
  tab.deg.acv[["art"]] = merge(tab.deg.art[[1]],tab.deg.art[[2]], by = "gene")
  tab.deg.acv[["art"]] = merge(tab.deg.acv[["art"]],tab.deg.art[[3]], by = "gene")
  tab.deg.acv[["art"]] = merge(tab.deg.acv[["art"]],tab.deg.art[[4]], by = "gene")
  tab.deg.acv[["art"]] = merge(tab.deg.acv[["art"]],tab.deg.art[[5]], by = "gene")
  tab.deg.acv[["art"]] = merge(tab.deg.acv[["art"]],tab.deg.art[[6]], by = "gene")
  tab.deg.acv[["art"]] = merge(tab.deg.acv[["art"]],tab.deg.art[[7]], by = "gene")

head(tab.deg.acv[["art"]])
ncol(tab.deg.acv[["art"]])

write.table(tab.deg.acv[["art"]], paste0("deg_art_3-4-6m_all_organs.txt"))

# combine capillary
tab.deg.acv[["cap"]] = merge(tab.deg.cap[[1]],tab.deg.cap[[2]], by = "gene")
tab.deg.acv[["cap"]] = merge(tab.deg.acv[["cap"]],tab.deg.cap[[3]], by = "gene")
tab.deg.acv[["cap"]] = merge(tab.deg.acv[["cap"]],tab.deg.cap[[4]], by = "gene")
tab.deg.acv[["cap"]] = merge(tab.deg.acv[["cap"]],tab.deg.cap[[5]], by = "gene")
tab.deg.acv[["cap"]] = merge(tab.deg.acv[["cap"]],tab.deg.cap[[6]], by = "gene")
tab.deg.acv[["cap"]] = merge(tab.deg.acv[["cap"]],tab.deg.cap[[7]], by = "gene")

head(tab.deg.acv[["cap"]])
ncol(tab.deg.acv[["cap"]])

write.table(tab.deg.acv[["cap"]], paste0("deg_cap_3-4-6m_all_organs.txt"))

# combine vein
tab.deg.acv[["ven"]] = merge(tab.deg.ven[[1]],tab.deg.ven[[2]], by = "gene")
tab.deg.acv[["ven"]] = merge(tab.deg.acv[["ven"]],tab.deg.ven[[3]], by = "gene")
tab.deg.acv[["ven"]] = merge(tab.deg.acv[["ven"]],tab.deg.ven[[4]], by = "gene")
tab.deg.acv[["ven"]] = merge(tab.deg.acv[["ven"]],tab.deg.ven[[5]], by = "gene")
tab.deg.acv[["ven"]] = merge(tab.deg.acv[["ven"]],tab.deg.ven[[6]], by = "gene")
tab.deg.acv[["ven"]] = merge(tab.deg.acv[["ven"]],tab.deg.ven[[7]], by = "gene")

head(tab.deg.acv[["ven"]])
ncol(tab.deg.acv[["ven"]])

write.table(tab.deg.acv[["ven"]], paste0("deg_ven_3-4-6m_all_organs.txt"))

tab.deg.acv = list()
tab.deg.acv[["art"]] = read.table("deg_art_3-4-6m_all_organs.txt")
tab.deg.acv[["cap"]] = read.table("deg_cap_3-4-6m_all_organs.txt")
tab.deg.acv[["ven"]] = read.table("deg_ven_3-4-6m_all_organs.txt")



klfs = c("Klf2","Klf3","Klf4","Klf6","Klf7","Klf9","Klf10","Klf12","Klf13","Klf16")
ap1 = c("Jun","Junb","Jund","Fos","Fosb","Egr1")

n = c("art", "cap", "ven")
# 6 months
for (i in c(1,2,3)) {
  
  # WD 6m c(6,16,26,36,46,54,64) ; Rev 3m c(10,20,30,40,50,58,68)
  # WD 4m c(4,14,24,34,44,54,64) ; Rev 1m c(8,18,28,38,48,58,68)
  # WD 3m c(2,12,22,32,42,52,62)
  b <- tab.deg.acv[[i]][,c(2,12,22,32,42,52,62)] 
  rownames(b) = tab.deg.acv[[i]]$gene
  b = b[rownames(b)%in%ap1,] #
  b <- b[order(match(rownames(b),ap1)),]
  
  # replase NA with 0
  b[is.na(b)] = 0
  # exclude -0.1<logFC<0.1 genes
  #b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]
  
  # make logFC max 0.5
  b[b>0.5] = 0.5
  b[b<(-0.5)] = -0.5
  
  mi <- min(b,na.rm=TRUE) #l
  ma <- max(b,na.rm=TRUE) #l
  b
  
  ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
  ColorLevels <- seq(mi, ma, length=length(ColorRamp))
  pdf(paste0("heatmap-",n[i],"-wd3m_ap1.pdf"),width=10,height=15,paper='special')
  heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = paste0("AP1 in ",n[i]," WD 3m, logFC(WD/chow)"), #PPAR signaling
            Rowv= F,Colv= F, 
            ColSideColors= my_color_palette,
            #RowSideColors= resp.chain$V3,
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(5,10),     # widens margins around plot
            keysize=0.75,
            dendrogram="none", cexRow=0.8, cexCol=0.8) 
  dev.off()
}
