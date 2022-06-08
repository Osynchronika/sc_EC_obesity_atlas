library(Seurat)
library(dplyr)
library(Matrix)
require(gplots)
require(ggplot2)
library(RColorBrewer)
library(cowplot)
require(scales)
library(Rfast)
library(corrplot)

setwd("~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/comparisons_new")

#############
# correlation
#############
# read in the DEG for capillaries from all organs
organ = c("brain", "lung", "heart", "liver", "kidney","lung","sc")

# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(organ))

tab.deg.cap = read.table("table_deg_capillary_all_organs.txt", sep='\t')
tab.deg.art = read.table("table_deg_artery_all_organs.txt", sep='\t')
tab.deg.ven = read.table("table_deg_vein_all_organs.txt", sep='\t')
tab.deg.all.logFC = read.table("tab_deg_all_logFC.txt", sep = '\t', header = T)

tab.deg.cap.logFC = tab.deg.cap[,seq(1,14,2)] # 
tab.deg.art.logFC = tab.deg.art[,seq(1,14,2)] # 
tab.deg.ven.logFC = tab.deg.ven[,seq(1,14,2)] # 
# brain heart lung kidney liver vis sc
tab.deg.cap.logFC = tab.deg.cap.logFC[,c(1,3,2,5,4,6,7)] # 
tab.deg.art.logFC = tab.deg.art.logFC[,c(1,3,2,5,4,6,7)] # 
tab.deg.ven.logFC = tab.deg.ven.logFC[,c(1,3,2,5,4,6,7)] # 
####
# 1. Capillaries

deg.capillaries = list()
deg.i = NULL
for (i in organ) {
  print(i)
  deg.i = read.table(paste("deg_gen/",i,"_EC-cap_Western_vs_chow.txt", sep = "", collapse = NULL), sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.capillaries[[i]] = deg.i
}

deg.capillaries$liver

# combine a table of logFC
tab.deg.cap = merge(deg.capillaries$brain[,c(2,5)], deg.capillaries$lung[,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.cap) = c("gene","brain_logFC","brain_Q-val","lung_logFC","lung_Q-val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries$heart[,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[6:7] = c("heart_logFC","heart_Q-val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries$liver[,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[8:9] = c("liver_logFC","liver_Q-val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries$kidney[,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[10:11] = c("kidney_logFC","kidney_Q-val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries$vis[,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[12:13] = c("vis_logFC","vis_Q-val")
tab.deg.cap = merge(tab.deg.cap, deg.capillaries$sc[,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.cap)[14:15] = c("sc_logFC","sc_Q-val")
head(tab.deg.cap)

rownames(tab.deg.cap) = tab.deg.cap$gene
tab.deg.cap = tab.deg.cap[,-1]
nrow(tab.deg.cap) # 8390 genes
write.table(tab.deg.cap,"table_deg_capillary_all_organs.txt", sep='\t')
tab.deg.cap.logFC = tab.deg.cap[,seq(1,14,2)] # tab.deg.cap[,seq(1,14,2)]

# Look for inflammation signature
nfkb = read.table("../NFkB.txt", sep='\t', header = T)
nfkb = nfkb$Symbol
nfkb = nfkb[-67]
library(stringr)
nfkb = str_to_sentence(nfkb, locale = "en")

TCA = read.table("../TCA_cycle.txt", sep ='\t')
TCA = TCA$V1
library(stringr)
TCA = str_to_sentence(TCA, locale = "en")
# order TCA cycle
TCA = TCA[c(18,19,20,21,22,23,5,6,1,4,2,3,9,10,11,12,13,7,16,17,28,29,30,31,24,25,26,27,8,14,15)]

chronic = read.table("../GO_0002544_chronic_inflammation.txt", sep='\t')
chronic = chronic$V1
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

# coagulation cascade
coagul = c("F3","F7","F8","F9","F10","F11","F12","Vwf","F2","F2r","F2rl2","F2rl3",
           "C3","C5","A2m","Serpina1","Serpina5","Serpinf2",
           "Serpine1","Serpinb2","Fga","Fgb","Fgg", # pro-coag
           "Proc","Serpina5","Serpinc1","Serpind1","Tfpi","Plg","Plau","Plat") # anti-coag

fatty.acid.trans = c("Fabp1","Fabp2","Fabp3","Fabp4","Fabp5","Fabp6","Fabp7","Cd36","Lpl","Slc27a1","Slc27a4","Dbi",
                     "Acsbg1","Acsbg2","Acsl1","Acsl2","Acsl3","Acsl4","Acsl5","Acsl6") # enzymes


# heatmap of the logFC WD/chow
b <- tab.deg.cap.logFC[rownames(tab.deg.cap.logFC) %in% nfkb,]
# replase NA with 0
b[is.na(b)] = 0
b

# exclude -0.1<logFC<0.1 genes
library(Rfast)
b = b[(rowMaxs(as.matrix(b),value = T)>=0.15)|(rowMins(as.matrix(b),value = T)<=(-0.15)),]

# make logFC max 0.5
b[b>0.5] = 0.5

mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = "NFkB targets in capillaries, logFC(WD/chow)", 
          Rowv= T,Colv= T, 
          ColSideColors= my_color_palette, # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="both", cexRow=0.8, cexCol=0.8) 



# heatmap of the logFC WD/chow
b <- tab.deg.cap.logFC[rownames(tab.deg.cap.logFC) %in% inflamm,]
# replase NA with 0
b[is.na(b)] = 0
b

# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.15)|(rowMins(as.matrix(b),value = T)<=(-0.15)),]

# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5

mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = "GO inflammation in capillaries, logFC(WD/chow)", 
          Rowv= T,Colv= T, 
          ColSideColors= my_color_palette, # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          dendrogram="both", cexRow=0.8, cexCol=0.8) 


# heatmap of the logFC WD/chow
b <- tab.deg.cap.logFC[rownames(tab.deg.cap.logFC) %in% glycolysis,] #FA.oxid
# replase NA with 0
b[is.na(b)] = 0
b

# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.15)|(rowMins(as.matrix(b),value = T)<=(-0.15)),]

# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5

mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = "Glycolysis in capillaries, logFC(WD/chow)", #FA oxidation
          Rowv= T,Colv= T, 
          ColSideColors= my_color_palette, # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="both", cexRow=0.8, cexCol=0.8) 


# heatmap of the logFC WD/chow
b <- tab.deg.cap.logFC[rownames(tab.deg.cap.logFC) %in% FA.trans,] #FA.oxid glycolysis
# replase NA with 0
b[is.na(b)] = 0
b
# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.15)|(rowMins(as.matrix(b),value = T)<=(-0.15)),]
# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5
mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = "FA transport in capillaries, logFC(WD/chow)", #FA oxidation Glycolysis
          Rowv= T,Colv= T, 
          ColSideColors= my_color_palette, # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          dendrogram="both", cexRow=0.8, cexCol=0.8) 

# heatmap of the logFC WD/chow
b <- tab.deg.cap.logFC[rownames(tab.deg.cap.logFC) %in% ecm,] #FA.oxid glycolysis
b <- tab.deg.art.logFC[rownames(tab.deg.art.logFC) %in% ecm,] #FA.oxid glycolysis
# replase NA with 0
b[is.na(b)] = 0
b

# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.2)|(rowMins(as.matrix(b),value = T)<=(-0.2)),]

# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5

mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = "ECM in capillaries, logFC(WD/chow)", #FA oxidation Glycolysis
          Rowv= T,Colv= F, 
          ColSideColors= my_color_palette, # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="row", cexRow=0.8, cexCol=0.8) 

# coagulation
# heatmap of the logFC WD/chow
b <- tab.deg.cap.logFC[rownames(tab.deg.cap.logFC) %in% platelet.act,] #coagul FA.oxid glycolysis
b <- tab.deg.art.logFC[rownames(tab.deg.art.logFC) %in% platelet.act,] #FA.oxid glycolysis
b <- tab.deg.ven.logFC[rownames(tab.deg.ven.logFC) %in% platelet.act,] #FA.oxid glycolysis
# replase NA with 0
b[is.na(b)] = 0
b
# match order - ordertable b rows in the marker list order
b <- b[order(match(rownames(b),coagul)),]

# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.2)|(rowMins(as.matrix(b),value = T)<=(-0.2)),]

# make logFC max 0.5
b[b>0.5] = 0.5

mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
          main = "Platelet activation in vein, logFC(WD/chow)", #Platelet activation Coagulation FA oxidation Glycolysis
          Rowv= T,Colv= T, 
          ColSideColors= my_color_palette, # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="both", cexRow=0.8, cexCol=0.8) 


# fatty acid transport
b <- tab.deg.cap.logFC[rownames(tab.deg.cap.logFC) %in% fatty.acid.trans,] #coagul FA.oxid glycolysis
b <- tab.deg.art.logFC[rownames(tab.deg.art.logFC) %in% fatty.acid.trans,] #FA.oxid glycolysis
b <- tab.deg.ven.logFC[rownames(tab.deg.ven.logFC) %in% fatty.acid.trans,] #FA.oxid glycolysis
# replase NA with 0
b[is.na(b)] = 0
b
# match order - ordertable b rows in the marker list order
b <- b[order(match(rownames(b),fatty.acid.trans)),]
# exclude -0.1<logFC<0.1 genes
#b = b[(rowMaxs(as.matrix(b),value = T)>=0.2)|(rowMins(as.matrix(b),value = T)<=(-0.2)),]
# make logFC max 0.5
b[b>0.5] = 0.5
mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma
ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
          main = "FA transport in vein, logFC(WD/chow)", #Platelet activation Coagulation FA oxidation Glycolysis
          Rowv= T,Colv= T, 
          ColSideColors= my_color_palette, # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="both", cexRow=0.8, cexCol=0.8) 

######################
# combine art-cap-vein into 1 table
colnames(tab.deg.art.logFC) = paste("art_", colnames(tab.deg.art.logFC),sep="")
tab.deg.art.logFC$gene = rownames(tab.deg.art.logFC)
colnames(tab.deg.cap.logFC) = paste("cap_", colnames(tab.deg.cap.logFC),sep="")
tab.deg.cap.logFC$gene = rownames(tab.deg.cap.logFC)
colnames(tab.deg.ven.logFC) = paste("ven_", colnames(tab.deg.ven.logFC),sep="")
tab.deg.ven.logFC$gene = rownames(tab.deg.ven.logFC)

tab.deg.all.logFC = merge(tab.deg.art.logFC,tab.deg.cap.logFC, by = "gene")
tab.deg.all.logFC = merge(tab.deg.all.logFC,tab.deg.ven.logFC, by = "gene")
# remove logFC from colnames
colnames(tab.deg.all.logFC) = gsub('.{6}$', '', colnames(tab.deg.all.logFC))
head(tab.deg.all.logFC)
rownames(tab.deg.all.logFC) = tab.deg.all.logFC$gene
tab.deg.all.logFC = tab.deg.all.logFC[,-1]
write.table(tab.deg.all.logFC, "tab_deg_all_logFC.txt", sep = '\t')

tab.deg.all.logFC = read.table("tab_deg_all_logFC.txt", sep = '\t', header = T)
# art cap ven order
#tab.deg.all.logFC = tab.deg.all.logFC[,c(1,8,15,2,9,16,3,10,17,4,11,18,5,12,19,6,13,20,7,14,21)]
# brain heart lung kidney liver vis sc
#tab.deg.all.logFC = tab.deg.all.logFC[,c(1,2,3,7,8,9,4,5,6,13,14,15,10,11,12,16,17,18,19,20,21)]
colnames(tab.deg.all.logFC)

beta.oxid = c("Cpt1", "Cpt2","Acadvl","Hadha","Hadhb","Acadm","Ehhadh", "Acat1", "Acaa1", "Acaa2")
sel.fa.trans = c("Lpl", "Cd36", "Fabp1", "Fabp4", "Fabp5","Dbi")
calcium = read.table("../GO_0019722_calcium_mediated_signal.txt", sep="\t")
calcium = unique(calcium$V1)
ca2.ec = c("Anxa2","Mylk","Camk2d","F2r","Trpc6","Icam1","Vcam1","Stim1")
leuk.recruit = c("Icam1","Icam2","Pecam1","Vcam1", "Sele", "Selp","Cd99","Hspg","Ptafr",
                 "Il1r","Tnfr1","Trap","Traf2","Traf6","Tirap","Tradd","Myd88","Irak1", "Irak4",
                 "Nfkb1","Nfbk2","Rela","Relb","Rel","Cox1","Cox2","Cxcl8","Cxcl10","Ccl26",
                 "Il6","Il1","Il4","Il7")
#c("Ace1","Agtr1","Agtr2")

ecm = c("Col4a1", "Col4a2", "Col13a1", "Col15a1", "Lama4", "Lama5", "Hspg2", "Lamb1", "Lamb2", "Slit3", "Thbs1")
infl = c("F2r","Cxcl9","Rela", "Relb", "Nfkb", "Il1a","Tgfb1", "Icam2", "Vegf")

b <- tab.deg.all.logFC[rownames(tab.deg.all.logFC) %in% infl,] #leuk.recruit ca2.ec sel.fa.trans TCA coagul nfkb inflamm glycolysis platelet.act
# replase NA with 0
b[is.na(b)] = 0
b
# match order - ordertable b rows in the marker list order
b <- b[order(match(rownames(b),junctions)),] # TCA
# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]
# make logFC max 0.5
b[b>0.3] = 0.3
b[b<(-0.3)] = -0.3
mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma
ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, 
          main = "Inflam selected, logFC(WD/chow)", #Leukocyte recruitment Ca2+ signaling Selected FA trans  Coagulation NFkB targets Inflammation Glycolysis TCA cycle
          Rowv= T,Colv= F, 
          ColSideColors= rep(my_color_palette,each=3), # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="row", cexRow=0.8, cexCol=0.8) 


###########################
# signif DEGs

sign.genes = rownames(tab.deg.cap[(tab.deg.cap$brain_Q.val<0.05)|(tab.deg.cap$lung_Q.val<0.05)|(tab.deg.cap$heart_Q.val<0.05)|(tab.deg.cap$liver_Q.val<0.05)|(tab.deg.cap$kidney_Q.val<0.05)|(tab.deg.cap$vis_Q.val<0.05)|(tab.deg.cap$sc_Q.val<0.05),])
sg = rownames(tab.deg.cap)[rownames(tab.deg.cap)%in%sign.genes]
tab.deg.cap.sign = tab.deg.cap[sg,]
nrow(tab.deg.cap.sign) # 2602 genes with signif Q-val
head(tab.deg.cap.sign)
write.table(tab.deg.cap.sign, "capillary_deg_sign.txt", sep = '\t')

# keep only logFC
tab.deg.cap.logFC = tab.deg.cap.sign[,seq(1,14,2)] # tab.deg.cap[,seq(1,14,2)]
head(tab.deg.cap.logFC)

# replace NA with 0s
tab.deg.cap.logFC[is.na(tab.deg.cap.logFC)] = 0
head(tab.deg.cap.logFC)

# remove rows with all 0s
tab.deg.cap.logFC = tab.deg.cap.logFC[!rowSums(tab.deg.cap.logFC)==0,]
head(tab.deg.cap.logFC)
nrow(tab.deg.cap.logFC) # 1657 genes with logFC

# plot correlations

pdf(paste0("corr-cap-brain.pdf"),width=11.5,height=6,paper='special')
p1 <- ggplot(tab.deg.cap.logFC, aes(brain_logFC, brain_logFC)) + geom_point() + ggtitle("Correlation brain vs brain") + theme_linedraw() + geom_smooth(method = "lm")
p2 <- ggplot(tab.deg.cap.logFC, aes(brain_logFC, heart_logFC)) + geom_point() + ggtitle("Correlation brain vs heart") + theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.cap.logFC, aes(brain_logFC, lung_logFC)) + geom_point() + ggtitle("Correlation brain vs lung") + theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.cap.logFC, aes(brain_logFC, liver_logFC)) + geom_point() + ggtitle("Correlation brain vs liver") + theme_linedraw() + geom_smooth(method = "lm")
p5 <- ggplot(tab.deg.cap.logFC, aes(brain_logFC, kidney_logFC)) + geom_point() + ggtitle("Correlation brain vs kidney") + theme_linedraw() + geom_smooth(method = "lm")
p6 <- ggplot(tab.deg.cap.logFC, aes(brain_logFC, vis_logFC)) + geom_point() + ggtitle("Correlation brain vs vis fat") + theme_linedraw() + geom_smooth(method = "lm")
p7 <- ggplot(tab.deg.cap.logFC, aes(brain_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation brain vs sc fat") + theme_linedraw() + geom_smooth(method = "lm")
print((p1 | p2 | p3 |p4) /
        (p5 | p6 | p7| p7)) # 11.5x6
dev.off()

pdf(paste0("corr-cap-heart.pdf"),width=8.5,height=6,paper='special')
p2 <- ggplot(tab.deg.cap.logFC, aes(heart_logFC, heart_logFC)) + geom_point() + ggtitle("Correlation heart vs heart")+ theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.cap.logFC, aes(heart_logFC, lung_logFC)) + geom_point() + ggtitle("Correlation heart vs lung")+ theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.cap.logFC, aes(heart_logFC, liver_logFC)) + geom_point() + ggtitle("Correlation heart vs liver")+ theme_linedraw() + geom_smooth(method = "lm")
p5 <- ggplot(tab.deg.cap.logFC, aes(heart_logFC, kidney_logFC)) + geom_point() + ggtitle("Correlation heart vs kidney")+ theme_linedraw() + geom_smooth(method = "lm")
p6 <- ggplot(tab.deg.cap.logFC, aes(heart_logFC, vis_logFC)) + geom_point() + ggtitle("Correlation heart vs vis fat")+ theme_linedraw() + geom_smooth(method = "lm")
p7 <- ggplot(tab.deg.cap.logFC, aes(heart_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation heart vs sc fat")+ theme_linedraw() + geom_smooth(method = "lm")
print(( p2 | p3 |p4) /
        (p5 | p6 | p7)) # 8.5x6
dev.off()

pdf(paste0("corr-cap-lung.pdf"),width=8.5,height=6,paper='special')
p3 <- ggplot(tab.deg.cap.logFC, aes(lung_logFC, lung_logFC)) + geom_point() + ggtitle("Correlation lung vs lung")+ theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.cap.logFC, aes(lung_logFC, liver_logFC)) + geom_point() + ggtitle("Correlation lung vs liver")+ theme_linedraw() + geom_smooth(method = "lm")
p5 <- ggplot(tab.deg.cap.logFC, aes(lung_logFC, kidney_logFC)) + geom_point() + ggtitle("Correlation lung vs kidney")+ theme_linedraw() + geom_smooth(method = "lm")
p6 <- ggplot(tab.deg.cap.logFC, aes(lung_logFC, vis_logFC)) + geom_point() + ggtitle("Correlation lung vs vis fat")+ theme_linedraw() + geom_smooth(method = "lm")
p7 <- ggplot(tab.deg.cap.logFC, aes(lung_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation lung vs sc fat")+ theme_linedraw() + geom_smooth(method = "lm")
print(( p3 | p3 |p4) /
        (p5 | p6 | p7)) # 8.5x6
dev.off()

pdf(paste0("corr-cap-liver.pdf"),width=5.8,height=6,paper='special')
p4 <- ggplot(tab.deg.cap.logFC, aes(liver_logFC, liver_logFC)) + geom_point() + ggtitle("Correlation liver vs liver")+ theme_linedraw() + geom_smooth(method = "lm")
p5 <- ggplot(tab.deg.cap.logFC, aes(liver_logFC, kidney_logFC)) + geom_point() + ggtitle("Correlation liver vs kidney")+ theme_linedraw() + geom_smooth(method = "lm")
p6 <- ggplot(tab.deg.cap.logFC, aes(liver_logFC, vis_logFC)) + geom_point() + ggtitle("Correlation liver vs vis fat")+ theme_linedraw() + geom_smooth(method = "lm")
p7 <- ggplot(tab.deg.cap.logFC, aes(liver_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation liver vs sc fat")+ theme_linedraw() + geom_smooth(method = "lm")
print(( p4 | p5 ) /
        ( p6 | p7))# 5.8x6
dev.off()

pdf(paste0("corr-cap-kidney.pdf"),width=8.3,height=3,paper='special')
p5 <- ggplot(tab.deg.cap.logFC, aes(kidney_logFC, kidney_logFC)) + geom_point() + ggtitle("Correlation kidney vs kidney")+ theme_linedraw() + geom_smooth(method = "lm")
p6 <- ggplot(tab.deg.cap.logFC, aes(kidney_logFC, vis_logFC)) + geom_point() + ggtitle("Correlation kidney vs vis fat")+ theme_linedraw() + geom_smooth(method = "lm")
p7 <- ggplot(tab.deg.cap.logFC, aes(kidney_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation kidney vs sc fat")+ theme_linedraw() + geom_smooth(method = "lm")
print(( p5 |p6 | p7)) # 8.3x3
dev.off()

pdf(paste0("corr-cap-fat.pdf"),width=8.3,height=3,paper='special')
p5 <- ggplot(tab.deg.cap.logFC, aes(vis_logFC, vis_logFC)) + geom_point() + ggtitle("Correlation vis vs vis fat")+ theme_linedraw() + geom_smooth(method = "lm")
p6 <- ggplot(tab.deg.cap.logFC, aes(vis_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation vis vs sc fat")+ theme_linedraw() + geom_smooth(method = "lm")
p7 <- ggplot(tab.deg.cap.logFC, aes(sc_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation sc vs sc fat")+ theme_linedraw() + geom_smooth(method = "lm")
print((p5 | p6 | p7)) # 8.3x3
dev.off()

# all clusters
library(corrplot)
m <- cor(tab.deg.cap.logFC)
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black",
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))



###########################
# DEGs abs(logFC)>=0.1

sign.genes = rownames(tab.deg.cap[(abs(tab.deg.cap$brain_logFC)>=0.1)|(abs(tab.deg.cap$lung_logFC)>=0.1)|
                                    (abs(tab.deg.cap$heart_logFC)>=0.1)|(abs(tab.deg.cap$liver_logFC)>=0.1)|
                                    (abs(tab.deg.cap$kidney_logFC)>=0.1)|(abs(tab.deg.cap$vis_logFC)>=0.1)|
                                    (abs(tab.deg.cap$sc_logFC)>=0.1),])
sg = rownames(tab.deg.cap)[rownames(tab.deg.cap)%in%sign.genes]
tab.deg.cap.sign = tab.deg.cap[sg,]
nrow(tab.deg.cap.sign) # 2242 genes with logFC>=0.1
head(tab.deg.cap.sign)
write.table(tab.deg.cap.sign, "capillary_deg_logFC0.1.txt", sep = '\t')

# keep only logFC
tab.deg.cap.logFC = tab.deg.cap.sign[,seq(1,14,2)] # tab.deg.cap[,seq(1,14,2)]
head(tab.deg.cap.logFC)

# replace NA with 0s
tab.deg.cap.logFC[is.na(tab.deg.cap.logFC)] = 0
head(tab.deg.cap.logFC)

# remove rows with all 0s
tab.deg.cap.logFC = tab.deg.cap.logFC[!rowSums(tab.deg.cap.logFC)==0,]
head(tab.deg.cap.logFC)
nrow(tab.deg.cap.logFC) # 1657 genes with logFC

# all clusters
library(corrplot)
m <- cor(tab.deg.cap.logFC)
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black",
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))

###########################
# DEGs abs(logFC)>=0.1

sign.genes = rownames(tab.deg.art[(abs(tab.deg.art$brain_logFC)>=0.1)|(abs(tab.deg.art$lung_logFC)>=0.1)|
                                    (abs(tab.deg.art$heart_logFC)>=0.1)|(abs(tab.deg.art$liver_logFC)>=0.1)|
                                    (abs(tab.deg.art$kidney_logFC)>=0.1)|(abs(tab.deg.art$vis_logFC)>=0.1)|
                                    (abs(tab.deg.art$sc_logFC)>=0.1),])
sg = rownames(tab.deg.art)[rownames(tab.deg.art)%in%sign.genes]
tab.deg.art.sign = tab.deg.art[sg,]
nrow(tab.deg.art.sign) # 6240 genes with logFC>=0.1
head(tab.deg.art.sign)
write.table(tab.deg.art.sign, "art_deg_logFC0.1.txt", sep = '\t')

# keep only logFC
tab.deg.art.logFC = tab.deg.art.sign[,seq(1,14,2)] # tab.deg.art[,seq(1,14,2)]
head(tab.deg.art.logFC)

# replace NA with 0s
tab.deg.art.logFC[is.na(tab.deg.art.logFC)] = 0
head(tab.deg.art.logFC)

# remove rows with all 0s
tab.deg.art.logFC = tab.deg.art.logFC[!rowSums(tab.deg.art.logFC)==0,]
head(tab.deg.art.logFC)
nrow(tab.deg.art.logFC) # 6240 genes with logFC

# all clusters
library(corrplot)
m <- cor(tab.deg.art.logFC)
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black",
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))
###########################
# DEGs abs(logFC)>=0.1

sign.genes = rownames(tab.deg.ven[(abs(tab.deg.ven$brain_logFC)>=0.1)|(abs(tab.deg.ven$lung_logFC)>=0.1)|
                                    (abs(tab.deg.ven$heart_logFC)>=0.1)|(abs(tab.deg.ven$liver_logFC)>=0.1)|
                                    (abs(tab.deg.ven$kidney_logFC)>=0.1)|(abs(tab.deg.ven$vis_logFC)>=0.1)|
                                    (abs(tab.deg.ven$sc_logFC)>=0.1),])
sg = rownames(tab.deg.ven)[rownames(tab.deg.ven)%in%sign.genes]
tab.deg.ven.sign = tab.deg.ven[sg,]
nrow(tab.deg.ven.sign) # 3542 genes with logFC>=0.1
head(tab.deg.ven.sign)
write.table(tab.deg.ven.sign, "ven_deg_logFC0.1.txt", sep = '\t')

# keep only logFC
tab.deg.ven.logFC = tab.deg.ven.sign[,seq(1,14,2)] # tab.deg.ven[,seq(1,14,2)]
head(tab.deg.ven.logFC)

# replace NA with 0s
tab.deg.ven.logFC[is.na(tab.deg.ven.logFC)] = 0
head(tab.deg.ven.logFC)

# remove rows with all 0s
tab.deg.ven.logFC = tab.deg.ven.logFC[!rowSums(tab.deg.ven.logFC)==0,]
head(tab.deg.ven.logFC)
nrow(tab.deg.ven.logFC) # 3542 genes with logFC

# all clusters
library(corrplot)
m <- cor(tab.deg.ven.logFC)
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black",
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))

##################
# correlate fats
tab.deg.cap.sign.fat = tab.deg.cap.sign[,c(11,12,13,14)] 
# signif DEGs
tab.deg.cap.sign.fat = tab.deg.cap.sign.fat[(tab.deg.cap.sign.fat$`vis_Q-val`<0.05)|(tab.deg.cap.sign.fat$`sc_Q-val`<0.05),]
# rm NAs
tab.deg.cap.sign.fat = tab.deg.cap.sign.fat[complete.cases(tab.deg.cap.sign.fat),]
nrow(tab.deg.cap.sign.fat) # 1140 genes with signif Q-val
head(tab.deg.cap.sign.fat)
tab.deg.cap.sign.fat = tab.deg.cap.sign.fat[,c(1,3)]
ggplot(tab.deg.cap.sign.fat, aes(vis_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation vis vs sc fat")+ theme_linedraw() + geom_smooth(method = "lm")

####
# 2. Arteries

deg.artery = list()
deg.i = NULL
for (i in organ) {
  print(i)
  deg.i = read.table(paste("deg_gen/",i,"_EC-art_Western_vs_chow.txt", sep = "", collapse = NULL), sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.artery[[i]] = deg.i
}

deg.artery$liver

# combine a table of logFC
tab.deg.art = merge(deg.artery$brain[,c(2,5)], deg.artery$lung[,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.art) = c("gene","brain_logFC","brain_Q-val","lung_logFC","lung_Q-val")
tab.deg.art = merge(tab.deg.art, deg.artery$heart[,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[6:7] = c("heart_logFC","heart_Q-val")
tab.deg.art = merge(tab.deg.art, deg.artery$liver[,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[8:9] = c("liver_logFC","liver_Q-val")
tab.deg.art = merge(tab.deg.art, deg.artery$kidney[,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[10:11] = c("kidney_logFC","kidney_Q-val")
tab.deg.art = merge(tab.deg.art, deg.artery$vis[,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[12:13] = c("vis_logFC","vis_Q-val")
tab.deg.art = merge(tab.deg.art, deg.artery$sc[,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.art)[14:15] = c("sc_logFC","sc_Q-val")
head(tab.deg.art)

rownames(tab.deg.art) = tab.deg.art$gene
tab.deg.art = tab.deg.art[,-1]
nrow(tab.deg.art) # 9392 genes
write.table(tab.deg.art,"table_deg_artery_all_organs.txt", sep='\t')

tab.deg.art.logFC = tab.deg.art[,seq(1,14,2)] # tab.deg.cap[,seq(1,14,2)]
# heatmap of the logFC WD/chow
b <- tab.deg.art.logFC[rownames(tab.deg.art.logFC) %in% nfkb,]
# replase NA with 0
b[is.na(b)] = 0
b

# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.15)|(rowMins(as.matrix(b),value = T)<=(-0.15)),]

# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5

mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = "NFkB targets in arteries, logFC(WD/chow)", 
          Rowv= T,Colv= T, 
          ColSideColors= my_color_palette, # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          dendrogram="both", cexRow=0.8, cexCol=0.8) 


# heatmap of the logFC WD/chow
b <- tab.deg.art.logFC[rownames(tab.deg.art.logFC) %in% inflamm,]
# replase NA with 0
b[is.na(b)] = 0
b

# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.2)|(rowMins(as.matrix(b),value = T)<=(-0.2)),]

# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5

mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = "GO inflammation in arteries, logFC(WD/chow)", 
          Rowv= T,Colv= T, 
          ColSideColors= my_color_palette, # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          dendrogram="both", cexRow=0.8, cexCol=0.8) 


# heatmap of the logFC WD/chow
b <- tab.deg.art.logFC[rownames(tab.deg.art.logFC) %in% glycolysis,] #FA.oxid glycolysis
# replase NA with 0
b[is.na(b)] = 0
b

# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.15)|(rowMins(as.matrix(b),value = T)<=(-0.15)),]

# make logFC max 0.5
b[b>0.5] = 0.5

mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = "Glycolysis in arteries, logFC(WD/chow)", #FA oxidation Glycolysis
          Rowv= T,Colv= T, 
          ColSideColors= my_color_palette, # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          dendrogram="both", cexRow=0.8, cexCol=0.8) 

# heatmap of the logFC WD/chow
b <- tab.deg.art.logFC[rownames(tab.deg.art.logFC) %in% FA.trans,] #FA.oxid glycolysis
# replase NA with 0
b[is.na(b)] = 0
b

# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.15)|(rowMins(as.matrix(b),value = T)<=(-0.15)),]

# make logFC max 0.5
b[b>0.5] = 0.5

mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = "FA transport in arteries, logFC(WD/chow)", #FA oxidation Glycolysis
          Rowv= T,Colv= T, 
          ColSideColors= my_color_palette, # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          dendrogram="both", cexRow=0.8, cexCol=0.8) 

# signif DEGs
sign.genes = rownames(tab.deg.art[(tab.deg.art$brain_Q.val<0.05)|(tab.deg.art$lung_Q.val<0.05)|(tab.deg.art$heart_Q.val<0.05)|(tab.deg.art$liver_Q.val<0.05)|(tab.deg.art$kidney_Q.val<0.05)|(tab.deg.art$vis_Q.val<0.05)|(tab.deg.art$sc_Q.val<0.05),])
sg = rownames(tab.deg.art)[rownames(tab.deg.art)%in%sign.genes]
tab.deg.art.sign = tab.deg.art[sg,]
nrow(tab.deg.art.sign) # 228 genes with signif Q-val
head(tab.deg.art.sign)
write.table(tab.deg.art.sign, "artery_deg_sign.txt", sep = '\t')

# keep only logFC
tab.deg.art.logFC = tab.deg.art.sign[,seq(1,14,2)] # tab.deg.art[,seq(1,14,2)]
head(tab.deg.art.logFC)

# replace NA with 0s
tab.deg.art.logFC[is.na(tab.deg.art.logFC)] = 0
head(tab.deg.art.logFC)

# remove rows with all 0s
tab.deg.art.logFC = tab.deg.art.logFC[!rowSums(tab.deg.art.logFC)==0,]
head(tab.deg.art.logFC)
nrow(tab.deg.art.logFC) # 177 genes with logFC
# plot correlations

pdf(paste0("corr-art-brain.pdf"),width=11.5,height=6,paper='special')
p1 <- ggplot(tab.deg.art.logFC, aes(brain_logFC, brain_logFC)) + geom_point() + ggtitle("Correlation brain vs brain") + theme_linedraw() + geom_smooth(method = "lm")
p2 <- ggplot(tab.deg.art.logFC, aes(brain_logFC, heart_logFC)) + geom_point() + ggtitle("Correlation brain vs heart") + theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.art.logFC, aes(brain_logFC, lung_logFC)) + geom_point() + ggtitle("Correlation brain vs lung") + theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.art.logFC, aes(brain_logFC, liver_logFC)) + geom_point() + ggtitle("Correlation brain vs liver") + theme_linedraw() + geom_smooth(method = "lm")
p5 <- ggplot(tab.deg.art.logFC, aes(brain_logFC, kidney_logFC)) + geom_point() + ggtitle("Correlation brain vs kidney") + theme_linedraw() + geom_smooth(method = "lm")
p6 <- ggplot(tab.deg.art.logFC, aes(brain_logFC, vis_logFC)) + geom_point() + ggtitle("Correlation brain vs vis fat") + theme_linedraw() + geom_smooth(method = "lm")
p7 <- ggplot(tab.deg.art.logFC, aes(brain_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation brain vs sc fat") + theme_linedraw() + geom_smooth(method = "lm")
print((p1 | p2 | p3 |p4) /
        (p5 | p6 | p7| p7)) # 11.5x6
dev.off()

pdf(paste0("corr-art-heart.pdf"),width=8.5,height=6,paper='special')
p2 <- ggplot(tab.deg.art.logFC, aes(heart_logFC, heart_logFC)) + geom_point() + ggtitle("Correlation heart vs heart")+ theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.art.logFC, aes(heart_logFC, lung_logFC)) + geom_point() + ggtitle("Correlation heart vs lung")+ theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.art.logFC, aes(heart_logFC, liver_logFC)) + geom_point() + ggtitle("Correlation heart vs liver")+ theme_linedraw() + geom_smooth(method = "lm")
p5 <- ggplot(tab.deg.art.logFC, aes(heart_logFC, kidney_logFC)) + geom_point() + ggtitle("Correlation heart vs kidney")+ theme_linedraw() + geom_smooth(method = "lm")
p6 <- ggplot(tab.deg.art.logFC, aes(heart_logFC, vis_logFC)) + geom_point() + ggtitle("Correlation heart vs vis fat")+ theme_linedraw() + geom_smooth(method = "lm")
p7 <- ggplot(tab.deg.art.logFC, aes(heart_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation heart vs sc fat")+ theme_linedraw() + geom_smooth(method = "lm")
print(( p2 | p3 |p4) /
        (p5 | p6 | p7)) # 8.5x6
dev.off()

pdf(paste0("corr-art-lung.pdf"),width=8.5,height=6,paper='special')
p3 <- ggplot(tab.deg.art.logFC, aes(lung_logFC, lung_logFC)) + geom_point() + ggtitle("Correlation lung vs lung")+ theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.art.logFC, aes(lung_logFC, liver_logFC)) + geom_point() + ggtitle("Correlation lung vs liver")+ theme_linedraw() + geom_smooth(method = "lm")
p5 <- ggplot(tab.deg.art.logFC, aes(lung_logFC, kidney_logFC)) + geom_point() + ggtitle("Correlation lung vs kidney")+ theme_linedraw() + geom_smooth(method = "lm")
p6 <- ggplot(tab.deg.art.logFC, aes(lung_logFC, vis_logFC)) + geom_point() + ggtitle("Correlation lung vs vis fat")+ theme_linedraw() + geom_smooth(method = "lm")
p7 <- ggplot(tab.deg.art.logFC, aes(lung_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation lung vs sc fat")+ theme_linedraw() + geom_smooth(method = "lm")
print(( p3 | p3 |p4) /
        (p5 | p6 | p7)) # 8.5x6
dev.off()

pdf(paste0("corr-art-liver.pdf"),width=5.8,height=6,paper='special')
p4 <- ggplot(tab.deg.art.logFC, aes(liver_logFC, liver_logFC)) + geom_point() + ggtitle("Correlation liver vs liver")+ theme_linedraw() + geom_smooth(method = "lm")
p5 <- ggplot(tab.deg.art.logFC, aes(liver_logFC, kidney_logFC)) + geom_point() + ggtitle("Correlation liver vs kidney")+ theme_linedraw() + geom_smooth(method = "lm")
p6 <- ggplot(tab.deg.art.logFC, aes(liver_logFC, vis_logFC)) + geom_point() + ggtitle("Correlation liver vs vis fat")+ theme_linedraw() + geom_smooth(method = "lm")
p7 <- ggplot(tab.deg.art.logFC, aes(liver_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation liver vs sc fat")+ theme_linedraw() + geom_smooth(method = "lm")
print(( p4 | p5 ) /
        ( p6 | p7))# 5.8x6
dev.off()

pdf(paste0("corr-art-kidney.pdf"),width=8.3,height=3,paper='special')
p5 <- ggplot(tab.deg.art.logFC, aes(kidney_logFC, kidney_logFC)) + geom_point() + ggtitle("Correlation kidney vs kidney")+ theme_linedraw() + geom_smooth(method = "lm")
p6 <- ggplot(tab.deg.art.logFC, aes(kidney_logFC, vis_logFC)) + geom_point() + ggtitle("Correlation kidney vs vis fat")+ theme_linedraw() + geom_smooth(method = "lm")
p7 <- ggplot(tab.deg.art.logFC, aes(kidney_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation kidney vs sc fat")+ theme_linedraw() + geom_smooth(method = "lm")
print(( p5 |p6 | p7)) # 8.3x3
dev.off()

pdf(paste0("corr-art-fat.pdf"),width=8.3,height=3,paper='special')
p5 <- ggplot(tab.deg.art.logFC, aes(vis_logFC, vis_logFC)) + geom_point() + ggtitle("Correlation vis vs vis fat")+ theme_linedraw() + geom_smooth(method = "lm")
p6 <- ggplot(tab.deg.art.logFC, aes(vis_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation vis vs sc fat")+ theme_linedraw() + geom_smooth(method = "lm")
p7 <- ggplot(tab.deg.art.logFC, aes(sc_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation sc vs sc fat")+ theme_linedraw() + geom_smooth(method = "lm")
print((p5 | p6 | p7)) # 8.3x3
dev.off()

# all clusters
m <- cor(tab.deg.art.logFC)
pdf(paste0("corr-all-artery-signif.pdf"),width=10,height=10,paper='special')
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black",
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))
dev.off()


####
# 3. Veins

deg.vein = list()
deg.i = NULL
for (i in organ) {
  print(i)
  deg.i = read.table(paste("deg_gen/",i,"_EC-ven_Western_vs_chow.txt", sep = "", collapse = NULL), sep = '\t', header = T)
  # sort rows by gene name
  deg.i = deg.i[order(rownames(deg.i),decreasing = F),]
  head(deg.i)
  deg.vein[[i]] = deg.i
}

deg.vein$liver

# combine a table of logFC
tab.deg.ven = merge(deg.vein$brain[,c(2,5)], deg.vein$lung[,c(2,5)], by="row.names", all=TRUE)
colnames(tab.deg.ven) = c("gene","brain_logFC","brain_Q-val","lung_logFC","lung_Q-val")
tab.deg.ven = merge(tab.deg.ven, deg.vein$heart[,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[6:7] = c("heart_logFC","heart_Q-val")
tab.deg.ven = merge(tab.deg.ven, deg.vein$liver[,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[8:9] = c("liver_logFC","liver_Q-val")
tab.deg.ven = merge(tab.deg.ven, deg.vein$kidney[,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[10:11] = c("kidney_logFC","kidney_Q-val")
tab.deg.ven = merge(tab.deg.ven, deg.vein$vis[,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[12:13] = c("vis_logFC","vis_Q-val")
tab.deg.ven = merge(tab.deg.ven, deg.vein$sc[,c(2,5)], by.x="gene",by.y = "row.names", all=TRUE)
colnames(tab.deg.ven)[14:15] = c("sc_logFC","sc_Q-val")
head(tab.deg.ven)

rownames(tab.deg.ven) = tab.deg.ven$gene
tab.deg.ven = tab.deg.ven[,-1]
nrow(tab.deg.ven) # 9146 genes
write.table(tab.deg.ven,"table_deg_vein_all_organs.txt", sep='\t')

tab.deg.ven.logFC = tab.deg.ven[,seq(1,14,2)] # 
# heatmap of the logFC WD/chow
b <- tab.deg.ven.logFC[rownames(tab.deg.ven.logFC) %in% nfkb,]
# replase NA with 0
b[is.na(b)] = 0
b

# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.15)|(rowMins(as.matrix(b),value = T)<=(-0.15)),]

# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5

mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = "NFkB targets in veins, logFC(WD/chow)", 
          Rowv= T,Colv= T, 
          ColSideColors= my_color_palette, # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          dendrogram="both", cexRow=0.8, cexCol=0.8) 


# heatmap of the logFC WD/chow
b <- tab.deg.ven.logFC[rownames(tab.deg.ven.logFC) %in% inflamm,]
# replase NA with 0
b[is.na(b)] = 0
b

# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.2)|(rowMins(as.matrix(b),value = T)<=(-0.2)),]

# make logFC max 0.5
b[b>0.5] = 0.5
b[b<(-0.5)] = -0.5

mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = "GO inflammation in veins, logFC(WD/chow)", 
          Rowv= T,Colv= T, 
          ColSideColors= my_color_palette, # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          dendrogram="both", cexRow=0.8, cexCol=0.8) 

# heatmap of the logFC WD/chow
b <- tab.deg.ven.logFC[rownames(tab.deg.ven.logFC) %in% FA.oxid,] #FA.oxid glycolysis
# replase NA with 0
b[is.na(b)] = 0
b

# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.15)|(rowMins(as.matrix(b),value = T)<=(-0.15)),]

# make logFC max 0.5
b[b>0.5] = 0.5

mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = "FA oxid in veins, logFC(WD/chow)", #FA oxidation Glycolysis
          Rowv= T,Colv= T, 
          ColSideColors= my_color_palette, # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          dendrogram="both", cexRow=0.8, cexCol=0.8) 

# heatmap of the logFC WD/chow
b <- tab.deg.ven.logFC[rownames(tab.deg.ven.logFC) %in% FA.trans,] #FA.oxid glycolysis
# replase NA with 0
b[is.na(b)] = 0
b

# exclude -0.1<logFC<0.1 genes
b = b[(rowMaxs(as.matrix(b),value = T)>=0.15)|(rowMins(as.matrix(b),value = T)<=(-0.15)),]

# make logFC max 0.5
b[b>0.5] = 0.5

mi <- min(b,na.rm=TRUE) #l
mi
ma <- max(b,na.rm=TRUE) #l
ma

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = "FA transport in veins, logFC(WD/chow)", #FA oxidation Glycolysis
          Rowv= T,Colv= T, 
          ColSideColors= my_color_palette, # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          dendrogram="both", cexRow=0.8, cexCol=0.8) 

# signif DEGs
sign.genes = rownames(tab.deg.ven[(tab.deg.ven$brain_Q.val<0.05)|(tab.deg.ven$lung_Q.val<0.05)|(tab.deg.ven$heart_Q.val<0.05)|(tab.deg.ven$liver_Q.val<0.05)|(tab.deg.ven$kidney_Q.val<0.05)|(tab.deg.ven$vis_Q.val<0.05)|(tab.deg.ven$sc_Q.val<0.05),])
sg = rownames(tab.deg.ven)[rownames(tab.deg.ven)%in%sign.genes]
tab.deg.ven.sign = tab.deg.ven[sg,]
nrow(tab.deg.ven.sign) # 648 genes with signif Q-val
head(tab.deg.ven.sign)
write.table(tab.deg.ven.sign, "vein_deg_sign.txt", sep = '\t')

# keep only logFC
tab.deg.ven.logFC = tab.deg.ven.sign[,seq(1,14,2)] # tab.deg.ven[,seq(1,14,2)]
head(tab.deg.ven.logFC)

# replace NA with 0s
tab.deg.ven.logFC[is.na(tab.deg.ven.logFC)] = 0
head(tab.deg.ven.logFC)

# remove rows with all 0s
tab.deg.ven.logFC = tab.deg.ven.logFC[!rowSums(tab.deg.ven.logFC)==0,]
head(tab.deg.ven.logFC)
nrow(tab.deg.ven.logFC) # 648 genes with logFC
# plot correlations

pdf(paste0("corr-ven-brain.pdf"),width=11.5,height=6,paper='special')
p1 <- ggplot(tab.deg.ven.logFC, aes(brain_logFC, brain_logFC)) + geom_point() + ggtitle("Correlation brain vs brain") + theme_linedraw() + geom_smooth(method = "lm")
p2 <- ggplot(tab.deg.ven.logFC, aes(brain_logFC, heart_logFC)) + geom_point() + ggtitle("Correlation brain vs heart") + theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.ven.logFC, aes(brain_logFC, lung_logFC)) + geom_point() + ggtitle("Correlation brain vs lung") + theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.ven.logFC, aes(brain_logFC, liver_logFC)) + geom_point() + ggtitle("Correlation brain vs liver") + theme_linedraw() + geom_smooth(method = "lm")
p5 <- ggplot(tab.deg.ven.logFC, aes(brain_logFC, kidney_logFC)) + geom_point() + ggtitle("Correlation brain vs kidney") + theme_linedraw() + geom_smooth(method = "lm")
p6 <- ggplot(tab.deg.ven.logFC, aes(brain_logFC, vis_logFC)) + geom_point() + ggtitle("Correlation brain vs vis fat") + theme_linedraw() + geom_smooth(method = "lm")
p7 <- ggplot(tab.deg.ven.logFC, aes(brain_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation brain vs sc fat") + theme_linedraw() + geom_smooth(method = "lm")
print((p1 | p2 | p3 |p4) /
        (p5 | p6 | p7| p7)) # 11.5x6
dev.off()

pdf(paste0("corr-ven-heart.pdf"),width=8.5,height=6,paper='special')
p2 <- ggplot(tab.deg.ven.logFC, aes(heart_logFC, heart_logFC)) + geom_point() + ggtitle("Correlation heart vs heart")+ theme_linedraw() + geom_smooth(method = "lm")
p3 <- ggplot(tab.deg.ven.logFC, aes(heart_logFC, lung_logFC)) + geom_point() + ggtitle("Correlation heart vs lung")+ theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.ven.logFC, aes(heart_logFC, liver_logFC)) + geom_point() + ggtitle("Correlation heart vs liver")+ theme_linedraw() + geom_smooth(method = "lm")
p5 <- ggplot(tab.deg.ven.logFC, aes(heart_logFC, kidney_logFC)) + geom_point() + ggtitle("Correlation heart vs kidney")+ theme_linedraw() + geom_smooth(method = "lm")
p6 <- ggplot(tab.deg.ven.logFC, aes(heart_logFC, vis_logFC)) + geom_point() + ggtitle("Correlation heart vs vis fat")+ theme_linedraw() + geom_smooth(method = "lm")
p7 <- ggplot(tab.deg.ven.logFC, aes(heart_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation heart vs sc fat")+ theme_linedraw() + geom_smooth(method = "lm")
print(( p2 | p3 |p4) /
        (p5 | p6 | p7)) # 8.5x6
dev.off()

pdf(paste0("corr-ven-lung.pdf"),width=8.5,height=6,paper='special')
p3 <- ggplot(tab.deg.ven.logFC, aes(lung_logFC, lung_logFC)) + geom_point() + ggtitle("Correlation lung vs lung")+ theme_linedraw() + geom_smooth(method = "lm")
p4 <- ggplot(tab.deg.ven.logFC, aes(lung_logFC, liver_logFC)) + geom_point() + ggtitle("Correlation lung vs liver")+ theme_linedraw() + geom_smooth(method = "lm")
p5 <- ggplot(tab.deg.ven.logFC, aes(lung_logFC, kidney_logFC)) + geom_point() + ggtitle("Correlation lung vs kidney")+ theme_linedraw() + geom_smooth(method = "lm")
p6 <- ggplot(tab.deg.ven.logFC, aes(lung_logFC, vis_logFC)) + geom_point() + ggtitle("Correlation lung vs vis fat")+ theme_linedraw() + geom_smooth(method = "lm")
p7 <- ggplot(tab.deg.ven.logFC, aes(lung_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation lung vs sc fat")+ theme_linedraw() + geom_smooth(method = "lm")
print(( p3 | p3 |p4) /
        (p5 | p6 | p7)) # 8.5x6
dev.off()

pdf(paste0("corr-ven-liver.pdf"),width=5.8,height=6,paper='special')
p4 <- ggplot(tab.deg.ven.logFC, aes(liver_logFC, liver_logFC)) + geom_point() + ggtitle("Correlation liver vs liver")+ theme_linedraw() + geom_smooth(method = "lm")
p5 <- ggplot(tab.deg.ven.logFC, aes(liver_logFC, kidney_logFC)) + geom_point() + ggtitle("Correlation liver vs kidney")+ theme_linedraw() + geom_smooth(method = "lm")
p6 <- ggplot(tab.deg.ven.logFC, aes(liver_logFC, vis_logFC)) + geom_point() + ggtitle("Correlation liver vs vis fat")+ theme_linedraw() + geom_smooth(method = "lm")
p7 <- ggplot(tab.deg.ven.logFC, aes(liver_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation liver vs sc fat")+ theme_linedraw() + geom_smooth(method = "lm")
print(( p4 | p5 ) /
        ( p6 | p7))# 5.8x6
dev.off()

pdf(paste0("corr-ven-kidney.pdf"),width=8.3,height=3,paper='special')
p5 <- ggplot(tab.deg.ven.logFC, aes(kidney_logFC, kidney_logFC)) + geom_point() + ggtitle("Correlation kidney vs kidney")+ theme_linedraw() + geom_smooth(method = "lm")
p6 <- ggplot(tab.deg.ven.logFC, aes(kidney_logFC, vis_logFC)) + geom_point() + ggtitle("Correlation kidney vs vis fat")+ theme_linedraw() + geom_smooth(method = "lm")
p7 <- ggplot(tab.deg.ven.logFC, aes(kidney_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation kidney vs sc fat")+ theme_linedraw() + geom_smooth(method = "lm")
print(( p5 |p6 | p7)) # 8.3x3
dev.off()

pdf(paste0("corr-ven-fat.pdf"),width=8.3,height=3,paper='special')
p5 <- ggplot(tab.deg.ven.logFC, aes(vis_logFC, vis_logFC)) + geom_point() + ggtitle("Correlation vis vs vis fat")+ theme_linedraw() + geom_smooth(method = "lm")
p6 <- ggplot(tab.deg.ven.logFC, aes(vis_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation vis vs sc fat")+ theme_linedraw() + geom_smooth(method = "lm")
p7 <- ggplot(tab.deg.ven.logFC, aes(sc_logFC, sc_logFC)) + geom_point() + ggtitle("Correlation sc vs sc fat")+ theme_linedraw() + geom_smooth(method = "lm")
print((p5 | p6 | p7)) # 8.3x3
dev.off()

# all clusters
m <- cor(tab.deg.ven.logFC)
pdf(paste0("corr-all-vein-signif.pdf"),width=10,height=10,paper='special')
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black",
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))
dev.off()

######
# Il6
b <- t(tab.deg.art.logFC[rownames(tab.deg.art.logFC)=="Il6",])
b = cbind(b,t(tab.deg.cap.logFC[rownames(tab.deg.cap.logFC)=="Il6",]))
b = cbind(b,t(tab.deg.ven.logFC[rownames(tab.deg.ven.logFC)=="Il6",]))
colnames(b) = c("Il6_art","Il6_ven")
b

# Junb
b <- t(tab.deg.art.logFC[rownames(tab.deg.art.logFC)=="Junb",])
b = cbind(b,t(tab.deg.cap.logFC[rownames(tab.deg.cap.logFC)=="Junb",]))
b = cbind(b,t(tab.deg.ven.logFC[rownames(tab.deg.ven.logFC)=="Junb",]))
colnames(b) = c("Junb_art","Junb_cap","Junb_ven")
b

ColorRamp <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100)) # 'RdGy' 'RdYlBu'
ColorLevels <- seq(mi, ma, length=length(ColorRamp))
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = "Junb, logFC(WD/chow)", 
          Rowv= T,Colv= T, 
          RowSideColors= my_color_palette, # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          dendrogram="both", cexRow=0.8, cexCol=0.8) 


###########################################################
# create tables for DEG logFC for all organs
folder = "~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/organs-3rd-processing/"
organ
dir = NULL
tabs = NULL
i=5
  organ[i]
  dir = paste0(folder,organ[i],"/deg/")
  tabs = list.files(dir)
  tabs
  deg.org = list()
  deg.i = NULL
  for (j in seq(1,length(tabs),1)) {
    print(tabs[j])
    deg.j = read.table(paste0(dir,"",tabs[j]), sep = '\t', header = T)
    # sort rows by gene name
    deg.j = deg.j[order(rownames(deg.j),decreasing = F),c(1,2,5)] # p-val, logFC, q-val
    head(deg.j)
    deg.org[[j]] = deg.j
  }
 

# merge
cl = c("EC-cap","EC-capA","EC-AP1","EC-art","EC-fenestr","EC-Hb","EC-platelet","EC-ven") # brain
cl = c("EC-anti-pathogen","EC-Aqp5a","EC-Aqp5b","EC-art","EC-cap1","EC-cap2",         
"EC-lymph","EC-platelet","EC-pulm1bc","EC-pulm1d","EC-pulm2a","EC-pulm2b","EC-ven","Prolif" ) # lung
cl = c("EC-ang","EC-AP1","EC-art","EC-arteriole","EC-cap1","EC-cap2","EC-Hb","EC-IFN",
       "EC-lymph","EC-ven","Prolif") # heart
cl = c("EC-cap1","EC-cap2","EC-cap3","EC-liver1","EC-liver2","EC-art","EC-ven") # liver
cl = c("EC-ang","EC-AP1","EC-Aqp1-arteriole","EC-art",           
       "EC-cap1","EC-cap2","EC-cap3","EC-glomeruli",     
       "EC-kidney1","EC-kidney2","EC-lymph","EC-ven","EC-venule") # kidney
cl = c("EC-ang","EC-antigen","EC-art","EC-cap1","EC-cap2","EC-ven","EC-venule", "EC-vis",
       "lymph1","lymph2","lymph3","lymph4","Prolif") # vis
cl = c("EC-ang","EC-antigen","EC-cap1","EC-cap2","EC-cap3","EC-cap4","EC-fenestr","EC-Hb","lymph",
       "EC-ven","EC-venule","Prolif","EC-art") # sc
length(cl)

deg.org.merge = merge(deg.org[[1]], deg.org[[2]], by="row.names", all=TRUE)
colnames(deg.org.merge) = c("gene",paste0(cl[1],"_p-val"),paste0(cl[1],"_logFC"),paste0(cl[1],"_Q-val"),
                          paste0(cl[2],"_p-val"),paste0(cl[2],"_logFC"),paste0(cl[2],"_Q-val"))
head(deg.org.merge)

deg.org.merge = merge(deg.org.merge, deg.org[[3]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[8:10] = c(paste0(cl[3],"_p-val"),paste0(cl[3],"_logFC"),paste0(cl[3],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[4]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[11:13] = c(paste0(cl[4],"_p-val"),paste0(cl[4],"_logFC"),paste0(cl[4],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[5]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[14:16] = c(paste0(cl[5],"_p-val"),paste0(cl[5],"_logFC"),paste0(cl[5],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[6]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[17:19] = c(paste0(cl[6],"_p-val"),paste0(cl[6],"_logFC"),paste0(cl[6],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[7]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[20:22] = c(paste0(cl[7],"_p-val"),paste0(cl[7],"_logFC"),paste0(cl[7],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[8]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[23:25] = c(paste0(cl[8],"_p-val"),paste0(cl[8],"_logFC"),paste0(cl[8],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[9]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[26:28] = c(paste0(cl[9],"_p-val"),paste0(cl[9],"_logFC"),paste0(cl[9],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[10]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[29:31] = c(paste0(cl[10],"_p-val"),paste0(cl[10],"_logFC"),paste0(cl[10],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[11]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[32:34] = c(paste0(cl[11],"_p-val"),paste0(cl[11],"_logFC"),paste0(cl[11],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[12]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[35:37] = c(paste0(cl[12],"_p-val"),paste0(cl[12],"_logFC"),paste0(cl[12],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[13]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[38:40] = c(paste0(cl[13],"_p-val"),paste0(cl[13],"_logFC"),paste0(cl[13],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[14]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[41:43] = c(paste0(cl[14],"_p-val"),paste0(cl[14],"_logFC"),paste0(cl[14],"_Q-val"))

head(deg.org.merge)
write.table(deg.org.merge, paste0("deg_clust_logFC_",organ[i],".txt"), sep = '\t')

##########
# Merged logFC tables for all clusters

deg.brain.merge = read.table("deg_organ_log/deg_clust_logFC_brain.txt",sep='\t', header = T)
rownames(deg.brain.merge) = deg.brain.merge$gene
deg.brain.merge = deg.brain.merge[,-1]

deg.lung.merge = read.table("deg_organ_log/deg_clust_logFC_lung.txt",sep='\t', header = T)
rownames(deg.lung.merge) = deg.lung.merge$gene
deg.lung.merge = deg.lung.merge[,-1]

deg.heart.merge = read.table("deg_organ_log/deg_clust_logFC_heart.txt",sep='\t', header = T)
rownames(deg.heart.merge) = deg.heart.merge$gene
deg.heart.merge = deg.heart.merge[,-1]

deg.liver.merge = read.table("deg_organ_log/deg_clust_logFC_liver.txt",sep='\t', header = T)
rownames(deg.liver.merge) = deg.liver.merge$gene
deg.liver.merge = deg.liver.merge[,-1]

#deg.kidney.merge = read.table("deg_organ_log/table_deg_kidney_clust_all_WD3m.txt",sep='\t', header = T)
deg.kidney.merge = read.table("deg_organ_log/deg_clust_logFC_kidney.txt",sep='\t', header = T)
rownames(deg.kidney.merge) = deg.kidney.merge$gene
deg.kidney.merge = deg.kidney.merge[,-1]

deg.vis.merge = read.table("deg_organ_log/deg_clust_logFC_vis.txt",sep='\t', header = T)
rownames(deg.vis.merge) = deg.vis.merge$gene
deg.vis.merge = deg.vis.merge[,-1]

deg.sc.merge = read.table("deg_organ_log/deg_clust_logFC_sc.txt",sep='\t', header = T)
rownames(deg.sc.merge) = deg.sc.merge$gene
deg.sc.merge = deg.sc.merge[,-1]

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

klfs = c("Klf2","Klf3","Klf4","Klf6","Klf7","Klf9","Klf10","Klf12","Klf13","Klf16")

### for lungs
h2 = c("H2-Aa","H2-Ab1","H2-D1","H2-DMa","H2-DMb1","H2-Eb1","H2-K1","H2-Ke6","H2-M3","H2-Q4","H2-Q6","H2-Q7","H2-T22","H2-T23")
rps = c("Rps9","Rpl3","Rps8","Rpl11","Rpl7","Rpl7a","Rps15a","Rpl18a","Rps19","Rpl24","Rps20","Rpl26","Rpl28",
        "Rps21","Rps13")


resp.chain = read.table("../resp_chain.txt", sep='\t')
compl.col = c("cyan2","coral","deeppink1","darkturquoise","red")
names(compl.col) = unique(resp.chain$V2)
resp.chain$V3 = 1
resp.chain$V3 = compl.col[match(resp.chain$V2, names(compl.col))]

ap1 = c("Jun","Junb","Jund","Fos","Fosb","Egr1")

mts = rownames(deg.brain.merge)[grep("mt-",rownames(deg.brain.merge))]

slc.brain = c("Slc16a1","Slc22a8","Slc38a5","Slc2a1","Slc40a1","Slc1a1","Slc6a6","Slc39a10",
              "Slc7a1","Slc9a8","Slc25a32","Slc12a6","Slc39a7","Slc27a3","Slc7a5","Slc3a2",
              "Slc25a37","Slc25a30","Slc24a5","Slc43a2","Slc29a1","Slc52a2",
              "Slc25a24","Slc9a1","Slc6a8","Slc30a1","Slc4a7","Slc35f5")
fat = c("Apoe","Apoc1","Abca1","Ttr","Alb","Apoa2")

# heatmap of the logFC WD/chow
b <- deg.liver.merge[rownames(deg.liver.merge) %in% fat,] # coagul junctions
write.table(b, "kidney_resp_chain_logFC_3m.txt")

b = b[,seq(2,ncol(b),3)] # leave only logFC
#b = b[,c(9,10)]
# replase NA with 0
b[is.na(b)] = 0
b
#b <- b[order(match(rownames(b),resp.chain$V1)),] #metabolism

# exclude -0.1<logFC<0.1 genes
#b = b[(rowMaxs(as.matrix(b),value = T)>=0.1)|(rowMins(as.matrix(b),value = T)<=(-0.1)),]

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
          main = "fat in liver, logFC(WD/chow)",  # Coagulation Junctions
          Rowv= T,Colv= T, 
          #ColSideColors= my_color_palette, # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(7,10),     # widens margins around plot
          keysize=0.8,
          dendrogram="both", cexRow=0.7, cexCol=0.8) 

clust.avg3 = read.table("../organs-4th-processing/clust_avgs/cluster_avg_expr_kidney.txt",sep=' ',header=T)
head(clust.avg3)

slc.all = rownames(clust.avg3)[grep("Slc",rownames(clust.avg3))]
  
b <- clust.avg3[rownames(clust.avg3) %in% slc.all,] 
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
heatmap.2(as.matrix(z), na.rm = TRUE, col= ColorRamp, 
          main = paste0("kidney SLCs, z-score of natural log expr"), #FA transport Platelet aggr Hemostasis
          Rowv= T,Colv= T, 
          #ColSideColors= my_color_palette, # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="both", cexRow=0.2, cexCol=0.8) 

##########################
# Correlation of clusters within each organ

deg.clust = list()
for (i in organ) {
  print(i)
  tbls = list.files(paste("../organs-3rd-processing/",i,"/deg/", sep = "", collapse = NULL))
  print(tbls)
  deg.j = NULL
  deg.i = data.frame()
  for (j in seq(1,length(tbls),1)) {
    deg.j = read.table(paste0("../organs-3rd-processing/",i,"/deg/",tbls[j]), sep = '\t', header = T)
    # sort rows by gene name
    deg.j = deg.j[order(rownames(deg.j),decreasing = F),]
    deg.j = deg.j[,c(2,5)] # logFC, q-val
    head(deg.j)
    deg.i = merge(deg.i,deg.j,  by="row.names", all = T)
    rownames(deg.i) = deg.i[,"Row.names"]
    deg.i[,"Row.names"] = NULL
  }
  head(deg.i)
  colnames(deg.i)[seq(1,2*length(tbls),2)] = tbls
  head(deg.i)
  deg.clust[[i]] = deg.i
  write.table(deg.clust[[i]], paste0(i, "_deg_all_clust.txt"), sep = '\t')
}

setwd("corrs_organ")
# read deg list
deg.clust = list()
for (i in organ) {
  deg.clust[[i]] = read.table(paste0(i, "_deg_all_clust.txt"), sep = '\t')
}

# keep only genes with at least one signif. q-val
i = "kidney"
#for (i in organ) {
  print(i)
  #deg.clust[[i]] = deg.clust[[i]][,-c(3,4)] # removing antigen=immune
  q.tab = deg.clust[[i]][,seq(2,ncol(deg.clust[[i]]),2)] # cols with q-val
  ncol(deg.clust[[i]])/2
  
  sign.genes = rownames(q.tab[(q.tab[,1]<0.05)|(q.tab[,2]<0.05)|(q.tab[,3]<0.05)|(q.tab[,4]<0.05)|(q.tab[,5]<0.05)|(q.tab[,6]<0.05)|(q.tab[,7]<0.05)|(q.tab[,8]<0.05)|(q.tab[,9]<0.05)|(q.tab[,10]<0.05)|(q.tab[,11]<0.05)|(q.tab[,12]<0.05)|(q.tab[,13]<0.05),]) #|(q.tab[,8]<0.05)|(q.tab[,9]<0.05)|(q.tab[,10]<0.05)|(q.tab[,11]<0.05)|(q.tab[,12]<0.05)|(q.tab[,13]<0.05)|(q.tab[,14]<0.05)
  sg = rownames(q.tab)[rownames(q.tab)%in%sign.genes]
  q.tab = q.tab[sg,] 
  nrow(q.tab) # 217 genes with signif Q-val
  head(q.tab)
  deg.clust[[i]] = deg.clust[[i]][rownames(q.tab),]
  head(deg.clust[[i]])
  write.table(deg.clust[[i]], paste0(i, "_deg_all_clust_signif.txt"), sep = '\t')
#}

  # keep only logFC colums
  organ = c("brain", "lung", "heart", "liver", "kidney","vis","sc")
  #deg.clust.sign = deg.clust
  deg.clust.sign = list()
  for (i in organ) {
    print(i)
    deg.clust.sign[[i]] = read.table(paste0(i, "_deg_all_clust_signif.txt"), sep = '\t') #"corrs_organ/",
    k = seq(1,ncol(deg.clust.sign[[i]]),2) # cols w/o q-val
    deg.clust.sign[[i]] = deg.clust.sign[[i]][,k]
    head(deg.clust.sign[[i]])
    # replace NA with 0
    deg.clust.sign[[i]][is.na(deg.clust.sign[[i]])] = 0
    
    # correlation
    # all clusters
    m <- cor(deg.clust.sign[[i]])
    res1 <- cor.mtest(deg.clust.sign[[i]], conf.level = .95) # significance
    hclust_avg <- hclust(dist(m), method = 'average')
   
    pdf(paste0("corr-all-clust-",i,"-hclust-tree.pdf"),width=10,height=7,paper='special')
    plot(hclust_avg)
    dev.off()
    
    pdf(paste0("corr-all-clust-",i,"-signif-hclust.pdf"),width=10,height=10,paper='special')
    #corrplot(m, method = "color", type = "lower", order = "AOE",col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))
    corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black", 
                  order = "hclust", hclust.method = "average", tl.cex = 0.6,
                 upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))
    dev.off()
  }

head(deg.clust.sign$brain)

### no cutoff for significance
deg.clust = list()
deg.clust.sign = list()
i="kidney"
for (i in organ) {
  print(i)
  deg.clust[[i]] = read.table(paste0(i, "_deg_all_clust.txt"), sep = '\t') #"corrs_organ/",
  head(deg.clust[[i]])
  #deg.clust[[i]] = deg.clust[[i]][,-c(3,4)] # removing antigen=immune
  lfc = deg.clust[[i]][,seq(1,ncol(deg.clust[[i]]),2)] # cols with logFC
  ncol(deg.clust[[i]])/2
  
  sign.genes = rownames(lfc[(abs(lfc[,1])>=0.1)|(abs(lfc[,2])>=0.1)|(abs(lfc[,3])>=0.1)|(abs(lfc[,4])>=0.1)|
                              (abs(lfc[,5])>=0.1)|(abs(lfc[,6])>=0.1)|(abs(lfc[,7])>=0.1)
                              ,]) #|(abs(lfc[,8])>=0.1)|(abs(lfc[,9])>=0.1)|(abs(lfc[,10])>=0.1)|(abs(lfc[,11])>=0.1)|(abs(lfc[,12])>=0.1)|(abs(lfc[,13])>=0.1)|(abs(lfc[,14])>=0.1)
  sg = rownames(lfc)[rownames(lfc)%in%sign.genes]
  lfc = lfc[sg,] 
  nrow(lfc) # 9850 genes with logFC 0.1
  head(lfc)
  deg.clust.sign[[i]] = deg.clust[[i]][rownames(lfc),]
  
  k = seq(1,ncol(deg.clust.sign[[i]]),2) # cols w/o q-val
  deg.clust.sign[[i]] = deg.clust.sign[[i]][,k]
  head(deg.clust.sign[[i]])
  # replace NA with 0
  deg.clust.sign[[i]][is.na(deg.clust.sign[[i]])] = 0
  
  # correlation
  # all clusters
  m <- cor(deg.clust.sign[[i]])
  res1 <- cor.mtest(deg.clust.sign[[i]], conf.level = .95) # significance
  hclust_avg <- hclust(dist(m), method = 'average')
  
  pdf(paste0("corr-all-clust-",i,"-lfc-cutoff-hclust-tree.pdf"),width=10,height=7,paper='special')
  plot(hclust_avg)
  dev.off()
  
  pdf(paste0("corr-all-clust-",i,"-lfc-cutoff-hclust.pdf"),width=20,height=20,paper='special')
  #corrplot(m, method = "color", type = "lower", order = "AOE",col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))
  corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black", 
                 order = "hclust", hclust.method = "average", tl.cex = 0.6,
                 upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))
  dev.off()
}

##########################
library(jtools) # for summ()
angiog = read.table("angiogen_organ.txt", sep='\t', header = T)
angiog$N = angiog$N_cells_chow + angiog$N_cells_WD
sign.angiog <- lm(chow ~ Western * N , data = angiog)
summ(sign.angiog)

##########################
# heatmaps for GOs
tab.deg.cap.sign = read.table("corrs/capillary_deg_sign.txt", sep='\t')
tab.deg.cap.logFC = tab.deg.cap.sign[,c(1,3,5,7,9,11,13)]

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

# lipid mobilization
lip.mob = c("Saa1","Saa2","Saa4","Pon1","Pon3","Lcat","Ttr","Apoe", "Apoa1","Apoa2", "Apoa4","Apoc1","Apoc2",
            "Apob","Apod","Serpina1","Alb","Abca1")

# GO full focal adhesion
go.focal.adh = read.table("../GO_0005925_focal_adhesion.txt", sep='\t')
go.focal.adh = unique(go.focal.adh$V1)

# heatmap of the logFC WD/chow
b <- tab.deg.cap.logFC[rownames(tab.deg.cap.logFC) %in% lip.mob,] #focal ppar integrin  c(lip.mob,lipid)

b <- tab.deg.all.logFC[rownames(tab.deg.all.logFC) %in% c(lip.mob,lipid),] #focal ppar integrin
# replase NA with 0
b[is.na(b)] = 0
b

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
heatmap.2(as.matrix(b), na.rm = TRUE, col= ColorRamp, main = "Lipid mobilization in capillaries, logFC(WD/chow)", #PPAR signaling
          Rowv= T,Colv= F, 
          ColSideColors= my_color_palette, # rep(my_color_palette, each=3)  my_color_palette
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,10),     # widens margins around plot
          keysize=0.75,
          dendrogram="row", cexRow=0.8, cexCol=0.8) 

####################################
####################################

###########################################################
# create tables for DEG logFC for all organs
folder = "~/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/analysis/organs-3rd-processing/"
organ = c("brain","heart","lung","liver","kidney","vis","sc")
dir = NULL
tabs = NULL
i=3
organ[i]
dir = paste0(folder,organ[i],"/deg/")
tabs = list.files(dir)
tabs
deg.org = list()
deg.i = NULL
for (j in seq(1,length(tabs),1)) {
  print(tabs[j])
  deg.j = read.table(paste0(dir,"",tabs[j]), sep = '\t', header = T)
  # sort rows by gene name
  deg.j = deg.j[order(rownames(deg.j),decreasing = F),c(1,2,5)] # p-val, logFC, q-val
  head(deg.j)
  deg.org[[j]] = deg.j
}


# merge
cl = c("EC-cap","EC-capA","EC-AP1","EC-art","EC-fenestr","EC-Hb","EC-platelet","EC-ven") # brain
cl = c("aEC","EC-Aqp5a","EC-Aqp5b","EC-art","EC-cap1","EC-cap2",         
       "EC-lymph","EC-platelet","EC-pneumocyte","pulmEC-a","pulmEC-b","EC-ven","Prolif" ) # lung
cl = c("EC-ang","EC-AP1","EC-art","EC-arteriole","EC-cap1","EC-cap2","EC-Hb","EC-IFN",
       "EC-lymph","EC-ven","Prolif") # heart
cl = c("EC-cap1","EC-cap2","EC-cap3","EC-liver1","EC-liver2","EC-art","EC-ven") # liver
cl = c("EC-ang","EC-AP1","EC-Aqp1-arteriole","EC-art",           
       "EC-cap1","EC-cap2","EC-cap3","EC-glomeruli",     
       "mEC1","EC-lymph","EC-ven","EC-venule") # kidney
cl = c("EC-ang","EC-art","EC-cap1","EC-cap2","EC-lymph1","EC-lymph2","EC-ven","EC-venule", "EC-vis",
       "Prolif") # vis
cl = c("EC-ang","EC-antigen","EC-cap1","EC-cap2","EC-cap3","EC-cap4","EC-fenestr","EC-Hb","lymph",
       "EC-ven","EC-venule","Prolif","EC-art") # sc
length(cl)

deg.org.merge = merge(deg.org[[1]], deg.org[[2]], by="row.names", all=TRUE)
colnames(deg.org.merge) = c("gene",paste0(cl[1],"_p-val"),paste0(cl[1],"_log2FC"),paste0(cl[1],"_Q-val"),
                            paste0(cl[2],"_p-val"),paste0(cl[2],"_logFC"),paste0(cl[2],"_Q-val"))
head(deg.org.merge)

deg.org.merge = merge(deg.org.merge, deg.org[[3]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[8:10] = c(paste0(cl[3],"_p-val"),paste0(cl[3],"_logFC"),paste0(cl[3],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[4]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[11:13] = c(paste0(cl[4],"_p-val"),paste0(cl[4],"_logFC"),paste0(cl[4],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[5]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[14:16] = c(paste0(cl[5],"_p-val"),paste0(cl[5],"_logFC"),paste0(cl[5],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[6]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[17:19] = c(paste0(cl[6],"_p-val"),paste0(cl[6],"_logFC"),paste0(cl[6],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[7]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[20:22] = c(paste0(cl[7],"_p-val"),paste0(cl[7],"_logFC"),paste0(cl[7],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[8]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[23:25] = c(paste0(cl[8],"_p-val"),paste0(cl[8],"_logFC"),paste0(cl[8],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[9]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[26:28] = c(paste0(cl[9],"_p-val"),paste0(cl[9],"_logFC"),paste0(cl[9],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[10]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[29:31] = c(paste0(cl[10],"_p-val"),paste0(cl[10],"_logFC"),paste0(cl[10],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[11]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[32:34] = c(paste0(cl[11],"_p-val"),paste0(cl[11],"_logFC"),paste0(cl[11],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[12]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[35:37] = c(paste0(cl[12],"_p-val"),paste0(cl[12],"_logFC"),paste0(cl[12],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[13]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[38:40] = c(paste0(cl[13],"_p-val"),paste0(cl[13],"_logFC"),paste0(cl[13],"_Q-val"))

deg.org.merge = merge(deg.org.merge, deg.org[[14]], by.x="gene",by.y = "row.names", all=TRUE)
colnames(deg.org.merge)[41:43] = c(paste0(cl[14],"_p-val"),paste0(cl[14],"_logFC"),paste0(cl[14],"_Q-val"))

head(deg.org.merge)
ncol(deg.org.merge)
# change the base
deg.org.merge$`EC-lymph1_logFC` = deg.org.merge$`EC-lymph1_log2FC`/log(exp(1),2)
deg.org.merge$`EC-lymph2_logFC` = deg.org.merge$`EC-lymph2_log2FC`/log(exp(1),2)

deg.org.merge$`aEC_logFC` = deg.org.merge$`aEC_log2FC`/log(exp(1),2)

# replace columns for lymph
deg.org.merge[,15] = deg.org.merge[,32]
deg.org.merge[,18] = deg.org.merge[,33]
deg.org.merge = deg.org.merge[,-c(32,33)]
names(deg.org.merge)[15] = "EC-lymph1_logFC"
names(deg.org.merge)[18] = "EC-lymph2_logFC"
head(deg.org.merge)

deg.org.merge[,3] = deg.org.merge[,41]
deg.org.merge = deg.org.merge[,-c(41)]
names(deg.org.merge)[3] = "aEC_logFC"

head(deg.org.merge)

write.table(deg.org.merge, paste0("deg_clust_logFC_",organ[i],".txt"), sep = '\t')

##########

deg.brain.merge = read.table("deg_organ_log/deg_clust_logFC_brain.txt",sep='\t', header = T)
rownames(deg.brain.merge) = deg.brain.merge$gene
deg.brain.merge = deg.brain.merge[,-1]

deg.lung.merge = read.table("deg_organ_log/deg_clust_logFC_lung.txt",sep='\t', header = T)
rownames(deg.lung.merge) = deg.lung.merge$gene
deg.lung.merge = deg.lung.merge[,-1]

deg.heart.merge = read.table("deg_organ_log/deg_clust_logFC_heart.txt",sep='\t', header = T)
rownames(deg.heart.merge) = deg.heart.merge$gene
deg.heart.merge = deg.heart.merge[,-1]

deg.liver.merge = read.table("deg_organ_log/deg_clust_logFC_liver.txt",sep='\t', header = T)
rownames(deg.liver.merge) = deg.liver.merge$gene
deg.liver.merge = deg.liver.merge[,-1]

deg.kidney.merge = read.table("deg_organ_log/deg_clust_logFC_kidney.txt",sep='\t', header = T)
rownames(deg.kidney.merge) = deg.kidney.merge$gene
deg.kidney.merge = deg.kidney.merge[,-1]

deg.vis.merge = read.table("deg_organ_log/deg_clust_logFC_vis.txt",sep='\t', header = T)
rownames(deg.vis.merge) = deg.vis.merge$gene
deg.vis.merge = deg.vis.merge[,-1]

deg.sc.merge = read.table("deg_organ_log/deg_clust_logFC_sc.txt",sep='\t', header = T)
rownames(deg.sc.merge) = deg.sc.merge$gene
deg.sc.merge = deg.sc.merge[,-1]

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

### for lungs
h2 = c("H2-Aa","H2-Ab1","H2-D1","H2-DMa","H2-DMb1","H2-Eb1","H2-K1","H2-Ke6","H2-M3","H2-Q4","H2-Q6","H2-Q7","H2-Q10","H2-T22","H2-T23","H2-T24")
rps = c("Rps9","Rpl3","Rps8","Rpl11","Rpl7","Rpl7a","Rps15a","Rpl18a","Rps19","Rpl24","Rps20","Rpl26","Rpl28",
        "Rps21","Rps13")

# heatmap of the logFC WD/chow
b <- deg.lung.merge[rownames(deg.lung.merge) %in% h2,] # coagul junctions leuk.recruit rps
b = b[,seq(2,ncol(b),3)] # leave only logFC
# replase NA with 0
b[is.na(b)] = 0
b

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
          main = "Histocomatiblility genes in lung, logFC(WD/chow)",  # Coagulation Junctions Leukocyte recruitment Ribosomal
          Rowv= T,Colv= T, 
          #ColSideColors= my_color_palette, # rep(cell_col.both, each=2)
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(7,10),     # widens margins around plot
          keysize=0.8,
          dendrogram="both", cexRow=0.8, cexCol=0.8) 

#####
deg.vis.merge = deg.org.merge
rownames(deg.vis.merge) = deg.vis.merge$gene
deg.vis.merge = deg.vis.merge[,-1]
# leave logFC, filter abs(logFC)>0.1
deg.vis.merge = deg.vis.merge[,c(2,5,8,11,14,17,20,23,26,29)]
head(deg.vis.merge)

sign.genes = rownames(deg.vis.merge[(abs(deg.vis.merge[,1])>=0.1)|(abs(deg.vis.merge[,2])>=0.1)|(abs(deg.vis.merge[,3])>=0.1)|(abs(deg.vis.merge[,4])>=0.1)|
                            (abs(deg.vis.merge[,5])>=0.1)|(abs(deg.vis.merge[,6])>=0.1)|(abs(deg.vis.merge[,7])>=0.1)|(abs(deg.vis.merge[,8])>=0.1)|(abs(deg.vis.merge[,9])>=0.1)|(abs(deg.vis.merge[,10])>=0.1)
                          ,]) #|(abs(lfc[,11])>=0.1)|(abs(lfc[,12])>=0.1)|(abs(lfc[,13])>=0.1)|(abs(lfc[,14])>=0.1)
sg = rownames(deg.vis.merge)[rownames(deg.vis.merge)%in%sign.genes]
deg.vis.merge = deg.vis.merge[sg,] 
nrow(deg.vis.merge) # 9750 genes with logFC 0.1
head(deg.vis.merge)

# replace NA with 0
deg.vis.merge[is.na(deg.vis.merge)] = 0

# correlation
# all clusters
m <- cor(deg.vis.merge)
res1 <- cor.mtest(deg.vis.merge, conf.level = .95) # significance
hclust_avg <- hclust(dist(m), method = 'average')

pdf(paste0("corr-all-clust-vis-lfc-cutoff-hclust-tree-new.pdf"),width=10,height=7,paper='special')
plot(hclust_avg)
dev.off()

pdf(paste0("corr-all-clust-vis-lfc-cutoff-hclust-new.pdf"),width=20,height=20,paper='special')
#corrplot(m, method = "color", type = "lower", order = "AOE",col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black", 
               order = "hclust", hclust.method = "average", tl.cex = 0.6,
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))
dev.off()

#####
deg.lung.merge = deg.org.merge
rownames(deg.lung.merge) = deg.lung.merge$gene
deg.lung.merge = deg.lung.merge[,-1]
# leave logFC, filter abs(logFC)>0.1
deg.lung.merge = deg.lung.merge[,c(2,5,8,11,14,17,20,23,26,29,32,35,38)]
head(deg.lung.merge)

sign.genes = rownames(deg.lung.merge[(abs(deg.lung.merge[,1])>=0.1)|(abs(deg.lung.merge[,2])>=0.1)|(abs(deg.lung.merge[,3])>=0.1)|(abs(deg.lung.merge[,4])>=0.1)|
                                      (abs(deg.lung.merge[,5])>=0.1)|(abs(deg.lung.merge[,6])>=0.1)|(abs(deg.lung.merge[,7])>=0.1)|(abs(deg.lung.merge[,8])>=0.1)|(abs(deg.lung.merge[,9])>=0.1)|(abs(deg.lung.merge[,10])>=0.1)|
                                     (abs(deg.lung.merge[,11])>=0.1)|(abs(deg.lung.merge[,12])>=0.1)|(abs(deg.lung.merge[,13])>=0.1),]) #|(abs(lfc[,14])>=0.1)
sg = rownames(deg.lung.merge)[rownames(deg.lung.merge)%in%sign.genes]
deg.lung.merge = deg.lung.merge[sg,] 
nrow(deg.lung.merge) # 8953 genes with logFC 0.1
head(deg.lung.merge)

# replace NA with 0
deg.lung.merge[is.na(deg.lung.merge)] = 0

# correlation
# all clusters
m <- cor(deg.lung.merge)
res1 <- cor.mtest(deg.lung.merge, conf.level = .95) # significance
hclust_avg <- hclust(dist(m), method = 'average')

pdf(paste0("corr-all-clust-lung-lfc-cutoff-hclust-tree-new.pdf"),width=10,height=7,paper='special')
plot(hclust_avg)
dev.off()

pdf(paste0("corr-all-clust-lung-lfc-cutoff-hclust-new.pdf"),width=20,height=20,paper='special')
#corrplot(m, method = "color", type = "lower", order = "AOE",col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))
corrplot.mixed(m, lower = "number", upper = "square",lower.col = "black", 
               order = "hclust", hclust.method = "average", tl.cex = 0.6,
               upper.col = rev(colorRampPalette(brewer.pal(5, 'RdBu'))(20)))
dev.off()
