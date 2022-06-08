# Convert Seurat object into Anndata -> .h5ad for cellxgene
#devtools::install_github("cellgeni/sceasy")
#BiocManager::install(c("LoomExperiment", "SingleCellExperiment"))
#install.packages('reticulate')
#loompy <- reticulate::import('loompy')


library(Seurat)
library(sceasy)
Sys.setenv(RETICULATE_PYTHON="/opt/miniconda3/envs/sceasy/bin/python3")
library(reticulate)
reticulate::use_condaenv('/opt/miniconda3/envs/sceasy/')
#reticulate::use_virtualenv('/opt/miniconda3/envs/sceasy/')
reticulate::py_discover_config()

setwd("/Users/olga.bondareva/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_02_23_M_EC_HFD_extra_seq/cellxgene")
#Seurat to AnnData
obj.dir = c("/Users/olga.bondareva/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_06_09_M_HFD_Rev_3m/analysis/organs/")
names = c("brain","heart","lung","liver","kidney") #"brain","heart","lung","liver","kidney","vis","sc"
obj.dir = c("/Users/olga.bondareva/Desktop/HI-MAG_OLGA/SeqData/Cellranger/2021_06_09_M_HFD_Rev_3m/analysis/NEW_SEQ/organs/")
names = c("vis","sc")
k=0
for (i in names) {
  k=k+1
  print(names[k])
  obj = readRDS(paste0(obj.dir,names[k],"/",names[k],"_subset_3m_4m_6m.rds",sep="")) #subset
  obj@meta.data[,grep("RNA_snn", colnames(obj@meta.data))] = NULL
  obj@meta.data[,grep("DF.", colnames(obj@meta.data))] = NULL
  obj@meta.data[,grep("pANN", colnames(obj@meta.data))] = NULL
  obj@meta.data$organ = NULL
  obj@meta.data$doublet = NULL
  obj@meta.data$ecs = NULL
  obj@meta.data$cdh5 = NULL
  obj@meta.data$pecam1 = NULL
  obj@meta.data$celltype = NULL
  obj@meta.data$time = NULL
  obj@meta.data$diet = NULL
  obj@meta.data$clust.diet = NULL
  obj@meta.data$orig.ident = NULL
  obj@meta.data$gen_celltype = NULL
  obj@meta.data$gen.diet = NULL
  obj@meta.data$seurat_clusters = NULL
  obj@meta.data$clust.diet_time = NULL
  obj@meta.data[,paste0(names[k],"_celltype")] = NULL
  obj@meta.data$repl = NULL
  obj@meta.data$S.Score = NULL
  obj@meta.data$G2M.Score = NULL
  obj@meta.data$Phase = NULL
  obj@meta.data[,paste0(names[k],"_celltype2.diet")] = NULL
  obj@meta.data[,paste0(names[k],"_celltype2.diet_time")] = NULL
  obj@meta.data$migration1 = NULL
  obj@meta.data$EndoMT1 = NULL
  obj@meta.data$old.ident = NULL
  obj@meta.data$vis_celltype2 = NULL
  
  print(colnames(obj@meta.data))
  
  sceasy::convertFormat(obj, from="seurat", to="anndata",
                       outFile=paste0(names[k],'_3m_4m_6m.h5ad', sep=""))
  
}


obj.dir = "/Users/olga.bondareva/Desktop/HI-MAG OLGA/SeqData/Cellranger/2021_06_09_M_HFD_Rev_3m/analysis/combined/combined_rev3m_merged_organs_ECs.rds"
combined = readRDS(obj.dir)
combined@meta.data[,grep("RNA_snn", colnames(combined@meta.data))] = NULL
combined@meta.data$cdh5 = NULL
combined@meta.data$pecam1 = NULL
combined@meta.data$celltype = NULL
combined@meta.data$ecs = NULL
combined@meta.data$doublet = NULL
sceasy::convertFormat(combined, from="seurat", to="anndata",
                      outFile='combined_ECs_3m_4m_6m.h5ad')

