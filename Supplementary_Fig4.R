#=================================================================================================================#
# Prepare R session
#=================================================================================================================#

library(Seurat)
library(ggplot2)
library(ggpubr)
library(pagoda2)
library(ggrepel)
library(clustree)
library(DESeq2)
library(SeuratWrappers)
library(dplyr)
library(stringr)
library(sctransform)
library(plyr)
library(ggcorrplot)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(scCustomize)
library(viridis)
library(RColorBrewer)
library(scplotter)
library(EnhancedVolcano)
library(patchwork)


#-----------------------------------------------------------------------------------------------------------------#
# !!!DEFINE ALL EXPERIMENT-SPECFIC VARIABLES!!!
#-----------------------------------------------------------------------------------------------------------------#
setwd("file")
analysis <- sub(".R", "", stringr::str_split(rstudioapi::getSourceEditorContext()$path, "/", simplify = T)[,10])
analysisid <- unlist(lapply(analysis, function(x){strsplit(as.character(x), '\\_')[[1]][1]}))
an.desc <- unlist(lapply(analysis, function(x){strsplit(as.character(x), '\\_')[[1]][2]}))
Exp_ID <- stringr::str_split(rstudioapi::getSourceEditorContext()$path, "/", simplify = T)[,8]

dir.create(paste0("plots"))
dir.create(paste0("files"))
dir.create(paste0("plots/", analysisid,"_",an.desc))
dir.create(paste0("files/", analysisid,"_",an.desc))
dir.figs <- paste0("plots/", analysisid, "_", an.desc, "/", analysisid, "_")
dir.files <- paste0("files/", analysisid, "_", an.desc, "/", analysisid, "_")

#-----------------------------------------------------------------------------------------------------------------#
# read variable needed
#-----------------------------------------------------------------------------------------------------------------#
lm2.mg <- readRDS("files/an0013.filtered.merge.sample2.3.rds")
mv.mg <- readRDS("files/an0026.filtered.merge.sample2.3.rds")
Br.mg <- readRDS("files/an0042.filtered.merge.batch2.rds")

# create metadata without sample number
c.sme <- unlist(lapply(lm2.mg$orig.ident, function(x){strsplit(as.character(x), '\\-')[[1]][1]}))
hx.sme <- unlist(lapply(lm2.mg$orig.ident, function(x){strsplit(as.character(x), '\\-')[[1]][3]}))
h.sme <- paste(c.sme,hx.sme,sep = "-")
lm2.mg$cell.type <- h.sme

#=================================================================================================================#
# Supplementary Fig4A
#=================================================================================================================#
#LM2
colo <- brewer.pal(9, "Purples")
FP.MKI67 <- FeaturePlot_scCustom(seurat_object = lm2.mg, features = "MKI67", split.by = "cell.type",num_columns = 4,colors_use = colo)+coord_fixed(ratio = 1)
ggsave(paste0(dir.figs, 'featureplot.MKI67.LM2', '.pdf'), FP.MKI67, height = 8, width = 12)
#=================================================================================================================#
#MVT1
FP.MKI67 <- FeaturePlot_scCustom(seurat_object = mv.mg, features = "Mki67", split.by = "cell.type",num_columns = 4,colors_use = colo)+coord_fixed(ratio = 1)
ggsave(paste0(dir.figs, 'featureplot.MKI67.MVT1', '.pdf'), FP.MKI67, height = 8, width = 12)
#=================================================================================================================#
# Br16
FP.MKI67 <- FeaturePlot_scCustom(seurat_object = Br.mg, features = "MKI67", split.by = "orig.ident",num_columns = 4,colors_use = colo)+coord_fixed(ratio = 1)
ggsave(paste0(dir.figs, 'featureplot.MKI67.Br16', '.pdf'), FP.MKI67, height = 8, width = 12)

#=================================================================================================================#
# sessioninfo
#=================================================================================================================#
#refer to files/sessioninfo_detailed.txt

