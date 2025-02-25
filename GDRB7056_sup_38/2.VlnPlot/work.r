#! /Bio/bin/Rscript-3.5.1_conda
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(RColorBrewer)
source("/public2/Bio/pipeline/SingleCell_Collections/SCellWare/v3.0/R/Seurat_lib.R", chdir = T)
source("box_violin.r")

load("/public2/Bio/Project/GDRB7056/GDRB7056_sup_28/1.pipe/pipe/2.Seurat/obj.Rda")
obj@meta.data$Groups[ obj@meta.data$orig.ident %in% c("WT-1","WT-2") ] = "WT"
obj@meta.data$Groups[ obj@meta.data$orig.ident %in% c("Aa-1","Aa-2") ] = "Aa"
obj@meta.data$Groups[ obj@meta.data$orig.ident %in% c("KO-1","KO-2") ] = "KO"
obj@meta.data$Groups = factor(obj@meta.data$Groups, levels = c("WT", "Aa", "KO"))
table(obj@meta.data$Groups)
obj@misc$color.group <- SetColor(obj@meta.data$Groups, "tsne", "set3")
names(obj@misc$color.group) = levels(obj@meta.data$Groups)
marker = readLines("marker.txt")

PlotVlnPlot(obj, marker, group.by = "Groups", outpref = "Point/All/ViolinPlot")
PlotVlnPlot2(obj, marker, group.by = "Groups", outpref = "Box/All/ViolinPlot")

obj_bak = obj
for (i in levels(obj_bak@meta.data$seurat_clusters)){
	dir.create(paste0("Point/Cluster", i), showWarnings = F, recursive = T)
	dir.create(paste0("Box/Cluster", i), showWarnings = F, recursive = T)
	obj = subset(obj_bak, idents = i)
	PlotVlnPlot(obj, marker, group.by = "Groups", outpref = paste0("Point/Cluster", i, "/ViolinPlot"))
	PlotVlnPlot2(obj, marker, group.by = "Groups", outpref = paste0("Box/Cluster", i, "/ViolinPlot"))
}
