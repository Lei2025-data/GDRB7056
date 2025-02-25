#! /Bio/bin/Rscript-3.5.1_conda
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(RColorBrewer)
source("/public2/Bio/pipeline/SingleCell_Collections/SCellWare/v3.0/R/Seurat_lib.R", chdir = T)
source("/public2/Bio/Project/GDRB7056/GDRB7056_sup_38/VlnPlot/box_violin.r")

load("/public2/Bio/Project/GDRB7056/GDRB7056_sup_28/1.pipe/pipe/2.Seurat/obj.Rda")
obj@meta.data$Groups[ obj@meta.data$orig.ident %in% c("WT-1","WT-2") ] = "WT"
obj@meta.data$Groups[ obj@meta.data$orig.ident %in% c("Aa-1","Aa-2") ] = "Aa"
obj@meta.data$Groups[ obj@meta.data$orig.ident %in% c("KO-1","KO-2") ] = "KO"
obj@meta.data$Groups = factor(obj@meta.data$Groups, levels = c("WT", "Aa", "KO"))
table(obj@meta.data$Groups)
obj@misc$color.group <- SetColor(obj@meta.data$Groups, "tsne", "set3")
names(obj@misc$color.group) = levels(obj@meta.data$Groups)
obj_bak = obj

for (dir in dir("/public2/Bio/Project/GDRB7056/GDRB7056_sup_41")){
	setwd("/public2/Bio/Project/GDRB7056/GDRB7056_sup_41")
	if (file.info(dir)$isdir){
		message("==> Processing ", dir, " <==")
		setwd(paste0("/public2/Bio/Project/GDRB7056/GDRB7056_sup_41/", dir))
		if (file.exists("marker.txt")){
			marker = read.table("marker.txt")
			marker = as.vector(marker$V2)
			dir.create("Point/All", recursive = T)
			dir.create("Box/All", recursive = T)
			for (id in marker){
				message("==> Gene: ", id, " <==")
				cell = Cells(obj_bak)[obj_bak@assays$RNA@counts[id,] > 0]
				obj_bak2 = subset(obj_bak, cells = cell)
				PlotVlnPlot(obj_bak2, id, group.by = "Groups", outpref = "Point/All/ViolinPlot")
				PlotVlnPlot2(obj_bak2, id, group.by = "Groups", outpref = "Box/All/ViolinPlot")
				
				for (i in levels(obj_bak@meta.data$seurat_clusters)){
					message("==> Cluster: ", i, " <==")
					dir.create(paste0("Point/Cluster", i), showWarnings = F, recursive = T)
					dir.create(paste0("Box/Cluster", i), showWarnings = F, recursive = T)
					if ( length(Cells(obj_bak2)[obj_bak2@meta.data$seurat_clusters == i]) == 0 ){
						obj = subset(obj_bak, idents = i)
					}else{
						obj = subset(obj_bak2, idents = i)
					}
					PlotVlnPlot(obj, id, group.by = "Groups", outpref = paste0("Point/Cluster", i, "/ViolinPlot"))
					PlotVlnPlot2(obj, id, group.by = "Groups", outpref = paste0("Box/Cluster", i, "/ViolinPlot"))
				}
			}
		}
	}
}
