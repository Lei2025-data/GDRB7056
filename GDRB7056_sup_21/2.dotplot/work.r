#! /state/partition1/WORK/lib-devel/working_library/anaconda2/envs/R-3.5.1/bin/Rscript
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
source("/public2/Bio/Project/GDRB7056/GDRB7056_sup_14/point1/Seurat_lib.R")
load("/public2/Bio/Project/GDRB7056/GDRB7056_sup_14/point1/Seurat/obj.Rda")
dat <- read.table ("marker.list", sep = "\t", header = F)
marker <- as.vector(dat$V2)

cell <- Cells(obj)[obj@meta.data$seurat_clusters %in% c(0,1,2,3,4,5,6,7,8,12,13,14,17,18,23)]
obj <- subset(obj, cells = cell)

for ( i in 0:8){
	obj@meta.data$seurat_clusters_new[obj@meta.data$seurat_clusters == i] <- i
}
obj@meta.data$seurat_clusters_new[obj@meta.data$seurat_clusters == 12] <- 9
obj@meta.data$seurat_clusters_new[obj@meta.data$seurat_clusters == 13] <- 10
obj@meta.data$seurat_clusters_new[obj@meta.data$seurat_clusters == 14] <- 11
obj@meta.data$seurat_clusters_new[obj@meta.data$seurat_clusters == 17] <- 12
obj@meta.data$seurat_clusters_new[obj@meta.data$seurat_clusters == 18] <- 13
obj@meta.data$seurat_clusters_new[obj@meta.data$seurat_clusters == 23] <- 14

obj@meta.data$seurat_clusters <- factor(obj@meta.data$seurat_clusters_new, levels = c(8,14,9,3,11,10,5,7,1,0,13,2,12,4,6))
obj@misc$color.cluster <- rainbow(length(levels(obj@meta.data$seurat_clusters)))
names(obj@misc$color.cluster) <- levels(obj@meta.data$seurat_clusters)
Idents(obj) <- "seurat_clusters"

PlotDotPlot(obj, features = rev(marker), outfile = "Epithelial.DotPlot.pdf")
