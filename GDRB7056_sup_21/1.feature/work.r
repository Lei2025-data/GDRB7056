#! /state/partition1/WORK/lib-devel/working_library/anaconda2/envs/R-3.5.1/bin/Rscript
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
source("Seurat_lib.R")
load("/public2/Bio/Project/GDRB7056/GDRB7056_sup_14/point1/Seurat/obj.Rda")
dat <- read.table ("marker.list", sep = "\t", header = F)
marker <- as.vector(dat$V1)
bak <- obj
for ( i in seq(marker)){
	for ( j in c("WT", "Aa", "KO")){
		sam1 <- paste0(j, "-1")
		sam2 <- paste0(j, "-2")
		cell <- Cells(bak)[bak@meta.data$orig.ident %in% c(sam1, sam2)]
		obj <- subset(bak, cells = cell)
		PlotFeaturePlot(obj, features = marker[i], outfile = paste0(dat$V2[i],"_",j,".Distribution.pdf"), reduction = "tsne")
	}
}
