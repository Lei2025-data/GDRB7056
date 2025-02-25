args <- commandArgs()
bin  <- dirname(normalizePath(sub('--file=', '',  args[grep('--file=', args)])))
args <- args[-seq(grep("--args", args))]


file   <- args[1]
outdir  <- args[2]
add_lib <- args[3]



library(yaml)
parameter <- yaml.load_file( file, handlers = list("bool#no" = function(x){if ( x %in% c("false", "FALSE") ) FALSE else x}, "bool#yes" = function(x){if ( x %in% c("true", "TRUE") ) TRUE else x}))

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

library(future)
options(future.globals.maxSize = 100 * 1024 * 1024^2)
plan("multicore", workers = 4)
#plan("sequential")
#reticulate::use_python("/state/partition1/WORK/lib-devel/working_library/anaconda2/envs/R-3.5.1/bin/python3.7", required = T)

source(paste0(bin, "/scDiffer_lib.R"))
if (! is.null(add_lib))
	source(add_lib)


setwd(outdir)

obj <- Load(parameter$data$object)

if ( !is.null(parameter$cluster.use) ) {
		obj <- subset(obj, idents = parameter$cluster.use)
}

group.data <- FindGroupIndex(obj, parameter)

marker <- DoFindGroupMarkers(obj, group.data, parameter)

save(marker, file = "marker.Rda")
WriteDifferMarker(marker, "all", object = obj)
StatDeGene(marker, group.by = "contrast")
PlotDeGeneVolcano(marker, parameter)

diff.marker <- lapply(marker, function(x) filter(x, sig != "nosig"))
WriteDifferMarker(diff.marker, "diff", object = obj)

PlotContrastHeatmap(obj, diff.marker, group.data, top.num = parameter$plots$top)
PlotContrastFeature(obj, diff.marker, group.data, top.num = parameter$plots$top, reduction = parameter$plots$reduce)
PlotContrastFeature(obj, diff.marker, group.data, top.num = parameter$plots$top, reduction = parameter$plots$reduce, is.demo = T)
PlotContrastDotPlot(obj, diff.marker, group.data, top.num = parameter$plots$top)

