
### Deal arguements
args <- commandArgs(T)
#bin  <- dirname(normalizePath(sub('--file=', '',  args[grep('--file=', args)])))
#args <- args[-seq(grep("--args", args))]

file <- args[1]
outdir <- args[2]
if ( is.null( file ) | is.na( file ) ){
		warning( "\n  Usage : Seurat.R <parameter.yaml> (<outdir>)\n" )
		quit()
}


### Loading Library
handlers <- list("bool#no" = function(x){if ( x %in% c("false", "FALSE") ) FALSE else x}, "bool#yes" = function(x){if ( x %in% c("true", "TRUE") ) TRUE else x})
parameter <- yaml::yaml.load_file( file, handlers = handlers)

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
#library(cowplot)

library(future)
options(future.globals.maxSize = 100 * 1024 * 1024^2)
plan("multiprocess", workers = 4)
#plan("sequential")

source("Seurat_lib.R")
#source("/Bio/Bin/pipeline/SingleCell_Collections/SCellWare/v1.0/Seurat_lib.R", chdir = T)

### Let's shake it
#setwd( parameter$outdir )
if ( ! is.na(outdir) ) setwd(outdir)

### Creat Seurat Object
message( "==>Reading 10x data<==" )
#obj <- MakeSeuratObj(parameter)
load("/public2/Bio/Project/GDRB7056/GDRB7056_sup_21/3.diff/obj.Rda")


### Add some check flag
message( "==>Adding MetaData<==" )
#if ( ! is.null(parameter$marker$mito_list) ) obj <- StatFeatures(obj, parameter$marker$mito_list, col.name = "percent.mito", stat_pct = T, add_to_pdata = T)
if ( ! is.null(parameter$marker$expected) )  obj <- StatFeatures(obj, parameter$marker$expected,  col.name = "expected.marker")
if ( ! is.null(parameter$marker$excluded) )  obj <- StatFeatures(obj, parameter$marker$excluded,  col.name = "exclude.marker")
if ( ! is.null(parameter$marker$more) )      obj@misc[["more.marker"]] <- FindFeaturesID(obj, parameter$marker$more, unlist = FALSE)
WriteTable(tibble::rownames_to_column(obj@meta.data, var = "Cells"), file = "metadata.xls")

if (F){
### Data Stat - before filter
message( "==>Stat before BasicInfo<==" )
PlotBasicStat(obj, "BasicInfo")
if ( ! is.null(parameter$Groups) ) {
		PlotBasicStat(obj, "BasicInfo.groups", group.by = "Groups")
}

### Filter
message( "==>Filter<==" )
obj <- FilterGenes(obj, parameter)
obj <- FilterCells(obj, parameter, do.stat = FALSE)
StatFilterCells(obj, group.by = "orig.ident", outfile = "Filter.stat.xls")
if ( ! is.null(parameter$Groups) ) {
		StatFilterCells(obj, group.by = "Groups", outfile = "Filter.stat.groups.xls")
}

### Data Stat - after filter
message( "==>Stat after BasicInfo<==" )
PlotBasicStat(obj, "AfterFilter.BasicInfo")
if ( ! is.null(parameter$Groups) ) {
		PlotBasicStat(obj, "AfterFilter.BasicInfo.groups", group.by = "Groups")
}


### Normalization Data
message( "==>Normalization Data<==" )
obj <- DoNormalization(obj, parameter, is_SCTransform = FALSE)

### Reduce dimension
message( "==>Reduce dimension<==" )
obj <- DoDimReduc(obj)

### Find clusters
message( "==>Find clusters<==" )
obj <- DoFindClusters(obj, reduction = "pca", dims = NULL, resolution = parameter$cluster_resolution)

if ( parameter$is.integration ) {
		if ( length(table(obj[["orig.ident"]])) > 1 ) {
				## check before integration visualization
				obj[["beforeInteg.cluster"]] <- Idents(object = obj)
				PlotCluster(obj, reduction = 'umap_RNA', outpref = "UMAP_before" )
				PlotCluster(obj, reduction = 'tsne_RNA', outpref = "tSNE_before" )
				if ( ! is.null(parameter$Groups) ) {
						PlotCluster(obj, reduction = 'umap_RNA', outpref = "UMAP_before.groups", split.by = "Groups", p1.group.by = "Groups" )
						PlotCluster(obj, reduction = 'tsne_RNA', outpref = "tSNE_before.groups", split.by = "Groups", p1.group.by = "Groups" )
				}

				### Integration
				message( "==> Do Integration <==" )
				obj <- DoIntegration(obj, split.by = "orig.ident", dims = 1:50)

				### Reduce dimension 
				obj <- DoDimReduc(obj)

				### Find clusters
				obj <- DoFindClusters(obj, reduction = "pca", dims = NULL, resolution = parameter$cluster_resolution)
		}
}
}

## Draw t-SNE plot
message( "==>Draw t-SNE plot<==" )
PlotCluster(obj, reduction = 'umap', outpref = "UMAP" )
PlotCluster(obj, reduction = 'tsne', outpref = "tSNE" )
if ( ! is.null(parameter$Groups) ) {
		PlotCluster(obj, reduction = 'umap', outpref = "UMAP.groups", split.by = "Groups", p1.group.by = "Groups" )
		PlotCluster(obj, reduction = 'tsne', outpref = "tSNE.groups", split.by = "Groups", p1.group.by = "Groups" )
}

### Save data object 
message( "==>Output obj.Rda<==" )
DefaultAssay(obj) <- "RNA"
#save(obj, file = "obj.Rda")

### stat table 
message( "==>Stat table<==" )
StatCluster(obj)
if ( ! is.null(parameter$Groups) ) {
		StatCluster(obj, "Groups")
}
CalAvgExp(obj)
ListCellCluster(obj)
PlotPresetMarker(obj)


### Find maker genes
message( "==>Find maker genes<==" )
obj.markers <- DoFindAllMarkers(obj, parameter)

message( "==>Output markers.Rda<==" )
save( obj.markers, file = "markers.Rda" )
#obj.markers$gene <- ChangeOUTName(obj.markers$gene, object@misc$fdata)

## stat marker
message( "==>stat marker<==" )
StatMarker(obj.markers, color = obj@misc$color.cluster)
ListMarker(obj, obj.markers)

### Top marker
message( "==>display top markers<==" )
top <- FindTopMarker(obj.markers, top_num = parameter$heatmap$top, object = obj)
PlotAboutFeatures(obj, features = unique(top$gene), outpref = "Top")

dir.create("DensityPlot/", showWarnings = F, recursive = T)
unlink("DensityPlot/*", recursive = T)
PlotDensityPlot(obj, unique(top$gene), reduction = 'tsne', outpref = "DensityPlot/")

dir.create("ExpPlot/", showWarnings = F, recursive = T)
unlink("ExpPlot/*", recursive = T)
PlotFeaturePlot(obj, unique(top$gene), reduction = 'tsne', outpref = "ExpPlot/ExpPlot", is.combine = FALSE)

dir.create("ViolinPlot/", showWarnings = F, recursive = T)
unlink("ViolinPlot/*", recursive = T)
PlotVlnPlot(obj, unique(top$gene), outpref = "ViolinPlot/ViolinPlot")

### Hasta la vista, baby
message( "==>All Done!<==" )


