
args <- commandArgs(T)

file    <- args[1]
outdir  <- args[2]
add_lib <- args[3]

handlers <- list("bool#no" = function(x){if ( x %in% c("false", "FALSE") ) FALSE else x}, "bool#yes" = function(x){if ( x %in% c("true", "TRUE") ) TRUE else x})
parameter <- yaml::yaml.load_file( file, handlers = handlers )
parameter$database <- yaml::yaml.load_file(parameter$database, handlers = handlers)

dir.create(outdir, showWarnings = F, recursive = T)
setwd(outdir)


library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

if ( !is.na(add_lib) )
	source(add_lib, chdir = T)
#source(paste0(bin, "/CellCycle_lib.R"), chdir = T)



obj     <- Load(parameter$obj_file)
cc.gene <- parameter$database$marker
cc.gene <- .FindCCGene(obj, cc.gene)
database2table(cc.gene, "CellCycle.genes.xls", obj)

obj <- DoClusterByMarker(obj, cc.gene, thres = parameter$thres, ident.name = "Phase", assign = parameter$database$assign, null = parameter$database$unassign)

color.cc <- if ( is.null(parameter$database$color.cc) ){
				c(scales::hue_pal()(length(levels(obj@meta.data$Phase)) - 1), "grey") 
		}else{
				parameter$database$color.cc
		}
names(color.cc) <- levels(obj@meta.data$Phase)

sample_name <- IfNull(parameter$sample_name, 'orig.ident')
cluster_name <- IfNull(parameter$cluster_name, 'seurat_clusters')
obj@meta.data$orig.ident <- obj@meta.data[[sample_name]]
obj@meta.data$seurat_clusters <- obj@meta.data[[cluster_name]]

vars <- if ( is.null(parameter$database$assign) ) names(cc.gene) else parameter$database$assign
cc_annot <- obj@meta.data %>%
	tibble::rownames_to_column("Cells") %>%
	mutate(CellCycle.Score = apply(.[,vars], 1, max)) %>%
	select(Cells, Samples = orig.ident, Clusters = seurat_clusters, !!!enquos(vars), CellCycle.Score, Phase)
WriteTable(cc_annot, file = "CellCycle.annot.xls")
saveRDS(cc_annot, file = "CellCycle.annot.Rds")

p1 <- ggplot(cc_annot, aes(x = Samples, y = CellCycle.Score)) + geom_boxplot() #+ mytheme_box 
p1 <- p1 + box_theme_default()
ggsave(p1, file = "CellCycle.boxplot.pdf")

stat_cc_samples <- cc_annot %>% group_by(Samples, Phase) %>% summarise("Number of cells" = n())
WriteTable(stat_cc_samples, file = "Phase.stat.Samples.xls")
stat_cc_cluster <- cc_annot %>% group_by(Clusters, Phase) %>% summarise("Number of cells" = n())
WriteTable(stat_cc_cluster, file = "Phase.stat.Cluster.xls")

.PlotClusterStat(obj, stat.what = "Phase", group.by = "orig.ident", color.st = color.cc)
.PlotClusterStat(obj, stat.what = "Phase", group.by = "seurat_clusters", color.st = color.cc)



PlotCellCycle(obj, color = color.cc, outpref = "CellCycle")	

pdf("CellCycle.Samples.circos.pdf")
.plotCircos(cc_annot[,c("Phase", "Samples")], color.1 = color.cc, color.2 = obj@misc$color.sample, show.axes = T)
dev.off()

pdf("CellCycle.Cluster.circos.pdf")
.plotCircos(cc_annot[,c("Phase", "Clusters")], color.1 = color.cc, color.2 = obj@misc$color.cluster, show.axes = T)
dev.off()

p <- ggplot(cc_annot, aes(x = Samples, fill = Phase)) +
	geom_bar(stat = "count", position = position_fill()) + 
	facet_grid(~Clusters) +
	scale_fill_manual(values = color.cc) +
	ylab("Fraction") + 
	bar_theme_default() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(p, file = "CellCycle.Samples_Cluster.barplot.pct.pdf", height = 8, width = 8 * length(unique(cc_annot$Samples)), limitsize = FALSE)


#### Plot cc genes
#top <- unique(unlist(lapply(cc.gene, head, 5)))
PlotAboutFeatures(obj, features = unlist(cc.gene), group.by = "Phase", group.colors = color.cc, outpref = "CellCycle")


#Idents(obj) <- "Phase"
if ( parameter$is.save_object ) {
	save(obj, file = "obj.Rda")
}


if ( ! is.null(parameter$FindMarkers$is.differ) && parameter$FindMarkers$is.differ ) {
	dir.create("differ/", showWarnings = F, recursive = T)
	setwd("differ/")

	CalAvgExp(obj, group.by = "Phase")

	### Find maker genes
	message( "==>Find maker genes<==" )
	obj.markers <- DoFindAllMarkers(obj, parameter, group.by = "Phase")
	print(dim(obj.markers))
	obj.markers <- obj.markers[!obj.markers$gene %in% unlist(cc.gene), ]
	print(dim(obj.markers))

	message( "==>Output markers.Rda<==" )
	save( obj.markers, file = "markers.Rda" )
	#obj.markers$gene <- ChangeOUTName(obj.markers$gene, object@misc$fdata)

	## stat marker
	message( "==>stat marker<==" )
	StatMarker(obj.markers, color = color.cc, Cluster_name = "Phase")
	ListMarker(obj, obj.markers, group.by = "Phase", group.by.name = "Phase")

	### Top marker
	message( "==>display top markers<==" )
	top <- FindTopMarker(obj.markers, top_num = parameter$FindMarkers$top, object = obj)
	PlotAboutFeatures(obj, features = unique(top$gene), outpref = "Top", group.by = "Phase", group.colors = color.cc)


	reduction <- if ( is.null(parameter$reduction) ) 'tsne' else parameter$reduction
	redct.name <- switch(reduction, tsne = 'TSNE', umap = 'UMAP')

	dir.create("DensityPlot/", showWarnings = F, recursive = T)
	unlink("DensityPlot/*", recursive = T)
	PlotDensityPlot(obj, unique(top$gene), reduction = reduction, outpref = "DensityPlot/")

	dir.create("ExpPlot/", showWarnings = F, recursive = T)
	unlink("ExpPlot/*", recursive = T)
	PlotFeaturePlot(obj, unique(top$gene), reduction = reduction, outpref = "ExpPlot/ExpPlot", is.combine = FALSE)

	dir.create("ViolinPlot/", showWarnings = F, recursive = T)
	unlink("ViolinPlot/*", recursive = T)
	PlotVlnPlot(obj, unique(top$gene), outpref = "ViolinPlot/ViolinPlot", group.by = "Phase", cols.use = color.cc)

	### Hasta la vista, baby
	message( "==>All Done!<==" )

}


