
DoFindClusters <- function(object, reduction = "pca", dims = NULL, resolution = 0.5){
		if ( is.null(dims) || max(dims) > length(object[[reduction]]) ) dims <- seq(object[[reduction]])
		object <- FindNeighbors(object, reduction = reduction, dims = dims, force.recalc = TRUE)
		object <- FindClusters(object, resolution = resolution, temp.file.location = getwd())

		color.cluster <- if ( exists("fetch_color") ) fetch_color(nlevels(object@meta.data$seurat_clusters), "tsne", "set1")
				else rainbow(nlevels(object@meta.data$seurat_clusters))
		names(color.cluster) <- levels(object@meta.data$seurat_clusters)
		object@misc[["color.cluster"]] <- color.cluster

		return(object)
}


StatCluster <- function(object, group.by = "orig.ident", outpref = "Cluster.stat", stat.what = "seurat_clusters", ...){
		.StatCluster(object, stat.what = stat.what, outpref = outpref)
		.StatCluster_by(object, stat.what = stat.what, group.by = group.by, outpref = outpref)
		.PlotClusterStat(object, stat.what = stat.what, group.by = group.by, outpref = outpref, ...)
}

.StatCluster <- function(object, outpref = "Cluster.stat", stat.what = "seurat_clusters" ){
		metadata <- if ( class(object) == "Seurat" ) object@meta.data else object
		name <- switch(stat.what, "seurat_clusters" = "Cluster", "orig.ident" = "Samples", stat.what)
		Cluster.stat <- metadata %>% group_by(!! name := !! as.name(stat.what)) %>% 
				summarise('Cells number' = n(),
						  'Median Genes per Cell' = median(nFeature_RNA),
						  'Median UMI Counts per Cell' = median(nCount_RNA))
		WriteTable(Cluster.stat, paste0(outpref, ".xls"))
}

.StatCluster_by <- function(object, group.by = "orig.ident", outpref = "Cluster.stat", stat.what = "seurat_clusters"){
		metadata <- if ( class(object) == "Seurat" ) object@meta.data else object
		name.stat.what <- switch(stat.what, "seurat_clusters" = "Cluster", "orig.ident" = "Samples", stat.what)
		Cluster.stat <- metadata %>% group_by(name = !! as.name(group.by), !! name.stat.what := !! as.name(stat.what)) %>%
				summarise(y = n()) %>% mutate(!! name.stat.what := factor(!! as.name(name.stat.what), levels = c("Total", levels(!! as.name(name.stat.what))))) %>% 
				full_join(x = metadata %>% group_by(name = !! as.name(group.by)) %>%
						summarise(!! name.stat.what := factor("Total", levels = c("Total", levels(!!as.name(stat.what)))), y = sum(n())) ) %>%
				reshape2::dcast(as.formula(paste0(name.stat.what," ~ name")), fill = 0) %>%
				mutate_if(is.numeric, list(~paste0(., " (", round(./.[1] * 100,2), "%)")))
#				mutate_if(is.numeric, funs(paste0(., " (", round(./.[1] * 100,2), "%)")))
		name <- switch(group.by, "orig.ident" = "Samples", "seurat_clusters" = "Cluster", group.by)
		WriteTable(Cluster.stat, paste0(outpref, ".", name, ".xls"))
}

ListCellCluster <- function(object, outfile = "Cells.cluster.list.xls"){
		data <- .GetMetaData(object, cols = c("Cluster" = "seurat_clusters", "Samples" = "orig.ident", "Groups"))
		WriteTable(data, file = outfile)
}


