DoDimReduc <- function(object, assay = NULL, pc.num = 50, is.checkpca = TRUE, n.components = 2L,
				reduction = NULL, reduction.surfix = NULL, dims = NULL, nn.name = NULL,
				is.tsne = TRUE, is.UMAP = TRUE, ...) {
		if ( is.null(assay) ) assay <- DefaultAssay(object)
		if ( is.null(reduction) ) reduction <- "pca"
		if ( is.null(reduction.surfix) ) reduction.surfix <- assay				
		if ( reduction == "pca" && is.null(nn.name) ) {
				message( "-->PCA<--" )
				pc.num <- min(pc.num, ncol(object) - 1)
				object <- RunPCA(object = object, assay = assay, npcs = pc.num, features = VariableFeatures(object), verbose = FALSE )
				object[[paste0("pca_", reduction.surfix)]] <- object[["pca"]]
				if ( is.checkpca ) CheckPCA(object)
				reduction <- "pca"
		}

		pc.num <- ncol(object[[reduction]])
		sig.PCs <- if ( is.null(dims) ) seq(pc.num) else dims
		if ( is.tsne && is.null(nn.name) ) {
				object <- DoRunTSNE(object, dims = sig.PCs, reduction = reduction, nn.name = nn.name, reduction.surfix = reduction.surfix, n.components = n.components, ...)
		}

		if ( is.UMAP ) {
				object <- DoRunUMAP(object, dims = sig.PCs, reduction = reduction, nn.name = nn.name, reduction.surfix = reduction.surfix, n.components = n.components)
		}

		return(object)
}

CheckPCA <- function(object, reduction = "pca"){
		p1 <- DimPlot(object, reduction = reduction, group.by = "orig.ident")
		w <- 6
		if ( exists("Phase", object@meta.data) ){
				p2 <- DimPlot(object, reduction = reduction, group.by = "Phase")
#				p1 <- plot_grid(p1, p2, nrow = 1)
				p1 <- p1 + p2
				w <- 12
		}
		ggsave(p1, file = "pcaPlot.pdf", width = w, height = 6)

		dims <- min(20, ncol(Reductions(object, reduction)))
		p3 <- DimHeatmap(object, dims = seq(dims), cells = 500, balanced = TRUE, ncol = 4, fast = FALSE, reduction = reduction)
		ggsave(p3, file = "pcaHeatmap.pdf", width = min(4, dims) * 4, height = ceiling(dims/4) * 4)
		if( !is.null(dev.list()) ) dev.off()

		p4 <- ElbowPlot(object, ndims = ncol(Reductions(object, reduction)), reduction = reduction)
		ggsave(p4, file = "pcaElbowPlot.pdf", width = 6, height = 6)
}


DoRunTSNE <- function(object, dims = 1:50, reduction = "pca", nn.name = NULL, reduction.surfix = NULL, n.components = 2L, ...) {
		message( "-->Run tSNE<--" )
		perplexity <- min(30, floor((ncol(object) - 1)/3))
		reduction.name <- if ( n.components == 2 ) 'tsne' else paste0('tsne', n.components)
#		if ( ! is.null(nn.name) ) {
#				object[["WNN"]] <- CreateDimReducObject(embeddings = object@neighbors[[nn.name]]@nn.dist, key = "WNN")
#				rownames(object[["WNN"]]@cell.embeddings) <- object@neighbors[[nn.name]]@cell.names
#				reduction <- "WNN"
#		}
		dims <- intersect(dims, seq(object[[reduction]]))
#		if ( is.null(nn.name) ) {
				object <- RunTSNE(object, dims = dims, dim.embed = n.components, perplexity = perplexity, 
						reduction.name = reduction.name, reduction = reduction, ...)
#		} else {
#				object <- RunTSNE(object, dim.embed = n.components, perplexity = perplexity,
#						reduction.name = reduction.name, distance.matrix = object@neighbors[[nn.name]]@nn.dist, ...)
#		}
		if ( ! is.null(reduction.surfix) ) {
				object[[paste0(reduction.name, "_", reduction.surfix)]] <- object[[reduction.name]]
		}
		return(object)
}


DoRunUMAP <- function(object, dims = 1:50, reduction = "pca", nn.name = NULL, reduction.surfix = NULL, n.components = 2L, umap.method = 'uwot', ...) {
		message( "-->Run UMAP<--" )
		n.neighbors <- min(30, length(dims))
		reduction.name <- if ( n.components == 2 ) 'umap' else paste0('umap', n.components)
		dims <- intersect(dims, seq(object[[reduction]]))
		if (  is.null(nn.name) ) {
				object <- RunUMAP(object, dims = dims, umap.method = umap.method, n.neighbors = n.neighbors,
						n.components = n.components, reduction.name = reduction.name, reduction = reduction, ...)
		} else {
				object <- RunUMAP(object, umap.method = umap.method, n.neighbors = n.neighbors,
						n.components = n.components, reduction.name = reduction.name, nn.name = nn.name, ...)
		}
		if ( ! is.null(reduction.surfix) ) {
				object[[paste0(reduction.name, "_", reduction.surfix)]] <- object[[reduction.name]]
		}
		return(object)
}
