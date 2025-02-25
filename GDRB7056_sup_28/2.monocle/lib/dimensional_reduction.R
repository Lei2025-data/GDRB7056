DoDimReduc <- function(object, assay = NULL, pc.num = 50, is.checkpca = TRUE, n.components = 2L,
				is.tsne = TRUE, is.UMAP = TRUE, ...) {
		if ( is.null(assay) ) assay <- DefaultAssay(object)				
		message( "-->PCA<--" )
		pc.num <- min(pc.num, ncol(object) - 1)
		object <- RunPCA(object = object, assay = assay, npcs = pc.num, features = VariableFeatures(object), verbose = FALSE )
		object[[paste0("pca_", assay)]] <- object[["pca"]]
		pc.num <- ncol(object[["pca"]])
		if ( is.checkpca ) CheckPCA(object)

		sig.PCs <- seq(pc.num)
		if ( is.tsne ) {
				message( "-->Run tSNE<--" )
				perplexity <- min(30, floor(pc.num - 1)/3)
				reduction.name <- if ( n.components == 2 ) 'tsne' else paste0('tsne', n.components)
				object <- RunTSNE(object, dims = sig.PCs, dim.embed = n.components, perplexity = perplexity, 
								reduction.name = reduction.name, ...)
				object[[paste0(reduction.name, "_", assay)]] <- object[[reduction.name]]
		}

		if ( is.UMAP ) {
				message( "-->Run UMAP<--" )
#				object <- RunUMAP(object, dims = sig.PCs, umap.method = "umap-learn")
				n.neighbors <- min(30, pc.num)
				reduction.name <- if ( n.components == 2 ) 'umap' else paste0('umap', n.components)
				object <- RunUMAP(object, dims = sig.PCs, umap.method = "uwot", n.neighbors = n.neighbors,
								n.components = n.components, reduction.name = reduction.name)
				object[[paste0(reduction.name, "_", assay)]] <- object[[reduction.name]]
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
