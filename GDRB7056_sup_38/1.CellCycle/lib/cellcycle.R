ClusterByMarker <- function (object, features, ident.name = "Cluster", assign = NULL, null = "unclutsered", set.ident = FALSE, thres = 0, ...) {
		name <- "Cell.Cycle"
		object.cc <- AddModuleScore(object = object, features = features, name = name,
						ctrl = min(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1))), ...)
		cc.columns <- grep(pattern = name, x = colnames(x = object.cc[[]]), value = TRUE)
		cc.scores <- object.cc[[cc.columns]]
		rm(object.cc)
		Seurat:::CheckGC()
		if ( is.null(assign) ) assign <- names(features)
		assignments <- .CC.assignments(cc.scores, thres = thres, assign = assign, null = null)
		cc.scores <- merge(x = cc.scores, y = data.frame(assignments), by = 0)
		colnames(x = cc.scores) <- c("rownames", assign, ident.name)
		rownames(x = cc.scores) <- cc.scores$rownames
		cc.scores <- cc.scores[, c(assign, ident.name)]
		object[[colnames(x = cc.scores)]] <- cc.scores
		object@meta.data[[ident.name]] <- factor(object@meta.data[[ident.name]], levels = c(assign, null))
		if (set.ident) {
				object[["old.ident"]] <- Idents(object = object)
				Idents(object = object) <- ident.name
		}
		return(object)
}

.CC.assignments <- function(cc.scores, thres = 0, assign = c("G1/S", "S", "G2/M", "M", "M/G1"), null = "non-cycling" ) {
		if ( ncol(cc.scores) != length(assign) ) 
				stop("the length of `assign` must equal ncol of `cc.scores`.")
		assignments <- apply(X = cc.scores, MARGIN = 1,
						FUN = function(scores) {
							if (all(scores < thres)) {
								return(null)
							} else {
								if (length(which(x = scores == max(scores))) > 1) {
#									return("Undecided")
									return(null)
								} else {
									return(assign[which(x = scores == max(scores))])
								}
							}
						})
		return(assignments)
}

.FindCCGene <- function(object, cc.genes) {
		for ( i in seq(cc.genes) ) {
			cc.genes[[i]] <- FindFeaturesID(object, cc.genes[[i]])
		}
		return(cc.genes)
}

DoClusterByMarker <- function(object, cc.genes, ...){
#		cc.genes <- lapply(cc.genes, FindFeaturesID, object = object)
		cc.genes <- .FindCCGene(object, cc.genes)
#		if ( any(sapply(list(cc.genes), length) < 2) ) {
#				warning("some cc.genes[[]] are least than 2. Return, don't run ClusterByMarker.")
#				return(object)
#		}
		object <- ClusterByMarker(object, cc.genes, set.ident = FALSE, ...)
		return(object)
}

PlotCellCycle <- function(object, group.by = "Phase", outpref = "CellCycle", color = NULL){
		.PlotCluster(object, reduction = 'tsne', p1.group.by = "orig.ident",      p2.group.by = group.by, outfile = paste0(outpref, ".Samples.tSNE.pdf"),  p2.color = color)
		.PlotCluster(object, reduction = 'tsne', p1.group.by = "seurat_clusters", p2.group.by = group.by, outfile = paste0(outpref, ".Cluster.tSNE.pdf"), p2.color = color)
		.PlotCluster(object, reduction = 'umap', p1.group.by = "orig.ident",      p2.group.by = group.by, outfile = paste0(outpref, ".Samples.UMAP.pdf"),  p2.color = color)
		.PlotCluster(object, reduction = 'umap', p1.group.by = "seurat_clusters", p2.group.by = group.by, outfile = paste0(outpref, ".Cluster.UMAP.pdf"), p2.color = color)
}

database2table <- function(cc, outfile, object = NULL){
		if ( ! is.null(object) ){
				for ( i in seq(cc) ) {
						cc[[i]] <- FindFeaturesName(object, cc[[i]], is.fast = TRUE)
				}
		}
		cc.len <- sapply(cc, length)
		dt <- do.call(cbind, cc)
		index <- unlist(sapply(seq(cc.len), function(i) 1:cc.len[i] + (i-1) * max(cc.len)))
		dt[-index] <- '-'

		write.table(dt, file = outfile, sep = "\t", quote = F, row.names = F)
}

