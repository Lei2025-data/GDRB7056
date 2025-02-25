
DoIntegration <- function(object, split.by = "orig.ident", dims = 1:50, nfeatures = 3000, is.SCT = FALSE, normalization.method = "LogNormalize"){
		old.assay <- DefaultAssay(object)
		object.list <- SplitObject(object, split.by = split.by)
		object.list <- SplitObject.Image(object.list)
		anchor.features <- nfeatures
		if ( is.SCT ) {
				for (i in seq(object.list)) {
						object.list[[i]] <- SCTransform(object.list[[i]], vars.to.regress = object@misc$vars.regress, verbose = FALSE, assay = old.assay)
				}
				anchor.features <- SelectIntegrationFeatures(object.list = object.list, nfeatures = nfeatures)
				object.list <- PrepSCTIntegration(object.list = object.list, anchor.features = anchor.features)
				normalization.method <- "SCT"
		}
		k.filter <- min(200, ceiling(min(sapply(object.list, ncol))/2))
		if ( any(sapply(object.list, ncol) <= max(dims)) ) {
				message(paste(sapply(object.list, ncol), collapse = " "))
				dims <- seq(min(sapply(object.list, ncol)) - 1)
		}
		anchors <- FindIntegrationAnchors(object.list = object.list, dims = dims, normalization.method = normalization.method, anchor.features = anchor.features, k.filter = k.filter)
		if ( nrow(anchors@anchors) == 0 ) {
				message("[Integrate] anchors is 0. No Integrating.")
				return(object)
		}
		integrated <- IntegrateData(anchorset = anchors, dims = dims, normalization.method = normalization.method, k.weight = k.filter)
		if ( ! is.SCT ) {
				integrated <- ScaleData(integrated, verbose = FALSE, vars.to.regress = object@misc$vars.regress)
		}
		integrated@misc <- object@misc
		integrated[[old.assay]] <- object[[old.assay]]
		integrated@reductions <- object@reductions
		integrated@meta.data <- object@meta.data
		return(integrated)
}

SplitObject.Image <- function(objects, names.in = "orig.ident"){
		if ( class(objects) != "list") {
				if ( "images" %in% slotNames(objects) && length(objects@images) > 1 ) {
						keep.image <- levels(droplevels(objects@meta.data)[[names.in]])
						objects@images <- objects@images[keep.image]
				}
		} else {
				for ( i in seq(objects) ) {
						if ( "images" %in% slotNames(objects[[i]]) && length(objects[[i]]@images) > 1 ) {
								keep.image <- levels(droplevels(objects[[i]]@meta.data)[[names.in]])
								objects[[i]]@images <- objects[[i]]@images[keep.image]
						}
				}
		}
		return(objects)
}


