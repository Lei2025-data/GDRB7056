
RenameFeatures <- function(object, new.names = NULL,
				form.type = c("id", "merge_name", "old_merge_name"),
				to.type = c("id", "merge_name", "old_merge_name"),
				features = NULL
				) {
		if ( is.null(new.names) ) {
				form.type <- match.arg(form.type)
				to.type <- match.arg(to.type)
				fdata <- object@misc$fdata
				fdata$id <- rownames(fdata)
				new.names <- fdata[,to.type]
				names(new.names) <- fdata[,form.type]
#				new.names <- data.frame(to.type = fdata[,to.type], row.names = fdata[,form.type])
		}
		if ( is.null(features) ){
				assays <- Seurat:::FilterObjects(object = object, classes.keep = "Assay")
				for (assay in assays) {
						slot(object = object, name = "assays")[[assay]] <- RenameFeatures.Assays(object = object[[assay]], new.names = new.names)
				}
				dimreducs <- Seurat:::FilterObjects(object = object, classes.keep = "DimReduc")
				for (dr in dimreducs) {
						object[[dr]] <- RenameFeatures.DimReduc(object = object[[dr]], new.names = new.names)
				}
				return(object)
		} else {
				return(as.vector(new.names[features]))
		}
}

RenameFeatures.Assays <- function(object, new.names = NULL ) {
		for (data.slot in c("counts", "data", "scale.data")) {
				old.data <- GetAssayData(object = object, slot = data.slot)
				if (nrow(x = old.data) <= 1) {
						next
				}
				old.name <- rownames(x = slot(object = object, name = data.slot))
				rownames(x = slot(object = object, name = data.slot)) <- new.names[old.name] ## as.vector(new.names[old.name]) ?? 
		}
		if ( length(slot(object = object, name = "var.features")) > 0 ) {
				old.name <- rownames(x = slot(object = object, name = "var.features"))
				slot(object = object, name = "var.features") <- new.names[old.name] ## as.vector(new.names[old.name]) ??
		}
#		old.name <- rownames(x = slot(object = object, name = "meta.features")) ## ??
#		rownames(x = slot(object = object, name = "meta.features")) <- as.vector(new.names[old.name])  ## ??
		return(object)
}

RenameFeatures.DimReduc <- function(object, new.names = NULL ) {
		for ( projected in c(TRUE, FALSE) ){
				data.slot <- ifelse(projected, "feature.loadings.projected", "feature.loadings")
				old.data <- Loadings(object = object, projected = projected)
				rownames(x = old.data) <- new.names[rownames(x = old.data)] ## as.vector(new.names[rownames(x = old.data)]) ??				
				slot(object = object, name = data.slot) <- old.data
		}
		return(object)
}



