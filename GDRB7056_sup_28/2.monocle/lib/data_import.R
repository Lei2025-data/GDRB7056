
MakeSeuratObj <- function(parameter = list(), assay = "RNA",
				data_name = parameter$data$name, data_dir = parameter$data$dir,
				name_list = parameter$name_list, group.use = parameter$Groups)
{
		object.list <- list()
		for ( i in seq(data_name) ){
				mat <- Read10X(data.dir = data_dir[i], gene.column = 1)
				object.list[[i]] <- CreateSeuratObject(counts = mat, project = data_name[i], assay = assay )
		}
		if ( length(object.list) == 1 ) {
				object <- object.list[[1]]
		}else{
				object <- merge(x = object.list[[1]], y = unlist(object.list[-1]), add.cell.ids = data_name)
		}

		object@meta.data$orig.ident <- factor(object@meta.data$orig.ident, levels = data_name)
		object <- SetIdent(object, value = "orig.ident")

		object@misc[["fdata"]] <- AddFData(object, name_list)
		object@misc[["pdata"]] <- FetchData(object, c("orig.ident", paste0("nFeature_", assay), paste0("nCount_", assay)))
		object@misc[["counts"]] <- GetAssayData(object, slot = "counts", assay = assay)

		color.sample <- fetch_color(nlevels(object@meta.data$orig.ident), "tsne", "set3")
		names(color.sample) <- levels(object@meta.data$orig.ident)
		object@misc[["color.sample"]] <- color.sample

		if ( ! is.null(group.use) ){
				object@meta.data$Groups <- object@meta.data$orig.ident
				group <- unlist(lapply(names(group.use), function(i) {
								x <- rep(i, length.out = length(group.use[[i]]))
								names(x) <- group.use[[i]]
								x
						}))
				levels(object@meta.data$Groups) <- group[levels(object@meta.data$Groups)]
				object@meta.data$Groups <- factor(object@meta.data$Groups, levels = names(group.use))
				object@misc[["pdata"]]$Groups <- object@meta.data$Groups

				color.group <- fetch_color(nlevels(object@meta.data$orig.ident), "tsne", "set2")
				names(color.group) <- levels(object@meta.data$Groups)
				object@misc[["color.group"]] <- color.group
		}

		return(object)
}


AddFData <- function(object, ref_name_file = NULL, col.name = NULL){
		if ( ! is.null(ref_name_file) && file.exists(ref_name_file) ) {
				fdata <- read.table(ref_name_file, row.names = 1, stringsAsFactors = F, sep = "\t")
				colnames(fdata) <- c("merge_name", "name", "type") ## old_merge_name
		}else{
				fdata <- data.frame(name = rownames(object), row.names = rownames(object), stringsAsFactors = F) # old_merge_name = rownames(object)
		}
		fdata$merge_name <- fdata$name
		fdata$merge_name[fdata$merge_name == "-"] <- rownames(fdata)[fdata$merge_name == "-"]
		index <- c(which(duplicated(fdata$merge_name, fromLast=T)), which(duplicated(fdata$merge_name, fromLast=F)))
		fdata$merge_name[index] <- paste0(fdata$merge_name[index], " (", rownames(fdata)[index], ")")
		fdata <- AddUnderscore(fdata)
		if ( is.null(col.name) ) {
				return(fdata)
		}else{
				object@misc[[col.name]] <- fdata
				return(object)
		}
}


StatFeatures <- function(object, features = NULL, col.name = NULL, stat_pct = FALSE, assay = NULL, add_to_pdata = FALSE){
		if ( file.exists(features[1]) )
				features <- readLines(con = features[1])
		features <- FindFeaturesID(object = object, features = features, unlist = FALSE)
		if ( is.null(assay) )
				assay <- DefaultAssay(object = object)
		metadata <- Matrix::colSums(x = GetAssayData(object = object, slot = "counts", assay = assay)[unlist(features), , drop = FALSE])
		if (stat_pct)
				metadata <- metadata / object@meta.data[[paste0("nCount_", assay)]] * 100

		if (!is.null(x = col.name)) {
				object@misc[[col.name]] <- features
				if (add_to_pdata)
						object@misc$pdata[[col.name]] <- metadata
				object <- AddMetaData(object = object, metadata = metadata, col.name = col.name)
				return(object)
		} else {
				return(metadata)
		}
}


RestoreObject <- function(object) {
		metadata <- if( exists("pdata", object@misc) ){ object@misc[["pdata"]] } else { object@meta.data }
		new.object <- CreateSeuratObject(counts = object@misc[["counts"]], meta.data = metadata, project = object@project.name, assay = "RNA" )
		new.object@misc <- object@misc
		return(new.object)
}


