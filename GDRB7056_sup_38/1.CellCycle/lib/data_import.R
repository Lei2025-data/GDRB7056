.Read10X <- function(data.path, use.names = FALSE, assay = NULL){
		data <- if ( dir.exists(data.path) ) {
				Seurat::Read10X(data.path, gene.column = ifelse(use.names, 2, 1))
			} else if ( grepl('\\.h5', data.path) ) {
					Seurat::Read10X_h5(data.path, use.names = use.names)
			} else {
					sep <- if ( grepl('\\.csv', data.path) ) ',' else '\t'
					read.table(data.path, header = T, row.names = 1, sep = sep)
			}
		if ( class(data) == "list" ) {
				names(data) <- sapply(names(data), function(x) switch(x, "Gene Expression" = "RNA", "Peaks" = "ATAC", x))
				if ( ! is.null(assay) ) {
						data <- data[assay]
				}
#				for ( i in names(data) ) {
#						if ( all(grepl(pattern = "-[0-9]+$", x = colnames(data[[i]]))) ){
#								colnames(data[[i]]) <- as.vector(x = as.character(x = sapply(X = colnames(data[[i]]), FUN = Seurat:::ExtractField, field = 1, delim = "-")))
#						}
#				}
				if ( length(data) == 1 ) {
						data <- data[[1]]
				}
		} else {
				if ( all(grepl(pattern = "-[0-9]+$", x = colnames(data))) ){
						colnames(data) <- as.vector(x = as.character(x = sapply(X = colnames(data), FUN = Seurat:::ExtractField, field = 1, delim = "-")))
				}
		}

		return(data)
}

MergeObject <- function(object.list, data_name = do.call(c, lapply(object.list, function(x) levels(x@meta.data$orig.ident))) ) {
		if ( class(object.list) == "list" ) {
				if ( length(object.list) == 1 ) {
						object <- RenameCells(object.list[[1]], add.cell.id = data_name)
						return(object)
				} else {
						object <- merge(object.list[[1]], object.list[-1], add.cell.ids = data_name)
						return(object)
				}
		} else { 
				return(object.list)
		}
}


MakeSeuratObj <- function(parameter = list(), assay = "RNA",
				data_name = parameter$data$name, data_dir = parameter$data$dir,
				name_list = parameter$name_list, group.use = parameter$Groups,
				use.names = FALSE)
{
		object.list <- list()
		for ( i in seq(data_name) ){
				mat <- .Read10X(data_dir[i], use.names = use.names, assay = assay)

				object.list[[i]] <- CreateSeuratObject(counts = mat, project = data_name[i], assay = assay )
		}
		object <- MergeObject(object.list)

		object <- SetSeuratInfo(object, data_name = data_name, name_list = name_list, group.use = group.use)

		return(object)
}

MakeSeuratObj_ARC <- function(parameter = list(), assay = "RNA", use.names = FALSE,
				assay.ATAC = "ATAC", do.merge_peak = TRUE,
				data_name = parameter$data$name, data_dir = parameter$data$dir,
				data_frag = parameter$data$fragment, data_meta = parameter$data$metadata,
				name_list = parameter$name_list, group.use = parameter$Groups,
				gtf_file  = parameter$ref$gtf)
{
		object.list <- list()
		for ( i in seq(data_name) ){
				counts <- .Read10X(data_dir[i], use.names = use.names)
				object.list[[i]] <- CreateSeuratObject(counts = counts[[assay]], project = data_name[i], assay = assay )
				object.list[[i]][[assay.ATAC]] <- CreateChromatinAssay(counts = counts[[assay.ATAC]], fragments = data_frag[i], sep = c(":", "-(?=\\d+$)"), min.cells = -1, min.features = -1)

		
				metadata <- read.csv(data_meta[i], header = TRUE, row.names = 1, stringsAsFactors = FALSE)
				metadata <- metadata[, c("is_cell", "atac_peak_region_fragments", "atac_fragments")]

				metadata <- metadata[colnames(object.list[[i]]), ]
				object.list[[i]] <- AddMetaData(object.list[[i]], metadata = metadata)
		}

		if ( do.merge_peak ) {
				object.list <- AddMergePeaks(object.list, assay.name = assay.ATAC)
		}

		object <- MergeObject(object.list)

		object <- SetSeuratInfo(object, data_name = data_name, name_list = name_list, group.use = group.use)

		if ( ! is.null(gtf_file) ) {
				gtf <- rtracklayer::import(con = gtf_file)
				if ( ! "gene_biotype" %in% names(gtf@elementMetadata@listData) ) {
						col <- grep("gene_type|biotype", names(gtf@elementMetadata@listData), value = T)[1]
						gtf$gene_biotype <- gtf@elementMetadata@listData[[col]]
				}
				Annotation(object[[assay.ATAC]]) <- gtf
		}

		DefaultAssay(object) <- assay

		return(object)
}

MakeSpatialObj <- function(parameter = list(), assay = "Spatial", filter.matrix = TRUE, 
				data_name = parameter$data$name, data_dir = parameter$data$dir,
				name_list = parameter$name_list, group.use = parameter$Groups,
				refdir = parameter$refdir, use.names = FALSE, image.dirname = "spatial",
				mat.dirname = "filtered_feature_bc_matrix", mat.h5name = "filtered_feature_bc_matrix.h5")
{
		object.list <- list()
		for ( i in seq(data_name) ){
				input <- file.path(data_dir[i], mat.dirname)
				if ( ! dir.exists(input) ) input <- file.path(data_dir[i], mat.h5name)
				mat <- .Read10X(input, use.names = use.names, assay = assay)
				object.list[[i]] <- CreateSeuratObject(counts = mat, project = data_name[i], assay = assay )

				image <- Read10X_Image(image.dir = file.path(data_dir[i], image.dirname), filter.matrix = filter.matrix)
				if ( all(grepl(pattern = "-[0-9]+$", x = rownames(image@coordinates))) ){
						rownames(image@coordinates) <- as.vector(x = as.character(x = sapply(X = rownames(image@coordinates), FUN = Seurat:::ExtractField, field = 1, delim = "-")))
				}
				image <- image[Cells(x = object.list[[i]])]
				DefaultAssay(image) <- assay

				object.list[[i]][[data_name[i]]] <- image
		}
		object <- MergeObject(object.list)
		names(object@images) <- data_name

		if ( is.null(name_list) && ! is.null(refdir) ) name_list <- file.path(refdir, "/genes/name_list.xls")
		object <- SetSeuratInfo(object, data_name = data_name, name_list = name_list, group.use = group.use, assay = assay)

		return(object)
}


SetSeuratInfo <- function(object, data_name = NULL, name_list = NULL,
				sample_col = "orig.ident", assay = "RNA",
				group.use = NULL, group_col = "Groups" ) {
		object@meta.data[[sample_col]] <- if ( is.null(data_name) ) {
				as.factor(object@meta.data[[sample_col]])
			} else { 
				factor(object@meta.data[[sample_col]], levels = data_name)
			}
		object <- SetIdent(object, value = sample_col)

		object@misc[["fdata"]] <- AddFData(object, name_list)
		object@misc[["pdata"]] <- FetchData(object, c(sample_col, paste0("nFeature_", assay), paste0("nCount_", assay)))
		object@misc[["counts"]] <- GetAssayData(object, slot = "counts", assay = assay)

		object@misc[["color.sample"]] <- SetColor(object@meta.data[[sample_col]], "tsne", "set3")

		if ( ! is.null(group.use) ){
				object@meta.data[[group_col]] <- object@meta.data[[sample_col]]
				group <- unlist(lapply(names(group.use), function(i) {
								x <- rep(i, length.out = length(group.use[[i]]))
								names(x) <- group.use[[i]]
								x
						}))
				levels(object@meta.data[[group_col]]) <- group[levels(object@meta.data[[group_col]])]
				object@meta.data[[group_col]] <- factor(object@meta.data[[group_col]], levels = names(group.use))
				object@misc[["pdata"]][[group_col]] <- object@meta.data[[group_col]]

				object@misc[["color.group"]] <- SetColor(object@meta.data[[group_col]], "tsne", "set2")
		}

		return(object)
}


AddFData <- function(object, ref_name_file = NULL, col.name = NULL){
		if ( ! is.null(ref_name_file) && file.exists(ref_name_file) ) {
				fdata <- read.table(ref_name_file, row.names = 1, stringsAsFactors = F, sep = "\t", quote = '"')
				if ( ncol(fdata) > 3 ) fdata <- fdata[, 1:3]
				colnames(fdata) <- c("merge_name", "name", "type")[1:ncol(fdata)] ## old_merge_name
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
		if ( is.null(features) || length(features) == 0 ) {
				warning("[features] is empty. return 'object' without any change.")
				return(object)
		}
		if ( file.exists(features[1]) )
				features <- readLines(con = features[1])
		features <- FindFeaturesID(object = object, features = features, unlist = FALSE)
		if ( is.null(assay) )
				assay <- DefaultAssay(object = object)
		features <- intersect(features, rownames(object[[assay]]))
		metadata <- Matrix::colSums(x = GetAssayData(object = object, slot = "counts", assay = assay)[features, , drop = FALSE])
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
		if ( ! is.null(object@misc[["counts"]]) && ncol(object@misc[["counts"]])!= ncol(object) ) {
				metadata <- if( exists("pdata", object@misc) ){ object@misc[["pdata"]] } else { object@meta.data }
				new.object <- CreateSeuratObject(counts = object@misc[["counts"]], meta.data = metadata, project = object@project.name, assay = "RNA" )
				new.object@misc <- object@misc
				return(new.object)
		} else {
				return(object)
		}
}


