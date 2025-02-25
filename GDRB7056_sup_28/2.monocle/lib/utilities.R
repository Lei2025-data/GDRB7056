readRDX <- function(file) {
	con <- gzfile(file)
	on.exit(close(con))
	magic <- readChar(con, 5L, useBytes = TRUE)
	if ( grepl("RD[ABX][2-9]\n", magic) ){
		object <- get(load(file))
	} else {
		object <- readRDS(file)
	}
	return(object)
}

Load <- function(file) {
	object <- readRDX(file)
	if ( "version" %in% slotNames(object) ) {
		if ( grepl('^2', object@version) ) {
			object <- Seurat::UpdateSeuratObject(object)
		}
	}
	return(object)
}


WriteTable <- function(x, file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, ...){
		write.table(x, file = file, quote = quote, sep = sep, row.names = row.names, col.names = col.names, ...  )
}


.GetMetaData <- function(object, cols = NULL) {
		name <- names(cols)		
		name[is.na(name) | name == ""] <- cols[is.na(name) | name == ""]
		if ( is.null(name) ) name <- cols
		names(cols) <- name

		cols <- cols[cols %in% colnames(obj@meta.data)]
		metadata <- object@meta.data[, cols]
		colnames(metadata) <- names(cols)
		metadata <- cbind(Cells = rownames(metadata), metadata)
		return(metadata)
}

CalAvgExp <- function(object, features = NULL, group.by = NULL, assay = NULL, slot = "data",
                      is.expm1 = ifelse(slot=="data", TRUE, FALSE), is.return = FALSE,
                      is.reverse = FALSE, is.bulk = FALSE, outfile = "AllGene.avg_exp.xls" ){
		data <- GetAssayData(object, assay = assay, slot = slot)
		if ( ! is.null(features) ) data <- data[features, ]
		if ( is.expm1 ) data <- expm1(data)
		if ( ! is.null(group.by) ) Idents(object) <- group.by
		mean_exp <- do.call(cbind, by(colnames(object), Idents(object), function(x)
							Matrix::rowMeans(data[, if(is.reverse){levels(x)[-as.numeric(x)]}else{x}, drop = F])))
		if ( is.bulk ) {
				mean_exp <- cbind(bulk = Matrix::rowMeans(data), mean_exp)
		}
		if ( is.return ) {
				return(mean_exp)
		}else{
				colnames(mean_exp) <- paste( "Cluster", colnames(mean_exp) )
				Gene_ID   <- ChangeOUTName(rownames(mean_exp), object@misc$fdata)
				Gene_name <- FindFeaturesName(object, rownames(mean_exp), "name")
				mean_exp  <- cbind(Gene_ID = Gene_ID, Gene_name = Gene_name, mean_exp)
				WriteTable(mean_exp, file = outfile)
		}
}

CalPctExp <- function(object, features = NULL, group.by = NULL, assay = NULL, slot = "counts",
					  is.return = FALSE, is.reverse = FALSE, is.bulk = FALSE,
					  outfile = "AllGene.pct_exp.xls"){
		data <- GetAssayData(object, assay = assay, slot = slot)
		if ( ! is.null(features) ) data <- data[features, ]
		if ( ! is.null(group.by) ) Idents(object) <- group.by
		mean_exp <- do.call(cbind, by(colnames(object), Idents(object), function(x) {
							y <- if(is.reverse){ levels(x)[-as.numeric(x)] } else { x }
							Matrix::rowSums(data[, y, drop = F] > 0 ) / length(y)
							}))
		if ( is.bulk ) {
				mean_exp <- cbind(bulk = Matrix::rowSums(data > 0) / ncol(data), mean_exp)
		}
		if ( is.return ) {
				return(mean_exp)
		}else{
				colnames(mean_exp) <- paste( "Cluster", colnames(mean_exp) )
				Gene_ID   <- ChangeOUTName(rownames(mean_exp), object@misc$fdata)
				Gene_name <- FindFeaturesName(object, rownames(mean_exp), "name")
				mean_exp  <- cbind(Gene_ID = Gene_ID, Gene_name = Gene_name, mean_exp)
				WriteTable(mean_exp, file = outfile)
		}
}

FindFeaturesID <- function(object, features, unlist = TRUE){
		object@misc$fdata <- AddUnderscore(object@misc$fdata)
		if ( all(rownames(object) %in% object@misc$fdata$dash) ) {
				rownames(object@misc$fdata) <- object@misc$fdata$dash
		}
		features <- sapply(X = features,
						FUN = function(x) {
							if ( !exists("fdata", object@misc) ) return(NULL)
							g1 <- toupper(rownames(object@misc$fdata)) %in% toupper(x)
							if ( sum(g1) > 0 ) return(rownames(object@misc$fdata)[g1])
							g2 <- toupper(object@misc$fdata$name) %in% toupper(x)
							if ( sum(g2) > 0 ) return(rownames(object@misc$fdata)[g2])
							g3 <- toupper(object@misc$fdata$merge_name) %in% toupper(x)
							if ( sum(g3) > 0 ) return(rownames(object@misc$fdata)[g3])
							g4 <- toupper(object@misc$fdata$underscore) %in% toupper(x)
							if ( sum(g4) > 0 ) return(rownames(object@misc$fdata)[g4])
							message("[WARNING] '", x, "' not found gene id.")
							return(NULL)
						})
		if ( unlist ) features <- unlist(features)
		return(features)
}

FindFeaturesName <- function(object, features, col = "merge_name", is.fast = FALSE) {
		if ( ! exists("fdata", object@misc) )
				return(features)
		object@misc$fdata <- AddUnderscore(object@misc$fdata)
		new <- gsub("_", "-", features)
		if ( all(new %in% object@misc$fdata$dash) ) {
				rownames(object@misc$fdata) <- object@misc$fdata$dash
				features <- new
		}
		if ( is.fast ) {
				Name <- object@misc$fdata[features, col]
				names(Name) <- features
		}else{
				Name <- sapply(features, function(x) {
								id <- object@misc$fdata[x, col]
								ifelse(is.null(id)||is.na(id), x, id)
						})
		}
		return(Name)
}

ChangeOUTName <- function(features, fdata) {
		features <- as.character(features)
		fdata <- AddUnderscore(fdata)
		if ( ! is.null(fdata) && all(features %in% fdata$dash) ) {
				underscore_id <- fdata$underscore
				names(underscore_id) <- fdata$dash
				features <- underscore_id[features]
		}
		return(features)
}

AddUnderscore <- function(data){
		if ( ! is.null(data) ) {
				if ( ! exists("underscore", data) || ! exists("dash", data) ) {
						data$underscore <- rownames(data)
						data$dash <- gsub("_", "-", rownames(data))
				}
		}
		return(data)
}

IfNull <- function(var, default = NULL){
		if ( is.null(var) ) {
				return(default)
		} else {
				return(var)
		}
}

