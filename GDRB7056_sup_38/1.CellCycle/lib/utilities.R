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

		cols <- cols[cols %in% colnames(object@meta.data)]
		metadata <- object@meta.data[, cols]
		colnames(metadata) <- names(cols)
		metadata <- cbind(Cells = rownames(metadata), metadata)
		return(metadata)
}

CalAvgExp <- function(object, features = NULL, group.by = NULL, assay = NULL, slot = "data",
                      is.expm1 = ifelse(slot=="data", TRUE, FALSE), is.return = FALSE,
                      is.reverse = FALSE, is.bulk = FALSE, outfile = "AllGene.avg_exp.xls" ){
		data <- GetAssayData(object, assay = assay, slot = slot)
		if ( ! is.null(features) ) data <- data[features, , drop = FALSE]
		if ( is.expm1 ) data <- expm1(data)
		if ( ! is.null(group.by) ) Idents(object) <- group.by
		mean_exp <- do.call(cbind, by(colnames(object), Idents(object), function(x) {
							y <- if(is.reverse){ setdiff(colnames(object), x) } else { x }
							Matrix::rowMeans(data[, y, drop = F])
							}, simplify = FALSE))
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
					  outfile = "AllGene.avg_pct.xls"){
		data <- GetAssayData(object, assay = assay, slot = slot)
		if ( ! is.null(features) ) data <- data[features, , drop = FALSE]
		if ( ! is.null(group.by) ) Idents(object) <- group.by
		mean_exp <- do.call(cbind, by(colnames(object), Idents(object), function(x) {
							y <- if(is.reverse){ setdiff(colnames(object), x) } else { x }
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


SetColor <- function(x, type = "tsne", tag = "set1", ...) {
		if ( ! is.factor(x) ) x <- as.factor(x)
		color <- if ( exists("fetch_color") ) fetch_color(n = nlevels(x), type = type, tag = tag, ...)
				else rainbow(n = nlevels(x))
		names(color) <- levels(x)
		return(color)
}


unlist.rev <- function(data) {
		data <- unlist(
				lapply(names(data), function(x)
						`names<-`(x = rep(x, length(data[[x]])), value = data[[x]])
				)
		)
		return(data)
}

SubsetObj <- function(object, cells = NULL, sample = NULL, cluster = NULL, sample.name = "orig.ident", cluster.name = "seurat_clusters", ...) {
		cells <- if ( ! is.null(cells) ) cells else Cells(object)
		if ( ! is.null(sample) ) {
				if ( is.list(sample) ) {
				## TODO : rename sample
				} else {
						cells.sample <- Cells(object)[object[[sample.name]][[1]] %in% sample]
				}
				cells <- intersect(cells, cells.sample)
		}
		if ( ! is.null(cluster) ) {
				if ( is.list(cluster) ) {
				## TODO : rename cluster
				} else {
						cells.cluster <- Cells(object)[object[[cluster.name]][[1]] %in% cluster]
				}
				cells <- intersect(cells, cells.cluster)
		}
		others <- list(...)
		if ( length(others) > 0 ) {

		}
		object <- object[, cells]
		object@meta.data <- droplevels(object@meta.data)
		if ( 'images' %in% slotNames(object) ) {
				## when object[, cells], if image's name is like 'A-B', it will be duplicated with name 'A.B',
				## here substracting right name images with 'sample.name' in object@meta.data
				object@images <- object@images[Images(object) %in% unique(object@meta.data[[sample.name]])]
				for ( image in Images(object) ) {
						## object[, cells] did not apply to object@images
						image.cells <- intersect(rownames(object@images[[image]]@coordinates), cells)
						object@images[[image]]@coordinates <- object@images[[image]]@coordinates[image.cells, ]
				}
		}
		return(object)
}

col2grey <- function(red, green, blue, algorithms = c("luminance", "luma", "average", "desaturation", "max", "min", "red", "blue", "green"), maxColorValue = 255){
		if (missing(green) && missing(blue)) {
				if (is.matrix(red) || is.data.frame(red)) {
						red <- data.matrix(red)
						if (ncol(red) < 3L) stop("at least 3 columns needed")
				} else {
						red <- t(col2rgb(red))
				}
				green <- red[, 2L]
				blue <- red[, 3L]
				red <- red[, 1L]
		}
		algorithms <- match.arg(algorithms)
		#https://tannerhelland.com/2011/10/01/grayscale-image-algorithm-vb6.html
		Y <- switch(algorithms,
				luminance = 0.299 * red + 0.587 * green + 0.114 * blue,
				luma = 0.2126 * red + 0.7152 * green + 0.0722 * blue,
				average = mean(c(red, green, blue)),
				desaturation = (max(red, green, blue) + min(red, green, blue))/2,
				max = max(red, green, blue),
				min = min(red, green, blue),
				red = red,
				green = green,
				blue = blue
				)
		color <- rgb(red = Y, green = Y, blue = Y, maxColorValue = maxColorValue)
		return(color)
}

FindWHNum <- function(data, ncol = NULL, nrow = NULL, by.col = TRUE ) {
		if ( length(data) == 1 ) {
				if ( ! is.numeric(data) ) stop()
				data <- round(data)
		} else {
				data <- length(data)
		}

		if ( ! is.null(ncol) && ! is.null(nrow) ) {
				if ( ncol * nrow < data ) stop()
		} else if ( ! is.null(ncol) ) {
				nrow <- ceiling(data / ncol)
		} else if ( ! is.null(nrow) ) {
				ncol <- ceiling(data / nrow)
		} else {
				if ( by.col ) {
						ncol <- ceiling(sqrt(data))
						nrow <- ceiling(data / ncol)
				} else {
						nrow <- ceiling(sqrt(data))
						ncol <- ceiling(data / nrow)
				}
		}
		return(c(nrow, ncol))
}

