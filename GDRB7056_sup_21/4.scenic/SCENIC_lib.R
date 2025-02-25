##### sub function  #####
library(SCENIC)
library(AUCell)

RunMessage <- function(strings)
{
	cat(paste0("[", as.character(Sys.time()) ,"] [Running]: ", strings),"\n")
}

WarnMessage <- function(strings)
{
	cat(paste0("[", as.character(Sys.time()) ,"] [Warning]: ", strings),"\n")
}

StopMessage <- function(strings)
{
	cat(paste0("[", as.character(Sys.time()) ,"] [Stopped]: ", strings),"\n")
	q(status=1)
}

Message <- function(strings)
{
	cat("\t\t\t\t   ", strings, "\n")
}



LibraryPackages <- function(packages, lib.loc)
{
	RunMessage('>>>> Start Load Rpackages <<<<')
	RunMessage('   ---> Check Rpackages <---  ')
	libs <- list.files(lib.loc)
	error.packages <- NULL
	error.packages <- packages[! packages %in% libs]
	if(length(error.packages)){
		StopMessage(paste0("   ---> Please install packages: [", paste(error.packages,collapse=", "), "]"))
	}
	RunMessage('   ---> Loading Rpackages <---  ')
	for (i in packages){
		Message(paste0("    ~~~> ", i, " <~~~"))
		library(i, character.only=T)
	}
}

CreateScenicOptions <- function (parameter){
	RunMessage(">>>> Initial scenicOptions <<<<")
	RunMessage("    ---> Import cisTarget database <---    ")
	org <- NULL
	cisTarget_database <- NULL
	if(parameter$species == "Homo sapiens"){
		org <- "hgnc"
		cisTarget_database <- "/public2/Bio/Project/yinquan/StoreHouse/cisTarget_database/Human"
	}else if(parameter$species == "Mus musculus"){
		org <- "mgi"
		cisTarget_database <- "/public2/Bio/Project/yinquan/StoreHouse/cisTarget_database/Mouse"
	}else if (parameter$species == "Drosophila melanogaster"){
		org <- "dmel"
		cisTarget_database <- "/public2/Bio/Project/yinquan/StoreHouse/cisTarget_database/Fly"
	}else{
		StopMessage("    ---> Species must be any one of [Homo sapiens, Mus musculus, Drosophila melanogaster ]")
	}
	data(defaultDbNames)
	dbs <- defaultDbNames[[org]]
	outdir <- dirname(parameter$outdir)
	scenicOptions <- initializeScenic(org=org, dbDir=cisTarget_database, dbs=dbs, datasetTitle=parameter$project, nCores=parameter$nCores, outdir = outdir)
	return(scenicOptions)
}

TransMetaData <- function(object, orig.labels = NULL, new.labels = NULL, orig.col.name = NULL, new.col.name = NULL, setIdent = F){
	if(length(orig.labels) != length(new.labels)){
		StopMessage("    ---> orig.labels must be equal to new.labels in length <---   ")
	}

	for( i in seq(length(new.labels))){
		object@meta.data[[new.col.name]][object@meta.data[[orig.col.name]] %in% orig.labels[[i]]] <- new.labels[i]
	}
	object@meta.data[[new.col.name]] <- factor(object@meta.data[[new.col.name]], levels = new.labels)
	if(setIdent){
		Idents(object) <- new.col.name
	}
	return(object)
}

## Add fData
AddFdata <- function(object, namelist)
{
        namelist <- read.table(file=namelist, sep="\t", header=F, row.names=1)
        colnames(namelist) <- c("merge_name","name","biotype")
        namelist$nounderscore <- gsub("_", "-", namelist$merge_name)
        namelist$dash <- gsub("_","-",rownames(namelist))
        namelist <- cbind(geneid=rownames(namelist), namelist)
        object@misc$fData <- namelist
        object
}

LoadObject <- function(parameter)
{
	RunMessage(">>>> Load Seurat object <<<<")
	
	tryCatch( { read.obj <- load(parameter$object) }, error = function(e) { read.obj <- readRDS(parameter$object) } )
	if ( class(read.obj) == "character" ) {
		object <- get(read.obj[[1]])
	}else {
		object <- read.obj
	}
	
	RunMessage("    ---> Checking Seurat object version <---    ")
	if(grepl('^2', object@version)){
		version <- "2x"
		object <- AddFdata(object, parameter$namelist)
		object <- UpdateSeuratObject(object)
	}else{
		version <- "3x"
#		if(!exists("fdata", object@misc)){
			object <- AddFdata(object, parameter$namelist)
			object@misc$fdata <- object@misc$fData
#		}
	}

	## rename cluster colname to seurat_clusters, and set ident
	if(parameter$metadata$cluster.col != "seurat_clusters"){
		colnames(object@meta.data)[which(colnames(object@meta.data)== "seurat_clusters")] <- "orig.seurat_clusters"
		colnames(object@meta.data)[which(colnames(object@meta.data)== parameter$metadata$cluster.col)] <- "seurat_clusters"
	}
	object@meta.data$seurat_clusters <- factor(object@meta.data$seurat_clusters)
	Idents(object) <- "seurat_clusters"
	## rename sample colname to orig.ident
	if(parameter$metadata$sample.col != "orig.ident"){
		colnames(object@meta.data)[which(colnames(object@meta.data)== "orig.ident")] <- "orig.ident.old"
		colnames(object@meta.data)[which(colnames(object@meta.data)== parameter$metadata$sample.col)] <- "orig.ident"
	}
	object@meta.data$orig.ident <- factor(object@meta.data$orig.ident)
	object@meta.data$sample <- object@meta.data$orig.ident

	RunMessage("    ---> Checking Input Samples <---    ")
	samples.use <- NULL
	if( is.null(parameter$samples.use) ){
		samples.use <- levels(object@meta.data$orig.ident)
	}else{
		samples.use <- parameter$samples.use
	}
	if(!is.null(parameter$metadata$groups$orig.label)){
		samples.use <- unlist(parameter$metadata$groups$orig.label)
	}
	
	RunMessage("    ---> Checking Input Clusters <---    ")
	clusters.use <- NULL
	if( is.null(parameter$clusters.use) ){
		clusters.use <- levels(object@meta.data$seurat_clusters)
	}else{
		clusters.use <- parameter$clusters.use
	}
	if(!is.null(parameter$metadata$clusters$orig.label)){
		clusters.use <- unlist(parameter$metadata$clusters$orig.label)
	}
	object@meta.data$cluster <- object@meta.data$seurat_clusters

	RunMessage("    ---> Subset Object <---    ")
	## Subset Object
	cells <- Cells(object)[object@meta.data$sample %in% samples.use & object@meta.data$cluster %in% clusters.use ]
	object <- object[, cells]
	object@meta.data <- droplevels(object@meta.data)
	object@active.ident <- droplevels(object@active.ident)
	samples.use <- samples.use[samples.use  %in% unique(object@meta.data$sample)]
	clusters.use <- clusters.use[clusters.use  %in% unique(object@meta.data$cluster)]
	## Sample and Cluster Info printed
	RunMessage(paste0("    --> Selected Sample: [ ", paste(samples.use, collapse = " | "), " ]"))
	RunMessage(paste0("    --> Selected Cluster: [ ", paste(clusters.use, collapse = " | "), " ]"))
	if(!is.null(parameter$metadata$groups$orig.label) | !is.null(parameter$samples.use)) object@meta.data$sample <- factor(object@meta.data$sample, levels = samples.use)
	
	
	## Add Group Info
	if( !is.null(parameter$metadata$groups$orig.label) ){
		RunMessage("    ---> Add Group Info <---    ")
		object@meta.data$group <- NULL
		object <- TransMetaData(object, orig.labels = parameter$metadata$groups$orig.label, new.labels = parameter$metadata$groups$new.label, orig.col.name = "orig.ident", new.col.name = "group")
	}

	## new Cluster was transformed
	if( !is.null(parameter$metadata$clusters$orig.label) ){
		RunMessage("    ---> Add new Cluster Info <---    ")
		object@meta.data$cluster <- NULL
		object <- TransMetaData(object, orig.labels = parameter$metadata$clusters$orig.label, new.labels = parameter$metadata$clusters$new.label, orig.col.name = "seurat_clusters", new.col.name = "cluster", setIdent = parameter$metadata$clusters$set.ident)
	}else{
		object@meta.data$seurat_clusters <- factor(object@meta.data$seurat_clusters, levels = clusters.use)
		object@meta.data$cluster <- object@meta.data$seurat_clusters
	}
	RunMessage("    ---> Get expession counts <---    ")
	object@misc$Expr <- object@assays$RNA@counts
	if(version == "2x"){
		rownames(object@misc$fData) <- as.vector(object@misc$fData$nounderscore)
		rownames(object@misc$Expr) <- as.vector(object@misc$fData[rownames(object@misc$Expr), "merge_name"])
	}else{
		if(! all(rownames(object@misc$Expr) %in% rownames(object@misc$fdata))){
			rownames(object@misc$fdata) <- object@misc$fdata$dash
		}
		nudup.id <- rownames(object@misc$fdata[!duplicated(object@misc$fdata$name), ])
		nudup.id <- intersect(nudup.id, rownames(object@misc$Expr))
		object@misc$Expr <- object@misc$Expr[nudup.id, ]
		rownames(object@misc$Expr) <- as.vector(object@misc$fdata[rownames(object@misc$Expr),"name"])
	}
	print(str(object@misc$Expr))
	return(object)
}

WriteTable <- function(data, outfile = NULL){
	write.table(data, sep="\t", quote=F, col.names=T, row.names=F, file= outfile)
}

Hclust <- function(data, method = "ave", labels.return = T){
	dat <- t(data)
	hc <- hclust(dist(dat),method=method)
	order.labels <- hc$labels[hc$order]
	data <- data[, order.labels]
	if(labels.return){
		return(order.labels)
	}else{
		return(data)
	}
}

rowmeans <- function(mat, names){
	data <- NULL
	if(length(names)==1){
		data <- mat[,names]
	}else{
		data <- rowMeans(mat[,names])
	}
	return(data)
}

AucHeatmap <- function(
		scenicOptions,
		plot.type = "all",
		group_by = NULL,
		group.use = NULL,
		TF.use = NULL,
		group.plot = FALSE,
		add.col.bar = FALSE,
		cluster_cols = FALSE,
		cluster_rows = TRUE,
		do.preclust = TRUE,
		row.cex = 8,
		col.cex = 12,
		cellwidth = NA,
		cellheight = NA,
		color.use = colorRampPalette(c("blue","white","red"))(100),
		border_color = NA,
		breaks = seq(-3, 3, length.out = 100),
		main = NA,
		show_rownames = T,
		show_colnames = T,
		border =F,
		use.log = F,
		pdf.width = NA,
		pdf.height = NA,
		outexp = F,
		out = NULL
){
	## Get TF-TARGET AUC (regulon => rownames, cell.names => colnames)
	## plot type => { all => all, centred = > centred, upstream => upstream}
	group.by <- NULL
	if(!is.null(group_by)){
		group.by <- group_by
		group_by <- group_by[1]
	}
	regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
	regulonActivity <- getAUC(regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ])
#	rownames(regulonActivity) <- gsub( "g)", " genes)", rownames(regulonActivity) )
	outtype <- "all"
	if(plot.type == "centred"){
		centred <- grep( pattern = 'extended', x = rownames(regulonActivity), invert = F )
		regulonActivity <- regulonActivity[centred, ]
		outtype <- "centred"
	}else if(plot.type == "upstream"){
		upstream <- grep( pattern = 'extended', x = rownames(regulonActivity), invert = T )
		regulonActivity <- regulonActivity[upstream, ]
		outtype <- "upstream"
	}else if(plot.type != "all"){
		StopMessage("plot type must be any one of [all | centred | upstream]")
	}
	
	cluster.type <- NULL
	if(cluster_cols){
		cluster.type <- c(cluster.type, "colcluster")
	}
	if(cluster_rows){
		cluster.type <- c(cluster.type, "rowcluster")
	}
	if( !is.null(cluster.type)){
		cluster.type <- paste(cluster.type, collapse=".")
	}else{
		cluster.type <- ""
	}

	if(use.log){
		regulonActivity <- log2(regulonActivity+1)
	}
	if(!is.null(TF.use)){
		TF.use <- TF.use[TF.use %in% rownames(regulonActivity)]
		regulonActivity <- regulonActivity[TF.use,]
	}
	
	cells.use <- colnames(regulonActivity)
	## check group_by and group.use
	if(is.null(group_by) & !is.null(group.use)){
		StopMessage("When using group.use parameter, group_by parameter can not be null!")
	}

	## load cellInfo and colVars
#	colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions,'colVars'))
	cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions,'cellInfo'))
	if(!is.null(group_by)){
		if(!is.factor(cellInfo[[group_by]])) cellInfo[[group_by]] <- factor(cellInfo[[group_by]])
		if(is.null(group.use)) group.use <- as.vector(unique(cellInfo[[group_by]]))
		cells <- NULL
		for (i in group.use){
			cells.tmp <- rownames(cellInfo)[cellInfo[[group_by]] == i]
			if(do.preclust){
				auc <- regulonActivity[,cells.tmp, drop=F]
				if (length(cells.tmp) >= 2) cells.tmp <- Hclust(auc)
			}
			cells <- c(cells, cells.tmp)
			rm(cells.tmp)
		}
		tmp <- intersect(cells, cells.use)
		if( length(tmp)  == 0){
			StopMessage("cells subsetted by group.use, have no common cells with cells.use\n")
		}
		cells.use <- tmp
		rm(tmp, i, cells)
		regulonActivity <- regulonActivity[, cells.use]
		cellInfo <- cellInfo[ cells.use, ]
		cellInfo <- droplevels(cellInfo)
		cellInfo[[group_by]] <- factor(cellInfo[[group_by]], levels=group.use)
		
#		color.group <- rainbow(length(group.use))
#		names(color.group) <- group.use
		annot_col <- cellInfo[, group.by, drop=FALSE]
#		if(group_by == "group"){
#			ann_colors = list(group = color.group)
#		}else if(group_by == "cluster"){
#			ann_colors = list(cluster = color.group)
#		}else if(group_by == "sample"){
#			ann_colors = list(sample = color.group)
#		}else{
#			ann_colors = list(type = color.group)
#		}
		pdf.width <- max(7, 0.1 * length(group.use))
		outprf <- NULL
		if(length(group.by) >1){
			outprf <- "RegulonActivity."
		}else{
			outprf <- paste0( group_by,".RegulonActivity.")
		}
		if(group.plot){
#			regulonActivity <- sapply(split(rownames(cellInfo), cellInfo[[group_by]]), function(cells) rowMeans(regulonActivity[ ,cells])) has bugs!
			regulonActivity <- sapply(split(rownames(cellInfo), cellInfo[[group_by]]), function(cells) rowmeans(regulonActivity, cells))
			pdf.height <- max(7, 0.1 * nrow(regulonActivity))
			p <- pheatmap(regulonActivity,
				scale="row",
				show_colnames=show_colnames,
				show_rownames = show_rownames,
				fontsize_row = row.cex,
				fontsize_col = col.cex,
				cluster_cols=cluster_cols,
				cluster_rows=cluster_rows, 
				width=pdf.width,
				height=pdf.height,
				filename=paste0(outprf, outtype, ".average.pdf"),
				breaks = breaks,
				border_color=border_color,
				border = border,
				color = color.use
				)
			if(outexp){
				if(cluster_cols & cluster_rows){
					regulonActivity <- regulonActivity[p$tree_row$order, p$tree_col$order]
				}else if(cluster_cols & (!cluster_rows)){
					regulonActivity <- regulonActivity[, p$tree_col$order]
				}else if(cluster_rows & (!cluster_cols)){
					regulonActivity <- regulonActivity[p$tree_row$order,]
				}
				WriteTable(cbind(id=rownames(regulonActivity), regulonActivity), outfile = paste0(outprf, outtype, ".average.xls" ))
			
			}
		}else{
			pdf.height <- max(7, 0.1 * nrow(regulonActivity))
			p <- pheatmap(regulonActivity,
				scale="row",
				show_colnames=F,
				show_rownames = T,
				annotation_col = annot_col,
				fontsize_row = row.cex,
				fontsize_col = col.cex,
				cluster_cols=cluster_cols,
				cluster_rows=cluster_rows,
				width=pdf.width,
				height=pdf.height,
				breaks = breaks,
				color = color.use,
				border_color = border_color,
				border = border,
				filename=paste0(outprf, outtype, ".pdf")
				)
				if(outexp){
					if(cluster_cols & cluster_rows){
						regulonActivity <- regulonActivity[p$tree_row$order, p$tree_col$order]
					}else if(cluster_cols & (!cluster_rows)){
						regulonActivity <- regulonActivity[, p$tree_col$order]
					}else if(cluster_rows & (!cluster_cols)){
						regulonActivity <- regulonActivity[p$tree_row$order,]
					}
					WriteTable(cbind(id=rownames(regulonActivity), regulonActivity), outfile = paste0( outprf, outtype, ".xls" ))
					
					if(outtype == "all"){
						regulonActivity.demo <- regulonActivity[, 1: min(10, ncol(regulonActivity))]
						WriteTable(cbind(id=rownames(regulonActivity.demo), regulonActivity.demo), outfile = paste0( outprf, outtype, ".demo.xls" ))
					}
				}
		}
	}else{
		pdf.height <- max(7, 0.1 * nrow(regulonActivity))
		p <- pheatmap(regulonActivity,
				scale="row",
				show_colnames=show_colnames,
				show_rownames = show_rownames,
				fontsize_row = row.cex,
				fontsize_col = col.cex,
				cluster_cols=cluster_cols,
				cluster_rows=cluster_rows,
				breaks = breaks,
				color = color.use,
				width=7,
				height=pdf.height,
				filename=paste0("RegulonActivity.", outtype, ".pdf"),
				border_color = border_color,
				border = border
				)
		if(outexp){
			if(cluster_cols & cluster_rows){
				regulonActivity <- regulonActivity[p$tree_row$order, p$tree_col$order]
			}else if(cluster_cols & (!cluster_rows)){
				regulonActivity <- regulonActivity[, p$tree_col$order]
			}else if(cluster_rows & (!cluster_cols)){
				regulonActivity <- regulonActivity[p$tree_row$order,]
			}
			WriteTable(cbind(id=rownames(regulonActivity), regulonActivity), outfile = paste0( "RegulonActivity.", outtype, ".xls" ))
			if(outtype == "all"){
				regulonActivity.demo <- regulonActivity[, 1: min(10, ncol(regulonActivity))]
				WriteTable(cbind(id=rownames(regulonActivity.demo), regulonActivity.demo), outfile = paste0( "RegulonActivity.", outtype, ".demo.xls" ))
			}
		}
	}
}

BinaryAucHeatmap <- function(
		scenicOptions,
		addthreshold=FALSE,
		minPerc = NULL,
		plot.type = "all",
		group_by = NULL,
		group.use = NULL,
		TF.use = NULL,
		group.plot = FALSE,
		add.col.bar = FALSE,
		cluster_cols = FALSE,
		cluster_rows = TRUE,
		do.preclust = TRUE,
		row.cex = 8,
		col.cex = 12,
		cellwidth = NA,
		cellheight = NA,
		color.use = ifelse(group.plot, colorRampPalette(c("white","pink","red"))(100), colorRampPalette(c("white","black"))(100)),
		border_color = NA,
		border = F,
		breaks = seq(0, 1, length.out = 100),
		main = NA,
		show_rownames = T,
		show_colnames = T,
		pdf.height = NA,
		outexp = F,
		out = NULL
){
	group.by <- NULL
	if(!is.null(group_by)){
		group.by <- group_by
		group_by <- group_by[1]
	}
	cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, 'cellInfo'))
	binaryregulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
	cellInfo <- cellInfo[which(rownames(cellInfo) %in% colnames(binaryregulonActivity)),, drop=FALSE]

	outtype <- "all"
	if(plot.type == "centred"){
		centred <- grep( pattern = 'extended', x = rownames(binaryregulonActivity), invert = F )
		binaryregulonActivity <- binaryregulonActivity[centred, ]
		outtype <- "centred"
	}else if(plot.type == "upstream"){
		upstream <- grep( pattern = 'extended', x = rownames(binaryregulonActivity), invert = T )
		binaryregulonActivity <- binaryregulonActivity[upstream, ]
		outtype <- "upstream"
	}else if(plot.type != "all"){
		StopMessage("plot type must be any one of [all | centred | upstream]")
	}
	
	if(!is.null(TF.use)){
		TF.use <- TF.use[TF.use %in% rownames(binaryregulonActivity)]
		binaryregulonActivity <- binaryregulonActivity[TF.use,]
	}
	
	cells.use <- colnames(binaryregulonActivity)
	## check group_by and group.use
	if(is.null(group_by) & !is.null(group.use)){
		StopMessage("When using group.use parameter, group_by parameter can not be null!")
	}
	if(!is.null(group_by)){
		if(!is.factor(cellInfo[[group_by]])) cellInfo[[group_by]] <- factor(cellInfo[[group_by]])
		if(is.null(group.use)) group.use <- as.vector(unique(cellInfo[[group_by]]))
		cells <- NULL
		for (i in group.use){
			cells.tmp <- rownames(cellInfo)[cellInfo[[group_by]] == i]
			if(do.preclust){
				auc <- binaryregulonActivity[,cells.tmp, drop=F]
				if (length(cells.tmp) >= 2) cells.tmp <- Hclust(auc)
			}
			cells <- c(cells, cells.tmp)
			rm(cells.tmp)
		}
		tmp <- intersect(cells, cells.use)
		if( length(tmp)  == 0){
			StopMessage("cells subsetted by group.use, have no common cells with cells.use\n")
		}
		cells.use <- tmp
		rm(tmp, i, cells)
		binaryregulonActivity <- binaryregulonActivity[, cells.use]
		cellInfo <- cellInfo[ cells.use, ]
		cellInfo <- droplevels(cellInfo)
		cellInfo[[group_by]] <- factor(cellInfo[[group_by]], levels=group.use)
#		color.group <- rainbow(length(group.use))
#		names(color.group) <- group.use
#		annot_col <- cellInfo[, c(group_by), drop=FALSE]
		annot_col <- cellInfo[, group.by, drop=FALSE]
#		if(group_by == "group"){
#			ann_colors = list(group = color.group)
#		}else if(group_by == "cluster"){
#			ann_colors = list(cluster = color.group)
#		}else if(group_by == "sample"){
#			ann_colors = list(sample = color.group)
#		}else{
#			ann_colors = list(type = color.group)
#		}


		pdf.width <- max(7, 0.1 * length(group.use))
		outprf <- NULL
		if(length(group.by)>1){
			outprf <- "BinaryRegulonActivity."
		}else{
			outprf <- paste0(group_by, ".BinaryRegulonActivity.")
		}
		if(group.plot){
#			binaryregulonActivity <- sapply(split(rownames(cellInfo), cellInfo[[group_by]]), function(cells) rowMeans(binaryregulonActivity[ ,cells])) has bugs
			binaryregulonActivity <- sapply(split(rownames(cellInfo), cellInfo[[group_by]]), function(cells) rowmeans(binaryregulonActivity ,cells))
			if(!is.null(minPerc)) binaryregulonActivity <- binaryregulonActivity[which(rowSums(binaryregulonActivity>minPerc)>0),]
			pdf.height <- max(7, 0.1 * nrow(binaryregulonActivity))
			p <- pheatmap(
					binaryregulonActivity,
					scale="none",
					show_colnames = show_colnames,
					show_rownames = show_rownames,
					fontsize_row = row.cex,
					fontsize_col = col.cex,
					cluster_cols=cluster_cols,
					cluster_rows=cluster_rows,
					width=pdf.width,
					height=pdf.height,
					filename=paste0(outprf, outtype, ".average.pdf"),
					breaks = breaks,
					color = colorRampPalette(c("white","pink","red"))(100),
					border_color=border_color,
					border=border
					)
			if(outexp){
				if(cluster_cols & cluster_rows){
					binaryregulonActivity <- binaryregulonActivity[p$tree_row$order, p$tree_col$order]
				}else if(cluster_cols & (!cluster_rows)){
					binaryregulonActivity <- binaryregulonActivity[, p$tree_col$order]
				}else if(cluster_rows & (!cluster_cols)){
					binaryregulonActivity <- binaryregulonActivity[p$tree_row$order, ]
				}
				WriteTable(cbind(id=rownames(binaryregulonActivity), binaryregulonActivity), outfile = paste0( outprf, outtype, ".average.xls" ))
			}
		}else{
			pdf.height <- max(7, 0.1 * nrow(binaryregulonActivity))
			p <- pheatmap(binaryregulonActivity,
					scale="none",
					show_colnames=F,
					show_rownames = T,
					annotation_col = annot_col,
					fontsize_row = row.cex,
					fontsize_col = col.cex,
					cluster_cols=cluster_cols,
					cluster_rows=cluster_rows,
					width=pdf.width,
					height=pdf.height,
					breaks = breaks,
					color = colorRampPalette(c("white","black"))(100),
					border_color=border_color,
					border=border,
					filename=paste0(outprf, outtype, ".pdf")
					)

			if(outexp){
				if(cluster_cols & cluster_rows){
					binaryregulonActivity <- binaryregulonActivity[p$tree_row$order, p$tree_col$order]
				}else if(cluster_cols & (!cluster_rows)){
					binaryregulonActivity <- binaryregulonActivity[, p$tree_col$order]
				}else if(cluster_rows & (!cluster_cols)){
					binaryregulonActivity <- binaryregulonActivity[p$tree_row$order,]
				}
				
				if(addthreshold){
					thresh_file <- getIntName(scenicOptions,'aucell_thresholdsTxt')
					theshold <- read.table(file=thresh_file, sep="\t", row.names=1, header=T, check.names=F)
					theshold <- theshold[rownames(binaryregulonActivity), ]
					binaryregulonActivity <- cbind(threshold=theshold$threshold, binaryregulonActivity)
				}
				WriteTable(cbind(id=rownames(binaryregulonActivity), binaryregulonActivity), outfile = paste0( outprf, outtype, ".xls" ))

				if( outtype == "all"){
					binaryregulonActivity.demo <- binaryregulonActivity[,1:min(10, ncol(binaryregulonActivity))]
					WriteTable(cbind(id=rownames(binaryregulonActivity.demo), binaryregulonActivity.demo), outfile = paste0( outprf, outtype, ".demo.xls" ))
				}
			}
		}
	}else{
		pdf.height <- max(7, 0.1 * nrow(binaryregulonActivity))
		p <- pheatmap(binaryregulonActivity,
				scale="none",
				show_colnames=show_colnames,
				show_rownames = show_rownames,
				fontsize_row = row.cex,
				fontsize_col = col.cex,
				cluster_cols=cluster_cols,
				cluster_rows=cluster_rows,
				breaks = breaks,
				color = colorRampPalette(c("white","black"))(100),
				border_color=border_color,
				border=border,
				width =7,
				height = pdf.height,
				filename=paste0("BinaryRegulonActivity.", outtype, ".pdf")
				)
		if(outexp){
			if(cluster_cols & cluster_rows){
				binaryregulonActivity <- binaryregulonActivity[p$tree_row$order, p$tree_col$order]
			}else if(cluster_cols & (!cluster_rows)){
				binaryregulonActivity <- binaryregulonActivity[, p$tree_col$order]
			}else if(cluster_rows & (!cluster_cols)){
				binaryregulonActivity <- binaryregulonActivity[p$tree_row$order,]
			}
			WriteTable(cbind(id=rownames(binaryregulonActivity), binaryregulonActivity), outfile = paste0("BinaryRegulonActivity.", outtype, ".xls" ))
			if( outtype == "all"){
				binaryregulonActivity.demo <- binaryregulonActivity[,1:min(10, ncol(binaryregulonActivity))]
				WriteTable(cbind(id=rownames(binaryregulonActivity.demo), binaryregulonActivity.demo), outfile = paste0( "BinaryRegulonActivity.", outtype, ".demo.xls" ))
			}
		}
	}
}



ExportRegulonsList <- function(regulonsRds, outfile){
	regulons <- readRDS(regulonsRds)
	names(regulons) <- gsub(" (.*)","",x=names(regulons))
	out <- NULL
	for(i in names(regulons)){
		genes <- as.vector(unlist(regulons[i]))
		for(j in genes){
			out <- rbind(out,c(i,j))
		}
	}
	colnames(out) <- c("Regulons","Targets")
	write.table(out,file=outfile,sep="\t",quote=F,col.names=T,row.names=F)
}




AttachScenic <- function(outdir = NULL){
	data <- list()
	scenicOptions <- if(file.exists(paste0(outdir, "/1.CreateScenicOptions/scenicOptions.Rds"))){
		Message("---> Attach scenicOption <---")
		readRDS(paste0(outdir, "/1.CreateScenicOptions/scenicOptions.Rds"))
	}else{
		NULL
	}

	data[["scenicOptions"]] <- scenicOptions

	object <- if(file.exists(paste0(outdir, "/1.CreateScenicOptions/obj.Rds"))){
		Message("---> Attach object <---")
		readRDS(paste0(outdir, "/1.CreateScenicOptions/obj.Rds"))
	}else{
		NULL
	}
	data[["object"]] <- object

	exprMat <- if(file.exists(paste0(outdir, "/2.FilterGene/exprMat.Rds"))){
		Message("---> Attach exprMat <---")
		readRDS(paste0(outdir, "/2.FilterGene/exprMat.Rds"))
	}else{
		NULL
	}
	data[["exprMat"]] <- exprMat
	data
}

RunGenie3 <- function (exprMat, scenicOptions, genes.use = NULL, index = NULL, nCores=1, ...)
{
	#nCores <- getSettings(scenicOptions, "nCores")
	nCores <- nCores
    if (is.data.frame(exprMat)) {
        supportedClasses <- paste(gsub("AUCell_buildRankings,", 
            "", methods("AUCell_buildRankings")), collapse = ", ")
        supportedClasses <- gsub("-method", "", supportedClasses)
        stop("'exprMat' should be one of the following classes: ", 
            supportedClasses, "(data.frames are not supported. Please, convert the expression matrix to one of these classes.)")
    }
    if (any(table(rownames(exprMat)) > 1)) 
        stop("The rownames (gene id/name) in the expression matrix should be unique.")
    allTFs <- getDbTfs(scenicOptions)
    inputTFs <- allTFs[allTFs %in% rownames(exprMat)]
    percMatched <- length(inputTFs)/length(allTFs)
    if (getSettings(scenicOptions, "verbose")) 
        message("Using ", length(inputTFs), " TFs as potential regulators...")
    if (percMatched < 0.4) 
        warning("Only ", round(percMatched * 100), "% of the ", 
            length(allTFs), " TFs in the database were found in the dataset. Do they use the same gene IDs?\n")
    weightMatrix <- NULL
    if (getSettings(scenicOptions, "verbose")) message("Running GENIE3")
    set.seed(getSettings(scenicOptions, "seed"))
    weightMatrix <- GENIE3::GENIE3(exprMat, regulators = inputTFs, 
        nCores = nCores, targets = genes.use, verbose=T, ...)
    fileName <- gsub(".Rds$", paste0("_part_", index, ".Rds"), 
        getIntName(scenicOptions, "genie3wm"))
    saveRDS(weightMatrix, file = fileName)
}
environment(RunGenie3) <- environment(runGenie3)

MergeGenie3 <- function(scenicOptions){
    intdir <- paste0(scenicOptions@fileNames$dir, "/3.Coexpression")
    weightMatrices <- list.files(path= intdir, pattern="weightMatrix")
    linkList_list <- list()
    for (i in 1:length(weightMatrices)) {
        weightMatrix <- readRDS(paste0(intdir, "/", weightMatrices[i]))
        linkList_list[[i]] <- GENIE3::getLinkList(weightMatrix, 
            threshold = getSettings(scenicOptions, "modules/weightThreshold"))
    }
    rm(weightMatrices)
    linkList <- do.call(rbind, linkList_list)
    colnames(linkList) <- c("TF", "Target", "weight")
    linkList <- linkList[order(linkList[, "weight"], decreasing = TRUE), 
        ]
    saveRDS(linkList, file = getIntName(scenicOptions, "genie3ll"))
    if (getSettings(scenicOptions, "verbose")) 
        message("Finished mergeing GENIE3.")
    invisible(linkList)
}

initializeScenic <- function (org = NULL, dbDir = "databases", dbs = NULL, datasetTitle = "", 
    nCores = 4, outdir = NULL) 
{
    inputDataset <- list(org = org, datasetTitle = datasetTitle, 
        cellInfo = c(paste0( outdir, "/1.CreateScenicOptions/cellInfo.Rds", NA)), colVars = c( paste0(outdir, "/1.CreateScenicOptions/colVars.Rds"), 
            NA), int_01 = c( paste0(outdir, "/1.CreateScenicOptions/cellColorNgenes.Rds", NA)))
    if (!org %in% c("mgi", "hgnc", "dmel")) 
        stop("'org' should be one of: mgi, hgnc, dmel.")
    defaultDBs <- FALSE
    if (is.null(dbs)) {
        data(defaultDbNames)
        dbs <- defaultDbNames[[org]]
    }
    dbsFound <- unlist(unname(lapply(dbs, function(x) setNames(file.exists(file.path(dbDir, 
        x)), unname(file.path(dbDir, x))))))
    if (any(!dbsFound)) {
        stop("The following RcisTarget databases were not found: ", 
            paste(paste0("\n- ", names(dbsFound[which(!dbsFound)])), 
                collapse = " "), "\nMake sure the arguments 'dbDir' and 'dbs' are correct.")
        dbs <- NULL
    }
    else {
        message("Motif databases selected: ", paste(paste0("\n  ", 
            dbs, collapse = " ")))
    }
    loadAttempt <- sapply(dbs, function(x) dbLoadingAttempt(file.path(dbDir, 
        x)))
    if (any(!loadAttempt)) 
        warning("It was not possible to load the following databses; check whether they are downloaded correctly: \n", 
            paste(dbs[which(!loadAttempt)], collapse = "\n"))
    db_mcVersion <- dbVersion(dbs)
    scenicSettings = list(dbs = dbs, dbDir = dbDir, db_mcVersion = db_mcVersion, 
        verbose = TRUE, nCores = nCores, seed = 123, devType = "pdf", 
        modules = list(weightThreshold = 0.001), regulons = list(), 
        aucell = list(smallestPopPercent = 0.25), defaultTsne = list(dims = 50, 
            perpl = 50, aucType = "AUC"), tSNE_filePrefix = "tSNE")
    scenicFiles <- list(output = c(s2_motifEnrichment = paste0( outdir, "/4.AUC/Step2_MotifEnrichment.tsv"), 
        s2_motifEnrichmentHtml = paste0( outdir, "/4.AUC/Step2_MotifEnrichment_preview.html"), 
        s2_regulonTargetsInfo = paste0( outdir, "/4.AUC/Step2_regulonTargetsInfo.tsv"), 
        s3_AUCheatmap = paste0( outdir, "/4.AUC/Step3_RegulonActivity_heatmap"), 
        s3_AUCtSNE_colAct = paste0( outdir, "/4.AUC/Step3_RegulonActivity_tSNE_colByActivity"), 
        s3_AUCtSNE_colProps = paste0( outdir, "/4.AUC/Step3_RegulonActivity_tSNE_colByCellProps"), 
        s4_boxplotBinaryActivity = paste0( outdir, "/4.AUC/Step4_BoxplotActiveCellsRegulon"), 
        s4_binaryActivityHeatmap = paste0( outdir, "/4.AUC/Step4_BinaryRegulonActivity_Heatmap_"), 
        s4_binarytSNE_colAct = paste0( outdir, "/4.AUC/Step4_BinaryRegulonActivity_tSNE_colByActivity"), 
        s4_binarytSNE_colProps = paste0( outdir, "/4.AUC/Step4_BinaryRegulonActivity_tSNE_colByCellProps"), 
        loomFile = paste0( outdir, "/4.AUC/scenic.loom")), int = list(genesKept = c(paste0(outdir, "/2.FilterGene/1.1_genesKept.Rds"), 
        TRUE), corrMat = c(paste0(outdir,"/2.FilterGene/1.2_corrMat.Rds"), TRUE), genie3wm = c(paste0(outdir, "/3.Coexpression/1.3_GENIE3_weightMatrix.Rds"), 
        FALSE), genie3ll = c(paste0(outdir, "/3.Coexpression/1.4_GENIE3_linkList.Rds"), TRUE), 
        genie3weighPlot = c(paste0(outdir, "/3.Coexpression/1.5_weightPlot"), TRUE), tfModules_asDF = c(paste0(outdir, "/4.AUC/1.6_tfModules_asDF.Rds"), 
            TRUE), tfModules_forEnrichment = c(paste0(outdir, "/4.AUC/2.1_tfModules_forMotifEnrichmet.Rds"), 
            FALSE), motifs_AUC = c(paste0(outdir , "/4.AUC/2.2_motifs_AUC.Rds"), 
            FALSE), motifEnrichment_full = c(paste0(outdir, "/4.AUC/2.3_motifEnrichment.Rds"), 
            FALSE), motifEnrichment_selfMotifs_wGenes = c( paste0(outdir, "/4.AUC/2.4_motifEnrichment_selfMotifs_wGenes.Rds"), 
            FALSE), regulonTargetsInfo = c( paste0(outdir, "/4.AUC/2.5_regulonTargetsInfo.Rds"), 
            NA), regulons = c( paste0(outdir, "/4.AUC/2.6_regulons_asGeneSet.Rds"), 
            NA), regulons_incidMat = c( paste0(outdir, "/4.AUC/2.6_regulons_asIncidMat.Rds"), 
            NA), aucell_regulons = c( paste0(outdir, "/4.AUC/3.1_regulons_forAUCell.Rds"), 
            NA), aucell_genesStatsPlot = c( paste0(outdir, "/4.AUC/3.2_aucellGenesStats"), 
            NA), aucell_rankings = c(paste0(outdir, "/4.AUC/3.3_aucellRankings.Rds"), 
            NA), aucell_regulonAUC = c(paste0(outdir, "/4.AUC/3.4_regulonAUC.Rds"), 
            NA), aucell_thresholds = c(paste0(outdir, "/4.AUC/3.5_AUCellThresholds.Rds"), 
            NA), aucell_thresholdsTxt = c(paste0(outdir, "/4.AUC/3.5_AUCellThresholds_Info.tsv"), 
            NA), aucell_binary_full = c(paste0(outdir, "/4.AUC/4.1_binaryRegulonActivity.Rds"), 
            NA), aucell_binary_nonDupl = c(paste0(outdir, "/4.AUC/4.2_binaryRegulonActivity_nonDupl.Rds"), 
            NA), aucell_regulonSelection = c(paste0(outdir, "/4.AUC/4.3_regulonSelections.Rds"), 
            NA), aucell_binaryRegulonOrder = c(paste0(outdir, "/4.AUC/4.4_binaryRegulonOrder.Rds"), 
            NA)))
    scenicFiles$output <- cbind(fileName = scenicFiles$output)
    scenicFiles$int <- do.call(rbind, scenicFiles$int)
    colnames(scenicFiles$int) <- c("fileName", "keep")
    scenicFiles$int <- scenicFiles$int[, "fileName", drop = FALSE]
#    dir.create("int", showWarnings = FALSE)
#    dir.create("output", showWarnings = FALSE)
    object <- new("ScenicOptions", inputDatasetInfo = inputDataset, 
        settings = scenicSettings, fileNames = scenicFiles)
    return(object)
}
environment(initializeScenic) <- asNamespace('SCENIC')
assignInNamespace("initializeScenic", initializeScenic, ns = "SCENIC")

AUCell_plotTSNE <- function (tSNE, exprMat = NULL, cellsAUC = NULL, thresholds = NULL, 
    reorderGeneSets = FALSE, cex = 1, alphaOn = 1, alphaOff = 0.2, 
    borderColor = adjustcolor("lightgray", alpha.f = 0.1), offColor = "lightgray", 
    plots = c("histogram", "binaryAUC", "AUC", "expression"), 
    exprCols = c("goldenrod1", "darkorange", "brown"), asPNG = FALSE, 
    ...) 
{
    if (is.null(rownames(tSNE))) 
        stop("Please, provide the cell rownames in the t-SNE")
    if (!is.matrix(tSNE) | ncol(tSNE) != 2) 
        stop("The t-SNE should be a matrix with 2 columns (cell coordinates)")
    if (any(grepl("binary", tolower(plots)))) {
        plots[grep("binary", tolower(plots))] <- "binaryAUC"
    }
    if (is.null(exprMat) && ("expression" %in% tolower(plots))) {
        plots <- plots[which(plots != "expression")]
        warning("Expression plot was requested, but no expression matrix provided.")
    }
    if (length(plots) == 0) 
        stop("Please, provide which plots to plot.")
    if (reorderGeneSets) {
        cellsAUC <- cellsAUC[orderAUC(cellsAUC), ]
    }
    if (is.logical(thresholds) && thresholds == FALSE) {
        thresholds <- FALSE
        if (any(grepl("binary", tolower(plots)))) 
            stop("Cannot plot binary AUC without calculating the thresholds.")
    }
    else {
        if (!is.null(thresholds)) {
            if (is.list(thresholds[1])) {
                if ("aucThr" %in% names(thresholds[[1]])) 
                  thresholds <- sapply(thresholds, function(x) unname(x$aucThr$selected))
                if ("threshold" %in% names(thresholds[[1]])) 
                  thresholds <- sapply(thresholds, function(x) unname(x$threshold))
            }
            if (!is.null(names(thresholds))) {
                geneSetNames <- rownames(cellsAUC)[which(rownames(cellsAUC) %in% 
                  names(thresholds))]
                cellsAUC <- cellsAUC[geneSetNames, ]
            }
            if (is.null(names(thresholds)) || length(thresholds) == 
                1) {
                thresholds <- setNames(rep(thresholds, nrow(cellsAUC)), 
                  rownames(cellsAUC))
            }
        }
    }
    if (!is.null(cellsAUC)) {
        selectedGeneSets <- rownames(cellsAUC)
    }
    else {
        selectedGeneSets <- rownames(exprMat)
        plots <- "expression"
    }
    cells_trhAssignment <- list()
    dirName <- "./"
    if (is.character(asPNG)) {
        if (!file.exists(asPNG)) 
            dir.create(asPNG)
        dirName <- paste0(asPNG, "/")
        asPNG <- TRUE
    }
    if (asPNG) {
        nCols <- length(plots)
        figsMatrix <- matrix(nrow = length(selectedGeneSets), 
            ncol = nCols)
        rownames(figsMatrix) <- selectedGeneSets
        colnames(figsMatrix) <- plots
    }
    for (geneSetName in selectedGeneSets) {
        if (is.null(thresholds) && any(c("histogram", "binaryAUC") %in% 
            plots)) {
            if (asPNG & ("histogram" %in% tolower(plots))) {
                imgFile <- paste0(geneSetName, "_histogram.png")
                pdfFile <- paste0(geneSetName, "_histogram.pdf")
                figsMatrix[geneSetName, "histogram"] <- imgFile
                pdf(paste0(dirName, pdfFile))
                a<-dev.cur()
                png(paste0(dirName, imgFile))
                dev.control("enable")
            }
            set.seed(123)
            cells_trhAssignment[[geneSetName]] <- AUCell_exploreThresholds(cellsAUC[geneSetName, 
                ], assignCells = TRUE, plotHist = ("histogram" %in% 
                tolower(plots)))[[geneSetName]]
            thisTrheshold <- cells_trhAssignment[[geneSetName]]$aucThr$selected
            thisAssignment <- cells_trhAssignment[[geneSetName]]$assignment
            if (asPNG & ("histogram" %in% tolower(plots))) 
                dev.copy(which=a)
                dev.off()
                dev.off()
        }
        else {
            if ("histogram" %in% tolower(plots)) {
                if (asPNG) {
                  imgFile <- paste0(geneSetName, "_histogram.png")
                  pdfFile <- paste0(geneSetName, "_histogram.pdf")
                  figsMatrix[geneSetName, "histogram"] <- imgFile
                  pdf(paste0(dirName, pdfFile))
                  a<-dev.cur()
                  png(paste0(dirName, imgFile))
                  dev.control("enable")
                }
                thisTrh <- as.vector(thresholds[geneSetName])
                tmp <- .auc_plotHist(auc = getAUC(cellsAUC)[geneSetName, 
                  ], gSetName = geneSetName, aucThr = min(thisTrh, 
                  1), nBreaks = 100, sub = "AUC")
                if (!is.null(thisTrh)) {
                  abline(v = thisTrh, lwd = 3, lty = 2, col = "darkorange")
                }
                if (asPNG)
                  dev.copy(which=a)
                  dev.off()
                  dev.off()
            }
            cells_trhAssignment[[geneSetName]] <- NULL
            if (!is.null(cellsAUC)) {
                thisTrheshold <- thresholds[geneSetName]
                if (is.matrix(cellsAUC)) {
                  matrixAUC <- cellsAUC
                }
                else {
                  matrixAUC <- getAUC(cellsAUC)
                }
                thisAssignment <- names(which(matrixAUC[geneSetName, 
                  ] > thisTrheshold))
                cells_trhAssignment[[geneSetName]] <- list(threshold = thisTrheshold, 
                  assignment = thisAssignment)
            }
        }
        if (any(grepl("binary", tolower(plots)))) {
            if (asPNG) {
                imgFile <- paste0(geneSetName, "_binaryAUC.png")
                pdfFile <- paste0(geneSetName, "_binaryAUC.pdf")
                figsMatrix[geneSetName, "binaryAUC"] <- imgFile
                pdf(paste0(dirName, pdfFile))
                a<-dev.cur()
                png(paste0(dirName, imgFile))
                dev.control("enable")
            }
            thrAsTxt <- ""
            if (is.numeric(thisTrheshold)) 
                thrAsTxt <- paste("Cells with AUC > ", signif(thisTrheshold, 
                  2), sep = "")
            .auc_plotBinaryTsne(tSNE, selectedCells = thisAssignment, 
                title = geneSetName, txt = thrAsTxt, cex = cex, 
                alphaOn = alphaOn, alphaOff = alphaOff, borderColor = borderColor, 
                offColor = offColor, ...)
            if (asPNG)
                dev.copy(which=a)
                dev.off()
                dev.off()
        }
        if ("auc" %in% tolower(plots)) {
            if (asPNG) {
                imgFile <- paste0(geneSetName, "_AUC.png")
                pdfFile <- paste0(geneSetName, "_AUC.pdf")
                figsMatrix[geneSetName, "AUC"] <- imgFile
                pdf(paste0(dirName, pdfFile))
                a<-dev.cur()
                png(paste0(dirName, imgFile))
                dev.control("enable")
            }
            .auc_plotGradientTsne(tSNE, cellProp = getAUC(cellsAUC)[geneSetName, 
                ], title = geneSetName, txt = "Gene set activity (AUC)", 
                cex = cex, alphaOn = alphaOn, alphaOff = alphaOff, 
                borderColor = borderColor, offColor = offColor, 
                ...)
            if (asPNG)
                dev.copy(which=a)
                dev.off()
                dev.off()
        }
        if ("expression" %in% tolower(plots)) {
            gene <- gsub("\\s\\(\\d+g)", "", geneSetName)
            gene <- gsub("_extended", "", gene)
            if (gene %in% rownames(exprMat)) {
                if (asPNG) {
                  imgFile <- paste0(geneSetName, "_expression.png")
                  pdfFile <- paste0(geneSetName, "_expression.pdf")
                  figsMatrix[geneSetName, "expression"] <- imgFile
                  pdf(paste0(dirName, pdfFile))
                  a<-dev.cur()
                  png(paste0(dirName, imgFile))
                  dev.control("enable")
                }
                .auc_plotGradientTsne(tSNE, cellProp = exprMat[gene, 
                  ], colorsForPal = exprCols, title = paste(gene, 
                  "expression"), txt = "", cex = cex, alphaOn = alphaOn, 
                  alphaOff = alphaOff, borderColor = borderColor, 
                  offColor = offColor, ...)
                if (asPNG)
                  dev.copy(which=a)
                  dev.off()
                  dev.off()
            }
        }
    }
    if (asPNG) {
        cells_trhAssignment <- c(cells_trhAssignment, figsMatrix = list(figsMatrix))
        if ("R2HTML" %in% rownames(installed.packages())) 
            asHTML(figsMatrix, dirName)
    }
    invisible(cells_trhAssignment)
}
environment(AUCell_plotTSNE) <- asNamespace('AUCell')
#environment("AUCell_plotTSNEs") <- environment("AUCell_plotTSNE")
assignInNamespace("AUCell_plotTSNE", AUCell_plotTSNE, ns = "AUCell")

plotTsne_AUCellHtml <- function (scenicOptions, exprMat, fileName, tSNE = NULL) 
{
    regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
    if (is.null(tSNE)) 
        tSNE <- readRDS(tsneFileName(scenicOptions))
    tSNE <- tSNE$Y
    thresholds <- loadInt(scenicOptions, "aucell_thresholds", 
        ifNotExists = "null")
    if (!is.null(thresholds)) 
#        plots <- c("histogram", "binaryAUC", "AUC", "expression")
        plots <- c("histogram")
    if (is.null(thresholds)) 
#        plots <- c("histogram", "AUC", "expression")
        plots <- c("histogram")
    plotsLoc <- dirname(fileName)
    plotsName <- basename(fileName)
    AUCell_plotTSNE(tSNE = tSNE, exprMat = exprMat, cellsAUC = regulonAUC, 
        thresholds = thresholds, plots = plots, asPNG = plotsName)
    AUCell_plotTSNEs(tSNE = tSNE, exprMat = exprMat, cellsAUC = regulonAUC,
          thresholds = thresholds, plots = plots, asPNG = plotsName)
    file.rename(plotsName, file.path(plotsLoc, plotsName))
    if (file.exists("index.html")) {
        file.rename("index.html", file.path(plotsLoc, paste0(plotsName, 
            ".html")))
    }
    else {
        if (!"R2HTML" %in% installed.packages()) 
            warning("R2HTML package is not installed. Cannot produce HTML overview.")
    }
    invisible(NULL)
}
environment(plotTsne_AUCellHtml) <- asNamespace('SCENIC')
assignInNamespace("plotTsne_AUCellHtml", plotTsne_AUCellHtml, ns = "SCENIC")

RunRegulonHeatmap <- function(scenicOptions, group_by=NULL, plot.type=c("all", "centred", "upstream"), addthreshold=FALSE){
	cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions,'cellInfo'))
	if(!is.null(group_by)) group.use <- levels(cellInfo[[group_by[1]]])
	Message('--->regulonActivtiy <---')


	for (i in plot.type){
		if(length(group_by) >1){
			AucHeatmap(scenicOptions, plot.type = i, group_by = group_by, show_rownames=T, show_colnames=F, group.plot=F, cluster_cols =F, cluster_rows=T, do.preclust =T, outexp=T, group.use = group.use)
		}else{
			if(length(group.use) > 1){
				Message(paste0('---> Draw ', group_by, " ", i, ' AUC <---'))
					AucHeatmap(scenicOptions, plot.type = i, group_by = group_by, group.plot=T, cluster_cols =F, cluster_rows=T, do.preclust =T, outexp=T, group.use = group.use)
					AucHeatmap(scenicOptions, plot.type = i, group_by = group_by, group.plot=F, cluster_cols =F, cluster_rows=T, do.preclust =T, outexp=T, group.use = group.use)
			}
		}
	}
	
	Message('---> binaryregulonActivtiy <---')
	
	for (i in plot.type){
		if(length(group_by)>1){
			BinaryAucHeatmap(scenicOptions, plot.type = i, group_by = group_by, show_rownames=T, show_colnames=F, group.plot=F, cluster_cols =F, cluster_rows=T, do.preclust =T, outexp=T, addthreshold=TRUE)
		}else{
			if(length(group.use) > 1){
				Message(paste0('---> Draw ', group_by, " ", i, ' binaryAUC <---'))
					BinaryAucHeatmap(scenicOptions, plot.type = i, group_by = group_by, group.plot=T, cluster_cols =F, cluster_rows=T, do.preclust =T, outexp=T, group.use = group.use)
					BinaryAucHeatmap(scenicOptions, plot.type = i, group_by = group_by, group.plot=F, cluster_cols =F, cluster_rows=T, do.preclust=T, outexp=T, group.use = group.use)
			}
		}
	}
}


AUCell_plotTSNEs <- function (tSNE, exprMat = NULL, cellsAUC = NULL, thresholds = NULL, 
    reorderGeneSets = FALSE, cex = 1, alphaOn = 1, alphaOff = 0.2, 
    borderColor = adjustcolor("lightgray", alpha.f = 0.1), offColor = "lightgray", 
    plots = c("histogram", "binaryAUC", "AUC", "expression"), 
    exprCols = c("goldenrod1", "darkorange", "brown"), asPNG = FALSE, 
    ...) 
{
    if (is.null(rownames(tSNE))) 
        stop("Please, provide the cell rownames in the t-SNE")
    if (!is.matrix(tSNE) | ncol(tSNE) != 2) 
        stop("The t-SNE should be a matrix with 2 columns (cell coordinates)")
    if (any(grepl("binary", tolower(plots)))) {
        plots[grep("binary", tolower(plots))] <- "binaryAUC"
    }
    if (is.null(exprMat) && ("expression" %in% tolower(plots))) {
        plots <- plots[which(plots != "expression")]
        warning("Expression plot was requested, but no expression matrix provided.")
    }
    if (length(plots) == 0) 
        stop("Please, provide which plots to plot.")
    if (reorderGeneSets) {
        cellsAUC <- cellsAUC[orderAUC(cellsAUC), ]
    }
    if (is.logical(thresholds) && thresholds == FALSE) {
        thresholds <- FALSE
        if (any(grepl("binary", tolower(plots)))) 
            stop("Cannot plot binary AUC without calculating the thresholds.")
    }
    else {
        if (!is.null(thresholds)) {
            if (is.list(thresholds[1])) {
                if ("aucThr" %in% names(thresholds[[1]])) 
                  thresholds <- sapply(thresholds, function(x) unname(x$aucThr$selected))
                if ("threshold" %in% names(thresholds[[1]])) 
                  thresholds <- sapply(thresholds, function(x) unname(x$threshold))
            }
            if (!is.null(names(thresholds))) {
                geneSetNames <- rownames(cellsAUC)[which(rownames(cellsAUC) %in% 
                  names(thresholds))]
                cellsAUC <- cellsAUC[geneSetNames, ]
            }
            if (is.null(names(thresholds)) || length(thresholds) == 
                1) {
                thresholds <- setNames(rep(thresholds, nrow(cellsAUC)), 
                  rownames(cellsAUC))
            }
        }
    }
    if (!is.null(cellsAUC)) {
        selectedGeneSets <- rownames(cellsAUC)
    }
    else {
        selectedGeneSets <- rownames(exprMat)
        plots <- "expression"
    }
    cells_trhAssignment <- list()
    dirName <- "./"
    if (is.character(asPNG)) {
        if (!file.exists(asPNG)) 
            dir.create(asPNG)
        dirName <- paste0(asPNG, "/")
        asPNG <- TRUE
    }
    if (asPNG) {
        nCols <- length(plots)
        figsMatrix <- matrix(nrow = length(selectedGeneSets), 
            ncol = nCols)
        rownames(figsMatrix) <- selectedGeneSets
        colnames(figsMatrix) <- plots
    }
    for (geneSetName in selectedGeneSets) {
        if (is.null(thresholds) && any(c("histogram", "binaryAUC") %in% 
            plots)) {
            if (asPNG & ("histogram" %in% tolower(plots))) {
                pdfFile <- paste0(geneSetName, "_histogram.pdf")
                figsMatrix[geneSetName, "histogram"] <- pdfFile
                pdf(paste0(dirName, pdfFile))
            }
            set.seed(123)
            cells_trhAssignment[[geneSetName]] <- AUCell_exploreThresholds(cellsAUC[geneSetName, 
                ], assignCells = TRUE, plotHist = ("histogram" %in% 
                tolower(plots)))[[geneSetName]]
            thisTrheshold <- cells_trhAssignment[[geneSetName]]$aucThr$selected
            thisAssignment <- cells_trhAssignment[[geneSetName]]$assignment
            if (asPNG & ("histogram" %in% tolower(plots))) 
                dev.off()
        }
        else {
            if ("histogram" %in% tolower(plots)) {
                if (asPNG) {
                  pdfFile <- paste0(geneSetName, "_histogram.pdf")
                  figsMatrix[geneSetName, "histogram"] <- pdfFile
                  pdf(paste0(dirName, pdfFile))
                }
                thisTrh <- as.vector(thresholds[geneSetName])
                tmp <- .auc_plotHist(auc = getAUC(cellsAUC)[geneSetName, 
                  ], gSetName = geneSetName, aucThr = min(thisTrh, 
                  1), nBreaks = 100, sub = "AUC")
                if (!is.null(thisTrh)) {
                  abline(v = thisTrh, lwd = 3, lty = 2, col = "darkorange")
                }
                if (asPNG)
                  dev.off()
            }
            cells_trhAssignment[[geneSetName]] <- NULL
            if (!is.null(cellsAUC)) {
                thisTrheshold <- thresholds[geneSetName]
                if (is.matrix(cellsAUC)) {
                  matrixAUC <- cellsAUC
                }
                else {
                  matrixAUC <- getAUC(cellsAUC)
                }
                thisAssignment <- names(which(matrixAUC[geneSetName, 
                  ] > thisTrheshold))
                cells_trhAssignment[[geneSetName]] <- list(threshold = thisTrheshold, 
                  assignment = thisAssignment)
            }
        }
        if (any(grepl("binary", tolower(plots)))) {
            if (asPNG) {
                pdfFile <- paste0(geneSetName, "_binaryAUC.pdf")
                figsMatrix[geneSetName, "binaryAUC"] <- pdfFile
                pdf(paste0(dirName, pdfFile))
            }
            thrAsTxt <- ""
            if (is.numeric(thisTrheshold)) 
                thrAsTxt <- paste("Cells with AUC > ", signif(thisTrheshold, 
                  2), sep = "")
            .auc_plotBinaryTsne(tSNE, selectedCells = thisAssignment, 
                title = geneSetName, txt = thrAsTxt, cex = cex, 
                alphaOn = alphaOn, alphaOff = alphaOff, borderColor = borderColor, 
                offColor = offColor, ...)
            if (asPNG)
                dev.off()
        }
        if ("auc" %in% tolower(plots)) {
            if (asPNG) {
                pdfFile <- paste0(geneSetName, "_AUC.pdf")
                figsMatrix[geneSetName, "AUC"] <- pdfFile
                pdf(paste0(dirName, pdfFile))
            }
            .auc_plotGradientTsne(tSNE, cellProp = getAUC(cellsAUC)[geneSetName, 
                ], title = geneSetName, txt = "Gene set activity (AUC)", 
                cex = cex, alphaOn = alphaOn, alphaOff = alphaOff, 
                borderColor = borderColor, offColor = offColor, 
                ...)
            if (asPNG)
                dev.off()
        }
        if ("expression" %in% tolower(plots)) {
            gene <- gsub("\\s\\(\\d+g)", "", geneSetName)
            gene <- gsub("_extended", "", gene)
            if (gene %in% rownames(exprMat)) {
                if (asPNG) {
                  pdfFile <- paste0(geneSetName, "_expression.pdf")
                  figsMatrix[geneSetName, "expression"] <- pdfFile
                  pdf(paste0(dirName, pdfFile))
                }
                .auc_plotGradientTsne(tSNE, cellProp = exprMat[gene, 
                  ], colorsForPal = exprCols, title = paste(gene, 
                  "expression"), txt = "", cex = cex, alphaOn = alphaOn, 
                  alphaOff = alphaOff, borderColor = borderColor, 
                  offColor = offColor, ...)
                if (asPNG)
                  dev.off()
            }
        }
    }
#    if (asPNG) {
#        cells_trhAssignment <- c(cells_trhAssignment, figsMatrix = list(figsMatrix))
#        if ("R2HTML" %in% rownames(installed.packages())) 
#            asHTML(figsMatrix, dirName)
#    }
    invisible(cells_trhAssignment)
}
environment(AUCell_plotTSNEs) <- environment(AUCell_plotTSNE)
#assignInNamespace("AUCell_plotTSNE", AUCell_plotTSNE, ns = "AUCell")
getRegulon2Targets <- function(scenicOptions){
	aucell_regulons <- loadInt(scenicOptions, 'aucell_regulons')
	aucell_regulons <- aucell_regulons[onlyNonDuplicatedExtended(names(aucell_regulons))]

	genie3ll <- loadInt(scenicOptions, 'genie3ll')
	genie3ll$TF <- as.character(genie3ll$TF)
	genie3ll$Target <- as.character(genie3ll$Target)
	tf2genes.coex <- split(genie3ll$Target, genie3ll$TF)
	
	TF2targets <-NULL
	tfs <- names(aucell_regulons)
	for(i in seq(tfs)){
		tf <- getTF(tfs[i])
		tf <- sub("_extended","", tf)
		targets.coex <- paste(tf2genes.coex[[tf]], collapse=";")
		targets.regulon <- paste(setdiff(aucell_regulons[[tfs[i]]], tf), collapse=";")
		TF2targets <- rbind(TF2targets, c( tfs[i], targets.coex, targets.regulon))
	}
	TF2targets <- as.data.frame(TF2targets)
	rownames(TF2targets) <- TF2targets$V1
	TF2targets <- TF2targets[,c(2,3)]
	colnames(TF2targets) <- c("tf_gene_coexpressed_module", "tf_targets")
	write.table(cbind(regulons=rownames(TF2targets), TF2targets), sep="\t", col.names=T, row.names=F, quote=F, file="Regulon_targets.xls")
	regulons_num <- length(aucell_regulons)
	TFs_num <- length(unique(sub("_extended","",  getTF(names(aucell_regulons)))))
	targets <- NULL
	for(i in names(aucell_regulons)){
		TF <- sub("_extended","",  getTF(i))
		tmp <- setdiff(aucell_regulons[[i]], TF)
		targets <- c(targets, tmp)
	}
	targets <- unique(targets)
	targets_num <- length(targets)
	stat <- data.frame(TF=TFs_num,target_gene=targets_num,regulons=regulons_num)
	write.table(stat,file="Regulon_stat.xls", sep="\t", row.names=F,col.names=T,quote=F)
}

NetPlot <- function(edge, node, outfile="TF-target.network.pdf"){
	node$shape <- "circle"
	node$color[node$annot == "tf"] <- "#781028"
	node$color[node$annot != "tf"] <- "#B6CECE"
	g <- graph_from_data_frame(edge, directed = F, vertices = node)
	layout <- layout_with_fr(g)
	plotIgraph <- function(g, layout, output) {
		if(grepl('.png$', output)) {
			png(file = output, width = 8, height = 6,units = "in", res = 300)
		} else if (grepl('.pdf$', output)) {
			pdf(file = output, width = 8, height = 6)
		} else {
			stop('Output not support!')
		}
		plot(g, layout = layout,vertex.label.cex = 1,vertex.shape = node$shape)
		legend( "topright",c("TF","Target"), 
				       pch = c(16, 16), col = c('#781028', '#B6CECE'), 
					          pt.cex = 2, cex = .8, bty = "n", ncol = 1)
		dev.off()
	}
	plotIgraph(g, layout, outfile)

}

DoNetwork <- function(scenicOptions, ntop = 5){
	aucell_regulons <- loadInt(scenicOptions, 'aucell_regulons')
	aucell_regulons <- aucell_regulons[onlyNonDuplicatedExtended(names(aucell_regulons))]

	for(i in names(aucell_regulons)){
		tf <- getTF(i)
		tf <- sub("_extended","", tf)
		aucell_regulons[[i]] <- setdiff(aucell_regulons[[i]], tf)
	}
	aucell_regulons <- reshape2::melt(aucell_regulons)
	aucell_regulons$L1 <- sub("_extended","", getTF(aucell_regulons$L1))
	aucell_regulons <- cbind(aucell_regulons, tf=aucell_regulons$L1, target= aucell_regulons$value)
	aucell_regulons <- aucell_regulons[,c("tf", "target")]

	edge <- aucell_regulons
	write.table(edge, file="tf_target.egde.tsv", sep="\t", col.names=T, row.names=F, quote=F)
	node <- data.frame(node=union(edge$tf, edge$target))
	node$annot <- ifelse(node$node %in% edge$tf, "tf", "target")
	write.table(node, file="tf_target.node.tsv", sep="\t", col.names=T, row.names=F, quote=F)
	tf_order_bytargets <- names(sort(table(edge$tf), decreasing=F))
	tfs <- tf_order_bytargets[1:ntop]
	edge <- edge[edge$tf %in% tfs,]
	node <- data.frame(node=union(edge$tf, edge$target))
	node$annot <- ifelse(node$node %in% edge$tf, "tf", "target")
	NetPlot(edge, node)

}


.TSNEplot <- function(
                scenicOptions,
                cells.use = NULL,
                group.by = "cluster",
                group.use = NULL,
                color.use = NULL,
                xlab = NULL,
                ylab = NULL,
                reduction = "tsne",
                pt.size=1,
                label.size=6,
                do.label=F,
                legend.label = NULL,
                legend.psize = 4,
                whitetheme = T,
                do.return = T
){
        meta.data <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")
        if(! group.by %in% colnames(meta.data)){
             return(0)
        }
        projection <- NULL
        reduction_key <- NULL
		reduction_dim <- scenicOptions@settings$defaultTsne$dims
		perpl <- scenicOptions@settings$defaultTsne$perpl
        if(reduction == "tsne"){
            projection <- readRDS(paste0(scenicOptions@fileNames$dir,"/4.AUC/tSNE_AUC_", reduction_dim, "pcs_", perpl,"perpl.Rds"))
            reduction_key <- "tSNE"
			projection <- projection$Y
        }else{
            projection <- readRDS(paste0(scenicOptions@fileNames$dir,"/4.AUC/UMAP_AUC_", reduction_dim, "pcs_", perpl,"perpl.Rds")) ## not supported temporary
            reduction_key <- "UMAP"
        }
        proj1.min <- min(projection[,1])
        proj2.max <- max(projection[,2])
		proj2.min <- min(projection[,2]) - 5
        colnames(projection) <- paste0(reduction_key, c("_1","_2"))
        meta.data <- cbind(meta.data, projection[rownames(meta.data),])
        colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists="null")
        
        dim1 <- paste0(reduction_key,"_1")
        dim2 <- paste0(reduction_key,"_2")
        meta.data$x <- meta.data[[dim1]]
        meta.data$y <- meta.data[[dim2]]
        meta.data$group <- meta.data[[group.by]]
        color.use <- colVars[[group.by]]
        if(!is.null(group.use)){
                meta.data <- meta.data[ meta.data[[group.by]] %in% group.use, ]
        }
        if(!is.null(cells.use)){
                meta.data <- meta.data[cells.use, ]
        }
        meta.data <- droplevels(meta.data)
        p <- ggplot(data = meta.data, mapping = aes(x = x, y = y))
        p <- p + geom_point(mapping = aes(colour = group), size = pt.size) + scale_y_continuous(limit=c(proj2.min, proj2.max))
        if(do.label){
                centers <- meta.data %>% dplyr::group_by(group) %>% summarize(x = median(x = x),y = median(x = y))
                p <- p + geom_point(data = centers, mapping = aes(x = x,y = y), size = 0, alpha = 0) +  geom_text(data = centers, mapping = aes(label = group), size = label.size)
        }
        
        if(is.null(legend.label)){
                p <- p +scale_colour_manual("",values=color.use)
        }else{
                p <- p + scale_colour_manual("",values=color.use,labels= as.vector(legend.label)) 
        }
		p <- p + labs(x=xlab,y=ylab)
#        p <- p + guides(color = guide_legend(override.aes = list(size = 4))) +theme(axis.line.x=element_line(size=0.6), axis.line.y=element_line(size=0.6), axis.ticks=element_line(size=0.6)) + labs(x=xlab,y=ylab)
        
        if(whitetheme){
            p <- p + theme_bw() + theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), panel.grid=element_blank(), rect=element_rect(size=1))
        }
		legend.row <- ceiling(length(color.use)/10 + 0.5)

		p <- p + theme(legend.position=c(0.5,0.06), legend.direction="horizontal") + guides(color = guide_legend(override.aes = list(size = 4), nrow=legend.row))
        if(do.return){
            return(p)
        }
}

TSNEplot <- function(scenicOptions, width=7, height = 7, reduction="tsne"){

	meta.data <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")
	
	outprf <- if(reduction == "tsne"){ "tSNE"}else{"UMAP"}
	if("group" %in% colnames(meta.data)){
		p<- .TSNEplot(scenicOptions, group.by="group", reduction=reduction)
		ggsave(p, file=paste0("group.", outprf, ".pdf"), width=width, height=height)

		for (i in unique(meta.data$group)){
			p <- .TSNEplot(scenicOptions, group.by = "group", group.use=i, reduction=reduction)
			ggsave(p, file=paste0(outprf,"_", i,".pdf"), width=width, height=height)
		}
	}

	if("cluster" %in% colnames(meta.data)){
		p <- .TSNEplot(scenicOptions, group.by="cluster", reduction=reduction)
		ggsave(p, file=paste0("cluster.", outprf, ".pdf"), width=width, height=height)

		for (i in unique(meta.data$cluster)){
			p <- .TSNEplot(scenicOptions, group.by = "cluster", group.use=i, reduction=reduction)
			ggsave(p, file=paste0(outprf,"_cluster", i,".pdf"), width=width, height=height)
		}
	}

	if("sample" %in% colnames(meta.data)){
		p <- .TSNEplot(scenicOptions, group.by="sample", reduction=reduction)
		ggsave(p, file=paste0("sample.", outprf, ".pdf"), width=width, height=height)

		for (i in unique(meta.data$sample)){
			p <- .TSNEplot(scenicOptions, group.by = "sample", group.use=i, reduction=reduction)
			ggsave(p, file=paste0(outprf,"_", i,".pdf"), width=width, height=height)
		}
	}
}

singleFeaturePlots <- function(
	data.use,
	feature,
	data.plot,
	pt.size,
	pch.use,
	cols.use,
	dim.codes,
	min.cutoff,
	max.cutoff,
	coord.fixed,
	no.axes,
	no.title = FALSE,
	no.legend,
	theme,
	features.blend,
	high.up= TRUE,
	title = NULL,
	breaks=NULL
){
	if(features.blend){
		data.gene <- as.data.frame(apply(data.use[feature,],2,function(x)mean(x)))
		colnames(data.gene) <- "gene"
		data.gene <- na.omit(object =data.gene)
		min.cutoff <- min(data.gene)
		max.cutoff <- max(data.gene)
	}else{
		data.gene <- na.omit(object = data.frame(data.use[feature,]))
		min.cutoff <-Seurat:::SetQuantile(cutoff = min.cutoff, data = data.gene)
		max.cutoff <- Seurat:::SetQuantile(cutoff = max.cutoff, data = data.gene)
	}
	data.gene <- sapply(
		X = data.gene,
		FUN = function(x){
			return(ifelse(test = x < min.cutoff, yes = min.cutoff,no = x))
		}
	)
	
	data.gene <- sapply(
		X = data.gene,
		FUN = function(x){
		return(ifelse(test = x > max.cutoff, yes = max.cutoff,no = x))
		}
	)

	data.plot$gene <- data.gene
	data.plot$cell.name <- rownames(data.plot)

	if(high.up){
		data.plot <- arrange(data.plot, gene)
		rownames(data.plot) <- data.plot$cell.name
	}
	
	p <- ggplot(data = data.plot, mapping = aes(x = x, y = y))
	
	if (all(data.plot$gene == data.plot$gene[1])) {
		warning(paste0("All cells have the same value of ",feature, "."))
		p <- p + geom_point(color = cols.use[1], size = pt.size,shape = pch.use)
	}else {
		p <- p + geom_point(mapping = aes(color = gene),size = pt.size, shape = pch.use)
		if(length(cols.use) == 2){
			p <- p + scale_color_gradientn(colors = cols.use, limits=breaks, guide = guide_colorbar(title = ""))
		}else if (length(cols.use) == 3){
			p <- p + scale_color_gradientn(limits=breaks[c(1,3)], guide = guide_colorbar(title = ""), colours = cols.use)
		}else{
			stop("cols.use muse be two or three colors")
		}
	}
	if (theme){
#		p <- p + theme_classic()
		p <- p + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title = element_blank())
	}
	
	if(features.blend){
		p <- p + labs(title = title, x = dim.codes[1],y = dim.codes[2])
	}else{
		if(!is.null(title)){
			p <- p + labs(title = title, x = dim.codes[1],y = dim.codes[2])
		}else{
			p <- p + labs(title = feature, x = dim.codes[1],y = dim.codes[2])
		}
	}
	p <- p +  theme(plot.title = element_text(hjust = 0.5))
	
	if (no.legend) {
		p <- p + theme(legend.position = "none")
	}
	if (coord.fixed) {
		p <- p + coord_fixed()
	}
	p <- p + theme(rect=element_rect(size=1))
	return(p)
}



.FeaturePlots <- function(
    scenicOptions,
    features.plot=NULL,
    features.type = c("AUC","expression", "BinaryAUC"),
    min.cutoff = NA,
    max.cutoff = NA,
    dim.1 = 1,
    dim.2 = 2,
    cells.use = NULL,
    pt.size = 3,
    cols.use= c("lightgrey", "blue"),
    pch.use = 20,
    reduction.use = "tsne",
    use.imputed = FALSE,
    nCol = NULL,
    no.axes = FALSE,
    no.legend = FALSE,
    coord.fixed = FALSE,
    theme= TRUE,
    features.blend = FALSE,
    high.up = TRUE,
    title=NULL,
    breaks=NULL,
    use.seurat.dimension=F,
    do.return =F
){
	setwd(paste0(scenicOptions@fileNames$dir, "/4.AUC/Step3_RegulonActivity_tSNE_colByActivity"))
	projection <- NULL
	reduction_key <- NULL
	meta.data <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")
	regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
	AUC <- getAUC(regulonAUC)
	reduction_dim <- scenicOptions@settings$defaultTsne$dims
	perpl <- scenicOptions@settings$defaultTsne$perpl
	if(reduction.use == "tsne"){
		projection <- readRDS(paste0(scenicOptions@fileNames$dir,"/4.AUC/tSNE_AUC_", reduction_dim, "pcs_", perpl,"perpl.Rds"))
		reduction_key <- "tSNE"
		projection <- projection$Y
	}else{
		projection <- readRDS(paste0(scenicOptions@fileNames$dir,"/4.AUC/UMAP_AUC_", reduction_dim, "pcs_", perpl,"perpl.Rds")) ## not supported temporary
		reduction_key <- "UMAP"
	}
	colnames(projection) <- paste0(reduction_key, c("_1","_2"))
	if(use.seurat.dimension){
		obj <- readRDS(paste0(scenicOptions@fileNames$dir,"/1.CreateScenicOptions/obj.Rds"))
		projection <- Embeddings(obj,reduction=reduction.use)
	}
	meta.data <- cbind(meta.data, projection[rownames(meta.data),])
	data.use <- NULL
	if (features.type == "expression"){
		data.use <- readRDS(paste0(scenicOptions@fileNames$dir,"/2.FilterGene/exprMat.Rds"))
		data.use <- log2(data.use + 1)
	}else if(features.type == "AUC"){
		data.use <- AUC
	}else{
		data.use <- loadInt(scenicOptions, "aucell_binary_full")
	}
#    cells.use <- Seurat:::SetIfNull(x = cells.use, default = colnames(x = data.use))
    cells.use <- if(is.null(cells.use)){colnames(x = data.use)}

    if (is.null(x = nCol)) {
        nCol <- 2
        if (length(x = features.plot) == 1) {
            nCol <- 1
        }
        if (length(x = features.plot) > 6) {
            nCol <- 3
        }
        if (length(x = features.plot) > 9) {
            nCol <- 4
        }
    }
    num.row <- floor(x = length(x = features.plot)/nCol - 1e-05) + 1
#    dim.code <- GetDimReduction(object = object, reduction.type = reduction.use, slot = "key")
    dim.code <- paste0(reduction_key, "_")
    dim.codes <- paste0(dim.code, c(dim.1, dim.2))
#    data.plot <- as.data.frame(GetCellEmbeddings(object = object, reduction.type = reduction.use, dims.use = c(dim.1, dim.2), cells.use = cells.use))
    data.plot <- as.data.frame(meta.data[cells.use, dim.codes])
    x1 <- paste0(dim.code, dim.1)
    x2 <- paste0(dim.code, dim.2)
    data.plot$x <- data.plot[, x1]
    data.plot$y <- data.plot[, x2]
    data.plot$pt.size <- pt.size
#    names(x = data.plot) <- c("x", "y")
#    data.use <- t(x = FetchData(object = object, vars.all = features.plot, cells.use = cells.use, use.imputed = use.imputed))
    min.cutoff <- mapply(FUN = function(cutoff, feature) {
        ifelse(test = is.na(x = cutoff), yes = min(data.use[feature, 
            ]), no = cutoff)
    }, cutoff = min.cutoff, feature = features.plot)
    max.cutoff <- mapply(FUN = function(cutoff, feature) {
        ifelse(test = is.na(x = cutoff), yes = max(data.use[feature, 
            ]), no = cutoff)
    }, cutoff = max.cutoff, feature = features.plot)
    check_lengths = unique(x = vapply(X = list(features.plot, 
        min.cutoff, max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
    if (length(x = check_lengths) != 1) {
        stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
		p <- singleFeaturePlots(
				feature = features.plot,
				min.cutoff = min.cutoff,
				max.cutoff = max.cutoff,
				coord.fixed = coord.fixed,
				data.use = data.use,
				data.plot = data.plot,
				pt.size = pt.size,
				pch.use = pch.use,
				cols.use = cols.use,
				dim.codes = dim.codes,
				no.axes = no.axes,
				no.legend = no.legend,
				theme = theme,
				high.up=high.up,
				features.blend = features.blend,
				title=title,
				breaks=breaks
				)
	if(features.type == "AUC" | features.type == "BinaryAUC") features.plot <- gsub(" ","", features.plot)
	if(do.return){
		return(p)
	}
	ggsave(p, file=paste0(features.plot, "_", features.type, ".pdf"))
	setwd(paste0(scenicOptions@fileNames$dir, "/4.AUC"))
}


FeaturePlots <- function(
    scenicOptions,
    features.type = c("AUC","expression", "BinaryAUC"),
    cols.use= c("lightgrey", "blue"),
    reduction.use = "tsne"
){
	regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
	AUC <- getAUC(regulonAUC)
	regulons <- rownames(AUC)
	regulons.NonDuplicatedExtended <- onlyNonDuplicatedExtended(regulons)
	features.plot <- NULL
	if (features.type == "expression"){
		features.plot <- names(regulons.NonDuplicatedExtended)
	}else if (features.type == "AUC"){
		features.plot <- regulons
	}else{
		BinaryAUC <- loadInt(scenicOptions, 'aucell_binary_full')
		features.plot <- rownames(BinaryAUC)
	}
	for( i in features.plot){
		.FeaturePlots(scenicOptions, features.plot=i, features.type=features.type, cols.use=cols.use, do.return=F, reduction.use=reduction.use, title=i)
	}
}




nodeColor <- function (sigmaObj, oneColor = "#E41A1C", colorAttr = NULL, colorPal = NULL)
{
  edges <- jsonlite::fromJSON(sigmaObj$x$data)$edges
  nodes <- jsonlite::fromJSON(sigmaObj$x$data)$nodes
  if (!is.null(colorPal)) {
    nodes$tempCol <- sigmaObj$x$graph$vertices[, colorAttr]
    name <- levels(as.factor(nodes$tempCol))
    names(colorPal) <- name
    pal <- colorPal[1:length(nodes$tempCol)]
    nodes$color <- pal[nodes$tempCol]
    nodes$tempCol <- NULL
  }
  else {
    nodes$color <- oneColor
  }
  graphOut <- list(nodes, edges)
  names(graphOut) <- c("nodes", "edges")
  sigmaObj$x$data <- jsonlite::toJSON(graphOut, pretty = TRUE)
  return(sigmaObj)
}

edgeColor <- function (sigmaObj, oneColor = "#E41A1C", colorAttr = NULL, colorPal = NULL) 
{
  edges <- jsonlite::fromJSON(sigmaObj$x$data)$edges
  nodes <- jsonlite::fromJSON(sigmaObj$x$data)$nodes
  if (!is.null(colorPal)) {
    edges$tempCol <- sigmaObj$x$graph$edges[, colorAttr]
    name <- levels(as.factor(edges$tempCol))
    names(colorPal) <- name
    pal <- colorPal[1:length(edges$tempCol)]
    edges$color <- pal[edges$tempCol]
    edges$tempCol <- NULL
  }
  else {
    edges$color <- oneColor
  }
  graphOut <- list(nodes, edges)
  names(graphOut) <- c("nodes", "edges")
  sigmaObj$x$data <- jsonlite::toJSON(graphOut, pretty = TRUE)
  return(sigmaObj)
}

GeneSplit <- function(features, splitnum=10){
	num <- floor(length(features)/splitnum)
	left <- length(features) %% splitnum
	for(i in seq(num)){
		index1 <- (i-1) * splitnum + 1
		index2 <- i * splitnum
		genes <- features[index1:index2]
		write.table(genes, sep="\t", quote=F, col.names=F, row.names=F, file= paste0( "genesplit-", i))
		}
	if(left){
		index1 <- num * splitnum + 1
		index2 <- length(features)
		genes <- features[index1:index2]
		write.table(genes, sep="\t", quote=F, col.names=F, row.names=F, file= paste0( "genesplit-", i))
	}
}
