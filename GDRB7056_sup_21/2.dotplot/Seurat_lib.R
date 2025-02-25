library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
		 
### subroutine
MakeSeuratObj <- function(parameter){
		object.list <- list()
		for ( i in seq(parameter$data$name) ){
				mat <- Read10X(data.dir = parameter$data$dir[i], gene.column = 1)
				object.list[[i]] <- CreateSeuratObject(counts = mat, project = parameter$data$name[i], assay = "RNA" )
		}
		if ( length(object.list) == 1 ) {
				object <- object.list[[1]]
		}else{
				object <- merge(x = object.list[[1]], y = unlist(object.list[-1]), add.cell.ids = parameter$data$name)
		}

		object@meta.data$orig.ident <- factor(object@meta.data$orig.ident, levels = parameter$data$name)
		object <- SetIdent(object, value = "orig.ident")

		object@misc[["fdata"]] <- AddFData(object, parameter$name_list)
		object@misc[["pdata"]] <- FetchData(object, "orig.ident")
		object@misc[["counts"]] <- GetAssayData(object, slot = "counts", assay = "RNA")

		color.sample <- scales::hue_pal()(length(levels(object@meta.data$orig.ident)))
		names(color.sample) <- levels(object@meta.data$orig.ident)
		object@misc[["color.sample"]] <- color.sample

		return(object)
}


AddFData <- function(object, ref_name_file = NULL, col.name = NULL){
		if ( ! is.null(ref_name_file) && file.exists(ref_name_file) ) {
				fdata <- read.table(ref_name_file, row.names = 1, stringsAsFactors = F)
				colnames(fdata) <- c("merge_name", "name", "type")
		}else{
				fdata <- data.frame(name = rownames(object), row.names = rownames(object), stringsAsFactors = F)
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


AddUnderscore <- function(data){
		if ( ! exists("underscore", data) || ! exists("dash", data) ) {
				data$underscore <- rownames(data)
				data$dash <- gsub("_", "-", rownames(data))
		}
		return(data)
}


RenameFeatures <- function(object, new.names = NULL, form.type = c("id", "merge_name"), to.type = c("id", "merge_name")) {
		if ( is.null(new.names) ) {
				form.type <- match.arg(form.type)
				to.type <- match.arg(to.type)
				fdata <- object@misc$fdata
				fdata$id <- rownames(fdata)
				new.names <- fdata[,to.type]
				names(new.names) <- fdata[,form.type]
#				new.names <- data.frame(to.type = fdata[,to.type], row.names = fdata[,form.type])
		}
		assays <- Seurat:::FilterObjects(object = object, classes.keep = "Assay")
		for (assay in assays) {
				slot(object = object, name = "assays")[[assay]] <- RenameFeatures.Assays(object = object[[assay]], new.names = new.names)
		}
		dimreducs <- Seurat:::FilterObjects(object = object, classes.keep = "DimReduc")
		for (dr in dimreducs) {
				object[[dr]] <- RenameFeatures.DimReduc(object = object[[dr]], new.names = new.names)
		}
		return(object)
}

RenameFeatures.Assays <- function(object, new.names = NULL ) {
		for (data.slot in c("counts", "data", "scale.data")) {
				old.data <- GetAssayData(object = object, slot = data.slot)
				if (nrow(x = old.data) <= 1) {
						next
				}
				old.name <- rownames(x = slot(object = object, name = data.slot))
				rownames(x = slot(object = object, name = data.slot)) <- new.names[old.name]
		}
		if ( length(slot(object = object, name = "var.features")) > 0 ) {
				old.name <- rownames(x = slot(object = object, name = "var.features"))
				slot(object = object, name = "var.features") <- new.names[old.name]
		}
		return(object)
}

RenameFeatures.DimReduc <- function(object, new.names = NULL ) {
		for ( projected in c(TRUE, FALSE) ){
				data.slot <- ifelse(projected, "feature.loadings.projected", "feature.loadings")
				old.data <- Loadings(object = object, projected = projected)
				rownames(x = old.data) <- new.names[rownames(x = old.data)]
				slot(object = object, name = data.slot) <- old.data
		}
		return(object)
}

StatFeatures <- function(object, features = NULL, col.name = NULL, stat_pct = FALSE, assay = NULL){
		if ( file.exists(features[1]) )
				features <- readLines(con = features[1])
		object@misc[[col.name]] <- FindFeaturesID(object = object, features = features, unlist = FALSE)
		features <- unlist(object@misc[[col.name]])
		if ( is.null(assay) )
				assay <- DefaultAssay(object = object)
		metadata <- Matrix::colSums(x = GetAssayData(object = object, slot = "counts", assay = assay)[features, , drop = FALSE])
		if (stat_pct)
				metadata <- metadata / object[[paste0("nCount_", assay)]] * 100
		if (!is.null(x = col.name)) {
				object <- AddMetaData(object = object, metadata = metadata, col.name = col.name)
				return(object)
		} else {
				return(metadata)
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
							message("[WARNING] '", x, "' not found gene id.")
							return(NULL)
						})
		if ( unlist ) features <- unlist(features)
		return(features)
}

FindFeaturesName <- function(object, features, col = "merge_name", is.fast = FALSE) {
		if ( ! exists("fdata", object@misc) ) return(features)
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
								ifelse(is.null(id), x, id)
								})
		}
		return(Name)
}

VlnplotShell <- function(object, features = NULL, labs.y = NULL, titles = NULL, outfile = NULL, group.by = "orig.ident", cols.use = NULL, do.split.save = FALSE, pref = NULL, ...){
		if ( is.null(features) ) return(1)
		nCol <- ifelse(length(features) == 3, 3, ceiling(sqrt(length(features))))
		nRow <- ceiling(length(features) / nCol)
		if ( "pt.size" %in% formalArgs(VlnPlot) ) {
		  plots <- VlnPlot(object = object, features = features, ncol = nCol, group.by = group.by, cols = cols.use, pt.size = 0.1, combine = FALSE, ...)
		}else{
		  plots <- VlnPlot(object = object, features = features, nCol = nCol, group.by = group.by, cols = cols.use, point.size.use = 0.1, combine = FALSE, ...)
		}
		for ( i in seq(features) ){
				if ( group.by == "orig.ident" )
						plots[[i]] <- plots[[i]] + xlab("Samples")
				if ( !is.null(labs.y[i]) && !is.na(labs.y[i]) )
						plots[[i]] <- plots[[i]] + ylab(labs.y[i])
				if ( !is.null(titles[i]) && !is.na(titles[i]) )
						plots[[i]] <- plots[[i]] + ggtitle(titles[i])
		}

		if ( do.split.save && !is.null(pref) ) {
				for ( i in seq(features) ){
						name <- ifelse( !is.null(names(features[i])), names(features[i]), features[i])
						outfile.tmp <- paste(pref, name, "pdf", sep = "." )
						ggsave(plots[[i]] + theme(legend.position = "none"), file = outfile.tmp, width = 6, height = 6)
				}
		}

		combine.args <- list(plots = plots, ncol = nCol)
		combine.args <- c(combine.args, list(...))
		if (!"legend" %in% names(x = combine.args)) {
				combine.args$legend <- "none"
		}
		plots.combined <- do.call(what = "CombinePlots", args = combine.args)
		if( !is.null(dev.list()) ) dev.off()
		
		if ( is.null(outfile) ){
				return(plots.combined)
		}else{
				ggsave(plots.combined, file = outfile, width = 6 * nCol, height = 6 * nRow, limitsize = F)
		}
}

StatPlot <- function(object, pref = NULL, color.sample = NULL, assay = "RNA" ){
		if ( is.null(color.sample) )
				color.sample <- object@misc[["color.sample"]]
		if ( is.null(assay) ) assay <- DefaultAssay(object)
		nFeature <- paste0("nFeature_", assay)
		nCount <- paste0("nCount_", assay)
		d <- data.frame("nGene" = c(nFeature, "Number"),
						"nUMI"  = c(nCount,   "Number"),
						"pMito" = c("percent.mito", "Percentage(%)"),
						stringsAsFactors = FALSE)
		for ( i in names(d) ) {
				if ( ! exists(d[[i]][1], object@meta.data) ) d[[i]] <- NULL
		}
		VlnplotShell(object, features = unlist(d[1,]), labs.y = unlist(d[2,]), titles = colnames(d), outfile = paste0(pref, ".merge.pdf"), cols.use = color.sample, do.split.save = TRUE, pref = pref )

		d <- data.frame(expected.marker = c("expected.marker", "Number"),
						exclude.marker = c("exclude.marker", "Number"),
						stringsAsFactors = FALSE)
		for ( i in names(d) ) {
				if ( ! exists(d[[i]][1], object@meta.data) ) d[[i]] <- NULL
		}
		VlnplotShell(object, features = unlist(d[1,]), labs.y = unlist(d[2,]), titles = colnames(d), outfile = paste0(pref, ".PresetMarker.pdf"), cols.use = color.sample)

		FeatureScatterShell(object, feature1 = nCount, feature2 = nFeature, outfile = paste0(pref, ".nUMI-nGene.pdf"), cols = color.sample)
		if ( exists("percent.mito", object@meta.data) )
				FeatureScatterShell(object, feature1 = nCount, feature2 = "percent.mito", outfile = paste0(pref, ".nUMI-pMito.pdf"), cols = color.sample)
}

FeatureScatterShell <- function(object, feature1 = NULL, feature2 = NULL, outfile = NULL, group.by = "orig.ident", cols = NULL){
		p <- FeatureScatter(object = object, feature1 = feature1, feature2 = feature2, span = 0.8, group.by = group.by, pt.size = 0.1, cols = cols) +
				guides(color = guide_legend(override.aes = list(size = 2), title = NULL))
		ggsave(p, file = outfile, width = 6, height = 6)
}

FilterGenes <- function(object, parameter) {
		message( "--->Filter Genes<---" ) 
		if ( exists( "min.cell", parameter$filter ) && is.numeric(parameter$filter$min.cell) && parameter$filter$min.cell > 0 ){
				if ( parameter$filter$min.cell < 1 )
						parameter$filter$min.cell <- parameter$filter$min.cell * length(object@cell.names)
				num.cells  <- Matrix::rowSums( object@assays$RNA@counts > 0 )
			
				genes.filter <- num.cells[which( num.cells < parameter$filter$min.cell )]
				genes.filter <- FindFeaturesName(object, genes.filter)
				write.table( genes.filter, file = "filtered_genes.xls", quote = F, sep = "\t", col.name = F )
				
				genes.use  <- names( num.cells[which( num.cells >= parameter$filter$min.cell )] )
				object     <- object[genes.use, ]

				## TODO : need to re-calculate 'percent.mito', 'expected.marker', ... value ?
		}
		return(object)
}


FilterCells <- function(object, parameter, do.stat = TRUE) {
		message( "--->Filter Cells<---" )
		cells.use <- Cells(object)
		if ( !is.null(parameter$filter$set.num) && parameter$filter$set.num != "none" ){
				if ( parameter$filter$set.num == "min" ){
						cell_num <- min(table(object@meta.data$orig.ident))
				} else if ( parameter$filter$set.num != "none" ) {
						cell_num <- min( as.integer(parameter$filter$set.num), max(table(object@meta.data$orig.ident)) )
				}
				seed <- ifelse(!is.null(parameter$filter$set.num.seed), parameter$filter$set.num.seed, 42)
				set.seed(seed)
				cells.use <- as.character(unlist(by(Cells(object), object@meta.data$orig.ident, function(x) sample(x, min(cell_num, length(x))))))
		}
		for ( i in names(parameter$filter$standard) ) {
				if ( exists(i, object@meta.data) ) {
				if ( length(parameter$filter$standard[[i]]) == 1 ) {
						value <- parameter$filter$standard[[i]]
						cells.use <- object@meta.data %>% tibble::rownames_to_column(var = "cells") %>%
								filter(.data[[i]] == value & cells %in% cells.use) %>%
								select(cells) %>% unlist()
				} else if ( length(parameter$filter$standard[[i]]) == 2 ) {
						lower <- parameter$filter$standard[[i]][1]
						upper <- parameter$filter$standard[[i]][2]
						cells.use <- object@meta.data %>% tibble::rownames_to_column(var = "cells") %>%
								filter(.data[[i]] >= lower & .data[[i]] <= upper & cells %in% cells.use) %>%
								select(cells) %>% unlist()
				}
				}
		}		
		object <- object[, cells.use]

		if ( do.stat ) {
				StatFilterCells(object)
		}

		return(object)
}


StatFilterCells <- function(object){
		filter_cells <- setdiff(rownames(object@misc$pdata), colnames(object))
		write.table( filter_cells, file = "filtered_cells.xls", quote = F, sep = "\t", col.name = F )

		filter_stat_table     <- as.data.frame( cbind( before_filter_num = table(object@misc$pdata$orig.ident), after_filter_num = table(object@meta.data$orig.ident) ) )
		filter_stat_table$pct <- paste( round( filter_stat_table$after_filter_num / filter_stat_table$before_filter_num * 100, 2) , "%", sep = '' )
		filter_stat_table     <- data.frame( Sample = row.names(filter_stat_table), filter_stat_table )
		filter_stat_table     <- object@meta.data %>% group_by(Sample = orig.ident) %>%
								 summarise(after_filter_median_UMI_per_cell = median(nCount_RNA), after_filter_median_genes_per_cell = median(nFeature_RNA)) %>%
								 left_join(x = filter_stat_table, by = "Sample")
		write.table( filter_stat_table, file = "Filter.stat.xls", quote = F, sep = "\t", row.names = F  )
}

RestoreObject <- function(object) {
		new.object <- CreateSeuratObject(counts = object@misc[["counts"]], meta.data = object@misc[["pdata"]], project = object@project.name, assay = "RNA" )
		new.object@misc <- object@misc
		return(new.object)
}

CheckVariableFeature <- function(object){
		top10 <- head(VariableFeatures(object), 10)
		plot1 <- VariableFeaturePlot(object)
		plot2 <- LabelPoints(plot = plot1, points = top10, labels = object@misc$fdata[top10,"merge_name"], repel = TRUE, xnudge = 0, ynudge = 0)
		plot2 <- plot2 + theme(legend.position = "top")
		ggsave(plot2, file = "Variable_gene.pdf", width = 6, height = 6)

		write.table(VariableFeatures(object), file = "var_gene.xls", quote = F, sep = "\t", col.name = F)
}


FindRegressVars <- function(object, parameter, vars.regress = NULL, force_recal = FALSE, ...){
		if ( is.null(vars.regress) ) {
#				vars.regress <- "nUMI"
				vars.regress <- paste0("nCount_", DefaultAssay(object)) # nCount_RNA
				if ( exists("percent.mito", object@meta.data) ) vars.regress <- c(vars.regress, "percent.mito")
				if ( exists("cell_cycle", parameter) && parameter$cell_cycle$is_remove ){
						message("-->CellCycle Scoring<--")
						if ( ! exists("CC.Difference", object@meta.data) || force_recal ) {
								object <- DoCellCycleScoring(object, ...)
						}
						if ( parameter$cell_cycle$is_rm_all_signal ) {
								vars.regress <- c(vars.regress, "S.Score", "G2M.Score")
						} else {
								vars.regress <- c(vars.regress, "CC.Difference")
						}
				}
		}else if ( vars.regress == "none" ){
				vars.regress <- NULL
		}
		object@misc$vars.regress <- vars.regress
		return(object)
}

DoCellCycleScoring <- function(object, ...){
		s.genes   <- FindFeaturesID(object, cc.genes$s.genes)
		g2m.genes <- FindFeaturesID(object, cc.genes$g2m.genes)
		if ( length(s.genes) < 2 || length(g2m.genes) < 2 ) return(object)
		object <- CellCycleScoring(object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, ...)
		object[["CC.Difference"]] <- object[["S.Score"]] - object[["G2M.Score"]]
		return(object)
}

DoNormalization <- function(object, parameter, assay = "RNA", is_SCTransform = FALSE, vars.regress = NULL, ...){
		DefaultAssay(object) <- assay
		### pre-deal with Cell Cycle
		object <- FindRegressVars(object, parameter, vars.regress = vars.regress, ...)

		if ( is_SCTransform ) {
				### Another workflow : SCTransform 
				object <- SCTransform(object, assay = assay, vars.to.regress = object@misc$vars.regress,
				                      verbose = FALSE, min_cells = 1, return.only.var.genes = F)
		}else{
				### Normalize Data
				message( "-->Normalize Data<--" )
				object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
				### Find Variable Features
				message( "-->Find Variable Genes<--" )
				object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
				CheckVariableFeature(object)
				### Confound factors and Scale Data
				message( "-->ScaleData<--" )
				object <- ScaleData(object = object, vars.to.regress = object@misc$vars.regress, features = rownames(object))
		}
		return(object)
}

DoDimReduc <- function(object, assay = NULL, pc.num = 50, ...) {
		if ( is.null(assay) ) assay <- DefaultAssay(object)				
		message( "-->PCA<--" )
		object <- RunPCA(object = object, assay = assay, npcs = pc.num, features = VariableFeatures(object), verbose = FALSE )
		object[[paste0("pca_", assay)]] <- object[["pca"]]
		CheckPCA(object)

		sig.PCs <- seq(pc.num)
		message( "-->Run tSNE<--" )
		object <- RunTSNE(object, dims = sig.PCs, ...)
		object[[paste0("tsne_", assay)]] <- object[["tsne"]]

		message( "-->Run UMAP<--" )
#		object <- RunUMAP(object, dims = sig.PCs, umap.method = "umap-learn")
		object <- RunUMAP(object, dims = sig.PCs, umap.method = "uwot")
		object[[paste0("umap_", assay)]] <- object[["umap"]]

		return(object)
}

DoIntegration <- function(object, split.by = "orig.ident", dims = 1:50, nfeatures = 3000, is.SCT = FALSE){
		old.assay <- DefaultAssay(object)
		object.list <- SplitObject(object, split.by = split.by)
		anchor.features <- nfeatures
		normalization.method <- "LogNormalize"
		if ( is.SCT ) {
#				for (i in seq(object.list)) {
#						object.list[[i]] <- SCTransform(object.list[[i]], vars.to.regress = object@misc$vars.regress, verbose = FALSE)
#				}
				anchor.features <- SelectIntegrationFeatures(object.list = object.list, nfeatures = nfeatures)
				object.list <- PrepSCTIntegration(object.list = object.list, anchor.features = anchor.features)
				normalization.method <- "SCT"
		}
#		k.filter <- min(200, min(sapply(object.list, ncol)))
		k.filter <- min(200, ceiling(min(sapply(object.list, ncol))/2))
		anchors <- FindIntegrationAnchors(object.list = object.list, dims = dims, normalization.method = normalization.method, anchor.features = anchor.features, k.filter = k.filter)
		integrated <- IntegrateData(anchorset = anchors, dims = dims, normalization.method = normalization.method)
		integrated <- ScaleData(integrated, verbose = FALSE)
		integrated@misc <- object@misc
		integrated[[old.assay]] <- object[[old.assay]]
		integrated@reductions <- object@reductions
		integrated@meta.data[[split.by]] <- factor(integrated@meta.data[[split.by]], levels = levels(object@meta.data[[split.by]]))
		return(integrated)
}

CheckPCA <- function(object){
		p1 <- DimPlot(object, reduction = "pca", group.by = "orig.ident")
		w <- 6
		if ( exists("Phase", object@meta.data) ){
				p2 <- DimPlot(object, reduction = "pca", group.by = "Phase")
				p1 <- plot_grid(p1, p2, nrow = 1)
				w <- 12
		}
		ggsave(p1, file = "pcaPlot.pdf", width = w, height = 6)

		dims <- min(20, ncol(Reductions(object, "pca")))
		p3 <- DimHeatmap(object, dims = seq(dims), cells = 500, balanced = TRUE, ncol = 4, fast = FALSE)
		ggsave(p3, file = "pcaHeatmap.pdf", width = min(4, dims) * 4, height = ceiling(dims/4) * 4)
		if( !is.null(dev.list()) ) dev.off()

		p4 <- ElbowPlot(object, ndims = ncol(Reductions(object, "pca")))
		ggsave(p4, file = "pcaElbowPlot.pdf", width = 6, height = 6)
}

DoFindClusters <- function(object, reduction = "pca", dims = NULL, resolution = 0.5){
		if ( is.null(dims) || max(dims) > length(object[[reduction]]) ) dims <- seq(object[[reduction]])
		object <- FindNeighbors(object, reduction = reduction, dims = dims, force.recalc = TRUE)
		object <- FindClusters(object, resolution = resolution, temp.file.location = getwd())

		color.cluster <- rainbow(length(levels(object@meta.data$seurat_clusters)))
		names(color.cluster) <- levels(object@meta.data$seurat_clusters)
		object@misc[["color.cluster"]] <- color.cluster

		return(object)
}


.PlotCluster <- function(object, reduction = NULL, cells = NULL, p1.group.by = "orig.ident", p2.group.by = "seurat_clusters", outfile = NULL, p1.color = object@misc$color.sample, p2.color = object@misc$color.cluster, p1.label = FALSE, p2.label = TRUE, ...){
		p1 <- DimPlot(object, reduction = reduction, cells = cells, group.by = p1.group.by, cols = p1.color, label = p1.label, ...)
		p2 <- DimPlot(object, reduction = reduction, cells = cells, group.by = p2.group.by, cols = p2.color, label = p2.label, ...)
		if ( exists("dot_theme_default") && "showtext" %in% installed.packages() ) {
				findplottitle <- function(object, reduction = NULL, group.by = NULL){
						if ( is.null(reduction) ) {
								reduction <- Seurat:::DefaultDimReduc(object = obj)
						}
						if ( grepl("umap", reduction, ignore.case = T) ) {
								reduction <- "UMAP"
						}else if ( grepl("tsne", reduction, ignore.case = T) ) {
								reduction <- "tSNE"
						}else if ( grepl("pca", reduction, ignore.case = T) ) {
								reduction <- "PCA"
						}
						title <- paste(reduction, "plot")
						if ( ! is.null(group.by) ) {
								if ( grepl("orig.ident|sample", group.by, ignore.case = T) ) {
										group.by <- "sample"
								}else if ( grepl("cluster", group.by, ignore.case = T) ) {
										group.by <- "cluster"
								}
								title <- paste(title, "for", group.by)
						}
						return(title)
				}
				p1 <- p1 + ggtitle(findplottitle(object, reduction, p1.group.by)) + dot_theme_default()
				p2 <- p2 + ggtitle(findplottitle(object, reduction, p1.group.by)) + dot_theme_default()
		}
#		p  <- plot_grid(p1, p2)
		p <- p1 + p2
		if ( is.null(outfile) ) {
				return(p)
		}else{
				ggsave(p, file = outfile, width = 12, height = 6 )
		}
}

PlotCluster <- function(object, reduction = 'umap', split.by = "orig.ident", outpref = NULL){
		.PlotCluster(object, reduction = reduction, outfile = paste0(outpref, ".pdf"))
		for ( i in unique(object@meta.data[[split.by]]) ){
				cells.use <- rownames(object@meta.data)[object@meta.data[[split.by]] == i]
				.PlotCluster(object, reduction = reduction, cells = cells.use, p1.group.by = split.by, outfile = paste0(outpref, "_", i, ".pdf"))
		}
}

WriteTable <- function(x, file){
		write.table(x, file = file, quote = F, sep = "\t", row.names = FALSE, col.names = TRUE  )
}

StatCluster <- function(object){
		.StatCluster_Num(object)
		.StatCluster_InSample(object)
		.PlotFraction_CLSP(object)
}

.StatCluster_Num <- function(object){
		Cluster.stat <- object@meta.data %>% group_by(Cluster = seurat_clusters) %>% 
				summarise('Cells number' = n(), 'Median Genes per Cell' = median(nFeature_RNA), 'Median UMI Counts per Cell' = median(nCount_RNA))
		WriteTable(Cluster.stat, "Cluster.stat.xls")
}

.StatCluster_InSample <- function(object){
		Cluster.sample.stat <- object@meta.data %>% group_by(Sample = orig.ident, Cluster = seurat_clusters) %>%
				summarise(y = n()) %>% mutate(Cluster = factor(Cluster, levels = c("Total", levels(Cluster)))) %>% 
				full_join(x = object@meta.data %>% group_by(Sample = orig.ident) %>% summarise(Cluster = factor("Total", levels = c("Total", levels(seurat_clusters))), y = sum(n())) ) %>%
				reshape2::dcast(Cluster ~ Sample, fill = 0) %>%
				mutate_if(is.numeric, list(~paste0(., " (", round(./.[1] * 100,2), "%)")))
#				mutate_if(is.numeric, funs(paste0(., " (", round(./.[1] * 100,2), "%)")))
		WriteTable(Cluster.sample.stat, "Cluster.sample.stat.xls")
}

.PlotFraction_CLSP <- function(object, color.cluster = NULL, color.sample = NULL){
		if ( is.null(color.cluster) ) color.cluster <- object@misc$color.cluster
		if ( is.null(color.sample) )  color.sample  <- object@misc$color.sample
		stat_sample <- object@meta.data %>% group_by(Sample = orig.ident, Cluster = seurat_clusters) %>% summarise("Number of cells" = n()) 
		p.sample  <- ggplot(stat_sample, aes(x = Sample,  y = `Number of cells`, fill = Cluster)) + scale_fill_manual(values = color.cluster) + theme_light()
		p.cluster <- ggplot(stat_sample, aes(x = Cluster, y = `Number of cells`, fill = Sample))  + scale_fill_manual(values = color.sample ) + theme_light()
		geom_stack <- geom_bar(stat = "identity", position = 'stack')
		geom_fill  <- geom_bar(stat = "identity", position = "fill" )
		ggsave( p.sample  + geom_stack, file = "Cluster.sample.stat.pdf",      height = 6, width = 8 )
		ggsave( p.sample  + geom_fill,  file = "Cluster.sample.stat.pct.pdf",  height = 6, width = 8 )
		ggsave( p.cluster + geom_stack, file = "Cluster.cluster.stat.pdf",     height = 6, width = 8 )
		ggsave( p.cluster + geom_fill,  file = "Cluster.cluster.stat.pct.pdf", height = 6, width = 8 )
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
#				rownames(mean_exp) <- ChangeOUTName(rownames(mean_exp), object@misc$fdata)
#				mean_exp <- cbind(Gene_ID = rownames(mean_exp), Gene_name = object@misc$fdata[rownames(mean_exp), "name"], mean_exp)
				Gene_name <- FindFeaturesName(object, rownames(mean_exp), "name")
				mean_exp <- cbind(Gene_ID = ChangeOUTName(rownames(mean_exp), object@misc$fdata), Gene_name = Gene_name, mean_exp)
				WriteTable(mean_exp, file = outfile)
		}
}

CalPctExp <- function(object, features = NULL, group.by = NULL, assay = NULL, slot = "counts", is.return = FALSE, is.reverse = FALSE, is.bulk = FALSE ){
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
				mean_exp <- cbind(Gene_ID = rownames(mean_exp), Gene_name = object@misc$fdata[rownames(mean_exp), "name"], mean_exp)
#				WriteTable(mean_exp, file = "AllGene.avg_exp.xls")
		}
}

.GetMetaData <- function(object, cols = NULL) {
		name <- names(cols)		
		name[is.na(name) | name == ""] <- cols[is.na(name) | name == ""]
		if ( is.null(name) ) name <- cols						
		cells_list <- object@meta.data[, cols]
		colnames(cells_list) <- name
		cells_list <- cbind(Cells = rownames(cells_list), cells_list)
		return(cells_list)
}
ListCellCluster <- function(object){
		data <- .GetMetaData(object, cols = c("Sample" = "orig.ident", "Cluster" = "seurat_clusters"))
		WriteTable(data, file = "Cells.cluster.list.xls")
}

PlotFeaturePlot <- function(object, features, outfile = NULL, reduction = NULL, is.use.name = TRUE, color.high = "blue", color.low = "lightgrey", show.cluster.label = FALSE, nCol = NULL, plot.basic.size = 4, group.by = "seurat_clusters" ) {
		if ( show.cluster.label ) Idents(object) <- group.by
		plots <- FeaturePlot(object, features = features, order = TRUE, reduction = reduction, combine = FALSE, label = show.cluster.label, cols = c(color.low, color.high))
		if ( is.use.name ) {
				name <- FindFeaturesName(object, features)
				for ( i in seq(plots) ){
						plots[[i]] <- plots[[i]] + ggtitle(name[i])
				}
		}

		if ( is.null(nCol) ) nCol <- ceiling(sqrt(length(features)))
		nRow <- ceiling(length(features) / nCol)

		p <- CombinePlots(plots = plots, ncol = nCol)
		if ( is.null(outfile) ) {
				return(p)
		} else {
				ggsave(p, file = outfile, width = plot.basic.size * (5/4) * nCol, height = plot.basic.size * nRow, limitsize = FALSE)
		}
}


PlotDotPlot <- function(object, features = NULL, outfile = NULL, group.by = "seurat_clusters", is.use.name = TRUE, ...){
		p <- DotPlot(object, features = features, group.by = group.by, ...) + RotatedAxis()
		if ( is.use.name ) {
				levels(p$data$features.plot) <- FindFeaturesName(object, levels(p$data$features.plot)) 
		}
		if ( is.null(outfile) ) {
				return(p)
		} else {
				w <- max(6, ceiling(length(features)) * 0.35 + 2)
				h <- max(6, length(unique(object@meta.data[[group.by]])) * 0.4)
				ggsave(p, file = outfile, width = w, height = h, limitsize = FALSE )
		}
}

PlotHeatmapPlot <- function(object, features = NULL, group.by = "seurat_clusters", is.use.name = TRUE, outfile = NULL) {
		group.colors <- if ( group.by == "seurat_clusters" ) {
								object@misc$color.cluster
						}else if ( group.by == "orig.ident" ) {
								object@misc$color.sample
						}else{
								NULL
						}
		p <- DoHeatmap( object = object, features = features, cells = NULL,
						group.by = group.by, group.colors = group.colors, 
						combine = FALSE, raster = FALSE)
		p <- p[[1]]
		p <- p + theme(legend.title = element_blank())
		p$layers[[2]] <- NULL

		if ( is.use.name ) levels(p$data$Feature) <- FindFeaturesName(object, levels(p$data$Feature))

#		p <- ggplotGrob(p)
#		for (i in grep('strip', p$layout$name)){
#				p$grobs[[i]]$layout$clip <- "off"
#		}

		if ( is.null(outfile) ) {
				return(p)
		} else {
				h <- max(7, length(unique(features)) * 0.11 + 2.5 )
				w <- h * 4 / 3
				ggsave(p, file = outfile, width = w, height = h, limitsize = FALSE )
		}
}

PlotAboutFeatures <- function(object, features = NULL, outpref = NULL ) {
		PlotDotPlot(object, features = features, outfile = paste0(outpref, ".DotPlot.pdf"))
		PlotFeaturePlot(object, features = features, outfile = paste0(outpref, ".Distribution.pdf"), reduction = "tsne")
		PlotHeatmapPlot(object, features = features, outfile = paste0(outpref, ".Heatmap.pdf"))
}

PlotPresetMarker <- function(object){
		PresetMarker <- union(object@misc[["expected.marker"]], object@misc[["more.marker"]])
		PresetMarker <- union(PresetMarker, object@misc[["exclude.marker"]])
		if ( length(PresetMarker) ) {
				VlnplotShell(object, features = PresetMarker, outfile = "PresetMarker.VlnPlot.pdf", titles = object@misc$fdata[PresetMarker, "merge_name"], cols.use = object@misc$color.cluster, group.by = "seurat_clusters")
				PlotAboutFeatures(object, features = PresetMarker, outpref = "PresetMarker")
		}
}

DoFindAllMarkers <- function(object, parameter, group.by = "seurat_clusters") {
		Idents(object) <- group.by
		object.markers <- FindAllMarkers(object = object, only.pos = TRUE,
						min.pct = parameter$FindMarkers$min_pct, logfc.threshold = parameter$FindMarkers$logfc,
						return.thresh = parameter$FindMarkers$pvalue, pseudocount.use = 0 )
		return(object.markers)
}


StatMarker <- function(object.markers, Cluster_name = "Up Gene Number", color = NULL, outpref = "DeGene.stat"){
		stat <- cbind( Cluster = Cluster_name, t(table(object.markers$cluster)) )
		WriteTable(stat, paste0(outpref, ".xls"))

		p <- ggplot(object.markers) + geom_bar(aes(x = cluster, fill = cluster), stat = "count") +
				theme_light() + labs(x = "Cluster", y = "Number of DE genes") 
		if ( !is.null(color) ) p <- p + scale_fill_manual(values = color)
		ggsave(p, file = paste0(outpref, ".pdf"), height = 6, width = 8 )
		return(stat)
}


ListMarker <- function(object, object.markers, outfile = "DeGene.list.xls", is.return = FALSE, is.fast = FALSE, group.by = "seurat_clusters", assay = DefaultAssay(object), slot = "data"){
		#assay <- DefaultAssay(object)
		Targets_mean <- CalAvgExp(object, unique(object.markers$gene), is.return = T, is.reverse = F, assay = assay, slot = slot, group.by = group.by) %>%
				reshape2::melt(varnames = c("gene", "cluster"), value.name = "Target_Cluster_mean") %>%
				mutate(cluster = factor(cluster))
		Others_mean  <- CalAvgExp(object, unique(object.markers$gene), is.return = T, is.reverse = T, assay = assay, slot = slot, group.by = group.by) %>%
				reshape2::melt(varnames = c("gene", "cluster"), value.name = "Other_Cluster_mean") %>%
				mutate(cluster = factor(cluster))
		Name <- FindFeaturesName(object, unique(object.markers$gene), "name", is.fast = is.fast)
		marker_list <- object.markers %>% left_join(y = Targets_mean) %>% left_join(y = Others_mean) %>% 
				mutate(Log2FC = log2(Target_Cluster_mean / Other_Cluster_mean), name = Name[gene]) %>%
				select("Target Cluster" = cluster, "Gene ID" = gene, "Gene Name" = name,
						Target_Cluster_mean, Other_Cluster_mean, Log2FC,
						Pvalue = p_val, Qvalue = p_val_adj)
		if ( is.null(outfile) || is.return ){
				return(marker_list)
		}else{
				marker_list[["Gene ID"]] <- ChangeOUTName(marker_list[["Gene ID"]], object@misc$fdata)
				WriteTable(marker_list, outfile)
		}
}

FindTopMarker <- function(object.markers, top_num = 20, object = NULL, outfile = "Top.avg_exp.xls.tmp"){
		top <- object.markers %>% group_by( cluster ) %>%
#				top_n(top_num, avg_logFC)
				arrange(desc(avg_logFC), p_val, p_val_adj, .by_group = TRUE) %>% filter(1:n() <= top_num)
		if ( ! is.null(outfile) ) {
				tmp <- top %>% select(Cluster = cluster, Gene_ID = gene)
				if ( ! is.null(object) ) {
						tmp$Gene_ID <- ChangeOUTName(tmp$Gene_ID, object@misc$fdata)
				}
				WriteTable(tmp, outfile )
		}
		return(top)
}
###############

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

