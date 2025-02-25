### function defined ###

## Message func
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


Version_check <- function(parameter, version = NULL){
	if(parameter$version != version) StopMessage("Use the correct version\n")
}

LibraryPackages <-function(packages, lib.loc = NULL)
{
	RunMessage('>>>> Start Load Rpackages <<<<')
	RunMessage('   ---> Loading Rpackages <---  ')
	for (i in packages){
		Message(paste0("    ~~~> ", i, " <~~~"))
		library(i, character.only=T)
	}
}


SetClusters <- function(object, ident.old, ident.new, col.names="single_clusters", set.ident=F)
{
	if(is.null(ident.old) | is.null(ident.new)){
		StopMessage("   ---> ident.old and ident.new must be defined <---   ")
	}else{
		for (i in seq(ident.old)){
			object@meta.data[[col.names]][object@active.ident %in% ident.old[[i]]] <- ident.new[i]
		}
	}
	object@meta.data[[col.names]] <- factor(object@meta.data[[col.names]])
	if(set.ident){
		Idents(object) <- col.names
	}
	return(object)
}


SetIdents <- function(object, col.name=NULL)
{
	if(! (col.name %in% colnames(object@meta.data))){
		StopMessage("   --> parameter col.names must be any one of object meatdata colnames <--  ")
	}

	Idents(object) <- col.name
	object@active.ident <- factor(object@active.ident)
	object@meta.data$Clusters <- NULL
	object <- AddMetaData(object=object, metadata = object@active.ident, col.name = 'Clusters')
	return(object)
}


Setlegend <- function(
	myPlot,
	pointSize = 1,
	textSize = 3,
	spaceLegend = 0.3,
	ncol = 2,
	legendtype = 0
){
	myPlot <-
		myPlot +
		theme(legend.text  = element_text(size = textSize),legend.key.size = unit(spaceLegend, "lines"))
	if(legendtype){
		myPlot <- myPlot + guides(colour = guide_legend(ncol = ncol,override.aes=list(size=pointSize)))
	}else{
		myPlot <- myPlot + guides(fill = guide_legend( ncol = ncol, byrow = TRUE))
	}
	return(myPlot)
}


ColorPalette <- function(object, color_by="orig.ident", cluster.colorpalette = FALSE, sample.colorpalette = FALSE){
	if(! color_by %in% colnames(object@meta.data)){
		cat("[FATAL ERROR]: color_by must be any one of ", "[ ", colnames(object@meta.data), " ]\n")
		quit(status=1)
	}

	if(! is.factor(object@meta.data[[color_by]])) object@meta.data[[color_by]] <- factor(object@meta.data[[color_by]])
	if(cluster.colorpalette & sample.colorpalette){
		cat("[FATAL ERROR]: Only one type color can be paletted at the same time\n")
		quit(status=1)
	}
	if(cluster.colorpalette){
		color <- fetch_color(nlevels(object@meta.data[[color_by]]),"tsne", "set1")
		object@misc[[color_by]] <- color
	}
	if(sample.colorpalette){
		color <- fetch_color(nlevels(object@meta.data[[color_by]]),"tsne", "set3")
		object@misc[[color_by]] <- color
	}
	return(object)
}

### Create Summarized Experiment for SingelR
CreateSummarizedExperiment <- function(parameter, sample_name = IfNull(parameter$sample_name, 'orig.ident'), cluster_name = IfNull(parameter$cluster_name, 'seurat_clusters'))
{
	RunMessage("   ---> Start Create Summarized Experiment <---   ")
	RunMessage("    --> Load object ... <--   ")
	tryCatch( { read.obj <- load(parameter$obj) }, error = function(e) { read.obj <- readRDS(parameter$obj) } )
	if ( class(read.obj) == "character" ) {
		object <- get(read.obj[[1]])
	}else {
		object <- read.obj
	}

	if(grepl('^2', object@version)){
		RunMessage("    --> Update objcet ... <--    ")
		object <- UpdateSeuratObject(object)
		object@meta.data$seurat_clusters <- object@meta.data$cluster
	}
	object@meta.data$orig.ident <- object@meta.data[[sample_name]]
	object@meta.data$seurat_clusters <- object@meta.data[[cluster_name]]


	RunMessage("    --> Check Input Samples ... <--    ")
	cells.use <- NULL
	if(!is.null(parameter$samples.use)){
		samples.use <- parameter$samples.use
		if(! all(samples.use %in% unique(object@meta.data$orig.ident))){
			StopMessage("    --> samples.use must be object meatdata orig.ident <--    ")
		}
		cells.use <- Cells(object)[object@meta.data$orig.ident %in% samples.use]
	}

	RunMessage("    --> Check Input idents ... <--    ")
	idents.use <- NULL
	if(!is.null(parameter$idents.use)){
		idents.use <- parameter$idents.use
		if(! all(idents.use %in% levels(object@active.ident))){
			StopMessage("    --> samples.use must be any one of object active ident or colnames of object metadata <--    ")
		}
	}
	
	if(! is.null(cells.use) | !is.null(idents.use)){
		RunMessage("    --> Subset object ... <--    ")
		object <- subset(object, cells=cells.use, idents=idents.use)
		object@meta.data$orig.ident <- factor(object@meta.data$orig.ident)
		object@meta.data$seurat_clusters <- factor(object@meta.data$seurat_clusters)
		object@active.ident <- factor(object@active.ident)
#		vars.regress <- FindVarsRegress(object,parameter)
#		object <-DoNormalization(object, parameter, vars.regress = vars.regress)
	}

	RunMessage("    --> Create Summarized Experiment ... <--    ")
	logcounts <- GetAssayData(object,slot='data')
	if (!is.null(parameter$name_aln)){
			name_aln <- read.table(parameter$name_aln, header = FALSE, row.names = 1, sep = "\t")
			nudup.id <- rownames(name_aln)[!duplicated(name_aln[[1]])]
			nudup.id <- intersect(rownames(logcounts), nudup.id)
			logcounts <- logcounts[nudup.id, ]
			rownames(logcounts) <- as.vector(name_aln[rownames(logcounts), 1])
	}else if(exists("fdata",object@misc)){
		if(! all(rownames(logcounts) %in% rownames(object@misc$fdata))){
			rownames(object@misc$fdata) <- object@misc$fdata$dash
		}
		nudup.id <- rownames(object@misc$fdata[!duplicated(object@misc$fdata$name), ])
		nudup.id <- intersect(rownames(logcounts), nudup.id)
		logcounts <- logcounts[nudup.id, ]
		rownames(logcounts) <- as.vector(object@misc$fdata[rownames(logcounts),"name"])
	}
	rownames(logcounts) <- toupper(rownames(logcounts))

	sc_data <- SummarizedExperiment(
			list(logcounts=logcounts),
			colData=DataFrame(cluster=object@meta.data$seurat_clusters)
			)

	object@assays$Single$data <- sc_data
	return(object)
}



AutoFindRef <- function(parameter = list(), ref.loc=parameter$reference$dir,
				spe = parameter$reference$species, tissue = parameter$reference$tissue){
	spe <- gsub(" ", "_", spe)
#	accept_species <- c("Homo_sapiens", "Mus_musculus", 'Oryza_sativa_indica', 'Oryza_sativa_japonica', 'Arabidopsis_thaliana')
#	if( ! (spe %in% accept_species))
#			StopMessage(paste0("---> Species must be [", accept_species, "]\n"))
	ref <- paste(c(ref.loc, spe, paste0(tissue,".rds")), collapse="/")
	if (!(file.exists(ref)))
			StopMessage(paste0(ref, " reference no exists !!!"))
	RunMessage(paste0("    [ Reference ] : ", ref))
	refdata <- readRDS(ref)
	return(refdata)
}


RunSingleR <- function(object, parameter = list(),
	label.col = parameter$reference$label.col, method = parameter$method,
	ref.loc=parameter$reference$dir, spe = parameter$reference$species, tissue = parameter$reference$tissue)
{
	RunMessage("   ---> Run SingleR <---   ")
	RunMessage("    --> Import Reference Data ... <--    ")
	ref.data <- AutoFindRef(ref.loc = ref.loc, spe = spe, tissue = tissue)
	rownames(ref.data) <- toupper(rownames(ref.data))
	object@assays$Single$refdata <- ref.data
	RunMessage("    --> Start SingleR ... <--    ")
	if(method == "single"){
		RunMessage("    --> Cell Annotated by Method~single ... <--    ")
		single.out <- NULL
		single.out <- SingleR(object@assays$Single$data, object@assays$Single$refdata, labels=object@assays$Single$refdata[[label.col]], method="single")
		object@meta.data$single_clusters <- single.out$labels
	}else if (method == "cluster"){
		RunMessage("    --> Cell Annotated by Method~cluster ... <--    ")
		single.out <- NULL
		single.out <- SingleR(object@assays$Single$data, object@assays$Single$refdata, labels=object@assays$Single$refdata[[label.col]], method=parameter$method, clusters=object@meta.data$seurat_clusters)
		single.out$cluster <- rownames(single.out)
		object <- SetClusters(object, single.out$cluster, single.out$labels, col.name="single_clusters", set.ident=F)
	}else{
		StopMessage("    --> Method must be one of [single | cluster] <--    ")
	}
	prob <- single.out$scores
	rownames(prob) <- rownames(single.out)
	object@misc$singleR_score <- prob
	return(object)
}



ClusterPlots <- function(object, reduction = 'tsne', cells = NULL, p1.group.by = "orig.ident", p2.group.by = "seurat_clusters", outpref = NULL, p1.color = NULL, p2.color = NULL, p1.label = FALSE, p2.label = TRUE, p1.label.size = 2, p2.label.size = 2, split.by = 'orig.ident', pdf.height = 7, pdf.width = pdf.height * 2 * 1.2, p1.legend.col = 1, p2.legend.col = 1, ...){
		p1 <- DimPlot(object, reduction = reduction, cells = cells, group.by = p1.group.by, cols = p1.color, label = p1.label, label.size = p1.label.size, ...)
		p1 <- Setlegend(p1, pointSize = 3, textSize = 6, ncol = p1.legend.col, spaceLegend = 0.3, legendtype=1)
		p2 <- DimPlot(object, reduction = reduction, cells = cells, group.by = p2.group.by, cols = p2.color, label = p2.label, label.size =p2.label.size, ...)
		p2 <- Setlegend(p2, pointSize = 3, textSize = 6, ncol = p2.legend.col, spaceLegend = 0.3, legendtype=1)
		p <- (p1 + p2) & dot_theme_default()
		ggsave(p, file = paste0(outpref, ".pdf"), width = pdf.width, height = pdf.height )
		if(!is.null(split.by)){
			for ( i in unique(object@meta.data[[split.by]]) ){
				cells.use <- rownames(object@meta.data)[object@meta.data[[split.by]] == i]
				p1 <- DimPlot(object, reduction = reduction, cells = cells.use, group.by = p1.group.by, cols = p1.color, label = p1.label, label.size = p1.label.size, ...)
				p1 <- Setlegend(p1, pointSize = 3, textSize = 6, ncol = p1.legend.col, spaceLegend = 0.3, legendtype=1)
				p2 <- DimPlot(object, reduction = reduction, cells = cells.use, group.by = p2.group.by, cols = p2.color, label = p2.label, label.size = p2.label.size, ...)
				p2 <- Setlegend(p2, pointSize = 3, textSize = 6, ncol = p2.legend.col, spaceLegend = 0.3, legendtype=1)
				p <- (p1 + p2) & dot_theme_default()
				ggsave(p, file = paste0(outpref, "_", i, ".pdf"), width = pdf.width, height = pdf.height )
			}
		}
}


ListMetaData <- function(object, cols=c("Sample" = "orig.ident","Seurat_Cluster" = "seurat_clusters"))
{
	data <- .GetMetaData(object, cols = cols)
	WriteTable(data, file = "Cells.metadata.list.xls")
}


corstat <- function(object, old_cluster = "seurat_clusters", new_cluster = "single_clusters", outfile = "Cluster.correlation.heatmap.xls"){
	exp.seurat <- CalAvgExp(object, group.by=old_cluster,is.return=T)
	exp.sing <- CalAvgExp(object, group.by=new_cluster,is.return=T)
	cor <- cor(exp.seurat, exp.sing)
	if ( is.null(outfile) ) {
			return(cor)
	}else{
			write.table(cbind(cluster=rownames(cor), cor), row.names=F, col.names=T, quote=F, sep="\t", file=outfile)
	}
}




