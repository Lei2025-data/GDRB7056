
######################
Read_yaml <- function(file, outdir_dir = NULL){
		handlers <- list("bool#no" = function(x){if ( x %in% c("false", "FALSE") ) FALSE else x}, "bool#yes" = function(x){if ( x %in% c("true", "TRUE") ) TRUE else x})
		parameter <- yaml::yaml.load_file( file, handlers = handlers )

		for ( i in names(parameter$infile) ) {
				parameter$infile[[i]] <- paste(parameter$project_dir, "indir", parameter$infile[[i]], sep = "/")
		}
		parameter$outdir_dir <- ifelse(is.null(outdir_dir), parameter$project_dir, outdir_dir)

		parameter$outdir <- list()
		parameter$outdir$indir <- paste(parameter$outdir_dir, "indir", sep = "/")
		parameter$outdir$list  <- paste(parameter$outdir_dir, "list", parameter$list_name, sep = "/")
		parameter$outdir$task <- list()
		for ( i in c("cluster", "marker", "differ", "group") ) {
				parameter$outdir$task[[i]] <- paste(parameter$outdir_dir, "task", parameter$task_name, i, sep = "/")
		}

		for ( i in names(parameter$samples) ) {
				parameter$new_sample[[i]] <- parameter$samples[[i]][1]
				parameter$new_group[[parameter$samples[[i]][1]]] <- parameter$samples[[i]][2]
		}

		i <- unlist(lapply(seq(parameter$two_groups_diff), function(x) if(parameter$two_groups_diff[[x]][2] == "others") x ))
		if ( length(i) ) parameter$two_groups_diff[i] <- NULL

		if ( parameter$filter$nUMI[1]     == "INF" ) parameter$filter$nUMI[1]     <- -Inf
		if ( parameter$filter$nGene[1]    == "INF" ) parameter$filter$nGene[1]    <- -Inf
		if ( parameter$filter$pct.mito[1] == "INF" ) parameter$filter$pct.mito[1] <- -Inf
		parameter$filter$nUMI     <- as.numeric(parameter$filter$nUMI)
		parameter$filter$nGene    <- as.numeric(parameter$filter$nGene)
		parameter$filter$pct.mito <- as.numeric(parameter$filter$pct.mito)
		parameter$filter$nCount_RNA   <- parameter$filter$nUMI
		parameter$filter$nFeature_RNA <- parameter$filter$nGene
		parameter$filter$percent.mito <- parameter$filter$pct.mito
		if ( ! exists("set.num", parameter$filter) ) parameter$filter$set.num <- "none"

		parameter$differ$logFC <- eval(parse(text = parameter$differ$logFC))
		parameter$group$logFC  <- eval(parse(text = parameter$group$logFC))

		if ( ! is.null(parameter$differ$qvalue) ) parameter$differ$pvalue <- parameter$differ$qvalue
		parameter$FindMarkers <- parameter$differ


		return(parameter)
}

Read_yaml_diff <- function(file){
		handlers <- list("bool#no" = function(x){if ( x %in% c("false", "FALSE") ) FALSE else x}, "bool#yes" = function(x){if ( x %in% c("true", "TRUE") ) TRUE else x})
		parameter <- yaml::yaml.load_file( file, handlers = handlers )

		for ( i in names(parameter$samples) ) {
				parameter$new_sample[[i]] <- parameter$samples[[i]][1]
				parameter$new_group[[parameter$samples[[i]][1]]] <- parameter$samples[[i]][2]
		}

		return(parameter)
}

LoadObject <- function(parameter, type = c("origin", "remake", "rename", "readlist", "cluster", "remaster")){
		type <- match.arg(type)
		if ( type == "origin" ) {
				object <- Load(parameter$infile$obj_file)
		}else if ( type == "rename" ) {
				object <- Load(parameter$infile$obj_file)
				object <- DoRenameCells(object, parameter)
		}else if ( type == "remake" ) {
				object <- Load(parameter$infile$obj_file)
				object <- RemakeObject(object)
		}else if ( type == "readlist" ) {
				object <- Load(parameter$infile$obj_file)
				object <- RemakeObject(object)
				dt <- read.table(paste(parameter$outdir$list, "cells.list.xls", sep = "/"),
						header = TRUE, stringsAsFactors = FALSE, sep = "\t")
				if ( ! all(colnames(object) %in% dt$oldname) )
						object <- object[, dt$oldname]
				object <- DoRenameCells(object, parameter)
		}else if ( type == "cluster" ) {
				object <- Load(paste(parameter$outdir$task$cluster, "obj.Rda", sep = "/"))
		}else if ( type == "remaster" ) {
				oldfile <- paste0(parameter$infile$obj_file, ".bak")
				if ( ! file.exists(oldfile) ) {
						if ( file.exists(parameter$infile$obj_file) ) {
								file.rename(parameter$infile$obj_file, oldfile)
						} else {
								stop("Can't not find file 'obj.Rda'.")
						}
				}else{
						if ( file.exists(parameter$infile$obj_file) )
								file.remove(parameter$infile$obj_file)
				}
				object <- Load(oldfile)
				if ( length(object@misc) == 0 ) {
						object <- RemasterObject(object, parameter = parameter)
				}
				save(object, file = parameter$infile$obj_file)
				file.remove(oldfile)
		}
		return(object)
}



RemasterObject <- function(object, name_list_file = NULL, parameter = NULL){
		if ( is.null(name_list_file) ) name_list_file <- parameter$infile$name_list_file
		object@misc[["fdata"]] <- AddFData(object, ref_name_file = name_list_file, is.oldstyle = F)
		object@misc[["fdata"]]$old_merge_name <- gsub(pattern = "_", replacement = "-", x = object@misc[["fdata"]]$old_merge_name)
		object <- RenameFeatures(object, from.type = "old_merge_name", to.type = "id")
		new.names <- gsub("_", "-", rownames(object))
		names(new.names) <- rownames(object)
		object <- RenameFeatures(object, new.names)
		object@misc[["pdata"]] <- object@meta.data
		object@misc[["counts"]] <- GetAssayData(object, slot = "counts", assay = "RNA")

		color.sample <- scales::hue_pal()(length(levels(object@meta.data$orig.ident)))
		names(color.sample) <- levels(object@meta.data$orig.ident)
		object@misc[["color.sample"]] <- color.sample
		
		if ( ! exists("seurat_clusters", object@meta.data) ) {
				if ( exists("cluster", object@meta.data) ) {
						object@meta.data$seurat_clusters <- object@meta.data$cluster
				}
		}

		if ( exists("pct.mito", object@meta.data) ) {
				object@meta.data$pct.mito <- object@meta.data$pct.mito * 100
				object@meta.data$percent.mito <- object@meta.data$pct.mito
		} else if ( exists("percent.mito", object@meta.data) ) {
				object@meta.data$percent.mito <- object@meta.data$percent.mito * 100
		}

		return(object)
}

UpdataParameter <- function(parameter){
		parameter$version <- 3.0
		parameter$cluster$is.integration <- TRUE
		parameter$infile <- list(
						obj_file       = "obj.Rda",
						differ_file    = "markers.Rda",
						name_list_file = "name_list.xls",
						seq_info_file  = "samples.sequence.stat.xls",
						ref_pref       = "annot/ref",
						violin_data    = "filter_violin_data.xls",
						matrix         = "expression.xls")
		parameter$outdir <- NULL
		return(parameter)
}
######################



######################
DoRenameCells <- function (object, parameter, new.cell.names = NULL) {
#		if ( is.null(new.cell.names) ) new.cell.names <- .find_new_name(colnames(object), parameter$new_sample)
#		object <- RenameCells(object, new.names = new.cell.names)

		object <- RenameGrouping(object, parameter)
#		object@meta.data$orig.ident <- as.factor(object@meta.data$orig.ident)
#		levels(object@meta.data$orig.ident) <- parameter$new_sample[levels(object@meta.data$orig.ident)]

#		object@meta.data$group <- object@meta.data$orig.ident
#		levels(object@meta.data$group) <- parameter$new_group[levels(object@meta.data$group)]

#		names(object@misc$color.sample) <- as.vector(parameter$new_sample[names(object@misc$color.sample)])

		return(object)
}

RenameGrouping <- function (object, parameter = list(), new_sample = parameter$new_sample, new_group = parameter$new_group) {

		if ( ! is.null(new_sample) ) {
				object@meta.data$orig.ident <- factor(new_sample[as.character(object@meta.data$orig.ident)], levels = new_sample)
				names(object@misc$color.sample) <- as.vector(new_sample[names(object@misc$color.sample)])
		}
		if ( ! class(object@meta.data$orig.ident) == "factor" ) 
				object@meta.data$orig.ident <- as.factor(object@meta.data$orig.ident)

		object@meta.data$group <- object@meta.data$orig.ident
		if ( ! is.null(new_group) ) {
				levels(object@meta.data$group) <- new_group[levels(object@meta.data$group)]
		}

		return(object)
}


RemakeObject <- function(object){
		object <- RestoreObject(object)
		if ( !is.null(object@misc$percent.mito) ) {
#				object <- StatFeatures(object, object@misc$percent.mito, col.name = "percent.mito", stat_pct = T)
		}else if ( exists("percent.mito", object@misc$pdata) ) {
				object@meta.data$percent.mito <- 0
				cells <- rownames(object@meta.data) %in% rownames(object@misc$pdata)
				object@meta.data$percent.mito[cells] <- object@misc$pdata[rownames(object@meta.data)[cells], "percent.mito"]
		}
		return(object)
}



ViolinData <- function(object, outfile = "filter_violin_data.xls") {
		col <- c("Sample" = "orig.ident", "nGene" = "nFeature_RNA", "nUMI" = "nCount_RNA")
		if ( exists("group", object@meta.data) ) col <- c(col, "Group" = "group")
		if ( exists("percent.mito", object@meta.data) ) col <- c(col, "percent.mito" = "percent.mito")
		dt <- object@meta.data[,col]
		colnames(dt) <- names(col)
		dt <- cbind(cell = rownames(dt), dt)
		write.table(dt, file = outfile, quote = F, sep = "\t", row.names = F)
}



AsMatrix <- function(object, outfile = "expression.xls", assay = "RNA", slot = "counts", max.row = 1000){
		DefaultAssay(object) <- assay
		header <- c("GeneID", "GeneName", colnames(object))
		WriteTable(t(header), outfile, col.names = FALSE)
#		k <- cut(1:ncol(object), breaks = 0:ceiling(ncol(object) / max.col) * ncol(object))
		k <- cut(1:nrow(object), breaks = 0:ceiling(nrow(object) / max.row) * max.row)
		for ( i in levels(k) ) {
				message("as.matrix : ", paste0(range(which(k == i)), collapse = "-") )
				dt <- as.matrix(GetAssayData(object[which(k == i), ], assay = assay, slot = slot))
				rownames(dt) <- ChangeOUTName(rownames(dt), object@misc$fdata)
				dt <- cbind(GeneID = rownames(dt), GeneName = FindFeaturesName(object, rownames(dt), "name"), dt)
				WriteTable(dt, outfile, append = TRUE, col.names = FALSE)
		}
}
############################

############################
DoFilterCells <- function(parameter, outtype = c("list", "violin")){
		dt <- read.table(parameter$infile$violin_data, header = TRUE, sep = "\t") %>%
			  filter(Sample %in% names(parameter$new_sample)) %>%
			  mutate(oldname = cell,
#					 cell    = .find_new_name(as.character(oldname), parameter$new_sample),	
					 cell    = cell,
					 Sample  = factor(parameter$new_sample[as.character(Sample)], levels = parameter$new_sample),
					 Group   = factor(parameter$new_group[as.character(Sample)], levels = unique(parameter$new_group)),
					 index   = 1:n() + 1) 
		cells.use <- FilterCells(dt, parameter)
		dt <- dt %>% filter(cell %in% cells.use)

		if ( "list" %in% outtype ) {
				cell.dt <- dt %>% select(oldname, cell, index)
				WriteTable(cell.dt, "cells.list.xls")
		}

		if ( "violin" %in% outtype ) {
				dt <- dt %>% select(-oldname, -index)
				WriteTable(dt, "filter_violin_data.xls")
		}
		
		return(dt)
}

FilterCells <- function(object, parameter, return.object = FALSE,
				set.num = IfNull(parameter$filter$set.num, "none"),
				set.num.seed = IfNull(parameter$filter$set.num.seed, 42),
				standard = parameter$filter
				) {
		message( "--->Filter Cells<---" )
		cells.use <- .FilterCells(object, set.num = set.num, set.num.seed = set.num.seed, standard  = standard)

		if ( class(object) == "Seurat" ) {
				if ( ! is.null(parameter$filter$marker$list) ){
						data <- GetAssayData(object, slot = "counts")
								parameter$filter$marker$list <- ChangeOUTName(parameter$filter$marker$list, object@misc$fdata)
								data <- data[parameter$filter$marker$list, ]
								index <- Matrix::colSums(data >= parameter$filter$marker$thres[1] & data <= parameter$filter$marker$thres[2])
						cells.use <- intersect(cells.use, colnames(object)[index == length(parameter$filter$marker$list)])
				}
				if ( return.object ) {
						object <- object[, cells.use]
						return(object)
				}
		}
		return(as.vector(cells.use))
}

StatBasicInfo <- function(metadata, parameter, outfile = "sample.details.xls") {
		dt2 <- read.table(parameter$infile$seq_info_file, header = T, row.names = 1, sep = "\t" ) %>%
				tibble::rownames_to_column(var = "sample") %>%
				mutate(sample = factor(parameter$new_sample[as.character(sample)], levels = parameter$new_sample))  %>% 
				select(sample, "样本饱和度" = Sequencing.Saturation)

		dt <- metadata %>% group_by(sample = Sample) %>%
				summarise(cells = n(),
						  UMI = sum(nUMI),
						  "UMI/cells" = round(UMI/cells, 2),
						  "线粒体UMI比例" = ifelse(exists("percent.mito", metadata), paste0(round(mean(percent.mito), 2), "%"), "0%")
						  ) %>%
				left_join(y = dt2, by = "sample")

		if ( ! is.null(outfile) ) {
				WriteTable(dt, outfile)
		} else {
				return(dt)
		}
}
############################


############################






stat_cluster <- function(object, slot = "data") {
#stat_cluster <- function(object, reduction = "tsne", slot = "data") {
#		embeddings <- object@reductions[[reduction]]@cell.embeddings

		types <- list("orig.ident" = "Sample", "seurat_clusters" = "Cluster","group" = "Group")
		metadata <- data.frame(count = object@meta.data$nCount_RNA)
		for ( type in names(types) ){
				if ( exists(type, object@meta.data) ) {
						metadata <- metadata %>% mutate( !! types[[type]] := object@meta.data[[type]] )

						dt <- CalAvgExp(object, group.by = type, slot = slot, is.return = T)
						colnames(dt) <- paste( types[[type]], colnames(dt), sep = " ")
						rownames(dt) <- ChangeOUTName(rownames(dt), object@misc$fdata)
						dt <- cbind(.change_name2(object, rownames(dt)), dt)
						WriteTable(dt, paste( "all.genes", types[[type]], "stat.xls", sep = "."))

						if ( type == "seurat_clusters" ) {
								pct.cell <- Matrix::rowSums(object@assays$RNA@counts > 0) / ncol(object)
								pct.cell <- paste0(round(pct.cell * 100, 2), "%")
								dt <- cbind(dt, "分布比例" = pct.cell)
								WriteTable(dt, "all.genes.marker.stat.xls")

								stat <- apply(dt[,-1:-2], 2, function(x) c(min = min(x), max = max(x)))
								stat <- cbind(term = rownames(stat), stat)
								WriteTable(stat, "all.genes.marker.stat.min_max.xls")
						}
				}
		}

#		plot.data <- cbind(cell = rownames(embeddings), embeddings, metadata)
#		WriteTable(plot.data, "tSNE.data.xls")

		cluster.sample <- as.data.frame(table(metadata[,c("Sample","Cluster")])) %>% reshape2::dcast(Cluster ~ Sample, value.var = "Freq" )
		WriteTable(cluster.sample, "Cluster.sample.stat.xls")

		cluster.group <- as.data.frame(table(metadata[,c("Group","Cluster")])) %>% reshape2::dcast(Cluster ~ Group, value.var = "Freq" )
		WriteTable(cluster.group, "Cluster.group.stat.xls")
}

get_embeddings <- function(object, reduction = "tsne" ){
		embeddings <- object@reductions[[reduction]]@cell.embeddings

		types <- list("orig.ident" = "Sample", "seurat_clusters" = "Cluster","group" = "Group")
		metadata <- data.frame(count = object@meta.data$nCount_RNA)
		for ( type in names(types) ){
				if ( exists(type, object@meta.data) ) {
						metadata <- metadata %>% mutate( !! types[[type]] := object@meta.data[[type]] )
				}
		}

		outpref <- if ( grepl("tsne", reduction, ignore.case = T) ) 'tSNE'
				else if ( grepl("umap", reduction, ignore.case = T) ) 'UMAP'
				else reduction
		if ( ncol(embeddings) > 2 ) {
				outpref <- paste0(outpref, ".", ncol(embeddings), "d" )
		}

		plot.data <- cbind(cell = rownames(embeddings), embeddings, metadata)
		WriteTable(plot.data, paste0(outpref, ".data.xls"))
}

#########################


stat_diff <- function(object, object.markers, do.plot = TRUE) {
		marker_list <- ListMarker(object, object.markers, is.return = TRUE)
		if ( ! is.null(parameter$differ$qvalue) )
				marker_list <- subset(marker_list, Qvalue < parameter$differ$qvalue )
		marker_list[["Gene ID"]] <- ChangeOUTName(marker_list[["Gene ID"]], object@misc$fdata)
		WriteTable(marker_list, "differ.genes.stat.xls")

		stat <- apply(marker_list[,4:8], 2, function(x) c( min = min(x), max = max(x)))
		stat <- cbind(term = rownames(stat), stat)
		WriteTable(stat, "differ.genes.stat.min_max.xls")
  
		if ( do.plot ) {
				top <- FindTopMarker(object.markers, top_num = 5, outfile = NULL)
				top_gene <- .change_name2(object, unique(top$gene))
				top_gene$GeneID <- ChangeOUTName(top_gene$GeneID, object@misc$fdata)
				WriteTable(top_gene, "Top.Heatmap.list.xls", col.names = FALSE)

				if ( length(unique(top$gene)) != 0 ){ 
						PlotHeatmapPlot(object, features = unique(top$gene), outfile = "Top.Heatmap.pdf")
						system( paste0( "convert -density 50 Top.Heatmap.pdf Top.Heatmap.png" ) )
				}
		} else {
				return(marker_list)
		}
}



#########################


#########################
SaveObj <- function(object, outfile = "obj.Rda", ...) {
		save(object, file = outfile, ...)
}

save_marker <- function(obj.markers, ...) {
		save(obj.markers, file = "markers.Rda", ...)
}

touch_success <- function(flag = NULL) {
		if ( !is.null(flag) && flag == "heatmap" ) {
				system("touch ./_complete_heatmap")
		} else {
				system("touch ./_complete")
		}
}
#########################



#########################

.change_name2 <- function(object, data, col = NULL ){
		if ( is.null(col) ){
				new.dt <- data.frame(GeneID = data, GeneName = FindFeaturesName(object, data, "name"))
		}else{
				tmp.dt <- data.frame(GeneID = data[[col]], GeneName = FindFeaturesName(object, data[[col]], "name"))
				nn <- which(colnames(data) == col)
				new.dt <- cbind(data[, 1 : nn, drop = F ], tmp.dt, data[, (nn + 1) : ncol(data), drop = F ])
				new.dt[[col]] <- NULL
		}
		return(new.dt)
}

.find_new_name <- function(x, new.list) {
		old.names <- strsplit(x = x, split = "_")
		new.names <- unlist(x = lapply(X = old.names, FUN = .rename, new.list = new.list))
#		names(new.names) <- x
		return(new.names)
}

.rename <- function(old.name, new.list) {
		old.sample <- paste(old.name[-length(old.name)], collapse = "_")
		if ( old.sample == '' ) return(old.name)
		new.sample <- new.list[[old.sample]]
		new.name   <- paste(c(new.sample, old.name[length(old.name)]), collapse = "_")
		return(new.name)
}

.default_heatmap <- function() {
  yaml <- '
data :
  use.genename : true
  scale.range : 2.5
## group.by : [ cluster | sample ]
  group.by : "cluster"
  filter :
    cluster : ~
    sample  : ~

hclust :
  col : false
  row : false
  show.tree :
    col : true
    row : true

font:
  family : Arial
  
axis : 
  x.text :
    show  : false
    size  : 10
    angle : 0
  y.text :
    show  : true
    size  : 10
    angle : 0
  group.text :
    show  : true
    size  : 15
    angle : 0
    ## loc : [ bottom | top ]
    loc : "bottom"

labs:
  title : ~
  xlab : ~
  ylab : ~
  size : 
    title : 18
    xlab  : 15
    ylab  : 15

legend :
## position : [ none | right | left | top | bottom ]
  position : "right"
  size : 
    title : 15
    label : 10
    width : 1.2
    height: 5
  color :
    high : "#FFFF00"
    mid  : "#000000"
    low  : "#FF00FF"
## direction : [ auto | h | v ]
  direction : "auto"
  title :
    name : ~
    ## loc : [ auto | left | right | top | bottom ]
    loc : "auto"
'
  pm <- yaml.load(yaml)
  return(pm)
}
#########################

########################
 
 



AddModuleScore <- function(
  object,
  features,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = NULL,
  name = 'Cluster',
  seed = 1,
  search = FALSE,
  ...
) {
## Note : 
##    alternate cut_number to Hmisc::cut2,
##    which will fail in some case
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  assay.old <- DefaultAssay(object = object)
  assay <- assay %||% assay.old
  DefaultAssay(object = object) <- assay
  assay.data <- GetAssayData(object = object)
  features.old <- features
  if (k) {
    .NotYetUsed(arg = 'k')
    features <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      features[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == i))
    }
    cluster.length <- length(x = features)
  } else {
    if (is.null(x = features)) {
      stop("Missing input feature list")
    }
    features <- lapply(
      X = features,
      FUN = function(x) {
        missing.features <- setdiff(x = x, y = rownames(x = object))
        if (length(x = missing.features) > 0) {
          warning(
            "The following features are not present in the object: ",
            paste(missing.features, collapse = ", "),
            ifelse(
              test = search,
              yes = ", attempting to find updated synonyms",
              no = ", not searching for symbol synonyms"
            ),
            call. = FALSE,
            immediate. = TRUE
          )
          if (search) {
            tryCatch(
              expr = {
                updated.features <- UpdateSymbolList(symbols = missing.features, ...)
                names(x = updated.features) <- missing.features
                for (miss in names(x = updated.features)) {
                  index <- which(x == miss)
                  x[index] <- updated.features[miss]
                }
              },
              error = function(...) {
                warning(
                  "Could not reach HGNC's gene names database",
                  call. = FALSE,
                  immediate. = TRUE
                )
              }
            )
            missing.features <- setdiff(x = x, y = rownames(x = object))
            if (length(x = missing.features) > 0) {
              warning(
                "The following features are still not present in the object: ",
                paste(missing.features, collapse = ", "),
                call. = FALSE,
                immediate. = TRUE
              )
            }
          }
        }
        return(intersect(x = x, y = rownames(x = object)))
      }
    )
    cluster.length <- length(x = features)
  }
  if (!all(LengthCheck(values = features))) {
    warning(paste(
      'Could not find enough features in the object from the following feature lists:',
      paste(names(x = which(x = !LengthCheck(values = features)))),
      'Attempting to match case...'
    ))
    features <- lapply(
      X = features.old,
      FUN = CaseMatch,
      match = rownames(x = object)
    )
  }
  if (!all(LengthCheck(values = features))) {
    stop(paste(
      'The following feature lists do not have enough features present in the object:',
      paste(names(x = which(x = !LengthCheck(values = features)))),
      'exiting...'
    ))
  }
  pool <- pool %||% rownames(x = object)
  data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
  data.avg <- data.avg[order(data.avg)]
#  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e30, n = nbin, labels = FALSE, right = FALSE)
  data.cut <- as.numeric(x = Hmisc::cut2(x = data.avg, m = round(x = length(x = data.avg) / (nbin + 1))))
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(x = sample(
          x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
          size = ctrl,
          replace = FALSE
        ))
      )
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(
    data = numeric(length = 1L),
    nrow = length(x = ctrl.use),
    ncol = ncol(x = object)
  )
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use, ])
  }
  features.scores <- matrix(
    data = numeric(length = 1L),
    nrow = cluster.length,
    ncol = ncol(x = object)
  )
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    data.use <- assay.data[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  features.scores.use <- features.scores - ctrl.scores
  rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  rownames(x = features.scores.use) <- colnames(x = object)
  object[[colnames(x = features.scores.use)]] <- features.scores.use
  CheckGC()
  DefaultAssay(object = object) <- assay.old
  return(object)
}
environment(AddModuleScore) <- asNamespace('Seurat')
assignInNamespace("AddModuleScore", AddModuleScore, ns = "Seurat")


CellSetFilter <- function(object, cell.set = NULL, col.name = "Clusters", is_renew_color = "color.cluster", new.name = NULL){
		cellset <- ReadCellSet(cell.set)
		if ( ! is.null(new.name) ) {
				cellset$cells <- new.name[cellset$cells]
		}
		if ( all(Cells(object) %in% cellset$cells) ) {
				object[[col.name]] <- cellset[Cells(object), "cluster"]
		} else {
				cell.use <- intersect(Cells(object), cellset$cells)
				if ( length(cell.use) == 0 ) 
						stop()
				object <- object[, cell.use]
				object[[col.name]] <- cellset[cell.use, "cluster"]
		}
		order <- names(sort(table(object[[col.name]]), decreasing = T))
		object@meta.data[[col.name]] <- factor(object@meta.data[[col.name]], levels = order)
		Idents(object) <- col.name

		if ( ! is.null(is_renew_color) ) {
				colors <- rainbow(length(levels(object@meta.data[[col.name]])))
				names(colors) <- levels(object@meta.data[[col.name]])
				object@misc[[is_renew_color]] <- colors
		}

		return(object)
}

ReadCellSet <- function(cell.set, header = FALSE, col.num = 1) {
		cellset <- NULL
		for ( name in names(cell.set) ){
				if ( ! file.exists(cell.set[[name]]) )
						stop()
				dt <- read.table(cell.set[[name]], header = header, sep = "\t")
				tmp <- data.frame(cells = dt[[col.num]]) %>% mutate(!!name := name)
				if ( is.null(cellset) ) {
						cellset <- tmp
				} else {
						cellset <- full_join(cellset, tmp, by = "cells")
				}
		}
		cellset <- tibble::column_to_rownames(cellset, "cells") %>%
				apply(1, function(x) {
						a <- na.omit(unlist(x))
						if (length(a) == 1) a else paste0("multi(", paste0(a, collapse = ","), ")")
				})
		cellset <- data.frame(cells = names(cellset), cluster = cellset, stringsAsFactors = FALSE)
		return(cellset)
}



### Dump
do_enrich <- function(parameter, do.run = TRUE){
		message("--> enrich : run <--")
		perl     <- "perl"
		script   <- paste0( bin, "/do_enrich.pl" )
		infile   <- paste0( parameter$outdir$task$differ, "/differ.genes.stat.xls" )
		indir    <- paste0( parameter$outdir$task$differ, "/enrich" )
		ref_pref <- parameter$infile$ref_pref
		flag     <- paste0( indir, "/_complete" )

		cmd <- paste( perl, script, infile, indir, ref_pref, sep = " ")
		if ( file.exists(flag) ) {
				message( "Enrich in ", indir, " has completed. Or delete ", flag, " first" )
		}else{
				print(cmd)
				error <- system(cmd)
				if (error != 0){
						stop("--> enrich : fail <--")
				}else{
						system(paste0("touch ", flag))
						message("--> enrich : DONE <--")
				}
		}
}

