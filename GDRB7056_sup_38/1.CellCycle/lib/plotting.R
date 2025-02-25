

VlnplotShell <- function(object, features = NULL, 
				group.by = "orig.ident", cols.use = NULL, 
				group.point.by = NULL, group.point.color = NULL,
				labs.y = NULL, titles = NULL, legend.position = "none",
				outfile = NULL, split.save.pref = NULL,
				pt.size = 0.1, alpha = 0.5, is.cline.x.text = TRUE,
				standard = NULL,
				...){
		if ( is.null(features) ) return(1)

		nCol <- ifelse(length(features) == 3, 3, ceiling(sqrt(length(features))))
		nRow <- ceiling(length(features) / nCol)

		plots <- VlnPlot(object = object, features = features, ncol = nCol, group.by = group.by, cols = cols.use, pt.size = pt.size, combine = FALSE, ...)
#		point.size.use = 0.1 ## Seurat 2.x, not supported anymore

		for ( i in seq(features) ) {
				if (!is.null(group.point.by) && pt.size > 0) {
						data <- cbind(plots[[i]]$data, object[[group.point.by]])
						plots[[i]]$layers[[2]] <- NULL
						plots[[i]] <- plots[[i]] + geom_jitter(aes_string(color = group.point.by), height = 0, size = pt.size, data = data)
						if ( !is.null(group.point.color) ) {
								plots[[i]] <- plots[[i]] + scale_color_manual(values = group.point.color)
						}
				}

				if ( group.by == "orig.ident" )
						plots[[i]] <- plots[[i]] + xlab("Samples")
				else if ( group.by == "Groups" )
						plots[[i]] <- plots[[i]] + xlab("Groups")
				else if ( group.by == "seurat_clusters" )
						plots[[i]] <- plots[[i]] + xlab("Clusters")

				if ( !is.null(labs.y[i]) && !is.na(labs.y[i]) )
						plots[[i]] <- plots[[i]] + ylab(labs.y[i])
				if ( !is.null(titles[i]) && !is.na(titles[i]) )
						plots[[i]] <- plots[[i]] + ggtitle(titles[i])
				# let point layer bottom
				plots[[i]]$layers[[1]]$aes_params$alpha <- alpha
				plots[[i]]$layers <- rev(plots[[i]]$layers)

				plots[[i]] <- plots[[i]] + box_theme_default()
				if ( is.cline.x.text ) 
						plots[[i]] <- plots[[i]] + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))

				if ( ! is.null(standard[[features[[i]]]]) ){ 
						plots[[i]] <- plots[[i]] + geom_hline(yintercept = standard[[features[[i]]]], color = "red")
				}
		}

		if ( !is.null(split.save.pref) ) {
				for ( i in seq(features) ){
						name <- ifelse( !is.null(names(features[i])), names(features[i]), features[i])
						outfile.tmp <- paste(split.save.pref, name, "pdf", sep = "." )
						ggsave(plots[[i]] + theme(legend.position = legend.position), file = outfile.tmp, width = 6, height = 6)
				}
		}

		plots.combined <- patchwork::wrap_plots(plots, ncol = nCol) & theme(legend.position = "none")
		if( !is.null(dev.list()) ) dev.off()

		if ( is.null(outfile) ){
				return(plots.combined)
		}else{
				ggsave(plots.combined, file = outfile, width = 6 * nCol, height = 6 * nRow, limitsize = F)
		}
}


FeatureScatterShell <- function(object, feature1 = NULL, feature2 = NULL, outfile = NULL, group.by = "orig.ident", cols = NULL, span = NULL, standard = NULL){
		p <- FeatureScatter(object = object, feature1 = feature1, feature2 = feature2, span = span, group.by = group.by, pt.size = 0.1, cols = cols) +
				guides(color = guide_legend(override.aes = list(size = 2), title = NULL))

		if ( ! is.null(standard) ) {
				if ( ! is.null(standard[[feature1]]) ){
						p <- p + geom_vline(xintercept = standard[[feature1]], color = "red")
				}

				if ( ! is.null(standard[[feature2]]) ){
						p <- p + geom_hline(yintercept = standard[[feature2]], color = "red")
				}
		}

		p <- p + dot_theme_default()
		ggsave(p, file = outfile, width = 7, height = 6)
}



PlotBasicStat <- function(object, outpref = NULL, color = NULL, assay = NULL, group.by = "orig.ident", span = NULL, stat_mito = TRUE, standard = NULL, ... ){
		if ( is.null(color) ) {
				color <- switch(group.by,
						"Groups" = object@misc[["color.group"]],
						"orig.ident" = object@misc[["color.sample"]],
						"seurat_clusters" = object@misc[["color.cluster"]])
		}

		if ( is.null(assay) ) {
				assay <- DefaultAssay(object)
		}

		nFeature <- paste0("nFeature_", assay)
		nCount   <- paste0("nCount_", assay)

		filter_feature <- function(features) unlist(sapply(features, function(x) if(x %in% colnames(object@meta.data)) x))
		get_y_name <- function(feature) sapply(features, switch, "percent.mito" = "Percentage(%)", "Number")

		features <- c(nFeature, nCount)
		if ( stat_mito )
				features <- c(features, "percent.mito")
		features <- filter_feature(features)
		VlnplotShell(object, features = features, labs.y = get_y_name(features), #titles = colnames(d),
						outfile = paste0(outpref, ".merge.pdf"), cols.use = color,
						group.by = group.by, assay = assay, standard = standard, ...)

		features <- c("expected.marker", "exclude.marker")
		features <- filter_feature(features)
		VlnplotShell(object, features = features, labs.y = get_y_name(features), #titles = colnames(d),
						outfile = paste0(outpref, ".PresetMarker.pdf"), cols.use = color,
						group.by = group.by, assay = assay, ...)

		FeatureScatterShell(object, feature1 = nCount, feature2 = nFeature, span = span,
						outfile = paste0(outpref, ".nCount-nFeature.pdf"), cols = color, group.by = group.by, standard = standard)
		if ( stat_mito && exists("percent.mito", object@meta.data) ) {
				FeatureScatterShell(object, feature1 = nCount, feature2 = "percent.mito", span = span,
						outfile = paste0(outpref, ".nCount-pMito.pdf"), cols = color, group.by = group.by, standard = standard)
		}
}


.PlotCluster <- function(object, reduction = NULL, cells = NULL, outfile = NULL,
				p1.group.by = "orig.ident",      p1.color = NULL, p1.label = FALSE,
				p2.group.by = "seurat_clusters", p2.color = NULL, p2.label = TRUE,
				plot.basic.size = 6, ...){
		if ( is.null(p1.color) && ! is.null(p1.group.by) ) {
				p1.color <- switch(p1.group.by,
						"Groups"          = object@misc[["color.group"]],
						"orig.ident"      = object@misc[["color.sample"]],
						"seurat_clusters" = object@misc[["color.cluster"]])
		}
		if ( is.null(p2.color) && ! is.null(p2.group.by) ) {
				p2.color <- switch(p2.group.by,
						"Groups"          = object@misc[["color.group"]],
						"orig.ident"      = object@misc[["color.sample"]],
						"seurat_clusters" = object@misc[["color.cluster"]])
		}
		p1 <- DimPlot(object, reduction = reduction, cells = cells, group.by = p1.group.by, cols = p1.color, label = p1.label, ...)
		p2 <- DimPlot(object, reduction = reduction, cells = cells, group.by = p2.group.by, cols = p2.color, label = p2.label, ...)
		p1 <- p1 + dot_theme_default() + ggtitle(NULL)
		p2 <- p2 + dot_theme_default() + ggtitle(NULL)
		if ( is.null(p2.group.by) ) {
				p <- p1
				width  <- plot.basic.size * 1.2
				height <- plot.basic.size
		} else if ( is.null(p1.group.by) ) {
				p <- p2
				width  <- plot.basic.size * 1.2
				height <- plot.basic.size
		} else {
				p <- p1 + p2
				width  <- plot.basic.size * 1.2 * 2
				height <- plot.basic.size
		}
		if ( is.null(outfile) ) {
				return(p)
		}else{
				ggsave(p, file = outfile, width = width, height = height, limitsize = FALSE )
		}
}


PlotCluster <- function(object, reduction = 'umap', p1.group.by = "orig.ident", split.by = p1.group.by, p2.group.by = "seurat_clusters", outpref = NULL, ...){
		.PlotCluster(object, reduction = reduction, p1.group.by = p1.group.by, p2.group.by = p2.group.by, outfile = paste0(outpref, ".pdf"), ...)
		for ( i in unique(object@meta.data[[split.by]]) ){
				cells.use <- rownames(object@meta.data)[object@meta.data[[split.by]] == i]
				.PlotCluster(object, reduction = reduction, cells = cells.use, p1.group.by = p1.group.by, p2.group.by = p2.group.by, outfile = paste0(outpref, ".", i, ".pdf"), ...)
		}
		data <- object[[reduction]]@cell.embeddings %>% as.data.frame() %>%
				tibble::rownames_to_column(var = "Cells") %>%
				left_join(.GetMetaData(object, cols = c("Samples" = p1.group.by, "Cluster" = p2.group.by, "Groups")))
		WriteTable(data, file = paste0(outpref, ".plot.data.tmp"))
}

.PlotClusterStat <- function(object, stat.what = "seurat_clusters", group.by = "orig.ident", color.st = NULL, color.gb = NULL, outpref = NULL, ...){
		if ( class(object) == "Seurat" ) {
				metadata <- object@meta.data
		} else {
				metadata <- object
		}
		if ( is.null(color.st) ) {
				if ( "misc" %in% slotNames(object) && exists(stat.what, object@misc) ) {
						color.st <- object@misc[[stat.what]]
				} else {
				color.st <- switch(stat.what,
						"seurat_clusters" = object@misc$color.cluster,
						"orig.ident" = object@misc$color.sample,
						"Groups" = object@misc$color.group
				)
				}
		}
		if ( is.null(color.gb) ) {
				if ( "misc" %in% slotNames(object) && exists(group.by, object@misc) ) {
						color.gb <- object@misc[[group.by]]
				} else {
				color.gb <- switch(group.by, 
						"seurat_clusters" = object@misc$color.cluster,
						"orig.ident" = object@misc$color.sample,
						"Groups" = object@misc$color.group
				)
				}
		}
		name.st <- switch(stat.what, "seurat_clusters" = "Cluster", "orig.ident" = "Samples", stat.what)
		name.gb <- switch(group.by,  "seurat_clusters" = "Cluster", "orig.ident" = "Samples", group.by)
		stat.what <- as.name(stat.what)
		group.by  <- as.name(group.by)

		stat_sample <- metadata %>%
				group_by(!! name.gb := !! group.by, !! name.st := !! stat.what) %>%
				summarise("Number of cells" = n())

		p <- list()
		p[["by"]] <- ggplot(stat_sample, aes_(x = as.name(name.gb), y = ~ `Number of cells`, fill = as.name(name.st)))
		p[["in"]] <- ggplot(stat_sample, aes_(x = as.name(name.st), y = ~ `Number of cells`, fill = as.name(name.gb)))
		if ( ! is.null(color.st) ) p[["by"]] <- p[["by"]] + scale_fill_manual(values = color.st)
		if ( ! is.null(color.gb) ) p[["in"]] <- p[["in"]] + scale_fill_manual(values = color.gb)
		geom_stack <- geom_bar(stat = "identity", position = 'stack')
		geom_fill  <- geom_bar(stat = "identity", position = "fill" )

		if ( is.null(outpref) ) {
				outpref <- paste0( name.st, ".stat")
		}
		for ( i in names(p) ) {
				p[[i]] <- p[[i]] + bar_theme_default()
				ggsave( p[[i]] + geom_stack, file = paste0( outpref, ".", i, name.gb, ".pdf"),     height = 6, width = 8 )
				ggsave( p[[i]] + geom_fill + ylab("Fraction of Cells"),  file = paste0( outpref, ".", i, name.gb, ".pct.pdf"), height = 6, width = 8 )
		}
}


.PlotFeaturePlot <- function(object, features, outfile = NULL, reduction = NULL, is.use.name = TRUE, color.high = "blue", color.low = "lightgrey", show.cluster.label = FALSE, nCol = NULL, plot.basic.size = 4, group.by = "seurat_clusters", is.combine = TRUE, cols = NULL, ... ) {
		if ( show.cluster.label ) Idents(object) <- group.by
		if ( is.null(cols) ) cols <- c(color.low, color.high)
		plots <- FeaturePlot(object, features = features, order = TRUE, reduction = reduction, combine = FALSE, label = show.cluster.label, cols = cols, ...)
		if ( is.use.name ) {
				name <- FindFeaturesName(object, features)
				for ( i in seq(plots) ){
						plots[[i]] <- plots[[i]] + ggtitle(name[i])
				}
				names(plots) <- name
		} else {
				names(plots) <- features
		}

		if ( is.null(nCol) ) nCol <- ceiling(sqrt(length(features)))
		nRow <- ceiling(length(features) / nCol)

		p <- wrap_plots(plots, ncol = nCol) & dot_theme_default()

		if ( is.null(outfile) ) {
				return(p)
		} else {
				ggsave(p, file = outfile, width = plot.basic.size * (6/5) * nCol, height = plot.basic.size * nRow, limitsize = FALSE)
		}
}

StatMarker <- function(object.markers, Cluster_name = "Cluster", Item_name = "Number of DE genes", color = NULL, outpref = "DeGene.stat"){
#		stat <- cbind( Cluster = Cluster_name, t(table(object.markers$cluster)) )
		stat <- cbind( tibble( !! Cluster_name := Item_name ) %>% as.matrix(), t(table(object.markers$cluster)) )
		WriteTable(stat, paste0(outpref, ".xls"))

		p <- ggplot(object.markers) + geom_bar(aes(x = cluster, fill = cluster), stat = "count") +
				theme_light() + labs(x = Cluster_name, y = Item_name)
		p <- p + bar_theme_default()
		p <- p + theme(legend.position = "none")
		if ( !is.null(color) ) p <- p + scale_fill_manual(values = color)
		ggsave(p, file = paste0(outpref, ".pdf"), height = 6, width = 8 )

		invisible(stat)
}



PlotDotPlot <- function(object, features = NULL, outfile = NULL, group.by = "seurat_clusters", is.use.name = TRUE, ...){
		p <- DotPlot(object, features = features, group.by = group.by, ...) + RotatedAxis()
		p <- p + dot_theme_default() +
				theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
				      axis.text.y = element_text(hjust = 1))
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


PlotHeatmapPlot <- function(object, features = NULL, group.by = "seurat_clusters", is.use.name = TRUE,
							outfile = NULL, group.colors = NULL, color = c("#FF00FF", "#000000", "#FFFF00")) {
		if ( is.null(group.colors) ) {
				group.colors <- switch(group.by,
								"seurat_clusters" = object@misc$color.cluster,
								"orig.ident" = object@misc$color.sample,
								"Groups" = object@misc$color.group
								)

		} else {
				group.colors <- group.colors[levels(object@meta.data[[group.by]])]
		}
		p <- DoHeatmap( object = object, features = features, cells = NULL,
						group.by = group.by, group.colors = group.colors,
						combine = FALSE, raster = FALSE)
		p <- p[[1]]
		p <- p + theme(legend.title = element_blank())
		p$layers[[2]] <- NULL

		if ( is.use.name ) {
				levels(p$data$Feature) <- FindFeaturesName(object, levels(p$data$Feature))
		}

		if (length(color) == 2) {
				p <- p + scale_fill_gradient(low = color[1], high = color[2], na.value = "white")
		} else if (length(color) == 3) {
				p <- p + scale_fill_gradient2(low = color[1], mid = color[2], high = color[3], na.value = "white")
		} else if (length(color) > 3) {
				p <- p + scale_fill_gradientn(colors = color, na.value = "white")
		}

#		p <- ggplotGrob(p)
#		for (i in grep('strip', p$layout$name)){
#				p$grobs[[i]]$layout$clip <- "off"
#		}

		if ( is.null(outfile) ) {
				return(p)
		} else {
				h <- max(7, length(unique(features)) * 0.11 + 2.5 )
				w <- h * 4 / 3
				p <- p + theme(plot.margin=margin(t=1,r=1,unit="lines"))
				ggsave(p, file = outfile, width = w, height = h, limitsize = FALSE )
		}
}

PlotAboutFeatures <- function(object, features = NULL, group.by = "seurat_clusters", outpref = NULL, group.colors = NULL, plot.feature = FALSE, is.use.name = TRUE, ... ) {
		object@meta.data <- droplevels(object@meta.data)
		PlotDotPlot(object, features = features, group.by = group.by, outfile = paste0(outpref, ".DotPlot.pdf"), is.use.name = is.use.name)
		PlotHeatmapPlot(object, features = features, group.by = group.by, outfile = paste0(outpref, ".Heatmap.pdf"), group.colors = group.colors, is.use.name = is.use.name, ...)
		if ( plot.feature ) {
				PlotFeaturePlot(object, features = features, outfile = paste0(outpref, ".Distribution.pdf"), reduction = "tsne", is.use.name = is.use.name)
		}
}

PlotPresetMarker <- function(object, cols.use = object@misc$color.cluster, group.by = "seurat_clusters", outpref = "PresetMarker"){
		PresetMarker <- union(object@misc[["expected.marker"]], object@misc[["more.marker"]])
		PresetMarker <- union(PresetMarker, object@misc[["exclude.marker"]])
		if ( length(PresetMarker) ) {
				VlnplotShell(object, features = PresetMarker, outfile = paste0(outpref, ".VlnPlot.pdf"), titles = object@misc$fdata[PresetMarker, "merge_name"], cols.use = cols.use, group.by = group.by)
				PlotAboutFeatures(object, features = PresetMarker, outpref = outpref, plot.feature = TRUE, group.by = group.by)
		}
}

PlotDensityPlot <- function(object, features = NULL, reduction = "umap", outpref = "DensityPlot/", ...){
		if ( is.null(features) ) {
				.PlotDensityPlot(object, reduction = reduction, outpref = outpref, ...)
		} else {
				for ( i in features ) {
						.PlotDensityPlot(object, i, reduction = reduction, outpref = outpref, ...)
				}
		}
}

PlotFeaturePlot <- function(object, features, reduction = "umap", is.combine = TRUE, outpref = NULL, outfile = NULL, is.use.name = TRUE, ...) {
		if ( is.combine ) {
				.PlotFeaturePlot(object, features = features, reduction = reduction, outfile = outfile, ...)
		}else{
				if ( is.null(outpref) ) outpref <- "ExpPlot"
				for ( i in features ) {
						if ( is.use.name ) {
								name <- FindFeaturesName(object, i)
								name <- gsub("[ /\\]", "_", name)
						} else {
								name = i
						}
						.PlotFeaturePlot(object, features = i, reduction = reduction, outfile = paste(c(outpref, name, "pdf"), collapse = "."), plot.basic.size = 6, ...)
				}
		}
}

PlotVlnPlot <- function(object, features, outpref = NULL, group.by = "seurat_clusters", cols.use = NULL, is.use.name = TRUE, ...) {
		name <- if ( is.use.name ) {
				FindFeaturesName(object, features)
		} else if ( is.null(names(features)) ) {
				features
		} else {
				names(features)
		}
		names(features) <- as.vector(name)
		if ( is.null(cols.use) ) {
				cols.use <- switch(group.by,
						seurat_clusters = object@misc$color.cluster,
						orig.ident = object@misc$color.sample)
#						NULL)
		}
		VlnplotShell(object, features, split.save.pref = outpref, group.by = group.by, titles = name, cols.use = cols.use, ...)
}




.PlotDensityPlot <- function(object, features = NULL, reduction = "umap", is.return = FALSE, outpref = NULL, is.consider.exp = FALSE, is.filter.noexp = TRUE) {
		dt <- as.data.frame(object[[reduction]]@cell.embeddings)
		if ( is.null(features) ) {
				if ( nrow(dt) < 2 ) return(NULL)
				dt$density <- KDE(x = dt[[1]], y = dt[[2]])
				outname <- paste0(outpref, "DensityPlot.pdf")
		}else{
				features <- features[1]
				exp <- t(as.data.frame(GetAssayData(object)[features,,drop = F]))
				dt <- cbind(dt, exp)
#				if ( is.filter.noexp ) {
						dt.filter <- dt[dt[[3]] > 0, ,drop = FALSE]
						if ( nrow(dt.filter) < 2 ) return(NULL)
						if ( is.consider.exp ) {
								dt.filter$density <- KDE(x = dt.filter[[1]], y = dt.filter[[2]], z = dt.filter[[3]])
						}else{
								dt.filter$density <- KDE(x = dt.filter[[1]], y = dt.filter[[2]])
						}
						dt <- full_join(dt, dt.filter)
#				} else {
#						dt$density <- KDE(x = dt[[1]], y = dt[[2]], z = dt[[3]])
#				}
				dt <- dt[order(-dt$density, na.last = F),,drop = FALSE]
				name <- FindFeaturesName(object, features)
				name <- gsub("[ /\\]", "_", name)
				outname <- paste0(outpref, "DensityPlot.", name, ".pdf")
		}
		p <- ggplot(dt, aes_string(x = colnames(dt)[1], y = colnames(dt)[2], color = "density")) +
			geom_point() + scale_color_viridis_c(option = "A", na.value = "grey90") #+ cowplot::theme_cowplot()
		if ( ! is.null(features) ) p <- p + ggtitle(name) + theme(plot.title = element_text(hjust = 0.5))
		
		p <- p + dot_theme_default()

		if ( is.return ) {
				return(list(p, dt))
		} else {
				ggsave(p, file = outname, width = 8, height = 7)
		}
}


kde3d <- function (x, y, z, h, n = 20, lims = c(range(x), range(y), range(z))) 
{
  nx <- length(x)
  if (length(y) != nx || length(z) != nx) 
    stop("data vectors must be the same length")
  if (missing(h)) 
    h <- c(MASS::bandwidth.nrd(x),
           MASS::bandwidth.nrd(y),
           MASS::bandwidth.nrd(z)) / 6
  else if (length(h) != 3)
    h <- rep(h, length = 3)
  if (length(n) != 3)
    n <- rep(n, length = 3)
  if (length(lims) == 2)
    lims <- rep(lims, length = 6)
  gx <- seq(lims[1], lims[2], length = n[1])
  gy <- seq(lims[3], lims[4], length = n[2])
  gz <- seq(lims[5], lims[6], length = n[3])
  mx <- matrix(outer(gx, x, dnorm, h[1]), n[1], nx)
  my <- matrix(outer(gy, y, dnorm, h[2]), n[2], nx)
  mz <- matrix(outer(gz, z, dnorm, h[3]), n[3], nx)
  v <- array(0, n)
  tmy.nx <- t(my) / nx
  for (k in 1:n[3]) {
    tmy.nz.zk <- tmy.nx * mz[k,] # uses recycling to scale the rows
    v[,,k] <- mx %*%  tmy.nz.zk
  }
  return(list(x = gx, y = gy, z = gz, d = v))
}


KDE <- function(x, y, z = NULL, h, n = 25, ...){
  hx <- MASS::bandwidth.nrd(x)
  hy <- MASS::bandwidth.nrd(y)
  if ( hx <= 0 ) hx <- diff(range(x)) / n
  if ( hy <= 0 ) hy <- diff(range(y)) / n
  if ( is.null(z) ) {
	if ( missing(h) ) h <- c(hx, hy)
    d <- MASS::kde2d(x = x, y = y, h = h, n = n, ...)
    gr <- data.frame(with(d, expand.grid(x,y)), as.vector(d$z))
    colnames(gr) <- c("xgr", "ygr", "zgr")
    mod <- loess(zgr~xgr*ygr, data=gr)
    dens <- predict(mod, newdata=data.frame(xgr = x, ygr = y))
  }else{
    hz <- MASS::bandwidth.nrd(z)
    if ( hz <= 0 ) hz <- diff(range(z)) / n
	if ( missing(h) ) h <- c(hx, hy, hz)
    d <- kde3d(x, y, z, h, n = n, ...)
    gr <- data.frame(with(d, expand.grid(x,y,z)), as.vector(d$d))
#    dens2 <- KDE(x = x, y = y, h = h, n = n, ...)
#    d <- MASS::kde2d(x = dens2, y = z, h = h, n = n, ...) 
#	gr <- data.frame(with(d, expand.grid(x,y,z)), as.vector(d$z))
    colnames(gr) <- c("xgr", "ygr", "zgr", "dgr")
    mod <- loess(dgr~xgr*ygr*zgr, data=gr)
    dens <- predict(mod, newdata=data.frame(xgr = x, ygr = y, zgr = z))
  }
  
  return(dens)
}

PlotHeatmap.online <- function(object, feature = NULL, parameter = list(), outfile = NULL,
		annot_with = parameter$annotation_col$with, annot_by = parameter$annotation_col$by,
		orders = parameter$annotation_col$order,
		filter_list = parameter$data$filter, assay = NULL, slot = "scale.data",		
		legend.position = IfNull(parameter$legend$position, "right"),
		color = IfNull(unlist(parameter$legend$color[c("low", "mid", "high")]), c("#A020F0", "#000000", "#FFFF00")),
		title = IfNull(parameter$labs$title, NA), xlab = parameter$labs$xlab, ylab = parameter$labs$ylab,
		scale.range = IfNull(parameter$data$scale.range, 3), scale.min = scale.range * -1, scale.max = scale.range,
		cluster_rows = IfNull(parameter$hclust$row, TRUE), cluster_cols = FALSE,
		show_rownames = IfNull(parameter$axis$y.text$show, TRUE), show_colnames = IfNull(parameter$axis$x.text$show, FALSE),
		is.useGeneName = IfNull(parameter$data$use.genename, TRUE),
		is.useIdName = IfNull(parameter$data$use.id_name, FALSE),
		fontsize = IfNull(parameter$font$size,10), font_family = parameter$font$family,
		...) {
		
		cells <- Cells(object)
		for ( i in names(filter_list) ) {
				colname <- switch(i, "cluster" = "seurat_clusters", "sample" = "orig.ident", "group" = "group")
				if ( ! is.null(filter_list[[i]]) && exists(colname, object@meta.data) ) {
						cells <- setdiff(cells, Cells(object)[object@meta.data[[colname]] %in% filter_list[[i]]])
				}
		}


		annotation_row <- NA
		annotation_col <- object@meta.data[cells, c("seurat_clusters", "orig.ident")]
		colnames(annotation_col) <- c("cluster", "sample")
		if ( exists("group", object@meta.data) )
				annotation_col$group <- object@meta.data[cells, "group"]
		if ( ! is.null(annot_with) )
				annotation_col <- annotation_col[intersect(annot_with, colnames(annotation_col))]
		order_col <- union(annot_by, annot_with)
		if ( ! is.null(orders) ) {
				for ( i in names(orders) ) {
						if ( i %in% colnames(annotation_col) ) {
								level <- levels(annotation_col[[i]])
								level <- c(intersect(orders[[i]], level), setdiff(level, orders[[i]]))
								annotation_col[[i]] <- factor(annotation_col[[i]], levels = level)
								if ( ! is.null(orders[[i]]) ) {
										annotation_col <- annotation_col[annotation_col[[i]] %in% orders[[i]],] ## show selected only
								}
						}
				}
		}
		annotation_col <- droplevels(annotation_col) ## show selected only
		annotation_col <- tibble::rownames_to_column(annotation_col, var = "cell") %>%
				arrange_at(vars(!!!lapply(order_col, as.name))) %>%
				tibble::column_to_rownames(var = "cell")

		annotation_colors <- sapply(colnames(annotation_col),
				function(x) switch(x,
						"sample"  = object@misc[["color.sample"]],
						"cluster" = object@misc[["color.cluster"]],
						"group"   = object@misc[["color.group"]]
				)[levels(annotation_col[[x]])],
				simplify = FALSE)
		annotation_colors <- annotation_colors[!sapply(annotation_colors, is.null)]

		is.legend <- if (legend.position == "none") FALSE else TRUE

		color <- colorRampPalette(color)(100)

		mat <- GetAssayData(object, assay = assay, slot = slot)
		if ( slot == "scale.data" && nrow(mat) == 0 ) {
				object <- ScaleData(object, assay = assay, features = rownames(object))
				mat <- GetAssayData(object, assay = assay, slot = slot)
		}
		mat <- mat[, rownames(annotation_col), drop = FALSE]
		if ( ! is.null(feature) ) {
				if ( ! is.null(object@misc$fdata) ) {
						feature <- FindFeaturesID(object, feature)
				}
				if ( ! all(feature %in% rownames(mat)) ) {
					object <- ScaleData(object, assay = assay, features = feature)
				}
				mat <- mat[as.character(feature),,drop = FALSE]
				if ( nrow(mat) < 2 )
						cluster_rows <- FALSE
		}
		mat <- MinMax(mat, scale.min, scale.max)
		breaks <- if ( length(unique(as.vector(mat))) == 1 ) {
				unique(as.vector(mat)) + c(-1, 0, 1)
		} else {
				NA
		}
		if ( is.useIdName ) {
				name <- FindFeaturesName(object, rownames(mat), col = "name")
				name <- paste0(rownames(mat), "(", name, ")" )
				rownames(mat) <- name
		} else {
		if ( is.useGeneName ) 
				rownames(mat) <- FindFeaturesName(object, rownames(mat), "name")
		}

		require(grid)
		ph <- pheatmap::pheatmap(mat, useRaster = T, border_color = NA, breaks = breaks,
						color = color, legend = is.legend,
						clustering_method = "ward.D2", treeheight_row = 20,
						cluster_rows = cluster_rows, cluster_cols = cluster_cols,
						annotation_row = annotation_row, annotation_col = annotation_col,
#						annotation_names_row = FALSE, annotation_names_col = FALSE,
						annotation_colors = annotation_colors,
						show_rownames = show_rownames, show_colnames = show_colnames,
						main = title, fontsize = fontsize, slient = TRUE)
		height <- max(7, convertHeight(sum(ph$gtable$heights), "inches", valueOnly = T))
		width  <- max(7, convertWidth( sum(ph$gtable$widths),  "inches", valueOnly = T))
		while(!is.null(dev.list())) { 
				dev.off()
		}

		if ( !is.null(outfile) ) {
				pdf(file = outfile, height = height, width = width)
		}
		if ( ! is.null(xlab) || ! is.null(ylab) ) {
				setHook("grid.newpage", 
						function() pushViewport(viewport(x = 1, y = 1, width = 0.9, height = 0.9, name = "vp", just = c("right","top"))),
						action = "prepend")
		}
		grid.draw(ph$gtable)
		if ( ! is.null(xlab) || ! is.null(ylab) ) {
				setHook("grid.newpage", NULL, "replace")
				grid.text(xlab, y = -0.05, gp = gpar(fontsize = fontsize))
				grid.text(ylab, x = -0.05, rot = 90, gp = gpar(fontsize = fontsize))
		}

		if ( !is.null(outfile) ) {
				dev.off()
		}
}

