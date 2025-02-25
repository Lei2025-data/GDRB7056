
PlotDeGeneVolcano <- function(marker, parameter = list(),
				log2fc = IfNull(parameter$FindMarkers$log2fc, 1),
				pvalue = IfNull(parameter$FindMarkers$pvalue, 0.05),
				use.q  = IfNull(parameter$FindMarkers$use.q, FALSE),
				show_genes = NULL
				){
		plist <- list()
		flist <- list(order = names(marker))
		for ( i in names(marker) ) {
				marker[[i]] <- marker[[i]] %>% filter(cells.1 * cells.2 != 0)
				marker[[i]]$cluster <- droplevels(marker[[i]]$cluster)
				flist[[i]] <- levels(marker[[i]]$cluster)
				for ( j in levels(marker[[i]]$cluster) ) {
						name <- paste0(i, "_", j)
						marker.tmp <- subset(marker[[i]], cluster == j)
						if ( ! is.null(show_genes[[i]][[j]]) ) {
								genes <- show_genes[[i]][[j]]
								marker.tmp$show_gene <- NA
								marker.tmp$show_gene[marker.tmp$gene %in% genes] <- if(!is.null(names(genes))) names(genes) else genes
						}
						yc  <- ifelse(use.q, "p_val_adj", "p_val")
						y   <- paste0("-log10(", yc, ")")
						ylb <- ifelse(use.q, expression(-log[10](P.adjust)), expression(-log[10](Pvalue)) )
						xc  <- "avg_logFC"
						marker.tmp[[yc]][marker.tmp[[yc]] < 10^-50] <- 10^-50
						marker.tmp[[yc]][marker.tmp[[yc]] > 1] <- 1
						marker.tmp[[xc]][marker.tmp[[xc]] > 10] <- 10
						marker.tmp[[xc]][marker.tmp[[xc]] < -10] <- -10
						plist[[name]] <- ggplot(marker.tmp, aes_string(x = xc, y = y)) +
								geom_point(aes(color = sig)) +
								scale_color_manual(values = c("up" = "#F8766D", "nosig" = "grey60", "down" = "#619CFF"), name = NULL) +
								geom_hline(yintercept = -log10(pvalue),   linetype = "dashed", color = "grey40") +
								geom_vline(xintercept = c(-1,1) * log2fc, linetype = "dashed", color = "grey40") +
								labs(x = expression(log[2](FC)), y = ylb) +
								xlim(-10,10) + ylim(0,50) + 
								theme_light() + 
								theme(legend.key = element_rect(color = "grey70")) 
						if ( ! is.null(show_genes[[i]]) ) {
								plist[[name]] <- plist[[name]] + ggrepel::geom_text_repel(aes(label = show_gene))
						}
						plist[[name]] <- plist[[name]] + dot_theme_default()
						name <- gsub("/", "_", name)
						ggsave(plist[[name]], file = paste0("Volcano.", name, ".pdf"), width = 7, height = 6)
				}
		}
		yaml::write_yaml(flist, file = "Volcano.list", indent.mapping.sequence = TRUE)
		invisible(plist)
}





PlotContrastHeatmap <- function(object, object.marker, group.data = NULL,
				top.num = 5, cluster.name = "seurat_clusters", do.return = FALSE,
				is.sort.by.contrast = TRUE,
				assay = NULL, slot = "scale.data", 
				scale.min = -2.5, scale.max = 2.5, 
				color = c("purple", "black", "yellow")
				) {
		for ( contrast_name in names(object.marker) ) {
				contrast <- strsplit(contrast_name, "-vs-")[[1]]

				annotation <- GetAnnotation(group.data, contrast, meta.data = object@meta.data, add.group = cluster.name, is.sort.by.contrast = is.sort.by.contrast)
				cells <- rownames(annotation)
				features <- object.marker[[contrast_name]] %>% group_by(cluster) %>% top_n(top.num, abs(avg_logFC))
				features <- unique(features$gene)

				mat <- GetAssayData(object, assay = assay, slot = slot)[features, cells]
				mat <- MinMax(mat, scale.min, scale.max)
				rownames(mat) <- FindFeaturesName(object, rownames(mat), "name")
				if ( do.return ){
				} else {
						annotation_colors <- list(switch(cluster.name,
								"Groups" = object@misc[["color.group"]],
								"orig.ident" = object@misc[["color.sample"]],
								"seurat_clusters" = object@misc[["color.cluster"]]))
						names(annotation_colors) <- switch(cluster.name, "seurat_clusters" = "cluster", "orig.ident" = "samples", cluster.name)
						if ( is.null(annotation_colors[[1]]) ) annotation_colors[[1]] <- SetColor(unique(annotation[[names(annotation_colors)]]))

						color <- colorRampPalette(color)(100)
						MAX <- max(abs(mat))
						breaks <- seq(-1 * MAX, MAX, length.out = length(color))
						pheatmap::pheatmap(mat, cluster_cols = F, show_colnames = F,
								annotation_col = annotation, annotation_colors = annotation_colors,
								scale = "none", color = color, breaks = breaks,
								cellheight = 10, #cellwidth = 10,
								filename = paste0("Heatmap.", contrast_name, ".pdf"))
				}

		}
}

GetAnnotation <- function(group.data, contrast, meta.data = NULL, add.group = "seurat_clusters", is.sort.by.contrast = TRUE, keep.group = FALSE) {
		if ( is.null(group.data) ) return(NULL)
		cell.1 <- group.data[[contrast[1]]]
		cell.2 <- group.data[[contrast[2]]]
		data <- data.frame(cell = rownames(group.data), row.names = rownames(group.data))
		if ( sum(cell.1&cell.2) == 0 ){
				data$contrast <- contrast[1]
				data$contrast[cell.2] <- contrast[2]
				data$contrast <- factor(data$contrast, levels = c(contrast[1], contrast[2]))
				data$group <- data$contrast
		} else {
				data[[contrast[1]]] <- sapply(cell.1, function(x) ifelse(x, contrast[1], NA))
				data[[contrast[2]]] <- sapply(cell.2, function(x) ifelse(x, contrast[2], NA))
				data$group <- "both"
				data$group[is.na(data[[contrast[1]]])] <- contrast[2]
				data$group[is.na(data[[contrast[2]]])] <- contrast[1]
				data$group <- factor(data$group, levels = c(contrast[1], "both", contrast[2]))
		}
		data <- data[cell.1|cell.2, , drop = F]
		if ( ! is.null(meta.data) && ! is.null(add.group) ) {
				name <- switch(add.group, "seurat_clusters" = "cluster", "orig.ident" = "samples", add.group)
				data[[name]] <- meta.data[rownames(data), add.group]
				if ( is.sort.by.contrast ) {
						data <- data %>% arrange(group, !! as.name(name))
				} else {
						data <- data %>% arrange(!! as.name(name), group)
				}
		} else {
				data <- data %>% arrange(group)
		}
		if (  packageVersion("dplyr") >= as.numeric_version("1.0.3") ) { ## damn it
				data <- data %>% select(-cell)
		} else {
				data <- tibble::column_to_rownames(data, var = "cell")
		}
		if ( ! keep.group ) {
				data$group <- NULL
		}
		return(data)
}


PlotContrastFeature <- function(object, object.marker, group.data, reduction = NULL, cluster.name = "seurat_clusters",
				top.num = 5, use.name = TRUE, color = NULL,
				legend = c("right", "bottom", "none"), is.combind = TRUE, do.return = FALSE, is.demo = FALSE, plot_size = 4,
				outpref = "FeaturePlot") {
		plots <- list()
		for ( contrast_name in names(object.marker) ) {
				contrast <- strsplit(contrast_name, "-vs-")[[1]]

				annotation <- GetAnnotation(group.data, contrast, meta.data = object@meta.data, add.group = cluster.name, keep.group = T)
				cells <- rownames(annotation)
				sub.object <- object[, cells]
				sub.object@meta.data$split.by <- annotation[colnames(sub.object),"group"]
				if ( is.demo ) {
						features <- as.character(object.marker[[contrast_name]]$gene[1])
				}else{
						features <- object.marker[[contrast_name]] %>% group_by(cluster) %>% top_n(top.num, abs(avg_logFC))
						features <- unique(features$gene)
				}


				data <- GetFeaturePlotData(sub.object, features = features, split.by = "split.by", reduction = reduction, max.cutoff = 3)
				dt <- reshape2::melt(data, measure.vars = features, variable.name = "features", value.name = "value")
				dt <- dt[order(dt[, "value"]), , drop = F]
				if ( use.name ) {
						new.names <- FindFeaturesName(object, features, "name")
				}else{
						new.names <- features
						names(new.names) <- features
				}
				my_color <- if ( is.null(color) ) {
								scale_color_viridis_c(NULL, option = "D")
						} else {
								if (length(x = color) == 1 && (is.numeric(x = color) || color %in% rownames(x = RColorBrewer::brewer.pal.info))) {
										scale_color_brewer(palette = color, na.value = "grey50")
#								} else if (length(x = color) == 1 && (color %in% c("alphabet", "alphabet2", "glasbey", "polychrome", "stepped"))) {
#										colors <- DiscretePalette(length(unique(data[[col.by]])), palette = color)
#										scale_color_manual(values = colors, na.value = "grey50")
								} else {
#										scale_color_manual(values = color, na.value = "grey50")
										scale_color_gradientn(colors = color, guide = "colorbar")
								}
						}
				ps <- lapply(levels(dt$features), function(i) {
								p <- ggplot(dt %>% filter(features == i), aes_string(x = colnames(dt)[1], y = colnames(dt)[2])) +
								      geom_point(aes(color = value), size = Seurat:::AutoPointSize(data = data)) +
									  facet_grid(features ~ split, labeller = labeller(features = as_labeller(new.names))) +
									  my_color +
									  dot_theme_default() +
									  theme(strip.background = element_blank(),
											strip.text = element_text(face = "bold", size = 15))
								p <- ggplotGrob(p)
								for (i in grep('strip', p$layout$name)){
										p$grobs[[i]]$layout$clip <- "off"
								}
								p <- cowplot::ggdraw(p)
								return(p)
							})
				if ( is.combind ) {
						ncol <- min(5, max(1, floor(sqrt(length(ps) * length(levels(data$split))) / length(levels(data$split))) ))
						nrow <- ceiling(length(ps) / ncol)
						ps <- cowplot::plot_grid(plotlist = ps, ncol = ncol, align = "hv", axis = "tblr")
						if ( ! do.return ) {
								if ( is.demo ) {
										outfile <- paste(outpref, contrast_name, reduction, "demo.pdf", sep = ".")
								} else {
										outfile <- paste(outpref, contrast_name, reduction, "pdf", sep = ".")
								}
								ggsave(ps, file = outfile, width = ncol * length(levels(data$split)) * plot_size, height = nrow * plot_size, limitsize = FALSE)
						}
				} else {
						if ( ! do.return ) {
								for ( feature in names(ps) ) {
										outfile <- paste(outpref, contrast_name, reduction, feature, "pdf", sep = ".")
										ggsave(ps, file = outfile, width = plot_size, height = plot_size, limitsize = FALSE)
								}
						}
				}
				plots[[contrast_name]] <- ps
		}
		if ( do.return ) {
				return(plots)
		}
}


PlotContrastDotPlot <- function(object, object.marker, group.data,
				color = NULL, use.name = TRUE,
				top.num = 5, cluster.name = "seurat_clusters", do.return = FALSE
				){
		if ( ! is.null(color) ){
				if ( length(color) < 2 ) stop("color number must be 2.")
				color <- color[1:2]
		} else {
				color <- c("blue", "red")
		}

		plots <- list()
		for ( contrast_name in names(object.marker) ) {
				contrast <- strsplit(contrast_name, "-vs-")[[1]]

				annotation <- GetAnnotation(group.data, contrast, object@meta.data, cluster.name, keep.group = T)
				cells <- rownames(annotation)
				Idents(object) <- cluster.name
				sub.object <- object[, cells]
				sub.object@meta.data$split.by <- annotation[colnames(sub.object),"group"]
				features <- object.marker[[contrast_name]] %>% group_by(cluster) %>% top_n(top.num, abs(avg_logFC))
				features <- unique(features$gene)

				names(color) <- contrast
				if ( "both" %in% levels(sub.object@meta.data$split.by) )
						color["both"] <- colorspace::hex(colorspace::mixcolor(alpha = 0.5,
												color1 = colorspace::RGB(t(col2rgb(color[1])/255)),
												color2 = colorspace::RGB(t(col2rgb(color[2])/255))))
				p <- DotPlot(sub.object, features = features, split.by = "split.by", cols = color) +
						dot_theme_default() +
						theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
							axis.text.y = element_text(hjust = 1),
							legend.title = element_text(hjust = 0))
				if ( use.name ) {
						levels(p$data$features.plot) <- FindFeaturesName(object, levels(p$data$features.plot), "name") 
				}
				if ( ! do.return ) {
						w <- max(6, ceiling(length(features)) * 0.35) #+ 2)
						h <- max(6, length(levels(p$data$id)) * 0.4)
#						ggsave(p, file = paste0("DotPlot.", contrast_name, ".pdf"), width = w, height = h, limitsize = FALSE)
						SavePlot(p, filename = paste0("DotPlot.", contrast_name, ".pdf"), plot.width = w, plot.height = h)
				}
				plots[[contrast_name]] <- p
		}
		if ( do.return ) {
				return(plots)
		}
}


GetFeaturePlotData <- function (object, features, dims = c(1, 2), cells = NULL, cols = if (blend) {
    c("lightgrey", "#ff0000", "#00ff00")
} else {
    c("lightgrey", "blue")
}, pt.size = NULL, order = FALSE, min.cutoff = NA, max.cutoff = NA, 
    reduction = NULL, split.by = NULL, shape.by = NULL, slot = "data", 
    blend = FALSE, blend.threshold = 0.5, label = FALSE, label.size = 4, 
    repel = FALSE, ncol = NULL, combine = TRUE, coord.fixed = FALSE, 
    by.col = TRUE, sort.cell = FALSE) 
{
	## copy from FeaturePlot and delete the part of plot data
    reduction <- reduction %||% DefaultDimReduc(object = object)
    if (length(x = dims) != 2 || !is.numeric(x = dims)) {
        stop("'dims' must be a two-length integer vector")
    }
    dims <- paste0(Key(object = object[[reduction]]), dims)
    cells <- cells %||% colnames(x = object)
    data <- FetchData(object = object, vars = c(dims, "ident", 
        features), cells = cells, slot = slot)
    if (ncol(x = data) < 4) {
        stop("None of the requested features were found: ", paste(features, 
            collapse = ", "), " in slot ", slot, call. = FALSE)
    }
    else if (!all(dims %in% colnames(x = data))) {
        stop("The dimensions requested were not found", call. = FALSE)
    }
    features <- colnames(x = data)[4:ncol(x = data)]
    min.cutoff <- mapply(FUN = function(cutoff, feature) {
        return(ifelse(test = is.na(x = cutoff), yes = min(data[, 
            feature]), no = cutoff))
    }, cutoff = min.cutoff, feature = features)
    max.cutoff <- mapply(FUN = function(cutoff, feature) {
        return(ifelse(test = is.na(x = cutoff), yes = max(data[, 
            feature]), no = cutoff))
    }, cutoff = max.cutoff, feature = features)
    check.lengths <- unique(x = vapply(X = list(features, min.cutoff, 
        max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
    if (length(x = check.lengths) != 1) {
        stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols, 
        ]$maxcolors, no = length(x = cols))
    data[, 4:ncol(x = data)] <- sapply(X = 4:ncol(x = data), 
        FUN = function(index) {
            data.feature <- as.vector(x = data[, index])
            min.use <- SetQuantile(cutoff = min.cutoff[index - 
                3], data.feature)
            max.use <- SetQuantile(cutoff = max.cutoff[index - 
                3], data.feature)
            data.feature[data.feature < min.use] <- min.use
            data.feature[data.feature > max.use] <- max.use
            if (brewer.gran == 2) {
                return(data.feature)
            }
            data.cut <- if (all(data.feature == 0)) {
                0
            }
            else {
                as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature), 
                  breaks = brewer.gran)))
            }
            return(data.cut)
        })
    colnames(x = data)[4:ncol(x = data)] <- features
    rownames(x = data) <- cells
    data$split <- if (is.null(x = split.by)) {
        RandomName()
    }
    else {
        switch(EXPR = split.by, ident = Idents(object = object)[cells], 
            object[[split.by, drop = TRUE]][cells])
    }
    if (!is.factor(x = data$split)) {
        data$split <- factor(x = data$split)
    }
    if (!is.null(x = shape.by)) {
        data[, shape.by] <- object[[shape.by, drop = TRUE]]
    }
	return(data)
}
environment(GetFeaturePlotData) <- asNamespace('Seurat')


