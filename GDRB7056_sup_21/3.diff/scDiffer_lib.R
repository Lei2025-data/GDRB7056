FindGroupIndex <- function(object, parameter, slot = "counts", assay = NULL, sample_col = "orig.ident"){
#		metadata <- object@meta.data[, sample_col, drop = F]
		metadata <- object@meta.data %>% select(sample = !! sample_col)
		for ( i in unique(metadata$sample)) {
				metadata <- metadata %>% mutate(!! i := sample %in% i )
		}
		for ( i in names(parameter$group$by_sample) ){
				metadata <- metadata %>% mutate(!! i := sample %in% parameter$group$by_sample[[i]])
		}
		metadata <- metadata %>% select(- sample)
		for ( i in names(parameter$group$by_gene) ) {
				gene_id <- parameter$group$by_gene[[i]][[1]]
				gene_id <- FindFeaturesID(object, gene_id)
				gene_thres <- parameter$group$by_gene[[i]][[2]]
				data <- GetAssayData(object, slot = slot, assay = assay)
				metadata <- metadata %>% mutate(!! i := data[gene_id, ] >= gene_thres[[1]] & data[gene_id, ] <= gene_thres[[2]])
		}
		rownames(metadata) <- colnames(object)
		return(metadata)
}

FindGroupMarkers <- function (object, group.data, contrast = NULL,
		return.thresh = 0.01, use.qvalue = FALSE,
		assay = NULL, features = NULL, logfc.threshold = 0.25, 
		test.use = "MAST", slot = "data", min.pct = 0.1, min.diff.pct = -Inf, 
		verbose = TRUE, only.pos = FALSE, max.cells.per.ident = Inf, 
		random.seed = 1, latent.vars = NULL, min.cells.feature = 3, 
		min.cells.group = 3, pseudocount.use = 1, ...) {
		
		if ((test.use == "roc") && (return.thresh == 0.01)) {
				return.thresh <- 0.7
		}

		cells.1 <- rownames(group.data)[group.data[[contrast[[1]]]]]
		cells.2 <- rownames(group.data)[group.data[[contrast[[2]]]]]
		if ( length(cells.1) == 0 || length(cells.2) == 0 ) {
				message( "===> cell num : \n",
						"        ", contrast[[1]], " : ", length(cells.1), "\n",
						"        ", contrast[[2]], " : ", length(cells.2), "\n",
						"        [SKIP]")
				return(NULL)
		}

		data.1 <- CalAvgExp(object[, cells.1], is.return = TRUE, is.bulk = TRUE) %>%
				reshape2::melt(varnames = c("gene", "cluster"), value.name = "mean.1")
		data.2 <- CalAvgExp(object[, cells.2], is.return = TRUE, is.bulk = TRUE) %>%
				reshape2::melt(varnames = c("gene", "cluster"), value.name = "mean.2")
		pct.1 <- CalPctExp(object[, cells.1], is.return = TRUE, is.bulk = TRUE) %>% round(digits = 3) %>%
				reshape2::melt(varnames = c("gene", "cluster"), value.name = "pct.1")
		pct.2 <- CalPctExp(object[, cells.2], is.return = TRUE, is.bulk = TRUE) %>% round(digits = 3) %>%
		        reshape2::melt(varnames = c("gene", "cluster"), value.name = "pct.2")				
#		cells.num <- rbind( bulk = c(length(cells.1), length(cells.2)),
#							cbind(cells.1 = table(Idents(object[,cells.1])),
#								  cells.2 = table(Idents(object[,cells.2])))) %>%
#					 as.data.frame() %>% tibble::rownames_to_column(var = "cluster")
		c1 <- c(bulk = length(cells.1), table(Idents(object[,cells.1])))
		c2 <- c(bulk = length(cells.2), table(Idents(object[,cells.2])))
		data <- full_join(x = data.1, y = data.2, by = c("gene", "cluster")) %>%
				left_join(y = pct.1, by = c("gene", "cluster") ) %>%
				left_join(y = pct.2, by = c("gene", "cluster") ) %>%
#				left_join(y = cells.num, by = "cluster") %>%
				mutate(cells.1 = c1[cluster], cells.2 = c2[cluster]) %>%
				mutate(cluster = factor(cluster, levels = c("bulk", levels(object)))) %>%
				mutate(avg_logFC = log2(mean.2/mean.1), p_val = 1, p_val_adj = 1, sig = "nosig") # switch order from T-vs-C to C-vs-T
#				mutate(avg_logFC = log2(mean.1/mean.2), p_val = 1, p_val_adj = 1, sig = "nosig") 
				
		data[is.na(data)] <- 0
		rm(data.1, data.2, pct.1, pct.2, c1, c2)

		idents.all <- c("bulk", levels(droplevels(object)))

		if ( is.null(features) ) {
				features.use <- data %>% filter(mean.1 * mean.2 != 0 &
												(pct.1 > min.pct | pct.2 > min.pct) &
												abs(avg_logFC) > logfc.threshold / log(2))
				features.use <- split(features.use, features.use$cluster)
				features.use <- lapply(features.use, function(x) as.character(x$gene))
				logfc.threshold <- 0
				pseudocount.use <- 0
		}else{
				features.use <- lapply(idents.all, function(x) features)
		}

		genes.de <- list()
		messages <- list()
		for (i in seq(idents.all)) {
				if (verbose) message("Calculating cluster ", idents.all[i])


				cells.idents <- if ( idents.all[[i]] == "bulk" ){ 
						WhichCells(object)
				} else {
						WhichCells(object, idents = idents.all[[i]])
				}
				cells.1.tmp <- intersect(cells.2, cells.idents) # switch order from T-vs-C to C-vs-T
				cells.2.tmp <- intersect(cells.1, cells.idents) # switch order from T-vs-C to C-vs-T
				genes.de[[i]] <- tryCatch(expr = {
								FindMarkers(object = object, ident.1 = cells.1.tmp, ident.2 = cells.2.tmp,
											assay = assay,
											features = features.use[[idents.all[[i]]]], logfc.threshold = logfc.threshold, 
											test.use = test.use, slot = slot, min.pct = min.pct, 
											min.diff.pct = min.diff.pct, verbose = verbose, 
											only.pos = only.pos, max.cells.per.ident = max.cells.per.ident, 
											random.seed = random.seed, latent.vars = latent.vars, 
											min.cells.feature = min.cells.feature, min.cells.group = min.cells.group, 
											pseudocount.use = pseudocount.use, ...)
								}, error = function(cond) { return(cond$message) })
				if (class(x = genes.de[[i]]) == "character") {
						messages[[i]] <- genes.de[[i]]
						genes.de[[i]] <- NULL
				}
		}

		gde.all <- data.frame()
		for (i in seq(idents.all)) {
				if (is.null(x = unlist(x = genes.de[i]))) {
						next
				}
				gde <- genes.de[[i]]
				if (nrow(x = gde) > 0) {
						if (test.use == "roc") {
								gde <- subset(x = gde, subset = (myAUC > return.thresh | myAUC < (1 - return.thresh)))
#						} else if (is.null(x = node) || test.use %in% c("bimod", "t")) {
						} else {
								gde <- gde[order(gde$p_val_adj, gde$p_val, -gde[, 2]), ]
								if ( use.qvalue ) {
										gde <- subset(x = gde, subset = p_val_adj < return.thresh)
								} else {
										gde <- subset(x = gde, subset = p_val < return.thresh)
								}
						}
				
#						if (nrow(x = gde) > 0) {
								genes.de[[i]]$sig <- "nosig"
								genes.de[[i]][rownames(gde)[gde$avg_logFC > 0], "sig"] <- "up"
								genes.de[[i]][rownames(gde)[gde$avg_logFC < 0], "sig"] <- "down"	
								genes.de[[i]]$cluster <- idents.all[i]
								genes.de[[i]]$gene <- rownames(x = genes.de[[i]])
								gde.all <- rbind(gde.all, genes.de[[i]])
#						}
				}
		}
		if ( nrow(x = gde.all) == 0) {
				warning("No DE genes identified", call. = FALSE, immediate. = TRUE)
				return(NULL)
		}
		gde.all <- gde.all %>%
				select(cluster, gene, p_val, p_val_adj, sig) %>%
				full_join(y = data, by = c("gene", "cluster")) %>%
				mutate( p_val     = if_else(is.na(p_val.x),     p_val.y,     p_val.x),
						p_val_adj = if_else(is.na(p_val_adj.x), p_val_adj.y, p_val_adj.x),
						sig       = if_else(is.na(sig.x),       sig.y,       sig.x)) %>%						
				mutate( sig     = factor(sig, levels = c("up", "nosig", "down")),
						cluster = factor(cluster, levels = levels(data$cluster)),
						contrast  = factor(paste0(contrast, collapse = "-vs-")) ) %>%
				select(-p_val.x, -p_val.y, -p_val_adj.x, -p_val_adj.y, -sig.x, -sig.y)

		if ((only.pos) && nrow(x = gde.all) > 0) {
				return(subset(x = gde.all, subset = gde.all[, "avg_logFC"] > 0))
		}
		if (nrow(x = gde.all) == 0) {
				warning("No DE genes identified", call. = FALSE, immediate. = TRUE)
		}
		if (length(x = messages) > 0) {
				warning("The following tests were not performed: ", call. = FALSE, immediate. = TRUE)
				for (i in 1:length(x = messages)) {
						if (!is.null(x = messages[[i]])) {
								warning("When testing cluster ", idents.all[i], " :\n\t", messages[[i]], call. = FALSE, immediate. = TRUE)
						}
				}
		}
		return(gde.all)
#		return(NULL)
}

DoFindGroupMarkers <- function(object, group.data, parameter, ...){
		if ( future::nbrOfWorkers() > 1 ) {
				marker <- future.apply::future_lapply(parameter$differ, function(i)
								FindGroupMarkers(object, group.data, contrast = i,
										return.thresh = parameter$FindMarkers$pvalue,
										use.qvalue = parameter$FindMarkers$use.q,
										logfc.threshold = parameter$FindMarkers$log2fc / log2(exp(1)),
										test.use = parameter$FindMarkers$method,
										min.pct = parameter$FindMarkers$min_pct, ...
										)
						)
				names(marker) <- sapply(parameter$differ, paste, collapse = "-vs-")
				for ( i in seq(marker) ) {
						if ( is.null(marker[[i]]) )  marker[[i]] <- NULL
				}
		}else{
				marker <- list()
				for ( i in parameter$differ ) {
						name <- paste(i, collapse = "-vs-")
						marker[[name]] <- FindGroupMarkers(object, group.data, contrast = i,
										return.thresh = parameter$FindMarkers$pvalue,
										use.qvalue = parameter$FindMarkers$use.q,
										logfc.threshold = parameter$FindMarkers$log2fc / log2(exp(1)),
										test.use = parameter$FindMarkers$method,
										min.pct = parameter$FindMarkers$min_pct, ...
										)
				}
		}
		return(marker)
}


WriteDifferMarker <- function(marker, add_name = NULL, object = NULL) {
		for ( contrast_name in names(marker) ) {
				contrast <- strsplit(contrast_name, "-vs-")[[1]]
				if ( ! is.null(object) ) marker[[contrast_name]]$gene <- ChangeOUTName(marker[[contrast_name]]$gene, object@misc$fdata)
				data <- marker[[contrast_name]] %>%
#						mutate(name = FindFeaturesName(obj, name, "name")) %>%
						select(Cluster = cluster, GeneID = gene, # "Gene Name" = name,
								!! contrast[[1]] := mean.1, !! contrast[[2]] := mean.2,
								cells.1, cells.2, pct.1, pct.2,
								log2FC = avg_logFC, p_value = p_val, p_val_adjust = p_val_adj, significance = sig)
				WriteTable(data, file = paste(c("DifferMarker", contrast_name, add_name, "xls"), collapse = "." ))
		}
}

PlotDeGeneVolcano <- function(marker, parameter){
		plist <- list()
		flist <- list(order = names(marker))
		for ( i in names(marker) ) {
				marker[[i]] <- marker[[i]] %>% filter(cells.1 * cells.2 != 0)
				marker[[i]]$cluster <- droplevels(marker[[i]]$cluster)
				flist[[i]] <- levels(marker[[i]]$cluster)
				for ( j in levels(marker[[i]]$cluster) ) {
						name <- paste0(i, "_", j)
						y <- ifelse(parameter$FindMarkers$use.q, "-log10(p_val_adj)", "-log10(p_val)")
						ylb <- ifelse(parameter$FindMarkers$use.q, expression(-log[10](P.adjust)), expression(-log[10](Pvalue)) )
						plist[[name]] <- ggplot(subset(marker[[i]], cluster == j), aes_string(x = "avg_logFC", y = y)) +
								geom_point(aes(color = sig)) +
								scale_color_manual(values = c("up" = "#F8766D", "nosig" = "grey60", "down" = "#619CFF"), name = NULL) +
								geom_hline(yintercept = -log10(parameter$FindMarkers$pvalue),   linetype = "dashed", color = "grey40") +
								geom_vline(xintercept = c(-1,1) * parameter$FindMarkers$log2fc, linetype = "dashed", color = "grey40") +
								labs(x = expression(log[2](FC)), y = ylb) +
								xlim(-10,10) + ylim(0,50) + 
								theme_light() + 
								theme(legend.key = element_rect(color = "grey70")) 
						ggsave(plist[[name]], file = paste0("Volcano.", name, ".pdf"), width = 7, height = 7)
				}
		}
		write_yaml(flist, file = "Volcano.list", indent.mapping.sequence = TRUE)
		invisible(plist)
}



StatDeGene <- function(marker, group.by = c("contrast", "cluster"), do.plot = TRUE){
		group.by <- match.arg(group.by)
		if ( is.null(marker) || length(marker) == 0 ) return(NULL)
		stat <- do.call(rbind, lapply(marker, function(x) x %>%
								group_by(cluster, contrast, cells.1, cells.2, sig) %>%
								summarise(count = n()) %>%
								as.data.frame() ) )
#		stat2 <- reshape2::dcast(stat, cluster + contrast ~ sig, value.var = "count", fill = 0, drop = F) %>% mutate(total = up + down)
		stat2 <- reshape2::dcast(stat, ... ~ sig, value.var = "count", fill = 0, drop = T) %>%
				mutate(up = if(exists("up")) up else 0, down = if(exists("down")) down else 0) %>% 
				select(-nosig) %>% mutate(total = up + down) %>% filter(cells.1 * cells.2 != 0)
		for ( i in levels(stat2[[group.by]]) ) {
				name <- i
				if ( group.by == "contrast" ) {
#						name <- paste0("Contrast.", name)
				} else if ( group.by == "cluster" ) {
#						name <- paste0("Cluster", name)
				}
				dt <- stat2 %>% filter(.data[[group.by]] == i) %>% select(- !! group.by)
				WriteTable(dt, file = paste0("Stat.", name, ".xls"))
		}

		if ( do.plot ) {
				stat <- reshape2::melt(stat2, id.var = c("cluster", "contrast"), measure.vars = c("up", "down"), variable.name = "sig", value.name = "count" )
				plist <- list()
				for ( i in levels(stat[[group.by]]) ) {
						name <- i
						if ( group.by == "contrast" ) {
#								name <- paste0("Contrast.", name)
								group.by.inv <- "cluster"
						} else if ( group.by == "cluster" ) {
#								name <- paste0("Cluster", name)
								group.by.inv <- "contrast"
						}
						dt <- stat %>% filter(.data[[group.by]] == i) %>%
								droplevels() %>% 
								filter(count != 0)
						plist[[i]] <- ggplot(dt, aes_string(x = group.by.inv, y = "count", fill = "sig")) + 
								geom_bar(stat = "identity", position = "dodge") +
								geom_text(aes(label = count), hjust = 0.5, vjust = -0.5, position = position_dodge(0.9)) + 
								scale_fill_discrete(NULL, drop = FALSE) +
								scale_x_discrete(NULL, drop = FALSE) + 
								ylab("Number of Genes") + theme_light() +
								theme(axis.text.x = element_text(angle=45, hjust=1))
#						x.text <- unique(p$data[[as.character(p$mapping$x)[2]]])
#						w <- grid::convertWidth(unit(1, "strwidth", x.text[which.max(nchar(x.text))]), "in", T) * length(x.text)
						w <- max(7, length(levels(stat[[group.by.inv]])) * 0.8)
						ggsave(plist[[i]], file = paste0("Stat.", name, ".pdf"), width = w, height = 7)
				}
		}
}



PlotContrastHeatmap <- function(object, object.marker, group.data = NULL,
				top.num = 5, cluster.name = "seurat_clusters", do.return = FALSE,
				assay = NULL, slot = "scale.data", 
				scale.min = -2.5, scale.max = 2.5, 
				color = c("purple", "black", "yellow")
				) {
		for ( contrast_name in names(object.marker) ) {
				contrast <- strsplit(contrast_name, "-vs-")[[1]]

				annotation <- GetAnnotation(group.data, contrast, object@meta.data, cluster.name)
				cells <- rownames(annotation)
				features <- object.marker[[contrast_name]] %>% group_by(cluster) %>% top_n(top.num, abs(avg_logFC))
				features <- unique(features$gene)

				mat <- GetAssayData(object, assay = assay, slot = slot)[features, cells]
				mat <- MinMax(mat, scale.min, scale.max)
#				rownames(mat) <- FindFeaturesName(object, features, "name")
				rownames(mat) <- FindFeaturesName(object, rownames(mat), "name")
				if ( do.return ){
				} else {
						if ("showtext" %in% installed.packages()) {
								library(showtext)
								showtext_auto(enable = TRUE)
						}
						annotation_colors <- list(cluster = object@misc[["color.cluster"]])
						pheatmap::pheatmap(mat, cluster_cols = F, show_colnames = F,
								annotation_col = annotation, annotation_colors = annotation_colors,
								scale = "none", color = colorRampPalette(color)(100),
								filename = paste0("Heatmap.", contrast_name, ".pdf"))
				}

		}
}



GetAnnotation <- function(group.data, contrast, meta.data = NULL, cluster.name = "seurat_clusters", keep.group = FALSE) {
		if ( is.null(group.data) ) return(NULL)
		cell.1 <- group.data[[contrast[1]]]
		cell.2 <- group.data[[contrast[2]]]
		data <- data.frame(cell = rownames(group.data), row.names = rownames(group.data))
		if ( sum(cell.1&cell.2) == 0 ){
				data$contrast <- contrast[1]
				data$contrast[cell.2] <- contrast[2]
				data$group <- data$contrast
		} else {
				data[[contrast[1]]] <- sapply(cell.1, function(x) ifelse(x, contrast[1], NA))
				data[[contrast[2]]] <- sapply(cell.2, function(x) ifelse(x, contrast[2], NA))
				data$group <- "both"
				data$group[is.na(data[[contrast[1]]])] <- contrast[2]
				data$group[is.na(data[[contrast[2]]])] <- contrast[1]
		}
		data <- data[cell.1|cell.2, , drop = F]
		if ( ! is.null(meta.data) && ! is.null(cluster.name) ) {
				data <- cbind(data, cluster = meta.data[rownames(data), cluster.name])
				data <- data %>% group_by(group) %>% arrange(cluster, .by_group = T) %>%
						tibble::column_to_rownames(var = "cell")						
		} else {
				data <- data %>% arrange(group) %>%
						tibble::column_to_rownames(var = "cell")
		}
		if ( ! keep.group ) {
				data$group <- NULL
		}
		return(data)
}


PlotContrastFeature <- function(object, object.marker, group.data, reduction = NULL, cluster.name = "seurat_clusters",
				top.num = 5, use.name = TRUE, color = NULL,
				legend = c("right", "bottom", "none"), is.combind = TRUE, do.return = FALSE, is.demo = FALSE) {
		plots <- list()
		for ( contrast_name in names(object.marker) ) {
				contrast <- strsplit(contrast_name, "-vs-")[[1]]

				annotation <- GetAnnotation(group.data, contrast, object@meta.data, cluster.name, keep.group = T)
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
#									  scale_color_viridis_c(NULL, option = "D") +
									  theme_cowplot() +
									  theme(strip.background = element_blank(),
											strip.text = element_text(face = "bold"))
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
										outfile <- paste0("FeaturePlot.", contrast_name, ".demo.pdf")
								} else {
										outfile <- paste0("FeaturePlot.", contrast_name, ".pdf")
								}
								ggsave(ps, file = outfile, width = ncol * length(levels(data$split)) * 4, height = nrow * 4, limitsize = FALSE)
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
						theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
				if ( use.name ) {
						levels(p$data$features.plot) <- FindFeaturesName(object, levels(p$data$features.plot), "name") 
				}
				if ( ! do.return ) {
						w <- max(6, ceiling(length(features)) * 0.35 + 2)
						h <- max(6, length(levels(p$data$id)) * 0.4)
						ggsave(p, file = paste0("DotPlot.", contrast_name, ".pdf"), width = w, height = h, limitsize = FALSE)
				}
				plots[[contrast_name]] <- p
		}
		if ( do.return ) {
				return(plots)
		}
}

DotPlot <- function (object, assay = NULL, features, cols = c("lightgrey", 
    "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
    group.by = NULL, split.by = NULL, scale.by = "radius", scale.min = NA, 
    scale.max = NA) 
{
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    scale.func <- switch(EXPR = scale.by, size = scale_size, 
        radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
    data.features <- FetchData(object = object, vars = features)
    data.features$id <- if (is.null(x = group.by)) {
        Idents(object = object)
    }
    else {
        object[[group.by, drop = TRUE]]
    }
    if (!is.factor(x = data.features$id)) {
        data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)
    if (!is.null(x = split.by)) {
        splits <- object[[split.by, drop = TRUE]]
        if (length(x = unique(x = splits)) > length(x = cols)) {
            stop("Not enought colors for the number of groups")
        }
        cols <- cols[1:length(x = unique(x = splits))]
        names(x = cols) <- unique(x = splits)
        data.features$id <- paste(data.features$id, splits, sep = "_")
        unique.splits <- unique(x = splits)
        id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
            "_", rep(x = unique(x = splits), times = length(x = id.levels)))

		tmp <- unique(data.frame(id = data.features$id, splits))
		cols <- cols[tmp$splits]
		names(cols) <- tmp$id
    }
    data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
        data.use <- data.features[data.features$id == ident, 
            1:(ncol(x = data.features) - 1), drop = FALSE]
        avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
            return(mean(x = expm1(x = x)))
        })
        pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
            threshold = 0)
        return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    })
    names(x = data.plot) <- unique(x = data.features$id)
    data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
        data.use <- as.data.frame(x = data.plot[[x]])
        data.use$features.plot <- rownames(x = data.use)
        data.use$id <- x
        return(data.use)
    })
    data.plot <- do.call(what = "rbind", args = data.plot)
    if (!is.null(x = id.levels)) {
        data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }
    avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
        FUN = function(x) {
            data.use <- data.plot[data.plot$features.plot == 
                x, "avg.exp"]
            data.use <- scale(x = data.use)
            data.use <- MinMax(data = data.use, min = col.min, 
                max = col.max)
            return(data.use)
        })
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    if (!is.null(x = split.by)) {
#        avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
        color.index <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 100))
    }
    data.plot$avg.exp.scaled <- avg.exp.scaled
    data.plot$features.plot <- factor(x = data.plot$features.plot, 
        levels = rev(x = features))
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (!is.null(x = split.by)) {
#        splits.use <- vapply(X = strsplit(x = as.character(x = data.plot$id), 
#            split = "_"), FUN = "[[", FUN.VALUE = character(length = 1L), 
#            2)
        data.plot$colors <- mapply(FUN = function(color, value) {
            return(colorRampPalette(colors = c("white", color))(100)[value])
        }, color = cols[as.character(data.plot$id)], value = color.index)
    }
    color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled", 
        no = "colors")
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }
    plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", 
        y = "id")) + geom_point(mapping = aes_string(size = "pct.exp", 
        color = color.by)) + scale.func(range = c(0, dot.scale), 
        limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) + 
        labs(x = "Features", y = ifelse(test = is.null(x = split.by), 
            yes = "Identity", no = "Split Identity")) + theme_cowplot()
    if (!is.null(x = split.by)) {
        plot <- plot + scale_color_identity() +
				geom_polygon(aes(x = 0, y = 0, fill = avg.exp.scaled)) + 
				scale_fill_distiller(name = "Average Expression", palette = "Greys", direction = 1)
    }
    else if (length(x = cols) == 1) {
        plot <- plot + scale_color_distiller(palette = cols)
    }
    else {
        plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
    }
    if (is.null(x = split.by)) {
        plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
    }
    return(plot)
}
environment(DotPlot) <- asNamespace('Seurat')
assignInNamespace("DotPlot", DotPlot, ns = "Seurat")

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



Load <- function(file) {
		object <- readRDX(file)
		if ( "version" %in% slotNames(object) ) {
				if ( grepl('^2', object@version) ) {
						object <- UpdateSeuratObject(object)
				}
		}
		return(object)
}


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

