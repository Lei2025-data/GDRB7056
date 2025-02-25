VlnplotShell2 <- function(object, features = NULL, 
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
#				if (!is.null(group.point.by) && pt.size > 0) {
						data <- cbind(plots[[i]]$data, object[[group.point.by]])
						plots[[i]]$layers[[2]] <- NULL
						plots[[i]] <- plots[[i]] + geom_boxplot( data = data, outlier.color = NA, width = 0.05, color = "black")
						if ( !is.null(group.point.color) ) {
								plots[[i]] <- plots[[i]] + scale_color_manual(values = group.point.color)
						}
#				}

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
#				plots[[i]]$layers[[1]]$aes_params$alpha <- alpha
#				plots[[i]]$layers <- rev(plots[[i]]$layers)

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

PlotVlnPlot2 <- function(object, features, outpref = NULL, group.by = "seurat_clusters", cols.use = NULL, is.use.name = TRUE, ...) {
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
		VlnplotShell2(object, features, split.save.pref = outpref, group.by = group.by, titles = name, cols.use = cols.use, ...)
}

