



CalSpatialAvgExp <- function(object, parameter = list(), groups = parameter$differ, sample.name = "orig.ident"){		
		for ( i in names(groups) ) {
				cells <- Cells(object)[object@meta.data[[sample.name]] %in% groups[[i]]]
				CalAvgExp(object[,cells], outfile = paste0("AllGene.", i, ".avg_exp.xls"))
		}
}

.Plot2Spatial <- function(object, image = NULL, cols = NULL, combine = TRUE, ...) {
  p1 <- SpatialDimPlot(object, images = image, pt.size.factor = 0, combine = FALSE)[[1]]
  p1 <- p1 + ggtitle(image) + NoLegend() + theme(plot.title = element_text(size = 18, hjust = 0.5))

  p2 <- .PlotSpatial(object, images = image, cols = cols, combine = FALSE, ...)[[1]]
  p2 <- p2 + ggtitle(NULL) + guides(fill = guide_legend(title = "cluster", override.aes = list(size = 3)))

  if ( combine ) {
		  p <- wrap_plots(p1, p2, ncol = 1)
		  return(p)
  } else {
		  return(list(p1,p2))
  }
}

PlotSpatial <- function(object, images = NULL, cols = NULL, outfile = "Spatial.cluster.pdf", scl = 6, legened.centered = TRUE, ...) {
		if ( is.null(images) ) images <- Images(object)
		plist <- lapply(images, .Plot2Spatial, object = object, cols = cols, combine = FALSE, ...)
#		p1 <- lapply(plist, `[[`, 1)
#		p2 <- lapply(plist, `[[`, 2)
#		p <- wrap_plots(c(p1, p2), nrow = 2) & theme(plot.margin = margin(0.5,0.5,0.5,0.5,"line"))
		if ( legened.centered ) {
				plist2 <- lapply(plist, function(x) wrap_plots(x, ncol = 1) + plot_layout(guides = "collect"))
		} else {
				plist2 <- lapply(plist, function(x) wrap_plots(x, ncol = 1))
		}
		w.list <- sapply(images, function(x) (diff(range(object@images[[x]]@coordinates$col)) / 127) / (diff(range(object@images[[x]]@coordinates$row) / 77)))
		p <- wrap_plots(plist2, nrow = 1, widths = w.list ) & theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "line"))
		if ( is.null(outfile) ) {
				return(p)
		} else {
				h <- scl * 2
				w <- scl * 1.2 * sum(w.list)
				ggsave(p, file = outfile, width = w, height = h, limitsize = FALSE)
		}
}

DoFindSpatialMarkers.cluster <- function(object, parameter = list(), groups = parameter$differ, sample.name = "orig.ident", ...){
  if ( future::nbrOfWorkers() > 1 ) {
    object.markers <- future.apply::future_lapply(names(groups), function(i){
      cells <- Cells(object)[object@meta.data[[sample.name]] %in% groups[[i]]]
      DoFindAllMarkers(object[, cells], parameter, ...)
    })
	names(object.markers) <- names(groups)
  } else {
    object.markers <- list()
    for ( i in names(groups) ){
      samples <- groups[[i]]
      cells <- Cells(object)[object@meta.data[[sample.name]] %in% samples]
      object.markers[[i]] <- DoFindAllMarkers(object[, cells], parameter, ...)
    } 
  }
  for ( i in names(object.markers) ){
		  object.markers[[i]] <- droplevels(object.markers[[i]])
  }
  return(object.markers)
}

DoFindSpatialMarkers.spatially <- function(object, assay = NULL, markvariogram.thres = 0.8, features = NULL, sample.name = "orig.ident", ...){
  if ( is.null(assay) ) assay <- DefaultAssay(object)
  marker <- list()
  for (image in Images(object)){
    cells <- Cells(object)[object@meta.data[[sample.name]] %in% image]
    mk <- FindVariableFeatures(object[,cells], nfeatures = 3000, assay = assay)
	select.features <- if ( is.null(features) ) VariableFeatures(mk) else features
	t1 <- Sys.time()
    mk <- FindSpatiallyVariableFeatures(mk, assay = assay, image = image,
                                        features = select.features,#VariableFeatures(object)[1:1000],
                                        selection.method = "markvariogram")
	t2 <- Sys.time()
	message("[", image, "] markvariogram time : ", format(t2 - t1) )
    vars <- grep(pattern = "r.metric", x = colnames(x = mk[[assay]]@meta.features), value = TRUE)
    tmp <- mk[[assay]]@meta.features[, vars, drop = FALSE]
	tmp <- na.omit(tmp)
    select.features <- rownames(tmp)[tmp[[vars]] <= markvariogram.thres]
    message("[", image, "] after filter(", markvariogram.thres, ") : ",
            length(select.features), " / ", nrow(tmp), " (", length(select.features)/nrow(tmp) * 100, "%)")
	t1 <- Sys.time()
    mk <- FindSpatiallyVariableFeatures(mk, assay = assay, image = image,
                                        features = select.features, #SpatiallyVariableFeatures(mk)[1:10],
                                        selection.method = "trendsceek", ...)
	t2 <- Sys.time()
	message("[", image, "] trendsceek time : ", t2 - t1 )
    vars <- c(vars, "p.val", "p_val_adjust")
    marker[[image]] <- mk[[assay]]@meta.features[, vars, drop = FALSE]
    colnames(marker[[image]]) <- c("r.metric", "p.val", "p_val_adjust")
    marker[[image]] <- marker[[image]][order(marker[[image]]$r.metric, marker[[image]]$p.val, marker[[image]]$p_val_adjust),]
  }
  return(marker)
}

GetDiffMarker <- function(marker, thres = 0.01, use.q = FALSE){
  for ( i in names(marker)) {
    if ( use.q ){
      marker[[i]] <- subset(marker[[i]], p_val_adjust < thres)
      marker[[i]] <- marker[[i]][order(marker[[i]]$p_val_adjust, marker[[i]]$r.metric),]
    }else{
      marker[[i]] <- subset(marker[[i]], p.val < thres)
      marker[[i]] <- marker[[i]][order(marker[[i]]$p.val, marker[[i]]$r.metric),]
    }
  }
  return(marker)
}

ListSpatialMarker.cluster <- function(object, object.markers, parameter = list(), is.return = F,
				assay = NULL, slot = "data", groups = parameter$differ, sample.name = "orig.ident",
				outpref = "DeGene", outsuf = "cluster.list", sep = "."){
  marker_list <- list()
  for ( differ in names(object.markers) ){
		  if ( differ %in% names(groups) ) {
    cells <- Cells(object)[object@meta.data[[sample.name]] %in% groups[[differ]]]
    marker_list[[differ]] <- ListMarker(object[, cells], object.markers[[differ]], is.return = T,
                                   assay = assay, slot = slot)
		  }
  }
  if ( is.return ){
    return(marker_list)
  } else {
    for ( differ in names(marker_list) ){
      marker_list[[differ]][["Gene ID"]] <- ChangeOUTName(marker_list[[differ]][["Gene ID"]], object@misc$fdata)
	  outfile <- paste(c(outpref, differ, outsuf), collapse = sep)
	  outfile <- paste0(outfile, ".xls")
      WriteTable(marker_list[[differ]], outfile)
    }
  }
}

ListSpatialMarker.spatial <- function(object, object.markers, is.return = F,
				assay = NULL, slot = "data",
				outpref = "DeGene", outsuf = "spatial.list", sep = "."){
  marker_list <- list()
  for ( differ in names(object.markers) ){
    Name <- FindFeaturesName(object, rownames(object.markers[[differ]]), "name", is.fast = T)
    marker_list[[differ]] <- object.markers[[differ]] %>%
      tibble::rownames_to_column(var = "gene") %>%
      mutate(name = Name[gene]) %>% 
      select("Gene ID" = gene, "Gene Name" = name, Pvalue = p.val, Qvalue = p_val_adjust)
  }
  if ( is.return ){
    return(marker_list)
  } else {
    for ( differ in names(marker_list) ){
      marker_list[[differ]][["Gene ID"]] <- ChangeOUTName(marker_list[[differ]][["Gene ID"]], object@misc$fdata)
	  outfile <- paste(c(outpref, differ, outsuf), collapse = sep)
	  outfile <- paste0(outfile, ".xls")
      WriteTable(marker_list[[differ]], outfile)
    }
  }
}

StatSpatialClusterMarker <- function(object.markers, color = NULL, outpref = "DeGene", outsuf = "cluster.stat", sep = "."){
  dt <- list()
  for ( differ in names(object.markers) ) {
    dt[[differ]] <- StatMarker(object.markers[[differ]], Cluster_name = differ, color = color, outpref = paste(c(outpref, differ, outsuf), collapse = sep))
	dt[[differ]] <- as.data.frame(dt[[differ]], stringsAsFactors = F)
  }
  dt <- lapply(dt, function(x) { x[[1]] <- names(x)[1]; names(x)[1] <- "Cluster"; x} )
  dt <- Reduce(full_join, dt)
  dt[is.na(dt)] <- '-'
  col_order <- colnames(dt)[c(1,order(as.numeric(colnames(dt)[-1]))+1)]
  dt <- dt[, col_order]
  outfile <- paste(c(outpref, outsuf), collapse = sep)
  outfile <- paste0(outfile, ".xls")
  WriteTable(dt, file = outfile)
}

PlotAboutSpatialFeatures <- function(object, features = NULL, outpref = NULL ) {
  PlotDotPlot(object, features = features, outfile = paste0(outpref, ".DotPlot.pdf"))
  #PlotFeaturePlot(object, features = features, outfile = paste0(outpref, ".Distribution.pdf"), reduction = "tsne")
  PlotHeatmapPlot(object, features = features, outfile = paste0(outpref, ".Heatmap.pdf"))
  
}



.PlotSpatial <- function(object, images = NULL, features = NULL, group.by = NULL, combine = TRUE, legend.title.position = "top", alpha = c(0.1,1), gap.between.images = TRUE, pt.size = 1.6, ... ){
	if ( is.null(images) ) images <- Images(object)

	pt.size.factor <- sapply(images, function(image) {
			cells <- intersect(Cells(object), rownames(object@images[[image]]@coordinates))
			radius.scale <- 0.01053 / object@images[[image]]@spot.radius
			diff.row <- diff(range(object@images[[image]]@coordinates[cells, "row"]))
			diff.col <- diff(range(object@images[[image]]@coordinates[cells, "col"]))
#			pt.size <- (78 * 240 - 100) / (diff.row * 240 + 140) * pt.size * radius.scale

			diff.imagerow <- diff(range(object@images[[image]]@coordinates[cells, "imagerow"]))
			diff.imagecol <- diff(range(object@images[[image]]@coordinates[cells, "imagecol"]))
			mr <- diff.imagerow / diff.row
			mc <- diff.imagecol / diff.col
#			print(c(mr,mc))
			if ( mr < 100 ) { ## is.rotated
				mr <- diff.imagecol / diff.row
				mc <- diff.imagerow / diff.col
			}
			pt.size.r <- (120         * 77)  / (mr * diff.row) * pt.size * radius.scale
			pt.size.c <- (120/sqrt(3) * 127) / (mc * diff.col) * pt.size * radius.scale
#			print(c(pt.size.r, pt.size.c))
			pt.size <- min(pt.size.r, pt.size.c)
			return(pt.size)
	})
	print(pt.size.factor)

	p.cls <- list()
	for ( i in images ) {
			p.cls <- c(p.cls, SpatialPlot(object, images = i, pt.size.factor = pt.size.factor[i],
					features = features, group.by = group.by,
					alpha = alpha, combine = F, ...))
	}

	if ( !is.null(features) ){
			p.cls <- p.cls[unlist(split(seq(p.cls), seq(features)))]
			for ( j in seq(p.cls)) {
					name <- FindFeaturesName(object, features = p.cls[[j]]$labels$fill)
					name <- stringr::str_wrap(name, width = 10)
					p.cls[[j]] <- p.cls[[j]] + guides(fill = guide_colorbar(title = name, title.position = legend.title.position))+
							theme(plot.margin = margin(t = 0.5, b = 0.5, unit = "line"), legend.title.align = 0.5)

					if ( gap.between.images && j %% length(images) == 0 ) {
						## add margin space between each features
						## feature1-image1|feature1-image2| [margin] |feature2-image1|feature2-image2
							p.cls[[j]] <- p.cls[[j]] + theme(plot.margin = margin(0.5,1,0.5,0,"line"))
					}
			}
	}

	if ( length(images) > 1 ) {
			for ( j in seq(p.cls)) {
					## add images name for each plots
					index <- j - floor((j - 0.0001)/length(images)) * length(images)
					p.cls[[j]] <- p.cls[[j]] + ggtitle(label = images[index]) +
							theme(plot.title = element_text(hjust = 0.5, margin = margin(0.5,0,0.5,0, unit = "line")))
			}
	}

	if ( combine ) {
			p <- wrap_plots(p.cls)
			return(p)
	} else {
			return(p.cls)
	}
}

WrapPlots <- function(object, plist, ncol = NULL, images = NULL, method = c("align", "byrow", "bycol") ){
		method <- match.arg(method)
		if ( is.null(ncol) ) {
				w <- ceiling(sqrt(length(plist)))
		} else {
				w <- ncol
		}
		h <- ceiling(length(plist) / w)
		if ( is.null(images) ) images <- Images(object)
		wh.ratio <- sapply(images, function(x) (diff(range(object@images[[x]]@coordinates$col)) / 127) / (diff(range(object@images[[x]]@coordinates$row) / 77)))
		is.rotated <- sapply(images, function(x) {
				row <- min(unlist(sapply(diff(sort(object@images[[x]]@coordinates$imagerow)), function(y) if (y > 10) y)))
				col <- min(unlist(sapply(diff(sort(object@images[[x]]@coordinates$imagecol)), function(y) if (y > 10) y)))
				col > row
		})
		wh.ratio[is.rotated] <- 1/wh.ratio[is.rotated]
		if ( method == "align" ) {
				p <- patchwork::wrap_plots(plist, nrow = h, ncol = w)
				max.w <- sapply(split(wh.ratio, rep(1:h, each = w)), sum)[1]
				scale.h <- h / sum(wh.ratio[1] / split(wh.ratio, rep(1:w, times = h))[[1]])
				max.w <- max.w * scale.h
		} else if ( method == "byrow" ) {
				p <- wrap_plots(lapply(split(plist, rep(1:h, each = w)), wrap_plots, nrow = 1), ncol = 1)
				max.w <- max(sapply(split(wh.ratio, rep(1:h, each = w)), sum))
		} else if ( method == "bycol" ) {
				p <- wrap_plots(lapply(split(plist, rep(1:w, times = h)), wrap_plots, ncol = 1), nrow = 1)
				scale.h <- h / max(sapply(split(wh.ratio, rep(1:w, times = h)), sum))
				max.w <- w * scale.h
		}
		return(list(p, max.w, h))
}

PlotSpatialMarkers <- function(object, features, outpref = "Marker.Spatial", images = NULL, each.plot.size = 4 ) {
		if ( is.null(images) ) {
				images <- Images(object)
		}
		if ( is.list(images) ) {
				if ( length(names(images)) != length(images) ){
						stop("if 'images' is.list, length(names(images)) should be equal to length(images).")
				}
				groups <- lapply(images, intersect, Images(object))
		} else {
				groups <- list(intersect(images, Images(object)))
		}
		return.p <- list()
		for ( i in seq(groups) ) {
				group <- names(groups)[[i]]
				if ( is.null(group) ) {
						images <- groups[[i]]
						surf <- "pdf"
				} else {
						images <- groups[[i]]
						surf <- paste0(group, ".pdf")
				}
				p.cls <- .PlotSpatial(object, images = images, features = features, combine = F)
				w.single <- ceiling(sqrt(length(p.cls)))
				n.single <- floor(w.single / length(images))
				w <- n.single * length(images)
				h <- ceiling(length(p.cls) / w)
				wh.ratio <- sapply(images, function(x) (diff(range(object@images[[x]]@coordinates$col)) / 127) / (diff(range(object@images[[x]]@coordinates$row) / 77)))
				max.w <- n.single * sum(wh.ratio)
				p <- patchwork::wrap_plots(p.cls, nrow = h, ncol = w)
				if ( ! is.null(outpref) ) {
						outfile <- paste0(outpref, ".", surf)
						ggsave(p, file = outfile, width = max.w * each.plot.size, height = h * each.plot.size, limitsize = FALSE)
				}
				return.p[[i]] <- p
		}
		if ( length(return.p) == 1 ) return.p <- return.p[[1]]
		return(invisible(return.p))
}

PlotSpatialMarkers.cluster <- function(object, object.markers, parameter = list(), outpref = NULL, groups = parameter$differ, top_num = parameter$top$cluster, sample.name = "orig.ident", ...){
  for (differ in names(object.markers) ){
    images <- groups[[differ]]
    cells <- Cells(object)[object@meta.data[[sample.name]] %in% images]
    object.tmp <- object[, cells]
	object.tmp@meta.data <- droplevels(object.tmp@meta.data)
    
    top <- FindTopMarker(object.markers[[differ]], top_num = top_num,
                         outfile = NULL)
    features <- unique(top$gene)
    if ( !is.null(outpref) ) {
      dt <- CalAvgExp(object.tmp, features = top$gene, is.return = TRUE)
      colnames(dt) <- paste0("Cluster ", colnames(dt))
      dt <- top %>% select(Cluster = cluster, Gene_ID = gene) %>%
            mutate(Gene_name = FindFeaturesName(object.tmp, Gene_ID, "name")) %>%
			mutate(Gene_ID = ChangeOUTName(Gene_ID, object.tmp@misc$fdata)) %>% 
			ungroup() %>% cbind(dt)
	  WriteTable(dt, paste0(outpref, ".", differ, ".cluster.avg_exp.xls"))

      PlotDotPlot(object.tmp, features = features, outfile = paste0(outpref, ".", differ, ".cluster.DotPlot.pdf"))
      PlotHeatmapPlot(object.tmp, features = features, outfile = paste0(outpref, ".", differ, ".cluster.Heatmap.pdf"))
    }
    
    plist <- list()
    for ( i in levels(top$cluster) ) {
      features <- top$gene[top$cluster == i]
      p.cls <- .PlotSpatial(object.tmp, images, features, combine = F, ...)
      p.cls[[1]] <- p.cls[[1]] + ylab(i) +
        theme(axis.title = element_text(margin = margin(0,1,0,0,"line"), size = 20),
              axis.title.x = element_blank())
      plist[[i]] <- wrap_plots(p.cls, nrow = 1)
    }
    p <- wrap_plots(plist, ncol = 1)
	w <- top_num * sum(sapply(images, function(x) (diff(range(object@images[[x]]@coordinates$col)) / 127) / (diff(range(object@images[[x]]@coordinates$row) / 77))))	
    h <- length(plist) * 1.2
    scl <- 4
    if ( !is.null(outpref) ) {
      ggsave(plist[[1]], file = paste0(outpref, ".", differ, ".cluster.SpatialPlot.demo.pdf"), width = scl * w, height = scl * 1.2, limitsize = F)
      ggsave(p, file = paste0(outpref, ".", differ, ".cluster.SpatialPlot.pdf"), width = scl * w, height = scl * h, limitsize = F)
    }
  }
}

PlotSpatialMarkers.spatially <- function(object, marker, outpref = NULL, sample.name = "orig.ident"){
  for ( image in names(marker) ){
    #    feature <- GetSigGene(marker)
    features <- rownames(marker[[image]])
	if ( length(features) > 0 ) {
    plist <- .PlotSpatial(object, image, features, combine = F)	
	p <- wrap_plots(plist)
    w <- ceiling(sqrt(length(features)))
    h <- ceiling(length(features) / w)
	w <- w * (diff(range(object@images[[image]]@coordinates$col)) / 127) / (diff(range(object@images[[image]]@coordinates$row) / 77))
    scl <- 4
    if ( !is.null(outpref) ) {
      ggsave(plist[[1]], file = paste0(outpref, ".", image, ".spatial.SpatialPlot.demo.pdf"), width = 7, height = 7)
      ggsave(p, file = paste0(outpref, ".", image, ".spatial.SpatialPlot.pdf"), width = scl * w, height = scl * h, limitsize = F)
      cells <- Cells(object)[object@meta.data[[sample.name]] %in% image]
      object.tmp <- object[, cells]
	  object.tmp@meta.data <- droplevels(object.tmp@meta.data)
      CalAvgExp(object.tmp, features = features, outfile = paste0(outpref, ".", image, ".spatial.avg_exp.xls"))
      PlotDotPlot(object.tmp, features = features, outfile = paste0(outpref, ".", image, ".spatial.DotPlot.pdf"))
      PlotHeatmapPlot(object.tmp, features = features, outfile = paste0(outpref, ".", image, ".spatial.Heatmap.pdf"))
    }
	}
  }
}

if ( FALSE && 'trendsceek' %in% installed.packages() ) {
trendsceek_test <- function (pp, nrand = 10000, ncores = NULL,
                             alpha_env = 0.1/ifelse(is.numeric(pp[["marks"]]), length(pp[["marks"]]), ncol(pp[["marks"]])),
                             alpha_nom_early = (alpha_bh * 4)/ifelse(ifelse(is.numeric(pp[["marks"]]), length(pp[["marks"]]), ncol(pp[["marks"]])) >= 500, 10, 1),
                             alpha_bh = 0.05, method = c("Emark", "markcorr", "markvario", "Vmark")) {
  bp_param = BiocParallel::MulticoreParam()
  if(!is.null(ncores)) bp_param = BiocParallel::MulticoreParam(workers = ncores)
  message("[trendsceek_test] use workers : ", bp_param$workers)
  marx = pp[["marks"]]
  if (is.numeric(marx)) {
    nfeats = 1
    feats = 1
  }
  else {
    nfeats = ncol(marx)
    feats = colnames(marx)
  }
  tstat_list = BiocParallel::bplapply(1:nfeats, calc_trendstats, 
                                      BPPARAM = bp_param, pp = pp, n.rand = nrand, alpha_env = alpha_env, 
                                      alpha_nom_early = alpha_nom_early, method = method)
  names(tstat_list) = feats
  supstats_list = tstat2supstat(tstat_list)
  trendstat_list = list(tstat = tstat_list, supstats = supstats_list) 
  
  return(trendstat_list)
}
environment(trendsceek_test) <- asNamespace('trendsceek')
assignInNamespace("trendsceek_test", trendsceek_test, ns = "trendsceek")


RunTrendsceek <- function(spatial.location, data, nrand = 100, ...) {
  pp <- trendsceek::pos2pp(pos_mat = spatial.location)
  pp <- trendsceek::set_marks(pp = pp, gene.marks = data)
  trendstat_list <- trendsceek::trendsceek_test(pp, nrand = nrand, method = "markvario", ...)
  #trendstat <- trendsceek::extract_sig_genes(trendstat_list, alpha = 1)
  svf.info <- trendstat_list$supstats$markvario[, c("min.pval", "p.bh"), drop = F]
  colnames(svf.info) <- c("p.val", "p_val_adjust")
  return(svf.info)
}

FindSpatiallyVariableFeatures.default <- function (object, spatial.location,
                                                   selection.method = c("markvariogram", "moransi", "trendsceek"),
                                                   r.metric = 5, x.cuts = NULL, y.cuts = NULL, verbose = TRUE, ...) {
  if (ncol(x = object) != nrow(x = spatial.location)) {
    stop("Please provide the same number of observations as spatial locations.")
  }
  if (!is.null(x = x.cuts) & !is.null(x = y.cuts)) {
    binned.data <- BinData(data = object, pos = spatial.location, 
                           x.cuts = x.cuts, y.cuts = y.cuts, verbose = verbose)
    object <- binned.data$data
    spatial.location <- binned.data$pos
  }
  svf.info <- switch(EXPR = selection.method,
                     markvariogram = RunMarkVario(spatial.location = spatial.location, data = object),
                     moransi = RunMoransI(data = object, pos = spatial.location, verbose = verbose),
                     trendsceek = RunTrendsceek(spatial.location = spatial.location, data = object, ...),
                     stop("Invalid selection method. Please choose one of: markvariogram, moransi."))
  return(svf.info)
}
environment(FindSpatiallyVariableFeatures.default) <- asNamespace('Seurat')
assignInNamespace("FindSpatiallyVariableFeatures.default", FindSpatiallyVariableFeatures.default, ns = "Seurat")

get_envstats <- function (tstat.df, nullstats.df, alpha_env) 
{
    null.mean = tstat.df[, "exp"]
    null.devs = apply(nullstats.df, 2, function(j.sim, null.mean) {
        abs(j.sim - null.mean)
    }, null.mean = null.mean)
	null.devs[is.nan(null.devs)] <- 0
    null.sorted = t(apply(null.devs, 1, sort))
    null.global = sort(apply(null.sorted, 2, max), decreasing = TRUE)
    nsim = ncol(nullstats.df)
    obs.dev = abs(tstat.df[, "obs"] - null.mean)
    hi.rank = unlist(lapply(obs.dev, function(j.r, null.global) {
        length(which(null.global >= j.r))
    }, null.global = null.global))
    hi.rank = hi.rank + 1
    p = hi.rank/(nsim + 1)
    kth.nullval = floor((nsim + 1) * alpha_env)
    if (kth.nullval == 0) {
        kth.nullval = 1
    }
    max.dev = null.global[kth.nullval]
    hi.global = null.mean + max.dev
    lo.global = null.mean - max.dev
    env.rel.dev = obs.dev/max.dev
    mean.rel.dev = obs.dev/null.mean
    envstats = cbind(obs.dev, lo.global, hi.global, env.rel.dev, 
        mean.rel.dev, p)
    return(envstats)
}
environment(get_envstats) <- asNamespace('trendsceek')
assignInNamespace("get_envstats", get_envstats, ns = "trendsceek")
}


SpatialPlot <- function (object, group.by = NULL, features = NULL, images = NULL, 
    cols = NULL, image.alpha = 1, crop = TRUE, slot = "data", assay = DefaultAssay(object), 
    min.cutoff = NA, max.cutoff = NA, cells.highlight = NULL, 
    cols.highlight = c("#DE2D26", "grey50"), facet.highlight = FALSE, 
    label = FALSE, label.size = 5, label.color = "white", label.box = TRUE, 
    repel = FALSE, ncol = NULL, combine = TRUE, pt.size.factor = 1.6, 
    alpha = c(1, 1), stroke = 0.25, interactive = FALSE, do.identify = FALSE, 
    identify.ident = NULL, do.hover = FALSE, information = NULL) 
{
    if (isTRUE(x = do.hover) || isTRUE(x = do.identify)) {
        warning("'do.hover' and 'do.identify' are deprecated as we are removing plotly-based interactive graphics, use 'interactive' instead for Shiny-based interactivity", 
            call. = FALSE, immediate. = TRUE)
        interactive <- TRUE
    }
    if (!is.null(x = group.by) & !is.null(x = features)) {
        stop("Please specific either group.by or features, not both.")
    }
    images <- images %||% Images(object = object, assay = assay)
    if (is.null(x = features)) {
        if (interactive) {
            return(ISpatialDimPlot(object = object, image = image, 
                group.by = group.by, alpha = alpha))
        }
        group.by <- group.by %||% "ident"
        object[["ident"]] <- Idents(object = object)
        data <- object[[group.by]]
        for (group in group.by) {
            if (!is.factor(x = data[, group])) {
                data[, group] <- factor(x = data[, group])
            }
        }
    }
    else {
        if (interactive) {
            return(ISpatialFeaturePlot(object = object, feature = features[1], 
                image = images[1], slot = slot, alpha = alpha))
        }
		DefaultAssay(object) <- assay ### add 
        data <- FetchData(object = object, vars = features, slot = slot)
        features <- colnames(x = data)
        min.cutoff <- mapply(FUN = function(cutoff, feature) {
            return(ifelse(test = is.na(x = cutoff), yes = min(data[, 
                feature]), no = cutoff))
        }, cutoff = min.cutoff, feature = features)
        max.cutoff <- mapply(FUN = function(cutoff, feature) {
            return(ifelse(test = is.na(x = cutoff), yes = max(data[, 
                feature]), no = cutoff))
        }, cutoff = max.cutoff, feature = features)
        check.lengths <- unique(x = vapply(X = list(features, 
            min.cutoff, max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
        if (length(x = check.lengths) != 1) {
            stop("There must be the same number of minimum and maximum cuttoffs as there are features")
        }
        data <- sapply(X = 1:ncol(x = data), FUN = function(index) {
            data.feature <- as.vector(x = data[, index])
            min.use <- SetQuantile(cutoff = min.cutoff[index], 
                data.feature)
            max.use <- SetQuantile(cutoff = max.cutoff[index], 
                data.feature)
            data.feature[data.feature < min.use] <- min.use
            data.feature[data.feature > max.use] <- max.use
            return(data.feature)
        })
        colnames(x = data) <- features
        rownames(x = data) <- Cells(x = object)
    }
    if (length(x = images) == 0) {
        images <- Images(object = object)
    }
    if (length(x = images) < 1) {
        stop("Could not find any spatial image information")
    }
    features <- colnames(x = data)
    colnames(x = data) <- features
    rownames(x = data) <- colnames(x = object)
    facet.highlight <- facet.highlight && (!is.null(x = cells.highlight) && 
        is.list(x = cells.highlight))
    if (do.hover) {
        if (length(x = images) > 1) {
            images <- images[1]
            warning("'do.hover' requires only one image, using image ", 
                images, call. = FALSE, immediate. = TRUE)
        }
        if (length(x = features) > 1) {
            features <- features[1]
            type <- ifelse(test = is.null(x = group.by), yes = "feature", 
                no = "grouping")
            warning("'do.hover' requires only one ", type, ", using ", 
                features, call. = FALSE, immediate. = TRUE)
        }
        if (facet.highlight) {
            warning("'do.hover' requires no faceting highlighted cells", 
                call. = FALSE, immediate. = TRUE)
            facet.highlight <- FALSE
        }
    }
    if (facet.highlight) {
        if (length(x = images) > 1) {
            images <- images[1]
            warning("Faceting the highlight only works with a single image, using image ", 
                images, call. = FALSE, immediate. = TRUE)
        }
        ncols <- length(x = cells.highlight)
    }
    else {
        ncols <- length(x = images)
    }
    plots <- vector(mode = "list", length = length(x = features) * 
        ncols)
    for (i in 1:ncols) {
        plot.idx <- i
        image.idx <- ifelse(test = facet.highlight, yes = 1, 
            no = i)
        image.use <- object[[images[[image.idx]]]]
        coordinates <- GetTissueCoordinates(object = image.use)
		coordinates <- coordinates[intersect(rownames(coordinates), rownames(data)),] ## added this
        highlight.use <- if (facet.highlight) {
            cells.highlight[i]
        }
        else {
            cells.highlight
        }
        for (j in 1:length(x = features)) {
            cols.unset <- is.factor(x = data[, features[j]]) && 
                is.null(x = cols)
            if (cols.unset) {
                cols <- scales::hue_pal()(n = length(x = levels(x = data[, 
                  features[j]])))
                names(x = cols) <- levels(x = data[, features[j]])
            }
            plot <- SingleSpatialPlot(data = cbind(coordinates, 
                data[rownames(x = coordinates), features[j], 
                  drop = FALSE]), image = image.use, image.alpha = image.alpha, 
                col.by = features[j], cols = cols, alpha.by = if (is.null(x = group.by)) {
                  features[j]
                }
                else {
                  NULL
                }, geom = if (inherits(x = image.use, what = "STARmap")) {
                  "poly"
                }
                else {
                  "spatial"
                }, cells.highlight = highlight.use, cols.highlight = cols.highlight, 
                pt.size.factor = pt.size.factor, stroke = stroke, 
                crop = crop)
            if (is.null(x = group.by)) {
                plot <- plot + scale_fill_gradientn(name = features[j], 
                  colours = SpatialColors(n = 100)) + theme(legend.position = "top") + 
                  scale_alpha(range = alpha) + guides(alpha = FALSE)
            }
            else if (label) {
                plot <- LabelClusters(plot = plot, id = ifelse(test = is.null(x = cells.highlight), 
                  yes = features[j], no = "highlight"), geom = if (inherits(x = image.use, 
                  what = "STARmap")) {
                  "GeomPolygon"
                }
                else {
                  "GeomSpatial"
                }, repel = repel, size = label.size, color = label.color, 
                  box = label.box, position = "nearest")
            }
            if (j == 1 && length(x = images) > 1 && !facet.highlight) {
                plot <- plot + ggtitle(label = images[[image.idx]]) + 
                  theme(plot.title = element_text(hjust = 0.5))
            }
            if (facet.highlight) {
                plot <- plot + ggtitle(label = names(x = cells.highlight)[i]) + 
                  theme(plot.title = element_text(hjust = 0.5)) + 
                  NoLegend()
            }
            plots[[plot.idx]] <- plot
            plot.idx <- plot.idx + ncols
            if (cols.unset) {
                cols <- NULL
            }
        }
    }
    if (length(x = images) > 1 && combine) {
        plots <- wrap_plots(plots = plots, ncol = length(x = images))
    }
    else if (length(x = images == 1) && combine) {
        plots <- wrap_plots(plots = plots, ncol = ncol)
    }
    return(plots)
}
environment(SpatialPlot) <- asNamespace("Seurat")
assignInNamespace("SpatialPlot", SpatialPlot, ns = "Seurat")


SpatialDimPlot <- function (object, group.by = NULL, images = NULL, cols = NULL, 
          crop = TRUE, cells.highlight = NULL, cols.highlight = c("#DE2D26", "grey50"),
          facet.highlight = FALSE, label = FALSE, label.size = 7, 
          label.color = "white", repel = FALSE, ncol = NULL, combine = TRUE, 
          pt.size.factor = 1.6, alpha = c(1, 1), stroke = 0.25, label.box = TRUE, 
          interactive = FALSE, information = NULL) 
{
  ## add "cols"
  return(SpatialPlot(object = object, group.by = group.by, cols = cols,
                     images = images, crop = crop, cells.highlight = cells.highlight, 
                     cols.highlight = cols.highlight, facet.highlight = facet.highlight, 
                     label = label, label.size = label.size, label.color = label.color, 
                     repel = repel, ncol = ncol, combine = combine, pt.size.factor = pt.size.factor, 
                     alpha = alpha, stroke = stroke, label.box = label.box, 
                     interactive = interactive, information = information))
}
environment(SpatialDimPlot) <- asNamespace("Seurat")
assignInNamespace("SpatialDimPlot", SpatialDimPlot, ns = "Seurat")


