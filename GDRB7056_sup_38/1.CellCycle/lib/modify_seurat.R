library(Seurat)

## modify for split color
DotPlot <- function (object, assay = NULL, features, cols = c("lightgrey", 
    "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
    group.by = NULL, split.by = NULL, split.fade.scale = 0.8,
    scale.by = "radius", scale.min = NA, scale.max = NA)
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
		if (!is.factor(x = splits)) {
			splits <- factor(splits)
		}
		if (nlevels(x = splits) > length(x = cols)) {
            stop("Not enought colors for the number of groups")
        }
		cols <- cols[seq(levels(x = splits))]
		names(x = cols) <- levels(x = splits)
		cols.split <- cols
        data.features$id <- paste(data.features$id, splits, sep = "_")
		id.levels <- paste0(rep(x = id.levels, each = nlevels(x = splits)), "_", rep(x = levels(x = splits), times = length(x = id.levels)))
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
	if ( nlevels(data.plot$id) == 1 ) {
		avg.exp.scaled <- scale(data.plot$avg.exp)
		avg.exp.scaled <- MinMax(data = avg.exp.scaled, min = col.min, max = col.max)
	} else {
    avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
        FUN = function(x) {
            data.use <- data.plot[data.plot$features.plot == 
                x, "avg.exp"]
            data.use <- scale(x = data.use)
            data.use <- MinMax(data = data.use, min = col.min, 
                max = col.max)
            return(data.use)
        })
	}
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
		split.col.function <- function(color, value) {
			color.hsv <- rgb2hsv(col2rgb(color))
			color.hsv[2,] <- color.hsv[2,] * (1 - split.fade.scale)
			min.color <- hsv(color.hsv[1,], color.hsv[2,], color.hsv[3,])
			return(colorRampPalette(colors = c(min.color, color))(100)[value])
		}
		data.plot$colors <- mapply(FUN = split.col.function, color = cols[as.character(data.plot$id)], value = color.index)

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
				scale_fill_gradientn(colours = sapply(split.col.function(cols.split[1]), col2grey))+
				geom_point(aes(x = 0, y = 0, alpha = avg.exp.scaled), size = 0)+
				scale_alpha("contrast", range = c(1, 1), limits = c(0, 1),
						breaks = seq(0, 1, length.out = length(cols.split)), labels = names(cols.split),
						guide = guide_legend(override.aes = list(color = cols.split, size = 3))
						)
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


## Note : 
##    alternate cut_number to Hmisc::cut2,
##    which will fail in some case
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
  data.avg <- Matrix::rowMeans(x = assay.data[pool, ,drop=F])
  data.avg <- data.avg[order(data.avg)]
#  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e30, n = nbin, labels = FALSE, right = FALSE)
  data.cut <- as.numeric(x = Hmisc::cut2(x = data.avg, m = round(x = length(x = data.avg) / (nbin + 1))))
  ctrl <- min(ctrl, min(table(data.cut)))
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
    ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use, , drop = F])
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


if ( packageVersion("Seurat") < as.numeric_version("4.0.0") ) {
FindMarkers.default <- function (object, slot = "data", counts = numeric(), cells.1 = NULL, 
    cells.2 = NULL, features = NULL, reduction = NULL, logfc.threshold = 0.25, 
    test.use = "wilcox", min.pct = 0.1, min.diff.pct = -Inf, 
    verbose = TRUE, only.pos = FALSE, max.cells.per.ident = Inf, 
    random.seed = 1, latent.vars = NULL, min.cells.feature = 3, 
    min.cells.group = 3, pseudocount.use = 1, ...) 
{
    features <- features %||% rownames(x = object)
    methods.noprefiliter <- c("DESeq2")
    if (test.use %in% methods.noprefiliter) {
        features <- rownames(x = object)
        min.diff.pct <- -Inf
        logfc.threshold <- 0
    }
    if (length(x = cells.1) == 0) {
        stop("Cell group 1 is empty - no cells with identity class ", 
            cells.1)
    }
    else if (length(x = cells.2) == 0) {
        stop("Cell group 2 is empty - no cells with identity class ", 
            cells.2)
        return(NULL)
    }
    else if (length(x = cells.1) < min.cells.group) {
        stop("Cell group 1 has fewer than ", min.cells.group, 
            " cells")
    }
    else if (length(x = cells.2) < min.cells.group) {
        stop("Cell group 2 has fewer than ", min.cells.group, 
            " cells")
    }
    else if (any(!cells.1 %in% colnames(x = object))) {
        bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.1) %in% 
            colnames(x = object))]
        stop("The following cell names provided to cells.1 are not present: ", 
            paste(bad.cells, collapse = ", "))
    }
    else if (any(!cells.2 %in% colnames(x = object))) {
        bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.2) %in% 
            colnames(x = object))]
        stop("The following cell names provided to cells.2 are not present: ", 
            paste(bad.cells, collapse = ", "))
    }
    data <- switch(EXPR = slot, scale.data = counts, object)
    if (is.null(x = reduction)) {
        thresh.min <- 0
        pct.1 <- round(x = rowSums(x = data[features, cells.1, 
            drop = FALSE] > thresh.min)/length(x = cells.1), 
            digits = 3)
        pct.2 <- round(x = rowSums(x = data[features, cells.2, 
            drop = FALSE] > thresh.min)/length(x = cells.2), 
            digits = 3)
        data.alpha <- cbind(pct.1, pct.2)
        colnames(x = data.alpha) <- c("pct.1", "pct.2")
        alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
        names(x = alpha.min) <- rownames(x = data.alpha)
        features <- names(x = which(x = alpha.min > min.pct))
        if (length(x = features) == 0) {
            stop("No features pass min.pct threshold")
        }
        alpha.diff <- alpha.min - apply(X = data.alpha, MARGIN = 1, 
            FUN = min)
        features <- names(x = which(x = alpha.min > min.pct & 
            alpha.diff > min.diff.pct))
        if (length(x = features) == 0) {
            stop("No features pass min.diff.pct threshold")
        }
    }
    else {
        data.alpha <- data.frame(pct.1 = rep(x = NA, times = length(x = features)), 
            pct.2 = rep(x = NA, times = length(x = features)))
    }
    mean.fxn <- if (is.null(x = reduction) && slot != "scale.data") {
        switch(EXPR = slot, data = function(x) {
            return(log(x = mean(x = expm1(x = x)) + pseudocount.use))
        }, function(x) {
            return(log(x = mean(x = x) + pseudocount.use))
        })
    }
    else {
        mean
    }
    data.1 <- apply(X = data[features, cells.1, drop = FALSE], 
        MARGIN = 1, FUN = mean.fxn)
    data.2 <- apply(X = data[features, cells.2, drop = FALSE], 
        MARGIN = 1, FUN = mean.fxn)
    total.diff <- (data.1 - data.2)
    if (is.null(x = reduction) && slot != "scale.data") {
        features.diff <- if (only.pos) {
            names(x = which(x = total.diff > logfc.threshold))
        }
        else {
            names(x = which(x = abs(x = total.diff) > logfc.threshold))
        }
        features <- intersect(x = features, y = features.diff)
        if (length(x = features) == 0) {
            stop("No features pass logfc.threshold threshold")
        }
    }
    if (max.cells.per.ident < Inf) {
        set.seed(seed = random.seed)
        if (length(x = cells.1) > max.cells.per.ident) {
            cells.1 <- sample(x = cells.1, size = max.cells.per.ident)
        }
        if (length(x = cells.2) > max.cells.per.ident) {
            cells.2 <- sample(x = cells.2, size = max.cells.per.ident)
        }
        if (!is.null(x = latent.vars)) {
            latent.vars <- latent.vars[c(cells.1, cells.2), , 
                drop = FALSE]
        }
    }
    if (!(test.use %in% c("negbinom", "poisson", "MAST", "LR")) && 
        !is.null(x = latent.vars)) {
        warning("'latent.vars' is only used for 'negbinom', 'poisson', 'LR', and 'MAST' tests", 
            call. = FALSE, immediate. = TRUE)
    }
    if (!test.use %in% c("wilcox", "MAST", "DESeq2")) {
        CheckDots(...)
    }
	if ( length(intersect(cells.1, cells.2)) > 0 ) {
			common.cells <- intersect(cells.1, cells.2)
			newname.cells <- paste0(common.cells, "_new")
			cells.2 <- setdiff(cells.2, cells.1)

			object <- object[, c(cells.1, cells.2, common.cells), drop = FALSE]
			colnames(object) <- c(cells.1, cells.2, newname.cells)
			cells.2 <- c(cells.2, newname.cells)
	}
    de.results <- switch(EXPR = test.use, wilcox = WilcoxDETest(data.use = object[features, 
        c(cells.1, cells.2), drop = FALSE], cells.1 = cells.1, 
        cells.2 = cells.2, verbose = verbose, ...), bimod = DiffExpTest(data.use = object[features, 
        c(cells.1, cells.2), drop = FALSE], cells.1 = cells.1, 
        cells.2 = cells.2, verbose = verbose), roc = MarkerTest(data.use = object[features, 
        c(cells.1, cells.2), drop = FALSE], cells.1 = cells.1, 
        cells.2 = cells.2, verbose = verbose), t = DiffTTest(data.use = object[features, 
        c(cells.1, cells.2), drop = FALSE], cells.1 = cells.1, 
        cells.2 = cells.2, verbose = verbose), negbinom = GLMDETest(data.use = object[features, 
        c(cells.1, cells.2), drop = FALSE], cells.1 = cells.1, 
        cells.2 = cells.2, min.cells = min.cells.feature, latent.vars = latent.vars, 
        test.use = test.use, verbose = verbose), poisson = GLMDETest(data.use = object[features, 
        c(cells.1, cells.2), drop = FALSE], cells.1 = cells.1, 
        cells.2 = cells.2, min.cells = min.cells.feature, latent.vars = latent.vars, 
        test.use = test.use, verbose = verbose), MAST = MASTDETest(data.use = object[features, 
        c(cells.1, cells.2), drop = FALSE], cells.1 = cells.1, 
        cells.2 = cells.2, latent.vars = latent.vars, verbose = verbose, 
        ...), DESeq2 = DESeq2DETest(data.use = object[features, 
        c(cells.1, cells.2), drop = FALSE], cells.1 = cells.1, 
        cells.2 = cells.2, verbose = verbose, ...), LR = LRDETest(data.use = object[features, 
        c(cells.1, cells.2), drop = FALSE], cells.1 = cells.1, 
        cells.2 = cells.2, latent.vars = latent.vars, verbose = verbose), 
        stop("Unknown test: ", test.use))
    if (is.null(x = reduction)) {
        diff.col <- ifelse(test = slot == "scale.data" || test.use == 
            "roc", yes = "avg_diff", no = "avg_logFC")
        de.results[, diff.col] <- total.diff[rownames(x = de.results)]
        de.results <- cbind(de.results, data.alpha[rownames(x = de.results), 
            , drop = FALSE])
    }
    else {
        diff.col <- "avg_diff"
        de.results[, diff.col] <- total.diff[rownames(x = de.results)]
    }
    if (only.pos) {
        de.results <- de.results[de.results[, diff.col] > 0, 
            , drop = FALSE]
    }
    if (test.use == "roc") {
        de.results <- de.results[order(-de.results$power, -de.results[, 
            diff.col]), ]
    }
    else {
        de.results <- de.results[order(de.results$p_val, -de.results[, 
            diff.col]), ]
        de.results$p_val_adj = p.adjust(p = de.results$p_val, 
            method = "bonferroni", n = nrow(x = object))
    }
    return(de.results)
}
environment(FindMarkers.default) <- asNamespace('Seurat')
assignInNamespace("FindMarkers.default", FindMarkers.default, ns = "Seurat")
#try(assignInNamespace("FindMarkers.default", FindMarkers.default, ns = "Seurat"), silent = TRUE)
## 'assignInNamespace' actually assigns internal function, but still raises an error.

} else {

PerformDE <- function (object, cells.1, cells.2, features, test.use, verbose, 
    min.cells.feature, latent.vars, densify, ...) 
{
    if (!(test.use %in% DEmethods_latent()) && !is.null(x = latent.vars)) {
        warning("'latent.vars' is only used for the following tests: ", 
            paste(DEmethods_latent(), collapse = ", "), call. = FALSE, 
            immediate. = TRUE)
    }
    if (!test.use %in% DEmethods_checkdots()) {
        CheckDots(...)
    }
	if ( length(intersect(cells.1, cells.2)) > 0 ) {
			common.cells <- intersect(cells.1, cells.2)
			newname.cells <- paste0(common.cells, "_new")
			cells.2 <- setdiff(cells.2, cells.1)

			object <- object[, c(cells.1, cells.2, common.cells), drop = FALSE]
			colnames(object) <- c(cells.1, cells.2, newname.cells)
			cells.2 <- c(cells.2, newname.cells)
	}
    data.use <- object[features, c(cells.1, cells.2), drop = FALSE]
    if (densify) {
        data.use <- as.matrix(x = data.use)
    }
    de.results <- switch(EXPR = test.use, wilcox = WilcoxDETest(data.use = data.use, 
        cells.1 = cells.1, cells.2 = cells.2, verbose = verbose, 
        ...), bimod = DiffExpTest(data.use = data.use, cells.1 = cells.1, 
        cells.2 = cells.2, verbose = verbose), roc = MarkerTest(data.use = data.use, 
        cells.1 = cells.1, cells.2 = cells.2, verbose = verbose), 
        t = DiffTTest(data.use = data.use, cells.1 = cells.1, 
            cells.2 = cells.2, verbose = verbose), negbinom = GLMDETest(data.use = data.use, 
            cells.1 = cells.1, cells.2 = cells.2, min.cells = min.cells.feature, 
            latent.vars = latent.vars, test.use = test.use, verbose = verbose), 
        poisson = GLMDETest(data.use = data.use, cells.1 = cells.1, 
            cells.2 = cells.2, min.cells = min.cells.feature, 
            latent.vars = latent.vars, test.use = test.use, verbose = verbose), 
        MAST = MASTDETest(data.use = data.use, cells.1 = cells.1, 
            cells.2 = cells.2, latent.vars = latent.vars, verbose = verbose, 
            ...), DESeq2 = DESeq2DETest(data.use = data.use, 
            cells.1 = cells.1, cells.2 = cells.2, verbose = verbose, 
            ...), LR = LRDETest(data.use = data.use, cells.1 = cells.1, 
            cells.2 = cells.2, latent.vars = latent.vars, verbose = verbose), 
        stop("Unknown test: ", test.use))
    return(de.results)
}
environment(PerformDE) <- asNamespace('Seurat')
assignInNamespace("PerformDE", PerformDE, ns = "Seurat")
}
