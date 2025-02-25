warning('    You may also need to load below packages : \n        "Seurat", "dplyr", "ggplot2", "patchwork"', call. = FALSE)
.GetMetaData <-
function (object, cols = NULL) 
{
    name <- names(cols)
    name[is.na(name) | name == ""] <- cols[is.na(name) | name == 
        ""]
    if (is.null(name)) 
        name <- cols
    names(cols) <- name
    cols <- cols[cols %in% colnames(obj@meta.data)]
    metadata <- object@meta.data[, cols]
    colnames(metadata) <- names(cols)
    metadata <- cbind(Cells = rownames(metadata), metadata)
    return(metadata)
}
.PlotCluster <-
function (object, reduction = NULL, cells = NULL, outfile = NULL, 
    p1.group.by = "orig.ident", p1.color = NULL, p1.label = FALSE, 
    p2.group.by = "seurat_clusters", p2.color = NULL, p2.label = TRUE, 
    width = 14, height = 6, ...) 
{
    if (is.null(p1.color)) {
        p1.color <- switch(p1.group.by, Groups = object@misc[["color.group"]], 
            orig.ident = object@misc[["color.sample"]], seurat_clusters = object@misc[["color.cluster"]])
    }
    if (is.null(p2.color)) {
        p2.color <- switch(p2.group.by, Groups = object@misc[["color.group"]], 
            orig.ident = object@misc[["color.sample"]], seurat_clusters = object@misc[["color.cluster"]])
    }
    p1 <- DimPlot(object, reduction = reduction, cells = cells, 
        group.by = p1.group.by, cols = p1.color, label = p1.label, 
        ...)
    p2 <- DimPlot(object, reduction = reduction, cells = cells, 
        group.by = p2.group.by, cols = p2.color, label = p2.label, 
        ...)
    p1 <- p1 + dot_theme_default()
    p2 <- p2 + dot_theme_default()
    p <- p1 + p2
    if (is.null(outfile)) {
        return(p)
    }
    else {
        ggsave(p, file = outfile, width = width, height = height, 
            limitsize = FALSE)
    }
}
.PlotClusterStat <-
function (object, stat.what = "seurat_clusters", group.by = "orig.ident", 
    color.st = NULL, color.gb = NULL, outpref = NULL) 
{
    if (is.null(color.st)) {
        color.st <- switch(stat.what, seurat_clusters = object@misc$color.cluster, 
            orig.ident = object@misc$color.sample, Groups = object@misc$color.group)
    }
    if (is.null(color.gb)) {
        color.gb <- switch(group.by, seurat_clusters = object@misc$color.cluster, 
            orig.ident = object@misc$color.sample, Groups = object@misc$color.group)
    }
    name.st <- switch(stat.what, seurat_clusters = "Cluster", 
        orig.ident = "Samples", stat.what)
    name.gb <- switch(group.by, seurat_clusters = "Cluster", 
        orig.ident = "Samples", group.by)
    stat.what <- as.name(stat.what)
    group.by <- as.name(group.by)
    stat_sample <- object@meta.data %>% group_by(`:=`(!!name.gb, 
        !!group.by), `:=`(!!name.st, !!stat.what)) %>% summarise(`Number of cells` = n())
    p <- list()
    p[["by"]] <- ggplot(stat_sample, aes_(x = as.name(name.gb), 
        y = ~`Number of cells`, fill = as.name(name.st))) + scale_fill_manual(values = color.st)
    p[["in"]] <- ggplot(stat_sample, aes_(x = as.name(name.st), 
        y = ~`Number of cells`, fill = as.name(name.gb))) + scale_fill_manual(values = color.gb)
    geom_stack <- geom_bar(stat = "identity", position = "stack")
    geom_fill <- geom_bar(stat = "identity", position = "fill")
    if (is.null(outpref)) {
        outpref <- paste0(name.st, ".stat")
    }
    for (i in names(p)) {
        p[[i]] <- p[[i]] + bar_theme_default()
        ggsave(p[[i]] + geom_stack, file = paste0(outpref, ".", 
            i, name.gb, ".pdf"), height = 6, width = 8)
        ggsave(p[[i]] + geom_fill, file = paste0(outpref, ".", 
            i, name.gb, ".pct.pdf"), height = 6, width = 8)
    }
}
.PlotDensityPlot <-
function (object, features = NULL, reduction = "umap", is.return = FALSE, 
    outpref = NULL, is.consider.exp = FALSE, is.filter.noexp = TRUE) 
{
    dt <- as.data.frame(object[[reduction]]@cell.embeddings)
    if (is.null(features)) {
        dt$density <- KDE(x = dt[[1]], y = dt[[2]])
        outname <- paste0(outpref, "DensityPlot.pdf")
    }
    else {
        features <- features[1]
        exp <- t(as.data.frame(GetAssayData(object)[features, 
            , drop = F]))
        dt <- cbind(dt, exp)
        dt.filter <- dt[dt[[3]] > 0, ]
        if (is.consider.exp) {
            dt.filter$density <- KDE(x = dt.filter[[1]], y = dt.filter[[2]], 
                z = dt.filter[[3]])
        }
        else {
            dt.filter$density <- KDE(x = dt.filter[[1]], y = dt.filter[[2]])
        }
        dt <- full_join(dt, dt.filter)
        dt <- dt[order(-dt$density, na.last = F), ]
        name <- FindFeaturesName(object, features)
        name <- gsub("[ /\\]", "_", name)
        outname <- paste0(outpref, "DensityPlot.", name, ".pdf")
    }
    p <- ggplot(dt, aes_string(x = colnames(dt)[1], y = colnames(dt)[2], 
        color = "density")) + geom_point() + scale_color_viridis_c(option = "A", 
        na.value = "grey90")
    if (!is.null(features)) 
        p <- p + ggtitle(name) + theme(plot.title = element_text(hjust = 0.5))
    p <- p + dot_theme_default()
    if (is.return) {
        return(list(p, dt))
    }
    else {
        ggsave(p, file = outname, width = 8, height = 7)
    }
}
.PlotFeaturePlot <-
function (object, features, outfile = NULL, reduction = NULL, 
    is.use.name = TRUE, color.high = "blue", color.low = "lightgrey", 
    show.cluster.label = FALSE, nCol = NULL, plot.basic.size = 4, 
    group.by = "seurat_clusters", is.combine = TRUE) 
{
    if (show.cluster.label) 
        Idents(object) <- group.by
    plots <- FeaturePlot(object, features = features, order = TRUE, 
        reduction = reduction, combine = FALSE, label = show.cluster.label, 
        cols = c(color.low, color.high))
    if (is.use.name) {
        name <- FindFeaturesName(object, features)
        for (i in seq(plots)) {
            plots[[i]] <- plots[[i]] + ggtitle(name[i])
        }
        names(plots) <- name
    }
    else {
        names(plots) <- features
    }
    if (is.null(nCol)) 
        nCol <- ceiling(sqrt(length(features)))
    nRow <- ceiling(length(features)/nCol)
    p <- wrap_plots(plots, ncol = nCol) & dot_theme_default()
    if (is.null(outfile)) {
        return(p)
    }
    else {
        ggsave(p, file = outfile, width = plot.basic.size * (6/5) * 
            nCol, height = plot.basic.size * nRow, limitsize = FALSE)
    }
}
.StatCluster <-
function (object, outpref = "Cluster.stat") 
{
    Cluster.stat <- object@meta.data %>% group_by(Cluster = seurat_clusters) %>% 
        summarise(`Cells number` = n(), `Median Genes per Cell` = median(nFeature_RNA), 
            `Median UMI Counts per Cell` = median(nCount_RNA))
    WriteTable(Cluster.stat, paste0(outpref, ".xls"))
}
.StatCluster_by <-
function (object, group.by = "orig.ident", outpref = "Cluster.stat") 
{
    name <- switch(group.by, orig.ident = "Samples", group.by)
    group.by <- as.name(group.by)
    Cluster.stat <- object@meta.data %>% group_by(name = !!group.by, 
        Cluster = seurat_clusters) %>% summarise(y = n()) %>% 
        mutate(Cluster = factor(Cluster, levels = c("Total", 
            levels(Cluster)))) %>% full_join(x = object@meta.data %>% 
        group_by(name = !!group.by) %>% summarise(Cluster = factor("Total", 
        levels = c("Total", levels(seurat_clusters))), y = sum(n()))) %>% 
        reshape2::dcast(Cluster ~ name, fill = 0) %>% mutate_if(is.numeric, 
        list(~paste0(., " (", round(./.[1] * 100, 2), "%)")))
    WriteTable(Cluster.stat, paste0(outpref, ".", name, ".xls"))
}
AddFData <-
function (object, ref_name_file = NULL, col.name = NULL) 
{
    if (!is.null(ref_name_file) && file.exists(ref_name_file)) {
        fdata <- read.table(ref_name_file, row.names = 1, stringsAsFactors = F, 
            sep = "\t")
        colnames(fdata) <- c("merge_name", "name", "type")
    }
    else {
        fdata <- data.frame(name = rownames(object), row.names = rownames(object), 
            stringsAsFactors = F)
    }
    fdata$merge_name <- fdata$name
    fdata$merge_name[fdata$merge_name == "-"] <- rownames(fdata)[fdata$merge_name == 
        "-"]
    index <- c(which(duplicated(fdata$merge_name, fromLast = T)), 
        which(duplicated(fdata$merge_name, fromLast = F)))
    fdata$merge_name[index] <- paste0(fdata$merge_name[index], 
        " (", rownames(fdata)[index], ")")
    fdata <- AddUnderscore(fdata)
    if (is.null(col.name)) {
        return(fdata)
    }
    else {
        object@misc[[col.name]] <- fdata
        return(object)
    }
}
AddUnderscore <-
function (data) 
{
    if (!exists("underscore", data) || !exists("dash", data)) {
        data$underscore <- rownames(data)
        data$dash <- gsub("_", "-", rownames(data))
    }
    return(data)
}
autothres <-
function (data, name = c("default", "nCount_RNA", "nFeature_RNA", 
    "percent.mito"), bin = 100, digits = 2) 
{
    thres <- Otsu(data = data, bin = bin)
    name <- match.arg(name)
    if (name == "nCount_RNA") {
        up <- signif(thres[2], digits)
        down <- -Inf
    }
    else if (name == "nFeature_RNA") {
        up <- signif(thres[2], digits)
        down <- signif(thres[1], digits)
        if (up > 200) {
            down <- max(down, 200)
        }
    }
    else if (name == "percent.mito") {
        if (all(data <= 1)) 
            up <- if (thres[2] <= 0.1) 
                0.1
            else if (thres[2] <= 0.25) 
                0.25
            else signif(thres[2], digits)
        else up <- if (thres[2] <= 10) 
            10
        else if (thres[2] <= 25) 
            25
        else signif(thres[2], digits)
        down <- -Inf
    }
    else {
        up <- signif(thres[2], digits)
        down <- signif(thres[1], digits)
    }
    return(c(down, up))
}
bar_theme_default <-
function () 
{
    library(ggplot2)
    options(scipen = -1)
    mytheme <- theme_bw() + theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.border = element_rect(color = "#000000", 
            size = 1), axis.text = element_text(color = "#000000", 
            size = 12), axis.text.x = element_text(angle = 45, 
            hjust = 1, vjust = 1), axis.text.y = element_text(hjust = 1, 
            vjust = 0.5), axis.title = element_text(color = "#000000", 
            size = 16, face = "bold"), axis.title.x = element_text(margin = margin(2.5, 
            0, 2.5, 0, "mm")), axis.title.y = element_text(margin = margin(0, 
            2.5, 0, 2.5, "mm")), axis.ticks = element_line(color = "#000000", 
            size = 0.6), axis.ticks.length = unit(0.13, "cm"), 
        legend.title = element_blank(), plot.title = element_text(size = 20, 
            face = "bold", hjust = 0.5), )
    if ("showtext" %in% installed.packages()) {
        library(showtext)
        showtext_auto(enable = TRUE)
        font_add("Arial", regular = "/Bio/Bin/pipeline/GeneralPlot/v1.0/fonts/tff/msttcore/arial.ttf")
        mytheme <- mytheme + theme(text = element_text(family = "Arial"))
    }
    mytheme
}
bin_stat <-
function (data, bin = 256) 
{
    data <- as.numeric(data)
    b <- floor((data - min(data))/diff(range(data)) * bin)
    b[b == bin] <- bin - 1
    names(b) <- data
    return(b)
}
box_theme_default <-
function () 
{
    library(ggplot2)
    mytheme <- theme_bw() + theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.border = element_rect(color = "#000000", 
            size = 1), axis.text = element_text(color = "#000000", 
            size = 12), axis.text.x = element_text(angle = 0, 
            hjust = 0.5, vjust = 0.5), axis.text.y = element_text(hjust = 0.5, 
            vjust = 0.5), axis.title = element_text(color = "#000000", 
            size = 16, face = "bold"), axis.title.x = element_text(margin = margin(2.5, 
            0, 2.5, 0, "mm")), axis.title.y = element_text(margin = margin(0, 
            2.5, 0, 2.5, "mm")), axis.ticks = element_line(color = "#000000", 
            size = 0.6), axis.ticks.length = unit(0.13, "cm"), 
        legend.title = element_blank(), legend.justification = "center", 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        )
    if ("showtext" %in% installed.packages()) {
        library(showtext)
        showtext_auto(enable = TRUE)
        font_add("Arial", regular = "/Bio/Bin/pipeline/GeneralPlot/v1.0/fonts/tff/msttcore/arial.ttf")
        mytheme <- mytheme + theme(text = element_text(family = "Arial"))
    }
    mytheme
}
CalAvgExp <-
function (object, features = NULL, group.by = NULL, assay = NULL, 
    slot = "data", is.expm1 = ifelse(slot == "data", TRUE, FALSE), 
    is.return = FALSE, is.reverse = FALSE, is.bulk = FALSE, outfile = "AllGene.avg_exp.xls") 
{
    data <- GetAssayData(object, assay = assay, slot = slot)
    if (!is.null(features)) 
        data <- data[features, ]
    if (is.expm1) 
        data <- expm1(data)
    if (!is.null(group.by)) 
        Idents(object) <- group.by
    mean_exp <- do.call(cbind, by(colnames(object), Idents(object), 
        function(x) Matrix::rowMeans(data[, if (is.reverse) {
            levels(x)[-as.numeric(x)]
        }
        else {
            x
        }, drop = F])))
    if (is.bulk) {
        mean_exp <- cbind(bulk = Matrix::rowMeans(data), mean_exp)
    }
    if (is.return) {
        return(mean_exp)
    }
    else {
        colnames(mean_exp) <- paste("Cluster", colnames(mean_exp))
        Gene_name <- FindFeaturesName(object, rownames(mean_exp), 
            "name")
        mean_exp <- cbind(Gene_ID = ChangeOUTName(rownames(mean_exp), 
            object@misc$fdata), Gene_name = Gene_name, mean_exp)
        WriteTable(mean_exp, file = outfile)
    }
}
CalPctExp <-
function (object, features = NULL, group.by = NULL, assay = NULL, 
    slot = "counts", is.return = FALSE, is.reverse = FALSE, is.bulk = FALSE, 
    outfile = "AllGene.pct_exp.xls") 
{
    data <- GetAssayData(object, assay = assay, slot = slot)
    if (!is.null(features)) 
        data <- data[features, ]
    if (!is.null(group.by)) 
        Idents(object) <- group.by
    mean_exp <- do.call(cbind, by(colnames(object), Idents(object), 
        function(x) {
            y <- if (is.reverse) {
                levels(x)[-as.numeric(x)]
            }
            else {
                x
            }
            Matrix::rowSums(data[, y, drop = F] > 0)/length(y)
        }))
    if (is.bulk) {
        mean_exp <- cbind(bulk = Matrix::rowSums(data > 0)/ncol(data), 
            mean_exp)
    }
    if (is.return) {
        return(mean_exp)
    }
    else {
        colnames(mean_exp) <- paste("Cluster", colnames(mean_exp))
        mean_exp <- cbind(Gene_ID = rownames(mean_exp), Gene_name = object@misc$fdata[rownames(mean_exp), 
            "name"], mean_exp)
        WriteTable(mean_exp, file = outfile)
    }
}
ChangeOUTName <-
function (features, fdata) 
{
    features <- as.character(features)
    fdata <- AddUnderscore(fdata)
    if (!is.null(fdata) && all(features %in% fdata$dash)) {
        underscore_id <- fdata$underscore
        names(underscore_id) <- fdata$dash
        features <- underscore_id[features]
    }
    return(features)
}
CheckPCA <-
function (object, reduction = "pca") 
{
    p1 <- DimPlot(object, reduction = reduction, group.by = "orig.ident")
    w <- 6
    if (exists("Phase", object@meta.data)) {
        p2 <- DimPlot(object, reduction = reduction, group.by = "Phase")
        p1 <- p1 + p2
        w <- 12
    }
    ggsave(p1, file = "pcaPlot.pdf", width = w, height = 6)
    dims <- min(20, ncol(Reductions(object, reduction)))
    p3 <- DimHeatmap(object, dims = seq(dims), cells = 500, balanced = TRUE, 
        ncol = 4, fast = FALSE, reduction = reduction)
    ggsave(p3, file = "pcaHeatmap.pdf", width = min(4, dims) * 
        4, height = ceiling(dims/4) * 4)
    if (!is.null(dev.list())) 
        dev.off()
    p4 <- ElbowPlot(object, ndims = ncol(Reductions(object, reduction)), 
        reduction = reduction)
    ggsave(p4, file = "pcaElbowPlot.pdf", width = 6, height = 6)
}
CheckVariableFeature <-
function (object) 
{
    top10 <- head(VariableFeatures(object), 10)
    plot1 <- VariableFeaturePlot(object)
    plot2 <- LabelPoints(plot = plot1, points = top10, labels = object@misc$fdata[top10, 
        "merge_name"], repel = TRUE, xnudge = 0, ynudge = 0)
    plot2 <- plot2 + theme(legend.position = "top")
    ggsave(plot2, file = "Variable_gene.pdf", width = 6, height = 6)
    write.table(VariableFeatures(object), file = "var_gene.xls", 
        quote = F, sep = "\t", col.name = F)
}
color.list <-
list(venn = list(default = list(`2` = c("#FFBDC0", "#C7D4EE"), 
    `3` = c("#69A4F9", "#FFCC66", "#FCB4B4"), `4` = c("#FFCC66", 
    "#BFE046", "#FCB4B4", "#69A4F9"), `5` = c("#BFE046", "#69A4F9", 
    "#ACB9EA", "#FCB4B4", "#FFCC66"), `6` = c("#FCB4B4", "#FFCC66", 
    "#BFE046", "#28C580", "#69A4F9", "#ACB9EA"), `7` = c("#FCB4B4", 
    "#FFCC66", "#BFE046", "#28C580", "#69A4F9", "#ACB9EA", "#C3C3C3"
    ), `0` = c("#FCB4B4", "#FFCC66", "#BFE046", "#28C580", "#69A4F9", 
    "#ACB9EA", "#C3C3C3"))), line = list(default = list(`1` = "#2771A7", 
    `2` = c("#D32421", "#2771A7"), `3` = c("#3A9736", "#2771A7", 
    "#D32421"), `4` = c("#3A9736", "#2771A7", "#C6AFD1", "#D32421"
    ), `5` = c("#3A9736", "#2771A7", "#C6AFD1", "#D32421", "#F3BB6F"
    ), `8` = c("#D32421", "#F09594", "#2771A7", "#3A9736", "#F3BB6F", 
    "#C6AFD1", "#831D20", "#A2C8DC"), `10` = c("#A2C8DC", "#F09594", 
    "#2771A7", "#C6AFD1", "#D32421", "#831D20", "#3A9736", "#F3BB6F", 
    "#A3A49E", "#5B4232"), `15` = c("#00468B", "#925E9F", "#759EDD", 
    "#0099B4", "#0A7C2E", "#B8D24D", "#EDE447", "#FAB158", "#FF7777", 
    "#FD0000", "#AD002A", "#AE8691", "#DEB8A1", "#CE9573", "#5B4232"
    ), `20` = c("#00468B", "#5377A7", "#6C6DA4", "#925E9F", "#759EDD", 
    "#0099B4", "#42C1BB", "#76D1B1", "#0A7C2E", "#B8D24D", "#EDE447", 
    "#FAB158", "#FDAF91", "#FF7777", "#FD0000", "#AD002A", "#AE8691", 
    "#DEB8A1", "#4C4E4E", "#5B4232"), `30` = c("#00468B", "#5377A7", 
    "#3B81AB", "#5C298F", "#6C6DA4", "#925E9F", "#759EDD", "#76C8DC", 
    "#0099B4", "#42C1BB", "#76D1B1", "#0F8074", "#0A7C2E", "#28AA6C", 
    "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", "#E67E74", "#FF7777", 
    "#FD0000", "#AD002A", "#792244", "#AD556B", "#AE8691", "#CE9573", 
    "#B09F91", "#DEB8A1", "#4C4E4E", "#5B4232"), `0` = c("#00468B", 
    "#5377A7", "#6C6DA4", "#925E9F", "#759EDD", "#0099B4", "#42C1BB", 
    "#76D1B1", "#0A7C2E", "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", 
    "#FF7777", "#FD0000", "#AD002A", "#AE8691", "#DEB8A1", "#4C4E4E", 
    "#5B4232"))), pie = list(default = list(`1` = "#42B540", 
    `2` = c("#42B540", "#EDE447"), `3` = c("#42B540", "#EDE447", 
    "#FF7777"), `4` = c("#00468B", "#42B540", "#EDE447", "#FF7777"
    ), `5` = c("#00468B", "#42B540", "#EDE447", "#759EDD", "#FF7777"
    ), `8` = c("#00468B", "#0099B4", "#76D1B1", "#42B540", "#EDE447", 
    "#FF7777", "#AD002A", "#759EDD"), `10` = c("#00468B", "#0099B4", 
    "#76D1B1", "#42B540", "#EDE447", "#FF7777", "#AD002A", "#759EDD", 
    "#DEB8A1", "#5B4232"), `20` = c("#00468B", "#5377A7", "#6C6DA4", 
    "#925E9F", "#759EDD", "#0099B4", "#42C1BB", "#76D1B1", "#42B540", 
    "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", "#FF7777", "#FD0000", 
    "#AD002A", "#AE8691", "#CE9573", "#DEB8A1", "#5B4232"), `30` = c("#00468B", 
    "#5377A7", "#3B81AB", "#5C298F", "#6C6DA4", "#925E9F", "#759EDD", 
    "#76C8DC", "#0099B4", "#42C1BB", "#76D1B1", "#0F8074", "#28AA6C", 
    "#42B540", "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", "#E67E74", 
    "#FF7777", "#FD0000", "#AD002A", "#792244", "#AD556B", "#AE8691", 
    "#CE9573", "#B09F91", "#756455", "#DEB8A1", "#5B4232"), `0` = c("#2e0a4a", 
    "#7a1b6c", "#15ad68", "#ded531", "#db9421", "#14b5b5", "#ede893", 
    "#76cfed", "#4599de", "#db2830", "#5F64AC", "#B271AD", "#eda4d3", 
    "#d33c67", "#EC6925", "#a155f9", "#70F2D3", "#6FBA33", "#EDAC2E", 
    "#096d42", "#4ec4b5", "#a36924", "#125fb2", "#7350EB", "#891a3a", 
    "#bf109a", "#E8851F", "#e77def", "#4ebee5", "#69a4f9", "#f9cfa5", 
    "#13D9B1", "#bfe046", "#DEB8A1", "#5B4232"))), twopie = list(
    default = list(`1` = "#2771A7", `2` = c("#D32421", "#2771A7"
    ), `3` = c("#3A9736", "#2771A7", "#D32421"), `4` = c("#3A9736", 
    "#2771A7", "#C6AFD1", "#D32421"), `5` = c("#3A9736", "#2771A7", 
    "#C6AFD1", "#D32421", "#F3BB6F"), `8` = c("#D32421", "#F09594", 
    "#2771A7", "#3A9736", "#F3BB6F", "#C6AFD1", "#831D20", "#A2C8DC"
    ), `10` = c("#A2C8DC", "#F09594", "#2771A7", "#C6AFD1", "#D32421", 
    "#831D20", "#3A9736", "#F3BB6F", "#A3A49E", "#5B4232"), `15` = c("#00468B", 
    "#925E9F", "#759EDD", "#0099B4", "#0A7C2E", "#B8D24D", "#EDE447", 
    "#FAB158", "#FF7777", "#FD0000", "#AD002A", "#AE8691", "#DEB8A1", 
    "#CE9573", "#5B4232"), `20` = c("#00468B", "#5377A7", "#6C6DA4", 
    "#925E9F", "#759EDD", "#0099B4", "#42C1BB", "#76D1B1", "#0A7C2E", 
    "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", "#FF7777", "#FD0000", 
    "#AD002A", "#AE8691", "#DEB8A1", "#4C4E4E", "#5B4232"), `30` = c("#00468B", 
    "#5377A7", "#3B81AB", "#5C298F", "#6C6DA4", "#925E9F", "#759EDD", 
    "#76C8DC", "#0099B4", "#42C1BB", "#76D1B1", "#0F8074", "#0A7C2E", 
    "#28AA6C", "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", "#E67E74", 
    "#FF7777", "#FD0000", "#AD002A", "#792244", "#AD556B", "#AE8691", 
    "#CE9573", "#B09F91", "#DEB8A1", "#4C4E4E", "#5B4232"), `0` = c("#00468B", 
    "#5377A7", "#6C6DA4", "#925E9F", "#759EDD", "#0099B4", "#42C1BB", 
    "#76D1B1", "#0A7C2E", "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", 
    "#FF7777", "#FD0000", "#AD002A", "#AE8691", "#DEB8A1", "#4C4E4E", 
    "#5B4232"))), bar = list(default = list(`1` = "#2771A7", 
    `2` = c("#D32421", "#2771A7"), `3` = c("#3A9736", "#2771A7", 
    "#D32421"), `4` = c("#3A9736", "#2771A7", "#C6AFD1", "#D32421"
    ), `5` = c("#3A9736", "#2771A7", "#C6AFD1", "#D32421", "#F3BB6F"
    ), `8` = c("#D32421", "#F09594", "#2771A7", "#3A9736", "#F3BB6F", 
    "#C6AFD1", "#831D20", "#A2C8DC"), `10` = c("#A2C8DC", "#F09594", 
    "#2771A7", "#C6AFD1", "#D32421", "#831D20", "#3A9736", "#F3BB6F", 
    "#A3A49E", "#5B4232"), `15` = c("#00468B", "#925E9F", "#759EDD", 
    "#0099B4", "#0A7C2E", "#B8D24D", "#EDE447", "#FAB158", "#FF7777", 
    "#FD0000", "#AD002A", "#AE8691", "#DEB8A1", "#CE9573", "#5B4232"
    ), `20` = c("#00468B", "#5377A7", "#6C6DA4", "#925E9F", "#759EDD", 
    "#0099B4", "#42C1BB", "#76D1B1", "#0A7C2E", "#B8D24D", "#EDE447", 
    "#FAB158", "#FDAF91", "#FF7777", "#FD0000", "#AD002A", "#AE8691", 
    "#DEB8A1", "#4C4E4E", "#5B4232"), `30` = c("#00468B", "#5377A7", 
    "#3B81AB", "#5C298F", "#6C6DA4", "#925E9F", "#759EDD", "#76C8DC", 
    "#0099B4", "#42C1BB", "#76D1B1", "#0F8074", "#0A7C2E", "#28AA6C", 
    "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", "#E67E74", "#FF7777", 
    "#FD0000", "#AD002A", "#792244", "#AD556B", "#AE8691", "#CE9573", 
    "#B09F91", "#DEB8A1", "#4C4E4E", "#5B4232"), `0` = c("#00468B", 
    "#5377A7", "#6C6DA4", "#925E9F", "#759EDD", "#0099B4", "#42C1BB", 
    "#76D1B1", "#0A7C2E", "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", 
    "#FF7777", "#FD0000", "#AD002A", "#AE8691", "#DEB8A1", "#4C4E4E", 
    "#5B4232")), stack = list(`1` = "#42B540", `2` = c("#42B540", 
"#EDE447"), `3` = c("#42B540", "#EDE447", "#FF7777"), `4` = c("#00468B", 
"#42B540", "#EDE447", "#FF7777"), `5` = c("#00468B", "#42B540", 
"#EDE447", "#759EDD", "#FF7777"), `8` = c("#00468B", "#0099B4", 
"#76D1B1", "#42B540", "#EDE447", "#FF7777", "#AD002A", "#759EDD"
), `10` = c("#00468B", "#0099B4", "#76D1B1", "#42B540", "#EDE447", 
"#FF7777", "#AD002A", "#759EDD", "#DEB8A1", "#5B4232"), `20` = c("#00468B", 
"#5377A7", "#6C6DA4", "#925E9F", "#759EDD", "#0099B4", "#42C1BB", 
"#76D1B1", "#42B540", "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", 
"#FF7777", "#FD0000", "#AD002A", "#AE8691", "#CE9573", "#DEB8A1", 
"#5B4232"), `30` = c("#00468B", "#5377A7", "#3B81AB", "#5C298F", 
"#6C6DA4", "#925E9F", "#759EDD", "#76C8DC", "#0099B4", "#42C1BB", 
"#76D1B1", "#0F8074", "#28AA6C", "#42B540", "#B8D24D", "#EDE447", 
"#FAB158", "#FDAF91", "#E67E74", "#FF7777", "#FD0000", "#AD002A", 
"#792244", "#AD556B", "#AE8691", "#CE9573", "#B09F91", "#756455", 
"#DEB8A1", "#5B4232"), `0` = c("#2e0a4a", "#7a1b6c", "#15ad68", 
"#ded531", "#db9421", "#14b5b5", "#ede893", "#76cfed", "#4599de", 
"#db2830", "#5F64AC", "#B271AD", "#eda4d3", "#d33c67", "#EC6925", 
"#a155f9", "#70F2D3", "#6FBA33", "#EDAC2E", "#096d42", "#4ec4b5", 
"#a36924", "#125fb2", "#7350EB", "#891a3a", "#bf109a", "#E8851F", 
"#e77def", "#4ebee5", "#69a4f9", "#f9cfa5", "#13D9B1", "#bfe046", 
"#DEB8A1", "#5B4232"))), hist = list(default = list(`1` = "#2771A7", 
    `2` = c("#D32421", "#2771A7"), `3` = c("#3A9736", "#2771A7", 
    "#D32421"), `4` = c("#3A9736", "#2771A7", "#C6AFD1", "#D32421"
    ), `5` = c("#3A9736", "#2771A7", "#C6AFD1", "#D32421", "#F3BB6F"
    ), `8` = c("#D32421", "#F09594", "#2771A7", "#3A9736", "#F3BB6F", 
    "#C6AFD1", "#831D20", "#A2C8DC"), `10` = c("#A2C8DC", "#F09594", 
    "#2771A7", "#C6AFD1", "#D32421", "#831D20", "#3A9736", "#F3BB6F", 
    "#A3A49E", "#5B4232"), `15` = c("#00468B", "#925E9F", "#759EDD", 
    "#0099B4", "#0A7C2E", "#B8D24D", "#EDE447", "#FAB158", "#FF7777", 
    "#FD0000", "#AD002A", "#AE8691", "#DEB8A1", "#CE9573", "#5B4232"
    ), `20` = c("#00468B", "#5377A7", "#6C6DA4", "#925E9F", "#759EDD", 
    "#0099B4", "#42C1BB", "#76D1B1", "#0A7C2E", "#B8D24D", "#EDE447", 
    "#FAB158", "#FDAF91", "#FF7777", "#FD0000", "#AD002A", "#AE8691", 
    "#DEB8A1", "#4C4E4E", "#5B4232"), `30` = c("#00468B", "#5377A7", 
    "#3B81AB", "#5C298F", "#6C6DA4", "#925E9F", "#759EDD", "#76C8DC", 
    "#0099B4", "#42C1BB", "#76D1B1", "#0F8074", "#0A7C2E", "#28AA6C", 
    "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", "#E67E74", "#FF7777", 
    "#FD0000", "#AD002A", "#792244", "#AD556B", "#AE8691", "#CE9573", 
    "#B09F91", "#DEB8A1", "#4C4E4E", "#5B4232"), `0` = c("#00468B", 
    "#5377A7", "#6C6DA4", "#925E9F", "#759EDD", "#0099B4", "#42C1BB", 
    "#76D1B1", "#0A7C2E", "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", 
    "#FF7777", "#FD0000", "#AD002A", "#AE8691", "#DEB8A1", "#4C4E4E", 
    "#5B4232"))), heatmap = list(default = list(`0` = c("#2F70AD", 
"#FFFFFF", "#BA2831"), `1` = c("#4785B6", "#FFFFFF", "#FF1717"
), `2` = c("#5B8089", "#FFFFFF", "#009933"), `3` = c("#5B8089", 
"#FDECBE", "#009933"), `4` = c("#3F99CB", "#FFFFFD", "#A1CD46"
), `5` = c("#FFFD36", "#FDFFFE", "#2C9B68"), `6` = c("#013F84", 
"#FDECBE", "#02908B")), corr = list(`0` = c("#FFFFFF", "#C4DEEC", 
"#2166AC"), `1` = c("#4785B6", "#FFFFFF", "#FF1717"), `2` = c("#5B8089", 
"#FFFFFF", "#009933"), `3` = c("#5B8089", "#FDECBE", "#009933"
), `4` = c("#3F99CB", "#FFFFFD", "#A1CD46"), `5` = c("#FFFD36", 
"#FDFFFE", "#2C9B68"), `6` = c("#013F84", "#FDECBE", "#02908B"
), `7` = c("#2F70AD", "#FFFFFF", "#BA2831"))), dot = list(default = list(
    `1` = "#2771A7", `2` = c("#D32421", "#2771A7"), `3` = c("#3A9736", 
    "#2771A7", "#D32421"), `4` = c("#3A9736", "#2771A7", "#C6AFD1", 
    "#D32421"), `5` = c("#3A9736", "#2771A7", "#C6AFD1", "#D32421", 
    "#F3BB6F"), `8` = c("#D32421", "#F09594", "#2771A7", "#3A9736", 
    "#F3BB6F", "#C6AFD1", "#831D20", "#A2C8DC"), `10` = c("#A2C8DC", 
    "#F09594", "#2771A7", "#C6AFD1", "#D32421", "#831D20", "#3A9736", 
    "#F3BB6F", "#A3A49E", "#5B4232"), `15` = c("#00468B", "#925E9F", 
    "#759EDD", "#0099B4", "#0A7C2E", "#B8D24D", "#EDE447", "#FAB158", 
    "#FF7777", "#FD0000", "#AD002A", "#AE8691", "#DEB8A1", "#CE9573", 
    "#5B4232"), `20` = c("#00468B", "#5377A7", "#6C6DA4", "#925E9F", 
    "#759EDD", "#0099B4", "#42C1BB", "#76D1B1", "#0A7C2E", "#B8D24D", 
    "#EDE447", "#FAB158", "#FDAF91", "#FF7777", "#FD0000", "#AD002A", 
    "#AE8691", "#DEB8A1", "#4C4E4E", "#5B4232"), `30` = c("#00468B", 
    "#5377A7", "#3B81AB", "#5C298F", "#6C6DA4", "#925E9F", "#759EDD", 
    "#76C8DC", "#0099B4", "#42C1BB", "#76D1B1", "#0F8074", "#0A7C2E", 
    "#28AA6C", "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", "#E67E74", 
    "#FF7777", "#FD0000", "#AD002A", "#792244", "#AD556B", "#AE8691", 
    "#CE9573", "#B09F91", "#DEB8A1", "#4C4E4E", "#5B4232"), `0` = c("#00468B", 
    "#5377A7", "#6C6DA4", "#925E9F", "#759EDD", "#0099B4", "#42C1BB", 
    "#76D1B1", "#0A7C2E", "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", 
    "#FF7777", "#FD0000", "#AD002A", "#AE8691", "#DEB8A1", "#4C4E4E", 
    "#5B4232")), volcano = list(`0` = c("#F9766D", "#00000032", 
"#609DFF"))), box = list(default = list(`1` = "#2771A7", `2` = c("#D32421", 
"#2771A7"), `3` = c("#3A9736", "#2771A7", "#D32421"), `4` = c("#3A9736", 
"#2771A7", "#C6AFD1", "#D32421"), `5` = c("#3A9736", "#2771A7", 
"#C6AFD1", "#D32421", "#F3BB6F"), `8` = c("#D32421", "#F09594", 
"#2771A7", "#3A9736", "#F3BB6F", "#C6AFD1", "#831D20", "#A2C8DC"
), `10` = c("#A2C8DC", "#F09594", "#2771A7", "#C6AFD1", "#D32421", 
"#831D20", "#3A9736", "#F3BB6F", "#A3A49E", "#5B4232"), `15` = c("#00468B", 
"#925E9F", "#759EDD", "#0099B4", "#0A7C2E", "#B8D24D", "#EDE447", 
"#FAB158", "#FF7777", "#FD0000", "#AD002A", "#AE8691", "#DEB8A1", 
"#CE9573", "#5B4232"), `20` = c("#00468B", "#5377A7", "#6C6DA4", 
"#925E9F", "#759EDD", "#0099B4", "#42C1BB", "#76D1B1", "#0A7C2E", 
"#B8D24D", "#EDE447", "#FAB158", "#FDAF91", "#FF7777", "#FD0000", 
"#AD002A", "#AE8691", "#DEB8A1", "#4C4E4E", "#5B4232"), `30` = c("#00468B", 
"#5377A7", "#3B81AB", "#5C298F", "#6C6DA4", "#925E9F", "#759EDD", 
"#76C8DC", "#0099B4", "#42C1BB", "#76D1B1", "#0F8074", "#0A7C2E", 
"#28AA6C", "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", "#E67E74", 
"#FF7777", "#FD0000", "#AD002A", "#792244", "#AD556B", "#AE8691", 
"#CE9573", "#B09F91", "#DEB8A1", "#4C4E4E", "#5B4232"), `0` = c("#00468B", 
"#5377A7", "#6C6DA4", "#925E9F", "#759EDD", "#0099B4", "#42C1BB", 
"#76D1B1", "#0A7C2E", "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", 
"#FF7777", "#FD0000", "#AD002A", "#AE8691", "#DEB8A1", "#4C4E4E", 
"#5B4232")), split = list(`0` = c("#F9766D", "#609DFF"), `1` = c("#56B4E9", 
"#E69F00"))), all = list(default = list(`1` = "#2771A7", `2` = c("#D32421", 
"#2771A7"), `3` = c("#3A9736", "#2771A7", "#D32421"), `4` = c("#3A9736", 
"#2771A7", "#C6AFD1", "#D32421"), `5` = c("#3A9736", "#2771A7", 
"#C6AFD1", "#D32421", "#F3BB6F"), `8` = c("#D32421", "#F09594", 
"#2771A7", "#3A9736", "#F3BB6F", "#C6AFD1", "#831D20", "#A2C8DC"
), `10` = c("#A2C8DC", "#F09594", "#2771A7", "#C6AFD1", "#D32421", 
"#831D20", "#3A9736", "#F3BB6F", "#A3A49E", "#5B4232"), `15` = c("#00468B", 
"#925E9F", "#759EDD", "#0099B4", "#0A7C2E", "#B8D24D", "#EDE447", 
"#FAB158", "#FF7777", "#FD0000", "#AD002A", "#AE8691", "#DEB8A1", 
"#CE9573", "#5B4232"), `20` = c("#00468B", "#5377A7", "#6C6DA4", 
"#925E9F", "#759EDD", "#0099B4", "#42C1BB", "#76D1B1", "#0A7C2E", 
"#B8D24D", "#EDE447", "#FAB158", "#FDAF91", "#FF7777", "#FD0000", 
"#AD002A", "#AE8691", "#DEB8A1", "#4C4E4E", "#5B4232"), `30` = c("#00468B", 
"#5377A7", "#3B81AB", "#5C298F", "#6C6DA4", "#925E9F", "#759EDD", 
"#76C8DC", "#0099B4", "#42C1BB", "#76D1B1", "#0F8074", "#0A7C2E", 
"#28AA6C", "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", "#E67E74", 
"#FF7777", "#FD0000", "#AD002A", "#792244", "#AD556B", "#AE8691", 
"#CE9573", "#B09F91", "#DEB8A1", "#4C4E4E", "#5B4232"), `0` = c("#00468B", 
"#5377A7", "#6C6DA4", "#925E9F", "#759EDD", "#0099B4", "#42C1BB", 
"#76D1B1", "#0A7C2E", "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", 
"#FF7777", "#FD0000", "#AD002A", "#AE8691", "#DEB8A1", "#4C4E4E", 
"#5B4232")), dark = list(`1` = "#0A7C2E", `2` = c("#D32421", 
"#00468B"), `3` = c("#0A7C2E", "#00468B", "#D32421"), `4` = c("#0A7C2E", 
"#00468B", "#925E9F", "#D32421"), `5` = c("#00468B", "#0A7C2E", 
"#0099B4", "#925E9F", "#AD002A"), `8` = c("#00468B", "#925E9F", 
"#0099B4", "#0A7C2E", "#FDAF91", "#FD0000", "#AD002A", "#DEB8A1"
), `10` = c("#00468B", "#925E9F", "#0099B4", "#3DB88C", "#0A7C2E", 
"#FDAF91", "#FD0000", "#AD002A", "#DEB8A1", "#5B4232"), `15` = c("#00468B", 
"#925E9F", "#759EDD", "#0099B4", "#0A7C2E", "#B8D24D", "#EDE447", 
"#FAB158", "#FF7777", "#FD0000", "#AD002A", "#AE8691", "#DEB8A1", 
"#CE9573", "#5B4232"), `20` = c("#00468B", "#5377A7", "#6C6DA4", 
"#925E9F", "#759EDD", "#0099B4", "#42C1BB", "#76D1B1", "#0A7C2E", 
"#B8D24D", "#EDE447", "#FAB158", "#FDAF91", "#FF7777", "#FD0000", 
"#AD002A", "#AE8691", "#DEB8A1", "#4C4E4E", "#5B4232"), `30` = c("#00468B", 
"#5377A7", "#3B81AB", "#5C298F", "#6C6DA4", "#925E9F", "#759EDD", 
"#76C8DC", "#0099B4", "#42C1BB", "#76D1B1", "#0F8074", "#0A7C2E", 
"#28AA6C", "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", "#E67E74", 
"#FF7777", "#FD0000", "#AD002A", "#792244", "#AD556B", "#AE8691", 
"#CE9573", "#B09F91", "#DEB8A1", "#4C4E4E", "#5B4232"), `0` = c("#00468B", 
"#5377A7", "#6C6DA4", "#925E9F", "#759EDD", "#0099B4", "#42C1BB", 
"#76D1B1", "#0A7C2E", "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", 
"#FF7777", "#FD0000", "#AD002A", "#AE8691", "#DEB8A1", "#4C4E4E", 
"#5B4232")), light = list(`1` = "#2771A7", `2` = c("#D32421", 
"#2771A7"), `3` = c("#3A9736", "#2771A7", "#D32421"), `4` = c("#3A9736", 
"#2771A7", "#C6AFD1", "#D32421"), `5` = c("#3A9736", "#2771A7", 
"#C6AFD1", "#D32421", "#F3BB6F"), `8` = c("#D32421", "#F09594", 
"#2771A7", "#3A9736", "#F3BB6F", "#C6AFD1", "#831D20", "#A2C8DC"
), `10` = c("#A2C8DC", "#F09594", "#2771A7", "#C6AFD1", "#D32421", 
"#831D20", "#3A9736", "#F3BB6F", "#A3A49E", "#5B4232"), `15` = c("#00468B", 
"#925E9F", "#759EDD", "#0099B4", "#0A7C2E", "#B8D24D", "#EDE447", 
"#FAB158", "#FF7777", "#FD0000", "#AD002A", "#AE8691", "#DEB8A1", 
"#CE9573", "#5B4232"), `20` = c("#00468B", "#5377A7", "#6C6DA4", 
"#925E9F", "#759EDD", "#0099B4", "#42C1BB", "#76D1B1", "#0A7C2E", 
"#B8D24D", "#EDE447", "#FAB158", "#FDAF91", "#FF7777", "#FD0000", 
"#AD002A", "#AE8691", "#DEB8A1", "#4C4E4E", "#5B4232"), `30` = c("#00468B", 
"#5377A7", "#3B81AB", "#5C298F", "#6C6DA4", "#925E9F", "#759EDD", 
"#76C8DC", "#0099B4", "#42C1BB", "#76D1B1", "#0F8074", "#0A7C2E", 
"#28AA6C", "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", "#E67E74", 
"#FF7777", "#FD0000", "#AD002A", "#792244", "#AD556B", "#AE8691", 
"#CE9573", "#B09F91", "#DEB8A1", "#4C4E4E", "#5B4232"), `0` = c("#00468B", 
"#5377A7", "#6C6DA4", "#925E9F", "#759EDD", "#0099B4", "#42C1BB", 
"#76D1B1", "#0A7C2E", "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", 
"#FF7777", "#FD0000", "#AD002A", "#AE8691", "#DEB8A1", "#4C4E4E", 
"#5B4232"))), gray = list(default = list(`2` = c("#D8DBE5", "#1B1919"
), `6` = c("#6C3C95", "#FDBF17", "#ABCF8E", "#EE2026", "#2E5696", 
"#C8C8CA"), `0` = c("#00468B", "#ED0000", "#42B540", "#0099B4", 
"#925E9F", "#FDAF91", "#AD002A", "#ADB6B6", "#1B1919"))), onlinereport = list(
    default = list(`0` = c("#CE2930", "#76C5E1", "#D28F26", "#7A1B6C", 
    "#21A465", "#D4CC38", "#2E0A4A", "#1FABAB", "#E6E191", "#458ECD", 
    "#891B20", "#3A322F", "#C8CAC9", "#5F64AB", "#B172AC", "#E09FC5", 
    "#D33C67", "#EB6924", "#A369F4", "#70F2D3", "#6EB933", "#EDAB2D", 
    "#096D42", "#4EC4B5", "#A36924", "#7350EB", "#891a3a", "#bf109a", 
    "#E8851F", "#e77def", "#4ebee5", "#69a4f9", "#f9cfa5", "#13D9B1", 
    "#bfe046", "#bfac88", "#d2f7ba", "#24e4ed", "#A6E2BE", "#c0e9fc", 
    "#f9efb4", "#fcb4b4", "#e5bcf7", "#acb9ea", "#a7dbd1", "#f2f225", 
    "#c67997", "#87a822", "#d39279", "#5DA6E8")), stack = list(
        `0` = c("#CE2930", "#76C5E1", "#D28F26", "#7A1B6C", "#21A465", 
        "#D4CC38", "#2E0A4A", "#1FABAB", "#E6E191", "#458ECD", 
        "#891B20", "#3A322F", "#C8CAC9", "#5F64AB", "#B172AC", 
        "#E09FC5", "#D33C67", "#EB6924", "#A369F4", "#70F2D3", 
        "#6EB933", "#EDAB2D", "#096D42", "#4EC4B5", "#A36924", 
        "#7350EB", "#891a3a", "#bf109a", "#E8851F", "#e77def", 
        "#4ebee5", "#69a4f9", "#f9cfa5", "#13D9B1", "#bfe046", 
        "#bfac88", "#d2f7ba", "#24e4ed", "#A6E2BE", "#c0e9fc", 
        "#f9efb4", "#fcb4b4", "#e5bcf7", "#acb9ea", "#a7dbd1", 
        "#f2f225", "#c67997", "#87a822", "#d39279", "#5DA6E8"
        )), set1 = list(`0` = c("#F0B142", "#62D6F6", "#F14149", 
    "#28C580", "#FBF73E", "#B02F96", "#1CE8BF", "#F3854F", "#561BB5", 
    "#C6C5E8", "#7ED8FA", "#FDF791", "#EFD53F", "#E46665", "#4846A6", 
    "#BEDF50", "#46BBA8", "#acb9ea", "#A17817", "#821762", "#A2E2B3", 
    "#406BDE", "#0CB1ED", "#51EBFC", "#c4f98c", "#FFDC48", "#FEA8C3", 
    "#F1835B", "#C4477A", "#821762", "#D544AE", "#50B5F5", "#90D1F4", 
    "#A0F5F1", "#D1FCE4", "#54C4B5", "#D4E746", "#EDB349", "#FEE05B", 
    "#F76C50", "#BA2D2D", "#3793DE", "#a7dbd1", "#ede893", "#76cfed", 
    "#c67997", "#bfac88", "#d39279", "#9e516b", "#666b96")), 
    set2 = list(`0` = c("#69a4f9", "#bfe046", "#EDAC2E", "#f2f225", 
    "#a155f9", "#db2830", "#4ec4b5", "#4599de", "#E8851F", "#6FBA33", 
    "#5F64AC", "#58B984", "#891a3a", "#bf109a", "#f7cb16", "#15ad68", 
    "#14b5b5", "#4ebee5", "#1a497c", "#EC6925", "#096d42", "#7a1b6c", 
    "#d33c67", "#db9421", "#ded531", "#521f84", "#87a822", "#43efce", 
    "#11C2CC", "#fc68c7", "#e77def", "#76cfed", "#125fb2", "#a36924", 
    "#c67997", "#B271AD", "#acb9ea", "#fcb4b4", "#f9cfa5", "#f9efb4", 
    "#d2f7ba", "#8bd6d6", "#c5fcdd", "#c0e9fc", "#e5bcf7", "#eda4d3", 
    "#d39279", "#d6beb2", "#e2a1ba", "#7350EB")), set3 = list(
        `0` = c("#bfe046", "#f7cb16", "#fcb4b4", "#e77def", "#acb9ea", 
        "#24e4ed", "#69a4f9", "#43efce", "#14b5b5", "#125fb2", 
        "#521f84", "#db2830", "#EC6925", "#f2f225", "#c0e9fc", 
        "#d2f7ba", "#f9cfa5", "#eda4d3", "#15ad68", "#6FBA33", 
        "#096d42", "#db9421", "#1a497c", "#7350EB", "#bf109a", 
        "#d33c67", "#fc68c7", "#EDAC2E", "#ded531", "#87a822", 
        "#4ec4b5", "#77D7F2", "#a155f9", "#B271AD", "#7a1b6c", 
        "#a36924", "#c67997", "#E8851F", "#891a3a", "#d39279", 
        "#e5bcf7", "#5F64AC", "#4599de", "#76cfed", "#e2a1ba", 
        "#f9efb4", "#d6beb2", "#58B984", "#c5fcdd", "#8bd6d6"
        ))), tsne = list(default = list(`1` = "#42B540", `2` = c("#42B540", 
"#EDE447"), `3` = c("#42B540", "#EDE447", "#FF7777"), `4` = c("#00468B", 
"#42B540", "#EDE447", "#FF7777"), `5` = c("#00468B", "#42B540", 
"#EDE447", "#759EDD", "#FF7777"), `10` = c("#00468B", "#925E9F", 
"#0099B4", "#76D1B1", "#42B540", "#EDE447", "#FAB158", "#FF7777", 
"#AD002A", "#759EDD"), `15` = c("#00468B", "#925E9F", "#759EDD", 
"#0099B4", "#76D1B1", "#42B540", "#B8D24D", "#EDE447", "#FAB158", 
"#FF7777", "#FD0000", "#AD002A", "#AE8691", "#CE9573", "#756455"
), `20` = c("#00468B", "#5377A7", "#6C6DA4", "#925E9F", "#759EDD", 
"#0099B4", "#42C1BB", "#76D1B1", "#42B540", "#B8D24D", "#EDE447", 
"#FAB158", "#FDAF91", "#FF7777", "#FD0000", "#AD002A", "#792244", 
"#AE8691", "#CE9573", "#756455"), `24` = c("#00468B", "#5377A7", 
"#3B81AB", "#5C298F", "#6C6DA4", "#925E9F", "#759EDD", "#76C8DC", 
"#0099B4", "#42C1BB", "#76D1B1", "#0F8074", "#42B540", "#B8D24D", 
"#EDE447", "#FAB158", "#FDAF91", "#FF7777", "#FD0000", "#AD002A", 
"#792244", "#AE8691", "#CE9573", "#756455"), `28` = c("#00468B", 
"#5377A7", "#3B81AB", "#5C298F", "#6C6DA4", "#925E9F", "#AD729B", 
"#759EDD", "#76C8DC", "#0099B4", "#42C1BB", "#76D1B1", "#0F8074", 
"#28AA6C", "#42B540", "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", 
"#E67E74", "#FF7777", "#FD0000", "#AD002A", "#792244", "#AE8691", 
"#CE9573", "#B09F91", "#756455"), `35` = c("#00468B", "#5377A7", 
"#3B81AB", "#5C298F", "#6C6DA4", "#925E9F", "#AD729B", "#BF8099", 
"#D18E96", "#759EDD", "#76C8DC", "#0099B4", "#42C1BB", "#0F8074", 
"#438424", "#28AA6C", "#42B540", "#B8D24D", "#EDE447", "#FCCD94", 
"#FAB158", "#DCA30C", "#DD7C06", "#E39C94", "#FDAF91", "#E67E74", 
"#FF7777", "#FD0000", "#AD002A", "#792244", "#AF556B", "#AE8691", 
"#CE9573", "#B09F91", "#756455"), `0` = c("#00468B", "#5377A7", 
"#3B81AB", "#5C298F", "#6C6DA4", "#925E9F", "#AD729B", "#BF8099", 
"#D18E96", "#759EDD", "#76C8DC", "#0099B4", "#42C1BB", "#0F8074", 
"#438424", "#28AA6C", "#42B540", "#B8D24D", "#EDE447", "#FCCD94", 
"#FAB158", "#DCA30C", "#DD7C06", "#E39C94", "#FDAF91", "#E67E74", 
"#FF7777", "#FD0000", "#AD002A", "#792244", "#AF556B", "#AE8691", 
"#CE9573", "#B09F91", "#756455")), set1 = list(`1` = "#42B540", 
    `2` = c("#42B540", "#EDE447"), `3` = c("#42B540", "#EDE447", 
    "#FF7777"), `4` = c("#00468B", "#42B540", "#EDE447", "#FF7777"
    ), `5` = c("#00468B", "#42B540", "#EDE447", "#759EDD", "#FF7777"
    ), `10` = c("#00468B", "#925E9F", "#0099B4", "#76D1B1", "#42B540", 
    "#EDE447", "#FAB158", "#FF7777", "#AD002A", "#759EDD"), `15` = c("#00468B", 
    "#925E9F", "#759EDD", "#0099B4", "#76D1B1", "#42B540", "#B8D24D", 
    "#EDE447", "#FAB158", "#FF7777", "#FD0000", "#AD002A", "#AE8691", 
    "#CE9573", "#756455"), `20` = c("#00468B", "#5377A7", "#6C6DA4", 
    "#925E9F", "#759EDD", "#0099B4", "#42C1BB", "#76D1B1", "#42B540", 
    "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", "#FF7777", "#FD0000", 
    "#AD002A", "#792244", "#AE8691", "#CE9573", "#756455"), `24` = c("#00468B", 
    "#5377A7", "#3B81AB", "#5C298F", "#6C6DA4", "#925E9F", "#759EDD", 
    "#76C8DC", "#0099B4", "#42C1BB", "#76D1B1", "#0F8074", "#42B540", 
    "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", "#FF7777", "#FD0000", 
    "#AD002A", "#792244", "#AE8691", "#CE9573", "#756455"), `28` = c("#00468B", 
    "#5377A7", "#3B81AB", "#5C298F", "#6C6DA4", "#925E9F", "#AD729B", 
    "#759EDD", "#76C8DC", "#0099B4", "#42C1BB", "#76D1B1", "#0F8074", 
    "#28AA6C", "#42B540", "#B8D24D", "#EDE447", "#FAB158", "#FDAF91", 
    "#E67E74", "#FF7777", "#FD0000", "#AD002A", "#792244", "#AE8691", 
    "#CE9573", "#B09F91", "#756455"), `35` = c("#00468B", "#5377A7", 
    "#3B81AB", "#5C298F", "#6C6DA4", "#925E9F", "#AD729B", "#BF8099", 
    "#D18E96", "#759EDD", "#76C8DC", "#0099B4", "#42C1BB", "#0F8074", 
    "#438424", "#28AA6C", "#42B540", "#B8D24D", "#EDE447", "#FCCD94", 
    "#FAB158", "#DCA30C", "#DD7C06", "#E39C94", "#FDAF91", "#E67E74", 
    "#FF7777", "#FD0000", "#AD002A", "#792244", "#AF556B", "#AE8691", 
    "#CE9573", "#B09F91", "#756455"), `0` = c("#00468B", "#5377A7", 
    "#3B81AB", "#5C298F", "#6C6DA4", "#925E9F", "#AD729B", "#BF8099", 
    "#D18E96", "#759EDD", "#76C8DC", "#0099B4", "#42C1BB", "#0F8074", 
    "#438424", "#28AA6C", "#42B540", "#B8D24D", "#EDE447", "#FCCD94", 
    "#FAB158", "#DCA30C", "#DD7C06", "#E39C94", "#FDAF91", "#E67E74", 
    "#FF7777", "#FD0000", "#AD002A", "#792244", "#AF556B", "#AE8691", 
    "#CE9573", "#B09F91", "#756455")), set2 = list(`1` = "#94B2CE", 
    `2` = c("#94B2CE", "#EBB196"), `3` = c("#94B2CE", "#EBB196", 
    "#AACDA2"), `4` = c("#B6D1FC", "#EBB196", "#AACDA2", "#E2BDD4"
    ), `5` = c("#94B2CE", "#EBB196", "#AACDA2", "#E2BDD4", "#AB91BB"
    ), `10` = c("#94B2CE", "#AEB7C7", "#AACDA2", "#DBB3B5", "#82BDE1", 
    "#E2BDD4", "#BBB978", "#CDBBDB", "#EBB196", "#AB91BB"), `16` = c("#94B2CE", 
    "#AEB7C7", "#AACDA2", "#DBB3B5", "#82BDE1", "#887474", "#E2BDD4", 
    "#BBB978", "#8AA5A2", "#DDB177", "#8892BF", "#DC7FAF", "#78AB91", 
    "#CDBBDB", "#EBB196", "#AB91BB"), `20` = c("#94B2CE", "#AEB7C7", 
    "#9B8EAC", "#AACDA2", "#C0767A", "#DBB3B5", "#82BDE1", "#778B8B", 
    "#887474", "#6D6D6D", "#BBB978", "#8AA5A2", "#DDB177", "#8892BF", 
    "#DC7FAF", "#78AB91", "#CDBBDB", "#86A0C7", "#EBB196", "#AB91BB"
    ), `25` = c("#94B2CE", "#AEB7C7", "#9B8EAC", "#AACDA2", "#C0767A", 
    "#DBB3B5", "#82BDE1", "#778B8B", "#887474", "#B0A9A5", "#E2BDD4", 
    "#6D6D6D", "#BBB978", "#8AA5A2", "#DDB177", "#8892BF", "#C48ABD", 
    "#DC7FAF", "#78AB91", "#919E73", "#CDBBDB", "#C08A67", "#86A0C7", 
    "#EBB196", "#AB91BB"), `30` = c("#94B2CE", "#AEB7C7", "#D0AE90", 
    "#9B8EAC", "#AACDA2", "#C0767A", "#DBB3B5", "#82BDE1", "#778B8B", 
    "#887474", "#B0A9A5", "#E2BDD4", "#FFADC9", "#6D6D6D", "#BBB978", 
    "#C5C6A7", "#8AA5A2", "#A9C3C7", "#DDB177", "#8892BF", "#C48ABD", 
    "#B6CB95", "#DC7FAF", "#78AB91", "#919E73", "#CDBBDB", "#C08A67", 
    "#86A0C7", "#EBB196", "#AB91BB"), `35` = c("#799FC2", "#AEB7C7", 
    "#C4B29A", "#D0AE90", "#9B8EAC", "#AACDA2", "#C0767A", "#D2A0A3", 
    "#63ACD9", "#B9B1C0", "#887474", "#B0A9A5", "#DBACC9", "#FFADC9", 
    "#6D6D6D", "#B7B7B9", "#BBB978", "#C5C6A7", "#8AA5A2", "#A9C3C7", 
    "#DDB177", "#8892BF", "#B1B076", "#C48ABD", "#B6CB95", "#DC7FAF", 
    "#78AB91", "#CB7C88", "#919E73", "#CDBBDB", "#C08A67", "#86A0C7", 
    "#EBB196", "#AB91BB", "#DAC097"), `0` = c("#799FC2", "#AEB7C7", 
    "#C4B29A", "#D0AE90", "#9B8EAC", "#AACDA2", "#C0767A", "#D2A0A3", 
    "#63ACD9", "#B9B1C0", "#887474", "#B0A9A5", "#DBACC9", "#FFADC9", 
    "#6D6D6D", "#B7B7B9", "#BBB978", "#C5C6A7", "#8AA5A2", "#A9C3C7", 
    "#DDB177", "#8892BF", "#B1B076", "#C48ABD", "#B6CB95", "#DC7FAF", 
    "#78AB91", "#CB7C88", "#919E73", "#CDBBDB", "#C08A67", "#86A0C7", 
    "#EBB196", "#AB91BB", "#DAC097")), set3 = list(`1` = "#E4A436", 
    `2` = c("#E4A436", "#48B9D0"), `3` = c("#E4A436", "#48B9D0", 
    "#A668AC"), `4` = c("#E4A436", "#48B9D0", "#A668AC", "#DFA7A4"
    ), `5` = c("#E4A436", "#48B9D0", "#A668AC", "#8EC5BF", "#DFA7A4"
    ), `10` = c("#E4A436", "#BD685F", "#48B9D0", "#83A5CD", "#A668AC", 
    "#815A6C", "#DFA7A4", "#8EC5BF", "#8E745D", "#CFC4B8"), `14` = c("#E4A436", 
    "#BD685F", "#48B9D0", "#94DBE0", "#83A5CD", "#A668AC", "#815A6C", 
    "#ED987F", "#8492B4", "#8EC5BF", "#D31715", "#835B43", "#A28C74", 
    "#D5CFAA"), `20` = c("#E4A436", "#BD685F", "#6CA3B4", "#48B9D0", 
    "#94DBE0", "#7CBA92", "#83A5CD", "#A668AC", "#815A6C", "#D6B5B1", 
    "#ED987F", "#703225", "#8492B4", "#8EC5BF", "#AB877D", "#D31715", 
    "#B5271D", "#835B43", "#A28C74", "#D5CFAA"), `25` = c("#E4A436", 
    "#BD685F", "#948689", "#6CA3B4", "#48B9D0", "#94DBE0", "#22ACA8", 
    "#66A39D", "#7CBA92", "#2073B7", "#83A5CD", "#A668AC", "#D6B5B1", 
    "#ED987F", "#703225", "#8492B4", "#87A3B8", "#8EC5BF", "#AB877D", 
    "#D31715", "#B5271D", "#835B43", "#CFC4B8", "#A28C74", "#D5CFAA"
    ), `30` = c("#E4A436", "#BD685F", "#948689", "#6CA3B4", "#48B9D0", 
    "#94DBE0", "#22ACA8", "#66A39D", "#7CBA92", "#2073B7", "#83A5CD", 
    "#A668AC", "#815A6C", "#D6B5B1", "#ED987F", "#703225", "#B37C79", 
    "#8492B4", "#87A3B8", "#8EC5BF", "#99BFB2", "#AB877D", "#DFA7A4", 
    "#CE0D09", "#B5271D", "#9C4130", "#8E745D", "#CFC4B8", "#A28C74", 
    "#D5CFAA"), `35` = c("#E4A436", "#BD685F", "#948689", "#6CA3B4", 
    "#48B9D0", "#94DBE0", "#22ACA8", "#66A39D", "#1B766C", "#7CBA92", 
    "#2073B7", "#83A5CD", "#A668AC", "#815A6C", "#D6B5B1", "#ED987F", 
    "#703225", "#B37C79", "#9E93A7", "#8492B4", "#87A3B8", "#8BB4BB", 
    "#8EC5BF", "#99BFB2", "#AB877D", "#DFA7A4", "#D31715", "#CE0D09", 
    "#B5271D", "#9C4130", "#835B43", "#8E745D", "#CFC4B8", "#A28C74", 
    "#D5CFAA"), `0` = c("#E4A436", "#BD685F", "#948689", "#6CA3B4", 
    "#48B9D0", "#94DBE0", "#22ACA8", "#66A39D", "#1B766C", "#7CBA92", 
    "#2073B7", "#83A5CD", "#A668AC", "#815A6C", "#D6B5B1", "#ED987F", 
    "#703225", "#B37C79", "#9E93A7", "#8492B4", "#87A3B8", "#8BB4BB", 
    "#8EC5BF", "#99BFB2", "#AB877D", "#DFA7A4", "#D31715", "#CE0D09", 
    "#B5271D", "#9C4130", "#835B43", "#8E745D", "#CFC4B8", "#A28C74", 
    "#D5CFAA"))), ggsci = list(npg = list(`0` = c("#E64B35", 
"#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", 
"#DC0000", "#7E6148", "#B09C85")), aaas = list(`0` = c("#3B4992", 
"#EE0000", "#008B45", "#631879", "#008280", "#BB0021", "#5F559B", 
"#A20056", "#808180", "#1B1919")), nejm = list(`0` = c("#BC3C29", 
"#0072B5", "#E18727", "#20854E", "#7876B1", "#6F99AD", "#FFDC91", 
"#EE4C97")), lancet = list(`0` = c("#00468B", "#ED0000", "#42B540", 
"#0099B4", "#925E9F", "#FDAF91", "#AD002A", "#ADB6B6", "#1B1919"
)), jama = list(`0` = c("#374E55", "#DF8F44", "#00A1D5", "#B24745", 
"#79AF97", "#6A6599", "#80796B")), jco = list(`0` = c("#0073C2", 
"#EFC000", "#868686", "#CD534C", "#7AA6DC", "#003C67", "#8F7700", 
"#3B3B3B", "#A73030", "#4A6990")), ucscgb = list(`0` = c("#FF0000", 
"#FF9900", "#FFCC00", "#00FF00", "#6699FF", "#CC33FF", "#99991E", 
"#999999", "#FF00CC", "#CC0000", "#FFCCCC", "#FFFF00", "#CCFF00", 
"#358000", "#0000CC", "#99CCFF", "#00FFFF", "#CCFFFF", "#9900CC", 
"#CC99FF", "#996600", "#666600", "#666666", "#CCCCCC", "#79CC3D", 
"#CCCC99")), d3_1 = list(`0` = c("#1F77B4", "#FF7F0E", "#2CA02C", 
"#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", 
"#17BECF")), d3_2 = list(`0` = c("#1F77B4", "#FF7F0E", "#2CA02C", 
"#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", 
"#17BECF", "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5", 
"#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5")), locuszoom = list(
    `0` = c("#D43F3A", "#EEA236", "#5CB85C", "#46B8DA", "#357EBD", 
    "#9632B8", "#B8B8B8")), igv = list(`0` = c("#5050FF", "#CE3D32", 
"#749B58", "#F0E685", "#466983", "#BA6338", "#5DB1DD", "#802268", 
"#6BD76B", "#D595A7", "#924822", "#837B8D", "#C75127", "#D58F5C", 
"#7A65A5", "#E4AF69", "#3B1B53", "#CDDEB7", "#612A79", "#AE1F63", 
"#E7C76F", "#5A655E", "#CC9900", "#99CC00", "#A9A9A9", "#CC9900", 
"#99CC00", "#33CC00", "#00CC33", "#00CC99", "#0099CC", "#0A47FF", 
"#4775FF", "#FFC20A", "#FFD147", "#990033", "#991A00", "#996600", 
"#809900", "#339900", "#00991A", "#009966", "#008099", "#003399", 
"#1A0099", "#660099", "#990080", "#D60047", "#FF1463", "#00D68F", 
"#14FFB1")), cosmic_1 = list(`0` = c("#2E2A2B", "#CF4E9C", "#8C57A2", 
"#358DB9", "#82581F", "#2F509E", "#E5614C", "#97A1A7", "#3DA873", 
"#DC9445")), cosmic_2 = list(`0` = c("#171717", "#7D0226", "#300049", 
"#165459", "#3F2327", "#0B1948", "#E71012", "#555555", "#193006", 
"#A8450C")), cosmic_3 = list(`0` = c("#5ABCEB", "#050708", "#D33C32", 
"#CBCACB", "#ABCD72", "#E7C9C6")), uchicago_1 = list(`0` = c("#800000", 
"#767676", "#FFA319", "#8A9045", "#155F83", "#C16622", "#8F3931", 
"#58593F", "#350E20")), uchicago_2 = list(`0` = c("#800000", 
"#D6D6CE", "#FFB547", "#ADB17D", "#5B8FA8", "#D49464", "#B1746F", 
"#8A8B79", "#725663")), uchicago_3 = list(`0` = c("#800000", 
"#767676", "#CC8214", "#616530", "#0F425C", "#9A5324", "#642822", 
"#3E3E23", "#350E20")), startrek = list(`0` = c("#CC0C00", "#5C88DA", 
"#84BD00", "#FFCD00", "#7C878E", "#00B5E2", "#00AF66")), tron = list(
    `0` = c("#FF410D", "#6EE2FF", "#F7C530", "#95CC5E", "#D0DFE6", 
    "#F79D1E", "#748AA6")), futurama = list(`0` = c("#FF6F00", 
"#C71000", "#008EA0", "#8A4198", "#5A9599", "#FF6348", "#84D7E1", 
"#FF95A8", "#3D3B25", "#ADE2D0", "#1A5354", "#3F4041")), rickandmorty = list(
    `0` = c("#FAFD7C", "#82491E", "#24325F", "#B7E4F9", "#FB6467", 
    "#526E2D", "#E762D7", "#E89242", "#FAE48B", "#A6EEE6", "#917C5D", 
    "#69C8EC")), simpsons = list(`0` = c("#FED439", "#709AE1", 
"#8A9197", "#D2AF81", "#FD7446", "#D5E4A2", "#197EC0", "#F05C3B", 
"#46732E", "#71D0F5", "#370335", "#075149", "#C80813", "#91331F", 
"#1A9993", "#FD8CC1")), gsea = list(`0` = c("#4500AC", "#2600D1", 
"#6B58EE", "#8787FF", "#C6C0FF", "#D4D4FF", "#FFBFE5", "#FF8888", 
"#FF707F", "#FF5959", "#EE3F3F", "#D60C00"))))
color.yaml <-
"Colors.yml"
core_otsu <-
function (dt) 
{
    weight.bg <- sapply(seq(0, max(dt)), function(i) sum(dt <= 
        i)/length(dt))
    weight.fg <- 1 - weight.bg
    mean.bg <- sapply(seq(0, max(dt)), function(i) mean(dt[dt <= 
        i]))
    mean.fg <- sapply(seq(0, max(dt)), function(i) mean(dt[dt > 
        i]))
    var.between <- weight.bg * weight.fg * (mean.bg - mean.fg)^2
    thres <- which.max(var.between)
    return(thres)
}
DoCellCycleScoring <-
function (object, s.genes = NULL, g2m.genes = NULL, ...) 
{
    if (is.null(s.genes)) 
        s.genes <- cc.genes$s.genes
    if (is.null(g2m.genes)) 
        g2m.genes <- cc.genes$g2m.genes
    s.genes <- FindFeaturesID(object, s.genes)
    g2m.genes <- FindFeaturesID(object, g2m.genes)
    if (length(s.genes) < 2 || length(g2m.genes) < 2) 
        return(object)
    object <- CellCycleScoring(object, s.features = s.genes, 
        g2m.features = g2m.genes, set.ident = FALSE, ...)
    object[["CC.Difference"]] <- object[["S.Score"]] - object[["G2M.Score"]]
    return(object)
}
DoDimReduc <-
function (object, assay = NULL, pc.num = 50, is.checkpca = TRUE, 
    ...) 
{
    if (is.null(assay)) 
        assay <- DefaultAssay(object)
    message("-->PCA<--")
    object <- RunPCA(object = object, assay = assay, npcs = pc.num, 
        features = VariableFeatures(object), verbose = FALSE)
    object[[paste0("pca_", assay)]] <- object[["pca"]]
    if (is.checkpca) 
        CheckPCA(object)
    sig.PCs <- seq(pc.num)
    message("-->Run tSNE<--")
    object <- RunTSNE(object, dims = sig.PCs, ...)
    object[[paste0("tsne_", assay)]] <- object[["tsne"]]
    message("-->Run UMAP<--")
    object <- RunUMAP(object, dims = sig.PCs, umap.method = "uwot")
    object[[paste0("umap_", assay)]] <- object[["umap"]]
    return(object)
}
DoFindAllMarkers <-
function (object, parameter = list(), group.by = "seurat_clusters", 
    min.pct = IfNull(parameter$FindMarkers$min_pct, 0.25), logfc.threshold = IfNull(parameter$FindMarkers$logfc, 
        0.25), return.thresh = IfNull(parameter$FindMarkers$pvalue, 
        0.01)) 
{
    Idents(object) <- group.by
    object.markers <- FindAllMarkers(object = object, only.pos = TRUE, 
        min.pct = min.pct, logfc.threshold = logfc.threshold, 
        return.thresh = return.thresh, pseudocount.use = 0)
    return(object.markers)
}
DoFindClusters <-
function (object, reduction = "pca", dims = NULL, resolution = 0.5) 
{
    if (is.null(dims) || max(dims) > length(object[[reduction]])) 
        dims <- seq(object[[reduction]])
    object <- FindNeighbors(object, reduction = reduction, dims = dims, 
        force.recalc = TRUE)
    object <- FindClusters(object, resolution = resolution, temp.file.location = getwd())
    color.cluster <- fetch_color(nlevels(object@meta.data$seurat_clusters), 
        "tsne", "set1")
    names(color.cluster) <- levels(object@meta.data$seurat_clusters)
    object@misc[["color.cluster"]] <- color.cluster
    return(object)
}
DoIntegration <-
function (object, split.by = "orig.ident", dims = 1:50, nfeatures = 3000, 
    is.SCT = FALSE) 
{
    old.assay <- DefaultAssay(object)
    object.list <- SplitObject(object, split.by = split.by)
    object.list <- SplitObject.Image(object.list)
    anchor.features <- nfeatures
    normalization.method <- "LogNormalize"
    if (is.SCT) {
        anchor.features <- SelectIntegrationFeatures(object.list = object.list, 
            nfeatures = nfeatures)
        object.list <- PrepSCTIntegration(object.list = object.list, 
            anchor.features = anchor.features)
        normalization.method <- "SCT"
    }
    k.filter <- min(200, ceiling(min(sapply(object.list, ncol))/2))
    anchors <- FindIntegrationAnchors(object.list = object.list, 
        dims = dims, normalization.method = normalization.method, 
        anchor.features = anchor.features, k.filter = k.filter)
    integrated <- IntegrateData(anchorset = anchors, dims = dims, 
        normalization.method = normalization.method)
    integrated <- ScaleData(integrated, verbose = FALSE)
    integrated@misc <- object@misc
    integrated[[old.assay]] <- object[[old.assay]]
    integrated@reductions <- object@reductions
    integrated@meta.data[[split.by]] <- factor(integrated@meta.data[[split.by]], 
        levels = levels(object@meta.data[[split.by]]))
    return(integrated)
}
DoNormalization <-
function (object, parameter = list(), assay = "RNA", is_SCTransform = FALSE, 
    vars.regress = NULL, is.check = TRUE, ...) 
{
    DefaultAssay(object) <- assay
    object <- FindRegressVars(object, parameter, vars.regress = vars.regress, 
        ...)
    if (is_SCTransform) {
        object <- SCTransform(object, assay = assay, vars.to.regress = object@misc$vars.regress, 
            verbose = FALSE, min_cells = 1, return.only.var.genes = F)
    }
    else {
        message("-->Normalize Data<--")
        object <- NormalizeData(object, normalization.method = "LogNormalize", 
            scale.factor = 10000)
        message("-->Find Variable Genes<--")
        object <- FindVariableFeatures(object, selection.method = "vst", 
            nfeatures = 2000)
        if (is.check) 
            CheckVariableFeature(object)
        message("-->ScaleData<--")
        object <- ScaleData(object = object, vars.to.regress = object@misc$vars.regress, 
            features = rownames(object))
    }
    return(object)
}
dot_theme_default <-
function () 
{
    library(ggplot2)
    mytheme <- theme_bw() + theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.border = element_rect(color = "#000000", 
            size = 1), axis.text = element_text(color = "#000000", 
            size = 14), axis.text.x = element_text(angle = 0, 
            hjust = 0.5, vjust = 0.5), axis.text.y = element_text(hjust = 0.5, 
            vjust = 0.5), axis.title = element_text(color = "#000000", 
            size = 16, face = "bold"), axis.title.x = element_text(margin = margin(2.5, 
            0, 2.5, 0, "mm")), axis.title.y = element_text(margin = margin(0, 
            2.5, 0, 2.5, "mm")), axis.ticks = element_line(color = "#000000", 
            size = 0.6), axis.ticks.length = unit(0.13, "cm"), 
        legend.title = element_blank(), legend.text = element_text(size = 12), 
        plot.title = element_text(color = "#000000", size = 20, 
            face = "bold", hjust = 0.5), )
    if ("showtext" %in% installed.packages()) {
        library(showtext)
        showtext_auto(enable = TRUE)
        font_add("Arial", regular = "/Bio/Bin/pipeline/GeneralPlot/v1.0/fonts/tff/msttcore/arial.ttf")
        mytheme <- mytheme + theme(text = element_text(family = "Arial"))
    }
    mytheme
}
dot_theme_manhattan <-
function () 
{
    library(ggplot2)
    mytheme <- theme_bw() + theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.border = element_rect(color = "#000000", 
            size = 1), axis.text = element_text(color = "#000000", 
            size = 12), axis.text.x = element_text(angle = 0, 
            hjust = 0.5, vjust = 0.5), axis.text.y = element_text(hjust = 1, 
            vjust = 0.5), axis.title = element_text(color = "#000000", 
            size = 16, face = "bold"), axis.title.x = element_text(margin = margin(2.5, 
            0, 2.5, 0, "mm")), axis.title.y = element_text(margin = margin(0, 
            2.5, 0, 2.5, "mm")), axis.ticks = element_line(color = "#000000", 
            size = 0.6), axis.ticks.length = unit(0.13, "cm"), 
        legend.title = element_blank(), legend.justification = "center", 
        plot.title = element_text(color = "#000000", size = 20, 
            face = "bold", hjust = 0.5), )
    mytheme
}
FeatureScatterShell <-
function (object, feature1 = NULL, feature2 = NULL, outfile = NULL, 
    group.by = "orig.ident", cols = NULL, span = NULL) 
{
    p <- FeatureScatter(object = object, feature1 = feature1, 
        feature2 = feature2, span = span, group.by = group.by, 
        pt.size = 0.1, cols = cols) + guides(color = guide_legend(override.aes = list(size = 2), 
        title = NULL))
    p <- p + dot_theme_default()
    ggsave(p, file = outfile, width = 7, height = 6)
}
fetch_color <-
function (n = 0, type = c(names(color.list), "random"), tag = NULL, 
    is.extend = TRUE, verbose = FALSE) 
{
    type <- match.arg(type)
    if (type == "random") {
        color.use <- fetch_random_color(n = n, usepalette = T)
        if (verbose) {
            message("color : ", type, "->", tag, "->", n)
        }
        return(color.use)
    }
    tag <- match.arg(tag, names(color.list[[type]]))
    if (!is.numeric(n)) 
        stop("'n' must be numeric.")
    n.available <- as.numeric(names(color.list[[type]][[tag]]))
    n.select <- n.available[n.available >= n]
    if (length(n.select)) {
        if (n == min(n.select)) {
            n.select <- as.character(min(n.select))
            color.use <- color.list[[type]][[tag]][[n.select]]
        }
        else {
            n.select <- as.character(min(n.select))
            color.use <- head(color.list[[type]][[tag]][[n.select]], 
                n)
        }
    }
    else {
        n.select <- as.character(0)
        color.use <- color.list[[type]][[tag]][[n.select]]
        if (is.extend) 
            color.use <- colorRampPalette(color.list[[type]][[tag]][[n.select]])(n)
    }
    if (verbose) {
        message("color : ", type, "->", tag, "->", n.select)
    }
    return(color.use)
}
fetch_random_color <-
function (n = 1, usepalette = FALSE, hue = " ", luminosity = " ") 
{
    library(randomcoloR)
    if (usepalette == TRUE) {
        set.seed(1)
        color.use <- distinctColorPalette(k = n)
    }
    else {
        color.use <- randomColor(count = n, hue = hue, luminosity = luminosity)
    }
    return(color.use)
}
FilterCells <-
function (object, parameter, do.stat = TRUE, set.num = IfNull(parameter$filter$set.num, 
    "none"), set.num.seed = IfNull(parameter$filter$set.num.seed, 
    42), standard = parameter$filter$standard) 
{
    message("--->Filter Cells<---")
    cells.use <- Cells(object)
    if (!is.null(set.num) && set.num != "none") {
        if (set.num == "min") {
            cell_num <- min(table(object@meta.data$orig.ident))
        }
        else if (set.num != "none") {
            cell_num <- min(as.integer(set.num), max(table(object@meta.data$orig.ident)))
        }
        seed <- ifelse(!is.null(set.num.seed), set.num.seed, 
            42)
        set.seed(seed)
        cells.use <- as.character(unlist(by(Cells(object), object@meta.data$orig.ident, 
            function(x) sample(x, min(cell_num, length(x))))))
    }
    pm.used <- list()
    for (i in names(standard)) {
        if (exists(i, object@meta.data)) {
            if (length(standard[[i]]) == 1) {
                value <- standard[[i]]
                if (value == "auto") {
                  standard[[i]] <- autothres(data = object@meta.data[[i]], 
                    name = i, bin = 100)
                  print(standard[[i]])
                }
                else {
                  cells.use <- object@meta.data %>% tibble::rownames_to_column(var = "cells") %>% 
                    filter(.data[[i]] == value & cells %in% cells.use) %>% 
                    select(cells) %>% unlist()
                }
            }
            if (length(standard[[i]]) == 2) {
                lower <- standard[[i]][1]
                upper <- standard[[i]][2]
                cells.use <- object@meta.data %>% tibble::rownames_to_column(var = "cells") %>% 
                  filter(.data[[i]] >= lower & .data[[i]] <= 
                    upper & cells %in% cells.use) %>% select(cells) %>% 
                  unlist()
            }
            pm.used[[i]] <- standard[[i]]
        }
    }
    yaml::write_yaml(pm.used, file = "filtered_used_parameter.yaml")
    object <- object[, cells.use]
    if (do.stat) {
        StatFilterCells(object)
    }
    return(object)
}
FilterGenes <-
function (object, parameter = list(), min.cell = IfNull(parameter$filter$min.cell, 
    0)) 
{
    message("--->Filter Genes<---")
    if (min.cell > 0) {
        if (min.cell < 1) {
            min.cell <- min.cell * length(object@cell.names)
        }
        num.cells <- Matrix::rowSums(object@assays$RNA@counts > 
            0)
        genes.filter <- num.cells[which(num.cells < min.cell)]
        genes.filter <- FindFeaturesName(object, genes.filter)
        write.table(genes.filter, file = "filtered_genes.xls", 
            quote = F, sep = "\t", col.name = F)
        genes.use <- names(num.cells[which(num.cells >= min.cell)])
        object <- object[genes.use, ]
    }
    return(object)
}
FindFeaturesID <-
function (object, features, unlist = TRUE) 
{
    object@misc$fdata <- AddUnderscore(object@misc$fdata)
    if (all(rownames(object) %in% object@misc$fdata$dash)) {
        rownames(object@misc$fdata) <- object@misc$fdata$dash
    }
    features <- sapply(X = features, FUN = function(x) {
        if (!exists("fdata", object@misc)) 
            return(NULL)
        g1 <- toupper(rownames(object@misc$fdata)) %in% toupper(x)
        if (sum(g1) > 0) 
            return(rownames(object@misc$fdata)[g1])
        g2 <- toupper(object@misc$fdata$name) %in% toupper(x)
        if (sum(g2) > 0) 
            return(rownames(object@misc$fdata)[g2])
        g3 <- toupper(object@misc$fdata$merge_name) %in% toupper(x)
        if (sum(g3) > 0) 
            return(rownames(object@misc$fdata)[g3])
        message("[WARNING] '", x, "' not found gene id.")
        return(NULL)
    })
    if (unlist) 
        features <- unlist(features)
    return(features)
}
FindFeaturesName <-
function (object, features, col = "merge_name", is.fast = FALSE) 
{
    if (!exists("fdata", object@misc)) 
        return(features)
    object@misc$fdata <- AddUnderscore(object@misc$fdata)
    new <- gsub("_", "-", features)
    if (all(new %in% object@misc$fdata$dash)) {
        rownames(object@misc$fdata) <- object@misc$fdata$dash
        features <- new
    }
    if (is.fast) {
        Name <- object@misc$fdata[features, col]
        names(Name) <- features
    }
    else {
        Name <- sapply(features, function(x) {
            id <- object@misc$fdata[x, col]
            ifelse(is.null(id) || is.na(id), x, id)
        })
    }
    return(Name)
}
FindRegressVars <-
function (object, parameter = list(), vars.regress = NULL, force_recal = FALSE, 
    is_rm_cc = IfNull(parameter$cell_cycle$is_remove, TRUE), 
    is_rm_all_cc_signal = IfNull(parameter$cell_cycle$is_rm_all_signal, 
        FALSE), ...) 
{
    if (is.null(vars.regress)) {
        vars.regress <- paste0("nCount_", DefaultAssay(object))
        if (exists("percent.mito", object@meta.data)) {
            vars.regress <- c(vars.regress, "percent.mito")
        }
        if (is_rm_cc) {
            message("-->CellCycle Scoring<--")
            if (!exists("CC.Difference", object@meta.data) || 
                force_recal) {
                object <- DoCellCycleScoring(object, ...)
            }
            if (is_rm_all_cc_signal) {
                vars.regress <- c(vars.regress, "S.Score", "G2M.Score")
            }
            else {
                vars.regress <- c(vars.regress, "CC.Difference")
            }
        }
    }
    else if (vars.regress == "none") {
        vars.regress <- NULL
    }
    object@misc$vars.regress <- vars.regress
    return(object)
}
FindTopMarker <-
function (object.markers, top_num = 20, object = NULL, outfile = "Top.avg_exp.xls.tmp") 
{
    top <- object.markers %>% group_by(cluster) %>% arrange(desc(avg_logFC), 
        p_val, p_val_adj, .by_group = TRUE) %>% filter(1:n() <= 
        top_num)
    if (!is.null(outfile)) {
        tmp <- top %>% select(Cluster = cluster, Gene_ID = gene)
        if (!is.null(object)) {
            tmp$Gene_ID <- ChangeOUTName(tmp$Gene_ID, object@misc$fdata)
        }
        WriteTable(tmp, outfile)
    }
    return(top)
}
heatmap_theme_default <-
function () 
{
    library(ggplot2)
    mytheme <- theme_void() + theme(axis.text = element_text(color = "#000000", 
        size = 12), axis.text.x = element_text(angle = 45, hjust = 1, 
        vjust = 1), axis.text.y = element_text(hjust = 0, vjust = 0.5), 
        legend.title = element_blank(), legend.justification = c(0.5, 
            0.5), legend.margin = margin(5, 5, 5, 5), legend.text = element_text(hjust = 0.5, 
            vjust = 0.5, margin = margin(2, 0, 0, 0)), plot.title = element_text(size = 20, 
            face = "bold", hjust = 0.5), plot.margin = unit(c(5, 
            5, 5, 5), "mm"))
    if ("showtext" %in% installed.packages()) {
        library(showtext)
        showtext_auto(enable = TRUE)
        font_add("Arial", regular = "/Bio/Bin/pipeline/GeneralPlot/v1.0/fonts/tff/msttcore/arial.ttf")
        mytheme <- mytheme + theme(text = element_text(family = "Arial"))
    }
    mytheme
}
heatmap_theme_ggcor <-
function () 
{
    library(ggplot2)
    library(ggcor)
    mytheme <- theme_cor() + theme(axis.text = element_text(color = "#000000", 
        size = 12), legend.title = element_blank(), legend.text = element_text(hjust = 0.5, 
        vjust = 0.5, margin = margin(2, 0, 0, 0)), plot.title = element_text(size = 20, 
        face = "bold", hjust = 0.5), )
    if ("showtext" %in% installed.packages()) {
        library(showtext)
        showtext_auto(enable = TRUE)
        font_add("Arial", regular = "/Bio/Bin/pipeline/GeneralPlot/v1.0/fonts/tff/msttcore/arial.ttf")
        mytheme <- mytheme + theme(text = element_text(family = "Arial"))
    }
    mytheme
}
hist_theme_default <-
function () 
{
    library(ggplot2)
    mytheme <- theme_bw() + theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.border = element_rect(color = "#000000", 
            size = 1), axis.text = element_text(color = "#000000", 
            size = 12), axis.text.x = element_text(angle = 0, 
            hjust = 0.5, vjust = 0.5), axis.text.y = element_text(hjust = 1, 
            vjust = 0.5), axis.title = element_text(color = "#000000", 
            size = 16, face = "bold"), axis.title.x = element_text(margin = margin(2.5, 
            0, 2.5, 0, "mm")), axis.title.y = element_text(margin = margin(0, 
            2.5, 0, 2.5, "mm")), axis.ticks = element_line(color = "#000000", 
            size = 0.6), axis.ticks.length = unit(0.13, "cm"), 
        legend.title = element_blank(), strip.text = element_blank(), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        )
    if ("showtext" %in% installed.packages()) {
        library(showtext)
        showtext_auto(enable = TRUE)
        font_add("Arial", regular = "/Bio/Bin/pipeline/GeneralPlot/v1.0/fonts/tff/msttcore/arial.ttf")
        mytheme <- mytheme + theme(text = element_text(family = "Arial"))
    }
    mytheme
}
IfNull <-
function (var, default = NULL) 
{
    if (is.null(var)) {
        return(default)
    }
    else {
        return(var)
    }
}
KDE <-
function (x, y, z = NULL, h, n = 25, ...) 
{
    hx <- MASS::bandwidth.nrd(x)
    hy <- MASS::bandwidth.nrd(y)
    if (hx <= 0) 
        hx <- diff(range(x))/n
    if (hy <= 0) 
        hy <- diff(range(y))/n
    if (is.null(z)) {
        if (missing(h)) 
            h <- c(hx, hy)
        d <- MASS::kde2d(x = x, y = y, h = h, n = n, ...)
        gr <- data.frame(with(d, expand.grid(x, y)), as.vector(d$z))
        colnames(gr) <- c("xgr", "ygr", "zgr")
        mod <- loess(zgr ~ xgr * ygr, data = gr)
        dens <- predict(mod, newdata = data.frame(xgr = x, ygr = y))
    }
    else {
        hz <- MASS::bandwidth.nrd(z)
        if (hz <= 0) 
            hz <- diff(range(z))/n
        if (missing(h)) 
            h <- c(hx, hy, hz)
        d <- kde3d(x, y, z, h, n = n, ...)
        gr <- data.frame(with(d, expand.grid(x, y, z)), as.vector(d$d))
        colnames(gr) <- c("xgr", "ygr", "zgr", "dgr")
        mod <- loess(dgr ~ xgr * ygr * zgr, data = gr)
        dens <- predict(mod, newdata = data.frame(xgr = x, ygr = y, 
            zgr = z))
    }
    return(dens)
}
kde3d <-
function (x, y, z, h, n = 20, lims = c(range(x), range(y), range(z))) 
{
    nx <- length(x)
    if (length(y) != nx || length(z) != nx) 
        stop("data vectors must be the same length")
    if (missing(h)) 
        h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y), 
            MASS::bandwidth.nrd(z))/6
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
    tmy.nx <- t(my)/nx
    for (k in 1:n[3]) {
        tmy.nz.zk <- tmy.nx * mz[k, ]
        v[, , k] <- mx %*% tmy.nz.zk
    }
    return(list(x = gx, y = gy, z = gz, d = v))
}
line_theme_default <-
function () 
{
    library(ggplot2)
    options(scipen = -1)
    mytheme <- theme_bw() + theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.border = element_rect(color = "#000000", 
            size = 1), axis.text = element_text(color = "#000000", 
            size = 12), axis.text.x = element_text(angle = 0, 
            hjust = 0.5, vjust = 0.5), axis.text.y = element_text(hjust = 1, 
            vjust = 0.5), axis.title = element_text(color = "#000000", 
            size = 16, face = "bold"), axis.title.x = element_text(margin = margin(2.5, 
            0, 2.5, 0, "mm")), axis.title.y = element_text(margin = margin(0, 
            2.5, 0, 2.5, "mm")), axis.ticks = element_line(color = "#000000", 
            size = 0.6), axis.ticks.length = unit(0.13, "cm"), 
        legend.title = element_blank(), plot.title = element_text(size = 20, 
            face = "bold", hjust = 0.5), )
    if ("showtext" %in% installed.packages()) {
        library(showtext)
        showtext_auto(enable = TRUE)
        font_add("Arial", regular = "/Bio/Bin/pipeline/GeneralPlot/v1.0/fonts/tff/msttcore/arial.ttf")
        mytheme <- mytheme + theme(text = element_text(family = "Arial"))
    }
    mytheme
}
ListCellCluster <-
function (object, outfile = "Cells.cluster.list.xls") 
{
    data <- .GetMetaData(object, cols = c(Cluster = "seurat_clusters", 
        Samples = "orig.ident", "Groups"))
    WriteTable(data, file = outfile)
}
ListMarker <-
function (object, object.markers, outfile = "DeGene.list.xls", 
    is.return = FALSE, is.fast = FALSE, group.by = "seurat_clusters", 
    assay = DefaultAssay(object), slot = "data", group.by.name = "Cluster") 
{
    Targets_name <- paste0("Target_", group.by.name)
    Others_name <- paste0("Other_", group.by.name)
    Targets_mean_name <- paste0(Targets_name, "_mean")
    Others_mean_name <- paste0(Others_name, "_mean")
    Targets_mean <- CalAvgExp(object, unique(object.markers$gene), 
        is.return = T, is.reverse = F, assay = assay, slot = slot, 
        group.by = group.by) %>% reshape2::melt(varnames = c("gene", 
        "cluster"), value.name = Targets_mean_name) %>% mutate(cluster = factor(cluster))
    Others_mean <- CalAvgExp(object, unique(object.markers$gene), 
        is.return = T, is.reverse = T, assay = assay, slot = slot, 
        group.by = group.by) %>% reshape2::melt(varnames = c("gene", 
        "cluster"), value.name = Others_mean_name) %>% mutate(cluster = factor(cluster))
    Name <- FindFeaturesName(object, unique(object.markers$gene), 
        "name", is.fast = is.fast)
    marker_list <- object.markers %>% left_join(y = Targets_mean) %>% 
        left_join(y = Others_mean) %>% mutate(Log2FC = log2(!!as.name(Targets_mean_name)/!!as.name(Others_mean_name)), 
        name = Name[gene]) %>% select(`:=`(!!Targets_name, cluster), 
        `Gene ID` = gene, `Gene Name` = name, !!Targets_mean_name, 
        !!Others_mean_name, Log2FC, Pvalue = p_val, Qvalue = p_val_adj)
    if (is.null(outfile) || is.return) {
        return(marker_list)
    }
    else {
        marker_list[["Gene ID"]] <- ChangeOUTName(marker_list[["Gene ID"]], 
            object@misc$fdata)
        WriteTable(marker_list, outfile)
    }
}
Load <-
function (file) 
{
    object <- readRDX(file)
    if ("version" %in% slotNames(object)) {
        if (grepl("^2", object@version)) {
            object <- Seurat::UpdateSeuratObject(object)
        }
    }
    return(object)
}
MakeSeuratObj <-
function (parameter = list(), assay = "RNA", data_name = parameter$data$name, 
    data_dir = parameter$data$dir, name_list = parameter$name_list, 
    group.use = parameter$Groups) 
{
    object.list <- list()
    for (i in seq(data_name)) {
        mat <- Read10X(data.dir = data_dir[i], gene.column = 1)
        object.list[[i]] <- CreateSeuratObject(counts = mat, 
            project = data_name[i], assay = assay)
    }
    if (length(object.list) == 1) {
        object <- object.list[[1]]
    }
    else {
        object <- merge(x = object.list[[1]], y = unlist(object.list[-1]), 
            add.cell.ids = data_name)
    }
    object@meta.data$orig.ident <- factor(object@meta.data$orig.ident, 
        levels = data_name)
    object <- SetIdent(object, value = "orig.ident")
    object@misc[["fdata"]] <- AddFData(object, name_list)
    object@misc[["pdata"]] <- FetchData(object, c("orig.ident", 
        paste0("nFeature_", assay), paste0("nCount_", assay)))
    object@misc[["counts"]] <- GetAssayData(object, slot = "counts", 
        assay = assay)
    color.sample <- fetch_color(nlevels(object@meta.data$orig.ident), 
        "tsne", "set3")
    names(color.sample) <- levels(object@meta.data$orig.ident)
    object@misc[["color.sample"]] <- color.sample
    if (!is.null(group.use)) {
        object@meta.data$Groups <- object@meta.data$orig.ident
        group <- unlist(lapply(names(group.use), function(i) {
            x <- rep(i, length.out = length(group.use[[i]]))
            names(x) <- group.use[[i]]
            x
        }))
        levels(object@meta.data$Groups) <- group[levels(object@meta.data$Groups)]
        object@meta.data$Groups <- factor(object@meta.data$Groups, 
            levels = names(group.use))
        object@misc[["pdata"]]$Groups <- object@meta.data$Groups
        color.group <- fetch_color(nlevels(object@meta.data$orig.ident), 
            "tsne", "set2")
        names(color.group) <- levels(object@meta.data$Groups)
        object@misc[["color.group"]] <- color.group
    }
    return(object)
}
Otsu <-
function (data, bin = 256, is.return.range = TRUE) 
{
    dt0 <- bin_stat(data, bin = bin)
    dt <- table(dt0)
    dt[dt > 100] <- 100
    thres <- core_otsu(dt)
    if (is.return.range) 
        return(range(as.numeric(names(dt0)[dt0 %in% names(dt)[dt >= 
            thres]])))
    else return(as.numeric(names(dt)[dt == thres - 1]))
}
Otsu.bak <-
function (data, digits = 0, is.return.range = TRUE) 
{
    dt0 <- table(round(data, digits = digits))
    dt <- dt0
    dt[dt > 100] <- 100
    thres <- core_otsu(dt)
    if (is.return.range) 
        return(range(as.numeric(names(dt0)[dt0 >= thres])))
    else return(thres)
}
pie_theme_default <-
function () 
{
    library(ggplot2)
    mytheme <- theme_void() + theme(legend.title = element_blank(), 
        legend.key.height = unit(5, "mm"), legend.key.width = unit(5, 
            "mm"), strip.text = element_text(size = 16, face = "bold", 
            hjust = 0.5), plot.title = element_text(size = 20, 
            face = "bold", hjust = 0.5), )
    if ("showtext" %in% installed.packages()) {
        library(showtext)
        showtext_auto(enable = TRUE)
        font_add("Arial", regular = "/Bio/Bin/pipeline/GeneralPlot/v1.0/fonts/tff/msttcore/arial.ttf")
        mytheme <- mytheme + theme(text = element_text(family = "Arial"))
    }
    mytheme
}
PlotAboutFeatures <-
function (object, features = NULL, group.by = "seurat_clusters", 
    outpref = NULL, group.colors = NULL, plot.feature = FALSE) 
{
    object@meta.data <- droplevels(object@meta.data)
    PlotDotPlot(object, features = features, group.by = group.by, 
        outfile = paste0(outpref, ".DotPlot.pdf"))
    PlotHeatmapPlot(object, features = features, group.by = group.by, 
        outfile = paste0(outpref, ".Heatmap.pdf"), group.colors = group.colors)
    if (plot.feature) {
        PlotFeaturePlot(object, features = features, outfile = paste0(outpref, 
            ".Distribution.pdf"), reduction = "tsne")
    }
}
PlotBasicStat <-
function (object, outpref = NULL, color = NULL, assay = NULL, 
    group.by = "orig.ident", span = NULL, ...) 
{
    if (is.null(color)) {
        color <- switch(group.by, Groups = object@misc[["color.group"]], 
            orig.ident = object@misc[["color.sample"]], seurat_clusters = object@misc[["color.cluster"]])
    }
    if (is.null(assay)) {
        assay <- DefaultAssay(object)
    }
    nFeature <- paste0("nFeature_", assay)
    nCount <- paste0("nCount_", assay)
    d <- data.frame(nGene = c(nFeature, "Number"), nUMI = c(nCount, 
        "Number"), pMito = c("percent.mito", "Percentage(%)"), 
        stringsAsFactors = FALSE)
    for (i in names(d)) {
        if (!exists(d[[i]][1], object@meta.data)) 
            d[[i]] <- NULL
    }
    VlnplotShell(object, features = unlist(d[1, ]), labs.y = unlist(d[2, 
        ]), titles = colnames(d), outfile = paste0(outpref, ".merge.pdf"), 
        cols.use = color, group.by = group.by, assay = assay, 
        ...)
    d <- data.frame(expected.marker = c("expected.marker", "Number"), 
        exclude.marker = c("exclude.marker", "Number"), stringsAsFactors = FALSE)
    for (i in names(d)) {
        if (!exists(d[[i]][1], object@meta.data)) 
            d[[i]] <- NULL
    }
    VlnplotShell(object, features = unlist(d[1, ]), labs.y = unlist(d[2, 
        ]), titles = colnames(d), outfile = paste0(outpref, ".PresetMarker.pdf"), 
        cols.use = color, group.by = group.by, assay = assay, 
        ...)
    FeatureScatterShell(object, feature1 = nCount, feature2 = nFeature, 
        span = span, outfile = paste0(outpref, ".nUMI-nGene.pdf"), 
        cols = color, group.by = group.by)
    if (exists("percent.mito", object@meta.data)) {
        FeatureScatterShell(object, feature1 = nCount, feature2 = "percent.mito", 
            span = span, outfile = paste0(outpref, ".nUMI-pMito.pdf"), 
            cols = color, group.by = group.by)
    }
}
PlotCluster <-
function (object, reduction = "umap", p1.group.by = "orig.ident", 
    split.by = p1.group.by, outpref = NULL, ...) 
{
    .PlotCluster(object, reduction = reduction, p1.group.by = p1.group.by, 
        outfile = paste0(outpref, ".pdf"), ...)
    for (i in unique(object@meta.data[[split.by]])) {
        cells.use <- rownames(object@meta.data)[object@meta.data[[split.by]] == 
            i]
        .PlotCluster(object, reduction = reduction, cells = cells.use, 
            p1.group.by = p1.group.by, outfile = paste0(outpref, 
                ".", i, ".pdf"), ...)
    }
    data <- object[[reduction]]@cell.embeddings %>% as.data.frame() %>% 
        tibble::rownames_to_column(var = "Cells") %>% left_join(.GetMetaData(object, 
        cols = c(Samples = "orig.ident", Cluster = "seurat_clusters", 
            "Groups")))
    WriteTable(data, file = paste0(outpref, ".plot.data.tmp"))
}
PlotDensityPlot <-
function (object, features = NULL, reduction = "umap", outpref = "DensityPlot/", 
    ...) 
{
    if (is.null(features)) {
        .PlotDensityPlot(object, reduction = reduction, outpref = outpref, 
            ...)
    }
    else {
        for (i in features) {
            .PlotDensityPlot(object, i, reduction = reduction, 
                outpref = outpref, ...)
        }
    }
}
PlotDotPlot <-
function (object, features = NULL, outfile = NULL, group.by = "seurat_clusters", 
    is.use.name = TRUE, ...) 
{
    p <- DotPlot(object, features = features, group.by = group.by, 
        ...) + RotatedAxis()
    p <- p + dot_theme_default() + theme(axis.text.x = element_text(angle = 45, 
        hjust = 1, vjust = 1))
    if (is.use.name) {
        levels(p$data$features.plot) <- FindFeaturesName(object, 
            levels(p$data$features.plot))
    }
    if (is.null(outfile)) {
        return(p)
    }
    else {
        w <- max(6, ceiling(length(features)) * 0.35 + 2)
        h <- max(6, length(unique(object@meta.data[[group.by]])) * 
            0.4)
        ggsave(p, file = outfile, width = w, height = h, limitsize = FALSE)
    }
}
PlotFeaturePlot <-
function (object, features, reduction = "umap", is.combine = TRUE, 
    outpref = NULL, outfile = NULL, ...) 
{
    if (is.combine) {
        .PlotFeaturePlot(object, features = features, reduction = reduction, 
            outfile = outfile, ...)
    }
    else {
        if (is.null(outpref)) 
            outpref <- "ExpPlot"
        for (i in features) {
            name <- FindFeaturesName(object, i)
            name <- gsub("[ /\\]", "_", name)
            .PlotFeaturePlot(object, features = i, reduction = reduction, 
                outfile = paste(c(outpref, name, "pdf"), collapse = "."), 
                plot.basic.size = 6, ...)
        }
    }
}
PlotHeatmapPlot <-
function (object, features = NULL, group.by = "seurat_clusters", 
    is.use.name = TRUE, outfile = NULL, group.colors = NULL) 
{
    if (is.null(group.colors)) {
        group.colors <- if (group.by == "seurat_clusters") 
            object@misc$color.cluster
        else if (group.by == "orig.ident") 
            object@misc$color.sample
        else NULL
    }
    else {
        group.colors <- group.colors[levels(object@meta.data[[group.by]])]
    }
    p <- DoHeatmap(object = object, features = features, cells = NULL, 
        group.by = group.by, group.colors = group.colors, combine = FALSE, 
        raster = FALSE)
    p <- p[[1]]
    p <- p + theme(legend.title = element_blank())
    p$layers[[2]] <- NULL
    if (is.use.name) {
        levels(p$data$Feature) <- FindFeaturesName(object, levels(p$data$Feature))
    }
    if (is.null(outfile)) {
        return(p)
    }
    else {
        h <- max(7, length(unique(features)) * 0.11 + 2.5)
        w <- h * 4/3
        ggsave(p, file = outfile, width = w, height = h, limitsize = FALSE)
    }
}
PlotPresetMarker <-
function (object, cols.use = object@misc$color.cluster, group.by = "seurat_clusters", 
    outpref = "PresetMarker") 
{
    PresetMarker <- union(object@misc[["expected.marker"]], object@misc[["more.marker"]])
    PresetMarker <- union(PresetMarker, object@misc[["exclude.marker"]])
    if (length(PresetMarker)) {
        VlnplotShell(object, features = PresetMarker, outfile = paste0(outpref, 
            ".VlnPlot.pdf"), titles = object@misc$fdata[PresetMarker, 
            "merge_name"], cols.use = cols.use, group.by = group.by)
        PlotAboutFeatures(object, features = PresetMarker, outpref = outpref, 
            plot.feature = TRUE, group.by = group.by)
    }
}
PlotVlnPlot <-
function (object, features, outpref = NULL, group.by = "seurat_clusters", 
    cols.use = NULL, ...) 
{
    name <- FindFeaturesName(object, features)
    names(features) <- as.vector(name)
    if (is.null(cols.use)) {
        cols.use <- switch(group.by, seurat_clusters = object@misc$color.cluster, 
            orig.ident = object@misc$color.sample)
    }
    VlnplotShell(object, features, do.split.save = TRUE, pref = outpref, 
        group.by = group.by, titles = name, cols.use = cols.use, 
        ...)
}
readRDX <-
function (file) 
{
    con <- gzfile(file)
    on.exit(close(con))
    magic <- readChar(con, 5L, useBytes = TRUE)
    if (grepl("RD[ABX][2-9]\n", magic)) {
        object <- get(load(file))
    }
    else {
        object <- readRDS(file)
    }
    return(object)
}
RenameFeatures <-
function (object, new.names = NULL, form.type = c("id", "merge_name"), 
    to.type = c("id", "merge_name")) 
{
    if (is.null(new.names)) {
        form.type <- match.arg(form.type)
        to.type <- match.arg(to.type)
        fdata <- object@misc$fdata
        fdata$id <- rownames(fdata)
        new.names <- fdata[, to.type]
        names(new.names) <- fdata[, form.type]
    }
    assays <- Seurat:::FilterObjects(object = object, classes.keep = "Assay")
    for (assay in assays) {
        slot(object = object, name = "assays")[[assay]] <- RenameFeatures.Assays(object = object[[assay]], 
            new.names = new.names)
    }
    dimreducs <- Seurat:::FilterObjects(object = object, classes.keep = "DimReduc")
    for (dr in dimreducs) {
        object[[dr]] <- RenameFeatures.DimReduc(object = object[[dr]], 
            new.names = new.names)
    }
    return(object)
}
RenameFeatures.Assays <-
function (object, new.names = NULL) 
{
    for (data.slot in c("counts", "data", "scale.data")) {
        old.data <- GetAssayData(object = object, slot = data.slot)
        if (nrow(x = old.data) <= 1) {
            next
        }
        old.name <- rownames(x = slot(object = object, name = data.slot))
        rownames(x = slot(object = object, name = data.slot)) <- new.names[old.name]
    }
    if (length(slot(object = object, name = "var.features")) > 
        0) {
        old.name <- rownames(x = slot(object = object, name = "var.features"))
        slot(object = object, name = "var.features") <- new.names[old.name]
    }
    return(object)
}
RenameFeatures.DimReduc <-
function (object, new.names = NULL) 
{
    for (projected in c(TRUE, FALSE)) {
        data.slot <- ifelse(projected, "feature.loadings.projected", 
            "feature.loadings")
        old.data <- Loadings(object = object, projected = projected)
        rownames(x = old.data) <- new.names[rownames(x = old.data)]
        slot(object = object, name = data.slot) <- old.data
    }
    return(object)
}
RestoreObject <-
function (object) 
{
    new.object <- CreateSeuratObject(counts = object@misc[["counts"]], 
        meta.data = object@misc[["pdata"]], project = object@project.name, 
        assay = "RNA")
    new.object@misc <- object@misc
    return(new.object)
}
show_all_color <-
function () 
{
    library(dplyr, warn.conflicts = F)
    dt <- reshape2::melt(color.list)
    dt <- dt %>% group_by(L1, L2, L3) %>% mutate(x = 1:n()) %>% 
        arrange(L1, L2, as.numeric(L3), x)
    dt$L3 <- factor(dt$L3, levels = as.character(0:max(as.numeric(dt$L3))))
    cc <- unique(as.character(dt$value))
    names(cc) <- cc
    library(ggplot2, warn.conflicts = F)
    p <- list()
    for (i in unique(dt$L1)) {
        p[[i]] <- ggplot(dt %>% filter(L1 == i), aes(x = x, y = L3)) + 
            geom_tile(aes(fill = value), color = "grey") + facet_grid(L2 ~ 
            ., scales = "free_y", switch = "y") + scale_fill_manual(values = cc) + 
            scale_x_continuous(expand = expand_scale()) + ylab(i) + 
            theme_minimal() + theme(legend.position = "none", 
            strip.placement = "outside", panel.grid = element_blank(), 
            axis.title.x = element_blank(), axis.text.x = element_blank())
    }
    cowplot::plot_grid(plotlist = p, ncol = 1, axis = "ltrb", 
        align = "hv")
}
SplitObject.Image <-
function (objects, names.in = "orig.ident") 
{
    if (class(objects) != "list") {
        if ("images" %in% slotNames(objects) && length(objects@images) > 
            1) {
            keep.image <- levels(droplevels(objects@meta.data)[[names.in]])
            objects@images <- objects@images[keep.image]
        }
    }
    else {
        for (i in seq(objects)) {
            if ("images" %in% slotNames(objects[[i]]) && length(objects[[i]]@images) > 
                1) {
                keep.image <- levels(droplevels(objects[[i]]@meta.data)[[names.in]])
                objects[[i]]@images <- objects[[i]]@images[keep.image]
            }
        }
    }
    return(objects)
}
StatCluster <-
function (object, group.by = "orig.ident", outpref = "Cluster.stat", 
    ...) 
{
    .StatCluster(object, outpref = outpref)
    .StatCluster_by(object, group.by = group.by, outpref = outpref)
    .PlotClusterStat(object, group.by = group.by, outpref = outpref, 
        ...)
}
StatFeatures <-
function (object, features = NULL, col.name = NULL, stat_pct = FALSE, 
    assay = NULL, add_to_pdata = FALSE) 
{
    if (file.exists(features[1])) 
        features <- readLines(con = features[1])
    features <- FindFeaturesID(object = object, features = features, 
        unlist = FALSE)
    if (is.null(assay)) 
        assay <- DefaultAssay(object = object)
    metadata <- Matrix::colSums(x = GetAssayData(object = object, 
        slot = "counts", assay = assay)[unlist(features), , drop = FALSE])
    if (stat_pct) 
        metadata <- metadata/object@meta.data[[paste0("nCount_", 
            assay)]] * 100
    if (!is.null(x = col.name)) {
        object@misc[[col.name]] <- features
        if (add_to_pdata) 
            object@misc$pdata[[col.name]] <- metadata
        object <- AddMetaData(object = object, metadata = metadata, 
            col.name = col.name)
        return(object)
    }
    else {
        return(metadata)
    }
}
StatFilterCells <-
function (object, group.by = "orig.ident", outfile = "Filter.stat.xls") 
{
    filter_cells <- setdiff(rownames(object@misc$pdata), colnames(object))
    write.table(filter_cells, file = "filtered_cells.xls", quote = F, 
        sep = "\t", col.name = F)
    name <- if (group.by == "orig.ident") 
        "Samples"
    else group.by
    group.by <- as.name(group.by)
    if (exists("percent.mito", object@meta.data) && exists("percent.mito", 
        object@misc$pdata)) {
        a <- object@meta.data %>% group_by(`:=`(!!name, !!group.by)) %>% 
            summarise(after_filter_num = n(), after_filter_median_UMI_per_cell = median(nCount_RNA), 
                after_filter_median_genes_per_cell = median(nFeature_RNA), 
                after_filter_median_MT_per_cell = median(percent.mito))
        b <- object@misc$pdata %>% group_by(`:=`(!!name, !!group.by)) %>% 
            summarise(before_filter_num = n(), before_filter_median_UMI_per_cell = median(nCount_RNA), 
                before_filter_median_genes_per_cell = median(nFeature_RNA), 
                before_filter_median_MT_per_cell = median(percent.mito))
        filter_stat_table <- left_join(a, b) %>% mutate(pct = paste0(round(after_filter_num/before_filter_num * 
            100, 2), "%")) %>% select(!!name, before_filter_num, 
            after_filter_num, pct, before_filter_median_UMI_per_cell, 
            after_filter_median_UMI_per_cell, before_filter_median_genes_per_cell, 
            after_filter_median_genes_per_cell, before_filter_median_MT_per_cell, 
            after_filter_median_MT_per_cell)
    }
    else {
        a <- object@meta.data %>% group_by(`:=`(!!name, !!group.by)) %>% 
            summarise(after_filter_num = n(), after_filter_median_UMI_per_cell = median(nCount_RNA), 
                after_filter_median_genes_per_cell = median(nFeature_RNA))
        b <- object@misc$pdata %>% group_by(`:=`(!!name, !!group.by)) %>% 
            summarise(before_filter_num = n(), before_filter_median_UMI_per_cell = median(nCount_RNA), 
                before_filter_median_genes_per_cell = median(nFeature_RNA))
        filter_stat_table <- left_join(a, b) %>% mutate(pct = paste0(round(after_filter_num/before_filter_num * 
            100, 2), "%")) %>% select(!!name, before_filter_num, 
            after_filter_num, pct, before_filter_median_UMI_per_cell, 
            after_filter_median_UMI_per_cell, before_filter_median_genes_per_cell, 
            after_filter_median_genes_per_cell)
    }
    write.table(filter_stat_table, file = outfile, quote = F, 
        sep = "\t", row.names = F)
}
StatMarker <-
function (object.markers, Cluster_name = "Cluster", Item_name = "Number of DE genes", 
    color = NULL, outpref = "DeGene.stat") 
{
    stat <- cbind(tibble(`:=`(!!Cluster_name, Item_name)) %>% 
        as.matrix(), t(table(object.markers$cluster)))
    WriteTable(stat, paste0(outpref, ".xls"))
    p <- ggplot(object.markers) + geom_bar(aes(x = cluster, fill = cluster), 
        stat = "count") + theme_light() + labs(x = Cluster_name, 
        y = Item_name)
    p <- p + bar_theme_default()
    p <- p + theme(legend.position = "none")
    if (!is.null(color)) 
        p <- p + scale_fill_manual(values = color)
    ggsave(p, file = paste0(outpref, ".pdf"), height = 6, width = 8)
    invisible(stat)
}
VlnplotShell <-
function (object, features = NULL, labs.y = NULL, titles = NULL, 
    outfile = NULL, group.by = "orig.ident", cols.use = NULL, 
    do.split.save = FALSE, pref = NULL, ...) 
{
    if (is.null(features)) 
        return(1)
    nCol <- ifelse(length(features) == 3, 3, ceiling(sqrt(length(features))))
    nRow <- ceiling(length(features)/nCol)
    plots <- VlnPlot(object = object, features = features, ncol = nCol, 
        group.by = group.by, cols = cols.use, pt.size = 0.1, 
        combine = FALSE, ...)
    for (i in seq(features)) {
        if (group.by == "orig.ident") 
            plots[[i]] <- plots[[i]] + xlab("Samples")
        else if (group.by == "Groups") 
            plots[[i]] <- plots[[i]] + xlab("Groups")
        else if (group.by == "seurat_clusters") 
            plots[[i]] <- plots[[i]] + xlab("Clusters")
        if (!is.null(labs.y[i]) && !is.na(labs.y[i])) 
            plots[[i]] <- plots[[i]] + ylab(labs.y[i])
        if (!is.null(titles[i]) && !is.na(titles[i])) 
            plots[[i]] <- plots[[i]] + ggtitle(titles[i])
        plots[[i]]$layers[[1]]$aes_params$alpha <- 0.5
        plots[[i]]$layers <- rev(plots[[i]]$layers)
        plots[[i]] <- plots[[i]] + box_theme_default() + theme(axis.text.x = element_text(angle = 45, 
            hjust = 1, vjust = 1))
    }
    if (do.split.save && !is.null(pref)) {
        for (i in seq(features)) {
            name <- ifelse(!is.null(names(features[i])), names(features[i]), 
                features[i])
            outfile.tmp <- paste(pref, name, "pdf", sep = ".")
            ggsave(plots[[i]] + theme(legend.position = "none"), 
                file = outfile.tmp, width = 6, height = 6)
        }
    }
    plots.combined <- patchwork::wrap_plots(plots, ncol = nCol) & 
        theme(legend.position = "none")
    if (!is.null(dev.list())) 
        dev.off()
    if (is.null(outfile)) {
        return(plots.combined)
    }
    else {
        ggsave(plots.combined, file = outfile, width = 6 * nCol, 
            height = 6 * nRow, limitsize = F)
    }
}
WriteTable <-
function (x, file) 
{
    write.table(x, file = file, quote = F, sep = "\t", row.names = FALSE, 
        col.names = TRUE)
}
