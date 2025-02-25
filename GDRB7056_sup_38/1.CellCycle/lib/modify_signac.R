#########################
GeneActivity <- function (object, assay = NULL, features = NULL, extend.upstream = 2000, 
    extend.downstream = 0, biotypes = "protein_coding", max.width = 5e+05, 
    verbose = TRUE, use.genename = FALSE, ...) 
{
    if (!is.null(x = features)) {
        if (length(x = features) == 0) {
            stop("Empty list of features provided")
        }
    }
    assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
    if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
        stop("The requested assay is not a ChromatinAssay.")
    }
    annotation <- Annotation(object = object[[assay]])
    if (length(x = annotation) == 0) {
        stop("No gene annotations present in object")
    }
    if (verbose) {
        message("Extracting gene coordinates")
    }
    transcripts <- CollapseToLongestTranscript(ranges = annotation)
    if (!is.null(x = biotypes)) {
        transcripts <- transcripts[transcripts$gene_biotype %in% 
            biotypes]
        if (length(x = transcripts) == 0) {
            stop("No genes remaining after filtering for requested biotypes")
        }
    }
    if (!is.null(x = features)) {
		if ( use.genename ) { ## here
				transcripts <- transcripts[transcripts$gene_name %in% features]
		} else {
				transcripts <- transcripts[transcripts$gene_id %in% features]
		}
        if (length(x = transcripts) == 0) {
            stop("None of the requested genes were found in the gene annotation")
        }
    }
    if (!is.null(x = max.width)) {
        transcript.keep <- which(x = width(x = transcripts) < 
            max.width)
        transcripts <- transcripts[transcript.keep]
        if (length(x = transcripts) == 0) {
            stop("No genes remaining after filtering for max.width")
        }
    }
    transcripts <- Extend(x = transcripts, upstream = extend.upstream, 
        downstream = extend.downstream)
    frags <- Fragments(object = object[[assay]])
    if (length(x = frags) == 0) {
        stop("No fragment information found for requested assay")
    }
    cells <- colnames(x = object[[assay]])
    counts <- FeatureMatrix(fragments = frags, features = transcripts, 
        cells = cells, verbose = verbose, ...)
	gene.key <- if ( use.genename ) transcripts$gene_name else transcripts$gene_id ## here
    names(x = gene.key) <- GRangesToString(grange = transcripts)
    rownames(x = counts) <- as.vector(x = gene.key[rownames(x = counts)])
    counts <- counts[rownames(x = counts) != "", ]
    return(counts)
}
environment(GeneActivity) <- asNamespace("Signac")
assignInNamespace("GeneActivity", GeneActivity, ns = "Signac")

FindMotifs <- function (object, features, background = 40000, assay = NULL, 
    verbose = TRUE, ...) 
{
    assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
    background <- SetIfNull(x = background, y = rownames(x = object))
    if (is(object = background, class2 = "numeric")) {
        if (verbose) {
            message("Selecting background regions to match input ", 
                "sequence characteristics")
        }
        meta.feature <- GetAssayData(object = object, assay = assay, 
            slot = "meta.features")
        mf.choose <- meta.feature[setdiff(x = rownames(x = meta.feature), 
            y = features), , drop = FALSE]
        mf.query <- meta.feature[features, , drop = FALSE]
        background <- MatchRegionStats(meta.feature = mf.choose, 
            query.feature = mf.query, regions = features, n = background, 
            verbose = verbose, ...)
    }
    if (verbose) {
        msg <- ifelse(test = length(x = features) > 1, yes = " regions", 
            no = " region")
        message("Testing motif enrichment in ", length(x = features), 
            msg)
    }
    if (length(x = features) < 10) {
        warning("Testing motif enrichment using a small number of regions is ", 
            "not recommended")
    }
    motif.all <- GetMotifData(object = object, assay = assay, 
        slot = "data")
    motif.names <- GetMotifData(object = object, assay = assay, 
        slot = "motif.names")
    query.motifs <- motif.all[features, , drop = FALSE]
    background.motifs <- motif.all[background, , drop = FALSE]
    query.counts <- Matrix::colSums(x = query.motifs) # here
    background.counts <- Matrix::colSums(x = background.motifs) # here
    percent.observed <- query.counts/length(x = features) * 100
    percent.background <- background.counts/length(x = background) * 
        100
    fold.enrichment <- percent.observed/percent.background
    p.list <- vector(mode = "numeric")
    for (i in seq_along(along.with = query.counts)) {
        p.list[[i]] <- phyper(q = query.counts[[i]] - 1, m = background.counts[[i]], 
            n = nrow(x = background.motifs) - background.counts[[i]], 
            k = length(x = features), lower.tail = FALSE)
    }
    results <- data.frame(motif = names(x = query.counts), observed = query.counts, 
        background = background.counts, percent.observed = percent.observed, 
        percent.background = percent.background, fold.enrichment = fold.enrichment, 
        pvalue = p.list, motif.name = as.vector(x = unlist(x = motif.names[names(x = query.counts)])), 
        stringsAsFactors = FALSE)
    if (nrow(x = results) == 0) {
        return(results)
    }
    else {
        return(results[order(results[, 7], -results[, 6]), ])
    }
}
environment(FindMotifs) <- asNamespace("Signac")
assignInNamespace("FindMotifs", FindMotifs, ns = "Signac")

