set_sep <- if ( packageVersion("Signac") > 1 ) c( "-", "-") else c( ":", "-" )
split_sep <- if ( packageVersion("Signac") > 1 ) c( "-(?=\\d+-\\d+$)", "-(?=\\d+$)" ) else c( ":", "-(?=\\d+$)" )



AddMergePeaks <- function(object.list, assay.name = "peaks"){
		if ( length(object.list) == 1 )				
				return(object.list)

		message("--> merge peaks <--")
		all.peaks   <- unlist(lapply(object.list, function(x) rownames(x[[assay.name]])))
		all.ranges  <- StringToGRanges(all.peaks, sep = split_sep)
		merge.peaks <- GenomicRanges::reduce(all.ranges)
		if ( identical(all.ranges, merge.peaks) )
				return(object.list)

		for ( i in seq(object.list) ) {
				DefaultAssay(object.list[[i]]) <- assay.name
				new.matrix <- FeatureMatrix(fragments = Fragments(object.list[[i]]), features = merge.peaks, cells = colnames(object.list[[i]]))
				object.list[[i]][[assay.name]] <- CreateChromatinAssay(counts = new.matrix, sep = set_sep, fragments = Fragments(object.list[[i]]), min.cells = -1, min.features = -1)
				DefaultAssay(object.list[[i]]) <- assay.name
		}
		return(object.list)
}


GetDims <- function (object, assay = NULL, reduction = "lsi", n = NULL, ...) {
		if ( is.null(assay) ) assay <- DefaultAssay(object = object)
		dr <- object[[reduction]]
		embed <- Embeddings(object = dr)
		counts <- object[[paste0("nCount_", assay)]]
		embed <- embed[rownames(x = counts), ]
		if ( is.null(n) ) n <- ncol(x = embed)
		embed <- embed[, seq_len(length.out = n)]
		depth.cor <- as.data.frame(cor(x = embed, y = counts, ...))
		depth.cor$counts <- depth.cor[, 1]
		depth.cor$Component <- seq_len(length.out = nrow(x = depth.cor))
		depth.cor <- subset(depth.cor, abs(counts) <= 0.5)

		return(depth.cor$Component)
}

stat.ClusterCells <- function(object, outdir,
				assay = DefaultAssay(object),
				name.ident  = "seurat_clusters",
				name.nft    = paste0("nFeature_", assay),
				name.ncnt   = paste0("nCount_", assay)
				) {
		name.ident <- as.name(name.ident)
		name.nft   <- as.name(name.nft)
		name.ncnt  <- as.name(name.ncnt)


		st <- object@meta.data %>%
				group_by(!! name.ident) %>%
				summarise(cell_num = n(), cell_pct = cell_num, median_peak = median(!! name.nft), median_fragment = median(!! name.ncnt)) %>%
				mutate(cell_pct = round(cell_num / sum(cell_num) * 100, 2)) %>%
				rename("Clusters" = !! name.ident, "Cell Numbers" = cell_num, "Percentage of Cells (%)" = cell_pct,
						"Median Peaks per Cells" = median_peak, "Median fragments per Cells" = median_fragment)
		write.table(st, file = paste0(outdir, "/Stat.ClusterCells.xls"), sep = "\t", quote = F, row.names = F )

		invisible(st)
}

stat.ClusterSamples <- function(object, outdir, 
				assay = DefaultAssay(object),
				name.ident  = "seurat_clusters",
				name.sample = "orig.ident",
				name.nft    = paste0("nFeature_", assay),
				name.ncnt   = paste0("nCount_", assay),
				do.plot = TRUE
				) {
		name.ident  <- as.name(name.ident)
		name.sample <- as.name(name.sample)
		name.nft    <- as.name(name.nft)
		name.ncnt   <- as.name(name.ncnt)

		st1 <- object@meta.data %>%
				group_by(!! name.ident, !! name.sample) %>%
				summarise(count = n()) %>%
				rename(Clusters = !! name.ident, Samples = !! name.sample)

		st2 <- st1 %>% group_by(Samples) %>%
				mutate(count = paste0(count, " (", round(count / sum(count) * 100, 2), "%)")) %>%
				reshape2::dcast(Clusters ~ Samples, value.var = "count", fill = "0 (0%)")
		WriteTable(st2, paste0(outdir, "/Stat.ClusterSamples.xls"))

		if ( do.plot ) {
				p1 <- ggplot(st1, aes(x = Samples, y = count, fill = Clusters)) +
						geom_bar(stat = "identity", position = "stack") +
						ylab("Number of cells") +
						theme_light()
				ggsave(p1, file = paste0(outdir, "/Stat.ClusterSamples.barplot.pdf"), width = 7, height = 7)

				p2 <- ggplot(st1 %>% filter_all(any_vars(count != 0)), aes(x = Samples, y = Clusters, color = Clusters, size = count / sum(count) )) +
						geom_point() +
						guides(size = guide_legend(title = "Proportion"), color = FALSE) +
						theme_grey()
				ggsave(p2, file = paste0(outdir, "/Stat.ClusterSamples.dotplot.pdf"), width = 7, height = 7)
		}

		invisible(st1)
}


FindDifferPeaks <- function(object, thresh = 0.01, use.q = FALSE, test.use = 'LR', logFC = log(2), min.pct = 0.1, pseudocount.use = 0, idents.name = NULL, ignore.idents = NULL, latent.vars = 'atac_peak_region_fragments'){
		pvalue <- thresh
		if ( ! is.null(idents.name) ) Idents(object) <- idents.name
		if ( ! is.null(ignore.idents) ) object <- subset(object, idents = setdiff(levels(object), ignore.idents))
		if ( packageVersion("Seurat") < as.numeric_version("4.0.0") ) {
		obj.marker <- FindAllMarkers(object, only.pos = T, return.thresh = thresh, test.use = test.use, latent.vars = latent.vars,
				min.pct = min.pct, pseudocount.use = pseudocount.use, logfc.threshold = logFC)
		} else {
				obj.marker <- FindAllMarkers(object, only.pos = T, return.thresh = thresh, test.use = test.use, latent.vars = latent.vars,
						min.pct = min.pct, pseudocount.use = pseudocount.use, logfc.threshold = logFC, base = exp(1))
		}
		if ( use.q ) {
				obj.marker <- subset(obj.marker, p_val_adj < thresh)
		}
		return(obj.marker)
}

GetAnnotatePeak <- function(object, parameter = list(), annotation.file = parameter$ref$gtf,
				upstream.distal     = IfNull(parameter$peak.annot$upstream.distal, 100000),
				upstream.proximal   = IfNull(parameter$peak.annot$upstream.proximal, 2000),
				include.body        = IfNull(parameter$peak.annot$include.genebody, TRUE),
				downstream.distal   = IfNull(parameter$peak.annot$downstream.distal, 0),
				downstream.proximal = IfNull(parameter$peak.annot$downstream.proximal, 0)
				) {
		distal   <- AnnotatePeak(object, annotation.file = annotation.file, upstream = upstream.distal, downstream = downstream.distal, include.body = include.body, annotation.name = "distal", only.annot = FALSE)
		proximal <- AnnotatePeak(object, annotation.file = annotation.file, upstream = upstream.proximal, downstream = downstream.proximal, include.body = include.body, annotation.name = "proximal")
		rownames(distal)   <- paste(distal[["Peak Name"]], distal[["Gene ID"]], sep = "_")
		rownames(proximal) <- paste(proximal[["Peak Name"]], proximal[["Gene ID"]], sep = "_")
		Annot <- distal
		Annot[rownames(proximal),] <- proximal
		return(Annot)
}


AnnotatePeak <- function(object, annotation.file = NULL, include.body = TRUE, upstream = 2000, downstream = 0, annotation.name = c("proximal", "distal"), only.annot = TRUE, method = c("overlap", "nearest"), assay = DefaultAssay(object)){
		annotation.name <- match.arg(annotation.name)
		method <- match.arg(method)
		peak.matrix <- object[[assay]]
		peaks.gr <- StringToGRanges(rownames(x = peak.matrix), split_sep, starts.in.df.are.0based = TRUE)
#		BiocGenerics::start(peaks.gr[BiocGenerics::start(peaks.gr) == 0, ]) <- 1

		if ( is.null(Annotation(object)) ){
				gtf <- rtracklayer::import(con = annotation.file)
		} else {
				gtf <- Annotation(object)
		}
		gtf.genes <- gtf[gtf$type == "gene"]
		if ( length(gtf.genes) == 0 ) {
				gtf.genes <- gtf[gtf$type == "transcript"]
				if ( length(gtf.genes) == 0 ) gtf.genes <- gtf[gtf$type == "exon"]
				df <- as.data.frame(gtf.genes) %>%
						group_by(gene_id, gene_name, seqnames, strand) %>%
						summarise(start = min(start), end = max(end), width = end - start + 1, type = "gene", score = NA, phase = NA)
				gtf.genes <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns = T)
		}
		gtf.genes$gene_id[grepl("^ENS", gtf.genes$gene_id)] <- gsub("\\.\\d+$", "", gtf.genes$gene_id[grepl("^ENS", gtf.genes$gene_id)], perl = T)
		GenomeInfoDb::seqlevels(gtf.genes) <- gsub("_", "-", GenomeInfoDb::seqlevels(gtf.genes))
		if (include.body) {
				gtf.body_prom <- Extend(x = gtf.genes, upstream = upstream, downstream = downstream)
		} else {
				gtf.body_prom <- SummarizedExperiment::promoters(x = gtf.genes, upstream = upstream, downstream = downstream)
		}
		if ( method == "nearest" ){
				gene.distances <- GenomicRanges::distanceToNearest(x = peaks.gr, subject = gtf.body_prom)
				keep.overlaps <- gene.distances
		}else if ( method == "overlap" ) {
				keep.overlaps <- GenomicRanges::findOverlaps(query = peaks.gr, subject = gtf.body_prom)
		}

		gene.ids <- gtf.genes[S4Vectors::subjectHits(x = keep.overlaps)]
		peak.ids <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps)]
		tss <- GenomicRanges::promoters(x = gene.ids, upstream = 0, downstream = 1)

		genes <- as.data.frame(gene.ids) %>%
				dplyr::mutate(d2tss = GenomicRanges::distance(x = peak.ids, y = tss)) %>%
				dplyr::select("Gene Chr" = seqnames, "Gene Start" = start, "Gene End" = end, "Gene Length" = width, "Gene Strand" = strand, "Dist to TSS" = d2tss, "Gene ID" = gene_id)
		peaks <- as.data.frame(x = peak.ids) %>%
				dplyr::mutate(peak = rownames(peak.matrix)[S4Vectors::queryHits(x = keep.overlaps)], Annotation = annotation.name, start = start + 1) %>%
				dplyr::select("Peak Name" = peak, Chr = seqnames, Start = start, End = end, Annotation)
		annotation <- cbind(peaks, genes)

		if ( ! only.annot ) {
				unannot <- as.data.frame(x = peaks.gr[-S4Vectors::queryHits(x = keep.overlaps)]) %>%
						dplyr::mutate(peak = rownames(peak.matrix)[-S4Vectors::queryHits(x = keep.overlaps)], Annotation = "none", start = start + 1) %>%
						dplyr::select("Peak Name" = peak, Chr = seqnames, Start = start, End = end, Annotation) %>%
						dplyr::mutate("Gene Chr" = "-", "Gene Start" = "-", "Gene End" = "-", "Gene Length" = "-", "Gene Strand" = "*", "Dist to TSS" = "-", "Gene ID" = "-")
				annotation <- rbind(annotation, unannot)
		}
		annotation <- droplevels(annotation)

		return(annotation)
}


AnnotMotifs <- function(object, genome.file = NULL, motif.matrix = JASPAR2020::JASPAR2020, tax_group = NULL){
		if ( grepl("JASPAR", class(motif.matrix)) ) {
				pfm <- TFBSTools::getMatrixSet(x = motif.matrix, opts = list(tax_group = tax_group, all_versions = FALSE))
		} else {
				pfm <- TFBSTools::readJASPARMatrix(motif.matrix)
		}

		genome <- rtracklayer::import(genome.file, format = "fasta")
		names(genome) <- gsub(" .+", "", names(genome))
		names(genome) <- gsub("_", "-", names(genome))

		assay <- DefaultAssay(object)
		if ( any(grepl("-0-", rownames(object[[assay]]))) ) {
				index <- grepl("-0-", rownames(object[[assay]]))
				range <- GetAssayData(object, slot = "ranges", assay = assay)
				range[index] <- IRanges::narrow(range[index], start = 2)

				motif <- AddMotifs(object = range, genome = genome, pfm = pfm)
				rownames(motif@data) <- rownames(object[[assay]])
				object[[assay]] <- SetAssayData(object[[assay]], slot = "motifs", new.data = motif)

				feature.metadata <- RegionStats(object = range, genome = genome)
				rownames(feature.metadata) <- rownames(object[[assay]])
				meta.data <- GetAssayData(object = object, slot = "meta.features", assay = assay)
				feature.metadata <- feature.metadata[rownames(meta.data), ]
				meta.data <- cbind(meta.data, feature.metadata)
				object[[assay]] <- SetAssayData(object[[assay]], slot = "meta.features", new.data = meta.data)
		} else {
				object <- AddMotifs(object = object, genome = genome, pfm = pfm)
		}

		return(object)
}

FindDifferMotifs <- function(object, object.marker){
		differ.motifs <- by(object.marker$gene, object.marker$cluster, function(features){
				enriched.motifs <- FindMotifs(object = object, features = unique(features), background = NULL)
				enriched.motifs <- enriched.motifs %>%
						mutate(adj_p_val = p.adjust(pvalue)) %>%
						select(motif_ID = motif, motif_alt_TD = motif.name, observed, background, percent.observed, percent.background, fold.enrichment, pvalue, adj_p_val)
				return(enriched.motifs)}
		)
		for ( i in names(differ.motifs) ) {
				if ( ! is.null(differ.motifs[[i]]) && nrow(differ.motifs[[i]]) != 0 )
				differ.motifs[[i]] <- data.frame(Cluster = i, differ.motifs[[i]])
		}
		differ.motifs <- do.call(rbind, differ.motifs)
		return(differ.motifs)
}

GetPeakMotifAln <- function(motif.object){
		pwm <- motif.object@motif.names
		peak_motif_aln <- reshape2::melt(as.matrix(motif.object@data)) %>%
				rename("Peak Name" = Var1, "motif_ID" = Var2) %>%
				dplyr::filter(value > 0) %>% select(-value) %>%
				mutate(motif_alt_TD = pwm[motif_ID])
		return(peak_motif_aln)
}

GetTFGeneLink <- function(peak.motif, peak.gene, peak.cluster, motif.clutser){
		x <- left_join(x = peak.gene, y = peak.motif, by = "Peak Name")
		y <- peak.cluster %>% select(Cluster = cluster, "Peak Name" = gene) %>% left_join(y = x, by = "Peak Name")
		z <- semi_join(x = y, y = motif.clutser, by = c("Cluster", "motif_ID")) %>% 
				rename("TFBS Motif ID" = motif_ID) %>%
				select(-motif_alt_TD) ## may be same as GeneName. May be duplicated after add annot data
		return(z)
}

GetPeakOneGeneAln <- function(peak.annot){
		peak.annot <- peak.annot %>% filter(Annotation != "none")
		q <- StringToGRanges(peak.annot[["Peak Name"]], sep = split_sep)
		s <- StringToGRanges(transmute(peak.annot, regions = paste(`Gene Chr`, `Gene Start`, `Gene End`, sep = ":"))[[1]], sep = c(":",":"))
		peak.annot$dist <- GenomicRanges::distance(x = q, y = s)
		peak_one_gene_aln <- peak.annot %>% group_by(`Peak Name`) %>%
				filter(dist == min(dist)) %>% filter(`Dist to TSS` == min(`Dist to TSS`)) %>%
				select(-dist) %>% 
				as.data.frame()
		return(peak_one_gene_aln)
}




plot.CellCalling <- function(meta, outpref){
		meta <- cbind(meta, transmute(meta,
					cells = ifelse(is_cell == 1, "Cells", "Non-cells"),
					frt = atac_peak_region_fragments / atac_fragments)
				)  
		meta[order(meta$atac_peak_region_fragments, decreasing = T), "order"] <- 1:nrow(meta)
		meta <- subset(meta, rownames(meta) != "NO_BARCODE" & passed_filters > 0)

		color <- c("Cells" = "#4DAF4A", "Non-cells" = "grey")
		p3 <- ggplot(meta, aes(x = passed_filters, y = frt, color = cells)) +
				geom_point(size = 1, alpha = 0.7, shape = 19) +
				scale_x_continuous(trans = "log10", labels = c(1, 100, "10K", "1M"), breaks = 100 ^ c(0:3)) + 
				labs(x = "Fragments per Barcode", y = "Faction Fragments Overlapping Called Peaks", title = "Singlecell Targeting(Peaks)") +
				scale_color_manual(NULL, values = color) + 
				theme_minimal() + 
				theme(axis.line = element_line(color = "black"),
						axis.title = element_text(size = 15),
						axis.text = element_text(size = 10),
						plot.title = element_text(hjust = 0.5, size = 20),
						legend.text = element_text(size = 15) ) +
				guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1)))

		ggsave(filename = paste0(outpref, ".Scatter.pdf"), plot = p3, width = 7, height = 6)
}

PlotExpCorr <- function(atac_exp, rna_exp){

		get_density <- function (x, y, ...) {
				dens <- MASS::kde2d(x, y, ...)
				ix <- findInterval(x, dens$x)
				iy <- findInterval(y, dens$y)
				ii <- cbind(ix, iy)
				return(dens$z[ii])
		}

		common_id <- intersect(names(atac_exp), names(rna_exp))

		data <- data.frame(atac_exp = atac_exp[common_id], rna_exp = rna_exp[common_id])
		data.cor <- cor(data$rna_exp, data$atac_exp)
		cor_anno <- paste0("R = ", sprintf("%0.2f", data.cor))

		data$atac_exp <- log2(data$atac_exp)
		data$rna_exp  <- log2(data$rna_exp)
		data <- data[is.finite(data$atac_exp) & is.finite(data$rna_exp),,drop = FALSE]
		data$density <- get_density(data$atac_exp, data$rna_exp, n = 100)

		p <- ggplot(data) + 
				geom_point(aes(atac_exp, rna_exp, color = density), alpha = 0.6, size = 0.5) + 
				scale_color_viridis_c(option = "D") +
				labs(x = expression(log[2] * "(Average Activity Gene expression)"),
					 y = expression(log[2] * "(Average RNA Gene expression)"),
					 color = paste0(cor_anno, "\n\n", "Density of dots")) +
				dot_theme_default() +
				theme(legend.title = element_text(size = 14))

		return(p)
}

