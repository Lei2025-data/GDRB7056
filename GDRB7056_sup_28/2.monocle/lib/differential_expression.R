
DoFindAllMarkers <- function(object, parameter = list(), group.by = "seurat_clusters",
				min.pct = IfNull(parameter$FindMarkers$min_pct, 0.25),
				logfc.threshold = IfNull(parameter$FindMarkers$logfc, 0.25),
				return.thresh = IfNull(parameter$FindMarkers$pvalue, 0.01),
				pseudocount.use = 0, only.pos = TRUE, ...)
{
		Idents(object) <- group.by
		object.markers <- FindAllMarkers(object = object, only.pos = only.pos,
						min.pct = min.pct, logfc.threshold = logfc.threshold, 
						return.thresh = return.thresh, pseudocount.use = pseudocount.use, ...)
		if ( nrow(object.markers) == 0 ) {
				object.markers <- data.frame("p_val" = 1, "avg_logFC" = 1, "pct.1" = 1, "pct.2" = 1, "p_val_adj" = 1, "cluster" = "1", "gene" = 1)
				object.markers <- object.markers[-1,]
		}
		return(object.markers)
}



ListMarker <- function(object, object.markers, outfile = "DeGene.list.xls", is.return = FALSE, is.fast = FALSE, group.by = "seurat_clusters",
				assay = DefaultAssay(object), slot = "data", group.by.name = "Cluster"){
		Targets_name <- paste0( "Target_", group.by.name )
		Others_name  <- paste0( "Other_",  group.by.name )
		Targets_mean_name <- paste0( Targets_name, "_mean" )
		Others_mean_name  <- paste0( Others_name,  "_mean" )		
		Targets_mean <- CalAvgExp(object, unique(object.markers$gene), is.return = T, is.reverse = F, assay = assay, slot = slot, group.by = group.by) %>%
				reshape2::melt(varnames = c("gene", "cluster"), value.name = Targets_mean_name) %>%
				mutate(cluster = factor(cluster))
		Others_mean  <- CalAvgExp(object, unique(object.markers$gene), is.return = T, is.reverse = T, assay = assay, slot = slot, group.by = group.by) %>%
				reshape2::melt(varnames = c("gene", "cluster"), value.name = Others_mean_name) %>%
				mutate(cluster = factor(cluster))
		Name <- FindFeaturesName(object, unique(object.markers$gene), "name", is.fast = is.fast)
		marker_list <- object.markers %>% left_join(y = Targets_mean) %>% left_join(y = Others_mean) %>% 
				mutate(Log2FC = log2(!! as.name(Targets_mean_name) / !! as.name(Others_mean_name)), name = Name[gene]) %>%
				select(!! Targets_name := cluster, "Gene ID" = gene, "Gene Name" = name,
						!! Targets_mean_name, !! Others_mean_name, Log2FC,
						Pvalue = p_val, Qvalue = p_val_adj)
		if ( is.null(outfile) || is.return ){
				return(marker_list)
		}else{
				marker_list[["Gene ID"]] <- ChangeOUTName(marker_list[["Gene ID"]], object@misc$fdata)
				WriteTable(marker_list, outfile)
		}
}

FindTopMarker <- function(object.markers, top_num = 20, object = NULL, outfile = "Top.avg_exp.xls.tmp"){
		top <- object.markers %>% group_by( cluster ) %>%
#				top_n(top_num, avg_logFC)
				arrange(desc(avg_logFC), p_val, p_val_adj, .by_group = TRUE) %>% filter(1:n() <= top_num)
		if ( ! is.null(outfile) ) {
				tmp <- top %>% select(Cluster = cluster, Gene_ID = gene)
				if ( ! is.null(object) ) {
						tmp$Gene_ID <- ChangeOUTName(tmp$Gene_ID, object@misc$fdata)
				}
				WriteTable(tmp, outfile )
		}
		return(top)
}

