
FilterGenes <- function(object, parameter = list(),
				min.cell = IfNull(parameter$filter$min.cell, 0)
				) {
		message( "--->Filter Genes<---" ) 
		if ( min.cell > 0 ){
				if ( min.cell < 1 ) {
						min.cell <- min.cell * length(object@cell.names)
				}
				num.cells  <- Matrix::rowSums( object@assays$RNA@counts > 0 )
			
				genes.filter <- num.cells[which( num.cells < min.cell )]
				genes.filter <- FindFeaturesName(object, genes.filter)				
				write.table( genes.filter, file = "filtered_genes.xls", quote = F, sep = "\t", col.name = F )
				
				genes.use  <- names( num.cells[which( num.cells >= min.cell )] )
				object     <- object[genes.use, ]

				## TODO : need to re-calculate 'percent.mito', 'expected.marker', ... value ?
		}
		return(object)
}

.FilterCells <- function(object, standard = NULL, set.num = "none", set.num.seed = 42) {
		if ( class(object) == "Seurat" ) {
				metadata <- object@meta.data
				metadata$cell <- rownames(metadata)
		}else{
				metadata <- object
				metadata$orig.ident <- metadata$Sample
		}
		cells.use <- metadata$cell
#		cells.use <- Cells(object)
		if ( !is.null(set.num) && set.num != "none" ){
				if ( set.num == "min" ){
						cell_num <- min(table(metadata$orig.ident))
				} else if ( set.num != "none" ) {
						cell_num <- min( as.integer(set.num), max(table(metadata$orig.ident)) )
				}
				seed <- set.num.seed
#				set.seed(seed)
				cells.use <- as.character(unlist(by(cells.use, metadata$orig.ident, function(x) sample(x, min(cell_num, length(x))))))
		}
		pm.used <- list()
		for ( i in names(standard) ) {
				if ( exists(i, metadata) ) {
						if ( length(standard[[i]]) == 1 ) {
								value <- standard[[i]] ## just keeping same form as below
								if ( value == "auto" ) {
										standard[[i]] <- autothres(data = metadata[[i]], name = i, bin = 100)
										print(standard[[i]])
								}else{
										cells.use <- metadata %>%
#												tibble::rownames_to_column(var = "cells") %>%
												filter(.data[[i]] == value & cell %in% cells.use) %>%
												select(cell) %>% unlist()
								}
						}
						if ( length(standard[[i]]) == 2 ) {
								lower <- standard[[i]][1]
								upper <- standard[[i]][2]
								cells.use <- metadata %>%
#										tibble::rownames_to_column(var = "cells") %>%
										filter(.data[[i]] >= lower & .data[[i]] <= upper & cell %in% cells.use) %>%
										select(cell) %>% unlist()
						}
						pm.used[[i]] <- standard[[i]]
				}
		}
		yaml::write_yaml(pm.used, file = "filtered_used_parameter.yaml")
		return(cells.use)
}


FilterCells <- function(object, parameter, do.stat = TRUE,
				set.num = IfNull(parameter$filter$set.num, "none"),
				set.num.seed = IfNull(parameter$filter$set.num.seed, 42),
				standard = parameter$filter$standard
				) {
		message( "--->Filter Cells<---" )
		cells.use <- .FilterCells(object, set.num = set.num, set.num.seed = set.num.seed, standard = standard )
		object <- object[, cells.use]

		if ( do.stat ) {
				StatFilterCells(object)
		}

		return(object)
}


StatFilterCells <- function(object, group.by = "orig.ident", outfile = "Filter.stat.xls"){
		filter_cells <- setdiff(rownames(object@misc$pdata), colnames(object))
		write.table( filter_cells, file = "filtered_cells.xls", quote = F, sep = "\t", col.name = F )

		name <- if ( group.by == "orig.ident" ) "Samples" else group.by
		group.by <- as.name(group.by)
		if ( exists("percent.mito", object@meta.data) && exists("percent.mito", object@misc$pdata) ) {
				a <- object@meta.data %>% group_by(!! name := !! group.by) %>%
						summarise(after_filter_num = n(),
								after_filter_median_UMI_per_cell   = median(nCount_RNA),
								after_filter_median_genes_per_cell = median(nFeature_RNA),
								after_filter_median_MT_per_cell    = median(percent.mito))
				b <- object@misc$pdata %>% group_by(!! name := !! group.by) %>%
						summarise(before_filter_num = n(),
								before_filter_median_UMI_per_cell   = median(nCount_RNA),
								before_filter_median_genes_per_cell = median(nFeature_RNA),
								before_filter_median_MT_per_cell    = median(percent.mito))
				filter_stat_table <- left_join(a,b) %>%
						mutate(pct = paste0( round( after_filter_num/before_filter_num * 100, 2) , "%")) %>%
						select(!!name, before_filter_num, after_filter_num, pct,
								before_filter_median_UMI_per_cell,   after_filter_median_UMI_per_cell,
								before_filter_median_genes_per_cell, after_filter_median_genes_per_cell,
								before_filter_median_MT_per_cell,    after_filter_median_MT_per_cell)
		}else{
				a <- object@meta.data %>% group_by(!! name := !! group.by) %>%
						summarise(after_filter_num = n(),
								after_filter_median_UMI_per_cell   = median(nCount_RNA),
								after_filter_median_genes_per_cell = median(nFeature_RNA))
				b <- object@misc$pdata %>% group_by(!! name := !! group.by) %>%
						summarise(before_filter_num = n(),
								before_filter_median_UMI_per_cell   = median(nCount_RNA),
								before_filter_median_genes_per_cell = median(nFeature_RNA))
				filter_stat_table <- left_join(a,b) %>%
						mutate(pct = paste0( round( after_filter_num/before_filter_num * 100, 2) , "%")) %>%
						select(!!name, before_filter_num, after_filter_num, pct,
								before_filter_median_UMI_per_cell,   after_filter_median_UMI_per_cell,
								before_filter_median_genes_per_cell, after_filter_median_genes_per_cell)	
		}
		write.table( filter_stat_table, file = outfile, quote = F, sep = "\t", row.names = F  )
}


