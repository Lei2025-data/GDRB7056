DoNormalization <- function(object, parameter = list(), assay = "RNA", is_SCTransform = FALSE, vars.regress = NULL, is.check = TRUE, nfeatures = 2000, normalization.method = "LogNormalize", ...){
		DefaultAssay(object) <- assay

		### pre-deal with Cell Cycle
		object <- FindRegressVars(object, parameter, vars.regress = vars.regress, ...)

		if ( is_SCTransform ) {
				### Another workflow : SCTransform 
				object <- SCTransform(object, assay = assay, vars.to.regress = object@misc$vars.regress,
				                      verbose = FALSE, min_cells = 1, return.only.var.genes = F)
		}else{
				### Normalize Data
				message( "-->Normalize Data<--" )
				object <- NormalizeData(object, normalization.method = normalization.method, scale.factor = 10000)

				### Find Variable Features
				message( "-->Find Variable Genes<--" )
				object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = nfeatures)
				if ( is.check ) CheckVariableFeature(object)

				### Confound factors and Scale Data
				message( "-->ScaleData<--" )
				object <- ScaleData(object = object, vars.to.regress = object@misc$vars.regress, features = rownames(object))
		}
		return(object)
}

CheckVariableFeature <- function(object){
		top10 <- head(VariableFeatures(object), 10)
		plot1 <- VariableFeaturePlot(object)
		plot2 <- LabelPoints(plot = plot1, points = top10, labels = object@misc$fdata[top10,"merge_name"], repel = TRUE, xnudge = 0, ynudge = 0)
		plot2 <- plot2 + theme(legend.position = "top")
		ggsave(plot2, file = "Variable_gene.pdf", width = 6, height = 6)

		write.table(VariableFeatures(object), file = "var_gene.xls", quote = F, sep = "\t", col.name = F)
}


FindRegressVars <- function(object, parameter = list(), vars.regress = NULL, force_recal = FALSE,
				is_rm_cc = IfNull(parameter$cell_cycle$is_remove, TRUE),
				is_rm_all_cc_signal = IfNull(parameter$cell_cycle$is_rm_all_signal, FALSE),
				...)
{
		if ( is.null(vars.regress) ) {
				vars.regress <- paste0("nCount_", DefaultAssay(object)) # nCount_RNA
				if ( exists("percent.mito", object@meta.data) ) {
						vars.regress <- c(vars.regress, "percent.mito")
				}
				if ( is_rm_cc ) {
						message("-->CellCycle Scoring<--")
						if ( ! exists("CC.Difference", object@meta.data) || force_recal ) {
								object <- DoCellCycleScoring(object, ...)
						}
						if ( is_rm_all_cc_signal ) {
								vars.regress <- c(vars.regress, "S.Score", "G2M.Score")
						} else {
								vars.regress <- c(vars.regress, "CC.Difference")
						}
				}
		}else if ( vars.regress == "none" ){
				vars.regress <- NULL
		}
		object@misc$vars.regress <- vars.regress
		return(object)
}

DoCellCycleScoring <- function(object, s.genes = NULL, g2m.genes = NULL, ...){
		if ( is.null(s.genes)   ) s.genes   <- cc.genes$s.genes
		if ( is.null(g2m.genes) ) g2m.genes <- cc.genes$g2m.genes
		s.genes <- FindFeaturesID(object, s.genes)
		g2m.genes <- FindFeaturesID(object, g2m.genes)
		if ( length(s.genes) < 2 || length(g2m.genes) < 2 ) return(object)
		object <- CellCycleScoring(object, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE, ...)
		object[["CC.Difference"]] <- object[["S.Score"]] - object[["G2M.Score"]]
		return(object)
}


