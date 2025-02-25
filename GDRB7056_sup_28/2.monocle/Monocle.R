cores <- 1
#if ( as.numeric(R.version$minor) >= 6 ){
#		.libPaths("/home/xushuyang/R/x86_64-unknown-linux-gnu-library/3.6")
#} else if ( as.numeric(R.version$minor) >= 5 ) {
#		.libPaths("/home/xushuyang/R/x86_64-unknown-linux-gnu-library/3.5")
#} else {
#		cores <- 4
#}

library(monocle)
library(igraph)
library(dplyr)

args <- commandArgs(T)

conf.file <- args[1]
outdir    <- args[2]
add_lib   <- args[3]
cls_col_name <- ifelse(!is.na(args[4]), args[4], "seurat_clusters")

source(add_lib, chdir = T)

setwd(outdir)


##########################################################
################### pipe line start here #################
handlers <- list("bool#no" = function(x){if ( x %in% c("false", "FALSE") ) FALSE else x}, "bool#yes" = function(x){if ( x %in% c("true", "TRUE") ) TRUE else x})
parameter <- yaml::yaml.load_file( conf.file, handlers = handlers)


### ready obj 
message( ">>> ready obj <<<" )
message( "==> Initial monocle <==" )
#mnc_obj <- MonocleObject(parameter, cls_col_name)
load(parameter$obj_file)
obj@meta.data$Samples <- obj@meta.data$orig.ident
obj@meta.data$Clusters <- obj@meta.data$seurat_clusters
cells.use   <- rownames(obj@meta.data)
cells.use <- sample(cells.use, size = 38733)
mnc_obj <- Initial( obj, cells.use)

if ( is.null(parameter$mnc_obj) || ! file.exists(parameter$mnc_obj) ) {
#if ( "run_monocle" %in% parameter$analyse ) {
### in my first design, when do second analysis, "run_monocle" will be turn off. 
### but product manager didn't think so. so it be this way.
		message( "==> run monocle <==" )
		mnc_obj <- RunMonocle(mnc_obj, parameter)
}

mnc_obj <- RenameState(mnc_obj, parameter)
saveRDS( mnc_obj, file = "Trajectory.obj.rds" )

message( "==> stat Trajectory <== " )
StatTrajectory(mnc_obj)
GetTrajectoryData(mnc_obj) # for online
StatCluster(pData(mnc_obj), group.by = "Samples", stat.what = "State", outpref = "State.stat" )
StatCluster(pData(mnc_obj), group.by = "Clusters", stat.what = "State", outpref = "State.stat" )
CDS_avg(mnc_obj)

### Diff
message( ">>> Differ <<<" )
if ( "diff_state" %in% parameter$analyse && nlevels(pData(mnc_obj)$State) > 1 ) {
		message( "==> state diff <== " )
		setwd(outdir)
		dir.create(path = "Diff_State", recursive = T, showWarnings = F)
		setwd("Diff_State")

		diff_state_res <- Diff.State(mnc_obj, cores = cores, thres = parameter$diff_state$thres, use.q = parameter$diff_state$use.q)

		GetGenesInPseudotime(mnc_obj[rownames(diff_state_res)[1], ], cores = cores)
		setwd(outdir)
}

if ( "diff_pseudotime" %in% parameter$analyse ) {
		message( "==> pseudotime diff <== " )
		setwd(outdir)
		dir.create(path = "Diff_Pseudotime", recursive = T, showWarnings = F)
		setwd("Diff_Pseudotime")

		diff_Pseudotime_sig <- Diff.Pseudotime(mnc_obj, cores = cores, thres = parameter$diff_pseudotime$thres, use.q = parameter$diff_pseudotime$use.q)

		GetGenesInPseudotime(mnc_obj[rownames(diff_Pseudotime_sig)[1], ], cores = cores)
		setwd(outdir)
}

if ( "diff_branch" %in% parameter$analyse && nlevels(pData(mnc_obj)$State) > 1) {
		message( "==> branch diff <== " )
		setwd(outdir)
		dir.create(path = "Diff_Branch", recursive = T, showWarnings = F)
		setwd("Diff_Branch")

		total_BEAM_sig <- Diff.Branch(mnc_obj, cores = cores, thres = parameter$diff_branch$thres, use.q = parameter$diff_branch$use.q)

		branch_name <- strsplit(as.character(total_BEAM_sig$branch[1]), " -vs- ")[[1]]
		GetGenesBranchedPseudotime(mnc_obj[as.character(total_BEAM_sig$GeneID[1]), ], branch_point = total_BEAM_sig$branch_node[1], branch_labels = branch_name, cores = cores)
		setwd(outdir)
}


message( ">>> __ALL DONE__ <<<" )
system("touch _complete")



#tip_leaves <- names(which(degree(minSpanningTree(cds)) == 1))


