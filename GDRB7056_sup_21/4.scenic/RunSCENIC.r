#! /Bio/User/yinquan/software/miniconda3/bin/Rscript

### Deal arguements
args <- commandArgs()
bin  <- dirname(normalizePath(sub('--file=', '',  args[grep('--file=', args)])))
args <- args[-seq(grep("--args", args))]

para <- args[1]
if ( is.null( para ) | is.null( para ) ){
        warning( "\n\tUsage : /Bio/User/yinquan/software/miniconda3/bin/Rscript runSCENIC.r runSCENIC.yaml \n" )
        quit()
}


## Source functions
source(paste0(bin, "/SCENIC_lib.R"))

## library packages
packages <- c("AUCell", "yaml","pheatmap", "dplyr", "reshape2", "ggplot2", "Seurat", "SCENIC", "RColorBrewer", "magrittr","NetPathMiner","sigmaNet")
LibraryPackages(packages, lib.loc="/Bio/User/yinquan/software/miniconda3/lib/R/library")

## parameter
RunMessage(">>>> scan parameter ... <<<<")
parameter <- yaml.load_file(para)

## setwd()
RunMessage(">>>> Set output directory ... <<<<")
setwd(parameter$outdir)
Message(paste0("outdir: ", parameter$outdir))

RunMessage(">>>> Attach SCENIC Data... <<<<")
dir <- dirname(parameter$outdir)
data <- AttachScenic(dir)
scenicOptions <- data[["scenicOptions"]]
exprMat <- data[["exprMat"]]
if(is.null(scenicOptions)) StopMessage(paste0(dir, "/1.CreateScenicOptions/scenicOptions.Rds isn't existed!"))
if(is.null(exprMat)) StopMessage(paste0(dir, "/2.FilterGene/exprMat.Rds isn't existed!"))
genesKept  <- loadInt(scenicOptions, "genesKept")
exprMat_filtered <- exprMat[genesKept, ]

RunMessage(">>>> Run SCENIC ... <<<<")
exprMat_filtered_log <- log2(exprMat_filtered+1)
RunMessage("    ---> Run SCENIC_1 coexNetwork2modules ... <---")
runSCENIC_1_coexNetwork2modules(scenicOptions)
RunMessage("    ---> Run SCENIC_2 createRegulons ... <---")
runSCENIC_2_createRegulons(scenicOptions, minGenes = parameter$minGenes)
RunMessage("    ---> Run SCENIC_3 scoreCells ... <---")
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log, skipBinaryThresholds = FALSE, skipHeatmap = FALSE, skipTsne = FALSE)
RunMessage("    ---> Binary regulon Activity ... <---")
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_filtered_log, skipBoxplot=FALSE, skipHeatmaps=FALSE, skipTsne=FALSE)

RunMessage(">>>> Do Plots ... <<<<")
RunMessage("    ---> Do Auc activity heatmap analysis  ... <---")
cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions,'cellInfo'))
group_by <- NULL
if(exists("cluster", cellInfo)) group_by <- c(group_by, "cluster")
if(exists("sample", cellInfo)) group_by <- c(group_by, "sample")
if(exists("group", cellInfo)) group_by <- c(group_by, "group")
RunRegulonHeatmap(scenicOptions, group_by="cluster")
RunRegulonHeatmap(scenicOptions, group_by="sample")
RunRegulonHeatmap(scenicOptions, group_by="group")
RunRegulonHeatmap(scenicOptions, group_by=group_by)

RunMessage("    ---> TSNE plot ... <---")
TSNEplot(scenicOptions)
RunMessage("    ---> Feature  plot ... <---")
FeaturePlots(scenicOptions, features.type = "AUC", cols.use= c("lightgrey", "red"))
FeaturePlots(scenicOptions, features.type = "expression", cols.use= c("lightgrey", "blue"))
FeaturePlots(scenicOptions, features.type = "BinaryAUC", cols.use= c("lightgrey", "#104E8B"))

RunMessage("    ---> Export regulon motif enrichment file  ... <---")
motifEnrichment <- read.table(getOutName(scenicOptions,'s2_motifEnrichment'), sep="\t", header=T, stringsAsFactors=F)
motifEnrichment <- motifEnrichment[,c("motif","motifDb","nEnrGenes","enrichedGenes","highlightedTFs","TFinDB","TF_highConf","TF_lowConf")]
colnames(motifEnrichment) <- c("Motif","motifDB","nEnrGenes","enrichedGenes","highlightedTFs","TFinDB","TF_highConf","TF_lowConf")
write.table(motifEnrichment, file="motifEnrichment.xls", sep="\t", quote=F, col.names=T, row.names=F)

RunMessage("    ---> Export regulon targets and generate TF-target network files  ... <---")
getRegulon2Targets(scenicOptions)
DoNetwork(scenicOptions, ntop=2)


unlink("output", recursive = T)
RunMessage(">>>> Done <<<<")
