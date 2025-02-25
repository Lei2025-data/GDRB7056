#! /Bio/User/yinquan/software/miniconda3/bin/Rscript

### Deal arguements
args <- commandArgs()
bin  <- dirname(normalizePath(sub('--file=', '',  args[grep('--file=', args)])))
args <- args[-seq(grep("--args", args))]

para <- args[1]
if ( is.null( para ) | is.na( para ) ){
        warning( "\n\tUsage : /Bio/User/yinquan/software/miniconda3/bin/Rscript CreateScenicOptions.r <scenicoption.yaml>\n" )
        quit()
}

## Source functions
source(paste0(bin, "/SCENIC_lib.R"))

## library packages
packages <- c("AUCell", "yaml","pheatmap", "dplyr", "reshape2", "ggplot2", "Seurat", "SCENIC", "RColorBrewer")
LibraryPackages(packages, lib.loc="/Bio/User/yinquan/software/miniconda3/lib/R/library")

## approach parameter
RunMessage(">>>> Scan parameter ... <<<<")
parameter <- yaml.load_file( para )

## setwd()
RunMessage(">>>> Set output directory ... <<<<")
setwd(parameter$outdir)
Message(paste0("outdir: ", parameter$outdir))
#dir.create("int", showWarnings = F)

## load object
obj <- LoadObject(parameter)
saveRDS(obj, file="obj.Rds")
head(obj@meta.data)
levels(obj@active.ident)
## initial scenicOptions
scenicOptions <- CreateScenicOptions(parameter)

## set colors
RunMessage(">>>> Set Senic Colors ... <<<<")
cellInfo <- NULL
## sample color
Message("     ---> sample color <---")
color.sample <- rainbow(length(levels(obj@meta.data$sample)))
names(color.sample) <- levels(obj@meta.data$sample)
## cluster color
Message("     ---> cluster color <---")
color.cluster <- rainbow(length(levels(obj@meta.data$cluster)))
names(color.cluster) <- levels(obj@meta.data$cluster)
cellInfo <- obj@meta.data[,c("sample", "cluster")]

RunMessage(">>>> Set tsne perplexity ... <<<<")
perplexity <- floor((length(Cells(obj)) - 1)/3)
if (perplexity >= 50){
	perplexity <- 50
}
scenicOptions@settings$defaultTsne$perpl <- perplexity

## group color
color.group <- NULL
if(! is.null(parameter$metadata$groups$orig.label)){
        Message("     ---> group color <---")
        color.group <- rainbow(length(levels(obj@meta.data$group)))
        names(color.group) <- levels(obj@meta.data$group)
        cellInfo <- cbind(cellInfo, group = obj@meta.data$group)
}
## save cellInfo, colVars and scenicOptions
RunMessage(">>>> Save cellInfo colVars scenicOptions ... <<<<")
saveRDS(cellInfo, file="cellInfo.Rds")
colVars <- NULL
if(!is.null(color.group)){
        colVars <- list(sample = color.sample, cluster = color.cluster, group = color.group)
}else{
        colVars <- list(sample = color.sample, cluster = color.cluster)
}
saveRDS(colVars, file="colVars.Rds")
scenicOptions@inputDatasetInfo$cellInfo <- paste0(parameter$outdir, "/cellInfo.Rds")
scenicOptions@inputDatasetInfo$colVars <- paste0(parameter$outdir, "/colVars.Rds")
scenicOptions@fileNames$dir <- dirname(parameter$outdir)
saveRDS(scenicOptions, file="scenicOptions.Rds")

## export group.list sample.list cluster.list
dir <- dirname(parameter$outdir)
dir <- paste0(dir, "/report_conf")
setwd(dir)
if(exists("group", cellInfo)){
	group.list <- unique(arrange(cellInfo, group)[, c("sample", "group")])
	write.table(group.list, file="group.list", sep="\t", col.names=F, row.names=F, quote=F)
}

if(length(levels(cellInfo$sample)) > 1) write.table(levels(cellInfo$sample), file="sample.list", quote=F, col.names=F, row.names=F)
if(length(levels(cellInfo$cluster)) > 1) write.table(levels(cellInfo$cluster), file="cluster.list", quote=F, col.names=F, row.names=F)

RunMessage(">>>> Done ... <<<<")
