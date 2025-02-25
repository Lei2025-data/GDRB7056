#! /Bio/User/yinquan/software/miniconda3/bin/Rscript

### Deal arguements
args <- commandArgs()
bin  <- dirname(normalizePath(sub('--file=', '',  args[grep('--file=', args)])))
args <- args[-seq(grep("--args", args))]

outdir <- args[1]
glist <- args[2]
index <- unlist(strsplit(basename(glist), split="-"))[2]
if ( is.na( outdir ) | is.na( glist ) ){
        warning( "\n\tUsage : /Bio/User/yinquan/software/miniconda3/bin/Rscript FilterGene.r <outdir> <gene list> <index>\n" )
        quit()
}


## Source functions
source(paste0(bin, "/SCENIC_lib.R"))

## library packages
packages <- c("AUCell", "yaml","pheatmap", "dplyr", "reshape2", "ggplot2", "Seurat", "SCENIC", "RColorBrewer")
LibraryPackages(packages, lib.loc="/Bio/User/yinquan/software/miniconda3/lib/R/library")

## setwd()
RunMessage(">>>> Set output directory ... <<<<")
setwd(outdir)
Message(paste0("outdir: ", outdir))

RunMessage(">>>> Attach SCENIC Data... <<<<")
dir <- dirname(outdir)
data <- AttachScenic(dir)
scenicOptions <- data[["scenicOptions"]]
exprMat <- data[["exprMat"]]
if(is.null(scenicOptions)) StopMessage(paste0(dir, "/1.CreateScenicOptions/scenicOptions.Rds isn't existed!"))
if(is.null(exprMat)) StopMessage(paste0(dir, "/2.FilterGene/exprMat.Rds isn't existed!"))
genesKept  <- loadInt(scenicOptions, "genesKept")
exprMat_filtered <- exprMat[genesKept, ]

RunMessage(">>>> RunGenie3 ... <<<<")
genes.use <- readLines(con=glist)
exprMat_filtered_log <- log2(exprMat_filtered+1)
RunGenie3(exprMat_filtered_log, scenicOptions, genes.use, index, nCores=2)

RunMessage(">>>> Done ... <<<<")
