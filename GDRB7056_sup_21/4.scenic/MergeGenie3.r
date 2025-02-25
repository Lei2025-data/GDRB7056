#! /Bio/User/yinquan/software/miniconda3/bin/Rscript

### Deal arguements
args <- commandArgs()
bin  <- dirname(normalizePath(sub('--file=', '',  args[grep('--file=', args)])))
args <- args[-seq(grep("--args", args))]

outdir <- args[1]
if ( is.na( outdir )){
        warning( "\n\tUsage : /Bio/User/yinquan/software/miniconda3/bin/Rscript mergeGenie3.r <outdir> \n" )
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
if(is.null(scenicOptions)) StopMessage(paste0(dir, "/1.CreateScenicOptions/scenicOptions.Rds isn't existed!"))

RunMessage(">>>> Merge Genie3 ... <<<<")
MergeGenie3(scenicOptions)

RunMessage(">>>> Done ... <<<<")
