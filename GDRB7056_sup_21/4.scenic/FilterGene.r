#! /Bio/User/yinquan/software/miniconda3/bin/Rscript

### Deal arguements
args <- commandArgs()
bin  <- dirname(normalizePath(sub('--file=', '',  args[grep('--file=', args)])))
args <- args[-seq(grep("--args", args))]

para <- args[1]
if ( is.null( para ) | is.na( para ) ){
        warning( "\n\tUsage : /Bio/User/yinquan/software/miniconda3/bin/Rscript FilterGene.r <FilterGene.yaml>\n" )
        quit()
}

## Source functions
source(paste0(bin, "/SCENIC_lib.R"))

## library packages
packages <- c("AUCell", "Rcpp","yaml","pheatmap", "dplyr", "reshape2", "ggplot2", "Seurat", "SCENIC", "RColorBrewer")
LibraryPackages(packages, lib.loc="/Bio/User/yinquan/software/miniconda3/lib/R/library")

#Rcpp::sourceCpp(code='
#include <Rcpp.h>
#using namespace Rcpp;
#// [[Rcpp::export]]
#NumericMatrix asMatrix(NumericVector rp, NumericVector cp,NumericVector z,int nrows,int ncols){
#        int k = z.size();
#        NumericMatrix  mat(nrows, ncols);
#         for (int i = 0; i < k; i++){
#                mat(rp[i],cp[i]) = z[i];
#        }
#        return mat;
#}
#' )
#as_matrix <- function(mat){
#        row_pos <- mat@i
#        col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])
#        tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,nrows =  mat@Dim[1], ncols = mat@Dim[2])
#        row.names(tmp) <- mat@Dimnames[[1]]
#        colnames(tmp) <- mat@Dimnames[[2]]
#        return(tmp)
#}




## approach parameter
RunMessage(">>>> Scan parameter ... <<<<")
parameter <- yaml.load_file( para )

## setwd()
RunMessage(">>>> Set output directory ... <<<<")
setwd(parameter$outdir)
Message(paste0("outdir: ", parameter$outdir))

RunMessage(">>>> Attach SCENIC Data... <<<<")
dir <- dirname(parameter$outdir)
data <- AttachScenic(dir)
scenicOptions <- data[["scenicOptions"]]
obj <- data[["object"]]
if(is.null(scenicOptions)) StopMessage(paste0(dir, "/1.CreateScenicOptions/scenicOptions.Rds isn't existed!"))
if(is.null(obj)) StopMessage(paste0(dir, "/1.CreateScenicOptions/obj.Rds isn't existed!"))

RunMessage(">>>> Filtering genes ... <<<<")
exprMat <- as.matrix(obj@misc$Expr)
#exprMat <- as_matrix(obj@misc$Expr)
saveRDS(exprMat,file="exprMat.Rds")

minCountsPerGene <- ifelse(is.null(parameter$CountsPerGene), parameter$CountsPerGene.pct * 3 * ncol(exprMat), parameter$CountsPerGene) 
genesKept <- geneFiltering(exprMat,
							scenicOptions = scenicOptions,
							minCountsPerGene = minCountsPerGene,
							minSamples= ncol(exprMat) * parameter$Samples.pct)

exprMat_filtered <- exprMat[genesKept, ]
### run Corr
RunMessage(">>>> Run Correlation ... <<<<")
runCorrelation(exprMat_filtered, scenicOptions)

RunMessage(">>>> split genes ... <<<<")
genesSplit <- suppressWarnings(split(sort(rownames(exprMat_filtered)), 1: parameter$split))
for (i in 1:length(genesSplit)) {
	genes <- genesSplit[[i]]
	write.table(genes, sep="\t", quote=F, col.names=F, row.names=F, file= paste0( "genesplit-", i))
}

#genes <- sort(rownames(exprMat_filtered))
#GeneSplit(genes, splitnum=20)

RunMessage(">>>> Done <<<<")
