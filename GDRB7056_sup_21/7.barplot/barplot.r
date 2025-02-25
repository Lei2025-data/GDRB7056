#! /state/partition1/WORK/lib-devel/working_library/anaconda2/envs/R-3.5.1/bin/Rscript
library(ggplot2)
args <- commandArgs(T)
if (length(args) < 6){
	print ("/state/partition1/WORK/lib-devel/working_library/anaconda2/envs/R-3.5.1/bin/Rscript barplot.r <infile> <outprefix> <xcol> <ycol> <xlab> <ylab>")
	q()
}

file <- args[1]
out <- args[2]
xcol <- as.numeric(args[3])
ycol <- as.numeric(args[4])
xlab <- args[5]
ylab <- args[6]
gcol <- as.numeric(args[7])
range <- args[8]

dat <- read.table (file, header = T, sep = "\t")
range <- strsplit(range, "-")[[1]]
dat <- dat[range[1]:range[2],]
colnames(dat)[xcol] <- "x"
colnames(dat)[ycol] <- "y"
#dat <- dat[order(dat$y, decreasing= T),]
#dat[,xcol] <- factor(dat[,xcol], levels = dat[,xcol])
dat[,gcol] <- factor(dat[,gcol], levels = c("WT", "Aa", "KO"))
if (length(rownames(dat))/4 < 7){
	pdf(paste0(out,".barplot.pdf"))
}else{
	pdf(paste0(out,".barplot.pdf"), width = length(rownames(dat))/5)
}
if (is.null(gcol)){
ggplot(dat, aes(x = x, y = -log10(y))) + geom_bar(stat = "identity", fill = "orange", alpha = -log10(dat$y)/14) + labs(x = xlab, y = ylab) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_discrete(labels = paste(dat$x, dat[,2], sep = ":"))
}else{
colnames(dat)[gcol] <- 'group'
ggplot(dat, aes(x = x, y = y, fill = group)) + geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.5) + labs(x = xlab, y = ylab) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}
