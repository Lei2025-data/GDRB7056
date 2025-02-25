#! /state/partition1/WORK/lib-devel/working_library/anaconda2/envs/R-3.5.1/bin/Rscript
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
if (F){
cc_annot <- readRDS("/public2/Bio/Project/GDRB7056/GDRB7056_sup_17/2.sup14/pipe/7.CellCycle/CellCycle.annot.Rds")
cc_annot <- cc_annot[cc_annot$Clusters %in% c(0:8,12,13,14,17,18,23),]
for (i in 0:8){
	cc_annot$Clusters_new[cc_annot$Clusters == i] <- i
}
cc_annot$Clusters_new[cc_annot$Clusters == 12] <- 9
cc_annot$Clusters_new[cc_annot$Clusters == 13] <- 10
cc_annot$Clusters_new[cc_annot$Clusters == 14] <- 11
cc_annot$Clusters_new[cc_annot$Clusters == 17] <- 12
cc_annot$Clusters_new[cc_annot$Clusters == 18] <- 13
cc_annot$Clusters_new[cc_annot$Clusters == 23] <- 14
cc_annot$Clusters <- factor(cc_annot$Clusters_new)
write.table(cc_annot[,c(1,2,3,10)], file = "cc.xls", quote = F, sep = "\t", row.names = F)
}

col = c("#FF7777", "#EDE447", "#42B540", "#0099B4", "#00468B", "grey")
dat = read.table("data.xls", sep = "\t", header = T)
dat$Samples_Clusters = factor(dat$Samples_Clusters, levels = c("WT-1_0", "WT-2_0", "Aa-1_0", "Aa-2_0", "KO-1_0", "KO-2_0", "WT-1_1", "WT-2_1", "Aa-1_1", "Aa-2_1", "KO-1_1", "KO-2_1", "WT-1_2", "WT-2_2", "Aa-1_2", "Aa-2_2", "KO-1_2", "KO-2_2", "WT-1_3", "WT-2_3", "Aa-1_3", "Aa-2_3", "KO-1_3", "KO-2_3", "WT-1_4", "WT-2_4", "Aa-1_4", "Aa-2_4", "KO-1_4", "KO-2_4", "WT-1_5", "WT-2_5", "Aa-1_5", "Aa-2_5", "KO-1_5", "KO-2_5", "WT-1_6", "WT-2_6", "Aa-1_6", "Aa-2_6", "KO-1_6", "KO-2_6", "WT-1_7", "WT-2_7", "Aa-1_7", "Aa-2_7", "KO-1_7", "KO-2_7", "WT-1_8", "WT-2_8", "Aa-1_8", "Aa-2_8", "KO-1_8", "KO-2_8", "WT-1_9", "WT-2_9", "Aa-1_9", "Aa-2_9", "KO-1_9", "KO-2_9", "WT-1_10", "WT-2_10", "Aa-1_10", "Aa-2_10", "KO-1_10", "KO-2_10", "WT-1_11", "WT-2_11", "Aa-1_11", "Aa-2_11", "KO-1_11", "KO-2_11", "WT-1_12", "WT-2_12", "Aa-1_12", "Aa-2_12", "KO-1_12", "KO-2_12", "WT-1_13", "WT-2_13", "Aa-1_13", "Aa-2_13", "KO-1_13", "KO-2_13", "WT-1_14", "WT-2_14", "Aa-1_14", "Aa-2_14", "KO-1_14", "KO-2_14"#, "WT-1_15", "WT-2_15", "Aa-1_15", "Aa-2_15", "KO-1_15", "KO-2_15", "WT-1_16", "WT-2_16", "Aa-1_16", "Aa-2_16", "KO-1_16", "KO-2_16", "WT-1_17", "WT-2_17", "Aa-1_17", "Aa-2_17", "KO-1_17", "KO-2_17", "WT-1_18", "WT-2_18", "Aa-1_18", "Aa-2_18", "KO-1_18", "KO-2_18", "WT-1_19", "WT-2_19", "Aa-1_19", "Aa-2_19", "KO-1_19", "KO-2_19", "WT-1_20", "WT-2_20", "Aa-1_20", "Aa-2_20", "KO-1_20", "KO-2_20", "WT-1_21", "WT-2_21", "Aa-1_21", "Aa-2_21", "KO-1_21", "KO-2_21", "WT-1_22", "WT-2_22", "Aa-1_22", "Aa-2_22", "KO-1_22", "KO-2_22", "WT-1_23", "WT-2_23", "Aa-1_23", "Aa-2_23", "KO-1_23", "KO-2_23"
))
#dat$Phase = factor(dat$Phase, levels = c("non-cycling", "M/G1", "M", "G2/M", "S", "G1/S"))
dat$Phase = factor(dat$Phase, levels = c("G1/S", "S", "G2/M", "M", "M/G1", "non-cycling"))
pdf("cellcycle.pdf", width = 30)
ggplot(dat, aes(x = Samples_Clusters, fill = Phase)) + geom_bar(stat = "count", position = "fill", width = 0.8) + labs(y = "Number of Cells") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5)) + scale_fill_manual(values = col)
