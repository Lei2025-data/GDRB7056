library(pheatmap)
args = commandArgs(T)
data = read.table (file = args[1], header = T, row.names = 1, sep = "\t", check.names = F)
data = data[,2:ncol(data)]
clean_data = data[which(rowSums(data)>0),]

main = "Significant Differential Gene Heatmap"
color = colorRampPalette(c("steelblue", "white", "red"))(256)
cluster_cols = T
cluster_rows = T
show_rownames = T

gene_num = nrow(clean_data)
if (gene_num <= 200){
	fontsize_row = 4
	cellheight = 4
}else if (gene_num > 200 && gene_num < 500){
	fontsize_row = 3
	cellheight = 3
}else if (gene_num > 500 && gene_num < 2000){
	fontsize_row = 2
	cellheight = 2
}else{
	show_rownames = F
	cellheight = 4000/gene_num
}

if (length(args) == 3){
	group_info = read.table (file = args[3], row.names = 2, check.names = F, sep = "\t")
	colnames(group_info) = c("Group")
	cluster_cols = F
	annotation_col = as.data.frame(group_info)

	scale = "row"
	filename = paste(args[2], ".heatmap.png", sep = "")
	pheatmap(clean_data, filename = filename, main = main, color = color, 
			cluster_cols = cluster_cols, cluster_rows = cluster_rows,
			show_rownames = show_rownames, fontsize_row = fontsize_row,
			cellheight = cellheight, annotation_col = annotation_col,
			scale = scale,
			)
	filename = paste(args[2], ".heatmap.pdf", sep = "")
	pheatmap(clean_data, filename = filename, main = main, color = color,
			cluster_cols = cluster_cols, cluster_rows = cluster_rows,
			show_rownames = show_rownames, fontsize_row = fontsize_row,
			cellheight = cellheight, annotation_col = annotation_col,
			scale = scale,
	)
}else{
	scale = "row"
	filename = paste(args[2], ".heatmap.png", sep = "")
	p = pheatmap(clean_data, filename = filename, main = main, color = color, 
			cluster_cols = cluster_cols, cluster_rows = cluster_rows,
			show_rownames = show_rownames, fontsize_row = fontsize_row,
			cellheight = cellheight, scale = scale,
			)
	print (str(p$tree_col))
	print (paste("Cluster",p$tree_col$order-1))
	filename = paste(args[2], ".heatmap.pdf", sep = "")
	pheatmap(clean_data, filename = filename, main = main, color = color, 
			cluster_cols = cluster_cols, cluster_rows = cluster_rows,
			show_rownames = show_rownames, fontsize_row = fontsize_row,
			cellheight = cellheight, scale = scale,
			)
}

