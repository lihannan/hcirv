library(pheatmap)
library(RColorBrewer)
library(gridExtra)
files = list.files(pattern="\\.csv$")

mat <- matrix(0, nrow = 0, ncol = 5)
df <- data.frame(mat, check.names = FALSE)
colnames(df) = c("A", "B", "C", "D","max")

plot_list = list()
for (file in files){
	data = read.csv(file,row.name = 1)
	name = strsplit(file,"\\.")[[1]][1]
	p = pheatmap(data, display_numbers = TRUE,legend_breaks=c(0.5,0.6,0.7,0.8,0.9,1), legend_labels=c(0.5,0.6,0.7,0.8,0.9,1), breaks = seq(0.5, 1, by = 0.01), color = c(colorRampPalette(c("blue", "red"))(50)), cutree_rows = 4, na_col = "#DDDDDD", treeheight_col = 0, number_color = "white",fontsize_number=16, cellwidth = 60, cellheight = 60, main = name, fontsize = 15, legend = F)
	plot_list[[file]] = p[[4]]
	df[name, "A"] = data[name, "A"]
	df[name, "B"] = data[name, "B"]
	df[name, "C"] = data[name, "C"]
	df[name, "D"] = data[name, "D"]
	data <- data[-which(names(data) == gsub("\\-",".",name)), -which(names(data) == gsub("\\-",".",name))]
	lower_tri_elements <- data[lower.tri(data)]
	df[name, "max"] = max(lower_tri_elements)
}

pdf("test.pdf",width = 40, height = 200)
do.call("grid.arrange", c(plot_list,  ncol = 5))
dev.off()
write.csv(df,"test.csv", row.names=T)