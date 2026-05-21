library("posterior")

setwd("/Users/nan/Desktop/科研/liver/liver.human.dna/write paper/Figure/Figure 5/香农指数")
data = read.table("all.align.fasta.seq",sep = "\t")

step = 20
column_names <- c("num", "shannong","average")
out_df <- data.frame(matrix(nrow = dim(data)[2], ncol = length(column_names)))
names(out_df) <- column_names
for (i in 1:dim(data)[2]) {
	temp = data[, i]
	temp = temp[temp %in% c("A","T","G","C")]
	out_df[i, 1] = length(temp)
	out_df[i, 2] = entropy(factor(temp,levels=c("A","T","G","C")))
}
for (l in 1: (dim(data)[2] - step)) {
	av = mean(out_df$shannong[l : (l + step)])
	out_df[l, 3] = av
}

#library(grDevices)
#my_red = adjustcolor("#845EC2", alpha.f = 1)
plot(out_df$shannong,type = "h",col="#DCDCDC", ylim=c(-0.7,1),axes = FALSE)
axis(side = 2, at = seq(-0.7, 1, 0.1))
axis(side = 1)
points(out_df$average,type = "l",col="#696969",lwd=0.5)
#library(grDevices)
#my_A = adjustcolor("#C2B95E", alpha.f = 0.8)
#points(A.out_df$average,type = "h",col=my_A,lwd=0.5)
#my_B = adjustcolor("#00B1FF", alpha.f = 0.8)
#points(B.out_df$average,type = "h",col="#00B1FF",lwd=0.5)


segments(0,-0.2,2045,-0.2,col = "#010604")
#polygon(c(0,0,2045,2045),c(-0.19,-0.21,-0.21,-0.19),col = "grey")
polygon(c(155,155,1116,1136,1116),c(-0.17,-0.23,-0.23,-0.2,-0.17),col = "#E1972E",lwd=0.1)
polygon(c(1402,1382,1402,2035,2035),c(-0.17,-0.2,-0.23,-0.23,-0.17),col = "#4FA3B8",lwd=0.1)

primer_col = "#F8600D"
t_col = "#2589D3"
u_col = "#4AB47D"
#tongyong
primer1 = c(574, 594)
loc = 0
polygon(c(primer1[1],primer1[1],primer1[2],primer1[2]),c(-0.3-loc,-0.33-loc,-0.33-loc,-0.3-loc),col = primer_col,lwd=0.1)
primer2 = c(1306, 1327)
polygon(c(primer2[1],primer2[1],primer2[2],primer2[2]),c(-0.3-loc,-0.33-loc,-0.33-loc,-0.3-loc),col = primer_col,lwd=0.1)
polygon(c(primer1[2],primer1[2],primer2[1],primer2[1]),c(-0.3-loc,-0.33-loc,-0.33-loc,-0.3-loc),col = t_col,lwd=0.1)


#A:
loc = 0.06
primer1 = c(697, 716)
polygon(c(primer1[1],primer1[1],primer1[2],primer1[2]),c(-0.3-loc,-0.33-loc,-0.33-loc,-0.3-loc),col = primer_col,lwd=0.1)
primer2 = c(1575, 1594)
polygon(c(primer2[1],primer2[1],primer2[2],primer2[2]),c(-0.3-loc,-0.33-loc,-0.33-loc,-0.3-loc),col = primer_col,lwd=0.1)
polygon(c(primer1[2],primer1[2],primer2[1],primer2[1]),c(-0.3-loc,-0.33-loc,-0.33-loc,-0.3-loc),col = u_col,lwd=0.1)

#B
loc = 0.12
primer1 = c(134, 153)
polygon(c(primer1[1],primer1[1],primer1[2],primer1[2]),c(-0.3-loc,-0.33-loc,-0.33-loc,-0.3-loc),col = primer_col,lwd=0.1)
primer2 = c(1875, 1894)
polygon(c(primer2[1],primer2[1],primer2[2],primer2[2]),c(-0.3-loc,-0.33-loc,-0.33-loc,-0.3-loc),col = primer_col,lwd=0.1)
polygon(c(primer1[2],primer1[2],primer2[1],primer2[1]),c(-0.3-loc,-0.33-loc,-0.33-loc,-0.3-loc),col = u_col,lwd=0.1)

#C
loc = 0.18
primer1 = c(1668, 1687)
polygon(c(primer1[1],primer1[1],primer1[2],primer1[2]),c(-0.3-loc,-0.33-loc,-0.33-loc,-0.3-loc),col = primer_col,lwd=0.1)
primer2 = c(1964, 1983)
polygon(c(primer2[1],primer2[1],primer2[2],primer2[2]),c(-0.3-loc,-0.33-loc,-0.33-loc,-0.3-loc),col = primer_col,lwd=0.1)
polygon(c(primer1[2],primer1[2],primer2[1],primer2[1]),c(-0.3-loc,-0.33-loc,-0.33-loc,-0.3-loc),col = u_col,lwd=0.1)

#D
loc = 0.24
primer1 = c(904, 924)
polygon(c(primer1[1],primer1[1],primer1[2],primer1[2]),c(-0.3-loc,-0.33-loc,-0.33-loc,-0.3-loc),col = primer_col,lwd=0.1)
primer2 = c(1849, 1874)
polygon(c(primer2[1],primer2[1],primer2[2],primer2[2]),c(-0.3-loc,-0.33-loc,-0.33-loc,-0.3-loc),col = primer_col,lwd=0.1)
polygon(c(primer1[2],primer1[2],primer2[1],primer2[1]),c(-0.3-loc,-0.33-loc,-0.33-loc,-0.3-loc),col = u_col,lwd=0.1)

