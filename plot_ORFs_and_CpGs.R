#Description : Visualize genome ORF architecture and CpG island annotations.
#               Each molecule (virus genome) is plotted as a separate panel.
#               Genes are shown as directional arrows; CpG islands as rectangles.
#               R version 4.3.3
# Input       : A tab-separated text file (.txt / .tsv) with header.
#               Expected columns (in order):
#                 molecule  - genome identifier (e.g. accession + virus name)
#                 source    - annotation source (Geneious, cpgplot, etc.)
#                 gene      - gene name or feature type (cap, rep, CpG, etc.)
#                 start     - feature start position (bp)
#                 end       - feature end position (bp)
#                 direction - strand: 1 = forward, -1 = reverse, 0 = unknown
#                 genome    - total genome length (bp)
#
# Output      : PDF file named "<input_file>.pdf"
#
# Usage       : Rscript plot_genome.R "ORFs and CpG structure data"
args = commandArgs(T)
in_file = args[1]

width = 0.18
data = read.table(in_file,header=T,sep='\t',stringsAsFactor=F)

molecules = data$molecule[!duplicated(data$molecule)]
molecule_num = length(data$molecule[!duplicated(data$molecule)])
max_length = max(data$end)
#plot(c(1,1),type='n',xlim=c(1,max_length),yaxt="n",xaxt="n",xlab="",ylab="",bty="n")
plot_genome=function(){
	par(mar=c(0,2,0,2))
	layout(matrix(1:molecule_num))
	for (molecule in molecules){
		plot(c(1,1),type='n',xlim=c(1,max_length),yaxt="n",xaxt="n",xlab="",ylab="",bty="n")
		temp_data = data[data$molecule == molecule,]
		genome_high = 1
		genome_name = temp_data$molecule[1]
		genome_length = temp_data$genome[1]
		orfs = temp_data[temp_data$molecule==temp_data$molecule[!duplicated(temp_data$molecule)],c(3,4,5,6)]
		rect(1,genome_high,genome_length,genome_high)
#		plot(c(1,genome_length),c(genome_high,genome_high),type='l',,col="black",xlim=c(1,max_length),lwd=3,yaxt="n",xaxt="n",xlab="",ylab="",bty="n")

#		text(1,genome_high-0.1,1,cex=0.75);text(genome_length,genome_high-0.1,genome_length,cex=0.75)
		orf_plot = orfs
#		orf_plot = orfs[!duplicated(orfs$ORF),]
		num = 0
		for (i in 1:dim(orf_plot)[1]){
			orf_name = orf_plot[i,][1]
			orf_start = orf_plot[i,][2]
			orf_end = orf_plot[i,][3]
			if (orf_plot[i,][4]==1) {
				polygon(c(orf_start,orf_end-15,orf_end,orf_end-15,orf_start),c(genome_high-width-num,genome_high-width-num,genome_high-num,genome_high+width-num,genome_high+width-num)-width,col="#F59319",lwd=0.2)
			}
			if (orf_plot[i,][4]==-1) {
				polygon(c(orf_start,orf_start+15,orf_end,orf_end,orf_start+15),c(genome_high-num,genome_high+width-num,genome_high+width-num,genome_high-width-num,genome_high-width-num)-width,col="#5AAED9",lwd=0.2)
			}
			if (orf_plot[i,][4]==0) {
				polygon(c(orf_start,orf_end,orf_end,orf_start),c(genome_high+width-num,genome_high+width-num,genome_high-width-num,genome_high-width-num) + +width,col="#8DF830",lwd=0.2)
			}
#			text((orf_start+orf_end)/2,genome_high-0.08-num,orf_name,cex=0.75)
			num = num + 0
#		text((1+genome_length)/2,genome_high,genome_name,cex=0.8)
		}
		genome_high = genome_high + 1
#		gene_plot = temp_data[!is.na(temp_data$GENES),]
#		for (i in 1:dim(gene_plot)[1]){
#			orf_start = gene_plot[i,][5]
#			gene_name = gene_plot[i,][7]
#			gene_start = orf_start + (gene_plot[i,][8]-1)*3
#			gene_end = orf_start + (gene_plot[i,][9]-1)*3
#			polygon(c(gene_start,gene_end-40,gene_end,gene_end-40,gene_start),c(genome_high-0.02,genome_high-0.02,genome_high,genome_high+0.02,genome_high+0.02),col="blue")
#			text((gene_start+gene_end)/2,genome_high+0.08,gene_name,cex=0.75)
#		}
	}
}

pdf(paste(in_file,'.pdf',sep=''),width=7,height=molecule_num/3)
plot_genome()
dev.off()
#emf(paste(in_file,'.emf',sep=''),width=7,height=genome_num)
#plot_genome()
#dev.off()
