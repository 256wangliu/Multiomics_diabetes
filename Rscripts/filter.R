setwd("/research/labs/bme/weiz/m238739/Paired-Tag/data/data16/mouse/bam")
dir()

library(Matrix)
plot_BC<-function(dna, rna, c1=1000, c2=1000, title=title){
layout(matrix(c(1,2,1,3),2,2,byrow=TRUE))
m<-merge(dna, rna, by="BC")
plot(m[,2:3],pch=19, cex=0.1, log="xy", xlab="DNA reads", ylab="RNA reads", col="grey", xlim=c(1,100000), ylim=c(1,100000), main=title)
lines(c(c1,c1),c(1,5000000), lty=2);lines(c(1,5000000),c(c2,c2), lty=2)
pf<-m[as.numeric(m[,2]) >= c1 & as.numeric(m[,3]) >= c2,]
points(x=pf[,2], y=pf[,3], pch=19, cex=0.1, col="firebrick")
legend("bottomright", legend=c( paste("DNA cutoff:", c1, sep=" "), paste("RNA cutoff:", c2, sep = " "), paste("# PF:", dim(pf)[1], sep=" ")), bty="n")
data<-na.omit(pf); val.cell<-data[,1]
write.table(val.cell, file=paste(title,"filt.xls", sep="_"), col.names=F, row.names=F, sep="\t", quote=F)
}

dna.path="/research/labs/bme/weiz/m238739/Paired-Tag/data/data16/mouse/bam/B14_D6_merge_mm10_filtered_10_mtx2"
rna.path="/research/labs/bme/weiz/m238739/Paired-Tag/data/data16/mouse/bam/RNA/B14_R6_sorted_rmdup_mtx2"

dna.mtx<-readMM(paste(dna.path, "/matrix.mtx", sep=""))
dna.bc<-read.table(paste(dna.path, "/barcodes.tsv", sep=""), head=F)
dna<-cbind(dna.bc[,1], colSums(dna.mtx))
rna.mtx<-readMM(paste(rna.path, "/matrix.mtx", sep=""))

rna.bc<-read.table(paste(rna.path, "/barcodes.tsv", sep=""), head=F)
rna<-cbind(rna.bc[,1], colSums(rna.mtx))

colnames(dna)<-c("BC", "n_reads_DNA")
colnames(rna)<-c("BC", "n_reads_RNA")

plot_BC(dna, rna,
c1=500,
c2=200,
title="B14_D6_merge_DNA_RNA")






















