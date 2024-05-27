library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(ChIPseeker)
library(clusterProfiler)
library(ggupset)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPpeakAnno)
library(Matrix)
library(GenomicRanges)
set.seed(1)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
dir()
proj_dir <-("/research/labs/bme/weiz/m238739/Paired-Tag/data/data16/mouse/")
out_dir <- (paste0(proj_dir,"figures/Manuscript_V1/Fig/figs3/data"))
setwd("/research/labs/bme/weiz/m238739/Paired-Tag/data/data16/mouse/bam/snap/splitbam/")
fig_dir <-("/research/labs/bme/weiz/m238739/Paired-Tag/data/data16/mouse/review1/fig/")

H3K27Ac <- readRDS(file = paste0(out_dir, "/","fig3B.h3k27Ac_integrated.rds"))
DefaultAssay(H3K27Ac) <- "peaks"
H3K27Ac[['peaks']]@counts
class(H3K27Ac[['peaks']]@counts)
count_matrix <- H3K27Ac[['peaks']]@counts

# Calculate row sums
row_sums <- rowSums(count_matrix)
# Convert row names into a GRanges object
row_names <- rownames(count_matrix)
seqnames <- sapply(row_names, function(x) strsplit(x, "-")[[1]][1])
starts <- sapply(row_names, function(x) as.numeric(strsplit(x, "-")[[1]][2]))
ends <- sapply(row_names, function(x) as.numeric(strsplit(x, "-")[[1]][3]))
# Create GRanges object
gr <- GRanges(seqnames = seqnames,
              ranges = IRanges(start = starts, end = ends),
              strand = rep("*", length(row_sums)))
# Add metadata columns
mcols(gr) <- DataFrame(V4 = rep("#", length(row_sums)), V5 = row_sums)
K27_peaks <-gr
K27_peaks <- keepStandardChromosomes(K27_peaks, pruning.mode = "coarse")
K27_peakAnno <- annotatePeak(K27_peaks, tssRegion=c(-1000, 1000),TxDb=txdb, annoDb="org.Mm.eg.db")
#----------------------------------------------------------------------------------------------------------------------
out_dir2 <- (paste0(proj_dir,"figures/Manuscript_V1/Fig/fig3/data"))
H3K4me1 <- readRDS(file = paste0(out_dir2, "/","fig3a.h3k4me1_integrated.rds"))
DefaultAssay(H3K4me1) <- "peaks"
H3K4me1[['peaks']]@counts
H3K4me1
count_matrix <- H3K4me1[['peaks']]@counts
# Calculate row sums
row_sums <- rowSums(count_matrix)
# Convert row names into a GRanges object
row_names <- rownames(count_matrix)
seqnames <- sapply(row_names, function(x) strsplit(x, "-")[[1]][1])
starts <- sapply(row_names, function(x) as.numeric(strsplit(x, "-")[[1]][2]))
ends <- sapply(row_names, function(x) as.numeric(strsplit(x, "-")[[1]][3]))
# Create GRanges object
gr <- GRanges(seqnames = seqnames,
              ranges = IRanges(start = starts, end = ends),
              strand = rep("*", length(row_sums)))
# Add metadata columns
mcols(gr) <- DataFrame(V4 = rep("#", length(row_sums)), V5 = row_sums)
K4_peaks <-gr
#K4_peaks <-granges(H3K4me1)
K4_peak <- keepStandardChromosomes(K4_peaks, pruning.mode = "coarse")
K4_peakAnno <- annotatePeak(K4_peak, tssRegion=c(-1000, 1000),TxDb=txdb, annoDb="org.Mm.eg.db")

peaks <- GRangesList(H3K4me1=K4_peaks,
                     H3K27ac=K27_peaks)
genomicElementDistribution(peaks,
                           TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                           promoterRegion=c(upstream=200, downstream=200),
                           geneDownstream=c(upstream=0, downstream=2000))

pdf(paste0(fig_dir,"H3K4me1_H3K27ac_peaks_001.pdf"), width = 10, height = 5, useDingbats = F)
p <- genomicElementDistribution(peaks,
                           TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                           promoterRegion=c(upstream=200, downstream=200),
                           geneDownstream=c(upstream=0, downstream=2000))
p
dev.off()


