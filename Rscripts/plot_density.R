library(Cairo)
library(Seurat)
library(dplyr)
library(magrittr)
library(patchwork)
library(monocle)
library(reshape2)
library(ggplot2)
library(DESeq2)
library(BiocParallel)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(magick)
library(readr)
library(Matrix)
library(Scillus)
library(clustree)
library(Nebulosa)
set.seed(123)

#======================================================
proj_dir <-("/research/labs/bme/weiz/m238739/Paired-Tag/data/data16/mouse/")
fig_dir <- (paste0(proj_dir,"figures/Manuscript_V1/Fig"))
out_dir <- (paste0(proj_dir,"figures/Manuscript_V1/Fig/fig1/data"))
dat_dir <- (paste0(proj_dir,"bam/RNA"))

ptag_list_integrated <- readRDS(file = paste0(out_dir, "/","ptag_list_integrated.rds"))
markers <- readRDS(file = paste0(out_dir, "/","ptag_list_integrated.markers.rds"))
paired = c("Beta-hi"="#34A047","Beta-low"="#1E78B4","Alpha"="#6A3E98", "Delta"="#B15928", "PP"="#F59899","Macrophage"="#E11E26",
           "Mki67-beta"="#FCBF6E","Endothelial"="#F47E1F","Duct"="#CAB2D6","Fibroblasts"="#FAF39B")


pdf(paste0(fig_dir, "/","fig1/fig4s.sex.UMAP.1.pdf"), height = 7, width = 10)
DimPlot(ptag_list_integrated, reduction = "umap", split.by = "Sex", label = TRUE, pt.size = 0.5,cols=paired)+NoLegend()+ NoAxes()
dev.off()
#------------------------------------------------------------------------extract barcode
class(ptag_list_integrated@meta.data$seurat_clusters)
sampledata <- c("B5","B6","B7","B8","B9","B10","B11","B12","B14","B15")
for (x in sampledata){
ptag_list_integrated.x <- subset(x = ptag_list_integrated, subset = Sample == x)
cellbarcode<- as.data.frame(ptag_list_integrated.x@meta.data[,c("CellType")])
rownames(cellbarcode)<- rownames(ptag_list_integrated.x@meta.data)
rownames(cellbarcode) <- gsub("\\_.*", replacement = "", rownames(cellbarcode))
rownames(cellbarcode) <- gsub("-", replacement = ":", rownames(cellbarcode))

colnames(cellbarcode)<- c("cluster")
cellbarcode<- tibble::rownames_to_column(cellbarcode, "Barcode")
head(cellbarcode)

write.csv(cellbarcode,paste0(out_dir, "/", x,"_rnaseq_clusters.csv"), row.names = FALSE)}
#-------------------------------------------------------------------------
H3K4me1_ptag_list_integrated <-subset(ptag_list_integrated, subset = Histone == 'H3K4me1')
class(H3K4me1_ptag_list_integrated@meta.data$seurat_clusters)
dim(H3K4me1_ptag_list_integrated@meta.data)
sampledata <- c("0","8","16")
for (x in sampledata){
  H3K4me1_ptag_list_integrated.x <- subset(x = H3K4me1_ptag_list_integrated, subset = Week == x)

  cellbarcode<- as.data.frame(H3K4me1_ptag_list_integrated.x@meta.data[,c("Week")])
  rownames(cellbarcode)<- rownames(H3K4me1_ptag_list_integrated.x@meta.data)

  colnames(cellbarcode)<- c("Week")
  cellbarcode<- tibble::rownames_to_column(cellbarcode, "Barcode")
  head(cellbarcode)

  cellbarcode$Barcode <- gsub("\\_.*", replacement = "", cellbarcode$Barcode)
  cellbarcode$Barcode <- gsub("-", replacement = ":", cellbarcode$Barcode)
write.csv(cellbarcode,paste0(out_dir, "/", x,"Week_rnaseq_H3K4me1.csv"), row.names = FALSE)}

#--------------------------------------------------------------------------------------------------------
H3K27ac_ptag_list_integrated <-subset(ptag_list_integrated, subset = Histone == 'H3K27ac')
class(H3K27ac_ptag_list_integrated@meta.data$seurat_clusters)
dim(H3K27ac_ptag_list_integrated@meta.data)
sampledata <- c("0","8","16")
for (x in sampledata){
  H3K27ac_ptag_list_integrated.x <- subset(x = H3K27ac_ptag_list_integrated, subset = Week == x)

  cellbarcode<- as.data.frame(H3K27ac_ptag_list_integrated.x@meta.data[,c("Week")])

  rownames(cellbarcode)<- rownames(H3K27ac_ptag_list_integrated.x@meta.data)

  colnames(cellbarcode)<- c("Week")


  cellbarcode<- tibble::rownames_to_column(cellbarcode, "Barcode")
  head(cellbarcode)

  cellbarcode$Barcode <- gsub("\\_.*", replacement = "", cellbarcode$Barcode)
  cellbarcode$Barcode <- gsub("-", replacement = ":", cellbarcode$Barcode)

  write.csv(cellbarcode,paste0(out_dir, "/", x,"Week_rnaseq_H3K27ac.csv"), row.names = FALSE)}
#-----------------------------------------------------------------------------------
sampledata <- c("B5","B6","B7","B8","B9","B10","B11","B12","B14","B15")
for (x in sampledata){
  ptag_list_integrated.x <- subset(x = ptag_list_integrated, subset = Sample == x)

  cellbarcode<- as.data.frame(ptag_list_integrated.x@meta.data[,c("Sample")])
  rownames(cellbarcode)<- rownames(ptag_list_integrated.x@meta.data)
  rownames(cellbarcode) <- gsub("\\_.*", replacement = "", rownames(cellbarcode))
  rownames(cellbarcode) <- gsub("-", replacement = ":", rownames(cellbarcode))

  colnames(cellbarcode)<- c("Sample")
  cellbarcode<- tibble::rownames_to_column(cellbarcode, "Barcode")
  head(cellbarcode)

  write.csv(cellbarcode,paste0(out_dir, "/", x,"_rnaseq_Sample.csv"), row.names = FALSE)
  write.table(cellbarcode[,1],paste0(out_dir, "/", x,"_rnaseq_Sample2.csv"), row.names = FALSE,col.names =FALSE,quote=F)
  }
#----------------------------------------------------------------------------
H3K27ac.RNA.metadata <-ptag_list_integrated@meta.data[ptag_list_integrated@meta.data$Histone=='H3K27ac',]
rownames(H3K27ac.RNA.metadata) <- with(H3K27ac.RNA.metadata, paste(orig.ident, rownames(H3K27ac.RNA.metadata), sep = "_"))
rownames(H3K27ac.RNA.metadata) <- sub("(_[^_]+)_.*", "\\1", rownames(H3K27ac.RNA.metadata))
rownames(H3K27ac.RNA.metadata) <- sub("_", "#", rownames(H3K27ac.RNA.metadata))
H3K27ac.beta.ArchR<- H3K27ac.RNA.metadata %>% dplyr::filter(CellType %in% c("Mki67-beta","Beta-hi","Beta-low"))
H3K27ac.beta.ArchR<-H3K27ac.beta.ArchR["CellType"]
class(H3K27ac.beta.ArchR)
write.csv(H3K27ac.beta.ArchR, paste0(out_dir, "/", "H3K27ac.beta.ArchR_clusters.csv"))
#---------------------------------------------------------------------------
H3K4me1.RNA.metadata <-ptag_list_integrated@meta.data[ptag_list_integrated@meta.data$Histone=='H3K4me1',]
rownames(H3K4me1.RNA.metadata) <- with(H3K4me1.RNA.metadata, paste(orig.ident, rownames(H3K4me1.RNA.metadata), sep = "_"))
rownames(H3K4me1.RNA.metadata) <- sub("(_[^_]+)_.*", "\\1", rownames(H3K4me1.RNA.metadata))
rownames(H3K4me1.RNA.metadata) <- sub("_", "#", rownames(H3K4me1.RNA.metadata))
H3K4me1.beta.ArchR<- H3K4me1.RNA.metadata %>% dplyr::filter(CellType %in% c("Mki67-beta","Beta-hi","Beta-low"))
H3K4me1.beta.ArchR<-H3K4me1.beta.ArchR["CellType"]
class(H3K4me1.beta.ArchR)
write.csv(H3K4me1.beta.ArchR, paste0(out_dir, "/", "H3K4me1.beta.ArchR_clusters.csv"))
#--------------------------------------------------------------------------------
seurat_obj <-ptag_list_integrated
seurat_obj$barcode <-paste0(seurat_obj@meta.data$orig.ident,"_", row.names(seurat_obj@meta.data))
seurat_obj$barcode <- sub("(_[^_]+)_.*", "\\1", seurat_obj$barcode)
seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
write.csv(seurat_obj@meta.data, file=paste0(out_dir,"/",'ptag_list_integrated_metadata.csv'), quote=F, row.names=F)
#-----------------------------------------------------------------------------
seurat_obj@meta.data$CellType2<-seurat_obj@meta.data$CellType
seurat_obj@meta.data$CellType2<-sub("Beta-hi", "Beta", seurat_obj@meta.data$CellType2)
seurat_obj@meta.data$CellType2<-sub("Beta-low", "Beta", seurat_obj@meta.data$CellType2)
write.csv(seurat_obj@meta.data, file=paste0(out_dir,"/",'ptag_list_integrated_metadata2.csv'), quote=F, row.names=F)
# write expression counts matrix
counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0(out_dir,"/", 'ptag_list_integrated_counts.mtx'))
#
write.csv(seurat_obj@reductions$pca@cell.embeddings, file=paste0(out_dir,"/",'ptag_list_integrated_pca.csv'), quote=F, row.names=F)
# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file=paste0(out_dir,"/",'ptag_list_integrated_gene_names.csv'),
  quote=F,row.names=F,col.names=F
)

#--------------------------------------------------------------------------------------------------------------
DefaultAssay(ptag_list_integrated) <- "SCT"
p4 <- plot_density(ptag_list_integrated, c("Ins1", 'Mki67'), joint = TRUE)
p4 + plot_layout(ncol = 1)
DimPlot(ptag_list_integrated, label = TRUE, repel = TRUE,  pt.size = 0.5, reduction = "umap", label.size = 8)

greyRed = c("1"="lightgrey", "2"="red")
pdf(paste0(fig_dir, "/","fig6/umap_3.pdf"))
plot_density(ptag_list_integrated, "Ngf",shape = 16,size = 2, pal = "viridis",adjust = 1,)+ scale_colour_gradientn(colours =brewer.pal(n = 11, name = "Reds"))
plot_density(ptag_list_integrated, "Manf",shape = 16,size = 2, pal = "viridis",adjust = 1,)+ scale_colour_gradientn(colours =greyRed)
dev.off()





















