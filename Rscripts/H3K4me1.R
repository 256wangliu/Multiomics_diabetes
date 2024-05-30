# load package
library(Seurat)
library(harmony)
library(clustree)
library(ggpubr)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(vegan)
library(future)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(flipPlots)
library(cowplot)
library(webshot)

plan(multicore)
options(future.globals.maxSize = 100000 * 1024^2)
set.seed(101)
#--------------------------------------------------------------------------------------------------
proj_dir <-("/research/labs/bme/weiz/m238739/Paired-Tag/data/data16/mouse/")
fig_dir <- (paste0(proj_dir,"figures/Manuscript_V1/Fig"))
out_dir <- (paste0(proj_dir,"figures/Manuscript_V1/Fig/fig3/data"))
out_dir3 <- (paste0(proj_dir,"figures/Manuscript_V1/Fig/fig1/data"))

integrated <- readRDS(file = paste0(out_dir, "/","fig3a.h3k4me1_integrated.rds"))
DefaultAssay(integrated) <- "SCT"
paired = c("Beta-hi"="#34A047","Beta-low"="#1E78B4","Alpha"="#6A3E98", "Delta"="#B15928", "PP"="#F59899","Macrophage"="#E11E26",
           "Mki67-beta"="#FCBF6E","Endothelial"="#F47E1F","Duct"="#CAB2D6","Fibroblasts"="#FAF39B")

pdf(paste0(fig_dir, "/","fig3/Fig3A.H3K4me1.UMAP.1.pdf"), height = 7, width = 7)
DimPlot(integrated, reduction = "integrated.wnn.umap", label = T, pt.size = 0.5,cols=paired )+NoLegend()+ NoAxes()
dev.off()


#-------------------------------------------------------------------------------------------------
h3k4me1_integrated <- readRDS(file = paste0(out_dir, "/","fig3a.h3k4me1_integrated.rds"))
H3K4ME1_markers <- readRDS(file = paste0(out_dir, "/","H3K4ME1_integrated.markers.rds"))
ptag_list_integrated <- readRDS(file = paste0(out_dir3, "/","ptag_list_integrated.rds"))

H3K4me1.RNA.metadata <-ptag_list_integrated@meta.data[ptag_list_integrated@meta.data$Histone=='H3K4me1',]
rownames(H3K4me1.RNA.metadata) <- with(H3K4me1.RNA.metadata, paste(orig.ident, rownames(H3K4me1.RNA.metadata), sep = "_"))
rownames(H3K4me1.RNA.metadata) <- sub("(_[^_]+)_.*", "\\1", rownames(H3K4me1.RNA.metadata))
subset_h3k4me1_integrated <- subset(h3k4me1_integrated, cells = rownames(H3K4me1.RNA.metadata))
subset_h3k4me1_integrated

classTab <- data.frame(row.names = rownames(H3K4me1.RNA.metadata),
                       RNA = H3K4me1.RNA.metadata$CellType,
                       combine = subset_h3k4me1_integrated@meta.data$CellType)

RNAClass <- classTab$RNA
combineClass <- classTab$combine
chisq.test(table(RNAClass, combineClass))

transMat <- as.data.frame(table(RNAClass,combineClass), stringsAsFactors = F)
transMat$RNAClass <- factor(transMat$RNAClass, levels = names(table(H3K4me1.RNA.metadata$CellType)))
transMat$combineClass <- factor(transMat$combineClass, levels = names(table(subset_h3k4me1_integrated@meta.data$CellType)))

options(repr.plot.height = 10, repr.plot.width = 2)
paired4 = c("#34A047","#1E78B4","#6A3E98", "#B15928","#F59899","#E11E26","#FCBF6E","#F47E1F","#CAB2D6","#FAF39B")

widget = SankeyDiagram(transMat[, -3],
                       max.categories = 40,
                       link.color = "Source",
                       label.show.varname = FALSE,
                       font.size = 12,
                       weights = transMat$Freq,node.padding = 30,
                       node.width = 100,
                       colors =paired4 ,
                       node.position.automatic = FALSE
)

htmlwidgets::saveWidget(widget, paste0(fig_dir, "/",'fig3/sankey_plot_extended3.html'))
webshot::webshot(paste0(fig_dir, "/",'fig3/sankey_plot_extended3.html'), file = paste0(fig_dir, "/",'fig3/sankey_plot_extended3.pdf'))
#saveRDS(transMat,paste0(out_dir, "/","h3k4me1.rna.peaks.transMat3.rds"))
#------------------------------------------------------------------------------------------------------------------------
source(file ="plot_dependencies.R")
h3k4me1_integrated@meta.data
h3k4me1_integrated$CellType2 <-h3k4me1_integrated$CellType
h3k4me1_integrated$CellType2 <- gsub("Beta-hi", replacement = "Beta", h3k4me1_integrated$CellType2)
h3k4me1_integrated$CellType2 <- gsub("Beta-low", replacement = "Beta", h3k4me1_integrated$CellType2)

data<- h3k4me1_integrated

DefaultAssay(data)<-"peaks"
annot<-Annotation(data)
a<-FindRegion(data, "Gcg", extend.downstream=100,extend.upstream=1000)
b<-FindRegion(data, "Ins1", extend.downstream=100,extend.upstream=1000)
c<-FindRegion(data, "Sst", extend.downstream=100,extend.upstream=1000)
e<-FindRegion(data, "Mki67", extend.downstream=100,extend.upstream=1000)
f<-FindRegion(data, "Ppy", extend.downstream=100,extend.upstream=1000)

paired2 = c("Alpha"="#6A3E98","Beta"="#34A047","Delta"="#B15928","Mki67-beta"="#FCBF6E","PP"="#F59899")

Annotation(data)<-annot[which(annot$gene_name=="Gcg"),]
p1<-CoveragePlot(data,window=600, idents=c("Alpha","Beta","Delta","PP","Mki67-beta"),group.by="CellType2",annotation=F,peaks=F, region =a,links = F) & scale_fill_manual(values=paired2)& theme(strip.text.y.left = element_blank(),strip.background = element_blank(),axis.text.x=element_blank(),plot.title = element_text(hjust = 0.5), axis.title.y=element_text(size=10))&ggtitle("Gcg")
anPlot<-AnnotationPlot(data, region=a)&scale_color_manual(values="darkgreen")& theme(strip.text.y.left = element_blank(),strip.background = element_blank(), axis.text.x=element_blank(),plot.title = element_text(hjust = 0.5))
p1<-CombineTracks(list(p1,anPlot),heights=c(5,1))

Annotation(data)<-annot[which(annot$gene_name=="Ins1"),]
p2<-CoveragePlot(data,window=100, idents=c("Alpha","Beta","Delta","PP","Mki67-beta"),group.by="CellType2",annotation=F,peaks=F, region =b,links = F)& scale_fill_manual(values=paired2)& theme(strip.text.y.left = element_blank(),strip.background = element_blank(),axis.text.x=element_blank(),plot.title = element_text(hjust = 0.5))&ggtitle("Ins1")
anPlot<-AnnotationPlot(data, region=b)& theme(strip.text.y.left = element_blank(),strip.background = element_blank(), axis.text.x=element_blank(),plot.title = element_text(hjust = 0.5))
p2<-CombineTracks(list(p2,anPlot),heights=c(5,1))

Annotation(data)<-annot[which(annot$gene_name=="Sst"),]
p3<-CoveragePlot(data,window=100, idents=c("Alpha","Beta","Delta","PP","Mki67-beta"),group.by="CellType2",annotation=F,peaks=F, region =c,links = F)& scale_fill_manual(values=paired2)& theme(strip.text.y.left = element_blank(),strip.background = element_blank(),axis.text.x=element_blank(),plot.title = element_text(hjust = 0.5))&ggtitle("Sst")
anPlot<-AnnotationPlot(data, region=c)&scale_color_manual(values="darkgreen")& theme(strip.text.y.left = element_blank(),strip.background = element_blank(), axis.text.x=element_blank(),plot.title = element_text(hjust = 0.5))
p3<-CombineTracks(list(p3,anPlot),heights=c(5,1))

Annotation(data)<-annot[which(annot$gene_name=="Mki67"),]
p5<-CoveragePlot(data,window=600, idents=c("Alpha","Beta","Delta","PP","Mki67-beta"),group.by="CellType2",annotation=F,peaks=F, region =e,links = F)& scale_fill_manual(values=paired2)& theme(strip.text.y.left = element_blank(),strip.background = element_blank(),axis.text.x=element_blank(),plot.title = element_text(hjust = 0.5))&ggtitle("Mki67")
anPlot<-AnnotationPlot(data, region=e)& theme(strip.text.y.left = element_blank(),strip.background = element_blank(), axis.text.x=element_blank(),plot.title = element_text(hjust = 0.5))
p5<-CombineTracks(list(p5,anPlot),heights=c(5,1))

Annotation(data)<-annot[which(annot$gene_name=="Ppy"),]
p6<-CoveragePlot(data,window=100, idents=c("Alpha","Beta","Delta","PP","Mki67-beta"), group.by="CellType2",annotation=F,peaks=F, region =f,links = F)& scale_fill_manual(values=paired2)& theme(strip.text.y.left = element_blank(),strip.background = element_blank(), axis.text.x=element_blank(),plot.title = element_text(hjust = 0.5))&ggtitle("Ppy")
anPlot<-AnnotationPlot(data, region=f)& theme(strip.text.y.left = element_blank(),strip.background = element_blank(), axis.text.x=element_blank(),plot.title = element_text(hjust = 0.5))
p6<-CombineTracks(list(p6,anPlot),heights=c(5,1))


pdf(paste0(fig_dir, "/","fig3/1.H3K4me1-track.03.pdf"), width=12, height=7)
plot_grid(p1,p2,p3,p5,p6, ncol=5)
dev.off()






















