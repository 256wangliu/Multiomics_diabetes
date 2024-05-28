library(Cairo)
library(Seurat)
library(dplyr)
library(magrittr)
library(patchwork)
library(reshape2)
library(BiocParallel)
library(RColorBrewer)
library(reshape2)
library(magick)
library(readr)
library(Matrix)
library(clustree)
library(Signac)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggrepel)
library(tidyverse)
library(ReactomePA)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(DOSE)
library(purrr)
library(AnnotationHub)
set.seed(1)
#---------------------------------------------------------------------------------------------------------------
out_dir <- (paste0(proj_dir,"figures/Manuscript_V1/Fig/fig1/data"))
fig_dir8 <- ("/research/labs/bme/weiz/m238739/Paired-Tag/data/data16/mouse/figures/Manuscript_V1/Fig/figs6/data")

ptag_list_integrated <- readRDS(file = paste0(out_dir, "/","ptag_list_integrated.rds"))
markers <- readRDS(file = paste0(out_dir, "/","ptag_list_integrated.markers.rds"))
# Get cell identity classes
Idents(object = ptag_list_integrated)
levels(x = ptag_list_integrated)
ptag_list_Beta <-subset(x = ptag_list_integrated, idents = c("Beta-hi", "Beta-low","Mki67-beta"))
Idents(object = ptag_list_Beta)
ptag_list_Beta.markers1 <- FindMarkers(ptag_list_Beta, ident.1 = "Mki67-beta", ident.2 = "Beta-hi")
ptag_list_Beta.markers2 <- FindMarkers(ptag_list_Beta, ident.1 = "Mki67-beta", ident.2 = "Beta-low")
#saveRDS(ptag_list_Beta.markers1,paste0(fig_dir8, "/","ptag_list_Beta.markers1.rds"))
#saveRDS(ptag_list_Beta.markers2,paste0(fig_dir8, "/","ptag_list_Beta.markers2.rds"))

ptag_list_Beta.markers1 <- readRDS(file = paste0(fig_dir8, "/","ptag_list_Beta.markers1.rds"))
ptag_list_Beta.markers2 <- readRDS(file = paste0(fig_dir8, "/","ptag_list_Beta.markers2.rds"))

pdf(paste0(fig_dir8, "/","mki67_topGENES001.pdf"), height = 10, width = 15)
VlnPlot(ptag_list_Beta,
        features = c("Cenpp", "Knl1","Smc2","Cep128","Kif11","Rad51b"),
        cols =c( "mediumorchid3", "coral3","darkgoldenrod1"),
        sort = F) + NoLegend()

dev.off()


Mki67vsBeta_HI_de<- ptag_list_Beta.markers1
Mki67vsBeta_LO_de<- ptag_list_Beta.markers2
#--------------------------------------------------------------------------------------------------------------------------
# add a column of NAs
Mki67vsBeta_HI_de$diffexpressed <- "NO"
Mki67vsBeta_HI_de$diffexpressed[Mki67vsBeta_HI_de$avg_log2FC > 1 & -log10(Mki67vsBeta_HI_de$p_val_adj) > 2] <- "UP"
Mki67vsBeta_HI_de$diffexpressed[Mki67vsBeta_HI_de$avg_log2FC < -1 & -log10(Mki67vsBeta_HI_de$p_val_adj) > 2] <- "DOWN"

Mki67vsBeta_HI_de$delabel <- NA
Mki67vsBeta_HI_de$delabel[Mki67vsBeta_HI_de$avg_log2FC >=13 | Mki67vsBeta_HI_de$avg_log2FC <=-20] <- rownames(Mki67vsBeta_HI_de)[Mki67vsBeta_HI_de$avg_log2FC >=13 | Mki67vsBeta_HI_de$avg_log2FC <=-20]

pdf(paste0(fig_dir8, "/","Mki67vsBeta_HI_topGENES_001.pdf"), height = 10, width = 15)
ggplot(data=Mki67vsBeta_HI_de, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed,label = delabel)) +
  geom_point(size = 2) +
  theme_minimal() +
  geom_text_repel(size = 8) +
  scale_color_manual(values=c("darkgreen", "gray", "darkred")) +
  geom_vline(xintercept=c(-1, 1), col="darkred",linetype = "longdash") +
  geom_hline(yintercept=1, col="darkred",linetype = "longdash")+ ylim(0, 300)+theme_classic() +
  ggtitle("Mki67vsBeta_HI")+coord_cartesian(clip = "off")+
labs(x = "Log2 Fold Change", y = "-Log10(Adjusted P-value)")

dev.off()
#----------------------------------------------------------------------------------------------------------------------------------
Mki67vsBeta_LO_de$diffexpressed <- "NO"
Mki67vsBeta_LO_de$diffexpressed[Mki67vsBeta_LO_de$avg_log2FC > 1 & -log10(Mki67vsBeta_LO_de$p_val_adj) > 2] <- "UP"
Mki67vsBeta_LO_de$diffexpressed[Mki67vsBeta_LO_de$avg_log2FC < -1 & -log10(Mki67vsBeta_LO_de$p_val_adj) > 2] <- "DOWN"
Mki67vsBeta_LO_de$delabel <- NA
Mki67vsBeta_LO_de$delabel[Mki67vsBeta_LO_de$avg_log2FC >=13 | Mki67vsBeta_LO_de$avg_log2FC <=-20] <- rownames(Mki67vsBeta_LO_de)[Mki67vsBeta_LO_de$avg_log2FC >=13 | Mki67vsBeta_LO_de$avg_log2FC <=-20]

pdf(paste0(fig_dir8, "/","Mki67vsBeta_LO_topGENES_001.pdf"), height = 10, width = 15)
ggplot(data=Mki67vsBeta_LO_de, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point(size = 2) +
  theme_minimal() +
  geom_text_repel(size = 8) +
  scale_color_manual(values=c("darkgreen", "gray", "darkred")) +
  geom_vline(xintercept=c(-1, 1), col="darkred",linetype = "longdash") +
  geom_hline(yintercept=1, col="darkred",linetype = "longdash")+ ylim(0, 300)+theme_classic() +
  ggtitle("Mki67vsBeta_LO")+coord_cartesian(clip = "off")+
  labs(x = "Log2 Fold Change", y = "-Log10(Adjusted P-value)")

dev.off()
#-------------------------------------------------------------------------------
#Get gene names per cluster
deg.ls <- split(rownames(markers), f = markers$cluster)
class(deg.ls)
#Transfer gene symbol into entrez id
library(AnnotationHub)
geneid.ls <- deg.ls %>% map( ~{

  gene.df <- AnnotationDbi::select(org.Mm.eg.db,keys = .x,columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")

  gene <- gene.df$ENTREZID
  gene <- gene[which(!is.na(gene))]
  gene <- unique(gene)

  return(gene)
} )

gene.ls <- geneid.ls[c("Beta-hi", "Beta-low","Mki67-beta")]

compKEGG <- compareCluster(geneCluster   = gene.ls,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH",
                           organism = "mmu")
compGO <- compareCluster(geneCluster   = gene.ls,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb = org.Mm.eg.db,
                         ont = 'BP')
## dot plot
pdf(paste0(fig_dir8, "/","GO1.pdf"))
g1 <- dotplot(compGO, showCategory = 5, title = "GO Enrichment Analysis")
g1
dev.off()



