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
library(ggrepel)
plan(multicore)
options(future.globals.maxSize = 100000 * 1024^2) # set 50G RAM
set.seed(1)
#===================================================================================================================
fig_dir <- ("/research/labs/bme/weiz/m238739/Paired-Tag/data/data16/mouse/figures")
out_dir <- ("/research/labs/bme/weiz/m238739/Paired-Tag/data/data16/mouse/bam/RNA/merged_RNA")
fig_dir10 <- ("/research/labs/bme/weiz/m238739/Paired-Tag/data/data16/mouse/figures/Manuscript_V1/Fig/figs5")
#==================================================================================================================
setwd(out_dir)
dir()
ptag_list_integrated <- readRDS(file = paste0(out_dir, "/","ptag_list_integrated.rds"))
markers <- readRDS(file = paste0(out_dir, "/","ptag_list_integrated.markers.rds"))
# Get cell identity classes
Idents(object = ptag_list_integrated)
levels(x = ptag_list_integrated)
ptag_list_Beta <-subset(x = ptag_list_integrated, idents = c("Beta-hi", "Beta-low","Mki67-beta"))
Idents(object = ptag_list_Beta)
ptag_list_Beta@meta.data
DefaultAssay(ptag_list_Beta) <- "RNA"
#------------------------------------------------------------------
ptag_list_Beta_0Week <-subset(x = ptag_list_Beta, subset = Week == "0")
ptag_list_Beta_8Week <-subset(x = ptag_list_Beta, subset = Week == "8")
ptag_list_Beta_16Week <-subset(x = ptag_list_Beta, subset = Week == "16")

ptag_list_Beta_0Week.MvsF <- FindMarkers(ptag_list_Beta_0Week, ident.1 = "M", ident.2 = "F",group.by = 'Sex',logfc.threshold = 0)
ptag_list_Beta_8Week.MvsF <- FindMarkers(ptag_list_Beta_8Week, ident.1 = "M", ident.2 = "F",group.by = 'Sex',logfc.threshold = 0)
ptag_list_Beta_16Week.MvsF <- FindMarkers(ptag_list_Beta_16Week, ident.1 = "M", ident.2 = "F",group.by = 'Sex',logfc.threshold = 0)

write.table(ptag_list_Beta_0Week.MvsF, paste0(fig_dir10, "/","ptag_list_Beta_0Week.MvsF2.csv"),col.names = TRUE, row.names = TRUE, sep = "\t", quote = TRUE)
write.table(ptag_list_Beta_8Week.MvsF, paste0(fig_dir10, "/","ptag_list_Beta_8Week.MvsF2.csv"),col.names = TRUE, row.names = TRUE, sep = "\t", quote = TRUE)
write.table(ptag_list_Beta_16Week.MvsF, paste0(fig_dir10, "/","ptag_list_Beta_16Week.MvsF2.csv"),col.names = TRUE, row.names = TRUE, sep = "\t", quote = TRUE)

#saveRDS(ptag_list_Beta_0Week.MvsF,paste0(fig_dir9, "/","ptag_list_Beta_0Week.MvsF2.rds"))
#saveRDS(ptag_list_Beta_8Week.MvsF,paste0(fig_dir9, "/","ptag_list_Beta_8Week.MvsF2.rds"))
#saveRDS(ptag_list_Beta_16Week.MvsF,paste0(fig_dir9, "/","ptag_list_Beta_16Week.MvsF2.rds"))
#--------------------------------------------------------
ptag_list_Beta_0Week.MvsF <- readRDS(file = paste0(fig_dir9, "/","ptag_list_Beta_0Week.MvsF2.rds"))
ptag_list_Beta_8Week.MvsF <- readRDS(file = paste0(fig_dir9, "/","ptag_list_Beta_8Week.MvsF2.rds"))
ptag_list_Beta_16Week.MvsF <- readRDS(file = paste0(fig_dir9, "/","ptag_list_Beta_16Week.MvsF2.rds"))
#---------------------------------------------------------
dim(ptag_list_Beta_0Week.MvsF[ptag_list_Beta_0Week.MvsF$p_val_adj<0.05,])
dim(ptag_list_Beta_8Week.MvsF[ptag_list_Beta_8Week.MvsF$p_val_adj<0.05,])
dim(ptag_list_Beta_16Week.MvsF[ptag_list_Beta_16Week.MvsF$p_val_adj<0.05,])

a1 <-ptag_list_Beta_0Week.MvsF[ptag_list_Beta_0Week.MvsF$p_val_adj<0.05,]
a2 <-ptag_list_Beta_8Week.MvsF[ptag_list_Beta_8Week.MvsF$p_val_adj<0.05,]
a3 <-ptag_list_Beta_16Week.MvsF[ptag_list_Beta_16Week.MvsF$p_val_adj<0.05,]

# Combine the dataframes into one
combined_df <- rbind(a1, a2, a3)
# Extract unique genes from the first column
unique_genes <- unique(rownames(combined_df))
# Count the unique genes
unique_gene_count <- length(unique_genes)

dim(ptag_list_Beta_0Week.MvsF[abs(ptag_list_Beta_0Week.MvsF$avg_log2FC)>0.3219281,])
dim(ptag_list_Beta_8Week.MvsF[abs(ptag_list_Beta_8Week.MvsF$avg_log2FC)>0.3219281,])
dim(ptag_list_Beta_16Week.MvsF[abs(ptag_list_Beta_16Week.MvsF$avg_log2FC)>0.3219281,])

a4 <-ptag_list_Beta_0Week.MvsF[abs(ptag_list_Beta_0Week.MvsF$avg_log2FC)>0.3219281,]
a5 <-ptag_list_Beta_8Week.MvsF[abs(ptag_list_Beta_8Week.MvsF$avg_log2FC)>0.3219281,]
a6 <-ptag_list_Beta_16Week.MvsF[abs(ptag_list_Beta_16Week.MvsF$avg_log2FC)>0.3219281,]

# Combine the dataframes into one
combined_df <- rbind(a4, a5, a6)
# Extract unique genes from the first column
unique_genes <- unique(rownames(combined_df))
# Count the unique genes
unique_gene_count <- length(unique_genes)##[1] 87

write.table(unique_genes, paste0(fig_dir10, "/","ptag_list_Beta_unique_genes.MvsF2.csv"),
            col.names = TRUE, row.names = TRUE, sep = "\t", quote = TRUE)
#--------------------------------------------------------------------------------------------------
dim(ptag_list_Beta_0Week.MvsF[abs(ptag_list_Beta_0Week.MvsF$avg_log2FC) > 0.3219281 & -log10(ptag_list_Beta_0Week.MvsF$p_val_adj) > 1.30103,])
dim(ptag_list_Beta_8Week.MvsF[abs(ptag_list_Beta_8Week.MvsF$avg_log2FC) > 0.3219281 & -log10(ptag_list_Beta_8Week.MvsF$p_val_adj) > 1.30103,])
dim(ptag_list_Beta_16Week.MvsF[abs(ptag_list_Beta_16Week.MvsF$avg_log2FC) > 0.3219281 & -log10(ptag_list_Beta_16Week.MvsF$p_val_adj) > 1.30103,])

a7 <-ptag_list_Beta_0Week.MvsF[abs(ptag_list_Beta_0Week.MvsF$avg_log2FC) > 0.3219281 & -log10(ptag_list_Beta_0Week.MvsF$p_val_adj) > 1.30103,]
a8 <-ptag_list_Beta_8Week.MvsF[abs(ptag_list_Beta_8Week.MvsF$avg_log2FC) > 0.3219281 & -log10(ptag_list_Beta_8Week.MvsF$p_val_adj) > 1.30103,]
a9 <-ptag_list_Beta_16Week.MvsF[abs(ptag_list_Beta_16Week.MvsF$avg_log2FC) > 0.3219281 & -log10(ptag_list_Beta_16Week.MvsF$p_val_adj) > 1.30103,]

# Combine the dataframes into one
combined_df <- rbind(a7, a8, a9)
# Extract unique genes from the first column
unique_genes <- unique(rownames(combined_df))
# Count the unique genes
unique_gene_count <- length(unique_genes)##[1] #84

write.table(unique_genes, paste0(fig_dir10, "/","ptag_list_Beta_unique_genes.MvsF_p0.05log2FC0.3.csv"),
            col.names = TRUE, row.names = TRUE, sep = "\t", quote = TRUE)
#----------------------------------------------------------------------------------------------------------
# Filter genes for each condition
male_upregulated_0Week <- ptag_list_Beta_0Week.MvsF[ptag_list_Beta_0Week.MvsF$avg_log2FC > 0.3219281 & -log10(ptag_list_Beta_0Week.MvsF$p_val_adj) > 1.30103,]
female_upregulated_0Week <- ptag_list_Beta_0Week.MvsF[ptag_list_Beta_0Week.MvsF$avg_log2FC < -0.3219281 & -log10(ptag_list_Beta_0Week.MvsF$p_val_adj) > 1.30103,]

male_upregulated_8Week <- ptag_list_Beta_8Week.MvsF[ptag_list_Beta_8Week.MvsF$avg_log2FC > 0.3219281 & -log10(ptag_list_Beta_8Week.MvsF$p_val_adj) > 1.30103,]
female_upregulated_8Week <- ptag_list_Beta_8Week.MvsF[ptag_list_Beta_8Week.MvsF$avg_log2FC < -0.3219281 & -log10(ptag_list_Beta_8Week.MvsF$p_val_adj) > 1.30103,]

male_upregulated_16Week <- ptag_list_Beta_16Week.MvsF[ptag_list_Beta_16Week.MvsF$avg_log2FC > 0.3219281 & -log10(ptag_list_Beta_16Week.MvsF$p_val_adj) > 1.30103,]
female_upregulated_16Week <- ptag_list_Beta_16Week.MvsF[ptag_list_Beta_16Week.MvsF$avg_log2FC < -0.3219281 & -log10(ptag_list_Beta_16Week.MvsF$p_val_adj) > 1.30103,]

# Combine the dataframes
male_combined <- rbind(male_upregulated_0Week, male_upregulated_8Week, male_upregulated_16Week)
female_combined <- rbind(female_upregulated_0Week, female_upregulated_8Week, female_upregulated_16Week)

# Extract unique genes
unique_male_genes <- unique(rownames(male_combined))
unique_female_genes <- unique(rownames(female_combined))

# Combine unique genes into a single dataframe
unique_genes_df <- data.frame(Gene = c(unique_male_genes, unique_female_genes),
                              Group = c(rep("Male", length(unique_male_genes)), rep("Female", length(unique_female_genes))))
# Write the result to a file
write.table(unique_genes_df, file = paste0(fig_dir10, "/", "ptag_list_Beta_unique_genes_grouped.MvsF_p0.05log2FC0.3_sex.csv"),
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = TRUE)

#---------------------------------------------------------------------------------------------------
Beta_0Week.MvsF<-ptag_list_Beta_0Week.MvsF
Beta_8Week.MvsF<-ptag_list_Beta_8Week.MvsF
Beta_16Week.MvsF<-ptag_list_Beta_16Week.MvsF

table(Beta_0Week.MvsF$diffexpressed)
table(Beta_8Week.MvsF$diffexpressed)
table(Beta_16Week.MvsF$diffexpressed)
# add a column of NAs
Beta_0Week.MvsF$diffexpressed <- "NO"
Beta_0Week.MvsF$diffexpressed[Beta_0Week.MvsF$avg_log2FC > 0.3219281 & -log10(Beta_0Week.MvsF$p_val_adj) > 1.30103] <- "UP"
Beta_0Week.MvsF$diffexpressed[Beta_0Week.MvsF$avg_log2FC < -0.3219281 & -log10(Beta_0Week.MvsF$p_val_adj) > 1.30103] <- "DOWN"

# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
Beta_0Week.MvsF$delabel <- NA
Beta_0Week.MvsF$delabel[Beta_0Week.MvsF$diffexpressed != "NO"] <- rownames(Beta_0Week.MvsF)[Beta_0Week.MvsF$diffexpressed != "NO"]
#-----------------------------------------------------------------------------------
pdf(paste0(fig_dir10, "/","Beta_0Week.MvsF_topGENES3_v01.pdf"), height = 8, width = 15)
ggplot(data=Beta_0Week.MvsF, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel(size = 6) +
  scale_color_manual(values=c("blue","gray", "red")) + ylim(0, 40)+theme_classic() + ggtitle("Beta_0Week.MvsF")+
  geom_vline(xintercept=c(-0.3219281, 0.3219281), col="red",linetype = "longdash")+geom_hline(yintercept=1.30103, col="red",linetype = "longdash") +
  coord_cartesian(clip = "off")+
  labs(x = "Log2 Fold Change", y = "-Log10(Adjusted P-value)")
dev.off()
#-----------------------------------------------------------------------------------
# add a column of NAs
Beta_8Week.MvsF$diffexpressed <- "NO"
Beta_8Week.MvsF$diffexpressed[Beta_8Week.MvsF$avg_log2FC > 0.3219281 & -log10(Beta_8Week.MvsF$p_val_adj) > 1.30103] <- "UP"
Beta_8Week.MvsF$diffexpressed[Beta_8Week.MvsF$avg_log2FC < -0.3219281 & -log10(Beta_8Week.MvsF$p_val_adj) > 1.30103] <- "DOWN"

# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
Beta_8Week.MvsF$delabel <- NA
Beta_8Week.MvsF$delabel[Beta_8Week.MvsF$diffexpressed != "NO"] <- rownames(Beta_8Week.MvsF)[Beta_8Week.MvsF$diffexpressed != "NO"]

#-------------------------------------------------------------------------------
pdf(paste0(fig_dir10, "/","Beta_8Week.MvsF_topGENES3_v01.pdf"), height = 8, width = 15)
ggplot(data=Beta_8Week.MvsF, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel(size = 6) +
  scale_color_manual(values=c("blue","gray", "red")) + ylim(0, 300)+theme_classic() + ggtitle("Beta_8Week.MvsF")+
  geom_vline(xintercept=c(-0.3219281, 0.3219281), col="red",linetype = "longdash") +geom_hline(yintercept=1.30103, col="red",linetype = "longdash")+
  coord_cartesian(clip = "off")+
  labs(x = "Log2 Fold Change", y = "-Log10(Adjusted P-value)")
dev.off()

#-------------------------------------------------------------------------------
# add a column of NAs
Beta_16Week.MvsF$diffexpressed <- "NO"
Beta_16Week.MvsF$diffexpressed[Beta_16Week.MvsF$avg_log2FC > 0.3219281 & -log10(Beta_16Week.MvsF$p_val_adj) > 1.30103] <- "UP"
Beta_16Week.MvsF$diffexpressed[Beta_16Week.MvsF$avg_log2FC < -0.3219281 & -log10(Beta_16Week.MvsF$p_val_adj) > 1.30103] <- "DOWN"

# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
Beta_16Week.MvsF$delabel <- NA
Beta_16Week.MvsF$delabel[Beta_16Week.MvsF$diffexpressed != "NO"] <- rownames(Beta_16Week.MvsF)[Beta_16Week.MvsF$diffexpressed != "NO"]
#------------------------------------------------------------------------
pdf(paste0(fig_dir10, "/","Beta_16Week.MvsF_topGENES3_v01.pdf"), height = 8, width = 15)
ggplot(data=Beta_16Week.MvsF, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel(size = 6) +
  scale_color_manual(values=c("blue","gray", "red")) + ylim(0, 300)+theme_classic() + ggtitle("Beta_16Week.MvsF")+
  geom_vline(xintercept=c(-0.3219281, 0.3219281), col="red",linetype = "longdash") +geom_hline(yintercept=1.30103, col="red",linetype = "longdash")+
  coord_cartesian(clip = "off")+
  labs(x = "Log2 Fold Change", y = "-Log10(Adjusted P-value)")
dev.off()
