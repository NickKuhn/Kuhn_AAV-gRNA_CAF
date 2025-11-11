# MANUSCRIPT TITLE: Discovery of a Linked Constellation of Gene Expression Revealed by Local Editing of Fibroblasts in Tumors
# Companion code: Figure S4
# Date: 2025-11-05

library(Seurat)
library(stringr)
library(ggplot2)
library(patchwork)
library(dplyr)
library(readxl)
library(EnhancedVolcano)
library(homologene)
library(reshape2)

# STEP 1: PROCESS RAW ----
# Processing and merging of demultiplexed sample_raw_feature_bc_matrix.h5 files to generate global Seurat object (shown below)
# Alternatively: pull .h5 files from CellRanger output 'per_sample_outs' folder
gTrac_A1.dataraw <- Read10X_h5("raw_feature_bc_matrix_h5_files_of_demultiplexed_samples/gTrac_A1_sample_raw_feature_bc_matrix.h5")
gTrac_A1 <- CreateSeuratObject(counts = gTrac_A1.dataraw, min.cells = 3, min.features = 150, project = 'gTrac_A1')
rm(gTrac_A1.dataraw)
gTrac_A2.dataraw <- Read10X_h5("raw_feature_bc_matrix_h5_files_of_demultiplexed_samples/gTrac_A2_sample_raw_feature_bc_matrix.h5")
gTrac_A2 <- CreateSeuratObject(counts = gTrac_A2.dataraw, min.cells = 3, min.features = 150, project = 'gTrac_A2')
rm(gTrac_A2.dataraw)
gTrac_B1.dataraw <- Read10X_h5("raw_feature_bc_matrix_h5_files_of_demultiplexed_samples/gTrac_B1_sample_raw_feature_bc_matrix.h5")
gTrac_B1 <- CreateSeuratObject(counts = gTrac_B1.dataraw, min.cells = 3, min.features = 150, project = 'gTrac_B1')
rm(gTrac_B1.dataraw)
gTrac_B2.dataraw <- Read10X_h5("raw_feature_bc_matrix_h5_files_of_demultiplexed_samples/gTrac_B2_sample_raw_feature_bc_matrix.h5")
gTrac_B2 <- CreateSeuratObject(counts = gTrac_B2.dataraw, min.cells = 3, min.features = 150, project = 'gTrac_B2')
rm(gTrac_B2.dataraw)
gOsmr_A3.dataraw <- Read10X_h5("raw_feature_bc_matrix_h5_files_of_demultiplexed_samples/gOsmr_A3_sample_raw_feature_bc_matrix.h5")
gOsmr_A3 <- CreateSeuratObject(counts = gOsmr_A3.dataraw, min.cells = 3, min.features = 150, project = 'gOsmr_A3')
rm(gOsmr_A3.dataraw)
gOsmr_A4.dataraw <- Read10X_h5("raw_feature_bc_matrix_h5_files_of_demultiplexed_samples/gOsmr_A4_sample_raw_feature_bc_matrix.h5")
gOsmr_A4 <- CreateSeuratObject(counts = gOsmr_A4.dataraw, min.cells = 3, min.features = 150, project = 'gOsmr_A4')
rm(gOsmr_A4.dataraw)
gOsmr_B3.dataraw <- Read10X_h5("raw_feature_bc_matrix_h5_files_of_demultiplexed_samples/gOsmr_B3_sample_raw_feature_bc_matrix.h5")
gOsmr_B3 <- CreateSeuratObject(counts = gOsmr_B3.dataraw, min.cells = 3, min.features = 150, project = 'gOsmr_B3')
rm(gOsmr_B3.dataraw)
gOsmr_B4.dataraw <- Read10X_h5("raw_feature_bc_matrix_h5_files_of_demultiplexed_samples/gOsmr_B4_sample_raw_feature_bc_matrix.h5")
gOsmr_B4 <- CreateSeuratObject(counts = gOsmr_B4.dataraw, min.cells = 3, min.features = 150, project = 'gOsmr_B4')
rm(gOsmr_B4.dataraw)
gTgfbr2_A5.dataraw <- Read10X_h5("raw_feature_bc_matrix_h5_files_of_demultiplexed_samples/gTgfbr2_A5_sample_raw_feature_bc_matrix.h5")
gTgfbr2_A5 <- CreateSeuratObject(counts = gTgfbr2_A5.dataraw, min.cells = 3, min.features = 150, project = 'gTgfbr2_A5')
rm(gTgfbr2_A5.dataraw)
gTgfbr2_A6.dataraw <- Read10X_h5("raw_feature_bc_matrix_h5_files_of_demultiplexed_samples/gTgfbr2_A6_sample_raw_feature_bc_matrix.h5")
gTgfbr2_A6 <- CreateSeuratObject(counts = gTgfbr2_A6.dataraw, min.cells = 3, min.features = 150, project = 'gTgfbr2_A6')
rm(gTgfbr2_A6.dataraw)
gTgfbr2_B5.dataraw <- Read10X_h5("raw_feature_bc_matrix_h5_files_of_demultiplexed_samples/gTgfbr2_B5_sample_raw_feature_bc_matrix.h5")
gTgfbr2_B5 <- CreateSeuratObject(counts = gTgfbr2_B5.dataraw, min.cells = 3, min.features = 150, project = 'gTgfbr2_B5')
rm(gTgfbr2_B5.dataraw)
gTgfbr2_B6.dataraw <- Read10X_h5("raw_feature_bc_matrix_h5_files_of_demultiplexed_samples/gTgfbr2_B6_sample_raw_feature_bc_matrix.h5")
gTgfbr2_B6 <- CreateSeuratObject(counts = gTgfbr2_B6.dataraw, min.cells = 3, min.features = 150, project = 'gTgfbr2_B6')
rm(gTgfbr2_B6.dataraw)
gIl1r1_A7.dataraw <- Read10X_h5("raw_feature_bc_matrix_h5_files_of_demultiplexed_samples/gIl1r1_A7_sample_raw_feature_bc_matrix.h5")
gIl1r1_A7 <- CreateSeuratObject(counts = gIl1r1_A7.dataraw, min.cells = 3, min.features = 150, project = 'gIl1r1_A7')
rm(gIl1r1_A7.dataraw)
gIl1r1_A8.dataraw <- Read10X_h5("raw_feature_bc_matrix_h5_files_of_demultiplexed_samples/gIl1r1_A8_sample_raw_feature_bc_matrix.h5")
gIl1r1_A8 <- CreateSeuratObject(counts = gIl1r1_A8.dataraw, min.cells = 3, min.features = 150, project = 'gIl1r1_A8')
rm(gIl1r1_A8.dataraw)
gIl1r1_B7.dataraw <- Read10X_h5("raw_feature_bc_matrix_h5_files_of_demultiplexed_samples/gIl1r1_B7_sample_raw_feature_bc_matrix.h5")
gIl1r1_B7 <- CreateSeuratObject(counts = gIl1r1_B7.dataraw, min.cells = 3, min.features = 150, project = 'gIl1r1_B7')
rm(gIl1r1_B7.dataraw)
gIl1r1_B8.dataraw <- Read10X_h5("raw_feature_bc_matrix_h5_files_of_demultiplexed_samples/gIl1r1_B8_sample_raw_feature_bc_matrix.h5")
gIl1r1_B8 <- CreateSeuratObject(counts = gIl1r1_B8.dataraw, min.cells = 3, min.features = 150, project = 'gIl1r1_B8')
rm(gIl1r1_B8.dataraw)

bigmerge <- merge(gTrac_A1, y = c(gTrac_A2, gTrac_B1, gTrac_B2,
                                  gOsmr_A3, gOsmr_A4, gOsmr_B3, gOsmr_B4,
                                  gTgfbr2_A5, gTgfbr2_A6, gTgfbr2_B5, gTgfbr2_B6,
                                  gIl1r1_A7, gIl1r1_A8, gIl1r1_B7, gIl1r1_B8),
                  project = 'bigmerge')
objects_to_remove <- c('gTrac_A1', 'gTrac_A2', 'gTrac_B1', 'gTrac_B2', 'gOsmr_A3', 'gOsmr_A4', 'gOsmr_B3', 'gOsmr_B4', 'gTgfbr2_A5', 'gTgfbr2_A6', 'gTgfbr2_B5', 'gTgfbr2_B6', 'gIl1r1_A7', 'gIl1r1_A8', 'gIl1r1_B7', 'gIl1r1_B8')
rm(list = objects_to_remove)

#add sampleID & batch number column
sampleID <- str_sub(bigmerge@meta.data$orig.ident, -2,-1)
bigmerge@meta.data$sampleID <- sampleID
bigmerge@meta.data$batch <- sampleID
bigmerge@meta.data$batch <- substr(bigmerge@meta.data$batch, 1, nchar(bigmerge@meta.data$batch) -1) #remove last character from column entry
#add AAV group
bigmerge@meta.data$aav <- bigmerge@meta.data$orig.ident
bigmerge@meta.data$aav <- as.character(bigmerge@meta.data$aav)
bigmerge@meta.data$aav <- substr(bigmerge@meta.data$aav, 1, nchar(bigmerge@meta.data$aav) -3)

bigmerge <- PercentageFeatureSet(bigmerge, pattern = "^mt-", col.name = "percent.mt")
bigmerge <- subset(bigmerge, subset = percent.mt <20)
bigmerge[["RNA"]] <- JoinLayers(bigmerge[["RNA"]])

#Visualize data:
bigmerge <- NormalizeData(bigmerge)
bigmerge <- FindVariableFeatures(bigmerge)
bigmerge <- ScaleData(bigmerge, vars.to.regress = 'percent.mt')
bigmerge <- RunPCA(bigmerge)
ElbowPlot(bigmerge, ndims = 40, reduction = 'pca')
bigmerge <- RunUMAP(bigmerge, dims=1:23, seed.use = 21212)
bigmerge <- FindNeighbors(bigmerge, dims=1:23, verbose = F)
bigmerge_0.2 <- FindClusters(bigmerge, verbose = TRUE, algorithm = 1, resolution = 0.2, random.seed = 21212)  
DimPlot(bigmerge_0.2, reduction='umap', label = T) + labs(color = "0.2")
FeaturePlot(bigmerge_0.2, features = c('Egfp', 'Pdgfra', 'Col1a2', 'Rgs5')) #Fibroblasts (Egfp, Pdgfra, Col1a2), pericytes (Rgs5)
FeaturePlot(bigmerge_0.2, features = c('Krt14', 'Krt5', 'Cdh1', 'Msln')) #epidermal basal cells (Krt14, Krt5), tumor cells (Cdh1, Msln)
FeaturePlot(bigmerge_0.2, features = c('Pecam1', 'Trac', 'Ncr1', 'S100a8')) #endothelial cells (Pecam1), T cells (Trac), NK cells (Ncr1), neutrophils (S100a8)
FeaturePlot(bigmerge_0.2, features = c('Mcpt4', 'H2-Ab1', 'Lyz2')) #mast cells (Mcpt4), myeloid cells (H2-Ab1, Lyz2)
annotations <- c(
  "tumor_1", #0
  "neutrophils", #1
  "myeloid_1", #2
  "tumor_2", #3
  "fibroblasts", #4
  "tumor_3", #5
  "TNK", #6
  "myeloid_2", #7 
  "endothelial cells", #8 
  "tumor_4", #9
  "myeloid_3", #10
  "mast cells", #11
  "pericytes", #12
  "epidermal basal cells" #13
)
names(annotations) <- levels(bigmerge_0.2)
bigmerge_0.2 <- RenameIdents(bigmerge_0.2, annotations)
bigmerge_0.2@meta.data$clustering_1 <- Idents(bigmerge_0.2)
DimPlot(bigmerge_0.2, group.by = 'clustering_1')

# STEP 2: Figure S4----
# import doublet-scrubbed subsets and re-merge to create global Seurat object
tumor_sub <- readRDS("R_objects/Seurat_objects/tumor_sub.rds") # tumor cells, downsampled to 10K total cells
tumor_sub@meta.data$cell.class <- 'tumor_cells'
tnk_clean_1 <- readRDS("R_objects/Seurat_objects/tnk_clean_1.rds") # T and NK cells
tnk_clean_1$cell.class <- 'TNK'
neuts_0.4 <- readRDS("R_objects/Seurat_objects/neuts_0.4.rds") # neutrophils
mmaav_0.4 <- readRDS("R_objects/Seurat_objects/mmaav_0.4.rds") # myeloid cells
mmaav_0.4@meta.data$cell.class <- 'MonoMac'
dcs <- readRDS("R_objects/Seurat_objects/dcs_0.2.RData") # dendritic cells
dcs@meta.data$cell.class <- 'DCs'
mastcells_clean_0.3 <- readRDS("R_objects/Seurat_objects/mastcells_clean_0.3.RData") # mast cells
mastcells_clean_0.3@meta.data$cell.class <- 'MastCells'
fbaav_0.4 <- readRDS("R_objects/Seurat_objects/fbaav_0.4.rds") # fibroblasts
fbaav_0.4@meta.data$cell.class <- 'fibroblasts'
endo_clean_0.3 <- readRDS("R_objects/Seurat_objects/endo_clean_0.3.RData") # endothelial cells
endo_clean_0.3@meta.data$cell.class <- 'endothelial_cells'
pericytes <- readRDS("R_objects/Seurat_objects/pericytes.rds") # pericytes
pericytes@meta.data$cell.class <- 'pericytes'

aav_object <- merge(tumor_sub, c(tnk_clean_1, neuts_0.4, mmaav_0.4, dcs, mastcells_clean_0.3, fbaav_0.4, endo_clean_0.3, pericytes), project = 'aav_object')
aav_object <- JoinLayers(aav_object)

rm(tumor_sub)
rm(tnk_clean_1)
rm(neuts_0.4)
rm(mmaav_0.4)
rm(dcs)
rm(mastcells_clean_0.3)
rm(fbaav_0.4)
rm(endo_clean_0.3)
rm(pericytes)

#Visualize data:
aav_object <- NormalizeData(aav_object)
aav_object <- FindVariableFeatures(aav_object)
aav_object <- ScaleData(aav_object, vars.to.regress = 'percent.mt')
aav_object <- RunPCA(aav_object)
ElbowPlot(aav_object, ndims = 40, reduction = 'pca')
aav_object <- RunUMAP(aav_object, dims=1:15, seed.use = 21212)
aav_object <- FindNeighbors(aav_object, dims=1:15, verbose = F) #k.param = 20 is default
Idents(aav_object) <- 'cell.class'
DimPlot(aav_object, label = T) + NoLegend() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

## Figure S4B----
# UMAP of all collected cells by scRNA-seq
my_levels <- c('fibroblasts', 'pericytes', 'endothelial_cells', 'MonoMac', 'neutrophils', 'DCs', 'MastCells', 'TNK', 'tumor_cells')
pdf("Figure_S4B.pdf", width = 10, height = 10)
DimPlot(aav_object, reduction = 'umap', label = T, label.size = 10.5, repel = F, raster = F, pt.size = 0.1) + NoLegend()  +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'), axis.line = element_line(linewidth = 0.1))
dev.off()

## Figure S4C----
# Bubble plot of scaled marker gene expression in main cell types
my_levels <- c('fibroblasts', 'pericytes', 'endothelial_cells', 'MonoMac', 'neutrophils', 'DCs', 'MastCells', 'TNK', 'tumor_cells')
my_levels <- rev(my_levels)
aav_object[['cell.class']] <- factor(x = aav_object@meta.data$cell.class, levels = my_levels)
Idents(aav_object) <- 'cell.class'

pdf("Figure_S4C.pdf", width = 4.2, height = 1.8)
DotPlot(aav_object,
        features = c('Egfp', 'Pdgfra', 'Col1a2', 'Rgs5', 'Cspg4', 'Pecam1', 'Csf1r', 'Lyz2', 'Itgam', 'Itgb2', 'S100a8', 'Cxcr2', 'H2-Ab1', 'Flt3', 'Mcpt4', 'Trac', 'Nkg7', 'Cdh1', 'Msln'),
        cols = 'RdBu', assay = "RNA") +
  theme(axis.text.x = element_text(angle=45, hjust = 1, size = 6), axis.text.y=element_text(size = 6), text = element_text(size = 6)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.ticks = element_line(linewidth = 0.1)) +
  theme(axis.line = element_line(linewidth = 0.1))
dev.off()

## Figure S4D----
# UMAP visualization of marker gene expression in all captured cells
pdf("Figure_S4D.pdf", width = 3.5, height = 2.5)
fp <- FeaturePlot(aav_object, c('Pdgfra', 'Rgs5', 'Pecam1', 'Lyz2', 'S100a8', 'H2-Ab1', 'Mcpt4', 'Trac', 'Msln'), order = T,
                  cols = c('lightgrey', 'red'))
fp <- lapply(fp, function(fp) fp +
               theme(plot.title = element_text(size = 6, face = 'italic'), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
                     legend.text = element_text(size = 6), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.1, 'lines'), legend.margin = margin(0.1, 0.1, 0.1, 0.1, 'mm'), legend.spacing = unit(0, 'mm'),
                     plot.margin = margin(0, 0, 0, 0, unit = 'pt'), axis.line = element_line(linewidth = 0.1)))
fp <- wrap_plots(fp, ncol = 3)
fp
dev.off()

## Figure S4E----
# Violin plot depicting mean expression of TGFbeta response signature from Mariathasan et al., 2018 in CAFs grouped by AAV-gRNA condition
fbaav_0.4 <- readRDS("R_objects/Seurat_objects/fbaav_0.4.rds") # fibroblasts
fbaav_0.4 <- AddModuleScore(fbaav_0.4, features = c('Acta2', 'Actg2', 'Adam12', 'Adam19', 'Cnn1', 'Col4a1',
                                                    'Ccn2', 'Ctps', 'Rflnb', 'Fstl3',
                                                    'Pxdc1', 'Sema7a', 'Sh3pxd2a', 'Tagln', 'Tgfbi', 'Tns1', 'Tpm1'),
                            assay = 'RNA', name = 'tgf.sig.Mariathasan')
pdf('Figure_S4E.pdf', width = 1.5, height = 1.2, bg = 'transparent')
sg1 <- VlnPlot(fbaav_0.4, 'tgf.sig.Mariathasan1', group.by = 'aav', pt.size = 0) + 
  ggtitle('TGF-beta response signature\n(Mariathasan et al., 2018)') +
  labs(y = 'arbitrary unit', title = 'TGF-beta response signature\n(Mariathasan et al., 2018)') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 6), axis.text = element_text(size = 6), axis.text.y = element_blank(),
        axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1),
        plot.title = element_text(size = 6, face = 'plain', margin = margin(0,0,0,0)),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA)) +
  NoLegend()
sg1$layers[[1]]$aes_params$size <- 0.1 #adjust linewidth of violins
sg1
dev.off()

## Figure S4F----
# UMAP visualization of CAFs from Krishnamurty et al. 2022
krishnamurty_tgfbr2ko <- readRDS("R_objects/public_datasets/krishnamurty_tgfbr2ko.rds")
fb_deg_markers <- read_excel("Supplementary_Tables/Supplementary Table 1 DEG scRNAseq.xlsx", sheet = 'fibroblasts')
fb.deg.col18a1 <- fb_deg_markers %>% filter(clustering_1 == 'Col18a1') %>% filter(p_val_adj < 0.05 & avg_log2FC > sqrt(2)) %>% arrange(desc(avg_log2FC)) %>% slice_head(n=50) %>% pull(gene)
genes_to_remove <- c('Lratd1', 'Gucy1a1', 'Smim41', 'Dlg2') # genes not found in Krishnamurty data set
#FeaturePlot(krishnamurty_tgfbr2ko, genes_to_remove)
fb.deg.col18a1 <- fb.deg.col18a1[!fb.deg.col18a1 %in% genes_to_remove]
krishnamurty_tgfbr2ko <- AddModuleScore(krishnamurty_tgfbr2ko, features = fb.deg.col18a1, assay = 'RNA', name = 'fb.deg.col18a1')

png("Figure_S4F.png", width = 1290, height = 3000, res = 600, bg = 'transparent')
gg1 <- DimPlot(krishnamurty_tgfbr2ko, group.by = 'seurat_clusters') + ggtitle(element_blank()) + theme(axis.text = element_blank(), axis.ticks = element_blank(), text = (element_text(size = 12)),legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.1, 'lines'), axis.line = element_line(linewidth = 0.1), plot.margin = margin(0, 0, 0, 0, unit = 'pt'), axis.title = element_blank())
gg2 <- DimPlot(krishnamurty_tgfbr2ko, pt.size = 0.1, cells.highlight = colnames(krishnamurty_tgfbr2ko)[krishnamurty_tgfbr2ko$condition=='WT'], sizes.highlight = 0.1) + labs(title='Dpt-wt/wt Tgfbr2-fl/fl') + NoLegend() + theme(axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size = 12, face = 'plain', margin=margin(0,0,0,0)), text = element_text(size = 12), axis.line = element_line(linewidth = 0.1), plot.margin = margin(0, 0, 0, 0, unit = 'pt'), axis.title = element_blank())
gg3 <- DimPlot(krishnamurty_tgfbr2ko, pt.size = 0.1, cells.highlight = colnames(krishnamurty_tgfbr2ko)[krishnamurty_tgfbr2ko$condition=='CKO'], sizes.highlight = 0.1) + labs(title='Dpt-ki/ki Tgfbr2-fl/fl') + NoLegend() + theme(axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size = 12, face = 'plain', margin=margin(0,0,0,0)), text = element_text(size = 12), axis.line = element_line(linewidth = 0.1), plot.margin = margin(0, 0, 0, 0, unit = 'pt'), axis.title = element_blank())
gg4 <- FeaturePlot(krishnamurty_tgfbr2ko, 'fb.deg.col18a11', order = T) + ggtitle('Col18a1 cluster sig.') + theme(axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size = 12, face = 'plain', margin=margin(0,0,0,0)), text = element_text(size = 12), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.1, 'lines'), axis.line = element_line(linewidth = 0.1), plot.margin = margin(0, 0, 0, 0, unit = 'pt'), axis.title = element_blank())
#gg4 <- FeaturePlot(fibro, 'Col18a1', order = T) + ggtitle('Col18a1') + theme(axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size = 12, face = 'plain', margin=margin(0,0,0,0)), text = element_text(size = 12), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.1, 'lines'), axis.line = element_line(linewidth = 0.1), plot.margin = margin(0, 0, 0, 0, unit = 'pt'), axis.title = element_blank())
krishnamurty_combined_plots <- (gg1 / gg2 / gg3 / gg4) + plot_annotation(title = 'Krishnamurty, 2022') & theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank())
krishnamurty_combined_plots
dev.off()

## Figure S4G----
# Volcano plots of differentially expressed genes (DEG, pseudobulk scRNA-Seq) in CAFs 
DEG_MAST_fb_gOsmr <- read_excel("Supplementary_Tables/Supplementary Table 2 DEG pseudobulk.xlsx", sheet = 'gOsmr')
DEG_MAST_fb_gTgfbr2 <- read_excel("Supplementary_Tables/Supplementary Table 2 DEG pseudobulk.xlsx", sheet = 'gTgfbr2')
DEG_MAST_fb_gIl1r1 <- read_excel("Supplementary_Tables/Supplementary Table 2 DEG pseudobulk.xlsx", sheet = 'gIl1r1')

# VOLCANO PLOT ---> gOsmr
lab_italics_gOsmr <- paste0("italic('", DEG_MAST_fb_gOsmr$primerid, "')")
selectLab_italics_gOsmr = paste0(
  "italic('", c('Tnn', 'Crispld2', 'Col23a1', 'Ltbp2', 'Sema3f', 'Antxr1', 'Sod2', 'Bid','Lmo4','Bmpr2',
                'Stard13', 'Klf5', 'Rarb', 'Pfkl', 'Osmr', 'Ccdc80', 'Gpi1', 'P4ha1'),"')")
vp1 <- EnhancedVolcano(DEG_MAST_fb_gOsmr,
                       lab = lab_italics_gOsmr,
                       labSize = 1.5,
                       selectLab = selectLab_italics_gOsmr,
                       x = 'coef',
                       y='fdr',
                       title = 'AAV-gOsmr',
                       subtitle = 'UP in ctrl    UP in gOsmr',
                       caption = element_blank(),
                       #legendLabels = element_blank(),
                       FCcutoff = 0,
                       pCutoff = 0.05,
                       cutoffLineWidth = 0.1,
                       cutoffLineCol = 'black',
                       col = c('grey', 'grey', 'grey', '#1A85FF'),
                       pointSize = 0.5,
                       colAlpha = 0.5,
                       drawConnectors = T,
                       min.segment.length = 0,
                       widthConnectors = 0.1,
                       max.overlaps = 20,
                       #boxedLabels = T,
                       parseLabels = T,
                       xlim = c(-1, 1),
                       ylim = c(0, -(log10(10e-7))),
                       #xlab = "log2FC",
                       #ylab = "log10(Padj)"
) +
  theme_classic() +
  ylab(expression(-log[10](Padj))) +
  xlab(expression(log[2](FC))) +
  theme(axis.text.x = element_text(size = 6, color = 'black'), axis.text.y = element_text(size = 6, color='black'),
        axis.title = element_text(size = 6),
        plot.title=element_text(face='plain', hjust=0.5, size = 6, margin = margin(b=1)), plot.subtitle = element_text(face='plain', hjust = 0.5, size=6, margin = margin(b=1)),
        axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1),
        plot.margin = margin(0, 1, 0, 0, unit = 'mm'),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA)) +
  NoLegend()
# VOLCANO PLOT ---> gTgfbr2
lab_italics_gTgfbr2 <- paste0("italic('", DEG_MAST_fb_gTgfbr2$primerid, "')")
selectLab_italics_gTgfbr2 = paste0(
  "italic('", c('Sav1','Cxcl5','Hgf', 'Dcn', 'Col18a1', 'Lum','P2rx4','Lgals3','S100a4', 'Cav1',
                'Lrrc15', 'Mmp13', 'Col8a1', 'Acta2', 'Dhrs3', 'Tagln', 'Fbn2', 'Loxl3', 'Pmepa1', 'Gask1b'),"')")
vp2 <- EnhancedVolcano(DEG_MAST_fb_gTgfbr2,
                       lab = lab_italics_gTgfbr2,
                       labSize = 1.5,
                       selectLab = selectLab_italics_gTgfbr2,
                       x = 'coef',
                       y='fdr',
                       title = 'AAV-gTgfbr2',
                       subtitle = 'UP in ctrl    UP in gTgfbr2',
                       caption = element_blank(),
                       #legendLabels = element_blank(),
                       FCcutoff = 0,
                       pCutoff = 0.05,
                       cutoffLineWidth = 0.1,
                       cutoffLineCol = 'black',
                       col = c('grey', 'grey', 'grey', '#D41159'),
                       pointSize = 0.5,
                       colAlpha = 0.5,
                       drawConnectors = T,
                       min.segment.length = 0,
                       widthConnectors = 0.1,
                       max.overlaps = 20,
                       #boxedLabels = T,
                       parseLabels = T,
                       xlim = c(-2, 2),
                       ylim = c(0, -(log10(10e-7))),
                       #xlab = "log2FC",
                       #ylab =
) +
  theme_classic() +
  ylab(expression(-log[10](Padj))) +
  xlab(expression(log[2](FC))) +
  theme(axis.text.x = element_text(size = 6, color = 'black'), axis.text.y = element_text(size = 6, color='black'),
        axis.title = element_text(size = 6),
        plot.title=element_text(face='plain', hjust=0.5, size = 6, margin = margin(b=1)), plot.subtitle = element_text(face='plain', hjust = 0.5, size=6, margin = margin(b=1)),
        axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1),
        plot.margin = margin(0, 1, 0, 0, unit = 'mm'),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA)) +
  NoLegend()
# VOLCANO PLOT ---> gIl1r1
lab_italics_gIl1r1 <- paste0("italic('", DEG_MAST_fb_gIl1r1$primerid, "')")
selectLab_italics_gIl1r1 = paste0(
  "italic('", c('Slc7a6', 'Utp11', 'Gramd1b', 'Mcpt2', 'Clic4', 'Ubxn6', 'Thra', 'Zfp395', 'Tspan8', 'Krt7',
                'Ptprj', 'Aggf1', 'Exoc5', 'Creb3', 'Chsy1', 'Ier3ip1', 'Psmd14'),"')")
vp3 <- EnhancedVolcano(DEG_MAST_fb_gIl1r1,
                       lab = lab_italics_gIl1r1,
                       labSize = 1.5,
                       selectLab = selectLab_italics_gIl1r1,
                       x = 'coef',
                       y='fdr',
                       title = 'AAV-gIl1r1',
                       subtitle = 'UP in ctrl    UP in gIl1r1',
                       caption = element_blank(),
                       #legendLabels = element_blank(),
                       FCcutoff = 0,
                       pCutoff = 0.05,
                       cutoffLineWidth = 0.1,
                       cutoffLineCol = 'black',
                       col = c('grey', 'grey', 'grey', '#40B0A6'),
                       pointSize = 0.5,
                       colAlpha = 0.5,
                       drawConnectors = T,
                       min.segment.length = 0,
                       widthConnectors = 0.1,
                       max.overlaps = 50,
                       #boxedLabels = T,
                       parseLabels = T,
                       xlim = c(-1, 1),
                       ylim = c(0, -(log10(10e-7))),
                       #xlab = "log2FC",
                       #ylab =
) +
  theme_classic() +
  ylab(expression(-log[10](Padj))) +
  xlab(expression(log[2](FC))) +
  theme(axis.text.x = element_text(size = 6, color = 'black'), axis.text.y = element_text(size = 6, color='black'),
        axis.title = element_text(size = 6),
        plot.title=element_text(face='plain', hjust=0.5, size = 6, margin = margin(b=1)), plot.subtitle = element_text(face='plain', hjust = 0.5, size=6, margin = margin(b=1)),
        axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1),
        plot.margin = margin(0, 1, 0, 0, unit = 'mm'),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA)) +
  NoLegend()
png("Figure_S4G.png", width = 4000, height = 1500, res = 300)
vp1 + vp2 + vp3
dev.off()

## Figure S4H----
# Venn diagram of up- and downregulated DEGs between control (ctrl) and AAV-gOsmr, AAV-gTgfbr2, or AAV-gIl1r1
DEG_MAST_fb_gOsmr_sig_up <- DEG_MAST_fb_gOsmr %>% filter(fdr < 0.05 & coef>=log2(1)) %>% pull(primerid)
DEG_MAST_fb_gOsmr_sig_down <- DEG_MAST_fb_gOsmr %>% filter(fdr < 0.05 & coef<log2(1)) %>% pull(primerid)
DEG_MAST_fb_gTgfbr2_sig_up <- DEG_MAST_fb_gTgfbr2 %>% filter(fdr < 0.05 & coef>=log2(1)) %>% pull(primerid)
DEG_MAST_fb_gTgfbr2_sig_down <- DEG_MAST_fb_gTgfbr2 %>% filter(fdr < 0.05 & coef<log2(1)) %>% pull(primerid)
DEG_MAST_fb_gIl1r1_sig_up <- DEG_MAST_fb_gIl1r1 %>% filter(fdr < 0.05 & coef>=log2(1)) %>% pull(primerid)
DEG_MAST_fb_gIl1r1_sig_down <- DEG_MAST_fb_gIl1r1 %>% filter(fdr < 0.05 & coef<log2(1)) %>% pull(primerid)

pdf('Figure_S4H_up.pdf', width = 1.25, height = 1.25, bg = 'transparent')
list_venn_up <- list('AAV-gOsmr' = DEG_MAST_fb_gOsmr_sig_up,
                     'AAV-gTgfbr2' = DEG_MAST_fb_gTgfbr2_sig_up,
                     'AAV-gIl1r1' = DEG_MAST_fb_gIl1r1_sig_up)
venn1 <- ggvenn::ggvenn(list_venn_up, c('AAV-gOsmr', 'AAV-gTgfbr2', 'AAV-gIl1r1'),
                        show_elements = F,
                        show_percentage = F,
                        #label_sep = "\n",
                        set_name_size = 2,
                        text_size = 2, #
                        stroke_size = 0.1173709, # ratio 2.13-to-1 in illustrator
                        fill_color = c('#1A85FF', '#D41159', '#40B0A6'),
                        digits = 0,
                        #auto_scale = T
) +
  ggtitle('upregulated vs. Ctrl') +
  theme(plot.title=element_text(hjust = 0.5, size = 6, margin = margin(b=0)),
        plot.margin = margin(1, 0, 1, 0, unit = 'mm')) +
  scale_x_continuous(expand = expansion(mult = .2)) + #expand ggplot plotting area to avoid cutting off of labels
  scale_y_continuous(expand = expansion(mult = .2))#expand ggplot plotting area to avoid cutting off of labels
venn1
dev.off()

pdf('Figure_S4H_down.pdf', width = 1.25, height = 1.25)
list_venn_down <- list('AAV-gOsmr' = DEG_MAST_fb_gOsmr_sig_down,
                     'AAV-gTgfbr2' = DEG_MAST_fb_gTgfbr2_sig_down,
                     'AAV-gIl1r1' = DEG_MAST_fb_gIl1r1_sig_down)
venn2 <- ggvenn::ggvenn(list_venn_down, c('AAV-gOsmr', 'AAV-gTgfbr2', 'AAV-gIl1r1'),
                        show_elements = F,
                        show_percentage = F,
                        #label_sep = "\n",
                        set_name_size = 2,
                        text_size = 2, #
                        stroke_size = 0.1173709, # ratio 2.13-to-1 in illustrator
                        fill_color = c('#1A85FF', '#D41159', '#40B0A6'),
                        digits = 0,
                        #auto_scale = T
) +
  ggtitle('downregulated vs. Ctrl') +
  theme(plot.title=element_text(hjust = 0.5, size = 6, margin = margin(b=0)),
        plot.margin = margin(1, 0, 1, 0, unit = 'mm')) +
  scale_x_continuous(expand = expansion(mult = .2)) + #expand ggplot plotting area to avoid cutting off of labels
  scale_y_continuous(expand = expansion(mult = .2))#expand ggplot plotting area to avoid cutting off of labels
venn2
dev.off()

## Figure S4I----
# Enriched gene sets (Gene Ontology hallmark) based on the up- and downregulated DEGs between control and AAV-gOsmr
#devtools::install_version("msigdbr", version = "7.2.1", repos = "http://cran.us.r-project.org")
library(msigdbr)
#packageVersion('msigdbr') # 7.2.1
#remotes::install_github("montilab/hypeR@devel")
library(hypeR)
#packageVersion('hypeR') # 2.0.0

HALLMARK <- msigdb_gsets(species = 'Mus musculus', category = 'H')

DEG_MAST_fb_gOsmr_sig_up <- DEG_MAST_fb_gOsmr %>% filter(fdr < 0.05 & coef>=log2(1)) %>% pull(primerid)
DEG_MAST_fb_gOsmr_sig_down <- DEG_MAST_fb_gOsmr %>% filter(fdr < 0.05 & coef<log2(1)) %>% pull(primerid)

hyp.gsea.hm.fb.deg.o.up <- hypeR(DEG_MAST_fb_gOsmr_sig_up, genesets = HALLMARK, test = 'hypergeometric', background=19405)
gsea.hm.fb.deg.o.up <- hyp_dots(hyp.gsea.hm.fb.deg.o.up, top=10, fdr=0.05, size_by = 'significance') #dots plot
gsea.hm.fb.deg.o.up
hyp.gsea.hm.fb.deg.o.down <- hypeR(DEG_MAST_fb_gOsmr_sig_down, genesets = HALLMARK, test = 'hypergeometric', background=19405)
gsea.hm.fb.deg.o.down <- hyp_dots(hyp.gsea.hm.fb.deg.o.down, top=10, fdr=0.05, size_by = 'significance') #dots plot
gsea.hm.fb.deg.o.down

gsea.hm.fb.deg.o.up.down <- list('osmr_up' = DEG_MAST_fb_gOsmr_sig_up,
                                 'osmr_down' = DEG_MAST_fb_gOsmr_sig_down)
hyp.gsea.hm.fb.deg.o.up.down <- hypeR(gsea.hm.fb.deg.o.up.down, genesets = HALLMARK, test = 'hypergeometric', background=19405)
gsea.hm.fb.deg.o.up.down <- hyp_dots(hyp.gsea.hm.fb.deg.o.up.down, top=10, fdr=0.05, size_by = 'significance') #dots plot

temp.osmr.up <- gsea.hm.fb.deg.o.up[['data']][, c('label', 'fdr')]
temp.osmr.up$neg_log10_adj_pval <- -log10(temp.osmr.up$fdr)
temp.osmr.down <- gsea.hm.fb.deg.o.down[['data']][, c('label', 'fdr')]
temp.osmr.down$neg_log10_adj_pval <- -log10(temp.osmr.down$fdr)
temp.osmr.down$neg_log10_adj_pval <- -temp.osmr.down$neg_log10_adj_pval
merge.osmr.up.down <- rbind(temp.osmr.up, temp.osmr.down)
merge.osmr.up.down$label <- gsub("HALLMARK_", "", merge.osmr.up.down$label)
merge.osmr.up.down$label <- factor(merge.osmr.up.down$label, levels = merge.osmr.up.down$label[order(merge.osmr.up.down$neg_log10_adj_pval)])
merge.osmr.up.down$sign <- ifelse(merge.osmr.up.down$neg_log10_adj_pval >= 0, "positive", "negative")

pdf('Figure_S4I.pdf', width = 2.75, height = 1.2, bg = 'transparent')
ggplot(merge.osmr.up.down, aes(x = label, y = neg_log10_adj_pval, fill = sign)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = 0, label = label, hjust = ifelse(neg_log10_adj_pval > 0, 1, 0)), size = 2) +
  scale_fill_manual(values = c('positive' = alpha('#1A85FF', 0.75), 'negative' = alpha('black', 0.5))) +
  coord_flip() +
  labs(title = "AAV-gOsmr (GO: hallmark)") +
  theme_classic() +
  xlab('') +
  ylab(expression(-log[10](Padj))) +
  theme(axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 6, color = 'black'), axis.text.y = element_blank(),
        axis.title = element_text(size = 6), title = element_text(size = 6),
        plot.title = element_text(size = 6, face = 'plain', margin = margin(0,0,0,0)),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA)) +
  NoLegend()
dev.off()

## Figure S4J----
# Enriched gene sets (Gene Ontology hallmark) based on the up- and downregulated DEGs between control and AAV-gTgfbr2
DEG_MAST_fb_gTgfbr2_sig_up <- DEG_MAST_fb_gTgfbr2 %>% filter(fdr < 0.05 & coef>=log2(1)) %>% pull(primerid)
DEG_MAST_fb_gTgfbr2_sig_down <- DEG_MAST_fb_gTgfbr2 %>% filter(fdr < 0.05 & coef<log2(1)) %>% pull(primerid)

hyp.gsea.hm.fb.deg.t2.up <- hypeR(DEG_MAST_fb_gTgfbr2_sig_up, genesets = HALLMARK, test = 'hypergeometric', background=19405)
gsea.hm.fb.deg.t2.up <- hyp_dots(hyp.gsea.hm.fb.deg.t2.up, top=10, fdr=0.05, size_by = 'significance', title='MAST DEGs fibro - upregulated gTgfbr2') #dots plot
gsea.hm.fb.deg.t2.up
hyp.gsea.hm.fb.deg.t2.down <- hypeR(DEG_MAST_fb_gTgfbr2_sig_down, genesets = HALLMARK, test = 'hypergeometric', background=19405)
gsea.hm.fb.deg.t2.down <- hyp_dots(hyp.gsea.hm.fb.deg.t2.down, top=10, fdr=0.05, size_by = 'significance', title='MAST DEGs fibro - downregulated gTgfbr2') #dots plot
gsea.hm.fb.deg.t2.down

temp.tgfbr2.up <- gsea.hm.fb.deg.t2.up[['data']][, c('label', 'fdr')]
temp.tgfbr2.up <- temp.tgfbr2.up %>% arrange(fdr) %>% slice_head(n = 6)
temp.tgfbr2.up$neg_log10_adj_pval <- -log10(temp.tgfbr2.up$fdr)
temp.tgfbr2.down <- gsea.hm.fb.deg.t2.down[['data']][, c('label', 'fdr')]
temp.tgfbr2.down <- temp.tgfbr2.down %>% arrange(fdr) %>% slice_head(n = 2)
temp.tgfbr2.down$neg_log10_adj_pval <- -log10(temp.tgfbr2.down$fdr)
temp.tgfbr2.down$neg_log10_adj_pval <- -temp.tgfbr2.down$neg_log10_adj_pval
merge.tgfbr2.up.down <- rbind(temp.tgfbr2.up, temp.tgfbr2.down)
merge.tgfbr2.up.down$label <- gsub("HALLMARK_", "", merge.tgfbr2.up.down$label)
merge.tgfbr2.up.down$label <- factor(merge.tgfbr2.up.down$label, levels = merge.tgfbr2.up.down$label[order(merge.tgfbr2.up.down$neg_log10_adj_pval)])
merge.tgfbr2.up.down$sign <- ifelse(merge.tgfbr2.up.down$neg_log10_adj_pval >= 0, "positive", "negative")

pdf('Figure_S4J.pdf', width = 2.75, height = 1.2, bg = 'transparent')
ggplot(merge.tgfbr2.up.down, aes(x = label, y = neg_log10_adj_pval, fill = sign)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = 0, label = label, hjust = ifelse(neg_log10_adj_pval > 0, 1, 0)), size = 2) +
  scale_fill_manual(values = c('positive' = alpha('#D41159', 0.75), 'negative' = alpha('black', 0.5))) +
  coord_flip() +
  labs(title = "AAV-gTgfbr2 (GO: hallmark)") +
  theme_classic() +
  xlab('') +
  ylab(expression(-log[10](Padj))) +
  theme(axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 6, color = 'black'), axis.text.y = element_blank(),
        axis.title = element_text(size = 6), title = element_text(size = 6),
        plot.title = element_text(size = 6, face = 'plain', margin = margin(0,0,0,0)),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA)) +
  NoLegend()
dev.off()

## Figure S4K----
# Similarity of DEGs between mouse CAF subsets (columns)
aav.fb.deg <- split(fb_deg_markers, fb_deg_markers$clustering_1)
# and human PDAC CAF subsets from Hwang et al., 2022 (rows) 
hwang.hu.adh <- read_excel("R_objects/public_datasets/cNMF_Hwang2022_CAFs.xlsx", sheet = 'Adhesive')
hwang.hu.immuno <- read_excel("R_objects/public_datasets/cNMF_Hwang2022_CAFs.xlsx", sheet = 'Immunomodulatory')
hwang.hu.myo <- read_excel("R_objects/public_datasets/cNMF_Hwang2022_CAFs.xlsx", sheet = 'Myofibroblastic')
hwang.hu.nrt <- read_excel("R_objects/public_datasets/cNMF_Hwang2022_CAFs.xlsx", sheet = 'Neurotropic')
hwang.hu.nmf <- rbind(hwang.hu.adh, hwang.hu.immuno, hwang.hu.myo, hwang.hu.nrt)
hwang.hu.nmf <- split(hwang.hu.nmf, hwang.hu.nmf$cNMF)
names(hwang.hu.nmf) <- c('0', '1', '2', '3')

# label transfer to mouse genes
for (i in 0:3) { # 0:3 is number of total clusters; adjust if necessary
  # Convert the loop counter to character to match the tibble names
  name <- as.character(i)
  
  # Extract gene names from the current tibble
  gene_list <- hwang.hu.nmf[[name]]$gene
  
  # Get the mouse orthologs using homologene
  gene_mapping <- homologene::human2mouse(gene_list)
  
  # Merge the mapping (selecting only humanGene and mouseGene) back into the tibble
  hwang.hu.nmf[[name]] <- merge(hwang.hu.nmf[[name]], 
                                gene_mapping[, c("humanGene", "mouseGene")],
                                by.x = "gene", 
                                by.y = "humanGene", 
                                all.x = TRUE)
}
names(hwang.hu.nmf) <- c('adhesive', 'immunomodulatory', 'myofibroblastic', 'neurotropic') # change back to original names

### HYPERGEOMETRIC TEST----
# Define background:
N <- 19405 # number of genes captured by FLEX kit
# read in lists
listA <- aav.fb.deg
listB <- hwang.hu.nmf

# Initialize a results data.frame
results <- data.frame(
  listA_name = character(),
  listB_name = character(),
  m = integer(),
  n = integer(),
  overlap = integer(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Nested loops for pairwise comparison
for(a in names(listA)) {
  for(b in names(listB)) {
    # select parameters for inclusion for listA 
    setA <- listA[[a]] %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>% pull(gene)

    # select parameters for inclusion for listB
    setB <- listB[[b]] %>% pull(mouseGene)

    m <- length(setA)
    n <- length(setB)
    # Count overlap between the two sets
    k <- length(intersect(setA, setB))
    
    # Calculate hypergeometric p-value
    p_val <- phyper(k - 1, m, N - m, n, lower.tail = FALSE)
    
    # Append the result
    results <- rbind(results, data.frame(
      listA_name = a,
      listB_name = b,
      m = m,
      n = n,
      overlap = k,
      p_value = p_val,
      stringsAsFactors = FALSE
    ))
  }
}

results$p_value_adj <- results$p_value*nrow(results)

deg.mat.hgt <- acast(results, listA_name ~ listB_name, value.var = "p_value_adj")

#create heatmap without hierarchical clustering
col_ord <- hclust(dist(t(deg.mat.hgt)))$order
y_levels <- colnames(deg.mat.hgt)[col_ord]

#create heatmap, first melt matrix
melt_hgt <- as.data.frame.table(deg.mat.hgt)
names(melt_hgt)[names(melt_hgt) == 'Freq'] <- 'P' 
head(melt_hgt)
melt_hgt$Var1 = as.character(melt_hgt$Var1)
melt_hgt$Var2 = as.character(melt_hgt$Var2)

#order x-axis
my_levels <- c('homeostatic', 'intermediate', 'infl_CAF', 'myCAF_1', 'myCAF_2', 'myCAF_3', 'myCAF_prolif', 'Col18a1')
melt_hgt[["Var1"]] <- factor(x = melt_hgt$Var1, levels = my_levels)
melt_hgt[["Var2"]] <- factor(x = melt_hgt$Var2, levels = y_levels)

pdf('Figure_S4K.pdf', width=2, height=0.8)
ggplot(data = melt_hgt, aes(x=Var1, y=Var2, fill=(-log10(P)))) + 
  geom_tile() +
  scale_fill_viridis_c(limits=c(0, 40), oob = scales::squish) +
  ggtitle('Hwang 2022, human PDAC CAFs') +
  theme_classic() +
  theme(plot.title = element_text(size =6, margin = margin(0,0,0,0,unit='pt')),
        axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.ticks = element_line(linewidth = 0.1),
        text = element_text(size = 6),
        axis.line = element_line(linewidth = 0.1),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'),
        legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.2, 'lines'), legend.margin = margin(0, 0, 0, 0, 'mm'), legend.spacing = unit(0, 'mm')
  )
dev.off()

## Figure S4L----
# Similarity of DEGs between mouse CAF subsets (columns) and
# human non-small cell lung carcinoma CAF subsets from Grout et al., 2022
# Gene lists from Grout et al., 2022, Supplemental Table S17
w <- 'CH25H,FRZB,TGM2,A2M,PPP1R14A,PTGDS,PLEKHH2,FMO2,PRELP,SLIT2,SCN7A,DST,SEPP1,LTBP4,ALDH1A1,ADH1B,SOD3,ABCA8,PLAC9,FBLN1,GSN,APOD,C7,C3,APOE'
adh1b <- strsplit(w, ",")[[1]]
x <- 'COL11A1,CTSK,CHN1,COL6A3,VCAN,SPARC,COL12A1,THBS2,COL1A2,CTHRC1,POSTN,DIO2,COL1A1,RIN2,MMP11,COL6A1,COL5A2,CRABP2,LRRC17,ADAM12,TWIST1,LRRC15,CILP,GJB2,TNFAIP6,INHBA,GREM1,CST1,IGFL2,MMP7,P4HA3,PLPP4'
fap <- strsplit(x, ",")[[1]]
y <- 'IGFBP2,ITGBL1,COL4A1,COL4A2,LTBP2,LTBP1,ELN,MYH11,WIF1,PCSK1N,WFDC1,MYLK,F2R,CITED2,ROBO2,LUZP2,TSLP,FAM150A,TYRP1,CST2,ANO4,ENTPD1-AS1,KRT17,WNT11,ATRNL1,COL9A1,BRF2'
myh11 <- strsplit(y, ",")[[1]]
z <- 'DERL3,CYP26A1,SULF1,C1QTNF3,CXCL14,SUGCT,PLAU,HOPX,HTRA3,EPYC,MMP1,SCUBE3'
fap.asma <- strsplit(z, ",")[[1]]
grout.deg <- list(adh1b, fap, myh11, fap.asma)

ms.grout.deg <- lapply(grout.deg, function(x) {
  v <- unlist(x, use.names = FALSE)
  m <- homologene::human2mouse(unique(v))
  lut <- setNames(m$mouseGene, m$humanGene)
  unname(lut[v])
})
names(ms.grout.deg) <- c('ADH1B+', 'FAP+', 'MYH11+ aSMA+', 'FAP+ aSMA+') # change back to original names

### HYPERGEOMETRIC TEST----
# Define background:
N <- 19405 # number of genes captured by FLEX kit
# read in lists
listA <- aav.fb.deg
listB <- ms.grout.deg

# Initialize a results data.frame
results <- data.frame(
  listA_name = character(),
  listB_name = character(),
  m = integer(),
  n = integer(),
  overlap = integer(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Nested loops for pairwise comparison
for(a in names(listA)) {
  for(b in names(listB)) {
    # select parameters for inclusion for listA
    setA <- listA[[a]] %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>% pull(gene)
   
    # select parameters for inclusion for listB
    setB <- listB[[b]]
    
    m <- length(setA)
    n <- length(setB)
    # Count overlap between the two sets
    k <- length(intersect(setA, setB))
    
    # Calculate hypergeometric p-value
    p_val <- phyper(k - 1, m, N - m, n, lower.tail = FALSE)
    
    # Append the result
    results <- rbind(results, data.frame(
      listA_name = a,
      listB_name = b,
      m = m,
      n = n,
      overlap = k,
      p_value = p_val,
      stringsAsFactors = FALSE
    ))
  }
}

# View results
print(results)

results$p_value_adj <- results$p_value*nrow(results)

deg.mat.hgt <- acast(results, listA_name ~ listB_name, value.var = "p_value_adj")

#create heatmap with manual ordering:
col_ord <- hclust(dist(t(deg.mat.hgt)))$order
y_levels <- c('ADH1B+', 'FAP+', 'FAP+ aSMA+', 'MYH11+ aSMA+') # Grout order
melt_hgt <- as.data.frame.table(deg.mat.hgt)
names(melt_hgt)[names(melt_hgt) == 'Freq'] <- 'P' 
#head(melt_hgt)
melt_hgt$Var1 = as.character(melt_hgt$Var1)
melt_hgt$Var2 = as.character(melt_hgt$Var2)

#order x-axis
my_levels <- c('homeostatic', 'intermediate', 'infl_CAF', 'myCAF_1', 'myCAF_2', 'myCAF_3', 'myCAF_prolif', 'Col18a1')
melt_hgt[["Var1"]] <- factor(x = melt_hgt$Var1, levels = my_levels)
melt_hgt[["Var2"]] <- factor(x = melt_hgt$Var2, levels = y_levels)

pdf('Figure_S4L.pdf', width=2, height=0.8)
ggplot(data = melt_hgt, aes(x=Var1, y=Var2, fill=(-log10(P)))) + 
  geom_tile() +
  scale_fill_viridis_c(limits=c(0, 10), oob = scales::squish) +
  ggtitle('Grout 2022, human myCAF, NSCLC') +
  theme_classic() +
  theme(plot.title = element_text(size =6, margin = margin(0,0,0,0,unit='pt')),
        axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.ticks = element_line(linewidth = 0.1),
        text = element_text(size = 6),
        axis.line = element_line(linewidth = 0.1),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'),
        legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.2, 'lines'), legend.margin = margin(0, 0, 0, 0, 'mm'), legend.spacing = unit(0, 'mm')
  ) + scale_y_discrete(limits = rev)
dev.off()

## Figure S4M----
# Similarity of DEGs between mouse CAF subsets (columns) and
# cross-tissue human fibroblast atlas subtypes from Gao et al, 2024
hs_fb_atlas_markers <- read_excel("R_objects/public_datasets/DEGs_Gao2024_human_PDAC.xlsx", sheet = 'Table S2A')
gao.deg <- split(hs_fb_atlas_markers, hs_fb_atlas_markers$cluster)
names(gao.deg) <- c('c01-MYH11', 'c02-ADAMDEC1', 'c03-COL15A1', 'c04-LRRC15', 'c05-PI16', 'c06-ADH1B', 'c07-IL6', 'c08-RGS5', 'c09-CTNNB1', 'c10-PRG4', 'c11-HOPX', 'c12-MSLN', 'c13-HGF', 'c14-HSPA6', 'c15-SOX6', 'c16-SFRP2', 'c17-STMN1', 'c18-HHIP', 'c19-MMP1', 'c20-CD74') # change names to annotation in paper

# label transfer to mouse genes
for (name in names(gao.deg)) {

  # Extract gene names from the current tibble
  gene_list <- gao.deg[[name]]$gene
  
  # Get the mouse orthologs using homologene
  gene_mapping <- homologene::human2mouse(gene_list)
  
  # Merge the mapping (selecting only humanGene and mouseGene) back into the tibble
  gao.deg[[name]] <- merge(gao.deg[[name]], 
                            gene_mapping[, c("humanGene", "mouseGene")],
                            by.x = "gene", 
                            by.y = "humanGene", 
                            all.x = TRUE)
}

### HYPERGEOMETRIC TEST----
# Define background:
N <- 19405 # number of genes captured by FLEX kit
# read in lists
listA <- aav.fb.deg
listB <- gao.deg

# Initialize a results data.frame
results <- data.frame(
  listA_name = character(),
  listB_name = character(),
  m = integer(),
  n = integer(),
  overlap = integer(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Nested loops for pairwise comparison
for(a in names(listA)) {
  for(b in names(listB)) {
    # select parameters for inclusion for listA
    setA <- listA[[a]] %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>% pull(gene)

    # select parameters for inclusion for listB
    setB <- listB[[b]] %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>% pull(mouseGene)

    m <- length(setA)
    n <- length(setB)
    # Count overlap between the two sets
    k <- length(intersect(setA, setB))
    
    # Calculate hypergeometric p-value
    p_val <- phyper(k - 1, m, N - m, n, lower.tail = FALSE)
    
    # Append the result
    results <- rbind(results, data.frame(
      listA_name = a,
      listB_name = b,
      m = m,
      n = n,
      overlap = k,
      p_value = p_val,
      stringsAsFactors = FALSE
    ))
  }
}

results$p_value_adj <- results$p_value*nrow(results)

deg.mat.hgt <- acast(results, listA_name ~ listB_name, value.var = "p_value_adj")

#create heatmap without hierarchical clustering
col_ord <- hclust(dist(t(deg.mat.hgt)))$order
y_levels <- c('c07-IL6', 'c13-HGF', 'c18-HHIP','c02-ADAMDEC1', 'c15-SOX6', 'c10-PRG4', 'c20-CD74', 'c17-STMN1', 'c11-HOPX', 'c12-MSLN', 'c06-ADH1B', 'c01-MYH11', 'c09-CTNNB1', 'c19-MMP1', 'c14-HSPA6', 'c08-RGS5', 'c04-LRRC15', 'c16-SFRP2', 'c03-COL15A1', 'c05-PI16') # Gao order

#create heatmap, first melt matrix
melt_hgt <- as.data.frame.table(deg.mat.hgt)
names(melt_hgt)[names(melt_hgt) == 'Freq'] <- 'P' 
melt_hgt$Var1 = as.character(melt_hgt$Var1)
melt_hgt$Var2 = as.character(melt_hgt$Var2)

#order x-axis
my_levels <- c('homeostatic', 'intermediate', 'infl_CAF', 'myCAF_1', 'myCAF_2', 'myCAF_3', 'myCAF_prolif', 'Col18a1')
melt_hgt[["Var1"]] <- factor(x = melt_hgt$Var1, levels = my_levels)
melt_hgt[["Var2"]] <- factor(x = melt_hgt$Var2, levels = y_levels)

pdf('Figure_S4M.pdf', width=1.8, height=1.8)
ggplot(data = melt_hgt, aes(x=Var1, y=Var2, fill=(-log10(P)))) + 
  geom_tile() +
  scale_fill_viridis_c(limits=c(0, 60), oob = scales::squish) +
  #scale_fill_viridis_c() +
  ggtitle('Gao 2024, human fibroblasts') +
  #scale_fill_gradient2(low = 'blue', high='red', mid = 'white', midpoint = 0) +
  theme_classic() +
  theme(plot.title = element_text(size =6, margin = margin(0,0,0,0,unit='pt')),
        axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.ticks = element_line(linewidth = 0.1),
        text = element_text(size = 6),
        axis.line = element_line(linewidth = 0.1),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'),
        legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.2, 'lines'), legend.margin = margin(0, 0, 0, 0, 'mm'), legend.spacing = unit(0, 'mm')
  )
dev.off()
