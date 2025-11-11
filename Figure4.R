# MANUSCRIPT TITLE: Discovery of a Linked Constellation of Gene Expression Revealed by Local Editing of Fibroblasts in Tumors
# Companion code: FIGURE 4
# Title: Tgfbr2 knockout-induced Col18a1hi CAFs are distinct from iCAFs and myCAFs and correlate with worse survival in human PDAC patients
# Date: 2025-11-06

library(monocle3)
library(Seurat)
library(dplyr)
library(SeuratWrappers)
library(corrplot)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(msigdbr)
library(hypeR)
library(data.table)
library(readxl)
library(survival)
library(survminer)

## Figure 4A----
# UMAP of CAF subsets without myCAF_prolif subset and pseudotime trajectory starting at the homeostatic fibroblast subset and ending at cluster myCAF_3 or Col18a1
fbaav_0.4 <- readRDS("R_objects/Seurat_objects/fbaav_0.4.rds") # fibroblast Seurat object

Idents(fbaav_0.4) <- 'clustering_1' #remove the proliferation clusters from fibroblast object
fb = subset(fbaav_0.4, idents = ('myCAF_prolif'), invert = T)

fb <- RunPCA(fb)
fb <- RunUMAP(fb, reduction = 'pca', dims=1:18, seed.use = 21212, n.neighbors = 50)

cds <- as.cell_data_set(fb)
cds <- cluster_cells(cds)

plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
cds <- learn_graph(cds,learn_graph_control = list(nn.k = 50))
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
cds <- order_cells(cds)
mp2 <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,graph_label_size = 5,trajectory_graph_segment_size = 2,cell_size = 0.8,
  trajectory_graph_color = 'black',
  label_branch_points = F,
  label_roots = T,
  label_leaves = F
)
mp2
Idents(fb) <- 'clustering_1'

png("Figure_4A_top.png", width = 800, height = 800, res = 300, bg='transparent')
mp1 <- DimPlot(fb, label = F, label.size = 2, cols = c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", "#FF61CC")) +
  NoLegend() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank(),
        plot.margin = margin(0,0,0,0, unit = 'pt'),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA))
mp1
dev.off()

png("Figure_4A_bottom.png", width = 800, height = 800, res = 300, bg='transparent')
mp2 <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE, graph_label_size = 5, trajectory_graph_segment_size = 1, cell_size = 0.5,
  trajectory_graph_color = 'white',
  label_branch_points = F,
  label_roots = T,
  label_leaves = F,
  alpha = 0.5) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.ticks = element_blank(), axis.line.y = element_blank(), axis.line.x = element_blank(),
        legend.position = 'none',
        #legend.text = element_text(size = 6), legend.title = element_text(size = 6), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.3, 'lines'), legend.margin = margin(0, 0, 0, 0, 'mm'), legend.spacing = unit(0, 'mm'),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA))
mp2
dev.off()

## Figure 4B----
# Violin plot of CAF subsets plotted according to their distribution in pseudotime
pseudotime <- cds@principal_graph_aux@listData$UMAP$pseudotime
fb <- AddMetaData(fb, pseudotime, col.name = 'pseudotime')
Idents(fb) <- 'clustering_1'
pdf('Figure_4B.pdf', width = 1.2, height = 2, bg = 'transparent')
gg1 <- VlnPlot(fb, 'pseudotime', pt.size = 0) + NoLegend() +
  labs(y='pseudotime', title = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y=element_text(size = 6), axis.text = element_text(size = 6),
        axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA))
gg1$layers[[1]]$aes_params$size <- 0.1 #adjust linewidth of violins
gg1
dev.off()

## Figure 4C----
# Pairwise Pearson correlation coefficient of CAF NMF factor expression across all CAFs. Rows and columns are ordered by hierarchical clustering
fb_18_programs_diag_L1_0_0 <- readRDS("Supplementary_Tables/fb_18_programs_diag_L1_0_0.rds")
df <- fb_18_programs_diag_L1_0_0[["h"]]
df_transposed <- t(df)
colnames(df_transposed) <- paste0("factor_", seq_len(ncol(df_transposed)))
correlation_matrix <- cor(df_transposed, method = "pearson")
correlation_df <- as.data.frame(correlation_matrix)
#print(correlation_df)
testRes = cor.mtest(correlation_df)
correlation_df <- as.matrix(correlation_df)
pdf('Figure_4C.pdf', width = 5, height = 5)
corrplot(correlation_df, p.mat = testRes$p, sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 1, pch.col = 'black', insig = 'label_sig', method = 'circle', order = 'hclust',
         type = 'lower', col = rev(brewer.pal(10 ,"RdBu")), tl.col = 'black', tl.cex = 1, col.lim = c(-1,1),
         title = 'Pearson correlation coefficient\nacross cells',
         mar = c(0,0,2,0))
dev.off()

## Figure 4D----
# UMAP visualization of select NMF factor expression in CAFs
h =fb_18_programs_diag_L1_0_0$h
colnames(h)<-colnames(fbaav_0.4)
row.names(h)<-paste('factor-',as.character(c(1:nrow(h))),sep ='')
fbaav_0.4[["NMF"]] <- CreateAssayObject(counts = h)

png("Figure_4D.png", width = 5000, height = 5000, res = 800, bg = 'transparent')
p5 <- FeaturePlot(fbaav_0.4, 'factor-5', max.cutoff = 5e-04, pt.size = 0.01) +
  scale_color_viridis_c(option = 'viridis', direction = 1) +
  labs(title = '', color = 'weight') +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        legend.text = element_blank(), legend.title = element_text(size = 6),
        plot.title = element_text(size = 6, face = 'plain'),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA)) +
  NoLegend()
p5
p6 <- FeaturePlot(fbaav_0.4, 'factor-6', max.cutoff = 5e-04, pt.size = 0.01) +
  scale_color_viridis_c(option = 'viridis', direction = 1) +
  labs(title = '', color = 'weight') +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        legend.text = element_blank(), legend.title = element_text(size = 6),
        plot.title = element_text(size = 6, face = 'plain'),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA)) +
  NoLegend()
p6
p7 <- FeaturePlot(fbaav_0.4, 'factor-7', max.cutoff = 5e-04, pt.size = 0.01) +
  scale_color_viridis_c(option = 'viridis', direction = 1) +
  labs(title = '', color = 'weight') +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        legend.text = element_blank(), legend.title = element_text(size = 6),
        plot.title = element_text(size = 6, face = 'plain'),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA)) +
  NoLegend()
p7
p9 <- FeaturePlot(fbaav_0.4, 'factor-9', max.cutoff = 5e-04, pt.size = 0.01) +
  scale_color_viridis_c(option = 'viridis', direction = 1) +
  labs(title = '', color = 'weight') +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        legend.text = element_blank(), legend.title = element_text(size = 6),
        plot.title = element_text(size = 6, face = 'plain'),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA)) +
  NoLegend()
p9
p10 <- FeaturePlot(fbaav_0.4, 'factor-10', max.cutoff = 5e-04, pt.size = 0.01) +
  scale_color_viridis_c(option = 'viridis', direction = 1) +
  labs(title = '', color = 'weight') +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        legend.text = element_blank(), legend.title = element_text(size = 6),
        plot.title = element_text(size = 6, face = 'plain'),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA)) +
  NoLegend()
p10
p11 <- FeaturePlot(fbaav_0.4, 'factor-11', max.cutoff = 5e-04, pt.size = 0.01) +
  scale_color_viridis_c(option = 'viridis', direction = 1) +
  labs(title = '', color = 'weight') +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        legend.text = element_blank(), legend.title = element_text(size = 6),
        plot.title = element_text(size = 6, face = 'plain'),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA)) +
  NoLegend()
p11
p13 <- FeaturePlot(fbaav_0.4, 'factor-13', max.cutoff = 5e-04, pt.size = 0.01) +
  scale_color_viridis_c(option = 'viridis', direction = 1) +
  labs(title = '', color = 'weight') +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        legend.text = element_blank(), legend.title = element_text(size = 6),
        plot.title = element_text(size = 6, face = 'plain'),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA)) +
  NoLegend()
p13
p14 <- FeaturePlot(fbaav_0.4, 'factor-14', max.cutoff = 5e-04, pt.size = 0.01) +
  scale_color_viridis_c(option = 'viridis', direction = 1) +
  labs(title = '', color = 'weight') +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        legend.text = element_blank(), legend.title = element_text(size = 6),
        plot.title = element_text(size = 6, face = 'plain'),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA)) +
  NoLegend()
p14
transparent_spacer <- plot_spacer() +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA)
  )
(p10 | p13 | p7) / (p5 | p11 | p14) / (p6 | p9 | transparent_spacer)
dev.off()


## Figure 4E----
# Violin plots of NMF factor_6 and NMF factor_9 expression in the CAF subsets
pdf('Figure_4E.pdf', width = 1.2, height = 2.5, bg = 'transparent')
gg1 <- VlnPlot(fbaav_0.4, 'factor-6', pt.size = 0) + NoLegend() +
  labs(y = 'factor\nweight', title = 'factor_6') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 6), axis.text = element_text(size = 6), axis.text.y = element_blank(),
        axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1),
        plot.title = element_text(size = 6, face = 'plain', margin = margin(0,0,0,0)),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA))
gg1$layers[[1]]$aes_params$size <- 0.1
gg2 <- VlnPlot(fbaav_0.4, 'factor-9', pt.size = 0) + NoLegend() +
  labs(y = 'factor\nweight', title = 'factor_9') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 6), axis.text = element_text(size = 6), axis.text.y = element_blank(),
        axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1),
        plot.title = element_text(size = 6, face = 'plain', margin = margin(0,0,0,0)),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA))
gg2$layers[[1]]$aes_params$size <- 0.1
gg1 / gg2
dev.off()

## Figure 4F----
# Enriched gene sets (Gene Ontology biological process) based on top contributing genes in NMF factor_6 and NMF factor_9
genesets <- msigdb_gsets(species = 'Mus musculus', category = 'C5', subcategory = 'BP') 

fb.factor.6.top50 <- rownames(fb_18_programs_diag_L1_0_0[["w"]][order(fb_18_programs_diag_L1_0_0[["w"]][, 6], decreasing = T), ]) [1:50]
fb.factor.9.top50 <- rownames(fb_18_programs_diag_L1_0_0[["w"]][order(fb_18_programs_diag_L1_0_0[["w"]][, 9], decreasing = T), ]) [1:50]

hyp.fb.factor.6.top50 <- hypeR(fb.factor.6.top50, genesets, test='hypergeometric', background=19405)
hyp.fb.factor.9.top50 <- hypeR(fb.factor.9.top50, genesets, test='hypergeometric', background=19405)

temp.nmf6.up <- hyp.fb.factor.6.top50[['data']][, c('label', 'fdr')]
temp.nmf6.up <- temp.nmf6.up %>% arrange(fdr) %>% slice(1:5)
temp.nmf6.up$neg_log10_adj_pval <- -log10(temp.nmf6.up$fdr)
temp.nmf6.up$label <- gsub("GO_", "", temp.nmf6.up$label)
temp.nmf6.up$label <- factor(temp.nmf6.up$label, levels = temp.nmf6.up$label[order(temp.nmf6.up$neg_log10_adj_pval)])
temp.nmf6.up$sign <- ifelse(temp.nmf6.up$neg_log10_adj_pval >= 0, "positive", "negative")

pdf('Figure_4F_top_factor6.pdf', width = 2.5, height = 1.1, bg = 'transparent')
ggplot(temp.nmf6.up, aes(x = label, y = neg_log10_adj_pval, fill = sign)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = 0, label = label, hjust = ifelse(neg_log10_adj_pval > 1, 0, 0)), size = 2) +
  scale_fill_manual(values = c('positive' = alpha('magenta', 0.5), 'negative' = alpha('black', 0.5))) +
  coord_flip() +
  labs(title = "factor_6 (GO: BP)", x = "", y = "-log10(adjusted P value)") +
  theme_classic() +
  ylim(0,36) +
  theme(axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 6, color = 'black'), axis.text.y = element_blank(),
        axis.title = element_text(size = 6), title = element_text(size = 6),
        plot.title = element_text(size = 6, face = 'plain', margin = margin(0,0,0,0)),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA)) +
  NoLegend()
dev.off()

temp.nmf9.up <- hyp.fb.factor.9.top50[['data']][, c('label', 'fdr')]
temp.nmf9.up <- temp.nmf9.up %>% arrange(fdr) %>% slice(1:5)
temp.nmf9.up$neg_log10_adj_pval <- -log10(temp.nmf9.up$fdr)
temp.nmf9.up$label <- gsub("GO_", "", temp.nmf9.up$label)
temp.nmf9.up$label <- factor(temp.nmf9.up$label, levels = temp.nmf9.up$label[order(temp.nmf9.up$neg_log10_adj_pval)])
temp.nmf9.up$sign <- ifelse(temp.nmf9.up$neg_log10_adj_pval >= 0, "positive", "negative")

pdf('Figure_4F_bottom_factor9.pdf', width = 2.5, height = 1.1, bg = 'transparent')
ggplot(temp.nmf9.up, aes(x = label, y = neg_log10_adj_pval, fill = sign)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = 0, label = label, hjust = ifelse(neg_log10_adj_pval > 1, 0, 0)), size = 2) +
  scale_fill_manual(values = c('positive' = alpha('magenta', 0.5), 'negative' = alpha('black', 0.5))) +
  coord_flip() +
  labs(title = "factor_9 (GO: BP)", x = "", y = "-log10(adjusted P value)") +
  theme_classic() +
  ylim(0,8) +
  theme(axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 6, color = 'black'), axis.text.y = element_blank(),
        axis.title = element_text(size = 6), title = element_text(size = 6),
        plot.title = element_text(size = 6, face = 'plain', margin = margin(0,0,0,0)),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA)) +
  NoLegend()
dev.off()

## Figure 4G----
# Kaplan-Meier survival plot of PAAD patients from TCGA excluding PNET (pancreatic neuroendocrine tumor) diagnosed patients categorized into upper and lower quartiles by the ratio of ‘AAV-gTgfbr2 up’ DEG signature to the fibroblast progenitor signature
tcga_paad_no_pnet <- fread('R_objects/public_datasets/PAADnoPNET_tcga_TPM_primary_annotated.tsv', sep='\t', header=T)
rownames(tcga_paad_no_pnet) <- tcga_paad_no_pnet$sample_name

signatures <- read_excel("Supplementary_Tables/Supplementary Table 4 signature genes.xlsx", sheet = 'signatures')
aav.tgfbr2.up.deg.sig <- signatures$`AAV-gTgfbr2 DEG signature`
progenitor.sig <- signatures$`progenitor signature (derived from Gao et al., 2024)`
progenitor.sig <- progenitor.sig[!is.na(progenitor.sig)]

# Subset data based on genes from signature
aav.tgfbr2.up.deg.sig_exp <- subset(tcga_paad_no_pnet, select = aav.tgfbr2.up.deg.sig)
aav.tgfbr2.up.deg.sig_exp_log <- as.data.frame(lapply(aav.tgfbr2.up.deg.sig_exp, function(x) log10(x + 1)))
# Calculate percentiles for each gene across samples
aav.tgfbr2.up.deg.sig_percentiles <- apply(aav.tgfbr2.up.deg.sig_exp_log, 2, function(x) ecdf(x)(x) * 100)
# Calculate the gene signature score for each sample (average of percentiles)
aav.tgfbr2.up.deg.sig_score_percentiles <- rowMeans(aav.tgfbr2.up.deg.sig_percentiles)
# Apply to original data table
tcga_paad_no_pnet$aav.tgfbr2.up.deg.sig_score_percentiles <- aav.tgfbr2.up.deg.sig_score_percentiles
# Survival code
tcga_paad_no_pnet$aav.tgfbr2.up.deg.sig_score_cat <- ifelse(tcga_paad_no_pnet$aav.tgfbr2.up.deg.sig_score_percentiles >= quantile(tcga_paad_no_pnet$aav.tgfbr2.up.deg.sig_score_percentiles, 0.75), "aav.tgfbr2.up.deg.sig_High",
                                    ifelse(tcga_paad_no_pnet$aav.tgfbr2.up.deg.sig_score_percentiles <= quantile(tcga_paad_no_pnet$aav.tgfbr2.up.deg.sig_score_percentiles, 0.25), "aav.tgfbr2.up.deg.sig_Low", "aav.tgfbr2.up.deg.sig_Intermediate"))
# Subset data based on genes from signature
progenitor.sig_exp <- subset(tcga_paad_no_pnet, select = progenitor.sig)
progenitor.sig_exp_log <- as.data.frame(lapply(progenitor.sig_exp, function(x) log10(x + 1)))
# Calculate percentiles for each gene across samples
progenitor.sig_percentiles <- apply(progenitor.sig_exp_log, 2, function(x) ecdf(x)(x) * 100)
# Calculate the gene signature score for each sample (average of percentiles)
progenitor.sig_score_percentiles <- rowMeans(progenitor.sig_percentiles)
# Apply to original data table
tcga_paad_no_pnet$progenitor.sig_score_percentiles <- progenitor.sig_score_percentiles
# Survival code
tcga_paad_no_pnet$progenitor.sig_score_cat <- ifelse(tcga_paad_no_pnet$progenitor.sig_score_percentiles >= quantile(tcga_paad_no_pnet$progenitor.sig_score_percentiles, 0.7), "progenitor.sig_High",
                                    ifelse(tcga_paad_no_pnet$progenitor.sig_score_percentiles <= quantile(tcga_paad_no_pnet$progenitor.sig_score_percentiles, 0.3), "progenitor.sig_Low", "progenitor.sig_Intermediate"))
tcga_paad_no_pnet$aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score <- tcga_paad_no_pnet$aav.tgfbr2.up.deg.sig_score_percentiles / tcga_paad_no_pnet$progenitor.sig_score_percentiles

tcga_paad_no_pnet$aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_cat <- ifelse(tcga_paad_no_pnet$aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score >= quantile(tcga_paad_no_pnet$aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score, 0.75), "aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_High",
                                                        ifelse(tcga_paad_no_pnet$aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score <= quantile(tcga_paad_no_pnet$aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score, 0.25), "aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_Low", "aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_Intermediate"))
table(tcga_paad_no_pnet$aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_cat)

ggplot(tcga_paad_no_pnet, aes(x=aav.tgfbr2.up.deg.sig_score_percentiles, y=progenitor.sig_score_percentiles, color = aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_cat)) + geom_point() + theme_bw()


Survival_param_aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_all <- c('aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_cat', 'OS.time', 'OS.status', 'indication', 'sex', 'age')
TCGA_aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_surv_all <- subset(tcga_paad_no_pnet, select = Survival_param_aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_all)
TCGA_aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_surv_all_sub <- subset(TCGA_aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_surv_all, aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_cat %in% c('aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_High', 'aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_Low'))

Survival_table <- survfit(Surv(OS.time, OS.status) ~ aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_cat, data = TCGA_aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_surv_all_sub)

# multivariate regression (i.e. sex & age)
do.cox.group <- function(df, group.order) {
  cox <- c()
  df$aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_cat <- factor(df$aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_cat, levels = group.order)
  if (length(unique(df$sex)) > 1) {
    cox <- coxph(Surv(OS.time, OS.status) ~ aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_cat + sex + age, data = df)
  } else {
    cox <- coxph(Surv(OS.time, OS.status) ~ aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_cat + age, data = df)
  }
  cox.summary <- summary(cox)
  nvar <- length(unique(df$aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_cat)) - 1
  HR <- round(cox.summary$conf.int[1:nvar, 1], 2)
  HR.lower <- round(cox.summary$conf.int[1:nvar, 3], 2)
  HR.upper <- round(cox.summary$conf.int[1:nvar, 4], 2)
  HR.range <- sprintf("(%.1f-%.1f)", HR.lower, HR.upper)
  coef <- cox.summary$coefficients[1:nvar, 1]
  coef.se <- cox.summary$coefficients[1:nvar, 3]
  Pval <- round(cox.summary$coefficients[1:nvar, 5], 4)
  aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_cat <- gsub("aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_cat", "", rownames(cox.summary$conf.int)[1:nvar])
  return(data.frame(aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_cat = aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_cat, HR = HR, HR.range = HR.range, coef = coef, coef.se = coef.se, Pval = Pval))
}

HR.pval.df <- do.cox.group(TCGA_aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_surv_all_sub, group.order = c("aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_Low", "aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_High"))
annotation_text <- paste("HR:", HR.pval.df$HR, '\n', "p =", HR.pval.df$Pval)

sorted_category_names <- sort(TCGA_aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_surv_all_sub$aav.tgfbr2.up.deg.sig_to_progenitor.sig_ratio_score_cat) #order the names of the category ('sig_High' & 'sig_Low') in alphabetical order, so the legend is correctly annotated!
legend_labels <- paste0(unique(sorted_category_names), " (n=", Survival_table$n, ")")
plot <- ggsurvplot(Survival_table, conf.int = FALSE, surv.size = 1.2, xlab = 'Days', ylab = 'Survival', legend.title = "", legend.labs = legend_labels,
                   ylim = c(0, NA), xlim = c(0, 2500), break.x.by = 1000, palette = c('red', 'black'), pval = F, size = 0.25,
                   censor.shape = 124, censor.size = 1) +
  ggtitle(paste("TCGA PAAD: AAV-gTgfbr2 up DEG / progenitor signature"))
plot$plot <- plot$plot + guides(color = guide_legend(ncol = 1)) + # make legend in vertical alignment
  annotate("text", x = Inf, y = Inf, label = annotation_text, hjust = 1.1, vjust = 1.1, size = 6*0.3527) # add HR and pval annotation
plot

# EXPORT FIGURE
pdf('Figure_4G.pdf', width = 1.6, height = 2.2, bg = 'transparent')
plot$plot + theme(axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1),
                  axis.text.x = element_text(size = 6, color = 'black'), axis.text.y = element_text(size = 6, color = 'black'),
                  axis.title.x = element_text(size = 6), axis.title.y = element_text(size = 6), title = element_text(size = 6),
                  legend.text = element_text(size = 6), 
                  plot.title = element_text(size = 6, face = 'plain', margin = margin(0,0,0,0)),
                  legend.background = element_rect(fill = 'transparent', color = NA),
                  plot.background = element_rect(fill = "transparent", color = NA),
                  panel.background = element_rect(fill = "transparent", color = NA))
dev.off()
