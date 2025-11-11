# MANUSCRIPT TITLE: Discovery of a Linked Constellation of Gene Expression Revealed by Local Editing of Fibroblasts in Tumors
# Companion code: FIGURE 6
# Title: Combinatorial knockout shows that emergence of Tgfbr2 knockout-induced Col18a1hi CAFs is dependent on TNFR1 and canonical Wnt signaling
# Date: 2025-11-06

library(khroma)
library(Seurat)
library(ggplot2)
library(readxl)
library(dplyr)
library(hypeR)
library(msigdbr)
library(CellChat)
library(decoupleR)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)

# Figure 6A----
# Chord diagrams of ligand-receptor pairs between TAM.Ndrg1 cells (senders) and Col18a1 CAFs (receivers) or neutrophil.4 subset (senders) and Col18a1 CAFs (receivers)
load("R_objects/CellChat_object/cellchat.RData")

source_cell <- 'TAM.Ndrg1'
target_cell <- "Col18a1"
cell.colors = c('TAM.Ndrg1' = 'green4', 'Col18a1' = 'black')
df_interaction <- subsetCommunication(cellchat, sources.use = source_cell, targets.use = target_cell)
df_ranked <- df_interaction[order(df_interaction$prob, decreasing = TRUE), ]
df_top15 <- df_ranked %>% arrange(desc(prob)) %>% slice_head(n = 15)
df_top15_lr <- df_top15 %>% select(interaction_name)
netVisual_chord_gene(cellchat, sources.use = source_cell, targets.use = target_cell, pairLR.use = df_top15_lr, color.use = cell.colors, small.gap = 3, big.gap = 10)

source_cell <- 'neutrophil.4'
target_cell <- "Col18a1"
cell.colors = c('neutrophil.4' = 'magenta3', 'Col18a1' = 'black')
df_interaction <- subsetCommunication(cellchat, sources.use = source_cell, targets.use = target_cell)
df_ranked <- df_interaction[order(df_interaction$prob, decreasing = TRUE), ]
df_top9 <- df_ranked %>% arrange(desc(prob)) %>% slice_head(n = 9)
df_top9_lr <- df_top9 %>% select(interaction_name)
netVisual_chord_gene(cellchat, sources.use = source_cell, targets.use = target_cell, pairLR.use = df_top9_lr, color.use = cell.colors, small.gap = 3, big.gap = 10)

# Figure 6B----
# Proportion of MonoMac subsets (scRNA-seq) in each AAV-gRNA group.
# Proportion of neutrophil subsets (scRNA-seq) in each AAV-gRNA group.
mmaav_0.4 <- readRDS("R_objects/Seurat_objects/mmaav_0.4.rds") # myeloid cells
neuts_0.4 <- readRDS("R_objects/Seurat_objects/neuts_0.4.rds") # neutrophils

discrete_rainbow <- color('discrete rainbow')
col_pal <- discrete_rainbow(11)

pdf('Figure_6B_myeloid.pdf', width = 1.8, height = 1.3, bg = 'transparent')
#Idents(neuts_0.4) <- 'annotation_1'
ggplot(mmaav_0.4@meta.data, aes(x=aav, fill = annotation_1)) +
  geom_bar(aes(fill=annotation_1),position='fill') +
  labs(title = 'MonoMac subset distribution', x='', y = 'Proportion (%)') +
  theme_classic() +
  scale_fill_manual(values = discrete_rainbow(11)) +
  theme(text = element_text(size = 6), axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5, size = 6, color = 'black'), axis.text.y=element_text(size = 6, color = 'black'), axis.title.x = element_blank(), axis.title.y = element_text(size = 6), plot.title = element_text(size = 6),
        axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1), plot.margin = margin(0,0,0,0),
        legend.key.size = unit(0.4, 'lines'), legend.key.width = unit(0.3, 'lines'), legend.title = element_blank(), legend.margin = margin(0,0,0,0), legend.text = element_text(size = 6), legend.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA))
dev.off()

pdf('Figure_6B_neutrophils.pdf', width = 1.6, height = 1.3, bg = 'transparent')
#Idents(neuts_0.4) <- 'annotation_1'
ggplot(neuts_0.4@meta.data, aes(x=aav, fill = annotation_1)) +
  geom_bar(aes(fill=annotation_1),position='fill') +
  labs(title = 'neutrophil subset distribution', x='', y = 'Proportion (%)') +
  theme_classic() +
  scale_fill_manual(values = discrete_rainbow(11)) +
  theme(text = element_text(size = 6), axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5, size = 6, color = 'black'), axis.text.y=element_text(size = 6, color = 'black'), axis.title.x = element_blank(), axis.title.y = element_text(size = 6), plot.title = element_text(size = 6),
        axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1), plot.margin = margin(0,0,0,0),
        legend.key.size = unit(0.4, 'lines'), legend.key.width = unit(0.3, 'lines'), legend.title = element_blank(), legend.margin = margin(0,0,0,0), legend.text = element_text(size = 6), legend.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA))
dev.off()

# Figure 6C----
# Enriched gene sets (Gene Ontology hallmark) based on DEGs from TAM.Ndrg1 cluster
# Enriched gene sets (Gene Ontology hallmark) based on DEGs from neutrophil.4 cluster
genesets <- msigdb_gsets(species = 'Mus musculus', category = 'H') # HALLMARK

mmaav_0.4_markers <- read_excel("Supplementary_Tables/Supplementary Table 1 DEG scRNAseq.xlsx", sheet = 'MonoMac')
mm.deg.ndrg1.signif <- mmaav_0.4_markers %>% filter(cluster == 'TAM.Ndrg1') %>% filter(p_val_adj < 0.05 & avg_log2FC > sqrt(0.25)) %>% pull(gene)

hyp.mm.deg.ndrg1.signif <- hypeR(mm.deg.ndrg1.signif, genesets, test='hypergeometric', background=19405)

temp.mm.deg.ndrg1.up <- hyp.mm.deg.ndrg1.signif[['data']][, c('label', 'fdr')]
temp.mm.deg.ndrg1.up <- temp.mm.deg.ndrg1.up %>% arrange(fdr) %>% slice_head(n = 5)
temp.mm.deg.ndrg1.up$neg_log10_adj_pval <- -log10(temp.mm.deg.ndrg1.up$fdr)
temp.mm.deg.ndrg1.up$label <- gsub("HALLMARK_", "", temp.mm.deg.ndrg1.up$label)
temp.mm.deg.ndrg1.up$label <- factor(temp.mm.deg.ndrg1.up$label, levels = temp.mm.deg.ndrg1.up$label[order(temp.mm.deg.ndrg1.up$neg_log10_adj_pval)])
temp.mm.deg.ndrg1.up$sign <- ifelse(temp.mm.deg.ndrg1.up$neg_log10_adj_pval >= 0, "positive", "negative")

pdf('Figure_6C_TAM.Ndrg1.pdf', width = 1.5, height = 1.1, bg = 'transparent')
ggplot(temp.mm.deg.ndrg1.up, aes(x = label, y = neg_log10_adj_pval, fill = sign)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = 0, label = label, hjust = ifelse(neg_log10_adj_pval > 1, 0, 0)), size = 2) +
  scale_fill_manual(values = c('positive' = alpha('black', 0.5), 'negative' = alpha('black', 0.5))) +
  coord_flip() +
  labs(title = "DEGs TAM.Ndrg1 (hallmark)") +
  theme_classic() +
  xlab('') +
  ylab(expression(-log[10](Padj))) +
  #ylim(0,20) +
  theme(axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 6, color = 'black'), axis.text.y = element_blank(),
        axis.title = element_text(size = 6), title = element_text(size = 6),
        plot.title = element_text(size = 6, face = 'plain', margin = margin(0,0,0,0)),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA)) +
  NoLegend()
dev.off()

neuts_0.4_markers <- read_excel("Supplementary_Tables/Supplementary Table 1 DEG scRNAseq.xlsx", sheet = 'Neutrophils')
neuts.deg.4.signif <- neuts_0.4_markers %>% filter(annotation_1 == 'neutrophil.4') %>% filter(p_val_adj < 0.05 & avg_log2FC > sqrt(0.25)) %>% pull(gene)

hyp.neuts.deg.4.signif <- hypeR(neuts.deg.4.signif, genesets, test='hypergeometric', background=19405)

temp.neuts.deg.4.up <- hyp.neuts.deg.4.signif[['data']][, c('label', 'fdr')]
temp.neuts.deg.4.up <- temp.neuts.deg.4.up %>% arrange(fdr) %>% slice_head(n = 5)
temp.neuts.deg.4.up$neg_log10_adj_pval <- -log10(temp.neuts.deg.4.up$fdr)
temp.neuts.deg.4.up$label <- gsub("HALLMARK_", "", temp.neuts.deg.4.up$label)
temp.neuts.deg.4.up$label <- factor(temp.neuts.deg.4.up$label, levels = temp.neuts.deg.4.up$label[order(temp.neuts.deg.4.up$neg_log10_adj_pval)])
temp.neuts.deg.4.up$sign <- ifelse(temp.neuts.deg.4.up$neg_log10_adj_pval >= 0, "positive", "negative")

pdf('Figure_6C_neutrophil.4.pdf', width = 1.5, height = 1.1, bg = 'transparent')
ggplot(temp.neuts.deg.4.up, aes(x = label, y = neg_log10_adj_pval, fill = sign)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = 0, label = label, hjust = ifelse(neg_log10_adj_pval > 1, 0, 0)), size = 2) +
  scale_fill_manual(values = c('positive' = alpha('black', 0.5), 'negative' = alpha('black', 0.5))) +
  coord_flip() +
  labs(title = "DEGs neutrophil.4 (hallmark)") +
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

# Figure 6D----
# Heatmap of average activity score calculated by PROGENy of osteoblast-associated transcription factors in CAF subsets
fbaav_0.4 <- readRDS("R_objects/Seurat_objects/fbaav_0.4.rds")
data <- fbaav_0.4

tf.net <- get_collectri(organism = 'mouse', split_complexes = F)
## if 'mouse' doesn't work, CLEAR CACHE:
#omnipath_cache_clean_db()

data@assays$RNA@layers$data@Dimnames[[1]] <- rownames(data@assays[["RNA"]]@features)
data@assays$RNA@layers$data@Dimnames[[2]] <- rownames(data@assays[["RNA"]]@cells)
mat <- as.matrix(data@assays$RNA@layers$data)

tf.acts <- run_ulm(mat=mat, net=tf.net, .source='source', .target='target',
                   .mor='mor', minsize = 5)

data[['tfsulm']] <- tf.acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

DefaultAssay(object = data) <- "tfsulm"
data <- ScaleData(data)
data@assays$tfsulm@data <- data@assays$tfsulm@scale.data

osteoblast_TFs <- c('Ctnnb1', 'Atf4', 'Cebpb', 'Cebpd', 'Dlx3', 'Dlx5', 'Ets1', 'Hes1', 'Hey1', 'Lef1', 'Men1', 'Msx2', 'Sp7', 'Runx1', 'Runx2', 'Runx3', 'Sall4', 'Satb2') # from https://www.nature.com/articles/s41421-024-00689-6/tables/1

DefaultAssay(object = data) <- "tfsulm"
Idents(data) <- 'clustering_1'

df <- t(as.matrix(data@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(data)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

tfs <- osteoblast_TFs

top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

col_fun <- colorRamp2(c(-1, 0, 1), c("darkblue", "white", 'darkred'))

pdf('Figure_6D.pdf', width = 3.2, height = 4.2, bg = 'transparent')
Heatmap(t(top_acts_mat),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_dend_reorder = TRUE,
        column_dend_reorder = TRUE,
        row_dend_gp = gpar(lwd = 1),
        column_dend_gp = gpar(lwd = 1),
        col = col_fun,
        heatmap_legend_param = list(title='z-score', title_gp = gpar(fontsize = 6, fontface = 'plain'), legend_width = unit(1, 'cm'), legend_height = unit(0.5, 'cm')),
        clustering_distance_rows = 'pearson',
        clustering_distance_columns = 'minkowski')
dev.off()
