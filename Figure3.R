# MANUSCRIPT TITLE: Discovery of a Linked Constellation of Gene Expression Revealed by Local Editing of Fibroblasts in Tumors
# Companion code: FIGURE 3
# Title: Local Tgfbr2 knockout in cancer-associated fibroblasts reveals emergence of unique Col18a1hi CAF cell state
# Date: 2025-11-06

library(Seurat)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(rstatix)
library(DescTools)
library(patchwork)
library(readxl)
library(dplyr)

## Figure 3B----
# Uniform manifold approximation and projection (UMAP) of cancer-associated fibroblast (CAF) subsets
fbaav_0.4 <- readRDS("R_objects/Seurat_objects/fbaav_0.4.rds") # fibroblast Seurat object
png("Figure_3B.png", width = 800, height = 800, res = 300, bg='transparent')
DimPlot(fbaav_0.4, reduction = 'umap', label = T, label.size = 3, repel = F, raster = F, pt.size = 0.1) + NoLegend()  &
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'), axis.line = element_line(linewidth = 0),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA))
dev.off()

# Figure 3C----
# Marker gene expression of identified CAF subsets
my_levels <- c('homeostatic', 'intermediate', 'infl_CAF', 'myCAF_1', 'myCAF_2', 'myCAF_3', 'myCAF_prolif', 'Col18a1')
fbaav_0.4[["clustering_1"]] <- factor(x = fbaav_0.4@meta.data$clustering_1, levels = my_levels)

Idents(fbaav_0.4) <- 'clustering_1'
pdf("Figure_3C.pdf", width = 3.2, height = 1.5)
DotPlot(fbaav_0.4,
        features = c('Pi16', 'Ly6c1', 'Ly6a', 'Il6', 'Lif', 'Saa3', 'Acta2', 'Lrrc15', 'Ncam1', 'Slit1', 'Fabp5', 'Wif1', 'Lamc3', 'Mki67', 'Col18a1', 'Lratd1'),
        cols = 'RdBu', assay = "RNA", dot.scale = 2) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, size = 6, face = 'oblique'), axis.text.y=element_text(size = 6), text = element_text(size = 6)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.ticks = element_line(linewidth = 0.1)) +
  theme(axis.line = element_line(linewidth = 0.1))
dev.off()

# Figure 3D----
# UMAP of CAFs highlighted by AAV-gRNA condition
png("Figure_3D.png", width = 5000, height = 5000, res = 300, bg = 'transparent')
p6i <- DimPlot(fbaav_0.4, pt.size = 2, cells.highlight = colnames(fbaav_0.4)[fbaav_0.4$aav=='gTrac'], sizes.highlight = 2) + NoLegend() & theme(title = element_text(size = 12.25), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), axis.line = element_line(linewidth = 0), panel.background = element_rect(fill = 'transparent', color = NA), plot.background = element_rect(fill = 'transparent', color = NA))
p6ii <- DimPlot(fbaav_0.4, pt.size = 2, cells.highlight = colnames(fbaav_0.4)[fbaav_0.4$aav=='gOsmr'], sizes.highlight = 2) + NoLegend() & theme(title = element_text(size = 12.25), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), axis.line = element_line(linewidth = 0), panel.background = element_rect(fill = 'transparent', color = NA), plot.background = element_rect(fill = 'transparent', color = NA))
p6iii <- DimPlot(fbaav_0.4, pt.size = 2, cells.highlight = colnames(fbaav_0.4)[fbaav_0.4$aav=='gTgfbr2'], sizes.highlight = 2) + NoLegend() & theme(title = element_text(size = 12.25), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), axis.line = element_line(linewidth = 0), panel.background = element_rect(fill = 'transparent', color = NA), plot.background = element_rect(fill = 'transparent', color = NA))
p6iv <- DimPlot(fbaav_0.4, pt.size = 2, cells.highlight = colnames(fbaav_0.4)[fbaav_0.4$aav=='gIl1r1'], sizes.highlight = 2) + NoLegend() & theme(title = element_text(size = 12.25), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), axis.line = element_line(linewidth = 0), panel.background = element_rect(fill = 'transparent', color = NA), plot.background = element_rect(fill = 'transparent', color = NA))
p6 <- ((p6i | p6ii) / (p6iii | p6iv))
p6
dev.off()

# Figure 3E----
# Bar charts of frequency of CAF subsets (scRNA-Seq) as a fraction of all fibroblasts by AAV-gRNA condition
cell_freq <- table(fbaav_0.4@meta.data$orig.ident, fbaav_0.4@meta.data$clustering_1)

cell_freq <- cell_freq/rowSums(cell_freq)
cell_freq <- reshape2::melt(cell_freq)

#add AAV annotation
cell_freq$AAV <- NA
cell_freq$AAV[which(cell_freq$Var1 %in% c('gTrac_A1', 'gTrac_A2', 'gTrac_B1', 'gTrac_B2'))] <- 'gTrac'
cell_freq$AAV[which(cell_freq$Var1 %in% c('gOsmr_A3', 'gOsmr_A4', 'gOsmr_B3', 'gOsmr_B4'))] <- 'gOsmr'
cell_freq$AAV[which(cell_freq$Var1 %in% c('gTgfbr2_A5', 'gTgfbr2_A6', 'gTgfbr2_B5', 'gTgfbr2_B6'))] <- 'gTgfbr2'
cell_freq$AAV[which(cell_freq$Var1 %in% c('gIl1r1_A7', 'gIl1r1_A8', 'gIl1r1_B7', 'gIl1r1_B8'))] <- 'gIl1r1'
#add replicate annotation
cell_freq$replicate <- NA
cell_freq$replicate[which(cell_freq$Var1 %in% c('gTrac_A1', 'gOsmr_A3', 'gTgfbr2_A5', 'gIl1r1_A7'))] <- '1'
cell_freq$replicate[which(cell_freq$Var1 %in% c('gTrac_A2', 'gOsmr_A4', 'gTgfbr2_A6', 'gIl1r1_A8'))] <- '2'
cell_freq$replicate[which(cell_freq$Var1 %in% c('gTrac_B1', 'gOsmr_B3', 'gTgfbr2_B5', 'gIl1r1_B7'))] <- '3'
cell_freq$replicate[which(cell_freq$Var1 %in% c('gTrac_B2', 'gOsmr_B4', 'gTgfbr2_B6', 'gIl1r1_B8'))] <- '4'
color_group <- c('lightgrey', '#1A85FF', '#D41159', '#40B0A6')
cell_freq$Var2 <- as.factor(cell_freq$Var2)
cell_freq$AAV <- factor(cell_freq$AAV, levels = c('gTrac', 'gOsmr', 'gTgfbr2', 'gIl1r1'))

#cast data.frame for stat analysis
df <- reshape2::dcast(cell_freq, Var1 ~ Var2, value = 'Freq')
#add AAV annotation
df$AAV <- NA
df$AAV[which(df$Var1 %in% c('gTrac_A1', 'gTrac_A2', 'gTrac_B1', 'gTrac_B2'))] <- 'gTrac'
df$AAV[which(df$Var1 %in% c('gOsmr_A3', 'gOsmr_A4', 'gOsmr_B3', 'gOsmr_B4'))] <- 'gOsmr'
df$AAV[which(df$Var1 %in% c('gTgfbr2_A5', 'gTgfbr2_A6', 'gTgfbr2_B5', 'gTgfbr2_B6'))] <- 'gTgfbr2'
df$AAV[which(df$Var1 %in% c('gIl1r1_A7', 'gIl1r1_A8', 'gIl1r1_B7', 'gIl1r1_B8'))] <- 'gIl1r1'
df$AAV <- factor(df$AAV, levels = c('gTrac', 'gOsmr', 'gTgfbr2', 'gIl1r1'))
#add replicate annotation
df$replicate <- NA
df$replicate[which(df$Var1 %in% c('gTrac_A1', 'gOsmr_A3', 'gTgfbr2_A5', 'gIl1r1_A7'))] <- '1'
df$replicate[which(df$Var1 %in% c('gTrac_A2', 'gOsmr_A4', 'gTgfbr2_A6', 'gIl1r1_A8'))] <- '2'
df$replicate[which(df$Var1 %in% c('gTrac_B1', 'gOsmr_B3', 'gTgfbr2_B5', 'gIl1r1_B7'))] <- '3'
df$replicate[which(df$Var1 %in% c('gTrac_B2', 'gOsmr_B4', 'gTgfbr2_B6', 'gIl1r1_B8'))] <- '4'

# Function to perform Dunnett's test with gTrac as reference and extract p-values
perform_dunnett_test <- function(column_name) {
  formula <- as.formula(paste(column_name, "~ AAV"))
  test_result <- DunnettTest(formula, data = df, control = 'gTrac')
  return(test_result)
}

# Columns to analyze
columns_to_analyze <- names(df)[2:(ncol(df)-2)]

# Perform Dunnett's test for each column and store the results
dunnett_test_results <- lapply(columns_to_analyze, perform_dunnett_test)

# Name the results for clarity
names(dunnett_test_results) <- columns_to_analyze

pval.0 <- data.frame(dunnett_test_results$myCAF_1$gTrac)
pval.0$signif <- NA
pval.0$signif[which(pval.0$pval <= 0.05)] <- '*'
pval.0$Var2 <- 'myCAF_1'
pval.0$Var2 <- as.factor(pval.0$Var2)
pval.0$group1 <- sub(".*-", "", rownames(pval.0))
pval.0$group2 <- sub("-.*", "", rownames(pval.0))

pval.1 <- data.frame(dunnett_test_results$infl_CAF$gTrac)
pval.1$signif <- NA
pval.1$signif[which(pval.1$pval <= 0.05)] <- '*'
pval.1$Var2 <- 'infl_CAF'
pval.1$Var2 <- as.factor(pval.1$Var2)
pval.1$group1 <- sub(".*-", "", rownames(pval.1))
pval.1$group2 <- sub("-.*", "", rownames(pval.1))

pval.2 <- data.frame(dunnett_test_results$myCAF_prolif$gTrac)
pval.2$signif <- NA
pval.2$signif[which(pval.2$pval <= 0.05)] <- '*'
pval.2$Var2 <- 'myCAF_prolif'
pval.2$Var2 <- as.factor(pval.2$Var2)
pval.2$group1 <- sub(".*-", "", rownames(pval.2))
pval.2$group2 <- sub("-.*", "", rownames(pval.2))

pval.3 <- data.frame(dunnett_test_results$Col18a1$gTrac)
pval.3$signif <- NA
pval.3$signif[which(pval.3$pval <= 0.01)] <- '**'
pval.3$Var2 <- 'Col18a1'
pval.3$Var2 <- as.factor(pval.3$Var2)
pval.3$group1 <- sub(".*-", "", rownames(pval.3))
pval.3$group2 <- sub("-.*", "", rownames(pval.3))

pval.4 <- data.frame(dunnett_test_results$myCAF_2$gTrac)
pval.4$signif <- NA
pval.4$signif[which(pval.4$pval <= 0.05)] <- '*'
pval.4$Var2 <- 'myCAF_2'
pval.4$Var2 <- as.factor(pval.4$Var2)
pval.4$group1 <- sub(".*-", "", rownames(pval.4))
pval.4$group2 <- sub("-.*", "", rownames(pval.4))

pval.5 <- data.frame(dunnett_test_results$homeostatic$gTrac)
pval.5$signif <- NA
pval.5$signif[which(pval.5$pval <= 0.05)] <- '*'
pval.5$Var2 <- 'homeostatic'
pval.5$Var2 <- as.factor(pval.5$Var2)
pval.5$group1 <- sub(".*-", "", rownames(pval.5))
pval.5$group2 <- sub("-.*", "", rownames(pval.5))

pval.6 <- data.frame(dunnett_test_results$intermediate$gTrac)
pval.6$signif <- NA
pval.6$signif[which(pval.6$pval <= 0.05)] <- '*'
pval.6$Var2 <- 'intermediate'
pval.6$Var2 <- as.factor(pval.6$Var2)
pval.6$group1 <- sub(".*-", "", rownames(pval.6))
pval.6$group2 <- sub("-.*", "", rownames(pval.6))

pval.7 <- data.frame(dunnett_test_results$myCAF_3$gTrac)
pval.7$signif <- NA
pval.7$signif[which(pval.7$pval <= 0.05)] <- '*'
pval.7$Var2 <- 'myCAF_3'
pval.7$Var2 <- as.factor(pval.7$Var2)
pval.7$group1 <- sub(".*-", "", rownames(pval.7))
pval.7$group2 <- sub("-.*", "", rownames(pval.7))

dunnett.pvals <- rbind(pval.5, pval.6, pval.1, pval.0, pval.4, pval.7, pval.2, pval.3)
df$AAV <- factor(df$AAV, levels = c('gOsmr', 'gTgfbr2', 'gIl1r1'))
cell_freq$Var2 <- factor(cell_freq$Var2, levels = c('homeostatic', 'intermediate', 'infl_CAF', 'myCAF_1', 'myCAF_2', 'myCAF_3', 'myCAF_prolif', 'Col18a1'))

dunnett.pvals$AAV <- dunnett.pvals$group2
dunnett.pvals$y.position <- 0.6

pdf("Figure_3E.pdf", width = 6, height = 3)
ggplot(data = cell_freq, mapping = aes(x = AAV,
                                       y = value,
                                       group = AAV)) +
  stat_summary(aes(fill = AAV),
               fun.data = mean_se, geom = 'col', color = 'black', position = position_dodge(1)) +
  stat_summary(fun.data = mean_se, geom = 'errorbar', width = 0.2, position = position_dodge(1)) +
  geom_jitter(aes(shape = replicate),
              position = position_jitterdodge(dodge.width = 1, jitter.width = 0.5),
              size = 1) +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  scale_fill_manual(values = color_group) +
  labs(y = "fraction of CAFs") +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size =14),
        legend.title = element_text(size = 14), legend.text = element_text(size = 14)) +
  stat_pvalue_manual(dunnett.pvals, label = "signif", tip.length = 0, step.group.by = 'Var2', step.increase = 0.1, hide.ns = F) +
  facet_wrap(~Var2, ncol = 4)
dev.off()

# Figure 3F----
# Similarity of DEGs between mouse CAF subsets (columns)
fb_deg_markers <- read_excel("Supplementary_Tables/Supplementary Table 1 DEG scRNAseq.xlsx", sheet = 'fibroblasts')
aav.fb.deg <- split(fb_deg_markers, fb_deg_markers$clustering_1)
# mouse fibroblasts in perturbed states from Buechler et al., 2021
buechler.deg.ms.ps <- readRDS("R_objects/public_datasets/buechler.deg.ms.ps.rds")

### HYPERGEOMETRIC TEST----
# Define background:
N <- 19405 # number of genes captured by FLEX kit
# read in lists
listA <- aav.fb.deg
listB <- buechler.deg.ms.ps

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
    setB <- listB[[b]] %>% filter(p_val_adj < 0.05 & avg_logFC > 0.25) %>% pull(Gene)
    
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
y_levels <- c('Pi16', 'Col15a1', 'Ccl19', 'Npnt', 'Cxcl12', 'Comp', 'Hhip', 'Adamdec1', 'Lrrc15', 'Cxcl5') # Buechler

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

pdf('Figure_3F.pdf', width=2, height=1.3)
ggplot(data = melt_hgt, aes(x=Var1, y=Var2, fill=(-log10(P)))) + 
  geom_tile() +
  scale_fill_viridis_c(limits=c(0, 100), oob = scales::squish) +
  ggtitle('Buechler 2021, mouse perturbed state') +
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

# Figure 3G----
# Similarity of DEGs between mouse CAF subsets (columns)
# and mouse PDAC CAF subsets from Elyada et al, 2019
elyada.ms.caf.markers <- read_excel("R_objects/public_datasets/DEGs_Elyada2019_mouse_KPC.xlsx")
elyada.ms.caf.markers$...1 <- NULL
elyada.ms.caf.deg <- split(elyada.ms.caf.markers, elyada.ms.caf.markers$cluster)

### HYPERGEOMETRIC TEST----
# Define background:
N <- 19405 # number of genes captured by FLEX kit
# read in lists
listA <- aav.fb.deg
listB <- elyada.ms.caf.deg

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
    setB <- listB[[b]] %>% pull(Associated.Gene.Name)
    
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

pdf('Figure_3G.pdf', width=2, height=0.8)
ggplot(data = melt_hgt, aes(x=Var1, y=Var2, fill=(-log10(P)))) + 
  geom_tile() +
  scale_fill_viridis_c(limits=c(0, 40), oob = scales::squish) +
  ggtitle('Elyada 2019, mouse PDAC') +
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

# Figure 3H----
# UMAP visualization of select genes in CAFs
png("Figure_3H.png", width = 3000, height = 1620, res = 300)
fp <- FeaturePlot(fbaav_0.4, c('Ly6c1', 'Acta2', 'Col18a1', 'Ly6a', 'Ncam1', 'Itga1'), order = T,
                  cols = c('lightgrey', 'red'), pt.size = 0.1)
fp <- lapply(fp, function(fp) fp +
               theme(plot.title = element_text(size = 24, face = 'italic', margin=margin(0,0,0,0)), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
                     legend.text = element_text(size = 24), legend.key.size = unit(1.5, 'lines'), legend.key.width = unit(0.1, 'lines'), legend.box.margin = margin(0, 0, 0, 0, 'mm'), legend.spacing = unit(0, 'mm'),
                     #legend.position = c(1, 0.5),
                     plot.margin = margin(0, 0, 0, 0, unit = 'pt'), axis.line = element_blank()))
fp <- wrap_plots(fp, ncol = 3)
fp
dev.off()

