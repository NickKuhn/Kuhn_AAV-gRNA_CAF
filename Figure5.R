# MANUSCRIPT TITLE: Discovery of a Linked Constellation of Gene Expression Revealed by Local Editing of Fibroblasts in Tumors
# Companion code: FIGURE 5
# Title: Tgfbr2 knockout-induced Col18a1hi CAFs recruit Siglec-Fhi neutrophils via high Cxcl5 expression
# Date: 2025-11-06

library(corrplot)
library(Seurat)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(DescTools)
library(RColorBrewer)

# Figure 5A----
# Heatmap showing the pairwise Pearson correlation of cell type frequencies among AAV-gRNA samples
fbaav_0.4 <- readRDS("R_objects/Seurat_objects/fbaav_0.4.rds") # fibroblasts
mmaav_0.4 <- readRDS("R_objects/Seurat_objects/mmaav_0.4.rds") # myeloid cells
tnk_clean <- readRDS("R_objects/Seurat_objects/tnk_clean_1.rds") # T and NK cells
dcs <- readRDS("R_objects/Seurat_objects/dcs_0.2.RData") # dendritic cells
neuts_0.4 <- readRDS("R_objects/Seurat_objects/neuts_0.4.rds") # neutrophils

cell_freq_fb <- table(fbaav_0.4@meta.data$orig.ident, fbaav_0.4@meta.data$clustering_1)
cell_freq_mm <- table(mmaav_0.4@meta.data$orig.ident, mmaav_0.4@meta.data$annotation_1)
cell_freq_tnk <- table(tnk_clean@meta.data$orig.ident, tnk_clean@meta.data$annotation_1)
cell_freq_dc <- table(dcs@meta.data$orig.ident, dcs@meta.data$clustering_1)
cell_freq_neuts <- table(neuts_0.4@meta.data$orig.ident, neuts_0.4@meta.data$annotation_1)

cell_freq_fb <- cell_freq_fb/rowSums(cell_freq_fb)
cell_freq_fb <- reshape2::melt(cell_freq_fb)
cell_freq_fb$Var2 <- paste0('FB.', cell_freq_fb$Var2)
cell_freq_mm <- cell_freq_mm/rowSums(cell_freq_mm)
cell_freq_mm <- reshape2::melt(cell_freq_mm)
cell_freq_tnk <- cell_freq_tnk/rowSums(cell_freq_tnk)
cell_freq_tnk <- reshape2::melt(cell_freq_tnk)
cell_freq_dc <- cell_freq_dc/rowSums(cell_freq_dc)
cell_freq_dc <- reshape2::melt(cell_freq_dc)
cell_freq_neuts <- cell_freq_neuts/rowSums(cell_freq_neuts)
cell_freq_neuts <- reshape2::melt(cell_freq_neuts)

mat.fb <- acast(cell_freq_fb, Var2 ~ Var1, value.var = 'value')
mat.mm <- acast(cell_freq_mm, Var2 ~ Var1, value.var = 'value')
mat.tnk <- acast(cell_freq_tnk, Var2 ~ Var1, value.var = 'value')
mat.dc <- acast(cell_freq_dc, Var2 ~ Var1, value.var = 'value')
mat.neuts <- acast(cell_freq_neuts, Var2 ~ Var1, value.var = 'value')

big.mat <- rbind(mat.fb, mat.mm, mat.tnk, mat.dc, mat.neuts) # no tumor, no endo

df_transposed <- t(big.mat)
correlation_matrix <- cor(df_transposed, method = "pearson")

correlation_df <- as.data.frame(correlation_matrix)
testRes = cor.mtest(correlation_df)
correlation_df <- as.matrix(correlation_df)

pdf('Figure_5A.pdf', width = 4.5, height = 4.5, bg = 'transparent')
corrplot(correlation_df, p.mat = testRes$p, sig.level = c(0.05),
         pch.cex = 0.5, tl.cex = 0.5, cl.cex = 0.5, insig = 'label_sig', method = 'circle', order = 'hclust', diag = T,
         type = 'full', col = rev(brewer.pal(10 ,"RdBu")), tl.col = 'black', col.lim = c(-1,1), addgrid.col = NA,
         mar = c(0,0,0,0))
dev.off()

# Figure 5B----
# Bar charts of frequency of neutrophil subsets (scRNA-Seq) as a fraction of all neutrophils by AAV-gRNA condition
cell_freq <- table(neuts_0.4@meta.data$orig.ident, neuts_0.4@meta.data$annotation_1)

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
df <- reshape2::dcast(cell_freq, Var1 ~ Var2, value = 'Freq') #the to-be-dcasted data.frame needs to be melted
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

pval.1 <- data.frame(dunnett_test_results$neutrophil.1$gTrac)
pval.1$signif <- NA
pval.1$signif[which(pval.1$pval <= 0.05)] <- '*'
pval.1$Var2 <- 'neutrophil.1'
pval.1$Var2 <- as.factor(pval.1$Var2)
pval.1$group1 <- sub(".*-", "", rownames(pval.1))
pval.1$group2 <- sub("-.*", "", rownames(pval.1))

pval.2 <- data.frame(dunnett_test_results$neutrophil.2$gTrac)
pval.2$signif <- NA
pval.2$signif[which(pval.2$pval <= 0.05)] <- '*'
pval.2$Var2 <- 'neutrophil.2'
pval.2$Var2 <- as.factor(pval.2$Var2)
pval.2$group1 <- sub(".*-", "", rownames(pval.2))
pval.2$group2 <- sub("-.*", "", rownames(pval.2))

pval.3 <- data.frame(dunnett_test_results$neutrophil.3$gTrac)
pval.3$signif <- NA
pval.3$signif[which(pval.3$pval <= 0.05)] <- '*'
pval.3$Var2 <- 'neutrophil.3'
pval.3$Var2 <- as.factor(pval.3$Var2)
pval.3$group1 <- sub(".*-", "", rownames(pval.3))
pval.3$group2 <- sub("-.*", "", rownames(pval.3))

pval.4 <- data.frame(dunnett_test_results$neutrophil.4$gTrac)
pval.4$signif <- NA
pval.4$signif[which(pval.4$pval <= 0.05)] <- '*'
pval.4$Var2 <- 'neutrophil.4'
pval.4$Var2 <- as.factor(pval.4$Var2)
pval.4$group1 <- sub(".*-", "", rownames(pval.4))
pval.4$group2 <- sub("-.*", "", rownames(pval.4))

pval.5 <- data.frame(dunnett_test_results$neutrophil.5$gTrac)
pval.5$signif <- NA
pval.5$signif[which(pval.5$pval <= 0.05)] <- '*'
pval.5$Var2 <- 'neutrophil.5'
pval.5$Var2 <- as.factor(pval.5$Var2)
pval.5$group1 <- sub(".*-", "", rownames(pval.5))
pval.5$group2 <- sub("-.*", "", rownames(pval.5))

pval.6 <- data.frame(dunnett_test_results$neutrophil.6$gTrac)
pval.6$signif <- NA
pval.6$signif[which(pval.6$pval <= 0.05)] <- '*'
pval.6$Var2 <- 'neutrophil.6'
pval.6$Var2 <- as.factor(pval.6$Var2)
pval.6$group1 <- sub(".*-", "", rownames(pval.6))
pval.6$group2 <- sub("-.*", "", rownames(pval.6))

pval.7 <- data.frame(dunnett_test_results$neutrophil.7$gTrac)
pval.7$signif <- NA
pval.7$signif[which(pval.7$pval <= 0.05)] <- '*'
pval.7$Var2 <- 'neutrophil.7'
pval.7$Var2 <- as.factor(pval.7$Var2)
pval.7$group1 <- sub(".*-", "", rownames(pval.7))
pval.7$group2 <- sub("-.*", "", rownames(pval.7))

pval.8 <- data.frame(dunnett_test_results$neutrophil.8$gTrac)
pval.8$signif <- NA
pval.8$signif[which(pval.8$pval <= 0.05)] <- '*'
pval.8$Var2 <- 'neutrophil.8'
pval.8$Var2 <- as.factor(pval.8$Var2)
pval.8$group1 <- sub(".*-", "", rownames(pval.8))
pval.8$group2 <- sub("-.*", "", rownames(pval.8))

pval.9 <- data.frame(dunnett_test_results$neutrophil.9$gTrac)
pval.9$signif <- NA
pval.9$signif[which(pval.9$pval <= 0.05)] <- '*'
pval.9$Var2 <- 'neutrophil.9'
pval.9$Var2 <- as.factor(pval.9$Var2)
pval.9$group1 <- sub(".*-", "", rownames(pval.9))
pval.9$group2 <- sub("-.*", "", rownames(pval.9))

pval.10 <- data.frame(dunnett_test_results$neutrophil.10$gTrac)
pval.10$signif <- NA
pval.10$signif[which(pval.10$pval <= 0.05)] <- '*'
pval.10$Var2 <- 'neutrophil.10'
pval.10$Var2 <- as.factor(pval.10$Var2)
pval.10$group1 <- sub(".*-", "", rownames(pval.10))
pval.10$group2 <- sub("-.*", "", rownames(pval.10))

pval.11 <- data.frame(dunnett_test_results$neutrophil.11$gTrac)
pval.11$signif <- NA
pval.11$signif[which(pval.11$pval <= 0.05)] <- '*'
pval.11$Var2 <- 'neutrophil.11'
pval.11$Var2 <- as.factor(pval.11$Var2)
pval.11$group1 <- sub(".*-", "", rownames(pval.11))
pval.11$group2 <- sub("-.*", "", rownames(pval.11))

dunnett.pvals <- rbind(pval.1, pval.2, pval.3, pval.4, pval.5, pval.6, pval.7, pval.8, pval.9, pval.10, pval.11)

dunnett.pvals$AAV <- dunnett.pvals$group2
dunnett.pvals$y.position <- 0.4

pdf('Figure_5B.pdf', width = 3, height = 2.7, bg = 'transparent')
ggplot(data = cell_freq, mapping = aes(x = AAV,
                                       y = value,
                                       group = AAV)) +
  stat_summary(aes(fill = AAV),
               fun.data = mean_se, geom = 'col', color = 'black', position = position_dodge(1), size = 0.1) +
  stat_summary(fun.data = mean_se, geom = 'errorbar', size = 0.1, width = 0.1, position = position_dodge(1)) +
  geom_jitter(aes(fill = AAV),
              position = position_jitterdodge(dodge.width = 1, jitter.width = 0.5),
              size = 0.1) +
  scale_fill_manual(values = color_group) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks = element_line(linewidth = 0.1), axis.line = element_line(linewidth = 0.1),
        axis.text.y = element_text(color = 'black', size = 6), axis.title.y = element_text(size = 6),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6),
        text = element_text(size = 6), strip.background = element_rect(linewidth = 0.1), strip.text = element_text(margin = margin(t=1,b=1)),
        legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.3, 'lines'), legend.margin = margin(0.1, 0.1, 0.1, 0.1, 'mm'), legend.spacing = unit(0, 'mm'))  +
  stat_pvalue_manual(dunnett.pvals, label = "signif", tip.length = 0, step.group.by = 'Var2', step.increase = 0.1, hide.ns = F, size = 2, bracket.size = 0.1) +
  facet_wrap(~Var2, ncol = 4)
dev.off()

# Figure 5C----
# Violin plot of Siglecf expression among all neutrophil subsets.
Idents(neuts_0.4) <- 'annotation_1'
pdf('Figure_5C.pdf', width = 1.6, height = 1.3, bg = 'transparent')
gg1 <- VlnPlot(neuts_0.4, 'Siglecf', pt.size = 0) + NoLegend() +
  labs(y = 'Expression level', title = 'Siglecf') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 6), axis.text = element_text(size = 6), axis.text.y = element_blank(),
        axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1),
        plot.title = element_text(size = 6, face = 'oblique', margin = margin(0,0,0,0)),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'), legend.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA))
gg1$layers[[1]]$aes_params$size <- 0.1 #adjust linewidth of violins
gg1
dev.off()

# Figure 5D----
# Bubble plot of scaled expression of all detected Ccl and Cxcl genes in all CAF subsets
ccl_genes <- grep('^Ccl', rownames(fbaav_0.4), value = T)
ccl_cxcl_genes <- grep('^(Ccl|Cxcl)', rownames(fbaav_0.4), value = TRUE)
ccl_cxcl_genes_ordered <- c('Ccl1', 'Ccl2', 'Ccl3', 'Ccl4', 'Ccl5', 'Ccl6', 'Ccl7', 'Ccl8', 'Ccl9', 'Ccl11', 'Ccl12', 'Ccl17', 'Ccl19', 'Ccl20', 'Ccl21a', 'Ccl21b', 'Ccl22', 'Ccl24', 'Ccl25', 'Ccl26', 'Ccl27b', 'Ccl28',
                            'Cxcl1', 'Cxcl2', 'Cxcl3', 'Cxcl5', 'Cxcl9', 'Cxcl10', 'Cxcl12', 'Cxcl13', 'Cxcl14', 'Cxcl15','Cxcl16', 'Cxcl17')

pdf('Figure_5D.pdf', width = 4, height = 1.1)
DotPlot(fbaav_0.4,
        features = ccl_cxcl_genes_ordered,
        cols = 'RdBu', assay = "RNA", dot.scale = 2) +
  theme(text = element_text(size = 6), axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5, size = 4),
        axis.text.y=element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1),
        legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.1, 'lines'))
dev.off()
