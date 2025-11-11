# MANUSCRIPT TITLE: Discovery of a Linked Constellation of Gene Expression Revealed by Local Editing of Fibroblasts in Tumors
# Companion code: FIGURE S5
# Date: 2025-11-06

library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(readxl)
library(data.table)
library(survival)
library(survminer)
library(plyr)
library(reshape2)
library(RcppML)
library(Matrix)
library(magrittr)
library(DescTools)

# Figure S5A----
# Gene expression of select genes (homeostatic Pi16; inflammatory Il6 and Lif; myofibroblast-associated Lrrc15, Acta2, and Ncam1; AAV-gTgfbr2-associated Col18a1) along pseudotime in CAFs
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

cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(fb[["RNA"]])

AFD_genes <- c('Pi16', 'Il6', 'Lif', 'Lrrc15', 'Acta2', 'Ncam1', 'Col18a1')
AFD_lineage_cds <- cds[rowData(cds)@rownames %in% AFD_genes,]

trajectories = plot_genes_in_pseudotime(AFD_lineage_cds,
                                        color_cells_by="clustering_1",
                                        trend_formula = "~ splines::ns(pseudotime, df=7)",
                                        min_expr=1,label_by_short_name = FALSE,cell_size = 0.1,vertical_jitter = 0.01,
                                        panel_order = AFD_genes)
pdf('Figure_S5A.pdf', width = 2, height = 3, bg = 'transparent')
trajectories +
  theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 6), axis.title.y = element_text(size = 6),
        axis.ticks = element_line(linewidth = 0.1), axis.line.y = element_line(linewidth = 0.1), element_line(linewidth = 0.1),
        legend.position = 'right', legend.title = element_blank(), legend.text = element_text(size = 6), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.3, 'lines'), legend.margin = margin(0, 0, 0, 0, 'mm'), legend.spacing = unit(0, 'mm'),
        strip.text = element_text(size = 6, face = 'oblique', margin = margin(t = 0, r = 0, b = 0, l = 0)),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'),
        panel.spacing = unit(0, "lines"),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA))
dev.off()

# Figure S5B----
# Pseudotime distribution along Monocle 2-generated trajectory of AAV-gTgfbr2 CAFs
library(monocle)
library(Seurat)
library(dplyr)
library(ggplot2)
library(grid)

seurat_obj <- readRDS("fbaav_0.4.rds")
create_monocle_cds <- function(seurat_subset, subset_name) {
  
  cat("Processing", subset_name, "subset...\n")
  expr_matrix <- as.matrix(GetAssayData(seurat_subset, slot = "counts"))
  
  cell_metadata <- seurat_subset@meta.data
  
  gene_metadata <- data.frame(
    gene_short_name = rownames(expr_matrix),
    row.names = rownames(expr_matrix)
  )
  
  cds <- newCellDataSet(expr_matrix,
                        phenoData = new("AnnotatedDataFrame", data = cell_metadata),
                        featureData = new("AnnotatedDataFrame", data = gene_metadata),
                        lowerDetectionLimit = 0.5,
                        expressionFamily = negbinomial.size())
  
  cds <- detectGenes(cds, min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
  cds <- cds[expressed_genes, ]
  
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  
  clustering_DEG_genes <- differentialGeneTest(cds,
                                               fullModelFormulaStr = '~ident',
                                               cores = 8)
  
  cat("Found", nrow(clustering_DEG_genes), "DEG genes for", subset_name, "\n")
  
  my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
  cds <- setOrderingFilter(cds, ordering_genes = my_ordering_genes)
  cds <- reduceDimension(cds, method = 'DDRTree')
  cds <- orderCells(cds)
  
  return(cds)
}

gtgfbr2_cells <- subset(seurat_obj, subset = group == "gTgfbr2")
gtgfbr2_cds <- create_monocle_cds(gtgfbr2_cells, "gTGFBR2")

plot_cell_trajectory(gtgfbr2_cds, color_by = "Pseudotime") + 
  ggtitle("gTGFBR2 Trajectory - by Pseudotime")

# Figure S5C----
# AAV-gTgfbr2 CAF subsets labelled on Monocle 2-generated pseudotime trajectory.
plot_cell_trajectory(gtgfbr2_cds, color_by = "ident") + 
  ggtitle("gTGFBR2 Trajectory - by Cell Type")

# Figure S5D----
# Cluster composition along pseudotime by branch
cell_data <- pData(gtgfbr2_cds)
analyze_pseudotime_bins <- function(cell_subset, branch_name, n_bins = 5) {
  
  cell_subset <- cell_subset %>%
    mutate(
      Pseudotime_bin = cut(Pseudotime, 
                           breaks = n_bins, 
                           labels = paste0("Bin", 1:n_bins),
                           include.lowest = TRUE)
    )
  
  bin_frequencies <- cell_subset %>%
    group_by(Pseudotime_bin, ident) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(Pseudotime_bin) %>%
    mutate(
      total = sum(count),
      frequency = count / total,
      branch = branch_name
    ) %>%
    ungroup()
  
  bin_ranges <- cell_subset %>%
    group_by(Pseudotime_bin) %>%
    summarise(
      min_pseudotime = min(Pseudotime),
      max_pseudotime = max(Pseudotime),
      mean_pseudotime = mean(Pseudotime),
      .groups = 'drop'
    )
  
  bin_frequencies <- bin_frequencies %>%
    left_join(bin_ranges, by = "Pseudotime_bin")
  
  return(bin_frequencies)
}

Branch1_freq <- analyze_pseudotime_bins(Branch1_cells, "Branch1", n_bins = 5)
Branch2_freq <- analyze_pseudotime_bins(Branch2_cells, "Branch2", n_bins = 5)
combined_freq <- bind_rows(Branch1_freq, Branch2_freq)
p1 <- ggplot(Branch1_freq, aes(x = Pseudotime_bin, y = frequency, fill = ident, group = ident)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  labs(title = "Branch1",
       x = "Pseudotime Bin",
       y = "Frequency",
       fill = "Cluster") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))
print(p1)
p2 <- ggplot(Branch2_freq, aes(x = Pseudotime_bin, y = frequency, fill = ident, group = ident)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  labs(title = "Branch2",
       x = "Pseudotime Bin",
       y = "Frequency",
       fill = "Cluster") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))
print(p2)

# Figure S5E----
# AAV-gTgfbr2 CAFs labelled by states that change as a function of pseudotime relative to the branchpoint.
plot_cell_trajectory(gtgfbr2_cds, color_by = "State") + 
  ggtitle("gTGFBR2 Trajectory - by State")

# Fgure S5F----
# Heatmap of differentially expressed genes between branch 1 (myCAF_1 trajectory) and branch 2 (Col18a1) trajectory
BEAM_res <- BEAM(
  gtgfbr2_cds,
  branch_point      = 1,
  progenitor_method = "duplicate",  
  cores             = 8
)

BEAM_res <- BEAM_res[order(BEAM_res$qval), ]
top_beam_genes_gtgfbr2 <- rownames(BEAM_res)[1:50]
gtgfbr2_branched_heatmap <- plot_genes_branched_heatmap(
  gtgfbr2_cds[top_beam_genes_gtgfbr2, ],
  branch_point        = 1,
  num_clusters        = 2,
  cores               = 8,
  use_gene_short_name = TRUE,
  show_rownames       = TRUE,
  return_heatmap      = TRUE
)

# Figure S5G----
# Gene expression of seven transcription factors differentially expressed between branch 1 (myCAF_1 trajectory) and branch 2 (Col18a1 trajectory) along pseudotime
tf_candidates <- c("Fos", "Fosb", "Jun", "Junb", "Egr1", "Nr4a1", "Zfp36")
# Generate the base plot
p <- plot_genes_branched_pseudotime(
  gtgfbr2_cds[tf_candidates, ],
  branch_point = 1,
  ncol = 4
)

# Figure S5H----
# Mean squared error (MSE) of an NMF model vs. model factorization rank (k = 1-50) in the fibroblast object
# Delta between MSE values from successive NMF model ranks to better visualize inflection point of top graph as defined by the point where the values plateau
fbaav_0.4<-FindVariableFeatures(fbaav_0.4,selection.method = 'vst',nfeatures = 3000,assay = 'RNA')
var.genes = VariableFeatures(fbaav_0.4,assay = 'RNA')
counts = GetAssayData(fbaav_0.4,assay = 'RNA',layer = "counts")
var.genes = row.names(counts)
counts = counts[match(var.genes,row.names(counts)),]
pass = Matrix::rowSums(counts>0)
pass = pass/ncol(fbaav_0.4)
pass = pass[pass>0.02]
length(pass)

bad.genes = c()
bad.genes = c(bad.genes,grep(pattern = "^Rps", x = names(pass)))
bad.genes = c(bad.genes,grep(pattern = "^Rpl", x = names(pass)))
bad.genes = c(bad.genes,grep(pattern = "^mt-", x = names(pass)))
pass.cleaned = pass[-bad.genes]
length(pass.cleaned)

rna = GetAssayData(fbaav_0.4,layer = 'data',assay = 'RNA')
var.genes.nmf = names(pass.cleaned)
var.norm.counts = rna[var.genes.nmf,] 
mse.array = c()
for (i in seq(1,50,1)){
  test = RcppML::nmf(var.norm.counts,k=i,diag = TRUE,L1 = c(0,0),verbose = FALSE,tol = 5e-3)
  mse.array = c(mse.array,RcppML::mse(var.norm.counts,w = test$w, d=test$d,h = test$h))
  print(paste('computing for k=',as.character(i)))
}
df.mse.array <- as.data.frame(mse.array)
mse.array2 = mse.array[2:length(mse.array)]-mse.array[1:length(mse.array)-1]
df.mse.array2 <- as.data.frame(mse.array2)

pdf('Figure_S5H.pdf', width = 1.75, height = 3.2, bg = 'transparent')
gg1 <- ggplot(df.mse.array, aes(x=c(1:50), y=mse.array)) +
  geom_point(size = 0.1) +
  labs(x = 'NMF factor rank', y = 'mean squared error (MSE)') +
  geom_vline(xintercept = 18, linetype = 'dashed', color = 'red', size = 0.3) +
  theme_classic() +
  theme(axis.text.x = element_text(color = 'black', size = 6), axis.title.x = element_text(color = 'black', size = 6), axis.ticks = element_line(linewidth = 0.1), axis.line = element_line(linewidth = 0.1),
        axis.text.y = element_text(color = 'black', size = 6), axis.title.y = element_text(color = 'black', size = 6),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6),
        text = element_text(size = 6), strip.background = element_rect(linewidth = 0.1), strip.text = element_text(margin = margin(t=1,b=1)),
        legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.3, 'lines'), legend.margin = margin(0, 0, 0, 0), legend.spacing = unit(0, 'mm'),
        plot.margin = margin(0,0,0,0))
gg2 <- ggplot(df.mse.array2, aes(x=c(1:49), y=mse.array2)) +
  geom_point(size = 0.1) +
  labs(x = 'NMF factor rank', y = '(MSE of rank N+1) - (MSE of rank N)') +
  geom_vline(xintercept = 18, linetype = 'dashed', color = 'red', size = 0.3) +
  geom_hline(yintercept = -4.338353e-04, linetype = 'dashed', color = 'grey', size = 0.3) +
  theme_classic() +
  theme(axis.text.x = element_text(color = 'black', size = 6), axis.title.x = element_text(color = 'black', size = 6), axis.ticks = element_line(linewidth = 0.1), axis.line = element_line(linewidth = 0.1),
        axis.text.y = element_text(color = 'black', size = 6), axis.title.y = element_text(color = 'black', size = 6),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6),
        text = element_text(size = 6), strip.background = element_rect(linewidth = 0.1), strip.text = element_text(margin = margin(t=1,b=1)),
        legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.3, 'lines'), legend.margin = margin(0, 0, 0, 0), legend.spacing = unit(0, 'mm'),
        plot.margin = margin(0,0,0,0))
gg1/gg2
dev.off()

# Figure S5I----
# Bar charts of average expression of NMF factors in CAFs grouped by AAV-gRNA condition 

# CALCULATE NMF OUTPUT
fb_18_programs_diag_L1_0_0 = RcppML::nmf(var.norm.counts,k=18,seed = 12345, diag = TRUE,L1 = c(0,0),verbose = T,maxit = 1000,tol = 1e-5)
h =fb_18_programs_diag_L1_0_0$h
colnames(h)<-colnames(fbaav_0.4)
row.names(h)<-paste('factor-',as.character(c(1:nrow(h))),sep ='')
fbaav_0.4[["NMF"]] <- CreateAssayObject(counts = h)

## PLOT NMF FACTOR AVERAGE LEVELS AS BAR PLOTS:
avg.nmf.fb = AverageExpression(fbaav_0.4, assays = 'NMF',slot = 'counts',features = rownames(fbaav_0.4@assays$NMF),group.by = 'orig.ident')
nmf_table <- as.data.frame(avg.nmf.fb)
nmf_table$NMF_factor <- row.names(nmf_table)
nmf_table <- melt(nmf_table)
nmf_table$variable <- as.character(nmf_table$variable)
nmf_table$variable <- substr(nmf_table$variable, 5, nchar(nmf_table$variable)) #remove first 4 entries
#add AAV annotation
nmf_table$AAV <- NA
nmf_table$AAV[which(nmf_table$variable %in% c('gTrac.A1', 'gTrac.A2', 'gTrac.B1', 'gTrac.B2'))] <- 'gTrac'
nmf_table$AAV[which(nmf_table$variable %in% c('gOsmr.A3', 'gOsmr.A4', 'gOsmr.B3', 'gOsmr.B4'))] <- 'gOsmr'
nmf_table$AAV[which(nmf_table$variable %in% c('gTgfbr2.A5', 'gTgfbr2.A6', 'gTgfbr2.B5', 'gTgfbr2.B6'))] <- 'gTgfbr2'
nmf_table$AAV[which(nmf_table$variable %in% c('gIl1r1.A7', 'gIl1r1.A8', 'gIl1r1.B7', 'gIl1r1.B8'))] <- 'gIl1r1'
#add replicate annotation
nmf_table$replicate <- NA
nmf_table$replicate[which(nmf_table$variable %in% c('gTrac.A1', 'gOsmr.A3', 'gTgfbr2.A5', 'gIl1r1.A7'))] <- '1'
nmf_table$replicate[which(nmf_table$variable %in% c('gTrac.A2', 'gOsmr.A4', 'gTgfbr2.A6', 'gIl1r1.A8'))] <- '2'
nmf_table$replicate[which(nmf_table$variable %in% c('gTrac.B1', 'gOsmr.B3', 'gTgfbr2.B5', 'gIl1r1.B7'))] <- '3'
nmf_table$replicate[which(nmf_table$variable %in% c('gTrac.B2', 'gOsmr.B4', 'gTgfbr2.B6', 'gIl1r1.B8'))] <- '4'
my_levels <- c('gTrac', 'gOsmr', 'gTgfbr2', 'gIl1r1')
nmf_table$AAV <- factor(nmf_table$AAV, levels = my_levels)
NMF_factor_levels <- c('factor-1','factor-2','factor-3','factor-4','factor-5','factor-6','factor-7',"factor-8","factor-9","factor-10",
                       "factor-11","factor-12","factor-13","factor-14","factor-15","factor-16", 'factor-17', 'factor-18')
nmf_table$NMF_factor <- factor(nmf_table$NMF_factor, levels = NMF_factor_levels)
nmf_table$replicate <- as.factor(nmf_table$replicate)

#cast data.frame for stat analysis
df <- reshape2::dcast(nmf_table, variable ~ NMF_factor, value = 'value')
#add AAV annotation
df$AAV <- NA
df$AAV[which(df$variable %in% c('gTrac.A1', 'gTrac.A2', 'gTrac.B1', 'gTrac.B2'))] <- 'gTrac'
df$AAV[which(df$variable %in% c('gOsmr.A3', 'gOsmr.A4', 'gOsmr.B3', 'gOsmr.B4'))] <- 'gOsmr'
df$AAV[which(df$variable %in% c('gTgfbr2.A5', 'gTgfbr2.A6', 'gTgfbr2.B5', 'gTgfbr2.B6'))] <- 'gTgfbr2'
df$AAV[which(df$variable %in% c('gIl1r1.A7', 'gIl1r1.A8', 'gIl1r1.B7', 'gIl1r1.B8'))] <- 'gIl1r1'
df$AAV <- factor(df$AAV, levels = c('gTrac', 'gOsmr', 'gTgfbr2', 'gIl1r1'))
#add replicate annotation
df$replicate <- NA
df$replicate[which(df$variable %in% c('gTrac.A1', 'gOsmr.A3', 'gTgfbr2.A5', 'gIl1r1.A7'))] <- '1'
df$replicate[which(df$variable %in% c('gTrac.A2', 'gOsmr.A4', 'gTgfbr2.A6', 'gIl1r1.A8'))] <- '2'
df$replicate[which(df$variable %in% c('gTrac.B1', 'gOsmr.B3', 'gTgfbr2.B5', 'gIl1r1.B7'))] <- '3'
df$replicate[which(df$variable %in% c('gTrac.B2', 'gOsmr.B4', 'gTgfbr2.B6', 'gIl1r1.B8'))] <- '4'
colnames(df) <- gsub("-", "_", colnames(df))

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

pval.1 <- data.frame(dunnett_test_results$factor_1$gTrac)
pval.1$signif <- NA
pval.1$signif[which(pval.1$pval <= 0.05)] <- '*'
pval.1$Var2 <- 'factor-1'
pval.1$Var2 <- as.factor(pval.1$Var2)
pval.1$group1 <- sub(".*-", "", rownames(pval.1))
pval.1$group2 <- sub("-.*", "", rownames(pval.1))

pval.2 <- data.frame(dunnett_test_results$factor_2$gTrac)
pval.2$signif <- NA
pval.2$signif[which(pval.2$pval <= 0.05)] <- '*'
pval.2$Var2 <- 'factor-2'
pval.2$Var2 <- as.factor(pval.2$Var2)
pval.2$group1 <- sub(".*-", "", rownames(pval.2))
pval.2$group2 <- sub("-.*", "", rownames(pval.2))

pval.3 <- data.frame(dunnett_test_results$factor_3$gTrac)
pval.3$signif <- NA
pval.3$signif[which(pval.3$pval <= 0.05)] <- '*'
pval.3$Var2 <- 'factor-3'
pval.3$Var2 <- as.factor(pval.3$Var2)
pval.3$group1 <- sub(".*-", "", rownames(pval.3))
pval.3$group2 <- sub("-.*", "", rownames(pval.3))

pval.4 <- data.frame(dunnett_test_results$factor_4$gTrac)
pval.4$signif <- NA
pval.4$signif[which(pval.4$pval <= 0.05)] <- '*'
pval.4$Var2 <- 'factor-4'
pval.4$Var2 <- as.factor(pval.4$Var2)
pval.4$group1 <- sub(".*-", "", rownames(pval.4))
pval.4$group2 <- sub("-.*", "", rownames(pval.4))

pval.5 <- data.frame(dunnett_test_results$factor_5$gTrac)
pval.5$signif <- NA
pval.5$signif[which(pval.5$pval <= 0.05)] <- '*'
pval.5$Var2 <- 'factor-5'
pval.5$Var2 <- as.factor(pval.5$Var2)
pval.5$group1 <- sub(".*-", "", rownames(pval.5))
pval.5$group2 <- sub("-.*", "", rownames(pval.5))

pval.6 <- data.frame(dunnett_test_results$factor_6$gTrac)
pval.6$signif <- NA
pval.6$signif[which(pval.6$pval <= 0.05)] <- '*'
pval.6$Var2 <- 'factor-6'
pval.6$Var2 <- as.factor(pval.6$Var2)
pval.6$group1 <- sub(".*-", "", rownames(pval.6))
pval.6$group2 <- sub("-.*", "", rownames(pval.6))

pval.7 <- data.frame(dunnett_test_results$factor_7$gTrac)
pval.7$signif <- NA
pval.7$signif[which(pval.7$pval <= 0.05)] <- '*'
pval.7$Var2 <- 'factor-7'
pval.7$Var2 <- as.factor(pval.7$Var2)
pval.7$group1 <- sub(".*-", "", rownames(pval.7))
pval.7$group2 <- sub("-.*", "", rownames(pval.7))

pval.8 <- data.frame(dunnett_test_results$factor_8$gTrac)
pval.8$signif <- NA
pval.8$signif[which(pval.8$pval <= 0.05)] <- '*'
pval.8$Var2 <- 'factor-8'
pval.8$Var2 <- as.factor(pval.8$Var2)
pval.8$group1 <- sub(".*-", "", rownames(pval.8))
pval.8$group2 <- sub("-.*", "", rownames(pval.8))

pval.9 <- data.frame(dunnett_test_results$factor_9$gTrac)
pval.9$signif <- NA
pval.9$signif[which(pval.9$pval <= 0.05)] <- '*'
pval.9$Var2 <- 'factor-9'
pval.9$Var2 <- as.factor(pval.9$Var2)
pval.9$group1 <- sub(".*-", "", rownames(pval.9))
pval.9$group2 <- sub("-.*", "", rownames(pval.9))

pval.10 <- data.frame(dunnett_test_results$factor_10$gTrac)
pval.10$signif <- NA
pval.10$signif[which(pval.10$pval <= 0.05)] <- '*'
pval.10$Var2 <- 'factor-10'
pval.10$Var2 <- as.factor(pval.10$Var2)
pval.10$group1 <- sub(".*-", "", rownames(pval.10))
pval.10$group2 <- sub("-.*", "", rownames(pval.10))

pval.11 <- data.frame(dunnett_test_results$factor_11$gTrac)
pval.11$signif <- NA
pval.11$signif[which(pval.11$pval <= 0.05)] <- '*'
pval.11$Var2 <- 'factor-11'
pval.11$Var2 <- as.factor(pval.11$Var2)
pval.11$group1 <- sub(".*-", "", rownames(pval.11))
pval.11$group2 <- sub("-.*", "", rownames(pval.11))

pval.12 <- data.frame(dunnett_test_results$factor_12$gTrac)
pval.12$signif <- NA
pval.12$signif[which(pval.12$pval <= 0.05)] <- '*'
pval.12$Var2 <- 'factor-12'
pval.12$Var2 <- as.factor(pval.12$Var2)
pval.12$group1 <- sub(".*-", "", rownames(pval.12))
pval.12$group2 <- sub("-.*", "", rownames(pval.12))

pval.13 <- data.frame(dunnett_test_results$factor_13$gTrac)
pval.13$signif <- NA
pval.13$signif[which(pval.13$pval <= 0.05)] <- '*'
pval.13$Var2 <- 'factor-13'
pval.13$Var2 <- as.factor(pval.13$Var2)
pval.13$group1 <- sub(".*-", "", rownames(pval.13))
pval.13$group2 <- sub("-.*", "", rownames(pval.13))

pval.14 <- data.frame(dunnett_test_results$factor_14$gTrac)
pval.14$signif <- NA
pval.14$signif[which(pval.14$pval <= 0.05)] <- '*'
pval.14$Var2 <- 'factor-14'
pval.14$Var2 <- as.factor(pval.14$Var2)
pval.14$group1 <- sub(".*-", "", rownames(pval.14))
pval.14$group2 <- sub("-.*", "", rownames(pval.14))

pval.15 <- data.frame(dunnett_test_results$factor_15$gTrac)
pval.15$signif <- NA
pval.15$signif[which(pval.15$pval <= 0.05)] <- '*'
pval.15$Var2 <- 'factor-15'
pval.15$Var2 <- as.factor(pval.15$Var2)
pval.15$group1 <- sub(".*-", "", rownames(pval.15))
pval.15$group2 <- sub("-.*", "", rownames(pval.15))

pval.16 <- data.frame(dunnett_test_results$factor_16$gTrac)
pval.16$signif <- NA
pval.16$signif[which(pval.16$pval <= 0.05)] <- '*'
pval.16$Var2 <- 'factor-16'
pval.16$Var2 <- as.factor(pval.16$Var2)
pval.16$group1 <- sub(".*-", "", rownames(pval.16))
pval.16$group2 <- sub("-.*", "", rownames(pval.16))

pval.17 <- data.frame(dunnett_test_results$factor_17$gTrac)
pval.17$signif <- NA
pval.17$signif[which(pval.17$pval <= 0.05)] <- '*'
pval.17$Var2 <- 'factor-17'
pval.17$Var2 <- as.factor(pval.17$Var2)
pval.17$group1 <- sub(".*-", "", rownames(pval.17))
pval.17$group2 <- sub("-.*", "", rownames(pval.17))

pval.18 <- data.frame(dunnett_test_results$factor_18$gTrac)
pval.18$signif <- NA
pval.18$signif[which(pval.18$pval <= 0.05)] <- '*'
pval.18$Var2 <- 'factor-18'
pval.18$Var2 <- as.factor(pval.18$Var2)
pval.18$group1 <- sub(".*-", "", rownames(pval.18))
pval.18$group2 <- sub("-.*", "", rownames(pval.18))

dunnett.pvals <- rbind(pval.1, pval.2, pval.3, pval.4, pval.5, pval.6, pval.7, pval.8, pval.9, pval.10, pval.11,pval.12, pval.13, pval.14, pval.15, pval.16, pval.17, pval.18)
dunnett.pvals$AAV <- dunnett.pvals$group2
dunnett.pvals$y.position <- 3.5e-04
colnames(dunnett.pvals)[colnames(dunnett.pvals) == 'Var2'] <- 'NMF_factor'
color_group <- c('lightgrey', '#1A85FF', '#D41159', '#40B0A6')

plot_factors <- nmf_table
plot_factors$NMF_factor <- gsub('factor-', '', plot_factors$NMF_factor)
plot_factors$NMF_factor <- as.factor(plot_factors$NMF_factor)
dunnett.pvals$NMF_factor <- gsub('factor-', '', dunnett.pvals$NMF_factor)
dunnett.pvals$NMF_factor <- as.factor(dunnett.pvals$NMF_factor)

NMF_factor_levels <- c('1','2','3','4','5','6','7',"8","9","10",
                       "11","12","13","14","15","16", '17', '18')
plot_factors$NMF_factor <- factor(plot_factors$NMF_factor, levels = NMF_factor_levels)
dunnett.pvals$NMF_factor <- factor(dunnett.pvals$NMF_factor, levels = NMF_factor_levels)

pdf('Figure_S5I.pdf', width = 2.3, height = 2.8, bg = 'transparent')
ggplot(data = plot_factors, mapping = aes(x = AAV,
                                          y = value,
                                          group = AAV)) +
  stat_summary(aes(fill = AAV),
               fun.data = mean_se, geom = 'col', color = 'black', position = position_dodge(1), size = 0.1, alpha = 0.75) +
  stat_summary(fun.data = mean_se, geom = 'errorbar', size = 0.1, width = 0.1, position = position_dodge(1)) +
  geom_jitter(aes(fill = AAV),
              position = position_jitterdodge(dodge.width = 1, jitter.width = 0.5),
              size = 0.1) +
  scale_fill_manual(values = color_group) +
  labs(y = 'average NMF factor weight') +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks = element_line(linewidth = 0.1), axis.line = element_line(linewidth = 0.1),
        axis.text.y = element_text(color = 'black', size = 6), axis.title.y = element_text(size = 6),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6),
        text = element_text(size = 6), strip.background = element_rect(linewidth = 0.1), strip.text = element_text(margin = margin(t=1,b=1)),
        legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.3, 'lines'), legend.margin = margin(0, 0, 0, 0), legend.spacing = unit(0, 'mm'),
        plot.margin = margin(0,0,0,0))  +
  stat_pvalue_manual(dunnett.pvals, label = "signif", tip.length = 0, step.group.by = 'NMF_factor', step.increase = 0.1, hide.ns = F, size = 2, bracket.size = 0.1) +
  facet_wrap(~NMF_factor, ncol = 6)
dev.off()

# Figure S5J----
# Similarity of DEGs between mouse CAF subsets (columns) and mouse bone marrow stroma cells from Baryawno et al., 2019
fb_deg_markers <- read_excel("Supplementary_Tables/Supplementary Table 1 DEG scRNAseq.xlsx", sheet = 'fibroblasts')
aav.fb.deg <- split(fb_deg_markers, fb_deg_markers$clustering_1)
# Bone Marrow, mouse, Baryawno 2019
bm.stroma.deg <- readRDS("R_objects/public_datasets/bm.stroma.deg.rds")

### HYPERGEOMETRIC TEST----
# Define background:
N <- 19405 # number of genes captured by FLEX kit
# read in lists
listA <- aav.fb.deg
listB <- bm.stroma.deg

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
    setB <- listB[[b]] %>% filter(p_val_adj < 0.05 & avg_logFC > 0.25) %>% pull(gene)
    
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
y_levels <- c('Fibro-1', 'Fibro-2', 'Fibro-3', 'Fibro-4', 'Fibro-5', 'Lepr-MSC', 'OLC-1', 'OLC-2', 'Pericytes', 'EC-arterial', 'EC-arteriolar', 'EC-sinusoidal', 'Chondro', 'Chondro-progen', 'Chondro-hyper', 'Chondro-prehyper-2', 'Chondro-prol-rest')

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

pdf('Figure_S5J.pdf', width=2, height=1.5)
ggplot(data = melt_hgt, aes(x=Var1, y=Var2, fill=(-log10(P)))) + 
  geom_tile() +
  scale_fill_viridis_c(limits=c(0, 100), oob = scales::squish) +
  #scale_fill_viridis_c() +
  ggtitle('Baryanow 2019, mouse bone marrow stroma') +
  #scale_fill_gradient2(low = 'blue', high='red', mid = 'white', midpoint = 0) +
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

# Figure S5K----
# Kaplan-Meier survival plots of TCGA patients categorized into upper and lower quartiles by the stroma score identified in Combes et al., 2022
# Forest plot of hazard ratios ± 95% confidence interval and p-values in Cox regression model for overall survival of TCGA patients separated by indication and categorized into upper and lower quartiles by stroma score

tcga <- fread('R_objects/public_datasets/Select_inds_tcga_TPM_primary_annotated.tsv',  sep='\t', header=T)
signatures <- read_excel("Supplementary_Tables/Supplementary Table 4 signature genes.xlsx", sheet = 'signatures')
stromal.sig <- signatures$`stromal signature (Combes et al., 2022)`
stromal.sig <- stromal.sig[!is.na(stromal.sig)]

# select which signature to plot
sig_oi <- stromal.sig #select signature of interest
length(sig_oi)
# Subset data based on genes from signature
stromal.sig_exp <- subset(tcga, select = sig_oi)
stromal.sig_exp_log <- as.data.frame(lapply(stromal.sig_exp, function(x) log10(x + 1)))
# Calculate percentiles for each gene across samples
stromal.sig_percentiles <- apply(stromal.sig_exp_log, 2, function(x) ecdf(x)(x) * 100)
# Calculate the gene signature score for each sample (average of percentiles)
stromal.sig_score_percentiles <- rowMeans(stromal.sig_percentiles)
# Apply to original data table
tcga$stromal.sig_score_percentiles <- stromal.sig_score_percentiles
tcga$stromal.sig_score_cat <- ifelse(tcga$stromal.sig_score_percentiles >= quantile(tcga$stromal.sig_score_percentiles, 0.75), "stromal.sig_High",
                                            ifelse(tcga$stromal.sig_score_percentiles <= quantile(tcga$stromal.sig_score_percentiles, 0.25), "stromal.sig_Low", "Intermediate")) #separate by quartiles
table(tcga$stromal.sig_score_cat)
Survival_param_stromal.sig_score_all <- c('stromal.sig_score_cat', 'OS.time', 'OS.status', 'indication', 'sex', 'age')
TCGA_stromal.sig_score_surv_all <- subset(tcga, select = Survival_param_stromal.sig_score_all)
TCGA_stromal.sig_score_surv_all_sub_stromal.sig <- subset(TCGA_stromal.sig_score_surv_all, stromal.sig_score_cat %in% c('stromal.sig_High', 'stromal.sig_Low'))
Survival_table <- survfit(Surv(OS.time, OS.status) ~ stromal.sig_score_cat, data = TCGA_stromal.sig_score_surv_all_sub_stromal.sig)

do.cox.group <- function(df, group.order) {
  cox <- c()
  df$stromal.sig_score_cat <- factor(df$stromal.sig_score_cat, levels = group.order)
  if (length(unique(df$sex)) > 1) {
    cox <- coxph(Surv(OS.time, OS.status) ~ stromal.sig_score_cat + sex + age + indication, data = df) # include indication
  } else {
    cox <- coxph(Surv(OS.time, OS.status) ~ stromal.sig_score_cat + age + indication, data = df) # include indication
  }
  cox.summary <- summary(cox)
  nvar <- length(unique(df$stromal.sig_score_cat)) - 1
  HR <- round(cox.summary$conf.int[1:nvar, 1], 2)
  HR.lower <- round(cox.summary$conf.int[1:nvar, 3], 2)
  HR.upper <- round(cox.summary$conf.int[1:nvar, 4], 2)
  HR.range <- sprintf("(%.1f-%.1f)", HR.lower, HR.upper)
  coef <- cox.summary$coefficients[1:nvar, 1]
  coef.se <- cox.summary$coefficients[1:nvar, 3]
  Pval <- round(cox.summary$coefficients[1:nvar, 5], 4)
  stromal.sig_score_cat <- gsub("stromal.sig_score_cat", "", rownames(cox.summary$conf.int)[1:nvar])
  return(data.frame(stromal.sig_score_cat = stromal.sig_score_cat, HR = HR, HR.range = HR.range, coef = coef, coef.se = coef.se, Pval = Pval))
}

HR.pval.df <- do.cox.group(TCGA_stromal.sig_score_surv_all_sub_stromal.sig, group.order = c("stromal.sig_Low", "stromal.sig_High"))
annotation_text <- paste("HR:", HR.pval.df$HR, '\n', "p =", HR.pval.df$Pval)

sorted_category_names <- sort(TCGA_stromal.sig_score_surv_all_sub_stromal.sig$stromal.sig_score_cat) #order the names of the category ('sig_High' & 'sig_Low') in alphabetical order, so the legend is correctly annotated!
legend_labels <- paste0(unique(sorted_category_names), " (n=", Survival_table$n, ")")
plot <- ggsurvplot(Survival_table, conf.int = FALSE, surv.size = 1.2, xlab = 'Days', ylab = 'Survival', legend.title = "", legend.labs = legend_labels,
                   ylim = c(0, NA), xlim = c(0, 3000), break.x.by = 1000, palette = c('red', 'black'), pval = F, size = 0.25,
                   censor.shape = 124, censor.size = 1) +
  ggtitle(paste("stromal signature (Combes et al., 2022)"))
plot$plot <- plot$plot + guides(color = guide_legend(ncol = 1)) + # make legend in vertical alignment
  annotate("text", x = Inf, y = Inf, label = annotation_text, hjust = 1.1, vjust = 1.1, size = 6*0.3527) # add HR and pval annotation
plot

pdf('Figure_S5K_KaplanMeier.pdf', width = 1.6, height = 2.2, bg = 'transparent')
plot$plot + theme(axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1),
                  axis.text.x = element_text(size = 6, color = 'black'), axis.text.y = element_text(size = 6, color = 'black'),
                  axis.title.x = element_text(size = 6), axis.title.y = element_text(size = 6), title = element_text(size = 6),
                  legend.text = element_text(size = 6), 
                  plot.title = element_text(size = 6, face = 'plain', margin = margin(0,0,0,0)),
                  legend.background = element_rect(fill = 'transparent', color = NA),
                  plot.background = element_rect(fill = "transparent", color = NA),
                  panel.background = element_rect(fill = "transparent", color = NA))
dev.off()

# FOREST PLOT:
indications <- sort(unique(tcga$indication))
plots_list <- list()

for (indic in indications) {
  indication_data <- subset(tcga, indication == indic)
  stromal.sig_exp <- subset(indication_data, select = stromal.sig)
  stromal.sig_exp_log <- as.data.frame(lapply(stromal.sig_exp, function(x) log10(x + 1)))
  stromal.sig_percentiles <- apply(stromal.sig_exp_log, 2, function(x) ecdf(x)(x) * 100)
  stromal.sig_score_percentiles <- rowMeans(stromal.sig_percentiles)
  indication_data$stromal.sig_score_perc <- stromal.sig_score_percentiles
  indication_data$stromal.sig_score_cat <- ifelse(indication_data$stromal.sig_score_perc >= quantile(indication_data$stromal.sig_score_perc, 0.75), "stromal.sig_High",
                                                         ifelse(indication_data$stromal.sig_score_perc <= quantile(indication_data$stromal.sig_score_perc, 0.25), "stromal.sig_Low", "stromal.sig_Intermediate"))
  Survival_param_stromal.sig_score_all <- c('sample_name', 'stromal.sig_score_cat', 'OS.time', 'OS.status', 'indication', 'sex', 'age')
  TCGA_stromal.sig_score_surv_all <- subset(indication_data, select = Survival_param_stromal.sig_score_all)
  TCGA_stromal.sig_score_surv_all_sub_stromal.sig <- subset(TCGA_stromal.sig_score_surv_all, stromal.sig_score_cat %in% c('stromal.sig_High', 'stromal.sig_Low'))
  plots_list[[indic]] <- TCGA_stromal.sig_score_surv_all_sub_stromal.sig
}

df <- do.call(rbind, lapply(plots_list, function(x) as.data.frame(x, stringsAsFactors = FALSE)))

do.cox.group <- function(df, group.order) {
  cox <- c()
  df$stromal.sig_score_cat <- factor(df$stromal.sig_score_cat, levels = group.order)
  if (length(unique(df$sex)) > 1) {
    cox <- coxph(Surv(OS.time, OS.status) ~ stromal.sig_score_cat + sex + age, data = df)
  } else {
    cox <- coxph(Surv(OS.time, OS.status) ~ stromal.sig_score_cat + age, data = df)
  }
  cox.summary <- summary(cox)
  nvar <- length(unique(df$stromal.sig_score_cat)) - 1
  HR <- round(cox.summary$conf.int[1:nvar, 1], 2)
  HR.lower <- round(cox.summary$conf.int[1:nvar, 3], 2)
  HR.upper <- round(cox.summary$conf.int[1:nvar, 4], 2)
  HR.range <- sprintf("(%.1f-%.1f)", HR.lower, HR.upper)
  coef <- cox.summary$coefficients[1:nvar, 1]
  coef.se <- cox.summary$coefficients[1:nvar, 3]
  Pval <- round(cox.summary$coefficients[1:nvar, 5], 4)
  stromal.sig_score_cat <- gsub("stromal.sig_score_cat", "", rownames(cox.summary$conf.int)[1:nvar])
  return(data.frame(stromal.sig_score_cat = stromal.sig_score_cat, HR = HR, HR.range = HR.range, coef = coef, coef.se = coef.se, Pval = Pval))
}

plot.data <- df[, c("sex", "age", "OS.time", "OS.status", "indication", "stromal.sig_score_cat")]
HR.df <- ddply(plot.data, c("indication"), do.cox.group, group.order = c("stromal.sig_Low", "stromal.sig_High"))
HR.df <- HR.df %>%
  group_by(stromal.sig_score_cat) %>%
  mutate(adj.Pval = p.adjust(Pval, method = "BH"))# %>%

meta <- meta::metagen(
  TE = HR.df$coef, seTE = HR.df$coef.se, studlab = HR.df$indication,
  common = F, random = T, prediction = F, sm = "HR"
)

pdf('Figure_S5K_ForestPlot.pdf', width = 5, height = 5, bg = 'transparent')
meta::forest(meta,
             test.overall.random = F, zero.pval = T, hetstat = F,
             leftcols = c('studlab'),
             leftlabs = c('indication'),
             rightcols = c('pval'),
             rightlabs = c('pval'),
             digits.pval = 3,
             digits = c(2),
             fontsize = 12,
             spacing = 0.7,
             lwd = 0.5)
dev.off()

# Figure S5L----
# Kaplan-Meier survival plots of TCGA patients categorized into upper and lower quartiles by the fibroblast TGFβ response signature from Mariathasan et al., 2018. Hazard ratio (HR) and p value of Cox regression fit are shown.
# Forest plot of hazard ratios ± 95% confidence interval and p-values in Cox regression model for overall survival of TCGA patients separated by indication and categorized into upper and lower quartiles by fibroblast TGFβ response signature.
tgfb.mariathasan.sig <- signatures$`fibroblast TGFbeta response signature (Mariathasan et al., 2018)`
tgfb.mariathasan.sig <- tgfb.mariathasan.sig[!is.na(tgfb.mariathasan.sig)]

# select which signature to plot
sig_oi <- tgfb.mariathasan.sig #select signature of interest
length(sig_oi)
# Subset data based on genes from signature
tgfb.mariathasan.sig_exp <- subset(tcga, select = sig_oi)
tgfb.mariathasan.sig_exp_log <- as.data.frame(lapply(tgfb.mariathasan.sig_exp, function(x) log10(x + 1)))
# Calculate percentiles for each gene across samples
tgfb.mariathasan.sig_percentiles <- apply(tgfb.mariathasan.sig_exp_log, 2, function(x) ecdf(x)(x) * 100)
# Calculate the gene signature score for each sample (average of percentiles)
tgfb.mariathasan.sig_score_percentiles <- rowMeans(tgfb.mariathasan.sig_percentiles)
# Apply to original data table
tcga$tgfb.mariathasan.sig_score_percentiles <- tgfb.mariathasan.sig_score_percentiles
tcga$tgfb.mariathasan.sig_score_cat <- ifelse(tcga$tgfb.mariathasan.sig_score_percentiles >= quantile(tcga$tgfb.mariathasan.sig_score_percentiles, 0.75), "tgfb.mariathasan.sig_High",
                                     ifelse(tcga$tgfb.mariathasan.sig_score_percentiles <= quantile(tcga$tgfb.mariathasan.sig_score_percentiles, 0.25), "tgfb.mariathasan.sig_Low", "Intermediate")) #separate by quartiles
table(tcga$tgfb.mariathasan.sig_score_cat)
Survival_param_tgfb.mariathasan.sig_score_all <- c('tgfb.mariathasan.sig_score_cat', 'OS.time', 'OS.status', 'indication', 'sex', 'age')
TCGA_tgfb.mariathasan.sig_score_surv_all <- subset(tcga, select = Survival_param_tgfb.mariathasan.sig_score_all)
TCGA_tgfb.mariathasan.sig_score_surv_all_sub_tgfb.mariathasan.sig <- subset(TCGA_tgfb.mariathasan.sig_score_surv_all, tgfb.mariathasan.sig_score_cat %in% c('tgfb.mariathasan.sig_High', 'tgfb.mariathasan.sig_Low'))
Survival_table <- survfit(Surv(OS.time, OS.status) ~ tgfb.mariathasan.sig_score_cat, data = TCGA_tgfb.mariathasan.sig_score_surv_all_sub_tgfb.mariathasan.sig)

do.cox.group <- function(df, group.order) {
  cox <- c()
  df$tgfb.mariathasan.sig_score_cat <- factor(df$tgfb.mariathasan.sig_score_cat, levels = group.order)
  if (length(unique(df$sex)) > 1) {
    cox <- coxph(Surv(OS.time, OS.status) ~ tgfb.mariathasan.sig_score_cat + sex + age + indication, data = df) # include indication
  } else {
    cox <- coxph(Surv(OS.time, OS.status) ~ tgfb.mariathasan.sig_score_cat + age + indication, data = df) # include indication
  }
  cox.summary <- summary(cox)
  nvar <- length(unique(df$tgfb.mariathasan.sig_score_cat)) - 1
  HR <- round(cox.summary$conf.int[1:nvar, 1], 2)
  HR.lower <- round(cox.summary$conf.int[1:nvar, 3], 2)
  HR.upper <- round(cox.summary$conf.int[1:nvar, 4], 2)
  HR.range <- sprintf("(%.1f-%.1f)", HR.lower, HR.upper)
  coef <- cox.summary$coefficients[1:nvar, 1]
  coef.se <- cox.summary$coefficients[1:nvar, 3]
  Pval <- round(cox.summary$coefficients[1:nvar, 5], 4)
  tgfb.mariathasan.sig_score_cat <- gsub("tgfb.mariathasan.sig_score_cat", "", rownames(cox.summary$conf.int)[1:nvar])
  return(data.frame(tgfb.mariathasan.sig_score_cat = tgfb.mariathasan.sig_score_cat, HR = HR, HR.range = HR.range, coef = coef, coef.se = coef.se, Pval = Pval))
}

HR.pval.df <- do.cox.group(TCGA_tgfb.mariathasan.sig_score_surv_all_sub_tgfb.mariathasan.sig, group.order = c("tgfb.mariathasan.sig_Low", "tgfb.mariathasan.sig_High"))
annotation_text <- paste("HR:", HR.pval.df$HR, '\n', "p =", HR.pval.df$Pval)

sorted_category_names <- sort(TCGA_tgfb.mariathasan.sig_score_surv_all_sub_tgfb.mariathasan.sig$tgfb.mariathasan.sig_score_cat) #order the names of the category ('sig_High' & 'sig_Low') in alphabetical order, so the legend is correctly annotated!
legend_labels <- paste0(unique(sorted_category_names), " (n=", Survival_table$n, ")")
plot <- ggsurvplot(Survival_table, conf.int = FALSE, surv.size = 1.2, xlab = 'Days', ylab = 'Survival', legend.title = "", legend.labs = legend_labels,
                   ylim = c(0, NA), xlim = c(0, 3000), break.x.by = 1000, palette = c('red', 'black'), pval = F, size = 0.25,
                   censor.shape = 124, censor.size = 1) +
  ggtitle(paste("fibroblast TGFb response sig. (Mariathasan et al., 2018)"))
plot$plot <- plot$plot + guides(color = guide_legend(ncol = 1)) + # make legend in vertical alignment
  annotate("text", x = Inf, y = Inf, label = annotation_text, hjust = 1.1, vjust = 1.1, size = 6*0.3527) # add HR and pval annotation
plot

pdf('Figure_S5L_KaplanMeier.pdf', width = 1.6, height = 2.2, bg = 'transparent')
plot$plot + theme(axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1),
                  axis.text.x = element_text(size = 6, color = 'black'), axis.text.y = element_text(size = 6, color = 'black'),
                  axis.title.x = element_text(size = 6), axis.title.y = element_text(size = 6), title = element_text(size = 6),
                  legend.text = element_text(size = 6), 
                  plot.title = element_text(size = 6, face = 'plain', margin = margin(0,0,0,0)),
                  legend.background = element_rect(fill = 'transparent', color = NA),
                  plot.background = element_rect(fill = "transparent", color = NA),
                  panel.background = element_rect(fill = "transparent", color = NA))
dev.off()

# FOREST PLOT:
indications <- sort(unique(tcga$indication))
plots_list <- list()

for (indic in indications) {
  indication_data <- subset(tcga, indication == indic)
  tgfb.mariathasan.sig_exp <- subset(indication_data, select = tgfb.mariathasan.sig)
  tgfb.mariathasan.sig_exp_log <- as.data.frame(lapply(tgfb.mariathasan.sig_exp, function(x) log10(x + 1)))
  tgfb.mariathasan.sig_percentiles <- apply(tgfb.mariathasan.sig_exp_log, 2, function(x) ecdf(x)(x) * 100)
  tgfb.mariathasan.sig_score_percentiles <- rowMeans(tgfb.mariathasan.sig_percentiles)
  indication_data$tgfb.mariathasan.sig_score_perc <- tgfb.mariathasan.sig_score_percentiles
  indication_data$tgfb.mariathasan.sig_score_cat <- ifelse(indication_data$tgfb.mariathasan.sig_score_perc >= quantile(indication_data$tgfb.mariathasan.sig_score_perc, 0.75), "tgfb.mariathasan.sig_High",
                                                  ifelse(indication_data$tgfb.mariathasan.sig_score_perc <= quantile(indication_data$tgfb.mariathasan.sig_score_perc, 0.25), "tgfb.mariathasan.sig_Low", "tgfb.mariathasan.sig_Intermediate"))
  Survival_param_tgfb.mariathasan.sig_score_all <- c('sample_name', 'tgfb.mariathasan.sig_score_cat', 'OS.time', 'OS.status', 'indication', 'sex', 'age')
  TCGA_tgfb.mariathasan.sig_score_surv_all <- subset(indication_data, select = Survival_param_tgfb.mariathasan.sig_score_all)
  TCGA_tgfb.mariathasan.sig_score_surv_all_sub_tgfb.mariathasan.sig <- subset(TCGA_tgfb.mariathasan.sig_score_surv_all, tgfb.mariathasan.sig_score_cat %in% c('tgfb.mariathasan.sig_High', 'tgfb.mariathasan.sig_Low'))
  plots_list[[indic]] <- TCGA_tgfb.mariathasan.sig_score_surv_all_sub_tgfb.mariathasan.sig
}

df <- do.call(rbind, lapply(plots_list, function(x) as.data.frame(x, stringsAsFactors = FALSE)))

do.cox.group <- function(df, group.order) {
  cox <- c()
  df$tgfb.mariathasan.sig_score_cat <- factor(df$tgfb.mariathasan.sig_score_cat, levels = group.order)
  if (length(unique(df$sex)) > 1) {
    cox <- coxph(Surv(OS.time, OS.status) ~ tgfb.mariathasan.sig_score_cat + sex + age, data = df)
  } else {
    cox <- coxph(Surv(OS.time, OS.status) ~ tgfb.mariathasan.sig_score_cat + age, data = df)
  }
  cox.summary <- summary(cox)
  nvar <- length(unique(df$tgfb.mariathasan.sig_score_cat)) - 1
  HR <- round(cox.summary$conf.int[1:nvar, 1], 2)
  HR.lower <- round(cox.summary$conf.int[1:nvar, 3], 2)
  HR.upper <- round(cox.summary$conf.int[1:nvar, 4], 2)
  HR.range <- sprintf("(%.1f-%.1f)", HR.lower, HR.upper)
  coef <- cox.summary$coefficients[1:nvar, 1]
  coef.se <- cox.summary$coefficients[1:nvar, 3]
  Pval <- round(cox.summary$coefficients[1:nvar, 5], 4)
  tgfb.mariathasan.sig_score_cat <- gsub("tgfb.mariathasan.sig_score_cat", "", rownames(cox.summary$conf.int)[1:nvar])
  return(data.frame(tgfb.mariathasan.sig_score_cat = tgfb.mariathasan.sig_score_cat, HR = HR, HR.range = HR.range, coef = coef, coef.se = coef.se, Pval = Pval))
}

plot.data <- df[, c("sex", "age", "OS.time", "OS.status", "indication", "tgfb.mariathasan.sig_score_cat")]
HR.df <- ddply(plot.data, c("indication"), do.cox.group, group.order = c("tgfb.mariathasan.sig_Low", "tgfb.mariathasan.sig_High"))
HR.df <- HR.df %>%
  group_by(tgfb.mariathasan.sig_score_cat) %>%
  mutate(adj.Pval = p.adjust(Pval, method = "BH"))# %>%

meta <- meta::metagen(
  TE = HR.df$coef, seTE = HR.df$coef.se, studlab = HR.df$indication,
  common = F, random = T, prediction = F, sm = "HR"
)

pdf('Figure_S5L_ForestPlot.pdf', width = 5, height = 5, bg = 'transparent')
meta::forest(meta,
             test.overall.random = F, zero.pval = T, hetstat = F,
             leftcols = c('studlab'),
             leftlabs = c('indication'),
             rightcols = c('pval'),
             rightlabs = c('pval'),
             digits.pval = 3,
             digits = c(2),
             fontsize = 12,
             spacing = 0.7,
             lwd = 0.5)
dev.off()

# Figure S5M----
# Kaplan-Meier survival plots of PAAD patients from TCGA excluding PNET (pancreatic neuroendocrine tumor) diagnosed patients categorized into upper and lower quartiles by the ratio of AAV-gOsmr DEG signature from this study, AAV-gIl1r1 DEG signature from this study, or TGFβ response signature from Mariathasan et al. to the fibroblast progenitor signature
tcga_paad_no_pnet <- fread('R_objects/public_datasets/PAADnoPNET_tcga_TPM_primary_annotated.tsv', sep='\t', header=T)
rownames(tcga_paad_no_pnet) <- tcga_paad_no_pnet$sample_name

signatures <- read_excel("Supplementary_Tables/Supplementary Table 4 signature genes.xlsx", sheet = 'signatures')
aav.osmr.up.deg.sig <- signatures$`AAV-gOsmr DEG signature`
aav.osmr.up.deg.sig <- aav.osmr.up.deg.sig[!is.na(aav.osmr.up.deg.sig)]
aav.il1r1.up.deg.sig <- signatures$`AAV-gIl1r1 DEG signature`
aav.il1r1.up.deg.sig <- aav.il1r1.up.deg.sig[!is.na(aav.il1r1.up.deg.sig)]
tgfb.mariathasan.sig <- signatures$`fibroblast TGFbeta response signature (Mariathasan et al., 2018)`
tgfb.mariathasan.sig <- tgfb.mariathasan.sig[!is.na(tgfb.mariathasan.sig)]
progenitor.sig <- signatures$`progenitor signature (derived from Gao et al., 2024)`
progenitor.sig <- progenitor.sig[!is.na(progenitor.sig)]

# select signature of interest
signatureA <- aav.osmr.up.deg.sig
#signatureA <- aav.il1r1.up.deg.sig
#signatureA <- tgfb.mariathasan.sig

# Subset data based on genes from signature
signatureA_exp <- subset(tcga_paad_no_pnet, select = signatureA)
signatureA_exp_log <- as.data.frame(lapply(signatureA_exp, function(x) log10(x + 1)))
# Calculate percentiles for each gene across samples
signatureA_percentiles <- apply(signatureA_exp_log, 2, function(x) ecdf(x)(x) * 100)
# Calculate the gene signature score for each sample (average of percentiles)
signatureA_score_percentiles <- rowMeans(signatureA_percentiles)
# Apply to original data table
tcga_paad_no_pnet$signatureA_score_percentiles <- signatureA_score_percentiles
# Survival code
tcga_paad_no_pnet$signatureA_score_cat <- ifelse(tcga_paad_no_pnet$signatureA_score_percentiles >= quantile(tcga_paad_no_pnet$signatureA_score_percentiles, 0.75), "signatureA_High",
                                                            ifelse(tcga_paad_no_pnet$signatureA_score_percentiles <= quantile(tcga_paad_no_pnet$signatureA_score_percentiles, 0.25), "signatureA_Low", "signatureA_Intermediate"))
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
tcga_paad_no_pnet$signatureA_to_progenitor.sig_ratio_score <- tcga_paad_no_pnet$signatureA_score_percentiles / tcga_paad_no_pnet$progenitor.sig_score_percentiles

tcga_paad_no_pnet$signatureA_to_progenitor.sig_ratio_score_cat <- ifelse(tcga_paad_no_pnet$signatureA_to_progenitor.sig_ratio_score >= quantile(tcga_paad_no_pnet$signatureA_to_progenitor.sig_ratio_score, 0.75), "signatureA_to_progenitor.sig_ratio_High",
                                                                                    ifelse(tcga_paad_no_pnet$signatureA_to_progenitor.sig_ratio_score <= quantile(tcga_paad_no_pnet$signatureA_to_progenitor.sig_ratio_score, 0.25), "signatureA_to_progenitor.sig_ratio_Low", "signatureA_to_progenitor.sig_ratio_Intermediate"))
table(tcga_paad_no_pnet$signatureA_to_progenitor.sig_ratio_score_cat)

ggplot(tcga_paad_no_pnet, aes(x=signatureA_score_percentiles, y=progenitor.sig_score_percentiles, color = signatureA_to_progenitor.sig_ratio_score_cat)) + geom_point() + theme_bw()


Survival_param_signatureA_to_progenitor.sig_ratio_score_all <- c('signatureA_to_progenitor.sig_ratio_score_cat', 'OS.time', 'OS.status', 'indication', 'sex', 'age')
TCGA_signatureA_to_progenitor.sig_ratio_score_surv_all <- subset(tcga_paad_no_pnet, select = Survival_param_signatureA_to_progenitor.sig_ratio_score_all)
TCGA_signatureA_to_progenitor.sig_ratio_score_surv_all_sub <- subset(TCGA_signatureA_to_progenitor.sig_ratio_score_surv_all, signatureA_to_progenitor.sig_ratio_score_cat %in% c('signatureA_to_progenitor.sig_ratio_High', 'signatureA_to_progenitor.sig_ratio_Low'))

Survival_table <- survfit(Surv(OS.time, OS.status) ~ signatureA_to_progenitor.sig_ratio_score_cat, data = TCGA_signatureA_to_progenitor.sig_ratio_score_surv_all_sub)

# multivariate regression (i.e. sex & age)
do.cox.group <- function(df, group.order) {
  cox <- c()
  df$signatureA_to_progenitor.sig_ratio_score_cat <- factor(df$signatureA_to_progenitor.sig_ratio_score_cat, levels = group.order)
  if (length(unique(df$sex)) > 1) {
    cox <- coxph(Surv(OS.time, OS.status) ~ signatureA_to_progenitor.sig_ratio_score_cat + sex + age, data = df)
  } else {
    cox <- coxph(Surv(OS.time, OS.status) ~ signatureA_to_progenitor.sig_ratio_score_cat + age, data = df)
  }
  cox.summary <- summary(cox)
  nvar <- length(unique(df$signatureA_to_progenitor.sig_ratio_score_cat)) - 1
  HR <- round(cox.summary$conf.int[1:nvar, 1], 2)
  HR.lower <- round(cox.summary$conf.int[1:nvar, 3], 2)
  HR.upper <- round(cox.summary$conf.int[1:nvar, 4], 2)
  HR.range <- sprintf("(%.1f-%.1f)", HR.lower, HR.upper)
  coef <- cox.summary$coefficients[1:nvar, 1]
  coef.se <- cox.summary$coefficients[1:nvar, 3]
  Pval <- round(cox.summary$coefficients[1:nvar, 5], 4)
  signatureA_to_progenitor.sig_ratio_score_cat <- gsub("signatureA_to_progenitor.sig_ratio_score_cat", "", rownames(cox.summary$conf.int)[1:nvar])
  return(data.frame(signatureA_to_progenitor.sig_ratio_score_cat = signatureA_to_progenitor.sig_ratio_score_cat, HR = HR, HR.range = HR.range, coef = coef, coef.se = coef.se, Pval = Pval))
}

HR.pval.df <- do.cox.group(TCGA_signatureA_to_progenitor.sig_ratio_score_surv_all_sub, group.order = c("signatureA_to_progenitor.sig_ratio_Low", "signatureA_to_progenitor.sig_ratio_High"))
annotation_text <- paste("HR:", HR.pval.df$HR, '\n', "p =", HR.pval.df$Pval)

sorted_category_names <- sort(TCGA_signatureA_to_progenitor.sig_ratio_score_surv_all_sub$signatureA_to_progenitor.sig_ratio_score_cat) #order the names of the category ('sig_High' & 'sig_Low') in alphabetical order, so the legend is correctly annotated!
legend_labels <- paste0(unique(sorted_category_names), " (n=", Survival_table$n, ")")
plot <- ggsurvplot(Survival_table, conf.int = FALSE, surv.size = 1.2, xlab = 'Days', ylab = 'Survival', legend.title = "", legend.labs = legend_labels,
                   ylim = c(0, NA), xlim = c(0, 2500), break.x.by = 1000, palette = c('red', 'black'), pval = F, size = 0.25,
                   censor.shape = 124, censor.size = 1) +
  ggtitle(paste("TCGA PAAD: signature A / progenitor signature"))
plot$plot <- plot$plot + guides(color = guide_legend(ncol = 1)) + # make legend in vertical alignment
  annotate("text", x = Inf, y = Inf, label = annotation_text, hjust = 1.1, vjust = 1.1, size = 6*0.3527) # add HR and pval annotation
plot

# EXPORT FIGURE
pdf('Figure_S5M.pdf', width = 1.6, height = 2.2, bg = 'transparent')
plot$plot + theme(axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1),
                  axis.text.x = element_text(size = 6, color = 'black'), axis.text.y = element_text(size = 6, color = 'black'),
                  axis.title.x = element_text(size = 6), axis.title.y = element_text(size = 6), title = element_text(size = 6),
                  legend.text = element_text(size = 6), 
                  plot.title = element_text(size = 6, face = 'plain', margin = margin(0,0,0,0)),
                  legend.background = element_rect(fill = 'transparent', color = NA),
                  plot.background = element_rect(fill = "transparent", color = NA),
                  panel.background = element_rect(fill = "transparent", color = NA))
dev.off()

# Figure S5N----
# UMAP of monocyte/macrophage (MonoMac) cell subsets from all AAV-gRNA tumors
# Bubble plot of scaled differentially expressed genes (DEGs) in MonoMac cell subsets
mmaav_0.4 <- readRDS("R_objects/Seurat_objects/mmaav_0.4.rds")
#Idents(mmaav_0.4) <- 'annotation_1'

png("Figure_S5N_UMAP.png", width = 900, height = 800, res = 300, bg='transparent')
DimPlot(mmaav_0.4, reduction = 'umap', label = T, label.size = 4, repel = F, raster = F, pt.size = 0.1) &
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(size = 12, face = 'plain'),
        plot.margin = margin(0, 2.5, 0, 0, unit = 'pt'), axis.line = element_line(linewidth = 0),
        legend.title = element_blank(), legend.text = element_blank(), legend.position = c(0.96, 0.5), legend.key.height = unit(0.7, 'cm'),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA))
dev.off()

mmaav_0.4_markers <- read_excel('Supplementary_Tables/Supplementary Table 1 DEG scRNAseq.xlsx', sheet = 'MonoMac')
mmaav_0.4_markers_padj0.1 <- mmaav_0.4_markers[which(mmaav_0.4_markers$p_val_adj<0.1),]
mmaav_0.4_markers_padj0.1 <- mmaav_0.4_markers_padj0.1[order(mmaav_0.4_markers_padj0.1$avg_log2FC,decreasing = TRUE),]
mmaav_0.4_markers_padj0.1 <- mmaav_0.4_markers_padj0.1[order(mmaav_0.4_markers_padj0.1$cluster,decreasing = FALSE),]
mmaav_0.4_markers_padj0.1_Top7 <- mmaav_0.4_markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC)

pdf('Figure_S5N_BubblePlot.pdf', width = 5.5, height = 1.4)
DotPlot(mmaav_0.4,
        features = unique(mmaav_0.4_markers_padj0.1_Top7$gene), dot.scale = 1,
        cols = 'RdBu', assay = "RNA") + theme(text = element_text(size = 6), axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5, size = 4), axis.text.y=element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank(),
                                              axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1),
                                              legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.1, 'lines'))
dev.off()

# Figure S5O----
# UMAP of dendritic cell (DC) subsets from all AAV-gRNA tumors. (Middle) UMAP visualization of select DC marker genes in the DC object
# Bubble plot of scaled DEGs in DC subsets
dcs <- readRDS("R_objects/Seurat_objects/dcs_0.2.RData") # dendritic cells

png("Figure_S5O_UMAP.png", width = 800, height = 600, res = 300, bg='transparent')
DimPlot(dcs, reduction = 'umap', label = T, label.size = 4, repel = F, raster = F, pt.size = 0.1) &
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(size = 12, face = 'plain'),
        plot.margin = margin(0, 2.5, 0, 0, unit = 'pt'), axis.line = element_line(linewidth = 0),
        legend.title = element_blank(), legend.text = element_blank(), legend.position = c(0.96, 0.5), legend.key.height = unit(0.7, 'cm'),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA))
dev.off()

dcs_markers <- read_excel('Supplementary_Tables/Supplementary Table 1 DEG scRNAseq.xlsx', sheet = 'Dendritic cells')
dcs_markers_padj0.1 <- dcs_markers[which(dcs_markers$p_val_adj<0.1),]
dcs_markers_padj0.1 <- dcs_markers_padj0.1[order(dcs_markers_padj0.1$avg_log2FC,decreasing = TRUE),]
dcs_markers_padj0.1 <- dcs_markers_padj0.1[order(dcs_markers_padj0.1$cluster,decreasing = FALSE),]
dcs_markers_padj0.1_Top5 <- dcs_markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

pdf('Figure_S5O_BubblePlot.pdf', width = 3.5, height = 1.1)
DotPlot(dcs,
        features = unique(dcs_markers_padj0.1_Top5$gene), dot.scale = 1,
        cols = 'RdBu', assay = "RNA") + theme(text = element_text(size = 6), axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5, size = 4), axis.text.y=element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank(),
                                              axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1),
                                              legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.1, 'lines'))
dev.off()

# Figure S5P----
# UMAP of neutrophil subsets from all AAV-gRNA tumors
# Bubble plot of scaled DEGs in neutrophil subsets
neuts_0.4 <- readRDS("R_objects/Seurat_objects/neuts_0.4.rds") # neutrophils

png("Figure_S5P_UMAP.png", width = 900, height = 800, res = 300, bg='transparent')
DimPlot(neuts_0.4, reduction = 'umap', label = T, label.size = 4, repel = F, raster = F, pt.size = 0.1) &
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(size = 12, face = 'plain'),
        plot.margin = margin(0, 2.5, 0, 0, unit = 'pt'), axis.line = element_line(linewidth = 0),
        legend.title = element_blank(), legend.text = element_blank(), legend.position = c(0.96, 0.5), legend.key.height = unit(0.6, 'cm'),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA))
dev.off()

neuts_0.4_markers <- read_excel('Supplementary_Tables/Supplementary Table 1 DEG scRNAseq.xlsx', sheet = 'Neutrophils')
neuts_0.4_markers_padj0.1 <- neuts_0.4_markers[which(neuts_0.4_markers$p_val_adj<0.1),]
neuts_0.4_markers_padj0.1 <- neuts_0.4_markers_padj0.1[order(neuts_0.4_markers_padj0.1$avg_log2FC,decreasing = TRUE),]
neuts_0.4_markers_padj0.1 <- neuts_0.4_markers_padj0.1[order(neuts_0.4_markers_padj0.1$cluster,decreasing = FALSE),]
neuts_0.4_markers_padj0.1_Top5 <- neuts_0.4_markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

pdf('Figure_S5P_BubblePlot.pdf', width = 5, height = 1.4)
DotPlot(neuts_0.4,
        features = unique(neuts_0.4_markers_padj0.1_Top5$gene), dot.scale = 1,
        cols = 'RdBu', assay = "RNA") + theme(text = element_text(size = 6), axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5, size = 4), axis.text.y=element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank(),
                                              axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1),
                                              legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.1, 'lines'))
dev.off()

# Figure S5Q----
# UMAP of lymphoid cell subsets from all AAV-gRNA tumors
# Bubble plot of scaled DEGs in lymphoid cell subsets
tnk_clean <- readRDS("R_objects/Seurat_objects/tnk_clean_1.rds") # T and NK cells

#Idents(tnk_clean) <- 'cell.type'
png("Figure_S5Q_UMAP.png", width = 900, height = 800, res = 300, bg='transparent')
DimPlot(tnk_clean, reduction = 'umap', label = T, label.size = 4, repel = F, raster = F, pt.size = 0.1) &
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(size = 12, face = 'plain'),
        plot.margin = margin(0, 2.5, 0, 0, unit = 'pt'), axis.line = element_line(linewidth = 0),
        legend.title = element_blank(), legend.text = element_blank(), legend.position = c(0.96, 0.5), legend.key.height = unit(0.8, 'cm'),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA))
dev.off()

tnk_clean_markers <- read_excel('Supplementary_Tables/Supplementary Table 1 DEG scRNAseq.xlsx', sheet = 'T NK cells')
tnk_clean_markers_padj0.1 <- tnk_clean_markers[which(tnk_clean_markers$p_val_adj<0.1),]
tnk_clean_markers_padj0.1 <- tnk_clean_markers_padj0.1[order(tnk_clean_markers_padj0.1$avg_log2FC,decreasing = TRUE),]
tnk_clean_markers_padj0.1 <- tnk_clean_markers_padj0.1[order(tnk_clean_markers_padj0.1$cluster,decreasing = FALSE),]
tnk_clean_markers_padj0.1_Top10 <- tnk_clean_markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

pdf('Figure_S5Q_BubblePlot.pdf', width = 5, height = 1)
DotPlot(tnk_clean,
        features = unique(tnk_clean_markers_padj0.1_Top10$gene), dot.scale = 1,
        cols = 'RdBu', assay = "RNA") + theme(text = element_text(size = 6), axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5, size = 4), axis.text.y=element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank(),
                                              axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1),
                                              legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.1, 'lines'))
dev.off()

# Figure S5R----
# Violin plots of H2-Ab1, H2-Aa, H2-Eb1 expression in the TAM.Ndrg1 subset, grouped by AAV-gRNA tumors
ndrg1 <- subset(mmaav_0.4, idents = 'TAM.Ndrg1')
pdf('Figure_S5R.pdf', width = 1, height = 2.3, bg = 'transparent')
gg1 <- VlnPlot(ndrg1, 'H2-Ab1', group.by = 'aav', pt.size = 0, ncol = 1) + NoLegend() +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 6), axis.text = element_text(size = 6), axis.text.y = element_blank(),
        axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1),
        plot.title = element_text(size = 6, face = 'plain', margin = margin(0,0,0,0)),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA))
gg1$layers[[1]]$aes_params$size <- 0.1 #adjust linewidth of violins
gg2 <- VlnPlot(ndrg1, 'H2-Aa', group.by = 'aav', pt.size = 0, ncol = 1) + NoLegend() +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 6), axis.text = element_text(size = 6), axis.text.y = element_blank(),
        axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1),
        plot.title = element_text(size = 6, face = 'plain', margin = margin(0,0,0,0)),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA))
gg2$layers[[1]]$aes_params$size <- 0.1 #adjust linewidth of violins
gg3 <- VlnPlot(ndrg1, 'H2-Eb1', group.by = 'aav', pt.size = 0, ncol = 1) + NoLegend() +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 6), axis.text = element_text(size = 6), axis.text.y = element_blank(),
        axis.line = element_line(linewidth = 0.1), axis.ticks = element_line(linewidth = 0.1),
        plot.title = element_text(size = 6, face = 'plain', margin = margin(0,0,0,0)),
        plot.margin = margin(0, 0, 0, 0, unit = 'pt'),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA))
gg3$layers[[1]]$aes_params$size <- 0.1 #adjust linewidth of violins
gg1 / gg2 / gg3
dev.off()
