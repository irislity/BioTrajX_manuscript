# =============================================================================
# GSE131847_preprocessing.R
#
# Preprocessing of the GSE131847 LCMV time-course CD8 T-cell dataset
# (naive, d3, d4, d5, d6, d7, d10, d14, d21, d32, d60, d90 post-infection):
# load raw per-day UMI tables, build/normalize a Seurat object, run
# PCA/UMAP, fit Slingshot trajectories, and reproduce the marker-gene,
# pseudotime-heatmap, and GSEA panels from the original analysis.
#
# SOURCE: This script is adapted from
# trajectory_analysis.R, published alongside the "Technology meets TILs:  Deciphering T cell function in the -omics era (Hudson, Wieland, 2023)" at
# Mendeley Data: https://doi.org/10.17632/3dvt79c7yt.1
#
# INPUT (expected in the working directory, as in the original analysis):
#   GSE131847_RAW/GSM3822179_Naive_cell-gene_UMI_table.tsv, ... (per-day UMI tables)
#   GSE41867_MEMORY_VS_EXHAUSTED_CD8_TCELL_DAY30_LCMV_DN.v7.5.1.grp
#   GSE41867_MEMORY_VS_EXHAUSTED_CD8_TCELL_DAY30_LCMV_UP.v7.5.1.grp
#
# OUTPUT:
#   dimplot.pdf, dimplot_no_legend.pdf, trajectories.pdf, selected_genes.pdf,
#   psuedotime_gene_expression.pdf, heatmap_legend.pdf,
#   realtime_gene_expression.pdf, geneset_plot.pdf, running_plot.pdf
#
# Usage:
#   Rscript manuscript/scripts/real/GSE131847_preprocessing.R
# =============================================================================

# ── Load raw per-day UMI count tables ───────────────────────────────────────
naive <- read.table(file = "GSE131847_RAW/GSM3822179_Naive_cell-gene_UMI_table.tsv", header = T, row.names = 1)
d3    <- read.table(file = "GSE131847_RAW/GSM3822180_D3_cell-gene_UMI_table.tsv",    header = T, row.names = 1)
d4    <- read.table(file = "GSE131847_RAW/GSM3822181_D4_cell-gene_UMI_table.tsv",    header = T, row.names = 1)
d5    <- read.table(file = "GSE131847_RAW/GSM3822183_D5_cell-gene_UMI_table.tsv",    header = T, row.names = 1)
d6    <- read.table(file = "GSE131847_RAW/GSM3822184_D6_cell-gene_UMI_table.tsv",    header = T, row.names = 1)
d7    <- read.table(file = "GSE131847_RAW/GSM3822185_D7_cell-gene_UMI_table.tsv",    header = T, row.names = 1)
d10   <- rbind(read.table(file = "GSE131847_RAW/GSM3822188_D10_cell-gene_UMI_table.tsv", header = T, row.names = 1),
               read.table(file = "GSE131847_RAW/GSM3822189_D10_2_cell-gene_UMI_table.tsv", header = T, row.names = 1))
d14   <- rbind(read.table(file = "GSE131847_RAW/GSM3822191_D14_cell-gene_UMI_table.tsv", header = T, row.names = 1),
               read.table(file = "GSE131847_RAW/GSM3822192_D14_2_cell-gene_UMI_table.tsv", header = T, row.names = 1))
d21   <- read.table(file = "GSE131847_RAW/GSM3822194_D21_cell-gene_UMI_table.tsv", header = T, row.names = 1)
d32   <- read.table(file = "GSE131847_RAW/GSM3822197_D32_cell-gene_UMI_table.tsv", header = T, row.names = 1)
d60   <- read.table(file = "GSE131847_RAW/GSM3822200_D60_cell-gene_UMI_table.tsv", header = T, row.names = 1)
d90   <- read.table(file = "GSE131847_RAW/GSM3822202_D90_cell-gene_UMI_table.tsv", header = T, row.names = 1)

# make a data frame with cell names and which day from which they originate
metadata <- data.frame(row.names = c(rownames(naive),
                                      rownames(d3),
                                      rownames(d4),
                                      rownames(d5),
                                      rownames(d6),
                                      rownames(d7),
                                      rownames(d10),
                                      rownames(d14),
                                      rownames(d21),
                                      rownames(d32),
                                      rownames(d60),
                                      rownames(d90)),
                        sample = c(rep("naive", length(rownames(naive))),
                                   rep("d3",    length(rownames(d3))),
                                   rep("d4",    length(rownames(d4))),
                                   rep("d5",    length(rownames(d5))),
                                   rep("d6",    length(rownames(d6))),
                                   rep("d7",    length(rownames(d7))),
                                   rep("d10",   length(rownames(d10))),
                                   rep("d14",   length(rownames(d14))),
                                   rep("d21",   length(rownames(d21))),
                                   rep("d32",   length(rownames(d32))),
                                   rep("d60",   length(rownames(d60))),
                                   rep("d90",   length(rownames(d90)))))

# combine counts
alldata <- rbind(naive, d3, d4, d5, d6, d7, d10, d14, d21, d32, d60, d90)
dim(alldata) # 26,261 cells; 9,941 genes

# ── Build Seurat object ──────────────────────────────────────────────────────
library(Seurat)
seurat_object <- CreateSeuratObject(t(alldata))
seurat_object <- SCTransform(seurat_object, vst.flavor = "v1") # normalize and scale

# add day of isolation to seurat object
seurat_object <- AddMetaData(seurat_object, factor(metadata$sample, levels = c("naive", "d3", "d4", "d5", "d6", "d7", "d10", "d14", "d21", "d32", "d60", "d90")), col.name = "cell_type")

# Dimensionality reduction
seurat_object <- RunPCA(seurat_object)          # 50 principal components is the default
seurat_object <- RunUMAP(seurat_object, dims = 1:50) # using default components

# Make colors for days
library(RColorBrewer)
library(ggplot2)
colors <- colorRampPalette(brewer.pal(8, "Spectral"))(12)
names(colors) <- unique(seurat_object$cell_type)
cell_colors <- colors[seurat_object$cell_type]
names(cell_colors) <- colnames(seurat_object)

DimPlot(seurat_object, reduction = "umap",
        group.by = "cell_type", pt.size = 0.5, label = TRUE, repel = TRUE, cols = colors)
ggsave("dimplot.pdf", height = 5, width = 5)

DimPlot(seurat_object, reduction = "umap",
        group.by = "cell_type", pt.size = 0.5, label = TRUE, repel = TRUE, cols = colors) +
  theme(legend.position = "None")
ggsave("dimplot_no_legend.pdf", height = 5, width = 5)

# ── Slingshot trajectories ───────────────────────────────────────────────────
library(slingshot)
library(cowplot)
library(SingleCellExperiment)

slingshot_output <- slingshot(Embeddings(seurat_object, "umap"), clusterLabels = seurat_object$cell_type)

slingshot_output_1 <- slingshot(Embeddings(seurat_object, "umap"), clusterLabels = seurat_object$cell_type, end.clus = "d90")

slingshot_output_2 <- slingshot(Embeddings(seurat_object, "umap"), clusterLabels = seurat_object$cell_type, start.clus = "naive", end.clus = "d90")

slingshot_output_pt <- slingPseudotime(slingshot_output)

slingshot_output # 2 curves

# make a data frame with all together
for_plot <- cbind(as.data.frame(seurat_object[["umap"]]@cell.embeddings),
                   as.data.frame(slingshot_output_pt[colnames(seurat_object), ]))

curve1 <- as.data.frame(slingCurves(slingshot_output)[[1]]$s[slingCurves(slingshot_output)[[1]]$ord, ])
curve2 <- as.data.frame(slingCurves(slingshot_output)[[2]]$s[slingCurves(slingshot_output)[[2]]$ord, ])

for_plot$cluster <- seurat_object$cell_type[rownames(for_plot)]

ggplot(for_plot) +
  geom_point(aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_path(data = curve1, aes(x = umap_1, y = umap_2), color = "black", size = 1.2,
            arrow = arrow(
              angle = 25,                # arrow head angle
              length = unit(0.25, "cm"), # arrow size
              ends = "last",             # put arrow at the end of the curve
              type = "open")) +
  # Curve 2
  geom_path(data = curve2, aes(x = umap_1, y = umap_2), color = "red", size = 1.2,
            arrow = arrow(
              angle = 25,                # arrow head angle
              length = unit(0.25, "cm"), # arrow size
              ends = "last",             # put arrow at the end of the curve
              type = "open")) +
  scale_color_manual(values = colors) +
  theme_cowplot() +
  theme(legend.position = "None")

ggsave("trajectories.pdf", height = 5, width = 5)

# ── Marker-gene violin plots and pseudotime/realtime heatmaps ───────────────
# Make violin plots of various genes over real time
plot_grid(
  VlnPlot(seurat_object, features = c("Sell"),  group.by = "cell_type", cols = colors, pt.size = 0) + theme(legend.position = "None") + theme(axis.title = element_blank()) + theme(plot.title = element_blank()),
  VlnPlot(seurat_object, features = c("Il7r"),  group.by = "cell_type", cols = colors, pt.size = 0) + theme(legend.position = "None") + theme(axis.title = element_blank()) + theme(plot.title = element_blank()),
  VlnPlot(seurat_object, features = c("Tcf7"),  group.by = "cell_type", cols = colors, pt.size = 0) + theme(legend.position = "None") + theme(axis.title = element_blank()) + theme(plot.title = element_blank()),
  VlnPlot(seurat_object, features = c("Mki67"), group.by = "cell_type", cols = colors, pt.size = 0) + theme(legend.position = "None") + theme(axis.title = element_blank()) + theme(plot.title = element_blank()),
  VlnPlot(seurat_object, features = c("Il2ra"), group.by = "cell_type", cols = colors, pt.size = 0) + theme(legend.position = "None") + theme(axis.title = element_blank()) + theme(plot.title = element_blank()),
  VlnPlot(seurat_object, features = c("Gzmb"),  group.by = "cell_type", cols = colors, pt.size = 0) + theme(legend.position = "None") + theme(axis.title = element_blank()) + theme(plot.title = element_blank()))

ggsave("selected_genes.pdf", height = 260 / 72, width = 580 / 72, units = "in")

# Show gene expression over pseudotime
orders <- for_plot[order(for_plot$curve1, decreasing = F), ]
genes  <- c("Sell", "Il7r", "Tcf7", "Mki67", "Il2ra", "Gzmb")
orders <- cbind(orders, t(GetAssayData(seurat_object, "scale.data")[genes, rownames(orders)]))
orders$cell_order <- 1:dim(orders)[1]
orders_melt <- reshape2::melt(orders, id = c("cell_order", "UMAP_1", "UMAP_2", "curve1", "#curve2", "curve1_x", "curve1_y", "#curve2_x", "#curve2_y", "cluster"))
orders_melt$value[orders_melt$value > 3] <- 3
orders_melt$value[orders_melt$value < -3] <- -3

plot_a <- ggplot(orders_melt, aes(x = cell_order, y = variable, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "black", high = "red", midpoint = 0) +
  theme_cowplot() +
  theme(legend.position = "None")

plot_b <- ggplot(orders_melt, aes(x = cell_order, y = 1, fill = cluster)) +
  geom_tile() +
  scale_fill_manual(values = colors) +
  theme_cowplot() +
  theme(legend.position = "None")

plot_grid(plot_a, plot_b, ncol = 1, align = "h")

ggsave("psuedotime_gene_expression.pdf", height = 8, width = 11, units = "in")

plot_a + theme(legend.position = "right")
ggsave("heatmap_legend.pdf")

# Show gene expression over realtime
orders_realtime <- for_plot
orders_realtime <- cbind(orders_realtime, t(GetAssayData(seurat_object, "scale.data")[genes, rownames(orders_realtime)]))
orders_realtime$cell_order <- 1:dim(orders_realtime)[1]
orders_realtime_melt <- reshape2::melt(orders_realtime, id = c("cell_order", "UMAP_1", "UMAP_2", "curve1", "#curve2", "curve1_x", "curve1_y", "#curve2_x", "#curve2_y", "cluster"))
orders_realtime_melt$value[orders_realtime_melt$value > 3] <- 3
orders_realtime_melt$value[orders_realtime_melt$value < -3] <- -3

plot_a <- ggplot(orders_realtime_melt, aes(x = cell_order, y = variable, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "black", high = "red", midpoint = 0) +
  theme_cowplot() +
  theme(legend.position = "None")

plot_b <- ggplot(orders_realtime_melt, aes(x = cell_order, y = 1, fill = cluster)) +
  geom_tile() +
  scale_fill_manual(values = colors) +
  theme_cowplot() +
  theme(legend.position = "None")

plot_grid(plot_a, plot_b, ncol = 1, align = "h")

ggsave("realtime_gene_expression.pdf", height = 8, width = 11, units = "in")

# ── GSEA ──────────────────────────────────────────────────────────────────────
# Geneset
genesets <- as.list(c(read.table("GSE41867_MEMORY_VS_EXHAUSTED_CD8_TCELL_DAY30_LCMV_DN.v7.5.1.grp", header = TRUE),
                       read.table("GSE41867_MEMORY_VS_EXHAUSTED_CD8_TCELL_DAY30_LCMV_UP.v7.5.1.grp", header = TRUE)))

# compare d6 with d90 - using log2FC as ranking statistic
d6vs90 <- FindMarkers(object = seurat_object, ident.1 = "d6", ident.2 = "d90", group.by = "cell_type", logfc.threshold = 0)
rank   <- data.frame(row.names = rownames(d6vs90), stat = d6vs90$avg_log2FC)
rank   <- rank[order(rank$stat, decreasing = T), , drop = F] # order from highest to lowest; positive = higher in d6
ranks  <- as.numeric(rank$stat)
names(ranks) <- toupper(row.names(rank))

# do GSEA analysis
library("fgsea")
set.seed(123)
gsea_results <- fgsea(genesets, ranks, 1000)
gsea_results # p=0.036 for GSE41867_MEMORY_VS_EXHAUSTED_CD8_TCELL_DAY30_LCMV_DN; will focus on this one

# plot of avg. expression
genes <- stringr::str_to_title(genesets[[1]]) # will plot genes higher in exhausted cells
genes <- genes[genes %in% rownames(seurat_object[["SCT"]])] # use only those genes expressed in this dataset
seurat_object <- AddMetaData(seurat_object, scale(colSums(seurat_object@assays$SCT@counts[genes, ])), "genelist") # scaled expression of sum of all genes in geneset

# values above/below 2.5/-2.5 capped
for_plot$genelist <- seurat_object[, rownames(for_plot)]$genelist
for_plot[for_plot$genelist > 2.5, ]$genelist <- 2.5
for_plot[for_plot$genelist < -2.5, ]$genelist <- -2.5 # none found; causes error

geneset_plot <- ggplot(for_plot) +
  geom_point(shape = 21, aes(x = UMAP_1, y = UMAP_2, fill = genelist), stroke = 0, size = 1.2) +
  theme_cowplot() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-2.5, 2.5), breaks = c(-2.5, -1.25, 0, 1.25, 2.5))
geneset_plot
ggsave("geneset_plot.pdf", height = 5, width = 5)

# make a running plot
pathway <- genesets[1]
to.plot <- data.frame(x = numeric(0), y = numeric(0))

gseaParam <- 1
rnk <- rank(-ranks) # make ranking list
ord <- order(rnk)   # make order
statsAdj <- ranks[ord] # order by rank
statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
statsAdj <- statsAdj / max(abs(statsAdj))

pathway2 <- unname(as.vector(na.omit(match(pathway[[1]], names(statsAdj)))))
pathway2 <- sort(pathway2) # where the genes are in the list

gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway2,
                         returnAllExtremes = TRUE)

bottoms <- gseaRes$bottoms
tops    <- gseaRes$tops

to.plot <- rbind(to.plot, data.frame(x = pathway2, y = tops))

running.plot <- ggplot(data = to.plot, aes(x = x, y = y)) +
  xlim(0, length(ranks)) +
  ylim(-1, 1) +
  theme_cowplot() +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = colors) +
  geom_rug(data = to.plot, sides = "b", color = colorRampPalette(c("#FCA85E", "#3288BD"))(length(ranks))[to.plot$x]) +
  scale_color_manual(values = colors) +
  coord_fixed(1500) +
  geom_point(alpha = 1, shape = 21, color = "black", stroke = .2, size = 3, fill = colorRampPalette(c("#FCA85E", "#3288BD"))(length(ranks))[to.plot$x]) +
  geom_vline(xintercept = length(ranks)) # to show cutoff
running.plot
ggsave("running_plot.pdf", height = 3, width = 4)
