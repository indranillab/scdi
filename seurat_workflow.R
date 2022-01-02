################################################################################
### Alignment workflow for Immune Dataset from Kang et al.
################################################################################

library(Seurat)
ctrl.data <- data1
stim.data <- data2

# Set up control object
ctrl <- CreateSeuratObject(ctrl.data, project = "Data 1", min.cells = 5)
ctrl@meta.data$stim <- "Data1"
ctrl <- FilterCells(ctrl, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
ctrl <- NormalizeData(ctrl)
ctrl <- ScaleData(ctrl)
# Set up stimulated object
stim <- CreateSeuratObject( stim.data, project = "Data 2", min.cells = 5)
stim@meta.data$stim <- "Data2"
stim <- FilterCells(stim, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
stim <- NormalizeData(stim)
stim <- ScaleData(stim)

# Gene selection for input to CCA
ctrl <- FindVariableGenes(ctrl, do.plot = F)
stim <- FindVariableGenes(stim, do.plot = F)
g.1 <- head(rownames(ctrl@hvg.info), 1000)
g.2 <- head(rownames(stim@hvg.info), 1000)
genes.use <- unique(c(g.1,g.2))
genes.use <- intersect(genes.use, rownames(ctrl@scale.data))
genes.use <- intersect(genes.use, rownames(stim@scale.data))

# Run CCA
immune.combined <- RunCCA(ctrl, stim, genes.use = genes.use, num.cc = 30)

# CC Selection
MetageneBicorPlot(immune.combined, grouping.var = "stim", dims.eval = 1:30)

# Run rare non-overlapping filtering (no cells detected here as expected)
immune.combined <- CalcVarExpRatio(immune.combined, reduction.type = "pca", grouping.var = "orig.ident", dims.use = 1:20)

# Alignment
immune.combined <- AlignSubspace(immune.combined, reduction.type = "cca", grouping.var = "orig.ident", dims.align = 1:20)

# t-SNE and Clustering
immune.combined <- RunTSNE(immune.combined, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T)
immune.combined <- FindClusters(immune.combined, reduction.type = "cca.aligned", resolution = 0.6,
                                dims.use = 1:20, k.param = 30, save.SNN = T, algorithm = 3, force.recalc = TRUE)

# Visualization
TSNEPlot(immune.combined, do.label = T)
