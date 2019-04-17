require(remotes)
remotes::install_github("UCSF-TI/fake-hdf5r") #fixes error in importing hdf5r package
require(Seurat)
library(mclust)
library(dplyr)

seurat <- CreateSeuratObject(
  raw.data = counts(sce_sc_10x_qc),
  min.cells = 3, 
  min.genes = 200
)
# violin plots -- useful for visualization in write up
VlnPlot(object = seurat, features.plot = c("nGene", "nUMI"), nCol = 2)
# gene plot
GenePlot(object = seurat, gene1 = "nUMI", gene2 = "nGene")
# no filtering needed
# PCA
seurat <- NormalizeData(seurat)
seurat <- FindVariableGenes(seurat,
                            mean.function = ExpMean,
                            dispersion.function = LogVMR, 
                            x.low.cutoff = 0.0125, 
                            x.high.cutoff = 3, 
                            y.cutoff = 0.5)


seurat <- ScaleData(object = seurat, vars.to.regress = c("nUMI"))
seurat <- RunPCA(object = seurat, pc.genes = seurat@var.genes,  do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCElbowPlot(seurat)

# clustering
