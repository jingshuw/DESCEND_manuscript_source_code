library(Seurat)

result.all <- readRDS("descend_result.rds")
data <- readRDS("sampled_combined_pure_data.rds")
Y <- as.matrix(data$Y)
colnames(Y) <- 1:ncol(Y)
library(descend)
hvg <- findHVG(result.all, spline.df = 15, quantile = 0.5, threshold = 8)
x.seurat <- CreateSeuratObject(raw.data = Y)
x.seurat <- NormalizeData(x.seurat)
x.seurat <- FindVariableGenes(x.seurat, do.plot = F)
genes1 <- hvg$HVG.genes
genes2 <- x.seurat@var.genes

require(VennDiagram)
venn.list <- list(genes1, genes2)
names(venn.list) <- c("DESCEND", 
                      "Seurat")
venn.plot <- venn.diagram(venn.list,
                          NULL, 
                          fill=c("red", "blue"), 
                          alpha=c(0.5,0.5), cex = 2.5,
                          cat.pos = 6, cat.cex = 2.5, 
                          cat.just = list(c(0.9, 0.1), c(0.2, 0.1))
                          )
grid.newpage()
grid.draw(venn.plot)

## Sample clustering analysis for number of PC = 30 ##
x.seurat@var.genes <- hvg$HVG.genes
meta <- data.frame(Cell_Type = data$labels)
rownames(meta) <- colnames(Y)
x.seurat <- AddMetaData(x.seurat, metadata = meta)

print(paste("# of HVGs: ", length(x.seurat@var.genes)))
x.seurat <- ScaleData(x.seurat)
x.seurat <- RunPCA(x.seurat, pcs.compute = 30)

## the resolution is chosen for Seurat to find the right number of clusters (here is 10)
## change the resolution and dims.use for a different chosen PC
resolution <- 2
x.seurat <- FindClusters(object = x.seurat, 
                         #    genes.use = x.seurat@var.genes, 
                         reduction.type = "pca", 
                         resolution = resolution, 
                         #    save.SNN = TRUE, 
                         dims.use = 1:30,
                         print.output = 0, save.SNN = T)
library(mclust)
adjustedRandIndex(x.seurat@meta.data[, 5], labels)



