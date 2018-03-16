results <- readRDS("all_cells_ests_forHVG.rds")
library(descend)
sel.gene.names <- findHVG(result, quantile = 0.5, threshold = 10, spline.df = 15)


library(Seurat)

## need to specifiy the folder for the downloaded data, trueMol counts for each ERCC spike-ins and calculate 
## cell specific efficiency beforehand. See functions_all.R file
data <- GetY(DATAFILEPATH, trueMol, alpha)

for (i in 1:10) gc()

rownames(data$Y)[rownames(data$Y) == "2-Mar"] <- c("2-Mar", "2-Mar_1") 

set.seed(1)
x <- data$Y
colnames(x) <- 1:ncol(x)

x.seurat <- CreateSeuratObject(raw.data = x)
x.seurat <- NormalizeData(x.seurat)
x.seurat <- FindVariableGenes(x.seurat, do.plot = F)

genes1 <- sel.gene.names$HVG.genes
genes2 <- x.seurat@var.genes

### plot the Venn diagram #########
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

## the clustering result is obtained by running Seurat and the code is omitted here
