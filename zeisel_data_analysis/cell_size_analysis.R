source("functions_all.R")
results <- readRDS("cell_size_multi_cell_type_descend_results.rds")
result.list <- results$descend.list.list
library(descend)

name.type <- c("Astrocytes-Ependymal", "Endothelial-Mural", "Interneurons",
               "Microglia", "Oligodendrocytes", "CA1 Pyramidal", "S1 Pyramidal", "")

col1 <- rgb(0, 0, 0, alpha = 0.6)
col2 <- rgb(1, 0, 0, alpha = 0.6)

all.results <- lapply(1:7, function(i) {
               result <- result.list[[i]]
               ests.all <- getEstimates(result)
               ests <- sapply(ests.all, function(mat)mat[, 1])
               sds <- sapply(ests.all, function(mat)mat[, 2])
               pvals <- getPval(result)
               return(list(ests = ests, sds = sds, pvals = pvals))
               })


names(all.results) <- name.type[1:7]

common.genes <- getCommonGenes(all.results)



############ Figure for cell size on Nonzero Mean for 6 cell types ######################
layout(matrix(1:6, nrow = 3))
par(mar = c(2, 2, 2, 1), mgp = c(2, 0.5, 0))
for (i in c(1, 3:7)) sizeEffect(i)

############ Figure for cell size on Nonzero Fraction for 6 cell types ######################
par(mar = c(2, 2, 2, 1), mgp = c(2, 0.5, 0))
for (i in c(1, 3:7)) freqEffect(i)

## how the GLM coefficients for the FISH are computed can be found in file ../RNA_FISH_analysis/compute_fish.R
fishes <- read.table("fish_coefs.txt", header = T)
### Effect of cell size on Nonzero Fraction for Endothelial-Mural cells with FISH ###############
par(mar = c(2, 2, 2, 0), mgp = c(2, 0.5, 0))
freqEffect(2, add.density = T, fish.values = fishes[, 1])
### Effect of cell size on Nonzero Fraction for Endothelial-Mural cells with FISH ###############
par(mar = c(2, 2, 2, 0), mgp = c(2, 0.5, 0))
sizeEffect(2, add.density = T, fish.values = fishes[, 2])


############ Violin plots for nonzero fraction before cell size adjustment ###############
temp.data1 <- data.frame(sapply(1:7, 
                                function(i) all.results[[i]]$ests[common.genes, 1]))
library(ggplot2)
temp.data1 <- reshape(temp.data1, 
                      direction = "long", varying = colnames(temp.data1),
                     v.name = "value")
temp.data1$time <- factor(temp.data1$time, labels = name.type[1:7])

temp.data1$Type <- as.factor(temp.data1$time %in% name.type[c(1, 2, 4:5)])
cols <- rep("Black", 7)
cols[c(2, 6)] <- "RED"
levels(temp.data1$Type) <- c("Neurons", "Non-neurons")
p1 <- ggplot(temp.data1, aes(x = time, y = value, fill = Type)) + 
     geom_violin(position = position_dodge(0.5), trim = T, 
                 draw_quantiles = 0.5, size = 0.3) + 
     theme_minimal() + 
     scale_y_continuous(breaks = c(0.3, 0.7, 1), limits = c(0.12, 1.03)) + 
     scale_fill_manual(values = c("#E69F00", "#56B4E9")) + 
     ggtitle("Active Fraction Estimated by DESCEND \n(before cell size adjustment)") +
     theme(axis.text.x = element_text(size = 13, angle = 60, 
                                    hjust = 1, vjust = 1, colour = cols), 
           axis.text.y = element_text(size = 14), 
           plot.title = element_text(size = 16, lineheight = 1),
           axis.title = element_blank(),
           legend.text = element_text(size = 13),
           legend.title = element_blank(), legend.position = "bottom",
        #   legend.margin = unit(-2, "line"),
           plot.margin = unit(x = c(0.2, 0.1, 0.1, 0), units = "cm")) + 
    geom_rect(aes(xmin = 2 - 0.5, xmax = 2 + 0.5, ymin = 0.15, ymax = 1.02),
          fill = "transparent", color = "red", size = 0.6) + 
    geom_rect(aes(xmin = 6 - 0.5, xmax = 6 + 0.5, ymin = 0.15, ymax = 1.02),
          fill = "transparent", color = "red", size = 0.6)
p1


############ Violin plots for nonzero fraction before cell size adjustment ############### 
temp.data2 <- data.frame(sapply(1:7, 
                                function(i) all.results[[i]]$ests[common.genes, 6]))
temp.data2 <- reshape(temp.data2, 
                      direction = "long", varying = colnames(temp.data2),
                     v.name = "value")
temp.data2$time <- factor(temp.data2$time, labels = name.type[1:7])

temp.data2$Type <- as.factor(temp.data2$time %in% name.type[c(1, 2, 4:5)])
cols <- rep("Black", 7)
cols[c(2, 6)] <- "RED"
levels(temp.data2$Type) <- c("Neurons", "Non-neurons")
p2 <- ggplot(temp.data2, aes(x = time, y = value, fill = Type)) + 
     geom_violin(position = position_dodge(0.5), trim = T,
                 draw_quantiles = 0.5, size = 0.3) + 
     theme_minimal() + 
     ggtitle("Active Fraction Estimated by DESCEND \n(after cell size adjustment)") +
   scale_y_continuous(breaks = c(0.3, 0.7, 1), limits = c(0.12, 1.03)) + 
     scale_fill_manual(values = c("#E69F00", "#56B4E9")) + 
     theme(axis.text.x = element_text(size = 13, angle = 60, 
                                    hjust = 1, vjust = 1, colour = cols), 
           axis.text.y = element_text(size = 14), 
           plot.title = element_text(size = 16, lineheight = 1),
           axis.title = element_blank(),
           legend.text = element_text(size = 13),
           legend.title = element_blank(), legend.position = "bottom",
plot.margin = unit(x = c(0.2, 0.1, 0.1, 0), units = "cm")) + 
    geom_rect(aes(xmin = 2 - 0.5, xmax = 2 + 0.5, ymin = 0.15, ymax = 1.02),
          fill = "transparent", color = "red", size = 0.6) + 
    geom_rect(aes(xmin = 6 - 0.5, xmax = 6 + 0.5, ymin = 0.15, ymax = 1.02),
          fill = "transparent", color = "red", size = 0.6)
p2

#################### Figures for Differential testing  ########################
DE.result <- readRDS("DEtest_celltypes_2_6_descend_results.rds")
diff.test <- getGeneList(DE.result)

common.genes26 <- getCommonGenes(all.results[c(2, 6)])

diff1 <- all.results[[2]]$ests[common.genes26, 1] - all.results[[6]]$ests[common.genes26, 1]
diff2<- all.results[[2]]$ests[common.genes26, 6] - all.results[[6]]$ests[common.genes26, 6]

### Figure of DETest on nonzero fraction ###
#pdf("~/Dropbox/sparse_factor_bic/notes/g_model/plots/zeisel_diff_new.pdf", width = 4.5, height = 4.5)
idx <- sample(1:length(diff1), 1500)
require(scales)
cols <- hue_pal()(4)[c(2, 1, 4)]
cols[1] <- hue_pal()(3)[2]
par(mar = c(3, 4, 3, 1), mgp = c(2, 0.5, 0))
plot(diff1[idx], diff2[idx], pch = 18, col = "gray",
     xlab = "BEFORE cell size adjustment",
     ylab = "AFTER cell size adjustment",
     main = "Difference in Nonzero Fraction", cex.lab = 1.4, cex.main = 1.5,
     xaxt = "n", yaxt = "n", font.main = 1, xlim = c(-0.8, 0.9), ylim = c(-0.8, 0.9))
abline(a = 0, b= 1, lty= 2)
idx1 <- intersect(names(which(diff.test$judge.comb[, 1])),  
                  names(which(!diff.test$judge.comb[, 2])))
points(diff1[idx1], diff2[idx1], col = cols[3], pch = 18)
idx1 <- intersect(names(which(diff.test$judge.comb[, 1])),  
                  names(which(diff.test$judge.comb[, 2])))
points(diff1[idx1], diff2[idx1], col = cols[1], pch = 18)
idx1 <- intersect(names(which(!diff.test$judge.comb[, 1])),  
                  names(which(diff.test$judge.comb[, 2])))
points(diff1[idx1], diff2[idx1], col = cols[2], pch = 18, cex = 1.2)
legend("topleft", bty = "n", pch = rep(18, 3), 
       col = cols[c(3, 2, 1)], pt.cex = 1.5,
       legend = c("Only before adjustment", "Only after adjustment",
                  "Significant for both"), cex = 1.15)
axis(2, at = c(-0.5, 0, 0.5), tck = -0.02, cex.axis = 1.4, las = 0)
axis(1, at = c(-0.5, 0, 0.5), tck = -0.02, cex.axis = 1.4, las = 1)


#### Venn diagram of selected genes before VS after cell size adjustment
require(VennDiagram)
venn.list <- apply(diff.test$judge.comb, 2, function(v) which(v))
names(venn.list) <- c("Before Cell Size Adjustment", 
                      "After Cell Size \n Adjustment")
venn.plot <- venn.diagram(venn.list,
                          NULL, 
                          fill=c("red", "blue"), 
                          alpha=c(0.25,0.25), cex = 2,
                          cat.pos = 6, cat.cex = 2.0, 
                          cat.just = list(c(0.5, 0.8), c(0.7, 0.2)))
grid.newpage()
grid.draw(venn.plot)


#### Figure for difference of Nonzero mean VS difference of Nonzero fraction
diff.gene.list.inten <- names(which(diff.test$judge.mean))


diff.genes <- rownames(diff.test$judge.comb)
temp1 <- log(all.results[[2]]$ests[diff.genes, 2]) 
temp2 <- log(all.results[[6]]$ests[diff.genes, 2])
sig.genes <- intersect(diff.genes, diff.gene.list.inten) 
fold.change <- temp1 - temp2
diff <- diff2[diff.genes]
avg <- (temp1 + temp2)/2


par(mar = c(4, 5.5, 1, 1), mgp = c(3, 0.5, 0))
smoothScatter(diff, fold.change, ylab = "logFC of\nNonzero Mean", 
              xlab = "Difference of\nNonzero Fraction",
              xaxt = "n", yaxt = "n",
     cex.axis = 1.2, cex.lab = 1.3)
axis(1, at = c(-0.5, 0, 0.5), tck = -0.02, cex.axis = 1.4)
axis(2, at = c(-2, 1, 4), tck = -0.02, cex.axis = 1.4, las = 1)


## Venn Diagram for nonzero fraction VS nonzero mean
require(VennDiagram)
venn.list2 <- list(which(diff.test$judge.comb[, 2]), diff.gene.list.inten)
names(venn.list2) <- c("Nonzero Fraction", "Nonzero Intensity")
venn.plot <- venn.diagram(venn.list2,
                          NULL, 
                          fill=c("red", "blue"), 
                          alpha=c(0.25,0.25), cex = 2.1,
                          cat.pos = 6, cat.cex = 2.3, 
                          cat.just = list(c(0.7, 0.5), c(0.5, 0)))
grid.newpage()
grid.draw(venn.plot)


## Nonzero mean MA plot
temp1 <- log(all.results[[2]]$ests[common.genes26, 2]) 
temp2 <- log(all.results[[6]]$ests[common.genes26, 2])
sig.genes <- intersect(common.genes26, diff.gene.list.inten) 
fold.change <- temp1 - temp2
avg <- (temp1 + temp2)/2
## down-sample the black dots to reduce figure size
idx <- sample.int(length(avg), 800)
par(mar = c(4, 4.5, 1, 1), mgp = c(3, 0.5, 0))

plot(exp(avg)[idx], fold.change[idx], log = "x", ylab = "logFC", 
     xlab = "Cell Size Adjusted\nNonzero Mean",
     xaxt = "n", yaxt = "n",
     pch = 20, cex.axis = 1.2, cex.lab = 1.3, col = col1)
points(exp(avg[sig.genes]), fold.change[sig.genes], col = col2, pch = 20)
axis(1, at = c(1, 5, 20), tck = -0.02, cex.axis = 1.4)
axis(2, at = c(-2, 0, 3), tck = -0.02, cex.axis = 1.4, las = 1)

abline(h = -1, col = "RED", lty= 2)
abline(h = 1, col = "RED", lty = 2)


## GO over-representation
require(clusterProfiler)
id11 <- getID(setdiff(names(which(diff.test$judge.comb[, 2])), diff.gene.list.inten))
id12 <- getID(intersect(names(which(diff.test$judge.comb[, 2])), diff.gene.list.inten))
yy1 <- enrichGO(id11, "org.Mm.eg.db", ont = "BP", readable = T,
               pAdjustMethod = "BH")
yy2 <- enrichGO(id12, "org.Mm.eg.db", ont = "BP", readable = T,
               pAdjustMethod = "BH")
p1 <- dotplot(yy1, showCategory = 15)
p2 <- dotplot(yy2, showCategory = 15)
library(grid)
library(gridExtra)
grid.arrange(p1, p2, ncol = 1)



################# Compare cell size among the cell types ####################################
## need to specify DATAFILEPATH, trueMol and alpha for the code to run, see the file functions_all.R
data <- GetY(DATAFILEPATH, trueMol, alpha)
alpha <- data$alpha
log.lib.size <- log(colSums(data$Y))

log.cell.size <- log.lib.size - log(alpha)
cell.size.avg <- sapply(1:7, function(i) 
                        mean(exp(log.cell.size[data$Ctype == levels(data$Ctype)[i]])))
names(cell.size.avg) <- name.type[1:7]
## bar-plot
size.data <- data.frame(value = cell.size.avg, cell.types = levels(data$Ctype), stringsAsFactors =F)
size.data$Type <- as.factor(size.data$cell.types %in% levels(data$Ctype)[c(1, 2, 4:5)])
levels(size.data$Type) <- c("Neurons", "Non-neurons")
p3 <- ggplot(size.data, aes(x = cell.types, y = value, fill = Type)) + 
     geom_bar(stat = "identity") + 
     theme_minimal() +  scale_fill_manual(values = c("#E69F00", "#56B4E9")) + 
     theme(axis.text.x = element_text(size = 13, angle = 60, 
                                    hjust = 1, vjust = 1), 
           axis.text.y = element_text(size = 14), 
           plot.title = element_text(size = 16, lineheight = 1),
           axis.title = element_blank(),
           legend.text = element_text(size = 13),
           legend.title = element_blank(), legend.position = "bottom",
plot.margin = unit(x = c(0.2, 0.1, 0.1, 0), units = "cm"))
p3

## violin plot
cell.size.list <- sapply(1:7, function(i) 
                        (log.cell.size)[data$Ctype == levels(data$Ctype)[i]])
names(cell.size.list) <- name.type[1:7]
size.data <- data.frame(value = unlist(cell.size.list), 
                        cell.types = unlist(lapply(levels(data$Ctype), 
                                                   function(str) rep(str, sum(data$Ctype == str)))), 
                        stringsAsFactors =F)
size.data$Type <- as.factor(size.data$cell.types %in% levels(data$Ctype)[c(1, 2, 4:5)])
levels(size.data$Type) <- c("Neurons", "Non-neurons")
p4 <- ggplot(size.data, aes(x = cell.types, y = value, fill = Type)) + 
#     geom_bar(stat = "identity") + 
   geom_violin(position = position_dodge(0.5), trim = T,
                 draw_quantiles = 0.5, size = 0.3) +
     theme_minimal() +  scale_fill_manual(values = c("#E69F00", "#56B4E9")) + 
     theme(axis.text.x = element_text(size = 13, angle = 60, 
                                    hjust = 1, vjust = 1), 
           axis.text.y = element_text(size = 14), 
           plot.title = element_text(size = 16, lineheight = 1),
           axis.title = element_blank(),
           legend.text = element_text(size = 13),
           legend.title = element_blank(), legend.position = "bottom",
plot.margin = unit(x = c(0.2, 0.1, 0.1, 0), units = "cm"))
p4

