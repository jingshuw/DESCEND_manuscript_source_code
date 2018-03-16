library(descend)
result.list <- lapply(c(0, 2, 4, 7), function(days) {
  result <- readRDS(paste("descend_result_day", days, ".rds", sep = ""))
  est <-  getEstimates(result)
  return(est[c("Mean", "Gini")])
})

gini.alldays <- sapply(result.list, 
                       function(ll) ll$Gini[, 1])
dropseq.raw <- readRDS("dropseq_raw_metrics.rds")
gini.dropseq <- sapply(1:4, function(i) {
                       dropseq.raw[[i]]$Gini})

library(colorspace)
col.seq <- 
sequential_hcl(15, power = 2.2, c= c(100, 30), l = c(40, 90), h = 130)[4:1 * 2 - 1]

thres <- 0.05
idx <- dropseq.raw[[1]]$pos.frac > thres & 
  dropseq.raw[[2]]$pos.frac > thres & dropseq.raw[[3]]$pos.frac > thres & dropseq.raw[[4]]$pos.frac > thres



library(ggplot2)
source("functions.R")

###################### Plot for the Gini of raw counts ##############################
gini.temp <- gini.dropseq
gini.temp[!idx, ] <- NA
data <- data.frame(gini.temp)
data <- reshape(data, direction = "long", varying = paste("X", 1:4, sep = ""),
                  v.names = "Gini")
data$time <- as.factor(data$time)
levels(data$time) <- c(0, 2, 4, 7)
p <- ggplot(data, aes(x = time, y = Gini, fill = time)) +
  geom_violin(trim=T, draw_quantiles = 0.5, size = 0.3)+
  theme_minimal() + 
  scale_y_continuous(breaks = c(0.3, 0.7, 1)) + 
  scale_fill_manual(values = col.seq) + 
  ggtitle("Gini coefficients \n of raw normalized counts") + 
  theme(axis.text=element_text(size=14),
        legend.position = "none",
        plot.title = element_text(size = 14, lineheight = 1),
        axis.title = element_blank())
p

###################### Plot for the Gini from DESCEND  ##############################
data <- data.frame(gini.alldays)
gini.temp <- gini.alldays
gini.temp[!idx, ] <- NA
data <- data.frame(gini.temp)
data <- reshape(data, direction = "long", varying = paste("X", 1:4, sep = ""),
                  v.names = "Gini")
data$time <- as.factor(data$time)
levels(data$time) <- c(0, 2, 4, 7)
p <- ggplot(data, aes(x = time, y = Gini, fill = time)) +
  geom_violin(trim=T, draw_quantiles = 0.5, size = 0.3)+
  theme_minimal() + 
  scale_y_continuous(breaks = c(0.2, 0.5, 0.8)) + 
  scale_fill_manual(values = col.seq) + 
  ggtitle("Gini coefficients \n estimated by DESCEND") + 
  theme(axis.text=element_text(size=14),
        legend.position = "none",
        plot.title = element_text(size = 14, lineheight = 1),
        axis.title = element_blank())
p

################### Library size barplot ###############################
par(mar = c(3, 3, 2, 1))
lib.sizes <- sapply(dropseq.raw, function(ll) mean(ll$lib.size))
names(lib.sizes) <- c(0, 2, 4, 7)
barplot(lib.sizes, main = "Average Library Size", cex.main = 1.4, font.main = 1,
        cex.names = 1.4, font.axis = 1, cex.axis = 1.4, yaxt = "n",
        col =col.seq)
axis(2, tck = -0.02, cex.axis = 1.3)

## compare Gini form raw counts at Day 2 VS Day 4
col <-rgb(0,0,0,alpha=0.3)
par(mar = c(4, 4, 1, 1))
plot(gini.dropseq[, 2:3], pch = 20, xlab = "Day 2", ylab = "Day 4",
     cex.lab = 1.4, font.lab = 1, col = col, cex.axis= 1.2)
abline(a = 0, b=  1, col = "RED", lwd = 2)


##############  Gini correlation plot #################################
cor.mat <- round(cor(gini.temp, use = "complete.obs"), 2)
par(mar = c(4, 4, 3, 2), mgp = c(0, 0.8, 0))
image(cor.mat, zlim = c(0.3, 1), 
      col = gray((32:10)/32), xaxt = "n", yaxt = "n", xlab = "", 
      ylab = "",
      main = "Correlation of \n DESCEND estimated Gini", 
      font.main = 1, cex.main = 1.3)
pos <- seq(0, 1, length.out = 4)
axis(1, pos, paste("", c(0, 2, 4, 7),sep= ""), cex.axis = 1.4, font.axis = 1, 
     tck = -0.02)
axis(2, pos, paste("", c(0, 2, 4, 7),sep= ""), cex.axis = 1.4, font.axis = 1,
     tck = -0.02, las = 1)
mat.pos <- matrix(rep(pos, 4), nrow = 4)
cor.mat <- matrix(as.character(cor.mat), nrow = 4)
diag(cor.mat) <- ""
text(as.vector(mat.pos), as.vector(t(mat.pos)), as.vector(cor.mat), col ="RED", font = 1, cex = 1.3)


############ Figures of the change of Gini for representative genes ############

### differentiation markers ###
par(mar = c(2.5, 2.7, 2.7, 1), mgp = c(1.3, 0.7, 0))
layout(matrix(c(1, 1, 2, 2, 3), nrow = 1, byrow = T))
gene.vec <- c("Krt8", "Krt18", "Tagln", "Cald1", "Tpm1", "Fxyd6")
cols <- rainbow_hcl(length(gene.vec),  c = 90, l = 65)

setPlot(gene.vec, result.list, 1, log = "y", add.legend = F, cols= cols)
setPlot(gene.vec, result.list, 3, log = "", add.legend= F, cols = cols)
par(mar = c(0, 0,0,0))
plot(0:1, type = "n", xaxt = "n", yaxt = "n", bty = "n")
legend("right", col = cols, lty = rep(1, length(gene.vec)),
         lwd = rep(2, length(gene.vec)),
         legend = gene.vec, cex = 1.2, text.font = 1, bty = "n")


### pluripotency factors ###
par(mar = c(3, 2.7, 2.7, 1), mgp = c(1.6, 0.7, 0))
layout(matrix(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4), nrow = 1, byrow = T))
gene.vec2 <- c("Pou5f1", "Dppa5a", "Sox2", "Nanog", "Zfp42", "Klf4", "Ccnd3")

cols <- rainbow_hcl(length(gene.vec2),  c = 90, l = 65)
setPlotSimple(gene.vec2, result.list, 1, log = "y", add.legend = F, cols = cols)
setPlotSimple(gene.vec2, result.list, 3, log = "", add.legend= F, cols = cols)
oriSetPlotSimple(gene.vec2, gini.dropseq, add.legend = F, cols = cols)
par(mar = c(0, 0,0,0))
plot(0:1, type = "n", xaxt = "n", yaxt = "n", bty = "n")
gene.vec2[1] <- "Pout5f1(Oct4)"
gene.vec2[5] <- "Zfp42(Rex1)"
legend("right", col = c(cols[1:3], cols[1], cols[4:6], cols[1], cols[7]), 
       lty = c(1, 1, 1, 0, 1,  1, 1, 0, 1),
       lwd = rep(2, length(gene.vec2)),
       legend = c(gene.vec2[1:3], 
                  "", gene.vec2[4:6], "", gene.vec2[7]), cex = 1, text.font = 1, bty = "n")

############## Supplementary figures of the comparison of mean and Gini for other genes #############
par(mar = c(2.5, 2.7, 2.7, 1), mgp = c(1.3, 0.7, 0))
layout(matrix(c(1, 1, 2, 2, 3), nrow = 1, byrow = T))

gene.vec <- c("Tagln", "Anxa2",  "H19",  "Sparc","Ccno")
cols <- rainbow_hcl(length(gene.vec),  c = 90, l = 65)

setPlot(gene.vec, result.list, 1, log = "y", add.legend = F, cols= cols)
setPlot(gene.vec, result.list, 3, log = "", add.legend= F, cols = cols)
par(mar = c(0, 0,0,0))
plot(0:1, type = "n", xaxt = "n", yaxt = "n", bty = "n")
legend("right", col = cols, lty = rep(1, length(gene.vec)),
         lwd = rep(2, length(gene.vec)),
         legend = gene.vec, cex = 1.2, text.font = 1, bty = "n")


par(mar = c(2.5, 2.7, 2.7, 1), mgp = c(1.3, 0.7, 0))
layout(matrix(c(1, 1, 2, 2, 3), nrow = 1, byrow = T))
gene.vec <- c("Jun", "Anxa3","Klf6", "Fos", "Dusp4")
cols <- rainbow_hcl(length(gene.vec),  c = 90, l = 65)

setPlot(gene.vec, result.list, 1, log = "y", add.legend = F, cols= cols)
setPlot(gene.vec, result.list, 3, log = "", add.legend= F, cols = cols)
par(mar = c(0, 0,0,0))
plot(0:1, type = "n", xaxt = "n", yaxt = "n", bty = "n")
legend("right", col = cols, lty = rep(1, length(gene.vec)),
         lwd = rep(2, length(gene.vec)),
         legend = gene.vec, cex = 1.2, text.font = 1, bty = "n")

################ Venn Diagrams in the supplementary figure#####################################
load("../deseq_result_summary.rda")
test.result <- readRDS("descend_test_result_Day02.rds")
v1 <- res$padj < 0.05
load("../deseq_result_summary.rda")
v0 <- res$padj < 0.05
pvals <- test.result$p.values[, "Gini"]
v2 <- p.adjust(pvals, "BH") < 0.05
v3 <- p.adjust(test.result$p.values[, "Mean"], "BH") < 0.05
v1[is.na(v2)] <- NA
require(VennDiagram)
venn.list <- list(which(v1), which(v2), which(v3))
names(venn.list) <- c("Mean (DESeq2)", "Gini", "Mean (DESCEND)")
## compare Mean (DESeq2) with Gini
venn.plot <- venn.diagram(venn.list[1:2],
                          NULL, 
                          height = 3000,
                          width = 3000,
                          fill=c("blue", "red"), 
                          alpha=c(0.3,0.3), 
                          main.cex = 4,cex = 4, lwd = 1,
                          cat.pos = 1.7, cat.cex = 4,
                          cat.just = list(c(0.4, 1), c(0.6, 0.5)))
grid.newpage()
grid.draw(venn.plot)

## compare Mean(DESCEND) with Gini
venn.plot <- venn.diagram(venn.list[c(3, 2)],
                          NULL, 
                          height = 3000,
                          width = 3000,
                          fill=c("green", "red"), 
                          alpha=c(0.3,0.3), lwd = 1, 
                          main.cex = 4,cex = 4, 
                          cat.pos = 1.7, cat.cex = 4,
                          cat.just = list(c(0.4, 1), c(0.6, 0.5)))
grid.newpage()
grid.draw(venn.plot)

## compare Mean(DESeq2) with Mean(DESCEND)
venn.plot <- venn.diagram(venn.list[c(1,3)],
                          NULL, 
                          height = 3000,
                          width = 3000,
                          fill=c("blue", "green"), 
                          alpha=c(0.3,0.3), 
                          main.cex = 3,cex = 3,
                          cat.pos = 1.7, cat.cex = 3,
                          cat.just = list(c(1.1, 1), c(0.15, 0.5)))
grid.newpage()
grid.draw(venn.plot)


