## get 3005 cells expression data
## original data of Zeisel et. al. can be downloaded from GSE60361
## computation of trueMol and alpha use function getData in the file functions_ercc.R in the folder ERCC_validation and the data in that folder
GetY <- function(data.file = "GSE60361_C1-3005-Expression.txt",
                 trueMol, alpha) {

  require(data.table)
  genestab<- fread(data.file, sep="\t", header=TRUE)
  headers<- read.table(paste(base, "zeisel_header_only.txt", 
                             sep = ""), sep="\t", header=TRUE)
  level1 <- paste(t(headers[8,-1]))
  gnames <- genestab[,1]
  Y <- as.matrix(genestab[,-1])
  rownames(Y) <- gnames

  ## remove Malat1 
  Y <- Y[-which(gnames == "Malat1"), ]
  
  return(list(Y = Y[rowSums(Y!=0) > 100, ], 
              trueMol = trueMol,
              Ctype = as.factor(level1),
              alpha = alpha))
}

getCommonGenes <- function(result.list) {
  gene.list <- lapply(result.list, 
                      function(ll) {
                        temp <- unique(rownames(ll$ests[!is.na(ll$ests[, 3]), ]))
                        temp[!is.na(temp)]
                      })
  common.genes <- Reduce(intersect, gene.list)
  return(common.genes)
}


## to plot the coefficient of cell size on nonzero Mean
sizeEffect <- function(i, selected.genes = NULL, fish.values = NULL, add.density = F) { 
  ylim  <- c(-1, 3)
  k <- 7
  col1 <- rgb(0, 0, 0, alpha = 0.5)
  col2 <- rgb(1, 0, 0, alpha = 0.5)

  if (is.null(selected.genes))
    selected.genes <- 1:nrow(all.results[[i]]$ests)

  ## 4th column is Mean
  idx <- !is.na(all.results[[i]]$ests[selected.genes, 3])
  ## genes whose gamma is significantly not 1
  sig.genes <- p.adjust(all.results[[i]]$pvals[selected.genes, 2], "BH") < 0.05
  xx <- sort(log(all.results[[i]]$ests[selected.genes, 3][idx]), index.return = T)
  
  yy <- all.results[[i]]$ests[selected.genes, 7][idx][xx$ix] 

  if (add.density)
    layout(matrix(c(1, 1, 2), nrow = 1))

  ### idx and idx2 are used to plot only a sample of points to make the figure smaller
  idx <- sample(1:length(yy), 1000)

  plot(xx$x[idx], yy[idx], ylim = ylim,
       pch = 20, main = name.type[i], font.main = 1, 
       col = col1, cex.main = 1.4,
       xlab = "", 
       ylab = "", xaxt = "n", yaxt = "n")
  axis(1, at = c(-1, 1, 3), tck = -0.02, cex.axis = 1.4)
  axis(2, at = c(-1, 1, 3), tck = -0.02, cex.axis = 1.4, las = 1)
  m <- sum(sig.genes, na.rm = T)
  idx2 <- sample.int(m, min(m, 500))
  points(log(all.results[[i]]$ests[selected.genes, 3])[sig.genes][idx2],
         all.results[[i]]$ests[selected.genes, 7][sig.genes][idx2],
         pch = 20, col = col2)

  abline(h = 1, col = "blue", lty = 2, lwd = 1.5)
  temp <- loess(yy ~ xx$x, family = "symmetric")

  if (add.density) {
    idx <- yy > ylim[1] & yy < ylim[2] & 
    !is.na(yy)

    temp <- density(yy)
     temp1 <- density(fish.values)
    mars <- par("mar")
    mars[2] <- 1
    par(mar = mars)
    plot(temp$y, temp$x, ylim = ylim, xlab= "", ylab= "",
         main = "", axes = F, type = "l", 
         xlim = range(c(temp$y, temp1$y)), lwd = 1.3)
       points(temp1$y, temp1$x, type = "l", col = "cyan3", lwd = 1.3)
    axis(2, fish.values, labels = F, col.ticks = "cyan3")
    abline(h = 1, col = "blue", lty = 2, lwd = 1.5)

    legend("bottomright", lty= rep(1, 2), col = c("black", "cyan3"),
           lwd = rep(2, 1.3),
           legend = c("Endothelial-\nMural", "RNA FISH"), bty = "n")
  }
}

## to plot the coefficient of cell size on nonzero Fraction
freqEffect <- function(i, add.density = F, 
                       fish.values= NULL) {
  results <- all.results[[i]]

  p0.pval <- results$pvals[, 1]
  beta.pval <- results$pval[, 3]
  ## only plot genes whose estimated nonzero fraction is significant less than 1 and the estimate less than 0.9 
  gene.sel <- p.adjust(p0.pval, "BH") < 0.05 
  gene.sel <- gene.sel & results$ests[, 1] < 0.9

  ## significant genes (will be plotted with red dots)
  beta.gene.sel <- p.adjust(beta.pval[gene.sel], "BH") < 0.05

  col1 <- rgb(0, 0, 0, alpha = 0.5)
  col2 <- rgb(1, 0, 0, alpha = 0.5)

  if (add.density)
    layout(matrix(c(1, 1, 2), nrow = 1))

  idx <- sample(1:sum(gene.sel, na.rm = T), min(sum(gene.sel, na.rm = T), 1000))

  plot(results$ests[gene.sel, 1][idx], -results$ests[gene.sel, 10][idx],
       pch = 20, ylim = c(-5, 6), col = col1, 
       xlim = c(0, 1), 
       xlab = "", ylab = "", xaxt = "n", yaxt = "n",
       main = name.type[i], font.main = 1, cex.main = 1.4)
  axis(1, at = c(0, 0.5, 1), tck = -0.02, cex.axis = 1.4)
  axis(2, at = c(-4, 0, 6), tck = -0.02, cex.axis = 1.4, las = 1)

  m <- sum(beta.gene.sel, na.rm  =T)
  idx2 <- c(sample.int(m, min(m, 400)), 
            intersect(idx, which(beta.gene.sel)))

  points(results$ests[gene.sel, 1][beta.gene.sel][idx2],
         -results$ests[gene.sel, 10][beta.gene.sel][idx2],
         pch = 20, col = col2)
  abline(h = 0, lty = 2, col ="blue")
  
  if (add.density) {
    idx <- results$ests[gene.sel, 10] < 5 & 
    results$ests[gene.sel, 10] > -6 & 
    !is.na(results$ests[gene.sel, 10])

    temp <- density(-results$ests[gene.sel, 10][idx])
    temp1 <- density(fish.values)
    mars <- par("mar")
    mars[2] <- 1
    par(mar = mars)
    plot(temp$y, temp$x, ylim = c(-5, 6), xlab= "", ylab= "",
         main = "", axes = F, type = "l", 
         xlim = range(c(temp$y, temp1$y)))
    points(temp1$y, temp1$x, type = "l", col = "cyan3", lwd = 1.3)
    axis(2, fish.values, labels = F, col.ticks = "cyan3", lwd = 1.3)
    abline(h = 0, lty = 2, col ="blue")

    legend("bottomright", lty= rep(1, 2), col = c("black", "cyan3"),
           lwd= rep(2, 1.3),
           legend = c("Endothelial-\nMural", "RNA FISH"), bty = "n")
  } 
}

### the two functions below are for GO analysis of the results
getGeneList <- function(DE.result, alpha = 0.05) {
  list.before <- p.adjust(DE.result$p.values[1, ], "BH") < alpha
  list.after <- p.adjust(DE.result$p.values[6, ], "BH") < alpha
  list.mean <- p.adjust(DE.result$p.values[2, ], "BH") < alpha


  require(gProfileR)
  require(Hmisc)
  enriched.sets <- gprofiler(names(which(list.after)), organism = "mmusculus", 
                             src_filter = "GO")
  enriched.sets$intersection <- sapply(enriched.sets$intersection, 
                                       function(str)paste(capitalize(tolower(strsplit(str, ",")[[1]])), 
                                                          sep = ","))
  judge.comb <- cbind(list.before, list.after)
  colnames(judge.comb) <- c("Before", "After")
  return(list(judge.comb = judge.comb, 
              judge.mean = list.mean,
              enriched.aftersets = enriched.sets))
  
}

getID <- function(gene.vec) {
  library(org.Mm.eg.db)
  xx <- as.list(org.Mm.egALIAS2EG)
  xx <- xx[!is.na(xx)]
  ids <- xx[gene.vec]
  ids <- lapply(ids, function(ll)ll[1])
  gene.list <- names(ids)
  gene.list <- gene.list[!is.na(gene.list)]
  ids <- as.numeric(unlist(ids))
  names(ids) <- gene.list
  return(ids)
}


### Figure for showing the down-sampling results ##############
PlotDownSample <- function(data, measure.name, col, logY = F, ylim = NULL,
                           legend.position = "none",
                           method = c("All", "Diff")) {
  method <- match.arg(method, c("All", "Diff"))
  if (method == "All")
    data[, -1] <- data[, -1] + data[, 1]
  temp.data <- reshape(data, varying = colnames(data), direction = "long", v.name = "value",
                       times = colnames(data), timevar = "Group")
  temp.data$Method <- rep("True", nrow(temp.data))
  temp.data$Method[grep("descend", temp.data$Group)] <- "DESCEND"
  temp.data$Method[grep("raw", temp.data$Group)] <- "Raw counts"
  temp.data$eff <- rep("True", nrow(temp.data))
  temp.data$eff[grep("20", temp.data$Group)] <- "20%"
  temp.data$eff[grep("10", temp.data$Group)] <- "10%"
  temp.data$eff[grep("5", temp.data$Group)] <- "5%"
  temp.data$eff <- factor(temp.data$eff, levels = c("True", "20%", "10%", "5%"), ordered = T)

  if (method == "Diff")
    temp.data <- temp.data[temp.data$Method != "True", ]
  p <- ggplot(temp.data, aes(x = eff, y = value, fill = Method)) + 
  geom_boxplot(position = position_dodge(0.7), width = 0.5, outlier.size = 0.2, outlier.shape = 15, lwd= 0.2) +
  xlab("") + ylab("") + ggtitle(measure.name) + 
  scale_fill_manual(values=col) + 
 # theme_minimal() + 
  theme(axis.text = element_text(size = 16), 
        title = element_text(size = 16),
        legend.position = legend.position,
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  guides(fill=guide_legend(
                         keywidth=0.3,
                         keyheight=0.3,
                         default.unit="inch")
      )
  if (method == "Diff")
    p <- p + geom_hline(yintercept = 0, linetype = 2)
  if (logY)
    p <- p + scale_y_continuous(trans = "log2", limits = ylim)
  if (!is.null(ylim) & !logY)
    p <- p + ylim(ylim) 

  return(p)
}


