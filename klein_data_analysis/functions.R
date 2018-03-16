setPlot <- function(gene.vec, ests.list, k,
                    add.legend = T, log = "y",
                    cols = NULL) {
  vals <- array(0, c(length(ests.list), length(gene.vec), 4))
  for (i in 1:length(ests.list)) {
    temp <- ests.list[[i]]
    vals[i, , ] <- cbind(temp$Mean[gene.vec, ], 
                         temp$Gini[gene.vec, ])
  }

  library(Hmisc)
  plot.names <- c("Mean", "",  
                  "Gini (DESCEND)")
  ylim <- range(vals[,,k] - vals[,,k + 1], vals[,,k] + vals[,,k+1],
                na.rm = T)
  require(colorspace)
  if (is.null(cols))
    cols <- rainbow_hcl(length(gene.vec),  c = 60, l = 75)

  x.temp <- c(0.4, 2:length(ests.list)) + rnorm(length(ests.list), sd = 0.07)
  plot(x.temp, vals[, 1, k], col = alpha(cols[1], 0.4), type = "l", 
       pch = 20, ylim = ylim, lwd = 2, yaxt = "n", xaxt = "n", xlim = c(0.2, 4.2),
       xlab = "", ylab = "", main = plot.names[k],
       cex.main = 1.5, font.main = 1, cex.axis = 1.4, font.axis = 1, log = log)
  points(x.temp[1:2], vals[1:2, 1, k], col = cols[1], type = "l", lwd = 2, 
       pch = 20)

  axis(1, c(0.4, 2:4), paste("", c(0, 2, 4, 7)), cex.axis = 1.4, tck = -0.02,  
       font.axis = 1)
  axis(2, tck = -0.02, cex.axis = 1.4)
  errbar(x.temp, vals[, 1, k], 
         vals[, 1, k] + vals[, 1, k + 1],
         vals[, 1, k] - vals[, 1, k + 1],
         errbar.col = alpha(cols[1], 0.3), add = T, pch = "", lwd = 2)
  errbar(x.temp[1:2], vals[1:2, 1, k], 
         vals[1:2, 1, k] + vals[1:2, 1, k + 1],
         vals[1:2, 1, k] - vals[1:2, 1, k + 1],
         errbar.col = cols[1], add = T, pch = "", lwd = 2)

  for (j in 2:length(gene.vec)) {
    x.temp <- c(0.4, 2:length(ests.list)) + rnorm(length(ests.list), sd = 0.07)
    points(x.temp, vals[, j, k], col = alpha(cols[j], 0.4), type = "l", 
           pch = 20, lwd = 2)
    points(x.temp[1:2], vals[1:2, j, k], col = cols[j], type = "l", 
           pch = 20, lwd = 2)

    errbar(x.temp, vals[, j, k], 
           vals[, j, k] + vals[, j, k + 1],
           vals[, j, k] - vals[, j, k + 1],
           errbar.col = alpha(cols[j], 0.3), add = T, pch = "", lwd = 2)
 errbar(x.temp[1:2], vals[1:2, j, k], 
           vals[1:2, j, k] + vals[1:2, j, k + 1],
           vals[1:2, j, k] - vals[1:2, j, k + 1],
           errbar.col = cols[j], add = T, pch = "", lwd = 2)
  }
  if (add.legend)
  legend("topleft", col = cols, lty = rep(1, length(gene.vec)),
         lwd = rep(2, length(gene.vec)),
         legend = gene.vec, cex = 1.1, text.font = 1)
  abline(v = 1.2, col = rgb(0, 0,0, 0.1), lty = 1, lwd = 90)

}


setPlotSimple <- function(gene.vec, ests.list, k,
                    add.legend = T, log = "y", cols = NULL) {
  vals <- array(0, c(length(ests.list), length(gene.vec), 4))
  for (i in 1:length(ests.list)) {
    temp <- ests.list[[i]]
    vals[i, , ] <- cbind(temp$Mean[gene.vec, ], 
                         temp$Gini[gene.vec, ])
  }

  library(Hmisc)
  plot.names <- c("Mean", "",  
                  "Gini (DESCEND)")
  ylim <- range(vals[,,k] - vals[,,k + 1], vals[,,k] + vals[,,k + 1],
                na.rm = T)
  require(colorspace)
  if (is.null(cols))
    cols <- rainbow_hcl(length(gene.vec),  c = 60, l = 75)

  x.temp <- 1:length(ests.list) + rnorm(length(ests.list), sd = 0.07)
  plot(x.temp, vals[, 1, k], col = cols[1], type = "l", 
       pch = 20, ylim = ylim, lwd = 2, xaxt = "n", yaxt = "n", xlim = c(0.8, 4.2),
       xlab = "", ylab = "", main = plot.names[k],
       cex.main = 1.5, font.main = 1, cex.axis = 1.3, font.axis = 1, log = log)
   axis(1, 1:4, paste("", c(0, 2, 4, 7)), cex.axis = 1.4, tck = -0.02,  
       font.axis = 1)
  axis(2, tck = -0.02, cex.axis = 1.4)
  errbar(x.temp, vals[, 1, k], 
         vals[, 1, k] + vals[, 1, k + 1],
         vals[, 1, k] - vals[, 1, k + 1],
         errbar.col = cols[1], add = T, pch = "", lwd = 2)

  for (j in 2:length(gene.vec)) {
    x.temp <- 1:length(ests.list) + rnorm(length(ests.list), sd = 0.07)
    points(x.temp, vals[, j, k], col = cols[j], type = "l", 
           pch = 20, lwd = 2)
      errbar(x.temp, vals[, j, k], 
           vals[, j, k] + vals[, j, k + 1],
           vals[, j, k] - vals[, j, k + 1],
           errbar.col = cols[j], add = T, pch = "", lwd = 2)

  }
  if (add.legend)
  legend("topleft", col = cols, lty = rep(1, length(gene.vec)),
         lwd = rep(2, length(gene.vec)),
         legend = gene.vec, cex = 1.1, text.font = 1)

}

oriSetPlotSimple <- function(gene.vec, gini.dropseq, 
                             add.legend = T, cols = NULL) {
  dropseq.gini <- gini.dropseq[gene.vec, ]

  ylim <- range(dropseq.gini,
                na.rm = T)
  require(colorspace)
  if (is.null(cols))
    cols <- rainbow_hcl(length(gene.vec),  c = 60, l = 75)
  x.temp <- 1:ncol(dropseq.gini) + rnorm(ncol(dropseq.gini), sd = 0.07)
  plot(x.temp, dropseq.gini[1, ], col = cols[1], type = "l", 
       pch = 20, ylim = ylim, lwd = 2, xaxt = "n", 
       yaxt = "n", xlim = c(0.8, 4.2),
       xlab = "", ylab = "", main = "Gini (InDrop)",
       cex.main = 1.5, font.main = 1, cex.axis = 1.3, font.axis = 1)
  axis(1, 1:4, paste("", c(0, 2, 4, 7)), cex.axis = 1.4, tck = -0.02,  
       font.axis = 1)
  axis(2, tck = -0.02, cex.axis = 1.4)

  for (j in 2:length(gene.vec)) {
    x.temp <- 1:ncol(dropseq.gini) + rnorm(ncol(dropseq.gini), sd = 0.07)
    points(dropseq.gini[j, ], col = cols[j], type = "l", 
           pch = 20, lwd = 2)

  }
  if (add.legend)
    legend("bottomright", col = cols, lty = rep(1, length(gene.vec)),
           lwd = rep(2, length(gene.vec)), bty = "n",
           legend = gene.vec, cex = 1, text.font = 1)

}
