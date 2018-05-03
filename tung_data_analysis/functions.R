plotGini <- function(ests, pvals, xlab= "", ylab = "") {
  temp <- sapply(ests, function(ll)ll[[5]][, 1])
  plot(sqrt(temp), 
       pch = 20, xlim = sqrt(c(0.005, 0.3)), ylim = sqrt(c(0.005, 0.3)), 
       xaxt = "n", yaxt = "n", col = rgb(0, 0, 0, 0.3), 
       xlab = xlab, ylab= ylab,
       cex.lab = 1.3)
  at.pt <- c(0.01, 0.1, 0.2)
  axis(1, at = sqrt(at.pt), labels = at.pt, cex.axis = 1.4, tck = -0.02)
  axis(2, at = sqrt(at.pt), labels = at.pt, cex.axis = 1.4, tck = -0.02, las = 1)
  abline(a = 0, b= 1, col = "blue", lwd = 1.2)
  points(sqrt(temp[p.adjust(pvals[, 5], "BH") < 0.05, ]), 
         pch = 20, col = rgb(1, 0, 0, 0.5))

}

plotCV <- function(ests, pvals, lims,
                   xlab = "", ylab = "") {
  temp <- sapply(ests, function(ll)ll[[4]][, 1])
  plot(sqrt(temp), 
       pch = 20, xlim = sqrt(lims), ylim = sqrt(lims), 
       xaxt = "n", yaxt = "n", col = rgb(0, 0, 0, 0.3), 
       xlab = xlab, ylab= ylab,
       cex.lab = 1.3)
  at.pt <- c(0.1, 0.2, 0.5)
  axis(1, at = sqrt(at.pt), labels = at.pt, cex.axis = 1.4, tck = -0.02)
  axis(2, at = sqrt(at.pt), labels = at.pt, cex.axis = 1.4, tck = -0.02, las = 1)
  abline(a = 0, b= 1, col = "blue", lwd = 1.2)
  points(sqrt(temp[p.adjust(pvals[, 4], "BH") < 0.05, ]), 
         pch = 20, col = rgb(1, 0, 0, 0.5))
}
