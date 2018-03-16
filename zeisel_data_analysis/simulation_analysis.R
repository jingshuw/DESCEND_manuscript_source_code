library(DESCEND)
source("functions_all.R")


####################### Sample-splitting Simulation ##########################

result <- readRDS("splitting_descend_results.rds")

ests <- getEstimates(result$result1)
est1 <- sapply(ests[c(1, 2, 4, 5)], function(mm)mm[, 1])
sd1 <- sapply(ests[c(1, 2, 4, 5)], function(mm)mm[, 2])

ests <- getEstimates(result$result2)
est2 <- sapply(ests[c(1, 2, 4, 5)], function(mm)mm[, 1])
sd2 <- sapply(ests[c(1, 2, 4, 5)], function(mm)mm[, 2])


gc()


bias1 <- sapply(result1, function(l) if(is.na(l)) rep(NA, 4) else l$estimates[c(1, 4:6), 2])
bias2 <- sapply(result2, function(l) if(is.na(l)) rep(NA, 4) else l$estimates[c(1, 4:6), 2])


score <- (est1 - est2) / sqrt(sd1^2 + sd2^2)

pval <- 2 * pt(-abs(score), df= 50)

est10 <- sapply(result10, function(l) if(is.na(l)) rep(NA, 6) else l$estimates[, 1])
est20 <- sapply(result20, function(l) if(is.na(l)) rep(NA, 6) else l$estimates[, 1])


idx <- !is.na(est1[, 2]) & !is.na(est2[, 2])

########### Figure for sample-splitting ##################
layout(matrix(1:3, nrow = 1, byrow = T))
par(mar = c(2, 2, 2, 1), mgp = c(0.5, 0.5, 0))

smoothScatter(log10(est1[idx, 2]), log10(est2[idx, 2]), xlab= "", ylab = "",
              main = "", pch = "",
              xaxt = "n", yaxt ="n")
axis(1, at = c(1, 2, 3), tck = -0.02, cex.axis = 1.4)
axis(2, at = c(1, 2, 3), tck = -0.02, cex.axis = 1.4, las = 1)
title("Nonzero Mean(log10)", line = 0.5, font.main = 1, cex.main = 1.4)
abline(a = 0, b= 1)


smoothScatter(est1[idx, 1], est2[idx, 1], xlab= "", ylab = "",
              cex.axis = 1.1,
              main = paste(""), pch  = "", 
              font.main  = 1, cex.main = 1.3, xaxt = "n", yaxt = "n")
axis(1, at = c(0.2, 0.6, 1.0), tck = -0.02, cex.axis = 1.4)
axis(2, at = c(0.2, 0.6, 1.0), tck = -0.02, cex.axis = 1.4, las = 1)
title("Nonzero Fraction", line = 0.5, font.main = 1, cex.main = 1.5)
abline(a = 0, b= 1)

smoothScatter(est1[idx, 4], est2[idx, 4], xlab= "", ylab = "",
              main = "", pch = "", xaxt = "n", yaxt = "n")
axis(1, at = c(0.2, 0.5, 0.8), tck = -0.02, cex.axis = 1.4)
axis(2, at = c(0.2, 0.5, 0.8), tck = -0.02, cex.axis = 1.4, las = 1)
title("Gini", line = 0.5, font.main = 1, cex.main = 1.5)
abline(a = 0, b= 1)





############################# Parametric Simulation ############################################################ Parametric Simulation ############################################################ Parametric Simulation ###############################

descend.result <- readRDS("parametric_simu_descend_results.rds")
ests <- getEstimates(descend.result)
est.par <- sapply(ests, function(mm)mm[, 1])
pval <- getPval(descend.result)

simu.data <- readRDS("parametric_simulated_data.rds")

p0.true <- 1 - simu.data$true.nonzero.frac
## beta is the negative of true
coefs.true <- simu.data$cell.size.effect.sizes
stats.true <- simu.data$stats



###### Figure for the parametric simulation ########
layout(matrix(c(rep(1, 7), rep(2, 7), rep(3, 8)), nrow = 1, byrow = T))
par(mar = c(2, 2, 3, 1), mgp = c(0.5, 0.5, 0))


idx1 <- coefs.true[, 1] > 0 & coefs.true[, 2] < 2
smoothScatter(coefs.true[idx1, 1], est.par[idx1, 6], xlim = c(0, 2),
              ylim = c(0, 2), pch = "", xaxt = "n", yaxt = "n",
              xlab = "", ylab = "")
axis(1, at = c(0, 1, 2), tck = -0.02, cex.axis = 1.4)
axis(2, at = c(0, 1, 2), tck = -0.02, cex.axis = 1.4, las = 1)
title(" coefficient on \n Nonzero Mean", 
      line = 0.5, font.main = 1, cex.main = 1.5)
abline(a = 0, b= 1)


idx2 <- p0.true != 0
smoothScatter(coefs.true[idx2, 2], -est.par[idx2, 8], 
              xlab= "", ylab = "",
              main = "", xlim = c(-1, 3), ylim = c(-1, 3), pch ="",
              xaxt = "n", yaxt = "n")

axis(1, at = c(-1, 1, 3), tck = -0.02, cex.axis = 1.4)
axis(2, at = c(-1, 1, 3), tck = -0.02, cex.axis = 1.4, las = 1)
title(" coefficient on \n Nonzero Fraction", 
      line = 0.5, font.main = 1, cex.main = 1.5)
abline(a = 0, b= 1)


par(mar = c(2, 4, 3, 2))
Fdp <- function(alpha.vec, pval, p0.true) {
  temp <- p.adjust(pval[, 1], "BH")
  aa <- sapply(alpha.vec, function(alpha)sum(temp <= alpha & p0.true == 0, na.rm = T)/
               sum(temp<=alpha, na.rm = T))
  plot(alpha.vec, aa, pch = 20, 
       xlab = "", ylab = "",
       xlim = c(0, 0.5), ylim = c(0, 0.5), 
       xaxt = "n", yaxt = "n",
       main = "")
  abline(a = 0, b= 1, col = "RED")
  axis(1, at = c(0, 0.2, 0.5), tck = -0.02, cex.axis = 1.4)
  axis(2, at = c(0, 0.2, 0.5), tck = -0.02, cex.axis = 1.4, las = 1)
  title("LRT Test on \n Active Fraction = 1", 
        line = 0.5, font.main = 1, cex.main = 1.5)

  return(aa)
}
Fdp(1:25 * 0.02, pval, p0.true)



###### Supplementary figure for the parametric simulation ########
layout(matrix(1:3, nrow = 1, byrow = T))
par(mar = c(2, 2, 3, 2))
smoothScatter(1-p0.true[p0.true!=0], est.par[p0.true!=0, 1], 
              xlab= "", ylab = "", xaxt = "n", yaxt = "n",
              cex.axis= 1.1, font.main = 1, cex.main = 1.3, 
              main = "", pch = "")
axis(1, at = c(0.2, 0.5, 0.8), tck = -0.02, cex.axis = 1.4)
axis(2, at = c(0.2, 0.6, 1.0), tck = -0.02, cex.axis = 1.4, las = 1)
title("Nonzero Fraction", line = 0.5, font.main = 1, cex.main = 1.5)
abline(a = 0, b= 1)

smoothScatter(stats.true[, 2], est.par[, 4], xlab= "", ylab = "", 
              xlim = c(0, 3), ylim = c(0, 3), xaxt = "n", yaxt = "n",
              cex.axis= 1.1, font.main = 1, cex.main = 1.3, 
              main = "", pch = "")
axis(1, at = c(0, 1.5, 3), tck = -0.02, cex.axis = 1.4)
axis(2, at = c(0, 1.5, 3), tck = -0.02, cex.axis = 1.4, las = 1)
title("CV", line = 0.5, font.main = 1, cex.main = 1.5)
abline(a = 0, b= 1)

smoothScatter(stats.true[, 3], est.par[, 5], xlab= "", ylab = "",
              cex.axis= 1.1, font.main = 1, cex.main = 1.3, 
              xaxt = "n", yaxt = "n",
              main = "", pch = "")
axis(1, at = c(0.2, 0.5, 0.8), tck = -0.02, cex.axis = 1.4)
axis(2, at = c(0.2, 0.5, 0.8), tck = -0.02, cex.axis = 1.4, las = 1)
title("Gini", line = 0.5, font.main = 1, cex.main = 1.5)
abline(a = 0, b= 1)



###################### Figure for down-sampling ###############################################
result <- readRDS("downsampling_simulation.rds")
library(reshape2)
library(ggplot2)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols <- gg_color_hue(3)

data.rr <- as.data.frame(result$Gini)
p1 <- PlotDownSample(data.rr, "Gini", ylim = c(0.2, 1), legend.position = "right",
                     col = cols)
data.rr <- as.data.frame(result$CV)
p2 <- PlotDownSample(data.rr, "CV", logY = T, ylim = c(0.5, 3.5),
                     col = cols)
data.rr <- as.data.frame(result$nonzero.frac)
p3 <- PlotDownSample(data.rr, "Nonzero Fraction", ylim = c(0, 1), col = cols)

library(grid)
library(gridExtra)
## boxplot for the distribution of estimated measures across genes 
grid.arrange(grobs = list(p3, p2, p1), widths = c(1,1, 1.6), nrow = 1)


data.rr <- as.data.frame(result$Gini)
p1 <- PlotDownSample(data.rr, "Gini", legend.position = "right", 
                     method = "Diff", col = cols)
data.rr <- as.data.frame(result$CV)
p2 <- PlotDownSample(data.rr, "CV", method = "Diff", col = cols)
data.rr <- as.data.frame(result$nonzero.frac)
p3 <- PlotDownSample(data.rr, "Nonzero Fraction", 
                     method = "Diff", col = cols)
## boxplot for the distribution of the difference between estimated measures and true values (supplementary figure)
grid.arrange(grobs = list(p3, p2, p1), widths = c(1,1, 1.5), nrow = 1)



