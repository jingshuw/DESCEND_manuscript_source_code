source("functions_fish.R")
library(descend)

## The RNA FISH data is included in the folder
fish <- read.table("fishSubset.txt", header = T)


###################### load Dropseq data ##############################3 
### UMI counts can be downloaded from GSE99330. The Dropseq data for the FISH genes are included
### folder. One can ignore this section if they are only interested in the data of these genes
library(data.table)
## the name of the downloaded txt file of the Dropseq data
### need to change this to your own directary 
PATHFORTORREDROPSEQ <- "GSE99330_dropseqUPM.txt" 
header <- read.table(PATHFORTORREDROPSEQ, header = T, nrow = 1)
dropseq.counts <- fread(PATHFORTORREDROPSEQ, skip = 1, header = F)
dropseq.counts[, lapply(.SD, sum, na.rm=TRUE), .SDcols=-c("Gene")]
lib.size <- colSums(dropseq.counts)
setnames(dropseq.counts, c("Gene", colnames(header)))

dropseq.fish <- dropseq.counts[dropseq.counts$Gene %in% colnames(fish), ]

rm(dropseq.counts)
for (i in 1:5)
  gc()

dropseq.fish <- as.data.frame(dropseq)

####################################################################

## The Dropseq data for the RNA FISH genes is included in the folder
dropseq.fish <- read.table("dropseqFish.txt", header = T)
lib.size <- read.table("libSize.txt", header = T)
lib.size <- unlist(lib.size)
## FISH filtering (based on GAPDH) and normalization
fish.idx <- (fish$GAPDH > 100) & (fish$GAPDH < 1000)
fish <- fish[fish.idx, ]/fish$GAPDH[fish.idx]
## DropSeq filtering
data.idx <- lib.size > 1000 
dropseq.fish <- dropseq.fish[, data.idx]
lib.size <- lib.size[data.idx]



result <- runDescend(dropseq.fish, 
                     scaling.consts = lib.size,
                     control = list(max.sparse = c(0.98, 100)),
                     n.cores = 3, verbose = F)
## Remove two problematic genes
result[["FOSL1"]] <- NA
result[["VGF"]] <- NA

ests <- getEstimates(result)

pp <- plotCVGini()

legend.idx <- 25
p.list <- lapply(1:26, function(i) {
                 if (i == legend.idx)
                   plot.legend <- T
                 else
                   plot.legend <- F
                 plotShapeSingle(i, add.legend = plot.legend)})


library(grid)
library(gridExtra)
## For Gini
grid.arrange(pp[[1]], pp[[2]], grob(), ncol = 1)
## For CV
grid.arrange(pp[[3]], pp[[4]], pp[[7]], ncol = 1)
## For Nonzero Fraction
grid.arrange(pp[[5]], pp[[6]], pp[[8]], ncol = 1)


## Plot the distribution
idx <- sapply(p.list, function(ll)length(ll) > 0)
idx[4] <- F
temp.fun <- function(...)grid.arrange(..., ncol = 3)
do.call(temp.fun, p.list[idx])



### Figure cell size effect on active fraction
library(Hmisc)
fish.coef <- sapply(1:26, function(i) {
  y <- fish[fish.idx, i] != 0
  temp <- glm(y~log(fish$GAPDH[fish.idx]), family = binomial(link="logit"))
  summary(temp)$coefficients[2, 1:2]
})
p0 <- colMeans(fish[fish.idx, ] == 0, na.rm = T)
rm.idx <- c(4, 11, 19)
ylim <- range(c(fish.coef[1, -rm.idx] - fish.coef[2, -rm.idx], fish.coef[1, -rm.idx] + fish.coef[2, -rm.idx]))
### Figure cell size effect on burst intensity
require(gamlss)
require(gamlss.tr)
require(Hmisc)
Trunzero.NB <- trun(par = 0, family = NBI, 
                    local = F)
fish.coef.size <- sapply(1:26, function(i) {
    print(i)
    y <- fish[fish.idx, i]
    x <- fish$GAPDH[fish.idx]
    idx.y <- which(y > 0)
    temp <- gamlss(y[idx.y] ~ log(x[idx.y]), family = Trunzero.NB,
                   control = gamlss.control(trace = F))
    summary(temp)[2, 1:2]
})
rm.idx <- c(4, 11, 19)
fish.mean <- colMeans(fish, na.rm = T)


## plot the effects of cell size on Nonzero Fraction/Mean
layout(matrix(1:2, nrow = 1))
plot(1 - p0[-rm.idx], fish.coef[1, -rm.idx], ylim = c(-2, ylim[2]),
     xlab = "Nonzero Fraction", 
     ylab = "Coefficient of Cell Size on Nonzero Fraction",
     col = "dodgerblue1", pch = 15, cex = 1.2, main = "RNA FISH")
errbar(1 - p0[-rm.idx], fish.coef[1, -rm.idx],
       fish.coef[1, -rm.idx] + fish.coef[2, -rm.idx],
       fish.coef[1, -rm.idx] - fish.coef[2, -rm.idx], add = T, pch = ".")
abline(h = 0, col = "blue", lty = 2)

plot(log(fish.mean[-rm.idx]), fish.coef.size[1, -rm.idx], ylim = c(-1, 3),
     xlab = "log of Nonzero Mean", 
     ylab = "Coefficient of Cell Size on Nonzero Mean",
     col = "dodgerblue1", pch = 15, cex = 1.2, main = "RNA FISH")
errbar(log(fish.mean[-rm.idx]), fish.coef.size[1, -rm.idx],
       fish.coef.size[1, -rm.idx] + fish.coef.size[2, -rm.idx],
       fish.coef.size[1, -rm.idx] - fish.coef.size[2, -rm.idx], add = T, pch = ".")
abline(h = 1, col = "blue", lty = 2)
dev.off()




