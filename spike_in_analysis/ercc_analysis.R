source("functions_ercc.R")
library(reldist)
library(descend)
library(Hmisc)

library(grid)
library(gridExtra)



ercc.info <- read.table("ercc-info.txt", header = T, sep = "\t")
data.info <- read.csv("ERCC_datasets.csv", stringsAsFactors = F)


p.list <- list()

## cells whose library sizes are at +- 5% tail are filtered out for Jaitin's data as these tail efficiencies either is almost 0 or is too high (low quality)
ercc.data <- getData(1, discard = 0.05)
result <- getDESCEND(ercc.data)
p.list <- c(p.list, plotResult(result, theta = 0))


ercc.data <- getData(2)
result <- getDESCEND(ercc.data)
p <- plotResult(result,, theta = 0)
p.list <- c(p.list, plotResult(result, theta = 0.015))


ercc.data <- getData(3)
result <- getDESCEND(ercc.data)
p.list <- c(p.list, plotResult(result, theta = 0))

ercc.data <- getData(4)
result <- getDESCEND(ercc.data)
p.list <- c(p.list, plotResult(result, theta = 0.015))

gc()
gc()

ercc.data <- getData(5)
result <- getDESCEND(ercc.data)
p.list <- c(p.list, plotResult(result, theta = 0))

ercc.data <- getData(6)
result <- getDESCEND(ercc.data)
p.list <- c(p.list, plotResult(result, theta = 0))

ercc.data <- getData(7)
result <- getDESCEND(ercc.data)
p.list <- c(p.list, plotResult(result, theta = 0))


ercc.data <- getData(8)
result <- getDESCEND(ercc.data)
p.list <- c(p.list, plotResult(result, theta = 0.015))



ercc.data <- getData(9)
p.list <- c(p.list, momentPlot(ercc.data))

ercc.data <- getData(10)
p.list <- c(p.list, momentPlot(ercc.data))


p.list <- c(p.list, list(grob()))



temp.fun <- function(...)grid.arrange(..., ncol = 9, as.table = F,
                                      padding = unit(1, "line"))
do.call(temp.fun, p.list)




