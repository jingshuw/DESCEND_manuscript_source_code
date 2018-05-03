library(descend)
source("functions.R")

## load DESCEND result
results.indiv <- readRDS("descend_results_between_individuals.rds")
results.rep <- readRDS("descend_results_between_replicates.rds")
ests.wi.indiv <- lapply(results.indiv$result.wi$descend.list.list, getEstimates)
ests.wo.rep <- lapply(results.rep$result.wo$descend.list.list, getEstimates)
ests.wi.rep <- lapply(results.rep$result.wi$descend.list.list, getEstimates)


######################## Plot for Gini Batch Effect Correction ###########################
layout(matrix(1:3, nrow = 1))
par(mar = c(2, 2.7, 1, 1), mgp = c(2, 0.5, 0))
plotGini(ests.wo.rep, results.rep$test.wo$p.values)
text(sqrt(0.002), sqrt(0.27), "Between Batches", cex = 1.4, pos = 4)
text(sqrt(0.002), sqrt(0.22), "(before correction)", cex = 1.3, pos = 4)

plotGini(ests.wi.rep, results.rep$test.wi$p.values)
text(sqrt(0.002), sqrt(0.27), "Between Batches", cex = 1.4, pos = 4)
text(sqrt(0.002), sqrt(0.22), "(after correction)", cex = 1.3, pos = 4)

par(mar = c(3, 3.5, 1, 1))
plotGini(ests.wi.indiv, results.indiv$test.wi$p.values,
         xlab = "NA19101", ylab = "NA19239") 
text(sqrt(0.002), sqrt(0.24), "Between\nIndividuals", cex = 1.4, pos = 4)


######################## Plot for CV Batch Effect Correction ###########################
layout(matrix(1:3, nrow = 1))
par(mar = c(2, 2.5, 1, 1), mgp = c(2, 0.5, 0))
plotCV(ests.wo.rep, results.rep$test.wo$p.values, c(0.02, 0.5))
text(sqrt(0.01), sqrt(0.47), "Between Batches", cex = 1.4, pos = 4)
text(sqrt(0.01), sqrt(0.4), "(before correction)", cex = 1.3, pos = 4)

plotCV(ests.wi.rep, results.rep$test.wi$p.values, c(0.01, 0.5))
text(sqrt(0.005), sqrt(0.47), "Between Batches", cex = 1.4, pos = 4)
text(sqrt(0.005), sqrt(0.4), "(after correction)", cex = 1.3, pos = 4)

par(mar = c(3, 3.5, 1, 1))
#at.pt <- c(0.01, 0.1, 0.2)
plotCV(ests.wi.indiv, results.indiv$test.wi$p.values, c(0.01, 0.5),
       xlab = "NA19101", ylab = "NA19239")
text(sqrt(0.003), sqrt(0.4), "Between\nIndividuals", cex = 1.4, pos = 4)




