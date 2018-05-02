bio.g <- readRDS("data.rds")
source("functions.R")
lib.size <- colSums(bio.g)

library(descend)

labels <- sapply(colnames(bio.g), 
                 function(v)paste(strsplit(v, "[.]")[[1]][1:2], 
                                  collapse = "."))
labels <- as.factor(labels)
sel.g <- rowMeans(bio.g) > 50
bio.g <- bio.g[sel.g, ]

require(parallel)
## the setup and generation of cl depend on the cluster setting
outfile <- paste("verbose_log_", as.numeric(Sys.time()), ".txt", sep = "")
cl <- makeCluster(10, outfile = outfile)

############### Differential Testing Between Batches ###########################################
set.seed(10)
select.idx <- unlist(lapply(levels(labels), function(ll) {
                                  n <- sum(labels == ll)
                                  idx <- rep(0, n)
                                  idx[sample(1:n, 50)] <- 1
                                  return(idx)
                                 }))
labels1 <- rep(0, length(labels))
labels1[labels %in% levels(labels)[c(1, 4, 6)] & select.idx] <- "Group1"
labels1[labels %in% levels(labels)[c(2, 3, 7)] & select.idx] <- "Group2"

Y <- bio.g[, labels1 != 0]
label.group <- labels1[labels1 != 0]
lib.size.group <- lib.size[labels1 != 0]
## There are in total 6 batches here, but as the batches are perfectly confounded with the groups
## Batch adjustment on the mean is impossible, but on CV/Gini is doable
batches1 <- model.matrix(~factor(labels[labels1 == "Group1"]))[, -1]
batches2 <- model.matrix(~factor(labels[labels1 == "Group2"]))[, -1]
## used for permutation
batches.comb <- model.matrix(~factor(labels[labels1 != 0]))[, -1]

batches <- batches.temp
batches[label.group == "Group1", ] <- batches.temp[1:150, ]
batches[label.group == "Group2", ] <- batches.temp[151:300, ]


## Without coviariates adjustment
result.wo <- descendMultiPop(Y, label.group, scaling.consts = lib.size.group,
                             cl = cl)
test.wo <- deTest(result.wo, unique(label.group), 
                  Y, label.group, N.genes.null = nrow(Y) * 5,
                  cl = cl)

result1 <- runDescend(Y[, label.group == "Group1"], 
                      scaling.consts = lib.size.group[label.group == "Group1"],
                      Z = batches1, cl = cl)
result2 <- runDescend(Y[, label.group == "Group2"], 
                      scaling.consts = lib.size.group[label.group == "Group2"],
                      Z = batches2, cl = cl)
result.wi <- list(descend.list.list = list(Group1 = result1, Group2 = result2))

model <- list(scaling.consts = lib.size.group,
              Z = batches.comb, 
              Z0 = NULL,
              control = list(), family = "Poisson", NB.size = 100)
test.wi <- deTest(result.wi, unique(label.group), Y,
                  label.group, N.genes.null = nrow(Y) * 5, params = model, cl = cl)

results <- list(result.wo = result.wo,
               result.wi = result.wi,
               test.wo = test.wo,
               test.wi = test.wi)

saveRDS(results, file = "descend_results_between_replicates.rds")


################## Differential testing between individuals ###################################
set.seed(10)
labels1 <- rep(0, length(labels))
labels1[labels %in% levels(labels)[3:5]] <- "NA19101"
labels1[labels %in% levels(labels)[6:8]] <- "NA19239"
Y <- bio.g[, labels1 != 0]
label.indiv <- labels1[labels1 != 0]
lib.size.indiv <- lib.size[labels1 != 0]
## There are in total 6 batches here, but as the batches are perfectly confounded 
## with the cell labels, we define the batch matrix this way simply to avoid computational issues
batches <- rbind(model.matrix(~factor(labels[labels1 == "NA19101"]))[, -1],
                 model.matrix(~factor(labels[labels1 == "NA19239"]))[, -1])

batches1 <- model.matrix(~factor(labels[labels1 == "NA19101"]))[, -1]
batches2 <- model.matrix(~factor(labels[labels1 == "NA19239"]))[, -1]
## used for permutation
batches.comb <- model.matrix(~factor(labels[labels1 != 0]))[, -1]




## Without coviariates adjustment
result.wo <- descendMultiPop(Y, label.indiv, scaling.consts = lib.size.indiv,
                             cl = cl)
test.wo <- deTest(result.wo, unique(label.indiv), 
                  Y, label.indiv, N.genes.null = nrow(Y) * 5,
                  cl = cl)


result1 <- runDescend(Y[, label.indiv == "NA19101"], 
                      scaling.consts = lib.size.indiv[label.indiv == "NA19101"],
                      Z = batches1, cl = cl)
result2 <- runDescend(Y[, label.indiv == "NA19239"], 
                      scaling.consts = lib.size.indiv[label.indiv == "NA19239"],
                      Z = batches2, cl = cl)
result.wi <- list(descend.list.list = list(NA19101 = result1, NA19239 = result2))

model <- list(scaling.consts = lib.size.indiv,
              Z = batches.comb, 
              Z0 = NULL,
              control = list(), family = "Poisson", NB.size = 100)
test.wi <- deTest(result.wi, unique(label.indiv), Y,
                  label.indiv, N.genes.null = nrow(Y) * 5, params = model, cl = cl)

results <- list(result.wo = result.wo,
               result.wi = result.wi,
               test.wo = test.wo,
               test.wi = test.wi)

saveRDS(results, file = "descend_results_between_individuals.rds")







