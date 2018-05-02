source("functions_all.R")
library(descend)

## need to specify DATAFILEPATH, trueMol and alpha for the code to run, see the file functions_all.R
## get the trueMol and alpha of Zeisel data
setwd("../spike_in_analysis")
source("functions_ercc.R")
library(Hmisc)

library(grid)
library(gridExtra)

ercc.info <- read.table("ercc-info.txt", header = T, sep = "\t")
data.info <- read.csv("ERCC_datasets.csv", stringsAsFactors = F)
ercc.data <- getData(2)
alpha <- ercc.data$alpha
trueMol <- ercc.data$trueMol
setwd("../zeisel_data_analysis/")

#DATAFILEPATH <- "GSE60361_C1-3005-Expression.txt"
data <- GetY(DATAFILEPATH, trueMol, alpha)
alpha <- data$alpha
log.lib.size <- log(colSums(data$Y))
Ctype <- levels(data$Ctype)

## instead of calculating cell size directly as lib.size/alpha, we put both lib.size and alpha (after log transformation) into the model and coefficient of lib.size then will be the same as coefficient of cell size.
Z <- cbind(log.lib.size, log(alpha))
## should center the data to make cell size adjusted values meaningful
Z <- apply(Z, 2, scale, scale = F)


require(parallel)
## the setup and generation of cl depend on the cluster setting
outfile <- paste("verbose_log_", as.numeric(Sys.time()), ".txt", sep = "")
cl <- makeCluster(100, outfile = outfile)

################ DESCEND result on 3005 cells together used for HVG selection #######################
computeForHVG <- function(cl) {

  Y <- data$Y
  ## This is the library size
  scaling <- colSums(Y)
  result <- runDescend(Y, scaling.consts = scaling,
                       cl = cl)
  ests <- getEstimates(result)

 saveRDS(ests,  
         file = "../result_data/all_cells_ests_forHVG.rds")
}


######## DESCEND results for each major cell type with cell size adjustment ########################
## we apply runDESCEND to each cell type separately, now you can use the function descendMultiPop to reproduce the result by setting control = list(max.sparse = c(0.95, 20))
computeMultiPop <- function(cl) {
  if (!file.exists("temp"))
    dir.create("temp")
  print("Directory ./temp will be used to store DESCEND result for each cell type!")

  result.list <- lapply(Ctype, function(str) {
                              print(paste("DESCEND is applied to cell population", 
                                          str, "!"))
                              idx <- which(data$Ctype == str)
                              sel.g <- (rowMeans(data$Y[, idx] == 0) < 0.95) & 
                              (rowMeans(data$Y[, idx]) > 0.3)

                              Z.temp <- Z[idx, ]
                              Z0.temp <- Z.temp

                           
                              result <- runDescend(data$Y[sel.g, idx],
                                                   ercc.matrix = NULL,
                                                   scaling.consts = rep(1, length(idx)),
                                                   Z = Z.temp,
                                                   Z0 = Z0.temp,
                                                   cl = cl, 
                                                   do.LRT.test = T,
                                                   control = list(LRT.Z.select = c(T, F),
                                                                  LRT.Z0.select = c(T, F),
                                                                  LRT.Z.values = 1))
                              save(result, file = paste("temp/DESCEND_result_", str, ".rda", sep = ""))
                              rm(result)
                              gc()
                              print(paste("Results for", str, "saved."))
                              return(list())
                            })

  result.list <- lapply(Ctype, function(str) {
                             load(paste("temp/DESCEND_result_", str, ".rda", sep = ""))
                             return(result)
                            }) 
  names(result.list) <- Ctype

  ## compute Z0 adjusted nonzero fraction
  result.list <- lapply(result.list, 
                        function(ll) {
                          result <- lapply(ll, 
                                           function(lll) {
                                             if (class(lll) != "DESCEND")
                                               return(NA)
                                             temp <- lll@estimates["Z0 effect: beta0", ]
                                             p0.adj <- exp(temp[1])/(1 + exp(temp[1]))
                                             p0.adj.bias <- temp[2] * p0.adj * 
                                             (1 - p0.adj)
                                             p0.adjsd <- temp[3] * p0.adj * (1 - p0.adj)
                                             temp <- c(1- p0.adj, -p0.adj.bias, p0.adjsd,
                                                       sqrt(p0.adj.bias^2 + p0.adjsd^2))

                                             lll@estimates <- rbind(lll@estimates[1:5,,drop = F],
                                                                    temp,
                                                                    lll@estimates[-(1:5),, drop = F])
                                             rownames(lll@estimates)[6] <- "Z0 Adjusted Nonzero Fraction"
                                             return(lll)

                                           })
                        })
  model <- list(scaling.consts = rep(1, length(ncol(data$Y))),
                Z = Z, Z0 = Z, control = list(), family = "Poisson", NB.size = NULL)
  saveRDS(list(descend.list.list = result.list,
               model = model), file = "cell_size_multi_cell_type_descend_results.rds")
}



####################### Differential testing between cell types with cell size adjustment ####################
DETestCompute <- function(cl) {
  result.multi <- readRDS("cell_size_multi_cell_type_descend_results.rds")
  multi.output <- result.multi$descend.list.list[c(2, 6)]
  common.genes.names <- intersect(names(multi.output[[1]])[sapply(multi.output[[1]], 
                                                                  function(ll) class(ll) == "DESCEND")], 
                                  names(multi.output[[2]])[sapply(multi.output[[2]], 
                                                                  function(ll) class(ll) == "DESCEND")])
  multi.output[[1]] <- multi.output[[1]][common.genes.names]
  multi.output[[2]] <- multi.output[[2]][common.genes.names]
  Y <- data$Y[common.genes.names, data$Ctype %in% names(multi.output)]
  model <- result.multi$model
  model$scaling.consts <- model$scaling.consts[data$Ctype %in% names(multi.output)]
  model$Z <- model$Z[data$Ctype %in% names(multi.output),]
  model$Z0 <- model$Z0[data$Ctype %in% names(multi.output),]

  DEresult <- deTest(list(descend.list.list = multi.output, model = model), 
                     names(multi.output),
                     Y, labels = data$Ctype[data$Ctype %in% names(multi.output)],
                     N.genes.null = length(common.genes.name), cl = cl)

  saveRDS(DEresult,
          file = "DEtest_celltypes_2_6_descend_results.rds")
}


################ get DESCEND results after sample splitting ############################
## i = 5: the Oligodendrocyte cell type
computeSplitSimu <- function(i = 5, cl) {
  idx <- (data$Ctype == Ctype[i]) & (alpha != 0)
  cell.sample <- which(idx)
  temp <- rbinom(sum(idx), 1, 0.5)
  idx1 <- cell.sample[temp == 0]
  idx2 <- cell.sample[temp == 1] 
  Y1 <- data$Y[, idx1]
  Y2 <- data$Y[, idx2]



  result1 <- runDescend(Y1, scaling.consts = alpha[idx1],
                        control = list(only.integer = T), cl = cl, control = list(max.sparse = c(0.98, 20)))

  result2 <- runDescend(Y2, scaling.consts = alpha[idx1],
                        control = list(only.integer = T), cl = cl, control = list(max.sparse = c(0.98, 20)))

  saveRDS(list(result1 = result1, result2 = result2),
          file = "splitting_descend_results.rds")
}


############## get DESCEND result for parametric simulation #################################
computeParametricSimu <- function(cl) {

  simu.data <- readRDS("parametric_simulated_data.rds")

  ## As both the library size and efficiency parameters are added as covariates 
  ## for Z0, the coefficient of the library size is the coefficient of cell size
  result <- runDescend(simu.data$Y[1:5, ],
                       scaling.consts = exp(simu.data$covariates[, 2]),
                       Z = simu.data$covariates[, 2],
                       Z0 = simu.data$covariates, do.LRT.test = T,
                       n.cores = 3)

  saveRDS(result, file = "parametric_simu_descend_results.rds")
}


############ get DESCEND result for Down-sampling experiment #########################
GetDownSampleData <- function(i = 5) {
  set.seed(1)
  idx <- data$Ctype == Ctype[i]
  Y <- data$Y[, idx]
  avg.umi <- apply(Y, 1, mean)
  sel.genes <- names(sort(avg.umi, decreasing = T)[1:150])
  CV <- apply(Y[sel.genes, ], 1, function(v)sd(v)/mean(v))
  sel.genes <- sel.genes[CV < 4]


  Y <- Y[sel.genes, ]
  avg.umi.genes <- avg.umi[sel.genes]

  require(reldist)
  Gini <- apply(Y, 1, gini)
  CV <- apply(Y, 1, function(v)sd(v)/mean(v))
  nonzero.frac <- apply(Y, 1, function(v)mean(v !=0))

  Y20 <- matrix(rpois(nrow(Y) * ncol(Y), 0.2 * as.vector(Y)), nrow = nrow(Y))

  result20 <- runDescend(Y20, 
                        scaling.consts = rep(0.2, ncol(Y)),
                        cl = cl)
  ests20 <- getEstimates(result20)

  Y10 <- matrix(rpois(nrow(Y) * ncol(Y), 0.1 * as.vector(Y)), nrow = nrow(Y))

  result10 <- runDescend(Y10, 
                        scaling.consts = rep(0.1, ncol(Y)),
                        cl = cl)
  ests10 <- getEstimates(result10)


  Y5 <- matrix(rpois(nrow(Y) * ncol(Y), 0.05 * as.vector(Y)), nrow = nrow(Y))

  result5 <- runDescend(Y5, 
                        scaling.consts = rep(0.05, ncol(Y)),
                        cl = cl)
  ests5 <- getEstimates(result5)

  Gini.diff <- cbind(Gini, 
                     ests20[[5]][, 1] - Gini, 
                     ests10[[5]][, 1] - Gini,
                     ests5[[5]][, 1] - Gini,
                     apply(Y20, 1, gini) - Gini,
                     apply(Y10, 1, gini) - Gini,
                     apply(Y5, 1, gini) - Gini)
  colnames(Gini.diff) <- c("true.gini", "20descend.err", "10descend.err", "5descend.err",
                           "20raw.err", "10raw.err", "5raw.err")

  CV.diff <- cbind(CV, 
                   ests20[[4]][, 1] - CV, 
                   ests10[[4]][, 1] - CV,
                   ests5[[4]][, 1] - CV,
                   apply(Y20, 1, function(v)sd(v)/mean(v)) - CV,
                   apply(Y10, 1, function(v)sd(v)/mean(v)) - CV,
                   apply(Y5, 1, function(v)sd(v)/mean(v)) - CV)
  colnames(CV.diff) <- c("true.CV", "20descend.err", "10descend.err", "5descend.err",
                           "20raw.err", "10raw.err", "5raw.err")

  nonzero.diff <- cbind(nonzero.frac, 
                        ests20[[1]][, 1] - nonzero.frac, 
                        ests10[[1]][, 1] - nonzero.frac,
                        ests5[[1]][, 1] - nonzero.frac,
                     apply(Y20, 1, function(v)mean(v !=0)) - nonzero.frac,
                     apply(Y10, 1, function(v)mean(v != 0)) - nonzero.frac,
                     apply(Y5, 1, function(v)mean(v != 0)) - nonzero.frac)
  colnames(nonzero.diff) <- c("true.nonzero", "20descend.err", "10descend.err", "5descend.err",
                              "20raw.err", "10raw.err", "5raw.err")

  result <- list(Gini = Gini.diff, CV = CV.diff, nonzero.fraction = nonzero.diff) 
  save(result,
       file = "downsampling_simulation.rds")
}


