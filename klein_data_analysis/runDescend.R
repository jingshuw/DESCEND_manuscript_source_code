require(parallel)

## the way to create cl depends on own your cluster parallization setup ##############
outfile <- paste("verbose_log_", as.numeric(Sys.time()), ".txt", sep = "")
cl <- makeCluster(100, outfile = outfile)

library(descend)

## This depend on your own data file folder and file name
Y0 <- read.csv("GSM1599494_ES_d0_main.csv",
              header = F, row.names = 1)
Y2 <- read.csv("GSM1599497_ES_d2_LIFminus.csv",
                  header = F, row.names = 1)
Y4 <- read.csv("GSM1599498_ES_d4_LIFminus.csv",
                  header = F, row.names = 1)
Y7 <- read.csv("GSM1599499_ES_d7_LIFminus.csv",
                  header = F, row.names = 1)


result1 <- runDescend(Y0, scaling.consts = colSums(Y0), cl = cl, 
                      control = list(max.sparse = c(0.95, 20)))

saveRDS(result1, file = "descend_result_day0.rds")

result2 <- runDescend(Y2, scaling.consts = colSums(Y2), cl = cl, 
                      control = list(max.sparse = c(0.95, 20)))

saveRDS(result1, file = "descend_result_day2.rds")

result3 <- runDescend(Y4, scaling.consts = colSums(Y4), cl = cl, 
                      control = list(max.sparse = c(0.95, 20)))

saveRDS(result1, file = "descend_result_day4.rds")

result4 <- runDescend(Y7, scaling.consts = colSums(Y7), cl = cl, 
                      control = list(max.sparse = c(0.95, 20)))

saveRDS(result1, file = "descend_result_day7.rds")

library(reldist)
metrics.dropseq <- lapply(list(Y0, Y2, Y4, Y7), function(Y) {
                    ginis <- apply(t(t(Y)/colSums(Y)), 1, gini)
                    pos.frac <- rowMeans(Y > 0)
                    lib.size <- colSums(Y)
                    return(list(lib.size = lib.size, pos.frac = pos.frac, Gini = ginis))
                  })
saveRDS(metrics.dropseq, file = "dropseq_raw_metrics.rds")


Y <- cbind(Y0, Y2, Y4, Y7)
lib.size <- colSums(Y)


result.list <- list(result1, result2, result3, result4)
names(result.list) <- c("Day0", "Day2", "Day4", "Day7")

labels <- c(rep("Day0", ncol(Y0)),
            rep("Day2", ncol(Y2)),
            rep("Day4", ncol(Y4)),
            rep("Day7", ncol(Y7)))

model <- list(scaling.consts = lib.size, Z= NULL, Z0 = NULL, control = list(max.sparse = c(0.95, 20)),
              family = "Poisson", NB.size = 100)

multi.result <- list(descend.list.list = result.list, model = model)

test.result <- deTest(multi.result, c("Day0", "Day2"),
                      count.matrix = Y, labels = labels,
                      N.genes.null = nrow(Y), cl = cl)

saveRDS(test.result, file = "descend_test_result_Day02.rds")



