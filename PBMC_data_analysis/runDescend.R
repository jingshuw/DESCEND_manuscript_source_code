data <- readRDS("sampled_combined_pure_data.rds")

## depend on your own setup of clusters ##############
cl <- makeCluster(100)


Y <- as.matrix(data$Y)
scaling <- colSums(Y)

## You may need to split the rows and run multiple times as this process on the whole Y matrix takes a lot of RAM
result <- runDescend(Y, scaling.consts = scaling,
                     cl = cl, control = list(max.sparse = c(0.999, 20)))


saveRDS(result, file = "descend_result.rds")

