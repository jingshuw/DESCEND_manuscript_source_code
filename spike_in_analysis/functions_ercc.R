getData <- function(i, top = 57, discard = 0) {
  print(paste("Top", top, "ERCC spike-in genes will be used!"))
  file.name <- paste(data.info[i, 1], ".txt", sep = "")
  print(paste("Analyze data from", file.name, "..."))
  Yspike <- read.table(file.name, header = T) 

  dilution <- as.numeric(data.info[i, 7])
  vol <- data.info[i, 8]

  trueMol <- ercc.info[, 3] *(10^(-18))*(6.0221417*(10^23)) / dilution / 1000 * vol

  names(trueMol) <- ercc.info$ERCC.ID
  trueMol <- trueMol[rownames(Yspike)]

  trueMol <- trueMol[order(names(trueMol))]
  Yspike <- Yspike[names(trueMol), ]

  idx <- trueMol >= sort(trueMol, decreasing = T)[min(top, length(trueMol))]

  Yspike <- Yspike[idx, ]
  trueMol <- trueMol[idx]

  alpha <- colSums(Yspike, na.rm = T)/sum(trueMol)

  thres1 <- sort(alpha)[max(1, round(length(alpha) * discard))]
  thres2 <- sort(alpha, decreasing = T)[max(1, round(length(alpha) * discard))]


  Yspike <- Yspike[, (alpha >= thres1) & (alpha <= thres2)]
  alpha <- alpha[(alpha >= thres1) & (alpha <= thres2)]

  return(list(Yspike = Yspike,
              trueMol = trueMol,
              alpha = alpha))
}


getDESCEND <- function(ercc.data) {

  results <- runDescend(ercc.data$Yspike, 
                        scaling.consts = ercc.data$alpha,
                        n.cores = 3, verbose = F)

  ests <- getEstimates(results)


  result.comb <- data.frame(trueMol = ercc.data$trueMol, 
                            p0 = 1 - ests[[1]][, 1], 
                            p0.sd = ests[[1]][, 2],
                            CV = ests[[4]][, 1],
                            CV.sd = ests[[4]][, 2],
                            gini = ests[[5]][, 1], 
                            gini.sd = ests[[5]][, 2])

  return(result.comb)
}


plotResult <- function(result.comb, theta  = 0.015, plim1 = 0.5) {



  library(ggplot2)
  library(scales)


#  par(mar = c(2,4,1, 1.2))

  result.comb$trueMolRand <- log10(result.comb$trueMol + rnorm(nrow(result.comb), mean = 0, 
                                   sd = result.comb$trueMol^1/10))
  result.comb$no <- 1 - result.comb$p0

  bg.pieces <- theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        text = element_text(size = 19),
        axis.text = element_text(size = 19),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        plot.margin = unit(c(4, 4, -18, -12), "pt"),
        plot.background=element_blank())

   fun.gini <- function(lambda) {
    pts <- rpois(50000, lambda)
    return(gini(pts))
  }
  fun.gini2 <- function(lambda) {
    pts <- rpois(50000, lambda)
    pts <- rgamma(50000, shape = 1/theta, scale = pts * theta)
    return(gini(pts))
  }

  fun3 <- function(lambda.seq) sapply(10^lambda.seq, fun.gini)
  fun4 <- function(lambda.seq) sapply(10^lambda.seq, fun.gini2)



  p <- ggplot(result.comb, aes(x = trueMolRand, y = gini)) +  
    geom_point(colour = rgb(0.3, 0.3, 0.3), size = 0.9) + 
    scale_x_continuous(breaks = c(0, 1, 3, 4, 6)) + 
       scale_y_continuous(breaks = c(0, 0.4, 0.8), limits = c(0, 0.9)) + 
        stat_function(fun = fun3, color = "red", size = 0.9) + 
        xlab("") + ylab("")

  if (theta > 0)
    p <- p + stat_function(fun = fun4, color = "blue", size= 0.9)

  p1 <- p + bg.pieces



  fun.CV <- function(lambda) sqrt(lambda)/lambda
  fun.CV2 <- function(lambda) sqrt(lambda + theta * lambda^2)/lambda

  fun1 <- function(lambda.seq) sapply(10^lambda.seq, fun.CV)
  fun2 <- function(lambda.seq) sapply(10^lambda.seq, fun.CV2)

  p <- ggplot(result.comb, aes(x = trueMolRand, y = CV)) +  
    geom_point(colour = rgb(0.2, 0.2, 0.2), size  = 0.9) + 
    scale_x_continuous(breaks = c(0, 1, 3, 4, 6)) +       
  scale_y_continuous(breaks = c(0, 0.9, 1.8), limits = c(0, 2)) + 
  stat_function(fun = fun1, color = "red", size = 0.9) + 
  xlab("") + ylab("")

if (theta > 0)
    p <- p + stat_function(fun = fun2, color = "blue", size = 0.9)

  p2 <- p + bg.pieces



 


  fun.p0 <- function(lambda) 1 - dpois(0, lambda)
  fun5 <- function(lambda.seq) sapply(10^lambda.seq, fun.p0)
  
  p <- ggplot(result.comb, aes(x = trueMolRand, y = no)) +  
    geom_point(colour = rgb(0.2, 0.2, 0.2), size=  0.9) + 
    scale_x_continuous(breaks = c(0, 1, 3, 4, 6)) + 
       scale_y_continuous(breaks = round(c(plim1, (1 + plim1)/2, 1.0), 1)
                          , limits = c(plim1, 1)) + 
        stat_function(fun = fun5, color = "red", size = 0.9) + 
        xlab("") + ylab("")

  p3 <- p + bg.pieces

  return(list(p1, p2, p3))


}


momentPlot <- function(ercc.data) {


  est.var <- sapply(1:nrow(ercc.data$Yspike), function(i) {
               sum((ercc.data$Yspike[i, ] - 
                    ercc.data$alpha * ercc.data$trueMol[i])^2)/sum(ercc.data$alpha)
                          })

  result.comb <- data.frame(trueMol = ercc.data$trueMol, 
                            CV = sqrt(est.var)/ercc.data$trueMol)

  result.comb$trueMolRand <- log10(result.comb$trueMol + rnorm(nrow(result.comb), mean = 0, 
                                   sd = result.comb$trueMol^1/10))

  fun.CV <- function(lambda) sqrt(lambda)/lambda

  fun1 <- function(lambda.seq) sapply(10^lambda.seq, fun.CV)


  p <- ggplot(result.comb, aes(x = trueMolRand, y = CV)) +  
    geom_point(colour = rgb(0.2, 0.2, 0.2), size  = 0.9) + 
    scale_x_continuous(breaks = c(0, 1, 3, 4, 6)) +       
  scale_y_continuous(breaks = c(0, 0.2, 0.4), limits = c(0, 0.5)) + 
  stat_function(fun = fun1, color = "red", size = 0.9) + 
  xlab("") + ylab("")


  bg.pieces <- theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        text = element_text(size = 19),
        axis.text = element_text(size = 19),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        plot.margin = unit(c(4, 4, -18, -6), "pt"),
        plot.background=element_blank())

  p <- p + bg.pieces
  return(list(p))
}


