plotShapeSingle <- function(i, add.legend = T) {
  if (class(result[[i]]) != "DESCEND")
    return(list())
  fish.temp <- t(t(fish)/colMeans(fish, na.rm = T))
  temp <- fish.temp[fish.temp[, i] < quantile(fish.temp[, i], 0.999, na.rm = T), i]
  dropseq.temp <- as.numeric(dropseq.fish[i, ])/lib.size
  dropseq.temp <- dropseq.temp[dropseq.temp < quantile(dropseq.temp, 0.98, na.rm = T)]
  dropseq.temp <- dropseq.temp / mean(dropseq.temp)

  bg.pieces <- theme_bw() + 
              theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              axis.line = element_line(colour = "black"),
              panel.border = element_blank(),
              plot.background=element_blank())


  data1 <- data.frame(x = temp)
  data2 <- data.frame(x = dropseq.temp)

  ## scale the DESCEND density to have the same mean as FISH for comparison
  scale.fac <- result[[i]]@estimates[3, 1]

  dens <- result[[i]]@density.points
  dens[, 1] <- dens[, 1]/scale.fac
  dens[, 2] <- dens[, 2] * scale.fac


  data3 <- data.frame(dens)


  tm <- data.frame(Marker = sort(rep(letters[1:3], 10)), x = runif(30), y = runif(30),
                   type = sort(rep(letters[c(1, 1, 2)], 10)))

  p1 <- ggplot() + geom_density(data = data1, aes(x), size = 0.9) + 
  geom_density(data = data2, aes(x), colour = "blue",
                                                                    lty = 2, size = 0.9) +  
      xlab("") + ylab("") + 
      xlim(min(temp, na.rm = T), max(temp, na.rm = T)) + 
      geom_line(data = data3, aes(x = theta, y= density), col = "red", size = 0.9) + bg.pieces 

  if (!add.legend)
    p <- p1 + ggtitle(colnames(fish.temp)[i]) +   
    theme(legend.position=c(1, 1), legend.justification = c(1, 1), legend.title = element_blank(),
          axis.text = element_text(size = 14), 
          plot.title = element_text(margin = margin(b = -15), hjust = 0.5, size = 18),
          plot.margin = unit(c(4, 5, -7, -7)/100, "npc"))
  else
    p <- p1 + geom_line(data = tm, aes(x = x, y = y, col = Marker, linetype = Marker), alpha = 0) + 
    ggtitle(colnames(fish.temp)[i]) +
    scale_linetype_manual(name = "Marker", values = c("solid", "solid", "dashed"),
                          labels = c("FISH", "DESCEND", "Drop-seq")) + 
scale_colour_manual(name = "Marker", values = c("BLACK", "RED", "BLUE"),
                    labels = 
                    c("FISH", "DESCEND", "Drop-seq")) + 
guides(linetype = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
theme(legend.position=c(0.9, 0.8), legend.justification = c(1, 1), legend.title = element_blank(),
      legend.text =element_text(size = 14),
      plot.title = element_text(margin = margin(b = -15), hjust = 0.5, size = 18),
      plot.margin = unit(c(4, 5, -7, -7)/100, "npc"),
      axis.text = element_text(size = 15))


  return(p)

  
}


plotCVGini <- function() {

  require(ggplot2)
  require(reldist)

  true.val <- apply(fish, 2, function(v)
                    c(sd(v, na.rm = T)/mean(v, na.rm = T), gini(v[!is.na(v)])))
  est.val <- sapply(ests, function(l) l[, 1])
  est.val <- t(est.val[, c(4, 5, 1)])

  est.val.sd <- sapply(ests, function(l) l[, 2])
  est.val.sd <- t(est.val.sd[, c(4, 5, 1)])



  normed.seq <- t(dropseq.fish)/lib.size
  seq.val1 <- apply(normed.seq, 2, function(v)
                    c(sd(v, na.rm = T)/mean(v, na.rm = T)))

  cv.lib <- sd(lib.size)/mean(lib.size)
  seq.val <- apply(dropseq.fish, 1, function(v) {
                     cv <- sqrt((sd(v)^2/mean(v)^2 - 1/mean(v) - cv.lib^2)/(1 + cv.lib^2))
                     return(c(cv, gini(v[!is.na(v)])))
                    })
  idx <- !is.na(est.val[1, ]) 
  idx[4] <- F

  #  library(cowplot)

  bg.pieces <- geom_abline(slope = 1, intercept = 0) + 
    theme_bw() + 
    theme(legend.position="none",
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border = element_blank(),
          plot.background=element_blank())




    temp <- data.frame(true = true.val[2, idx], est = est.val[2, idx], 
                       sd = est.val.sd[2, idx])
    p1 <- ggplot(temp, aes(x = true, y = est)) + 
      geom_errorbar(aes(ymin = est -sd, ymax = est+sd), colour = "dodgerblue1", 
                    width = .04, size = 1) + 
    geom_point(colour = "dodgerblue1", pch = 17, cex = 3) +  
    scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
    scale_y_continuous(breaks= c(0, 0.5, 1), limits = c(0, 1)) + 
    #xlim(0, 1) + ylim(0, 1) + 
    xlab("") + ylab("Estimated") +  
    ggtitle("DESCEND") + geom_abline(slope = 1, intercept = 0) +
    theme_bw() + 
    theme(legend.position="none",
          text = element_text(size = 19),
          axis.text = element_text(size = 19),
          plot.title = element_text(margin = margin(b=-10), hjust = 0.07, vjust = -0.5),
          plot.margin = unit(c(1, 5, 0, 1)/100, "npc"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          plot.background=element_blank())



    temp <- data.frame(true = true.val[2, idx], est = seq.val[2, idx])
    p2 <- ggplot(temp, aes(x = true, y = est)) + 
    geom_point(colour = "dodgerblue1", pch = 17, cex = 3) +  
    scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
    scale_y_continuous(breaks= c(0, 0.5, 1), limits = c(0, 1)) + 
    #xlim(0, 1) + ylim(0, 1) + 
    xlab("FISH") + ylab("Estimated") +  
    ggtitle("Drop-seq") + geom_abline(slope = 1, intercept = 0) + 
    theme_bw() + 
    theme(legend.position="none",
          text = element_text(size = 19),
          axis.text = element_text(size = 19),
          plot.title = element_text(margin = margin(b=-10), hjust = 0.07),
          plot.margin = unit(c(1, 5, 0, 1)/100, "npc"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          plot.background=element_blank())




    temp <- data.frame(true = true.val[1, idx], est = est.val[1, idx], 
                   sd = est.val.sd[1, idx])
    p3 <- ggplot(temp, aes(x = true, y = est)) + 
    geom_errorbar(aes(ymin = est -sd, ymax = est+sd), colour = "dodgerblue1", 
                  width = .1, size = 1) + 
    geom_point(colour = "dodgerblue1", pch = 17, cex = 3) +  
    scale_x_continuous(breaks = c(0, 1.2, 2.5), limits = c(0, 2.5)) +
    scale_y_continuous(breaks= c(0, 1.2, 2.5), limits = c(0, 2.5)) + 
    #xlim(0, 1) + ylim(0, 1) + 
    xlab("") + ylab("Estimated") +  
    ggtitle("DESCEND") + geom_abline(slope = 1, intercept = 0) + 
    theme_bw() + 
    theme(legend.position="none",
          text = element_text(size = 19),
          axis.text = element_text(size = 19),
          plot.title = element_text(margin = margin(b=-10), hjust = 0.07),
          plot.margin = unit(c(1, 5, 0, 1)/100, "npc"),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())



    temp2 <- data.frame(x = true.val[1, idx], y = seq.val[1, idx])
    temp <- data.frame(true = true.val[1, idx], est = seq.val1[idx])

    p4 <- ggplot(temp, aes(x = true, y = est)) + 
    geom_point(colour = "dodgerblue1", pch = 17, cex = 3) +  
    scale_x_continuous(breaks = c(0, 1.2, 2.5), limits = c(0, 2.5)) +
    scale_y_continuous(breaks= c(0, 4.0, 8), limits = c(0, 8.5)) + 
    #xlim(0, 1) + ylim(0, 1) + 
    xlab("") + ylab("Estimated") +  
    ggtitle("Drop-seq") + geom_abline(slope = 1, intercept = 0) + 
    theme_bw() + 
    theme(legend.position="none",
          text=  element_text(size = 19),
          axis.text = element_text(size = 19),
          plot.title = element_text(margin = margin(b=-10), hjust = 0.07),
          plot.margin = unit(c(1, 5, 0, 1)/100, "npc"),
          panel.grid.major=element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())

    p7 <- ggplot(temp, aes(x = true, y = est)) + 
    geom_point(data = temp2, aes(x = x, y = y), colour = "dodgerblue1", pch = 17, cex = 3) +
    scale_x_continuous(breaks = c(0, 1.2, 2.5), limits = c(0, 2.5)) +
    scale_y_continuous(breaks= c(0, 1.2, 2.5), limits = c(0, 2.5)) + 
    #xlim(0, 1) + ylim(0, 1) + 
    xlab("FISH") + ylab("Estimated") +  
    ggtitle("Variance \n Decomposition") + geom_abline(slope = 1, intercept = 0) + 
    theme_bw() + 
    theme(legend.position="none",
          text=  element_text(size = 19),
          axis.text = element_text(size = 19),
          plot.title = element_text(margin = margin(b=-10), hjust = 0.07),
          plot.margin = unit(c(1, 5, 0, 1)/100, "npc"),
          panel.grid.major=element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())





    p0 <- colMeans(fish == 0, na.rm = T)

    temp <- data.frame(true = 1 - p0[idx], est = est.val[3, idx], sd = est.val.sd[3, idx])
    p5 <- ggplot(temp, aes(x = true, y = est)) + 
    geom_errorbar(aes(ymin = est -sd, ymax = est+sd), colour = "dodgerblue1", 
                  width = .05, size = 1.2) + 
    geom_point(colour = "dodgerblue1", pch = 17, cex = 3) +  
    scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
    scale_y_continuous(breaks= c(0, 0.5, 1), limits = c(0, 1)) + 
    #xlim(0, 1) + ylim(0, 1) + 
    xlab("") + ylab("Estimated") +  
    ggtitle("DESCEND") + geom_abline(slope = 1, intercept = 0) + 
    theme_bw() + 
    theme(legend.position="none",
          text = element_text(size = 19),
          axis.text = element_text(size = 19),
          plot.title = element_text(margin = margin(b=-10), hjust = 0.07),
          plot.margin = unit(c(1, 5, 0, 1)/100, "npc"),
          panel.grid.major=element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())


    result.qvark <- read.table("QVARKS_results.txt", header = T, stringsAsFactor = F)

    p0 <- colMeans(fish == 0, na.rm = T)

    p0.seq <- rowMeans(dropseq.fish == 0)
    p0.seq[!idx] <- NA

    temp <- data.frame(true = 1 - p0, est = 1 - p0.seq)
    temp2 <- data.frame(true = 1 - p0[result.qvark$gene], est = result.qvark$ppi.c1, 
                        sd = result.qvark$ppi.c1.sd)
    p6 <- ggplot(temp2, aes(x = true, y = est)) + 
    geom_point(data = temp, aes(x = true, y = est), colour = "dodgerblue1", pch = 17, cex = 3) +  
    scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
    scale_y_continuous(breaks= c(0, 0.5, 1), limits = c(0, 1)) + 
    #xlim(0, 1) + ylim(0, 1) + 
    xlab("") + ylab("Estimated") +  
    ggtitle("Drop-seq") + geom_abline(slope = 1, intercept = 0) + 
    theme_bw() + 
    theme(legend.position="none",
          text = element_text(size = 19),
          axis.text = element_text(size = 19),
          plot.title = element_text(margin = margin(b=-10), hjust = 0.07),
          plot.margin = unit(c(1, 5, 0, 1)/100, "npc"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          plot.background=element_blank())

  p8 <- ggplot(temp2, aes(x = true, y = est)) + 
  geom_point(colour = "dodgerblue1", cex = 3, pch = 17) + 
  geom_errorbar(aes(ymin = est -sd, ymax = est+sd), colour = "dodgerblue1", 
                width = .05, size = 1.2, linetype = 1) + 
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  scale_y_continuous(breaks= c(0, 0.5, 1), limits = c(0, 1)) + 
  #xlim(0, 1) + ylim(0, 1) + 
  xlab("FISH") + ylab("Estimated") +  
  ggtitle("QVARKS") + geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + 
  theme(legend.position="none",
        text = element_text(size = 19),
        axis.text = element_text(size = 19),
        plot.title = element_text(margin = margin(b=-10), hjust = 0.07),
        plot.margin = unit(c(1, 5, 0, 1)/100, "npc"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        plot.background=element_blank())


  return(list(p1, p2, p3, p4, p5, p6, p7, p8))

}

