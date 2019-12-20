# visulization functions

#' visulize stability measurements of leaf nodes
#' @export
plotStability <- function(stability.measure){

  p.fai <- ggplot2::ggplot(data.frame(aRI=stability.measure$flat$ari), ggplot2::aes(x=1,y=aRI)) +
    ggplot2::geom_boxplot(notch=T,outlier.shape=NA) +
    ggplot2::geom_point(shape=16, position = ggplot2::position_jitter(), alpha=0.3) +
    ggplot2::guides(color=FALSE) +
    ggplot2::geom_hline(yintercept=1, linetype="dashed", alpha=0.2) +
    ggplot2::ylim(c(0,1)) + ggplot2::labs(x=" ", y="adjusted Rand Index") +
    ggplot2::theme(legend.position="none", axis.ticks.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank()) +
    ggplot2::theme_bw()

  df <- reshape2::melt(stability.measure$flat$jc)
  colnames(df) <- c('rep','cluster','jc')
  df$cluster <- factor(colnames(stability.measure$flat$jc)[df$cluster], levels=unique(as.numeric(colnames(stability.measure$flat$jc)[df$cluster])))

  nclusters <- length(unique(df$cluster))
  p.fjc <- ggplot2::ggplot(df, ggplot2::aes(x=cluster,y=jc,color=cluster)) +
    ggplot2::geom_boxplot(ggplot2::aes(color=cluster),notch=T,outlier.shape=NA) +
    ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.2),alpha=0.3) +
    ggplot2::guides(color=FALSE) +
    ggplot2::geom_hline(yintercept=1, linetype="dashed", alpha=0.2) +
    ggplot2::ylab("Jaccard coefficient (flat)") + ggplot2::ylim(c(0,1)) + ggplot2::theme_bw()
  cowplot::plot_grid(plotlist=list(p.fai,p.fjc),nrow=1, rel_widths=c(4,nclusters,nclusters))

}

#' plot 3-factor legend
#' @export
plot3factorLengend <- function(rate=15, base=0.001){
  sm <- 500
  hc <- tan(pi/3)/2 # sqrt(3)/2
  x <- do.call(c, sapply(1:(sm*sqrt(3)/2)/sm,  function(i) (i*sm/sqrt(3)):(sm-i*sm/sqrt(3))/sm))
  y <- do.call(c, sapply(1:(sm*sqrt(3)/2)/sm,  function(i) rep(i, length((i*sm/sqrt(3)):(sm-i*sm/sqrt(3))))))
  d.red = y/hc
  d.green= pmax(0,sin(pi/3)* (x-y/tan(pi/3))/hc)
  d.blue= pmax(0,sin(-pi/3)* (x-(y+tan(-pi/3))/tan(-pi/3))/hc)

  # normalize in the same way
  cm <- cbind(d.red,d.green,d.blue)

  cm <- dexp(cm, rate)
  #cm <- cbind(dexp(cm[,1], 0.3), dexp(cm[,2], 0.4), dexp(cm[,3], 0.3))
  cm <- cm/rate * (1-base)

  library(Cairo)
  #CairoPNG(file='tri.legend.png',width=300,height=300)
  plot(NA,NA,xlim=c(0,1),ylim=c(0,1),asp=1,bty="n",axes=F,xlab="",ylab="")
  lines(c(0,1,0.5,0),c(0,0,hc,0),col=1,lwd=8)
  col <- rgb(cm[,1], cm[,2], cm[,3])
  col <- adjustcolor(col, offset = c(0.5, 0.5, 0.5, 0))
  points(x, y, col=col, pch=19, cex=0.5)
  #dev.off()
}
