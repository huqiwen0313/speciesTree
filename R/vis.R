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
