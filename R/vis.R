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
  cm <- cbind(d.red, d.green, d.blue)

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

#' plot stability distribution based on descretized height in the tree
#' @param dend dendrogram obj
#' @export
plotStabilityDistr <- function(dend){
  height <- get_nodes_attr(dend, "height")
  #height.discretized <- discretize(height, breaks=c(0, 0.1, 0.2, 0.3, 1, max(height)), include.lowest=TRUE, method="fixed")
  height.discretized <- discretize(height, breaks=seq(0, max_depth(dend), 1), include.lowest=TRUE, method="fixed")
  stability <- get_nodes_attr(dend, "stability")
  stability.info <- data.frame(height=height.discretized, stability=stability)
  ggplot(stability.info, aes(x=height, y=stability)) +
    geom_boxplot() + theme_classic() + ggtitle("Distribution of Stability scores") +
    scale_color_grey()
}

# functions for mapping summary of homologous clusters
# from Allen
reorder_matrix <- function(matrix1, by.rows = TRUE){
  if (by.rows == TRUE) {
    conf.order <- order(apply(matrix1, 1, which.max))
    matrix1.reordered <- matrix1[conf.order, ]
  } else {
    conf.order <- order(apply(matrix1, 2, which.max))
    matrix1.reordered <- matrix1[, conf.order]
  }
}

compare_cl <- function(cl, ref.cl,
                       plot.title = NA, plot.silent = TRUE,
                       heat.colors = colorRampPalette(c("white", "grey70", "black"))(100),
                       row.cl.num = min(length(unique(cl)),
                                        length(unique(ref.cl)))) {
  library(grid)
  library(pheatmap)

  conf1 <- table(cl, ref.cl)
  conf1 <- sweep(conf1, 1, rowSums(conf1), "/")
  conf2 <- reorder_matrix(conf1)

  # Cluster co-occurence
  cl.prop.cocl <- apply(conf1, 2, function(x) {
    grid1 <- expand.grid(x, x)
    min.prop <- apply(grid1, 1, min)
  })

  cl.prop.cocl.total <- apply(cl.prop.cocl, 1, sum)
  cl.prop.cocl.m <- matrix(cl.prop.cocl.total, nrow(conf1), nrow(conf1),
                           dimnames = list(rownames(conf1), rownames(conf1)))

  ph1 <- pheatmap(conf2, cutree_rows = row.cl.num, clustering_method = "ward.D2",
                  color = heat.colors, fontsize = 6,
                  main = plot.title, silent = plot.silent)
  return(list(conf = conf2, cocl = cl.prop.cocl.m, ph = ph1))
}

#' visulization of homologous mapping - pairwise
#' @param ref.cl vector contains cell annotation and cell ids (e.g. human_cellID),
#' @param integrated.cl vector contatins integrated cluster and cell ids -
#'                        cellIDs should be matching to ref.cl
#' @param compareFac vector showing which species to compare: e.g. c("human", "marmoset")
#' @param cellorder vector contains ordered cells
#' @param upperlevelannot cells with upper-level annotation
#' @return homologous heatmap
#' @export
MappingHeatmap <- function(ref.cl, integrated.cl, compareFac, cellorder=NULL, upperlevelannot=NULL){
  # color platte
  clpalette <- function(n){
    all.colors <- colors()
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    cols <- c(cols, sample(cols))
    cols[1:n]
  }

  # maching labels
  integrated.cl <- integrated.cl[match(names(ref.cl), names(integrated.cl))]
  upperlevelannot <- upperlevelannot[names(upperlevelannot) %in% names(ref.cl)]
  upperlevelannot <- upperlevelannot[match(names(ref.cl), names(upperlevelannot))]


  cl.conf <- compare_cl(as.factor(ref.cl), as.factor(integrated.cl))
  cocl <- cl.conf$cocl
  cocl.subset <- cocl[grepl(compareFac[1], row.names(cocl)),
                      grepl(compareFac[2], row.names(cocl))]

  row.names(cocl.subset) <- sub(paste0(compareFac[1], "_"), "",
                                row.names(cocl.subset))

  colnames(cocl.subset) <- sub(paste0(compareFac[2], "_"), "",
                               colnames(cocl.subset))
  heat.colors <- colorRampPalette(c("white", "grey70", "black"))(100)

  if(!is.null(cellorder)){
    cl.order <- intersect(cellorder, row.names(cocl.subset))
    cocl.subset2 <- reorder_matrix(cocl.subset[cl.order, ], by.rows = FALSE)

  } else{
    cocl.subset2 <- reorder_matrix(cocl.subset, by.rows = FALSE)
  }

  library(fpc)
  clus.method <- "single"
  clus.num <- pamk(cocl, 1:(min(nrow(cocl), ncol(cocl)) - 1))$nc

  if(!is.null(upperlevelannot)){
    refCluster <- gsub(paste0(compareFac[1], "_"), "", ref.cl)
    refCluster <- gsub(paste0(compareFac[2], "_"), "", refCluster)
    annotation_row = data.frame(Subclass = factor(upperlevelannot[match(rownames(cocl.subset2), refCluster)],
                                                  levels=names(table(upperlevelannot))))
    annotation_col = data.frame(Subclass = factor(upperlevelannot[match(colnames(cocl.subset2), refCluster)],
                                                  levels=names(table(upperlevelannot))))
    rownames(annotation_row) <- rownames(cocl.subset2)
    rownames(annotation_col) <- colnames(cocl.subset2)

    subclassColor <- clpalette(length(table(upperlevelannot)))
    names(subclassColor) <- names(table(upperlevelannot))
    ann_colors <- list(
      Subclass = subclassColor
    )
    pheatmap(cocl.subset2, cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = annotation_col,
             annotation_row = annotation_row, color = heat.colors, annotation_colors = ann_colors,
             fontsize = 6, cellwidth = 5, cellheight = 5)
  } else{
    pheatmap(cocl.subset2, cluster_rows = FALSE, cluster_cols = FALSE, color = heat.colors,
             fontsize = 6, cellwidth = 5, cellheight = 5)
  }
}

#' plot heatmap of expression of marker genes that differentiate two children from given a dendrogram
#' @param d dendrogram obj
#' @param count gene expression matrix
#' @export
expHeatMap <- function(d, count, prefix=c("human", "mouse", "marmoset"), title=NULL){
  leftChildCells <- attr(d[[1]], "nodesinfo")
  rightChildCells <- attr(d[[2]], "nodesinfo")
  markers <- c(as.character(attr(d[[1]], "marker")$Gene), as.character(attr(d[[2]], "marker")$Gene))

  leftChildCellExp <- as.data.frame(dplyr::bind_cols(lapply(prefix, function(r){
    apply(count[, colnames(count) %in% leftChildCells[grep(r, leftChildCells)]], 1, mean)
  })))
  names(leftChildCellExp) <- paste(prefix, "1", sep="_")
  rownames(leftChildCellExp) <- rownames(count)
  leftChildCellExp[is.na(leftChildCellExp)] <- 0

  rightChildCellExp <- as.data.frame(dplyr::bind_cols(lapply(prefix, function(r){
    apply(count[, colnames(count) %in% rightChildCells[grep(r, rightChildCells)]], 1, mean)
  })))
  names(rightChildCellExp) <- paste(prefix, "2", sep="_")
  rownames(rightChildCellExp) <- rownames(count)
  rightChildCellExp[is.na(rightChildCellExp)] <- 0

  geneExp <- cbind(leftChildCellExp, rightChildCellExp)
  geneExp <- geneExp[rownames(geneExp) %in% markers, ]

  col.annot <- data.frame(group=c(rep("group1", length(prefix)), rep("group2", length(prefix))))
  rownames(col.annot) <- colnames(geneExp)
  pheatmap::pheatmap(geneExp, cluster_rows=FALSE, cluster_cols=FALSE, scale="none",
                     color=colorRampPalette(c("white", "red"))(100), annotation_col=col.annot,
                     main=title)
}

