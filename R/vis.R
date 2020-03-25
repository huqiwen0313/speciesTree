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

#' visualize expression heatmap across species at specific node in dendrogram
#' dendrogram object
#' @param d dendrogram obj
#' @param count raw count matrix - rows are cells and columns are genes
#' @param prefix factors(species) to visulize - currently only support 2-pairwise comparison
#' @return gene expression heatmaps across factors in a specific node
#' @export
HeatmapSpecies <- function(d, count, prefix=c("human", "mouse", "marmoset")){
  cells <- attr(d, "nodesinfo")
  genes <- attr(d, "assGenes")$Gene
  group1Cells <- cells[grep(prefix[1], cells)]
  group2Cells <- cells[grep(prefix[2], cells)]
  group3Cells <- cells[grep(prefix[3], cells)]
  groups <- setNames(c(rep(prefix[1], length(group1Cells)), rep(prefix[2], length(group2Cells)),
                       rep(prefix[3], length(group3Cells))), c(group1Cells, group2Cells, group3Cells))
  #count <- count[, colnames(count) %in% genes]
  countByGroup <- conos:::collapseCellsByType(count, groups) %>% apply(2, function(r){
    r/sum(r)
  })
  countByGroup <- as.data.frame(t(countByGroup[, genes]))
  
  # classify genesets
  expDiff <- t(abs(apply(countByGroup, 1, diff)))
  geneClassCode <- apply(expDiff, 1, function(r){
    if(min(r)<0.1 & max(r)<0.1){
      "Conserved in human, mouse and marmoset"
    }else if(min(r) < 0.1){
      "Conserved in two species"
    } else{
      "Divergent gene"
    }
  })
  geneClassCode <- geneClassCode[match(rownames(countByGroup), names(geneClassCode))]
  
  countByGroup$classCode <- geneClassCode
  #ordering
  countByGroup <- countByGroup[order(countByGroup$classCode), ]
  rowAnnot <- data.frame(status=countByGroup$classCode)
  rownames(rowAnnot) <- rownames(countByGroup)
  
  pheatmap(countByGroup[, -4], cluster_rows=FALSE, cluster_cols=FALSE,
           color=colorRampPalette(c("white", "brown2"))(100), annotation_row=rowAnnot)
}

#' add summary pie chart at each node to show the conservation
heatmapPropCells <- function(d, col="blue"){
  visulize_ratio <- function(d){
    for(i in seq_len(length(d))){
      if(!is.leaf(d[[i]])){
        divergentRatio <- attr(d[[i]], "divergentRatio")
        conservationRatio <- (1-divergentRatio[1])/(1-div)
        d[[i]] <- calculate_divergence(d[[i]])
      }
    }
    return(d)
  }
  
  cbm <- function(dend){
    if(is.leaf(dend)) {
      nodeCells <- attr(dend, 'nodesinfo')
      targetCells <- nodeCells[nodeCells %in% cells]
      purity <- length(targetCells)/length(nodeCells)
      
      if(purity > pure.cutoff){
        attr(dend, "nodePar")$pch <- 19
        attr(dend, "nodePar")$col <- col
        attr(dend, "nodePar")$cex <- purity * 1.8
        #attr(dend, "nodePar")$cex <- 1
      }
      return(dend)
    } else{
      oa <- attributes(dend)
      dend <- lapply(dend, cbm)
      attributes(dend) <- oa
      return(dend)
    }
  }
  cbm(dend)
}

#' add summary pie chart at each node to show the conservative/divergent patterns
#' @param d dendrogram with conservation attributes
#' @import plotrix
#' @import RColorBrewer
#' @export
PlotTreeConservationPie <- function(d){
  require("RColorBrewer")
  colorpallete <- colorRampPalette(c("grey", "black"))(10)
  
  # plotting
  plot(d)
  conservation <- get_nodes_attr(d, "conservation")
  
  xy <- get_nodes_xy(d)
  lapply(1:length(conservation), function(r){
    if(!is.na(conservation[[r]])){
      # set conservation color scale
      conservationRatio <- round(conservation[[r]][1]/sum(conservation[[r]])*10, 0)
      plotrix::floating.pie(xy[r, 1], xy[r, 2], conservation[[r]],radius=0.5,
                            col=c("#ff0000","#80ff00","#00ffff"), border=colorpallete[conservationRatio])
    }
  })
}
