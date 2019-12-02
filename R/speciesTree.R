# tree functions

#' get leaf information from igraph obj
getLeafcontent <- function(cls){
  leafcontent <- list()
  for(i in unique(cls$membership)){
    leafcontent[[i]] <- cls$names[which(cls$membership==i)]
  }
  return(leafcontent)
}

#' get leaf information from cell annotation
#' @param cls.groups factor contains celltype annotation for each cell
#' @ return a list contains cells from every leaf
#' @export
getLeafcontentFlat <- function(cls.groups){
  leafcontent <- list()
  for(i in 1:length(unique(as.numeric(cls.groups)))){
    leafcontent[[i]] <- names(cls.groups[which(as.numeric(cls.groups)==i)])
  }
  return(leafcontent)
}

#'
transferLeaflabel <- function(leafcontent, cellannot){
  leafcontentAnnot <- lapply(leafcontent, function(x){cellannot[names(cellannot) %in% x]})
  return(leafcontentAnnot)
}

#' mapping leaf node according to cell annotation
#' @param d dendrogram
#' @param cellannot factor contains celltype annotation for each cell
#' @return leaf annotation and purity measurement
#' @export
mappingLabel <- function(d, leafcontent, cellannot, humanAnnot=T){
  # mapping celltype labels onto the tree
  leaftransfer <- transferLeaflabel(leafcontent, cellannot)

  totalCells <- table(cellannot)
  # label with largest cell annotations
  if(humanAnnot == T){
    leaftransfer <- lapply(leaftransfer, function(r){
      #map human annotation
      if(length(grep("bi", names(r))) == 0){
        human.annot <- r
      } else{
        human.annot <- r[-1*grep("bi", names(r))]
      }

      totalPopCells <- table(human.annot)
      if(length(human.annot) > 0){
        #percentage <- round(table(human.annot)/length(human.annot), 2)
        maxProbLabel <- names(table(human.annot))[which(table(human.annot) == max(table(human.annot)))][1]
        percentage <- table(human.annot)[names(table(human.annot)) == maxProbLabel] / totalCells[names(totalCells) == maxProbLabel]
        percentage <- round(percentage, 2)

        # calculate purity of branch
        percentage.Pop <- table(human.annot)[names(table(human.annot)) == maxProbLabel] / sum(totalPopCells)
        percentage.Pop <- round(percentage.Pop, 2)

        if(percentage > 0.05 | length(grep("bi", names(r))) == 0){
          paste(maxProbLabel, "(", percentage.Pop, "/", percentage, ")", sep = "")
        } else{
          marmo.annot <- r[grep("bi", names(r))]
          totalPopCells <- table(marmo.annot)
          maxProbLabel <- names(table(marmo.annot))[which(table(marmo.annot) == max(table(marmo.annot)))][1]
          percentage <- table(marmo.annot)[names(table(marmo.annot)) == maxProbLabel] / totalCells[names(totalCells) == maxProbLabel]
          percentage <- round(percentage, 2)

          # calculate purity of branch
          percentage.Pop <- table(marmo.annot)[names(table(marmo.annot)) == maxProbLabel] / sum(totalPopCells)
          percentage.Pop <- round(percentage.Pop, 2)

          paste(maxProbLabel, "(", percentage.Pop, "/", percentage, ")", sep = "")
        }

      } else{
        #percentage <- round(table(r)/length(r), 2)
        maxProbLabel <- names(table(r))[which(table(r) == max(table(r)))][1]
        percentage <- table(r)[names(table(r)) == maxProbLabel] / totalCells[names(totalCells) == maxProbLabel]
        percentage <- round(percentage, 2)

        # calculate purity of branch
        percentage.Pop <- table(r)[names(table(r)) == maxProbLabel] / sum(totalPopCells)
        percentage.Pop <- round(percentage.Pop, 2)

        paste(maxProbLabel, "(", percentage.Pop, "/", percentage, ")", sep = "")
      }
    })
  } else{
    leaftransfer <- lapply(leaftransfer, function(r){
      #map human annotation
      maxProbLabel <- names(table(r))[which(table(r) == max(table(r)))][1]
      percentage <- table(r)[names(table(r)) == maxProbLabel] / totalCells[names(totalCells) == maxProbLabel]
      percentage <- round(percentage, 3)

      # calculate purity of branch
      percentage.Pop <- table(r)[names(table(r)) == maxProbLabel] / sum(table(r))
      percentage.Pop <- round(percentage.Pop, 2)

      paste(maxProbLabel, "(", percentage.Pop, "/", percentage, ")", sep = "")
    })
  }

  leaftransfer <- unlist(lapply(leaftransfer, `[[`, 1))
  clusters <- as.numeric(d %>% labels())
  leaftransfer <- leaftransfer[clusters]
  return(leaftransfer)
}

#' Change the format of tree
#' @param dend dendrogram object
#' @param renameCluster if rename the leafnode
#' @param cls.groups factor contains celltype annotation
#' @return list contains dendrogram obj, cells per leaf node and transfered table
#' @export
TransferDend <- function(dend, renameCluster=TRUE, cls.groups){

  cls.groups <- as.factor(cls.groups)
  cls.levs <- levels(cls.groups)

  require(dendextend)
  if(renameCluster){
    labels <- dend %>% labels()
    labels.new <- 1:length(labels)
    dend <- set_labels(dend, labels.new)

    label.table <- data.frame(old=as.character(labels), new=as.character(labels.new), stringsAsFactors = FALSE)

    # reassign label to groups
    cls.groups.new <- as.character(cls.groups)

    cls.groups.new <- label.table[match(cls.groups.new, label.table$old), ]$new
    #for(i in 1:nrow(label.table)){
    #  cls.groups.new[cls.groups.new == label.table[i, ]$old] <- label.table[i, ]$new
    #}
    names(cls.groups.new) <- names(cls.groups)

    leafcontent <- getLeafcontentFlat(cls.groups.new)
  } else {
    leafcontent <- getLeafcontentFlat(cls.groups)
    label.table <- NULL
  }
  return(list(dendrogram=dend, leafcontent=leafcontent, label.table=label.table))
}

#' subsampling graph
#' @param g igraph obj
#' @param stability.subsample # of subsampling
#' @param method clustering algorithm (e.g. walktrap, leiden etc)
#' @param saveGraph if save the susample result
#'
#' @export
subSamplingGraph <- function(g, method=rleiden.detection, stability.subsamples=10,
                             stability.subsampling.fraction=0.95, saveGraph=T, prefix=NULL){

  cls <- rleiden.detection(g, K=3, resolution=c(0.5, 0.3, 0.3), min.community.size = 10)

  leafcontent <- getLeafcontent(cls)

  cls.mem <- membership(cls)
  cls.groups <- as.factor(cls.mem)
  cls.levs <- levels(cls.groups)

  sr <- NULL
  subset.clustering <- function(g, f=stability.subsampling.fraction, seed=NULL, cut=NULL){
    if(!is.null(seed)) { set.seed(seed) }
    vi <- sample(1:length(V(g)), ceiling(length(V(g))*(f)))
    sg <- induced_subgraph(g,vi)
    t <- method(sg,  K=3, resolution=c(0.5, 0.3, 0.3), min.community.size = 10)
    if(is.null(cut)){
      return(t)
    } else{
      membership <- cut_at(t, cut)
      t$membership <- membership
      return(t)
    }
  }

  if(is.null(sr)){
    sr <- conos:::papply(1:stability.subsamples, function(i) subset.clustering(g,f=stability.subsampling.fraction,seed=i), n.cores=20)
  }

  # save subsampled graph
  if(saveGraph == T){
    if(is.null(prefix)){
      saveRDS(list(mem=cls, subsample.mem=sr), "graph.membership.rds")
    } else{
      saveRDS(list(mem=cls, subsample.mem=sr), paste(prefix, "graph.membership.rds", sep="."))
    }
  }
  return(list(mem=cls, subsample.mem=sr))
}

#' return stability score based on walktrap in hierachical structure
#' @export
TreeStability <- function(g, dend, cls.groups, cls.subsamples, stability.subsamples=10,
                          stability.subsampling.fraction=0.95, min.group.size=30, n.cores=10){

  cls.groups <- as.factor(cls.groups)
  cls.levs <- levels(cls.groups)

  hc <- as.hclust(dend)
  clm <- hc$merge

  # transform to igraph format
  nleafs <- nrow(clm) + 1; clm[clm>0] <- clm[clm>0] + nleafs; clm[clm<=nleafs] <- -1*clm[clm<=nleafs]

  jc.hstats <- do.call(rbind, mclapply(cls.subsamples, function(st1) {
    mf <- as.factor(st1)
    st1g <- conos:::getClusterGraph(g, mf, plot=F, normalize=T)
    st1w <- walktrap.community(st1g, steps=8)

    x <- conos:::bestClusterTreeThresholds(st1w, mf, cls.groups, clm)
    x$threshold
  }, mc.cores=n.cores))

  jc.hstats <- apply(jc.hstats, 1, function(x){
    if(length(grep("Error", x)) == 0){
      as.numeric(x)
    }
  })
  jc.hstats <- do.call(rbind, jc.hstats)

  stability <- list()
  stability$upper.tree <- clm
  stability$sr <- cls.groups
  stability$hierarchical <- list(jc=jc.hstats);
  #if(verbose) cat("done\n");

  require(dendextend)
  # depth-first traversal of a merge matrix
  t.dfirst <- function(m,i=nrow(m)) {
    rl <- m[i,1]; if(rl<0) { rl <- abs(rl) } else { rl <- t.dfirst(m,rl) }
    rr <- m[i,2]; if(rr<0) { rr <- abs(rr) } else { rr <- t.dfirst(m,rr) }
    c(i+nrow(m)+1,rl,rr)
  }
  xy <- get_nodes_xy(dend)
  to <- t.dfirst(hc$merge)
  x <- apply(jc.hstats, 2, median)

  return(list(dendrogram=dend, stability=xy, labels=round(x[to], 2), leafcontent=leafcontent))
}

# Caculate flat stability score based subsampling clusters
#' @export
TreeStabilityFlat <- function(cls.groups, cls.subsamples, n.cores=10){

  jc.stats <- do.call(rbind,conos:::papply(cls.subsamples, function(o) {
    p1 <- membership(o);
    p2 <- cls.groups[names(p1)]; p1 <- as.character(p1)
    x <- tapply(1:length(p2),p2,function(i1) {
      i2 <- which(p1==p1[i1[[1]]])
      length(intersect(i1, i2))/length(unique(c(i1,i2)))
    })
  },n.cores=n.cores,mc.preschedule=T))

  # Adjusted rand index
  ari <- unlist(conos:::papply(sr,function(o) { ol <- membership(o); adjustedRand(as.integer(ol),
                                                                                  as.integer(cls.groups[names(ol)]),randMethod='HA') },n.cores=n.cores))
  stability <- list(flat=list(jc=jc.stats,ari=ari))
  return(stability)
}

#' Add attribute into the tree
#' @param d dendrogram obj
#' @param fac factor contains cell information and its correspondent species
#' @param leafContent cells in leaf nodes
#' @return dendrogram with related attributes added into the tree
#' @export
AddTreeAttribute <- function(d, fac, leafContent){
    fac <- as.factor(fac);
    if(length(levels(fac))>3) stop("factor with more than 3 levels are not supported")
    if(length(levels(fac))<2) stop("factor with less than 2 levels are not supported")

    totalCells <- table(fac)

    cbm <- function(d,fac) {
      if(is.leaf(d)) {
        lc <- fac[leafContent[[as.numeric(attr(d,'label'))]]]
        cc <- c(sum(is.na(lc)),table(lc))
        attr(d,'cc') <- cc
        size <- length(lc)
        attr(d, 'size') <- size
        return(d)
      } else {
        oa <- attributes(d)
        d <- lapply(d,cbm,fac=fac);
        attributes(d) <- oa
        cc <- attr(d[[1]],'cc')+attr(d[[2]],'cc')
        attr(d,'cc') <- cc
        size <- attr(d[[1]],'size')+attr(d[[2]],'size')
        attr(d, 'size') <- size
        return(d)
      }
    }
    cbm(d, fac)
}

#' add upperlevel attributes into the tree
#' @export
UpperLevelInfo <- function(d, cellannot, leafContent, propCutoff=0.15){
  # d: dendrogram
  # cellannot: factor contains upperlevel annotation for each cell
  # leafContent: leafcontent structure for each cluster

  cbm <- function(d, cellannot){
    if(is.leaf(d)) {
      label <- attr(d, 'label')
      label <- as.numeric(gsub(" .*", "", label))
      leafs <- leafContent[[label]]
      attr(d,"nodesinfo") <- leafs
      lc <- cellannot[leafs]

      Cellpercentage <- table(lc)/sum(table(lc))

      Upperlabel <- names(Cellpercentage[which(Cellpercentage == max(Cellpercentage))])[1]
      attr(d, "upperlabel") <- Upperlabel
      attr(d, "percentage") <- Cellpercentage[which(Cellpercentage == max(Cellpercentage))][1]
      return(d)
    } else {
      oa <- attributes(d)
      d <- lapply(d, cbm, cellannot=cellannot);
      attributes(d) <- oa
      leafs <- c(attr(d[[1]],'nodesinfo'), attr(d[[2]],'nodesinfo'))
      attr(d,"nodesinfo") <- leafs
      lc <- cellannot[leafs]

      Cellpercentage <- table(lc)/sum(table(lc))

      Upperlabel <- names(Cellpercentage[which(Cellpercentage == max(Cellpercentage))])[1]
      attr(d, "upperlabel") <- Upperlabel
      attr(d, "percentage") <- Cellpercentage[which(Cellpercentage == max(Cellpercentage))][1]
      return(d)
    }
  }
  cbm(d, cellannot)
}

#' Get all upperlevel annotations
#' @export
AllUpperLevelInfo <- function(d, cellannot, leafContent, propCutoff=0.15){
  # Add all upper-level annotations in each branch concated with "+"
  # d: dendrogram
  # cellannot: factor contains upperlevel annotation for each cell
  # leafContent: leafcontent structure for each cluster

  cbm <- function(d, cellannot){
    if(is.leaf(d)) {
      label <- attr(d, 'label')
      label <- as.numeric(gsub(" .*", "", label))
      leafs <- leafContent[[label]]
      attr(d,"nodesinfo") <- leafs
      lc <- cellannot[leafs]

      Cellpercentage <- table(lc)/sum(table(lc))
      Cellpercentage <- Cellpercentage[Cellpercentage > propCutoff]
      Upperlabel <- paste(names(table(lc))[names(table(lc)) %in% names(Cellpercentage)], collapse = "+")
      attr(d, "upperlabelAll") <- Upperlabel
      return(d)
    } else {
      oa <- attributes(d)
      d <- lapply(d, cbm, cellannot=cellannot);
      attributes(d) <- oa
      leafs <- c(attr(d[[1]],'nodesinfo'), attr(d[[2]],'nodesinfo'))
      attr(d, "nodesinfo") <- leafs
      lc <- cellannot[leafs]

      Cellpercentage <- table(lc)/sum(table(lc))
      Cellpercentage <- Cellpercentage[Cellpercentage > propCutoff]
      Upperlabel <- paste(names(table(lc))[names(table(lc)) %in% names(Cellpercentage)], collapse = "+")
      attr(d, "upperlabelAll") <- Upperlabel
      return(d)
    }
  }
  cbm(d, cellannot)
}


#' get highest node with upperlevel information
#' @export
getUpperLevelNode <- function(d, cutoff=0.7){
  # d: dendrogram with upperlevel and percentage features

  upperlabel <- dend %>% get_nodes_attr("upperlabel")
  percentage <- round(dend %>% get_nodes_attr("percentage"), 2)
  xy <- get_nodes_xy(dend)
  leaf <- dend %>% get_nodes_attr("leaf")
  height <- dend %>% get_nodes_attr("height")

  LabelInfo <- data.frame(height=height, upperlabel=upperlabel, percentage=percentage, leaf=leaf)
  LabelInfo <- cbind(xy, LabelInfo)

  # cutoff leaves
  LabelInfo <- LabelInfo[-1*which(LabelInfo$leaf==TRUE), ]
  LabelInfo <- LabelInfo[LabelInfo$percentage > cutoff, ]

  uppernodes <- lapply(unique(LabelInfo$upperlabel), function(r){
    nodes <- LabelInfo[LabelInfo$upperlabel == r, ]
    nodes[which(nodes$height == max(nodes$height)), ]
  })
  uppernodes <- dplyr::bind_rows(uppernodes)
  return(list(upperlabel=uppernodes$upperlabel, xy=uppernodes[, 1:2], percentage=uppernodes$percentage, height=uppernodes$height))
}

#' normalize to upperlevel
#'
#' @param d: dendrogram
#' @param cellannot: factor contains upperlevel annotation for each cell
#' @param leafContent: leafcontent structure for each cluster
#' @return dendrogram with normalized attributes
#' @export
NormTree <- function(d, upperLevelnodes, cellannot, fac, cutoff=0.7){
  fac <- as.factor(fac)
  if(length(levels(fac))>3) stop("factor with more than 3 levels are not supported")
  if(length(levels(fac))<2) stop("factor with less than 2 levels are not supported")

  Cellinfo <- data.frame(barcodes=names(cellannot), celltype=cellannot)
  Cellinfo <- merge(Cellinfo, data.frame(barcodes=names(fac), species=fac), by=c("barcodes"))
  Cellprop <- table(Cellinfo$celltype, Cellinfo$species)
  facTotal <- table(fac)

  upperLevelnode <- data.frame(upperLevelnodes$xy, upperLabel=upperLevelnodes$upperlabel, height=upperLevelnodes$height)

  cbm <- function(d, cellannot){
    if(is.leaf(d)) {
      upperlabel <- attr(d, 'upperlabel')
      height <- attr(d, 'height')
      percentage <- attr(d, 'percentage')
      cc <- attr(d, 'cc')

      if(percentage > cutoff){
        prop <- Cellprop[rownames(Cellprop) == upperlabel]
        if(length(cc) == 3){
          facs <- cc[2:3]
        } else{
          facs <- cc[2:4]
        }

        facs <- facs/prop
        normPercentage <- round((facs/sum(facs))[1], 2)
        attr(d, "normPercentage") <- normPercentage
      } else{
        if(length(cc) == 3){
          facs <- cc[2:3]
        } else{
          facs <- cc[2:4]
        }
        facs <- facs/facTotal
        normPercentage <- round((facs/sum(facs))[1], 2)
        attr(d, "normPercentage") <- normPercentage
      }
      return(d)
    } else {
      oa <- attributes(d)
      d <- lapply(d, cbm, cellannot=cellannot);
      attributes(d) <- oa
      percentage <- attr(d, 'percentage')
      upperlabel <- attr(d, 'upperlabel')
      cc <- attr(d, 'cc')
      if(percentage > cutoff){
        prop <- Cellprop[rownames(Cellprop) == upperlabel]
        if(length(cc) == 3){
          facs <- cc[2:3]
        } else{
          facs <- cc[2:4]
        }
        facs <- facs/prop
        normPercentage <- round((facs/sum(facs))[1], 2)
        attr(d, "normPercentage") <- normPercentage
      } else{
        if(length(cc) == 3){
          facs <- cc[2:3]
        } else{
          facs <- cc[2:4]
        }
        facs <- facs/facTotal
        normPercentage <- round((facs/sum(facs))[1], 2)
        attr(d, "normPercentage") <- normPercentage
      }
      return(d)
    }
  }
  cbm(d, cellannot)
}

#' Set width of tree base on size of clusters
#'
#' @param d dendrogram obj
#' @param scale scales to set branch size
#' @return dendrogram with the width of branches indicates cluster size
#' @export
dendSetWidthBysize <- function(d, scale=10){
  cbm <- function(d, fac) {
    if(is.leaf(d)){
      size <- attr(d,'size')
      lwd <- log(size + 1)/log(scale)
      attr(d,"edgePar") <- c(attr(d,"edgePar"),list(lwd=lwd))
      return(d)
    } else {
      oa <- attributes(d)
      d <- lapply(d,cbm)
      attributes(d) <- oa
      size <- attr(d,'size')
      lwd <- log(size + 1)/log(scale)
      attr(d,"edgePar") <- c(attr(d,"edgePar"),list(lwd=lwd))
      return(d)
    }
  }
  cbm(d,fac)
}

#' Color tree based on mixing of species
#'
#' @param d dendrogram obj
#' @param fac factor contains species and barcodes information
#' @param colorpallete color pallete to color the tree
#'          e.g. colorRampPalette(c("blue", "grey", "grey", "grey", "red"))(101)
#' @return colored dendrogram
#' @export
dendSetColorByMixing <- function(d,fac,leafContent, colorpallete){
  fac <- as.factor(fac);
  if(length(levels(fac))>3) stop("factor with more than 3 levels are not supported")
  if(length(levels(fac))<2) stop("factor with less than 2 levels are not supported")

  totalCells <- table(fac)

  cc2col <- function(cc,base=0.1){
    if(sum(cc)==0) {
      cc <- rep(1, length(cc))
    } else {
      # normalized by total number of cells
      if(length(cc) == 3){
        cc <- round((cc[2:3]/totalCells)/sum(cc[2:3]/totalCells), 2)
      } else{
        cc <- round((cc[2:4]/totalCells)/sum(cc[2:4]/totalCells), 2)
      }

    }

    cv <- round(cc[1]*100, 0)
    return(colorpallete[cv + 1])
  }

  cbm <- function(d,fac) {
    if(is.leaf(d)) {
      lc <- fac[leafContent[[as.numeric(attr(d,'label'))]]]
      cc <- c(sum(is.na(lc)),table(lc));
      col <- cc2col(cc)
      attr(d,"edgePar") <- c(attr(d,"edgePar"),list(col=col))
      attr(d,'cc') <- cc;
      return(d);
    } else {
      oa <- attributes(d);
      d <- lapply(d,cbm,fac=fac);
      attributes(d) <- oa;
      cc <- attr(d[[1]],'cc')+attr(d[[2]],'cc')
      col <- cc2col(cc)
      attr(d,"edgePar") <- c(attr(d,"edgePar"),list(col=col))
      attr(d,'cc') <- cc;
      return(d);
    }
  }
  cbm(d,fac);
}

#' colored tree by species mixing based on normalized value
#' @export
dendSetColorByNormMixing <- function(d, colorpallete){
  cc2col <- function(Normpercent, base=0.1){
    cv <- round(Normpercent*100, 0)
    return(colorpallete[cv + 1])
  }

  cbm <- function(d,fac) {
    if(is.leaf(d)) {
      normpercent <- attr(d, "normPercentage")
      col <- cc2col(normpercent)
      attr(d,"edgePar") <- c(attr(d,"edgePar"),list(col=col))
      return(d);
    } else {
      oa <- attributes(d);
      d <- lapply(d,cbm);
      attributes(d) <- oa;
      normpercent <- attr(d, "normPercentage")
      col <- cc2col(normpercent)
      attr(d,"edgePar") <- c(attr(d,"edgePar"),list(col=col))
      return(d);
    }
  }
  cbm(d,colorpallete);
}

#' recursive prunning tree
#' @export
prunningTree <- function(d, fac, leafContent, minSize=20){
  cbm <- function(d,fac){
    if(is.leaf(d)) {
      label <- attr(d,'label')
      if(length(grep(" ", label))>0){
        # if tranfer leaf label before
        label <- unlist(strsplit(label, " "))[1]
      }
      lc <- fac[leafContent[[as.numeric(label)]]]
      cc <- length(lc)
      if(cc > minSize){
        return(d);
      }
    } else {
      oa <- attributes(d);
      d <- lapply(d,cbm,fac=fac);
      if(is.null(d[[1]])){
        d <- d[[2]]
        oa <- attributes(d)
        attributes(d) <- oa
      } else if(is.null(d[[2]])){
        d <- d[[1]]
        oa <- attributes(d)
        attributes(d) <- oa
      } else{
        attributes(d) <- oa
      }
      return(d);
    }
  }
  cbm(d,fac);
}

#' cross-species mapping summary: 1:1, or many:many
#' @export
clusterMappingSummary <- function(leafannot, prefix="bi", speciesNames=c("marmo", "human"), cutoff=0.01){
  # leafannot: list contains celltype annotation for each leafnode
  # prefix: prefix in cells that distinguish two species
  # cutoff: cutoff to filter small percentage cells
  # returns a summary of mapping statistics
  MappingStat <- lapply(leafannot, function(r){
    r[grep(prefix, names(r))] <-  paste(prefix, r[grep(prefix, names(r))])
    summary <- table(r)
    # filtering small proportion cells
    cellProp <- summary/sum(summary)
    summary <- summary[names(summary) %in% names(cellProp[cellProp > cutoff])]

    celltypeCount <- c(length(grep(prefix, names(summary))), length(summary) - length(grep(prefix, names(summary))))
    names(celltypeCount) <- speciesNames
    return(celltypeCount)
  })
  MappingStat <- do.call(rbind, MappingStat)
  return(MappingStat)
}

#' Run recursive leiden clustering algorithm - modified from conos (https://github.com/hms-dbmi/conos)
#'
#' @param graph igraph object
#' @export

rleiden.detection <- function(graph, K=2, renameCluter=TRUE, n.cores=parallel::detectCores(logical=F),
                              min.community.size=10, verbose=FALSE, resolution=1, K.current=1, hierarchical=FALSE, ...){

  if(verbose & K.current==1) cat(paste0("running ",K,"-recursive Leiden clustering: "));
  if(length(resolution)>1) {
    if(length(resolution)!=K) { stop("resolution value must be either a single number or a vector of length K")}
    res <- resolution[K.current]
  } else { res <- resolution }
  mt <- leiden.community(graph, resolution=res, ...);

  mem <- membership(mt);
  tx <- table(mem)
  ivn <- names(tx)[tx<min.community.size]
  if(length(ivn)>1) {
    mem[mem %in% ivn] <- as.integer(ivn[1]); # collapse into one group
  }
  if(verbose) cat(length(unique(mem)),' ');

  if(K.current<K) {
    # start recursive run
    if(n.cores>1) {
      wtl <- mclapply(conos:::sn(unique(mem)), function(cluster) {
        cn <- names(mem)[which(mem==cluster)]
        sg <- induced.subgraph(graph,cn)
        rleiden.detection(induced.subgraph(graph,cn), K=K, resolution=resolution, K.current=K.current+1,
                          min.community.size=min.community.size, hierarchical=hierarchical, verbose=verbose, n.cores=n.cores)
      },mc.cores=n.cores,mc.allow.recursive = FALSE)
    } else {
      wtl <- lapply(conos:::sn(unique(mem)), function(cluster) {
        cn <- names(mem)[which(mem==cluster)]
        sg <- induced.subgraph(graph,cn)
        rleiden.detection(induced.subgraph(graph,cn), K=K, resolution=resolution, K.current=K.current+1,
                          min.community.size=min.community.size, hierarchical=hierarchical, verbose=verbose, n.cores=n.cores)
      })
    }
    # merge clusters, cleanup
    mbl <- lapply(wtl,membership);
    # combined clustering factor
    fv <- unlist(lapply(names(wtl),function(cn) {
      paste(cn,as.character(mbl[[cn]]),sep='-')
    }))
    names(fv) <- unlist(lapply(mbl,names))
  } else {
    fv <- mem;
    if(hierarchical) {
      # use walktrap on the last level
      wtl <- conos:::papply(conos:::sn(unique(mem)), function(cluster) {
        cn <- names(mem)[which(mem==cluster)]
        sg <- induced.subgraph(graph,cn)
        res <- walktrap.community(induced.subgraph(graph,cn))
        res$merges <- igraph:::complete.dend(res,FALSE)
        res
      },n.cores=n.cores)
    }
  }

  if(K.current==1) {
    if(verbose) {
      cat(paste0(' detected a total of ',length(unique(fv)),' clusters '));
      cat("done\n");
    }
    # rename clusters
    if(renameCluter){
      tmp <- fv
      label.table <- data.frame(old=as.character(unique(fv)), new=as.character(seq(1, length(unique(fv)), 1)),
                                stringsAsFactors = FALSE)

      # reassign label to groups
      fv.new <- as.character(fv)

      fv.new <- as.character(match(fv.new, label.table$old))
      fv <- fv.new
      names(fv) <- names(tmp)
    }

  }

  # enclose in a masquerading class
  combd <- NULL
  res <- list(membership=fv,dendrogram=combd,algorithm='rleiden');
  class(res) <- rev("fakeCommunities")
  return(res)
}

