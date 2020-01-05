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
#' @return a list contains cells from every leaf
#' @export
getLeafcontentFlat <- function(cls.groups){
  leafcontent <- list()
  if(class(cls.groups) == "factor"){
    tmp <- as.character(cls.groups)
    names(tmp) <- names(cls.groups)
    cls.groups <- tmp
  }
  for(i in 1:length(unique(as.numeric(cls.groups)))){
    leafcontent[[i]] <- names(cls.groups[which(as.numeric(cls.groups)==i)])
  }
  return(leafcontent)
}

#' mapping leaf cells to cell annotation
transferLeaflabel <- function(leafcontent, cellannot){
  leafcontentAnnot <- lapply(leafcontent, function(x){cellannot[names(cellannot) %in% x]})
  return(leafcontentAnnot)
}

#' caculate jacard coefficiency from best maching tree
#' @param res either subsampled hierarchical graphs or a list contains merge matrix that have the similar format as igraph:::complete.dend
#' @param leaf.factor: leaf factor showing the assignment of cell ID in leaf nodes
#' @param clusters: cluster factor
bestClusterTreeThresholds <- function (res, leaf.factor, clusters, clmerges = NULL){
  clusters <- as.factor(clusters)
  #cl <- as.integer(as.character(clusters[names(leaf.factor)]))
  cl <- as.integer(clusters[names(leaf.factor)])
  clT <- tabulate(cl, nbins = length(levels(clusters)))
  mt <- table(cl, leaf.factor)
  if(class(res) != "list"){
    merges <- igraph:::complete.dend(res, FALSE)
  }

  if (is.null(clmerges)) {
    x <- conos:::treeJaccard(res$merges - 1L, as.matrix(mt),
                             clT)
    names(x$threshold) <- levels(clusters)
  }
  else {
    x <- conos:::treeJaccard(res$merges - 1L, as.matrix(mt),
                             clT, clmerges - 1L)
  }
  x
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
        percentage <- table(human.annot)[names(table(human.annot)) == maxProbLabel] / (totalCells[names(totalCells) == maxProbLabel] + 1)
        percentage <- round(percentage, 2)

        # calculate purity of branch
        percentage.Pop <- table(human.annot)[names(table(human.annot)) == maxProbLabel] / (sum(totalPopCells) + 1)
        percentage.Pop <- round(percentage.Pop, 2)

        if(percentage > 0.05 | length(grep("bi", names(r))) == 0){
          paste(maxProbLabel, "(", percentage.Pop, "/", percentage, ")", sep = "")
        } else{
          marmo.annot <- r[grep("bi", names(r))]
          totalPopCells <- table(marmo.annot)
          maxProbLabel <- names(table(marmo.annot))[which(table(marmo.annot) == max(table(marmo.annot)))][1]
          percentage <- table(marmo.annot)[names(table(marmo.annot)) == maxProbLabel] / (totalCells[names(totalCells) == maxProbLabel] + 1)
          percentage <- round(percentage, 2)

          # calculate purity of branch
          percentage.Pop <- table(marmo.annot)[names(table(marmo.annot)) == maxProbLabel] / (sum(totalPopCells) + 1)
          percentage.Pop <- round(percentage.Pop, 2)

          paste(maxProbLabel, "(", percentage.Pop, "/", percentage, ")", sep = "")
        }

      } else{
        #percentage <- round(table(r)/length(r), 2)
        totalPopCells <- table(r)
        maxProbLabel <- names(table(r))[which(table(r) == max(table(r)))][1]
        percentage <- table(r)[names(table(r)) == maxProbLabel] / (totalCells[names(totalCells) == maxProbLabel] + 1)
        percentage <- round(percentage, 2)

        # calculate purity of branch
        percentage.Pop <- table(r)[names(table(r)) == maxProbLabel] / (sum(totalPopCells) + 1)
        percentage.Pop <- round(percentage.Pop, 2)

        paste(maxProbLabel, "(", percentage.Pop, "/", percentage, ")", sep = "")
      }
    })
  } else{
    leaftransfer <- lapply(leaftransfer, function(r){
      #map human annotation
      maxProbLabel <- names(table(r))[which(table(r) == max(table(r)))][1]
      percentage <- table(r)[names(table(r)) == maxProbLabel] / (totalCells[names(totalCells) == maxProbLabel] + 1)
      percentage <- round(percentage, 3)

      # calculate purity of branch
      percentage.Pop <- table(r)[names(table(r)) == maxProbLabel] / (sum(table(r)) + 1)
      percentage.Pop <- round(percentage.Pop, 2)

      paste(maxProbLabel, "(", percentage.Pop, "/", percentage, ")", sep = "")
    })
  }

  leaftransfer <- unlist(lapply(leaftransfer, `[[`, 1))
  clusters <- as.numeric(d %>% labels())
  leaftransfer <- leaftransfer[clusters]
  return(leaftransfer)
}

#' assign attribute value to dendrogram
#' @param dend dendrogram obj
#' @param attribute attribute name
#' @param value vector contains attribute values
#' @return dendrogram obj with new attributes added to each node
assign_values_to_nodes <- function(dend, attribute, value){
  if (!is.dendrogram(dend)) stop("'dend' should be a dendrogram.")

  if (missing(value)) {
    warning("value is missing, returning the dendrogram as is.")
    return(dend)
  }

  nodes_length <- nnodes(dend) # length(labels(dend)) # it will be faster to use order.dendrogram than labels...
  if (nodes_length > length(value)) {
    if (warn) warning("Length of value vector was shorter than the number of nodes - vector value recycled")
    value <- rep(value, length.out = nodes_length)
  }

  set_value_to_node <- function(dend_node) {
    i_node_number <<- i_node_number + 1

    to_update_attr <- !is.infinite(value[i_node_number]) | is.na(value[i_node_number])
    if (to_update_attr) {
      # if(!is.infinite2(value[i_node_number])) {
      attr(dend_node, attribute) <- value[i_node_number]
    }

    if (length(attr(dend_node, attribute)) == 0) {
      attr(dend_node, attribute) <- NULL # remove attribute if it is empty
    }
    return(unclass(dend_node))
  }

  i_node_number <- 0
  new_dend <- dendrapply(dend, set_value_to_node)
  class(new_dend) <- "dendrogram"
  return(new_dend)
}

#' Change the format of tree
#' @param dend dendrogram object
#' @param renameCluster if rename the leafnode
#' @param cls.groups factor contains celltype annotation
#' @return list contains dendrogram obj, cells per leaf node and transfered groups
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
    cls.groups.new <- NULL
  }
  return(list(dendrogram=dend, leafcontent=leafcontent, new.groups=cls.groups.new))
}

#' subsampling graph
#' @param g igraph obj
#' @param stability.subsample # of subsampling
#' @param method clustering algorithm (e.g. walktrap, leiden etc)
#' @param saveGraph if save the susample result
#' @return original and subsampled graphs
#' @export
subSamplingGraph <- function(g, method=rleiden.detection, stability.subsamples=10,
                             stability.subsampling.fraction=0.95, saveGraph=T, prefix=NULL){

  cls <- rleiden.detection(g, K=3, resolution=c(0.7, 0.3, 0.3), min.community.size = 10)

  leafcontent <- getLeafcontent(cls)

  cls.mem <- membership(cls)
  cls.groups <- as.factor(cls.mem)
  cls.levs <- levels(cls.groups)

  sr <- NULL
  subset.clustering <- function(g, f=stability.subsampling.fraction, seed=NULL, cut=NULL){
    if(!is.null(seed)) { set.seed(seed) }
    vi <- sample(1:length(V(g)), ceiling(length(V(g))*(f)))
    sg <- induced_subgraph(g,vi)
    t <- method(sg,  K=4, resolution=c(0.7, 0.3, 0.3, 0.3), min.community.size = 10)
    if(is.null(cut)){
      return(t)
    } else{
      membership <- cut_at(t, cut)
      t$membership <- membership
      return(t)
    }
  }

  if(is.null(sr)){
    sr <- conos:::papply(1:stability.subsamples, function(i) subset.clustering(g, f=stability.subsampling.fraction,seed=i), n.cores=20)
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

#' generate subsampled dendrograms from gene expression matrix
#' @param counts scaled count matrix - rows are cells, colums are genes
#' @param cls.group original group to subsample
#' @param subsample.group list contains subsampled groups
#' @export
subSampleTree <- function(counts, cls.groups=NULL, subsample.groups=NULL, stability.subsamples=10, Variablegenes=NULL,
                             stability.subsampling.fraction=0.95, renameCluster=FALSE){
  if(is.null(cls.groups) & is.null(subsample.groups)){
    stop('need to specify either orignal groups or subsample groups')
  } else if(!is.null(cls.groups)){
    groups <- as.factor(cls.groups)
  }

  useVariablegenes = TRUE
  if(is.null(Variablegenes)){
    useVariablegenes = FALSE
  }

  sr <- NULL
  if(!is.null(cls.groups)){
    subsampling.cells <- function(counts, f=stability.subsampling.fraction, seed=NULL){
      if(!is.null(seed)) { set.seed(seed) }
      samples.loc <- sample(1:nrow(counts), ceiling(nrow(counts)*(f)))
      sampled.counts <- counts[samples.loc, ]
      sub.cls.groups <- cls.groups[names(cls.groups) %in% rownames(sampled.counts)]

      d <- cluster.matrix.expression.distances(sampled.counts, groups=cls.groups, dist="cor", useVariablegenes=useVariablegenes, variableGenes=Variablegenes,
                                               use.single.cell.comparisons=TRUE, use.scaled.data=TRUE)
      dendr <- hclust(as.dist(d), method='ward.D2')
      dend <- as.dendrogram(dendr)
      dendr <- TransferDend(dend, renameCluster=renameCluster, cls.groups = sub.cls.groups)
      if(renameCluster){
        sub.cls.groups <- dendr$new.groups
      }
      return(list(dend=dendr$dendrogram, clusters=sub.cls.groups))
    }

    if(is.null(sr)){
      sr <- parallel::mclapply(1:stability.subsamples, function(i) subsampling.cells(counts, f=stability.subsampling.fraction, seed=i), mc.cores=2)
    }
  } else{
    # have subsampled groups
    subsampling.cells <- function(counts, subsampled.groups, seed=NULL){
      sampled.counts <- counts[rownames(counts) %in% names(subsampled.groups), ]
      d <- cluster.matrix.expression.distances(sampled.counts, groups=subsampled.groups, dist="cor", useVariablegenes=useVariablegenes,
                                               variableGenes=Variablegenes, use.single.cell.comparisons=FALSE, use.scaled.data=FALSE)
      if(any(is.na(d))){
        d[is.na(d)] <- 0
      }
      dendr <- hclust(as.dist(d), method='ward.D2')
      dend <- as.dendrogram(dendr)
      dendr <- TransferDend(dend, renameCluster=renameCluster, cls.groups = subsampled.groups)
      if(renameCluster){
        subsampled.groups <- dendr$new.groups
      }
      return(list(dend=dendr$dendrogram, clusters=subsampled.groups))
    }
    if(is.null(sr)){
      sr <- parallel::mclapply(1:length(subsample.groups), function(i) subsampling.cells(counts, subsample.groups[[i]], seed=i), mc.cores=5)
    }
  }
  return(sr)
}

#' return stability scores for each node by searching for optimal matching tree
#' @param dataobj igraph object if algorithm=walktrap, conos or gene expression matrix if algorithm=expression
#' @param dend original dendrogram
#' @param algorithm algorithms for contruct hierachical tree in subsamples: walktrap or based on expression (expression)
#' @export
TreeStability <- function(dataobj, dend, algorithm="expression", cls.groups, cls.subsamples, stability.subsamples=10,
                          stability.subsampling.fraction=0.95, min.group.size=30, n.cores=10){

  supported.algorithms <- c("walktrap", "expression")
  if(!algorithm %in% supported.algorithms) {
    stop(paste0("only the following algorithm are currently supported: [",paste(supported.algorithms, collapse=' '),"]"))
  }

  if(algorithm == "walktrap"){
    if(!class(dataobj) == "igraph"){
      stop("Please provide an igraph object for walktrap algorithm")
    }
  } else if(algorithm == "expression"){
    if(!class(dataobj) == "Conos" & !class(dataobj) == "matrix"){
      stop("Please provide an conos object or matrix for expression algorithm")
    }
  }

  cls.groups <- as.factor(cls.groups)
  cls.levs <- levels(cls.groups)

  hc <- as.hclust(dend)
  clm <- hc$merge

  # transform to igraph format
  nleafs <- nrow(clm) + 1; clm[clm>0] <- clm[clm>0] + nleafs; clm[clm<=nleafs] <- -1*clm[clm<=nleafs]

  jc.hstats <- do.call(rbind, parallel::mclapply(cls.subsamples, function(st1){
    mf <- as.factor(membership(st1))
    if(algorithm == "walktrap"){
      st1g <- conos:::getClusterGraph(g, mf, plot=F, normalize=T)
      st1w <- walktrap.community(st1g, steps=8)
    } else if(algorithm == "expression"){

      # contruct tree structure
      ## for testing
      #source("~pkharchenko/m/pavan/DLI/conp2.r")
      d <- cluster.expression.distances(dataobj, groups=mf, dist='JS', min.cluster.size=1, min.samples=1, n.cores=1)
      subtreedend <- hclust(cluster::daisy(d), method='average')
      merges <- subtreedend$merge
      # transform to igraph format
      nleafs <- nrow(merges) + 1; merges[merges>0] <- merges[merges>0] + nleafs;
      merges[merges<=nleafs] <- -1*merges[merges<=nleafs]
      ###
      st1w = list(merges=merges)
    }
    x <- bestClusterTreeThresholds(st1w, mf, cls.groups, clm)
    x$threshold
  }, mc.cores=n.cores))

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

  return(list(dendrogram=dend, stability.loc=xy, stability.labels=round(x[to], 2)))
}

#' caculate stability only based on dendrogram and subsampled dendrograms
#' @param dend original dendrogram
#' @param cls.groups original cell clusters
#' @param subsamples list contains subsampled dendrograms and correspondent cell clusters
#' @import dendextend
#' @export
TreeStabilityDend <- function(dend, cls.groups, subsamples, assignValuestoNode=TRUE, n.cores=5){

  cls.groups <- as.factor(cls.groups)
  cls.levs <- levels(cls.groups)

  hc <- as.hclust(dend)
  clm <- hc$merge

  # getting into the same order
  cf <- factor(setNames(as.character(cls.groups), names(cls.groups)), levels=hc$labels)
  # transform to igraph format
  nleafs <- nrow(clm) + 1; clm[clm>0] <- clm[clm>0] + nleafs; clm[clm<=nleafs] <- -1*clm[clm<=nleafs]

  jc.hstats <- do.call(rbind, parallel::mclapply(subsamples, function(st1){
    mf <- as.factor(st1$clusters)
    # remove non-numeric annotation
    if(length(grep("singleton", mf)) > 0){
      mf <- mf[-1*grep("singleton", mf)]
    }

    subdendr <- TransferDend(st1$dend, renameCluster=FALSE, mf)
    subdend <- subdendr$dendrogram

    subtreedend <- as.hclust(subdend)
    # getting into the same order
    mf <- factor(setNames(as.character(mf), names(mf)), levels=subtreedend$labels)

    merges <- subtreedend$merge

    # transform to igraph format
    nleafs <- nrow(merges) + 1; merges[merges>0] <- merges[merges>0] + nleafs;
    merges[merges<=nleafs] <- -1*merges[merges<=nleafs]
    ###
    st1w = list(merges=merges)

    x <- bestClusterTreeThresholds(st1w, mf, cf, clm)
    x$threshold
  }, mc.cores=n.cores))

  if(length(grep("Error", jc.hstats)) > 0){
    loc <- unlist(lapply(1:nrow(jc.hstats), function(r){
      if(length(grep("Error", jc.hstats[r, ])) > 0){
        return(r)
      }
    }))
    jc.hstats <- t(apply(jc.hstats[-loc, ], 1, as.numeric))
  }

  stability <- list()
  stability$upper.tree <- clm
  stability$sr <- cls.groups
  stability$hierarchical <- list(jc=jc.hstats)

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

  # assign stability attribute to dendrogram
  if(assignValuestoNode){
    dend <- assign_values_to_nodes(dend, "stability", round(x[to], 2))
  }
  return(list(dendrogram=dend, stability.loc=xy, stability.labels=round(x[to], 2)))
}

#' Caculate flat stability score based subsampling clusters
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
#' @param d: dendrogram
#' @param cellannot: factor contains upperlevel annotation for each cell
#' @param leafContent: cells for each leaf node
#' @export
UpperLevelInfo <- function(d, cellannot, leafContent, propCutoff=0.15){

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

#' Get all upperlevel annotations concanated with "+"
#' @param d: dendrogram object
#' @param cellannot: factor contains upperlevel annotation for each cell
#' @param leafContent: leafcontent structure for each cluster
#' @export
AllUpperLevelInfo <- function(d, cellannot, leafContent, propCutoff=0.15){

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
#' @param dend dendrogram obj
#' @param cutoff purity cutoff
#' @export
getUpperLevelNode <- function(dend, cutoff=0.7){
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
    if(is.leaf(d)){
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

        facs <- facs/(prop + 1)
        normPercentage <- round((facs/sum(facs)), 2)
        attr(d, "normPercentage") <- normPercentage
        attr(d, "normFac") <- facs
      } else{
        if(length(cc) == 3){
          facs <- cc[2:3]
        } else{
          facs <- cc[2:4]
        }
        facs <- facs/(facTotal + 1)
        normPercentage <- round((facs/sum(facs)), 2)
        attr(d, "normPercentage") <- normPercentage
        attr(d, "normFac") <- facs
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
        facs <- facs/(prop + 1)
        normPercentage <- round((facs/sum(facs)), 2)
        attr(d, "normPercentage") <- normPercentage
        attr(d, "normFac") <- facs
      } else{
        if(length(cc) == 3){
          facs <- cc[2:3]
        } else{
          facs <- cc[2:4]
        }
        facs <- facs/(facTotal + 1)
        normPercentage <- round((facs/sum(facs)), 2)
        attr(d, "normPercentage") <- normPercentage
        attr(d, "normFac") <- facs
      }
      return(d)
    }
  }
  cbm(d, cellannot)
}

#' Normalize tree based on within-factors
#' @param dend dendrogram obj
#' @param fac factor contain all cells and their correspondent annotation
#' @return dendrogram with normalized values
#' @export
NormTreeWithinFactor <- function(d, fac, facprefix=c("human", "marmoset", "mouse")){
  celltable <- table(fac)
  cbm <- function(d){
    if(is.leaf(d)){
      # mapping cell to within-factor annotation table
      nodeinfo <- attr(d, "nodesinfo")
      nodeAnnot <- fac[names(fac) %in% nodeinfo]
      # get node cluster count distribution
      nodedistr <- table(nodeAnnot)
      # mapping node cells into total cell distr
      mappedAnnot <- celltable[match(names(nodedistr), names(celltable))]
      normalizedFac <- nodedistr/mappedAnnot
      prefix <- gsub(" .*", "", names(normalizedFac))
      withinFacCount <- unlist(lapply(facprefix, function(f){
        if(length(which(prefix == f)) == 0){
          return(0)
        } else{
          return(sum(normalizedFac[which(prefix==f)]))
        }
      }))
      names(withinFacCount) <- facprefix
      attr(d, "withinFacCount") <- withinFacCount
      attr(d, "withinFacPercent") <- round(withinFacCount/sum(withinFacCount), 2)
      return(d)
    } else{
      oa <- attributes(d)
      d <- lapply(d, cbm);
      attributes(d) <- oa
      nodeinfo <- attr(d, "nodesinfo")
      nodeAnnot <- fac[names(fac) %in% nodeinfo]
      # get node cluster count distribution
      nodedistr <- table(nodeAnnot)
      # mapping node cells into total cell distr
      mappedAnnot <- celltable[match(names(nodedistr), names(celltable))]
      normalizedFac <- nodedistr/mappedAnnot
      prefix <- gsub(" .*", "", names(normalizedFac))
      withinFacCount <- unlist(lapply(facprefix, function(f){
        if(length(which(prefix == f)) == 0){
          return(0)
        } else{
          return(sum(normalizedFac[which(prefix==f)]))
        }
      }))
      names(withinFacCount) <- facprefix
      attr(d, "withinFacCount") <- withinFacCount
      attr(d, "withinFacPercent") <- round(withinFacCount/sum(withinFacCount), 2)
    }
    return(d)
  }
  cbm(d)
}

#' caculate branch entropy attribute for normalized value
#' @param d dendrogram
#' @return add entropy attributes into the dendrogram
#' @import entropy
#' @export
TreeEntropy <- function(d, entropy.cutoff=2.0, withnFacNorm=FALSE){
  cbm <- function(d){
    if(is.leaf(d)){
      if(withnFacNorm){
        normpercentage <- attr(d, "withinFacPercent")
      } else{
        normpercentage <- attr(d, 'normPercentage')
      }

      if(is.null(normpercentage)){
        stop("please normalize tree first..........")
      }
      normentropy <- entropy::entropy(normpercentage, method='MM', unit='log2')
      attr(d, "entropy") <- normentropy
      if(normentropy > entropy.cutoff){
        attr(d, "nodePar")$pch <- 19
        attr(d, "nodePar")$col <- "blue"
      }
      return(d)
    } else{
      oa <- attributes(d)
      d <- lapply(d, cbm)
      attributes(d) <- oa
      if(withnFacNorm){
        normpercentage <- attr(d, "withinFacPercent")
      } else{
        normpercentage <- attr(d, 'normPercentage')
      }

      normentropy <- entropy::entropy(normpercentage, method='MM', unit='log2')
      attr(d, "entropy") <- normentropy
      if(normentropy > entropy.cutoff){
        attr(d, "nodePar")$pch <- 19
        attr(d, "nodePar")$col <- "blue"
      }
      return(d)
    }
  }
  cbm(d)
}

#' Mark nodes that mix well
#' @param d dendrogram
#' @return dendragram marked with nodes that have good mixing between factors
#' @import entropy
#' @export
MarkMixedNode <- function(d, mixing.cutoff=0.3, max.cutoff=0.6){
  cbm <- function(d){
    if(is.leaf(d)) {
      normpercentage <- attr(d, 'normPercentage')
      if(is.null(normpercentage)){
        stop("please normalize tree first..........")
      }

      if(normpercentage[1]>mixing.cutoff & normpercentage[1]<max.cutoff & !is.na(normpercentage[1])){
        attr(d, "nodePar")$pch <- 19
        attr(d, "nodePar")$col <- "blue"
      }
      return(d)
    } else{
      oa <- attributes(d)
      d <- lapply(d, cbm)
      attributes(d) <- oa
      normpercentage <- attr(d, 'normPercentage')
      if(normpercentage[1]>mixing.cutoff & normpercentage[1]<max.cutoff & !is.na(normpercentage[1])){
        attr(d, "nodePar")$pch <- 19
        attr(d, "nodePar")$col <- "blue"
      }
      return(d)
    }
  }
  cbm(d)
}

#' get mixed nodes by pairwise factor comparison
#' @param d: dendrogram ojb
#' @param mixing.cutoff: minimal cutoff value
#' @param max.cutoff: maximum cutoff value
#' @return dendrogram marked with well mixed nodes
MarkMixedNodePairwise <- function(d, mixing.cutoff=0.3, max.cutoff=0.6){

  nodeMixing <- function(d, p, mixing.cutoff=0.3, max.cutoff=0.6){
    if(p>mixing.cutoff & p<max.cutoff & !is.na(p)){
      return(TRUE)
    } else{
      return(FALSE)
    }
  }

  cbm <- function(d){
    if(is.leaf(d)) {
      normfac <- attr(d, 'normFac')
      if(is.null(normfac)){
        stop("please normalize tree first..........")
      }

      facname <- names(normfac)
      if(length(facname) > 3){
        stop("factor size > 3 is not supported")
      }

      if(length(facname) == 3){
        p1 <- normfac[1]/sum(normfac[1:2])
        p2 <- normfac[1]/sum(normfac[1:3])
        p3 <- normfac[2]/sum(normfac[2:3])

        # mixing across 3 species
        if(nodeMixing(d, p1) & nodeMixing(d, p2) & nodeMixing(d, p3)){
          attr(d, "nodePar")$pch <- 19
          attr(d, "nodePar")$col <- "bisque4"
        } else if(nodeMixing(d, p1)){
          attr(d, "nodePar")$pch <- 19
          attr(d, "nodePar")$col <- "blue"
        } else if(nodeMixing(d, p2)){
          attr(d, "nodePar")$pch <- 19
          attr(d, "nodePar")$col <- "green"
        } else if(nodeMixing(d, p3)){
          attr(d, "nodePar")$pch <- 19
          attr(d, "nodePar")$col <- "coral"
        }
      }

      return(d)
    } else{
      oa <- attributes(d)
      d <- lapply(d, cbm)
      attributes(d) <- oa
      normfac <- attr(d, 'normFac')
      if(is.null(normfac)){
        stop("please normalize tree first..........")
      }

      facname <- names(normfac)
      if(length(facname) > 3){
        stop("factor size > 3 is not supported")
      }

      if(length(facname) == 3){
        p1 <- normfac[1]/sum(normfac[1:2])
        p2 <- normfac[1]/sum(normfac[1:3])
        p3 <- normfac[2]/sum(normfac[2:3])

        # mixing across 3 species
        if(nodeMixing(d, p1) & nodeMixing(d, p2) & nodeMixing(d, p3)){
          attr(d, "nodePar")$pch <- 19
          attr(d, "nodePar")$col <- "bisque4"
        } else if(nodeMixing(d, p1)){
          attr(d, "nodePar")$pch <- 19
          attr(d, "nodePar")$col <- "blue"
        } else if(nodeMixing(d, p2)){
          attr(d, "nodePar")$pch <- 19
          attr(d, "nodePar")$col <- "green"
        } else if(nodeMixing(d, p3)){
          attr(d, "nodePar")$pch <- 19
          attr(d, "nodePar")$col <- "coral"
        }
      }
      return(d)
    }
  }
  cbm(d)
}

#' Set width of tree base on size of clusters
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

#' Color tree based on mixing of species - unnormalized
#' @param d dendrogram obj
#' @param fac factor contains species and barcodes information
#' @param colorpallete color pallete to color the tree
#'          e.g. colorRampPalette(c("blue", "grey", "grey", "grey", "red"))(101)
#' @return colored dendrogram
#' @export
dendSetColorByMixing <- function(d, fac, leafContent, normTofac=TRUE){
  fac <- as.factor(fac);
  if(length(levels(fac))>3) stop("factor with more than 3 levels are not supported")
  if(length(levels(fac))<2) stop("factor with less than 2 levels are not supported")

  totalCells <- table(fac)

  cc2col <- function(cc, rate=15, base=0.001){
    if(sum(cc)==0) {
      cc <- rep(1, length(cc))
    } else {
      # normalized by total number of cells
      if(length(cc) == 3){
        if(normTofac){
          cc <- round((cc[2:3]/totalCells)/sum(cc[2:3]/totalCells), 2)
        } else{
          cc <- round((cc[2:3])/sum(cc[2:3]), 2)
        }
        cv <- round(cc[1]*100, 0)
        colorpallete <- colorRampPalette(c("blue", "grey", "grey", "grey", "red"))(101)
        col <- colorpallete[cv + 1]

      } else{
        if(normTofac){
          cc <- round((cc[2:4]/totalCells)/sum(cc[2:4]/totalCells), 2)
        } else{
          cc <- round((cc[2:4])/sum(cc[2:4]), 2)
        }
        cv <- cc
        cv <- dexp(cv, rate)
        cv <- cv/rate * (1-base)
        col <- adjustcolor(rgb(cv[1],cv[2],cv[3], 1), offset = c(0.5, 0.5, 0.5, 0.1))
      }
    }
    return(col)
  }

  cbm <- function(d,fac) {
    if(is.leaf(d)) {
      #lc <- fac[leafContent[[as.numeric(attr(d,'label'))]]]
      cc <- attr(d, "cc")
      col <- cc2col(cc)
      attr(d,"edgePar") <- c(attr(d,"edgePar"),list(col=col))
      return(d);
    } else {
      oa <- attributes(d);
      d <- lapply(d,cbm,fac=fac);
      attributes(d) <- oa;
      cc <- attr(d, "cc")
      col <- cc2col(cc)
      attr(d,"edgePar") <- c(attr(d,"edgePar"),list(col=col))
      return(d);
    }
  }
  cbm(d,fac);
}

#' colored tree by species mixing based on normalized value - 2-species only
#' @param d dendrogram obj
#' @param colorpallete vector contains colors
#' @param facPrefix factor prefix selected to color - e.g. c("human", "mouse")
#' @param reNormalize re-caculate mixing
#' @param addValueToNode add re-caculated value to node
#' @import RColorBrewer
#' @export
dendSetColor2factorNormMixing <- function(d, reNormalize=TRUE, facPrefix, addValueToNode=TRUE, colorpallete){
  cc2col <- function(Normpercent, base=0.1){
    cv <- round(Normpercent*100, 0)
    return(colorpallete[cv + 1])
  }

  cbm <- function(d) {
    if(is.leaf(d)) {
      if(reNormalize){
        normfac <- attr(d, "normFac")
        if(is.null(normfac)){
          stop("please normalize the tree first......")
        } else{
          normfac <- normfac[match(facPrefix, names(normfac))]
          normpercent <- normfac/sum(normfac)
          if(addValueToNode){
            attr(d, "normPercentage") <- normpercent
          }
        }
      } else{
        normpercent <- attr(d, "normPercentage")
      }
      col <- cc2col(normpercent[1])
      attr(d,"edgePar") <- c(attr(d,"edgePar"),list(col=col))
      return(d);
    } else {
      oa <- attributes(d);
      d <- lapply(d,cbm);
      attributes(d) <- oa;
      if(reNormalize){
        normfac <- attr(d, "normFac")
        if(is.null(normfac)){
          stop("please normalize the tree first......")
        } else{
          normfac <- normfac[match(facPrefix, names(normfac))]
          normpercent <- normfac/sum(normfac)
          if(addValueToNode){
            attr(d, "normPercentage") <- normpercent
          }
        }
      } else{
        normpercent <- attr(d, "normPercentage")
      }
      col <- cc2col(normpercent[1])
      attr(d,"edgePar") <- c(attr(d,"edgePar"),list(col=col))
      return(d);
    }
  }
  cbm(d);
}

#' colored tree by species mixing based on normalized value
#' @param d dendrogram
#' @param withinFacNorm if colored the tree based on wihinfactor normalization
#' @return colored dendrogram
#' @export
dendSetColorByNormMixing <- function(d, withinFacNorm=FALSE){
  cc2col <- function(cc, rate=15, base=0.001) {
    if(length(cc)==2) { # 2-color
      cv <- cc
      cv <- dexp(c(cv[1], 0, cv[2]), rate) / base * (1-0.001)
      #cv <- c(cc[1],0,cc[2])+base; cv <- cv/max(cv) * (1-base)
      #rgb(base+cc[2],base,base+cc[3],1)
      adjustcolor(rgb(cv[1], cv[2], cv[3], 1), offset = c(0.5, 0.5, 0.5, 0.1))
    } else if(length(cc)==3) { # 3-color
      #cv <- 1 - cc
      cv <- cc
      cv <- dexp(cv, rate)
      #cv <- cv/max(cv) * (1-0.2)
      cv <- cv/rate * (1-base)
      adjustcolor(rgb(cv[1],cv[2],cv[3], 1), offset = c(0.5, 0.5, 0.5, 0.1))
    }
  }

  cbm <- function(d,fac){
    if(is.leaf(d)) {
      if(withinFacNorm){
        normpercent <- attr(d, "withinFacPercent")
      } else{
        normpercent <- attr(d, "normPercentage")
      }

      col <- cc2col(normpercent)
      attr(d,"edgePar") <- c(attr(d,"edgePar"),list(col=col))
      return(d);
    } else {
      oa <- attributes(d);
      d <- lapply(d,cbm);
      attributes(d) <- oa;
      if(withinFacNorm){
        normpercent <- attr(d, "withinFacPercent")
      } else{
        normpercent <- attr(d, "normPercentage")
      }

      col <- cc2col(normpercent)
      attr(d,"edgePar") <- c(attr(d,"edgePar"),list(col=col))
      return(d);
    }
  }
  cbm(d)
}

#' prune tree based on mixing
#' @param dend dendrogram obj
#' @param minCutoff minimal value for mixing cutoff
#' @param maxCutoff max value for mixing cutoff
#' @export
pruneTree <- function(dend, minCutoff=0.3, maxCutoff=0.4){

  is.branch.mixed <- function(dend, minCutoff=0.3, maxCutoff=0.4){
    is.mixed <- FALSE
    normPercentage <- attr(dend, "normPercentage")
    if(normPercentage[1] > minCutoff & normPercentage[1] < maxCutoff & !is.na(normPercentage[1])){
      is.mixed <- TRUE
    }
    return(is.mixed)
  }

  is.father.of.subtree.to.merge <- function(dend, minCutoff, maxCutoff) {
    # this function checks if the subtree we wish to merge is the direct child of the current branch (dend) we entered the function
    is.father <- FALSE
    for (i in seq_len(length(dend))){

      if(is.branch.mixed(dend[[i]], minCutoff, maxCutoff) == FALSE & !is.leaf(dend[[i]])) is.father <- TRUE
    }
    return(is.father)
  }

  search_mixed_subtree <- function(dend, minCutoff, maxCutoff){
    if (!is.father.of.subtree.to.merge(dend, minCutoff=minCutoff, maxCutoff=maxCutoff)){
      for (i in seq_len(length(dend))){
        if(is.leaf(dend[[i]])){
          next()
        }
        dend[[i]] <- search_mixed_subtree(dend[[i]], minCutoff, maxCutoff)
      }
    } else { # we'll merge
      subtree_loc <- 1
      if(!is.branch.mixed(dend[[subtree_loc]], minCutoff, maxCutoff)){
        # achieve attributes
        branch.attributes <- attributes(dend[[subtree_loc]])
        # get all the leaf labels in this subtree
        leaf.labels <- labels(dend[[subtree_loc]])
        # prune subtree
        if(length(leaf.labels) > 1){
          dend <- prune(dend, leaf.labels[-1])
          # assign attributes
          attr(dend[[subtree_loc]], "cc") <- branch.attributes$cc
          attr(dend[[subtree_loc]], "size") <- branch.attributes$size
          attr(dend[[subtree_loc]], "edgePar") <- branch.attributes$edgePar
          attr(dend[[subtree_loc]], "nodesinfo") <- branch.attributes$nodesinfo
          attr(dend[[subtree_loc]], "percentage") <- branch.attributes$percentage
          attr(dend[[subtree_loc]], "normPercentage") <- branch.attributes$normPercentage
        }
      }

      if(!is.branch.mixed(dend[[subtree_loc+1]], minCutoff, maxCutoff)){
        # achieve attributes
        branch.attributes <- attributes(dend[[subtree_loc+1]])
        # get all the leaf labels in this subtree
        leaf.labels <- labels(dend[[subtree_loc+1]])
        # prune subtree
        if(length(leaf.labels) > 1){
          dend <- prune(dend, leaf.labels[-1])
          # assign attributes
          attr(dend[[subtree_loc+1]], "cc") <- branch.attributes$cc
          attr(dend[[subtree_loc+1]], "size") <- branch.attributes$size
          attr(dend[[subtree_loc+1]], "edgePar") <- branch.attributes$edgePar
          attr(dend[[subtree_loc+1]], "nodesinfo") <- branch.attributes$nodesinfo
          attr(dend[[subtree_loc+1]], "percentage") <- branch.attributes$percentage
          attr(dend[[subtree_loc+1]], "normPercentage") <- branch.attributes$normPercentage
        }
      }
    }
    return(dend)
  }

  prune_mixed_subtree <- function(dend, minCutoff, maxCutoff){
    dend.label <- labels(dend)
    dend <- search_mixed_subtree(dend, minCutoff, maxCutoff)
    label.updated <- labels(dend)

    if(length(dend.label)==length(label.updated)){
      return(dend)
    } else{
      #dend.label <- label.updated
      dend <- prune_mixed_subtree(dend, minCutoff, maxCutoff)
      #label.updated <- labels(dend)
    }
  }

  new_dend <- prune_mixed_subtree(dend, minCutoff, maxCutoff)
  # re-assign space
  new_dend <- ladderize(new_dend, right=FALSE)
  new_dend <- ladderize(fix_members_attr.dendrogram(new_dend), right=FALSE)
  return(new_dend)
}

#' prune tree based on entropy
#' @param dend dendrogram obj
#' @param cutoff cutoff value
#' @export
pruneTreeEntropy <- function(dend, cutoff=2.9){

  is.branch.mixed <- function(dend, cutoff){
    is.mixed <- FALSE
    entropy <- attr(dend, "entropy")
    if(entropy > cutoff){
      is.mixed <- TRUE
    }
    return(is.mixed)
  }

  is.father.of.subtree.to.merge <- function(dend, cutoff) {
    # this function checks if the subtree we wish to merge is the direct child of the current branch (dend) we entered the function
    is.father <- FALSE
    for (i in seq_len(length(dend))){
      if(is.branch.mixed(dend[[i]], cutoff) == FALSE & !is.leaf(dend[[i]])) is.father <- TRUE
    }
    return(is.father)
  }

  search_mixed_subtree <- function(dend, cutoff){
    if (!is.father.of.subtree.to.merge(dend, cutoff)){
      for (i in seq_len(length(dend))){
        if(is.leaf(dend[[i]])){
          next()
        }
        dend[[i]] <- search_mixed_subtree(dend[[i]], cutoff)
      }
    } else { # we'll merge
      subtree_loc <- 1
      if(!is.branch.mixed(dend[[subtree_loc]], cutoff)){
        # achieve attributes
        branch.attributes <- attributes(dend[[subtree_loc]])
        # get all the leaf labels in this subtree
        leaf.labels <- labels(dend[[subtree_loc]])
        # prune subtree
        if(length(leaf.labels) > 1){
          dend <- prune(dend, leaf.labels[-1])
          # assign attributes
          attr(dend[[subtree_loc]], "cc") <- branch.attributes$cc
          attr(dend[[subtree_loc]], "size") <- branch.attributes$size
          attr(dend[[subtree_loc]], "edgePar") <- branch.attributes$edgePar
          attr(dend[[subtree_loc]], "nodesinfo") <- branch.attributes$nodesinfo
          attr(dend[[subtree_loc]], "percentage") <- branch.attributes$percentage
          attr(dend[[subtree_loc]], "normPercentage") <- branch.attributes$normPercentage
          attr(dend[[subtree_loc]], "stability") <- branch.attributes$stability
          attr(dend[[subtree_loc]], "entropy") <- branch.attributes$entropy
        }
      }

      if(!is.branch.mixed(dend[[subtree_loc+1]], cutoff)){
        # achieve attributes
        branch.attributes <- attributes(dend[[subtree_loc+1]])
        # get all the leaf labels in this subtree
        leaf.labels <- labels(dend[[subtree_loc+1]])
        # prune subtree
        if(length(leaf.labels) > 1){
          dend <- prune(dend, leaf.labels[-1])
          # assign attributes
          attr(dend[[subtree_loc+1]], "cc") <- branch.attributes$cc
          attr(dend[[subtree_loc+1]], "size") <- branch.attributes$size
          attr(dend[[subtree_loc+1]], "edgePar") <- branch.attributes$edgePar
          attr(dend[[subtree_loc+1]], "nodesinfo") <- branch.attributes$nodesinfo
          attr(dend[[subtree_loc+1]], "percentage") <- branch.attributes$percentage
          attr(dend[[subtree_loc+1]], "normPercentage") <- branch.attributes$normPercentage
          attr(dend[[subtree_loc+1]], "stability") <- branch.attributes$stability
          attr(dend[[subtree_loc+1]], "entropy") <- branch.attributes$entropy
        }
      }
    }
    return(dend)
  }

  prune_mixed_subtree <- function(dend, cutoff){
    dend.label <- labels(dend)
    dend <- search_mixed_subtree(dend, cutoff)
    label.updated <- labels(dend)

    if(length(dend.label)==length(label.updated)){
      return(dend)
    } else{
      #dend.label <- label.updated
      dend <- prune_mixed_subtree(dend, cutoff)
      #label.updated <- labels(dend)
    }
  }

  new_dend <- prune_mixed_subtree(dend, cutoff)
  # re-assign space
  new_dend <- ladderize(new_dend, right=FALSE)
  new_dend <- ladderize(fix_members_attr.dendrogram(new_dend), right=FALSE)
  return(new_dend)
}

#' remove leaves with low stability score
#' @param d dendrogram obj
#' @param cutoff stability cutoff to prune the tree
#' @return trimmed dendrogram
removeLowStabilityLeaf <- function(d, cutoff=0.5){
  if(is.null(attr(d, "stability"))){
    stop("Please measure stability first ...")
  }
  leaf.stability <- get_leaves_attr(d, "stability")
  leaf.labels <- get_leaves_attr(d, "label")
  leaf.info <- data.frame(labels=leaf.labels, stability=leaf.stability)
  leaf.info.filtered <- leaf.info[leaf.info$stability<=cutoff, ]$labels
  d <- prune(d, leaves=as.character(leaf.info.filtered))
  return(d)
}

#' Get robust homologous clusters from the tree
#'  extract homologous clusters at the brach node where species begin to differentiate
#' @param d pruned dendrogram
#' @return factor contains homologous clusters
getClusters <- function(d){

  is.father.of.leafnodes <- function(d){
    is.father <- FALSE
    for (i in seq_len(length(d))){
      if(is.leaf(d[[i]])) is.father <- TRUE
    }
    return(is.father)
  }

  search_tree <- function(d){
    if (!is.father.of.leafnodes(d)){
      for (i in seq_len(length(d))){
        d[[i]] <- search_tree(d[[i]])
      }
    } else { # get node information
      attr(d, "leaf") <- TRUE
    }
    #d <- ladderize(d, right=FALSE)
    return(d)
  }

  d.pruned <- search_tree(d)
  # get leaf clusters
  cells <- get_leaves_attr(d.pruned, "nodesinfo")
  sizes <- get_leaves_attr(d.pruned, "size")

  clusters <- unlist(lapply(1:length(sizes), function(r){
    rep(r, sizes[r])
  }))
  names(clusters) <- cells
  clusters <- as.factor(clusters)

  return(clusters)
}

#' build entire tree based on dendrogram and subsampled dendrograms automatically
#' @param dend dendrogram obj of original tree
#' @param subsamples list contains subsampled dendrograms and correspondent cell clusters
#' @param cls.groups clusters of original dendrograms
#' @param cellannot  cell annotations
#' @param species factor annoted which cell belong to which species
#' @param upperlevelannot higher-level cell annotations
#' @param plot plot the built tree or not
#' @export
buildSpeciesTree <- function(dend, subsamples=NULL, cls.groups, cellannot, species, upperlevelannot=NULL, renameCluster=TRUE, plot=TRUE){
  dendr <- TransferDend(dend, renameCluster=renameCluster, cls.groups = cls.groups)
  cls.groups <- dendr$new.groups

  dend <- dendr$dendrogram
  leafcontent <- dendr$leafcontent

  if(!is.null(subsamples)){
    stability.measurements <- TreeStabilityDend(dend, cls.groups, subsamples, n.cores=10)
  } else{
    stability.measurements = NULL
  }

  # add cluster attribute to dendrogram
  dend <- AddTreeAttribute(dend, species, leafcontent)

  dend <- dendSetWidthBysize(dend, scale=8)
  leaflabels <- mappingLabel(dend, leafcontent, cellannot, humanAnnot=T)

  # set labels
  dend <- set_labels(dend, paste(dend %>% labels(), leaflabels, sep=" "))

  # add upperlevel info to each nodes
  if(!is.null(upperlevelannot[1])){
    dend <- UpperLevelInfo(dend, cellannot=upperlevelannot, leafcontent, propCutoff = 0.1)
    upperLevelnodes <- getUpperLevelNode(dend, cutoff=0.65)

    # normalize Tree
    dend <- NormTree(dend, upperLevelnodes, upperlevelannot, species)
    dend <- dendSetColorByNormMixing(dend)
  } else{
    # need to modify
    colorpallete <- colorRampPalette(c("blue", "grey", "grey", "grey", "red"))(101)
    dend <- dendSetColorByMixing(dend, species, leafContent, colorpallete)
    upperLevelnodes = NULL
    ##
  }
  if(plot){
    if(!is.null(upperlevelannot) & !is.null(subsamples)){
      par(cex=1, mar=c(20, 5, 0, 10))
      plot(dend)
      text(upperLevelnodes$xy, labels=upperLevelnodes$upperlabel, adj=c(0.5, -1.2), cex=0.8, col="red")
      text(stability.measurements$stability.loc, labels=stability.measurements$stability.labels,
           adj=c(0.4, 0.1), cex=0.35, col="red")
    } else if(!is.null(upperlevelannot)){
      par(cex=1, mar=c(20, 5, 0, 10))
      plot(dend)
      text(upperLevelnodes$xy, labels=upperLevelnodes$upperlabel, adj=c(0.5, -1.2), cex=0.8, col="red")
    } else if(!is.null(subsamples)){
      par(cex=1, mar=c(20, 5, 0, 10))
      plot(dend)
      text(stability.measurements$stability.loc, labels=stability.measurements$stability.labels,
           adj=c(0.4, 0.1), cex=0.35, col="red")
    }
    else{
      plot(dend)
    }
  }
  return(list(dend=dend, stability=stability.measurements, upperLevelnodes=upperLevelnodes))
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

