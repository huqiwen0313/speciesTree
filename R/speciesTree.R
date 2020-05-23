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
bestClusterTreeThresholds <- function(res, leaf.factor, clusters, clmerges = NULL){
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
#' @param purity TRUE/FLASE if TRUE, return purity measurements
#' @return leaf annotation and purity measurement (if purity=TRUE)
#' @export
mappingLabel <- function(d, leafcontent, cellannot, humanAnnot=FALSE, purity=FALSE){
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

          if(purity){
            paste(maxProbLabel, "(", percentage.Pop, "/", percentage, ")", sep = "")
          } else{
            maxProbLabel
          }
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

        if(purity){
          paste(maxProbLabel, "(", percentage.Pop, "/", percentage, ")", sep = "")
        } else{
          maxProbLabel
        }
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

      if(purity){
        paste(maxProbLabel, "(", percentage.Pop, "/", percentage, ")", sep = "")
      } else{
        maxProbLabel
      }
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
#' @export
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

  #cls.groups <- as.factor(cls.groups)
  #cls.levs <- levels(cls.groups)

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
                             stability.subsampling.fraction=0.95, renameCluster=FALSE, use.scaled.data=TRUE){
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
                                               use.single.cell.comparisons=FALSE, use.scaled.data=use.scaled.data)
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
                                               variableGenes=Variablegenes, use.single.cell.comparisons=FALSE, use.scaled.data=use.scaled.data)
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
      #sr <- lapply(1:length(subsample.groups), function(i) subsampling.cells(counts, subsample.groups[[i]], seed=i))
    }
  }
  return(sr)
}

#' Generate subsampled tree from conos joint graph based on distance on the joint graph
#' @param graph conos joint graph
#' @param subsample.group list contains subsampled groups
#' @import sccore
subsampleTree.conos <- function(graph, subsample.groups=NULL, stability.subsamples=10,
                                stability.subsampling.fraction=0.95, renameCluster=FALSE){
  if(is.null(subsample.groups)){
    stop('please provide subsampled groups')
  }

  # based on conos graph function
  computeGraphDistance <- function(graph){
    conn.comps <- igraph::components(graph)
    #min.visited.verts=1000
    #min.visited.verts = min(min.visited.verts, min(conn.comps$csize) - 1)
    min.visited.verts=1000
    min.visited.verts = min(min.visited.verts, min(conn.comps$csize) - 1)
    max.hitting.nn.num=0
    max.commute.nn.num=0
    min.prob.lower=1e-7
    min.prob=1e-3

    if(max.hitting.nn.num == 0) {
      max.hitting.nn.num <- length(igraph::V(graph)) - 1
    }
    adj.info <- graphToAdjList(graph)

    commute.times <- get_nearest_neighbors(adj.info$idx, adj.info$probabilities, min_prob=min.prob,
                                           min_visited_verts=min.visited.verts, n_cores=1, max_hitting_nn_num=max.hitting.nn.num,
                                           max_commute_nn_num=max.commute.nn.num, min_prob_lower=min.prob.lower, verbose=TRUE)
    n.neighbors <- length(V(graph))
    ct.top <- sapply(commute.times$dist, `[`, 1:(n.neighbors-1)) %>% t() + 1
    ct.top.ids <- sapply(commute.times$idx, `[`, 1:(n.neighbors-1)) %>% t() + 1

    ct.top.ids <- cbind(1:nrow(ct.top.ids), ct.top.ids)
    ct.top <- cbind(rep(0, nrow(ct.top)), ct.top)
    wins <- quantile(ct.top, 0.75)
    ct.top[ct.top > wins] <- wins

    graph.distance.matrix <- ct.top.ids
    for(i in 1:nrow(graph.distance.matrix)){
      graph.distance.matrix[i, ct.top.ids [i,]] <- ct.top[i, ]
    }
    return(graph.distance.matrix)
  }

  subsampling.cells <- function(graph, subsampled.groups, seed=NULL){
     subsampledGraph <- sccore::getClusterGraph(graph, as.factor(subsampled.groups), method="sum")
     E(subsampledGraph)$weight %<>% log1p() %>% range()
     graph.distance.matrix <- computeGraphDistance(subsampledGraph)
     dendr <- graph.distance.matrix %>% log1p() %>% as.dist() %>% hclust(method='ward.D')
     dend <- as.dendrogram(dendr)
     dendr <- TransferDend(dend, renameCluster=renameCluster, cls.groups = subsampled.groups)
     if(renameCluster){
      subsampled.groups <- dendr$new.groups
    }
    return(list(dend=dendr$dendrogram, clusters=subsampled.groups))
  }

  sr <- NULL
  if(is.null(sr)){
    sr <- parallel::mclapply(1:length(subsample.groups), function(i) subsampling.cells(graph, subsample.groups[[i]], seed=i), mc.cores=5)
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
        d <- lapply(d,cbm,fac=fac)
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
#' @return dendrogram
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

#' normalize tree based upperlevel annotation
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

#' Normalize tree based on within-factor mixing
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

#' pairwise normalization - normalize nodes based on pairwise-species (factor) comparison
#' @param d dendrogram
#' @param facPrefix factor prefix selected to color - e.g. c("human", "mouse")
#' @return dendrogram with normalized value
#' @export
pairwiseNormTree <- function(d, facPrefix, reNormalize=TRUE){
  cbm <- function(d){
    if(is.leaf(d)){
      if(reNormalize){
        normfac <- attr(d, "normFac")
        if(is.null(normfac)){
          stop("please normalize the tree first......")
        } else{
          normfac <- normfac[match(facPrefix, names(normfac))]
          normpercent <- normfac/sum(normfac)
          attr(d, "normPercentage") <- normpercent
        }
      } else{
        normpercent <- attr(d, "normPercentage")
      }
      return(d);
    } else{
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
          attr(d, "normPercentage") <- normpercent
        }
      } else{
        normpercent <- attr(d, "normPercentage")
      }
      return(d);
    }
  }
  cbm(d)
}

#' caculate branch entropy for normalized value
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
#' @param heuristicStop TRUE: stop at the nodes where the mixing cutoff criteria is met
#' @param FurthurStop: TRUE: one-stop furthur where the mixing cutoff criteria is met
#' @export
pruneTree <- function(dend, minCutoff=0.3, maxCutoff=0.4, heuristicStop=TRUE, FurthurStop=FALSE){

  is.branch.mixed <- function(dend, minCutoff=0.3, maxCutoff=0.4){
    is.mixed <- FALSE
    normPercentage <- attr(dend, "normPercentage")
    if(normPercentage[1] >= minCutoff & normPercentage[1] < maxCutoff & !is.na(normPercentage[1])){
      is.mixed <- TRUE
    }
    return(is.mixed)
  }

  is.father.of.subtree.to.merge <- function(dend, minCutoff, maxCutoff) {
    # this function checks if the subtree we wish to merge is the direct child of the current branch (dend) we entered the function
    is.father <- FALSE
    for (i in seq_len(length(dend))){
      if(is.branch.mixed(dend[[i]], minCutoff, maxCutoff) == FALSE) is.father <- TRUE
    }
    return(is.father)
  }

  search_mixed_subtree <- function(dend, minCutoff, maxCutoff, heuristicStop){
    if(!is.father.of.subtree.to.merge(dend, minCutoff=minCutoff, maxCutoff=maxCutoff)){
      for (i in seq_len(length(dend))){
        if(is.leaf(dend[[i]])){
          next()
        }
        dend[[i]] <- search_mixed_subtree(dend[[i]], minCutoff, maxCutoff, heuristicStop = heuristicStop)
      }
    } else{ # we'll merge
      subtree_loc <- 1
      if(heuristicStop==TRUE){
        branch.attributes <- attributes(dend)
        dend <- prune(dend, labels(dend)[-1])
        attr(dend, "cc") <- branch.attributes$cc
        attr(dend, "size") <- branch.attributes$size
        attr(dend, "edgePar") <- branch.attributes$edgePar
        attr(dend, "nodesinfo") <- branch.attributes$nodesinfo
        attr(dend, "percentage") <- branch.attributes$percentage
        attr(dend, "normPercentage") <- branch.attributes$normPercentage
        attr(dend, "stability") <- branch.attributes$stability
        attr(dend, "height_original") <- branch.attributes$height
        attr(dend, "leaf") <- TRUE
        attr(dend, "members") <- 1
        attr(dend, "label") <- labels(dend)[1]
        attr(dend, "height") <- 0
      } else{
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
          attr(dend[[subtree_loc]], "stability") <- branch.attributes$stability
          attr(dend[[subtree_loc]], "height_original") <- branch.attributes$height
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
          attr(dend[[subtree_loc+1]], "stability") <- branch.attributes$stability
          attr(dend[[subtree_loc+1]], "height_original") <- branch.attributes$height
        }
       }
      }
    }
    return(dend)
  }

  prune_mixed_subtree <- function(dend, minCutoff, maxCutoff){
    dend.label <- labels(dend)
    dend <- search_mixed_subtree(dend, minCutoff, maxCutoff, heuristicStop)
    label.updated <- labels(dend)

    if(length(dend.label)==length(label.updated)){
      return(dend)
    } else{
      dend <- prune_mixed_subtree(dend, minCutoff, maxCutoff)
    }
  }

  if(heuristicStop){
    new_dend <- search_mixed_subtree(dend, minCutoff, maxCutoff, heuristicStop)
  } else{
    new_dend <- prune_mixed_subtree(dend, minCutoff, maxCutoff)
  }

  # re-assign space
  #new_dend <- ladderize(new_dend, right=FALSE)
  new_dend <- ladderize(fix_members_attr.dendrogram(new_dend), right=FALSE)
  return(new_dend)
}

#' prune tree based on entropy and other criteria
#' @param dend dendrogram obj
#' @param cutoff cutoff value
#' @param mixing if turn on factor mixing criteria - at least
#' @export
pruneTreeEntropy <- function(dend, cutoff=2.9, mixing=FALSE, minCutoff=NULL){

  is.branch.mixed <- function(dend, cutoff){
    is.mixed <- FALSE

    if(mixing){
      if(is.null(minCutoff)){
        stop("please provide cutoff values")
      }
    }
    entropy <- attr(dend, "entropy")
    normPercentage <- attr(dend, "normPercentage")
    if(mixing){
      if(entropy > cutoff | length(normPercentage[normPercentage >= minCutoff]) > 2){
        is.mixed <- TRUE
      }
    } else{
      if(entropy > cutoff){
        is.mixed <- TRUE
      }
    }
    return(is.mixed)
  }

  is.father.of.subtree.to.merge <- function(dend, cutoff){
    # this function checks if the subtree we wish to merge is the direct child of the current branch (dend)
    #  we entered the function
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
          attr(dend[[subtree_loc]], "nodePar") <- branch.attributes$nodePar
          attr(dend[[subtree_loc]], "nodesinfo") <- branch.attributes$nodesinfo
          attr(dend[[subtree_loc]], "percentage") <- branch.attributes$percentage
          attr(dend[[subtree_loc]], "normPercentage") <- branch.attributes$normPercentage
          attr(dend[[subtree_loc]], "stability") <- branch.attributes$stability
          attr(dend[[subtree_loc]], "entropy") <- branch.attributes$entropy
          attr(dend[[subtree_loc]], "normFac") <- branch.attributes$normFac
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
          attr(dend[[subtree_loc+1]], "nodePar") <- branch.attributes$nodePar
          attr(dend[[subtree_loc+1]], "nodesinfo") <- branch.attributes$nodesinfo
          attr(dend[[subtree_loc+1]], "percentage") <- branch.attributes$percentage
          attr(dend[[subtree_loc+1]], "normPercentage") <- branch.attributes$normPercentage
          attr(dend[[subtree_loc+1]], "stability") <- branch.attributes$stability
          attr(dend[[subtree_loc+1]], "entropy") <- branch.attributes$entropy
          attr(dend[[subtree_loc+1]], "normFac") <- branch.attributes$normFac
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

#' Prune tree based on within species cluster purity
#'  search the nodes until it reach to a cutoff value for at lease two species
#' @param dend dendrogram obj
#' @param cutoff cutoff value
#' @param mixing if turn on factor mixing criteria - at least
#' @export
pruneTreeWithSpecies <- function(d, purityCutoff=0.8){

  is.branch.mixed <- function(d, purityCutoff){
    is.mixed <- FALSE
    withinSpeciesComp <- attr(d, "withinSpeciesComp")
    if(is.null(withinSpeciesComp)){
      stop("Please calculate within-species composition first")
    }

    if(length(withinSpeciesComp$human) == 0){
      human.withinSpeciesComp <- 0
    } else{
      human.withinSpeciesComp <- max(withinSpeciesComp$human)[1]
    }

    if(length(withinSpeciesComp$marmo) == 0){
      marmo.withinSpeciesComp <- 0
    } else{
      marmo.withinSpeciesComp <- max(withinSpeciesComp$marmo)[1]
    }

    if(length(withinSpeciesComp$mouse) == 0){
      mouse.withinSpeciesComp <- 0
    } else{
      mouse.withinSpeciesComp <- max(withinSpeciesComp$mouse)[1]
    }
    composition <- c(human.withinSpeciesComp, marmo.withinSpeciesComp, mouse.withinSpeciesComp)

    if(length(which(composition > purityCutoff)) > 2){
      is.mixed <- TRUE
    }
    return(is.mixed)
  }

  is.father.of.subtree.to.merge <- function(d, purityCutoff){
    # this function checks if the subtree we wish to merge is the direct child of the current branch (dend)
    #  we entered the function
    is.father <- FALSE
    for (i in seq_len(length(d))){
      if(is.branch.mixed(d[[i]], purityCutoff) == FALSE & !is.leaf(d[[i]])) is.father <- TRUE
    }
    return(is.father)
  }

  search_mixed_subtree <- function(d, purityCutoff){
    if (!is.father.of.subtree.to.merge(d, purityCutoff)){
      for (i in seq_len(length(d))){
        if(is.leaf(d[[i]])){
          next()
        }
        d[[i]] <- search_mixed_subtree(d[[i]], purityCutoff)
      }
    } else { # we'll merge
      subtree_loc <- 1
      if(!is.branch.mixed(d[[subtree_loc]], purityCutoff)){
        # achieve attributes
        branch.attributes <- attributes(d[[subtree_loc]])
        # get all the leaf labels in this subtree
        leaf.labels <- labels(d[[subtree_loc]])
        # prune subtree
        if(length(leaf.labels) > 1){
          d <- prune(d, leaf.labels[-1])
          # assign attributes
          attr(d[[subtree_loc]], "cc") <- branch.attributes$cc
          attr(d[[subtree_loc]], "size") <- branch.attributes$size
          attr(d[[subtree_loc]], "edgePar") <- branch.attributes$edgePar
          attr(d[[subtree_loc]], "nodePar") <- branch.attributes$nodePar
          attr(d[[subtree_loc]], "nodesinfo") <- branch.attributes$nodesinfo
          attr(d[[subtree_loc]], "percentage") <- branch.attributes$percentage
          attr(d[[subtree_loc]], "normPercentage") <- branch.attributes$normPercentage
          attr(d[[subtree_loc]], "stability") <- branch.attributes$stability
          attr(d[[subtree_loc]], "entropy") <- branch.attributes$entropy
          attr(d[[subtree_loc]], "normFac") <- branch.attributes$normFac
        }
      }

      if(!is.branch.mixed(d[[subtree_loc+1]], purityCutoff)){
        # achieve attributes
        branch.attributes <- attributes(d[[subtree_loc+1]])
        # get all the leaf labels in this subtree
        leaf.labels <- labels(d[[subtree_loc+1]])
        # prune subtree
        if(length(leaf.labels) > 1){
          d <- prune(d, leaf.labels[-1])
          # assign attributes
          attr(d[[subtree_loc+1]], "cc") <- branch.attributes$cc
          attr(d[[subtree_loc+1]], "size") <- branch.attributes$size
          attr(d[[subtree_loc+1]], "edgePar") <- branch.attributes$edgePar
          attr(d[[subtree_loc+1]], "nodePar") <- branch.attributes$nodePar
          attr(d[[subtree_loc+1]], "nodesinfo") <- branch.attributes$nodesinfo
          attr(d[[subtree_loc+1]], "percentage") <- branch.attributes$percentage
          attr(d[[subtree_loc+1]], "normPercentage") <- branch.attributes$normPercentage
          attr(d[[subtree_loc+1]], "stability") <- branch.attributes$stability
          attr(d[[subtree_loc+1]], "entropy") <- branch.attributes$entropy
          attr(d[[subtree_loc+1]], "normFac") <- branch.attributes$normFac
        }
      }
    }
    return(d)
  }

  prune_mixed_subtree <- function(d, purityCutoff){
    d.label <- labels(d)
    d <- search_mixed_subtree(d, purityCutoff)
    label.updated <- labels(d)

    if(length(d.label)==length(label.updated)){
      return(d)
    } else{
      #dend.label <- label.updated
      d <- prune_mixed_subtree(d, purityCutoff)
      #label.updated <- labels(dend)
    }
  }

  new_dend <- prune_mixed_subtree(d, purityCutoff)
  # re-assign space
  new_dend <- ladderize(new_dend, right=FALSE)
  new_dend <- ladderize(fix_members_attr.dendrogram(new_dend), right=FALSE)
  return(new_dend)
}

#' remove leaves with low stability score
#' @param d dendrogram obj
#' @param cutoff stability cutoff to prune the tree
#' @param sizeCutoff cut leaf nodes that have size below sizeCutoff
#' @param sizeOnly only cut tree based on the size of leaf nodes
#' @return trimmed dendrogram
removeLowStabilityLeaf <- function(d, cutoff=0.5, sizeOnly=FALSE, sizeCutoff=1){
  if(is.null(attr(d, "stability")) & sizeOnly==FALSE){
    stop("Please measure stability first ...")
  }
  if(sizeOnly){
    leaf.labels <- get_leaves_attr(d, "label")
    size <- get_leaves_attr(d, "size")
    leaf.info <- data.frame(labels=leaf.labels, size=size)
    leaf.info.filtered <- leaf.info[leaf.info$size<=sizeCutoff, ]$labels
    d <- prune(d, leaves=as.character(leaf.info.filtered))
  } else{
    leaf.stability <- get_leaves_attr(d, "stability")
    leaf.labels <- get_leaves_attr(d, "label")
    size <- get_leaves_attr(d, "size")
    leaf.info <- data.frame(labels=leaf.labels, stability=leaf.stability, size=size)
    leaf.info.filtered <- leaf.info[leaf.info$stability<=cutoff | leaf.info$size<=sizeCutoff, ]$labels
    d <- prune(d, leaves=as.character(leaf.info.filtered))
  }
  return(d)
}

#' Get robust homologous clusters from the tree based on species (factor) mixing.
#' Clusters are selected at the nodes where one species (factor) begin to differentiate from the other species (factors)
#' @param d pruned dendrogram
#' @param plotTree if TRUE, plot the dendrogram and its correspondent integrative clusters
#' @param plotleafLabel TRUE/FALSE, if TRUE visualize leaf labels
#' @param upperlevelannot upperlevel annotation, plot top level annotation in the tree
#' @param cell.group list contains specific cell group to visualize at the bottom of the tree
#' @param scale True: scale the cell groups according to their purity in the node
#' @param furtherStop if TRUE it will go one-step further beyond the current homologous nodes. In this case,
#'              it may get one/two species-specific cluster(s) and one/0 well-mixed clusters.
#' @param withinSpeciesStop fine tuning for the within-Species composition for the integrative clusters
#'                          in the furtherStop mode. Stop splitting of a node if the split of the node
#'                          will seperate the within-species clusters below certain percentage
#' @param withinSpeciesStopCutoff cutoff value for the purity of withinSpecies clusters
#' @param plotClusterTree TRUE: plot tree which each leaf node indicates the integrative clusters and directly return
#'                              dendrogram
#' @return factor contains robust cluster IDs and their correspondent cell IDs
#' @import RColorBrewer
#' @import colorspace
#' @export
getClusters <- function(d, plotTree=TRUE, plotleafLabel=FALSE, upperlevelannot=NULL, cell.group=NULL,
                        scale=TRUE, furtherStop=FALSE, withinSpeciesStop=FALSE, withinSpeciesStopCutoff=0.8,
                        plotClusterTree=FALSE){
  require("RColorBrewer")
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
    } else{ # get node information
      if(furtherStop){
        if(withinSpeciesStop){
          withinSpeciesAttr <- attr(d, "withinSpeciesComp")
          if(is.null(withinSpeciesAttr)){
            stop("Please provide within-species cluster information")
          }

          branch1WithinSpeciesAttr <- attr(d[[1]], "withinSpeciesComp")
          branch2WithinSpeciesAttr <- attr(d[[2]], "withinSpeciesComp")
          branch1WithinSpeciesPercent <- max(branch1WithinSpeciesAttr$human, branch1WithinSpeciesAttr$marmo,
                                          branch1WithinSpeciesAttr$mouse)
          branch2WithinSpeciesPercent <- max(branch2WithinSpeciesAttr$human, branch2WithinSpeciesAttr$marmo,
                                             branch2WithinSpeciesAttr$mouse)
          # if the composition of any of the leaf branch larger than cutoff, split
          if(branch1WithinSpciesPercent >= withinSpeciesStopCutoff | branch2WithinSpeciesPercent >= withinSpeciesStopCutoff){
            if(!is.leaf(d[[1]])){
              branch.attributes <- attributes(d[[1]])
              leafsize <- length(labels(d[[1]]))
              d[[1]] <- prune(d[[1]], labels(d[[1]])[-1])
              attributes(d[[1]]) <- branch.attributes[-2]
              attr(d[[1]], "leaf") <- TRUE
              attr(d[[1]], "members") <- 1
              attr(d[[1]], "label") <- labels(d)[1]
              attr(d[[1]], "height") <- 0
              attr(d[[1]], "leaflabels") <- labels(d[[1]])
              attr(d[[1]], "leafSize") <- leafsize
              attr(d[[2]], "leaflabels") <- labels(d[[2]])
              attr(d[[2]], "leafSize") <- length(labels(d[[2]]))
            } else if(!is.leaf(d[[2]])){
              branch.attributes <- attributes(d[[2]])
              leafsize <- length(labels(d[[2]]))
              d[[2]] <- prune(d[[2]], labels(d[[2]])[-1])
              attributes(d[[2]]) <- branch.attributes[-2]
              attr(d[[2]], "leaf") <- TRUE
              attr(d[[2]], "members") <- 1
              attr(d[[2]], "label") <- labels(d)[1]
              attr(d[[2]], "height") <- 0
              attr(d[[2]], "leaflabels") <- labels(d[[2]])
              attr(d[[2]], "leafSize") <- leafsize
              attr(d[[1]], "leaflabels") <- labels(d[[1]])
              attr(d[[1]], "leafSize") <- length(labels(d[[1]]))
            } else{
              attr(d[[1]], "leaflabels") <- labels(d[[1]])
              attr(d[[1]], "leafSize") <- length(labels(d[[1]]))
              attr(d[[2]], "leaflabels") <- labels(d[[2]])
              attr(d[[2]], "leafSize") <- length(labels(d[[2]]))
            }
          } else {
            #attr(d, "leaf") <- TRUE
            branch.attributes <- attributes(d)
            leafsize <- length(labels(d))
            d <- prune(d, labels(d)[-1])
            attributes(d) <- branch.attributes[-2]
            attr(d, "leaf") <- TRUE
            attr(d, "members") <- 1
            attr(d, "label") <- labels(d)[1]
            attr(d, "height") <- 0
            attr(d, "leaflabels") <- labels(d)
            attr(d, "leafSize") <- leafsize
          }

        } else{
          if(!is.leaf(d[[1]])){
            branch.attributes <- attributes(d[[1]])
            leafsize <- length(labels(d[[1]]))
            d[[1]] <- prune(d[[1]], labels(d[[1]])[-1])
            attributes(d[[1]]) <- branch.attributes[-2]
            attr(d[[1]], "leaf") <- TRUE
            attr(d[[1]], "leaflabels") <- labels(d[[1]])
            attr(d[[1]], "leafSize") <-  leafsize
            attr(d[[2]], "leaflabels") <- labels(d[[2]])
            attr(d[[2]], "leafSize") <- length(labels(d[[2]]))
          } else if(!is.leaf(d[[2]])){
            branch.attributes <- attributes(d[[2]])
            leafsize <- length(labels(d[[2]]))
            d[[2]] <- prune(d[[2]], labels(d[[2]])[-1])
            attributes(d[[2]]) <- branch.attributes[-2]
            attr(d[[2]], "leaf") <- TRUE
            attr(d[[2]], "leaflabels") <- labels(d[[2]])
            attr(d[[2]], "leafSize") <- leafsize
            attr(d[[1]], "leaflabels") <- labels(d[[1]])
            attr(d[[1]], "leafSize") <- length(labels(d[[1]]))
          } else{
            attr(d[[1]], "leaflabels") <- labels(d[[1]])
            attr(d[[1]], "leafSize") <- length(labels(d[[1]]))
            attr(d[[2]], "leaflabels") <- labels(d[[2]])
            attr(d[[2]], "leafSize") <- length(labels(d[[2]]))
          }
        }
      } else{
        attr(d, "leaf") <- TRUE
        attr(d, "leaflabels") <- labels(d)
        attr(d, "leafSize") <- length(labels(d))
      }
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
  labels(d.pruned) <- as.character(unique(clusters))

  if(plotClusterTree){
    d.pruned.ladder <- ladderize(fix_members_attr.dendrogram(d.pruned), right=FALSE)
    upperLevelnodes <- getUpperLevelNode(d.pruned.ladder, cutoff=0.65)
    plot(d.pruned.ladder)
    text(upperLevelnodes$xy, labels=upperLevelnodes$upperlabel, adj=c(0.5, -0.5), cex=1, col="red")
    return(d.pruned.ladder)
  }

  if(plotTree){
    clpalette <- function(n){
      all.colors <- colors()
      qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
      cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      cols <- c(cols, sample(cols))
      cols[1:n]
    }
    bars <- get_leaves_attr(d.pruned, "leafSize")
    colorbars <- unlist(lapply(1:length(bars), function(r){
      rep(r, bars[r])
    }))
    clusterColorbars <- clpalette(max(colorbars))
    clusterColorbars <- clusterColorbars[colorbars]

    if(plotleafLabel){
      plot(d)
    } else{
      plot(d, leaflab="none")
    }
    if(!is.null(upperlevelannot)){
      text(upperlevelannot$xy, labels=upperlevelannot$upperlabel, adj=c(0.5, -0.5), cex=1, col="red")
    }

    # visualize cell groups
    cellcolorbar <- NULL
    if(!is.null(cell.group)){
      # check the location of cell groups in leafnode
      getCellpos <- function(cells){
        leafcluster <- getLeafClusters(d)
        cluster.name <- names(table(leafcluster))
        cell.color.bar <- unlist(lapply(1:length(cluster.name), function(r){
          leafcells <- names(leafcluster[leafcluster==cluster.name[r]])
          cellProp <- intersect(cells, leafcells)
          if(length(cellProp)/min(length(cells), length(leafcells)) > 0.1){
            #1
            round(length(cellProp)/min(length(cells), length(leafcells)), 1)
          } else{
            #8
            0
          }
        }))
        return(cell.color.bar)
      }
      cellcolorbar <- dplyr::bind_cols(lapply(cell.group, getCellpos)) * 10
      # generate color pallete
      colorPallete <- colorspace::sequential_hcl(11, h = 260, c = c(100, 0), l = c(25, 100), rev = TRUE, power = 1)
      cellcolorbar <- apply(cellcolorbar, 2, function(r){colorPallete[r+1]})
    }
    if(!is.null(cellcolorbar)){
      thebar <- cbind(cellcolorbar, as.data.frame(clusterColorbars))
      names(thebar)[ncol(thebar)]  <- "Clusters"
    } else{
      thebar <- as.data.frame(clusterColorbars)
      names(thebar)[ncol(thebar)]  <- "Clusters"
    }
    colored_bars(colors=thebar, dend=d, y_shift=-0.06, sort_by_labels_order = FALSE, cex.rowLabels = 0.5)
  }
  return(clusters)
}

#' testing function: get homologous clusters for heuristic criteria 3
#' @param d pruned tree
#' @param plot TRUE plot the tree that correspondent to final homologous clusters
#' @return homologous clusters
get_divergent_heuristic_clusters <- function(d, plot=TRUE){
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
    } else{
      if(is.leaf(d[[1]]) & is.leaf(d[[2]])){
        branch.attributes <- attributes(d)
        d.label <- labels(d)[1]
        d <- prune(d, labels(d)[-1])
        attributes(d) <- branch.attributes[-2]
        attr(d, "leaf") <- TRUE
        attr(d, "members") <- 1
        attr(d, "label") <- d.label
        attr(d, "height") <- 0
      } else{
        if(!is.leaf(d[[1]])){
          d[[1]]<- search_tree(d[[1]])
        }
        if(!is.leaf(d[[2]])){
          d[[2]] <- search_tree(d[[2]])
        }
      }
    }
    return(d)
  }

  d.pruned <- search_tree(d)
  d.pruned <- ladderize(fix_members_attr.dendrogram(d.pruned), right=FALSE)
  if(plot){
    upperLevelnodes <- getUpperLevelNode(d.pruned, cutoff=0.65)
    plot(d.pruned)
    text(upperLevelnodes$xy, labels=upperLevelnodes$upperlabel, adj=c(0.5, -0.5), cex=1, col="red")
  }

  homo.clusters <- getLeafClusters(d.pruned)
  return(homo.clusters)
}


#' get leaf clusters
#' @param dend dendrogram object
#' @param renameClusterlabel if TRUE, rename the cluster labels at leaf nodes
#' @return cell clusters at leaf nodes
#' @export
getLeafClusters <- function(dend, renameClusterlabel=FALSE){
  if(renameClusterlabel){
    dend <- assign_values_to_nodes(dend, "label", seq(1, nnodes(dend), 1))
    clusters <- labels(dend)
  } else{
    clusters <- labels(dend)
  }
  cells <- get_leaves_attr(dend, "nodesinfo")
  sizes <- get_leaves_attr(dend, "size")
  clusters <- unlist(lapply(1:length(sizes), function(r){ rep(clusters[r], sizes[r])}))
  names(clusters) <- cells
  return(clusters)
}

#' build entire tree based on dendrogram and subsampled dendrograms automatically
#' @param dend dendrogram obj of original tree
#' @param expMatrix gene expression matrix - rows are cells
#' @param subsampleClusters list contains subsampled cell clusters
#' @param cls.groups clusters of original dendrograms
#' @param cellannot  cell annotations
#' @param species factor annoted which cell belong to which species
#' @param upperlevelannot higher-level cell annotations
#' @param plot plot the built tree or not
#' @export
buildSpeciesTree <- function(dend, expMatrix, subsampleClusters=NULL, cls.groups, mappingcells=FALSE, cellannot=NULL, species,
                             upperlevelannot=NULL, renameCluster=TRUE, plot=TRUE){
  dendr <- TransferDend(dend, renameCluster=renameCluster, cls.groups = cls.groups)
  cls.groups <- dendr$new.groups

  dend <- dendr$dendrogram
  leafcontent <- dendr$leafcontent

  if(!is.null(subsampleClusters)){
    subsampled.dend <- subSampleTree(expMatrix, subsample.groups=subsampleClusters)
    stability.measurements <- TreeStabilityDend(dend, cls.groups=cls.groups, subsampled.dend, n.cores=10)
    dend <- stability.measurements$dendrogram
    #stability.measurements <- TreeStabilityDend(dend, cls.groups, subsamples, n.cores=10)
  } else{
    stability.measurements = NULL
  }

  # add cluster attribute to dendrogram
  dend <- AddTreeAttribute(dend, species, leafcontent)
  dend <- dendSetWidthBysize(dend, scale=8)

  if(mappingcells){
    if(is.null(cellannot)){
      stop("Please provide cell annotations")
    }
    leaflabels <- mappingLabel(dend, leafcontent, cellannot, humanAnnot=T)
    # set labels
    dend <- set_labels(dend, paste(dend %>% labels(), leaflabels, sep=" "))
  }

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
  }
  if(plot){
    if(!is.null(upperlevelannot) & !is.null(subsamples)){
      par(cex=1, mar=c(20, 5, 0, 10))
      plot(dend)
      text(upperLevelnodes$xy, labels=upperLevelnodes$upperlabel, adj=c(0.5, -1.2), cex=0.8, col="red")
      text(get_nodes_xy(dend), labels=get_nodes_attr(dend, "stability"), adj=c(0.4,0.4), cex=0.3, col="red")
    } else if(!is.null(upperlevelannot)){
      par(cex=1, mar=c(20, 5, 0, 10))
      plot(dend)
      text(get_nodes_xy(dend), labels=get_nodes_attr(dend, "stability"), adj=c(0.4,0.4), cex=0.3, col="red")
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
  return(list(dend=dend, upperLevelnodes=upperLevelnodes, leafcontent=leafcontent))
}

#' adjust height for the tree based on branch level,
#' @param d dendrogram obj
#' @param heightDiff differences in height between two branches
#' @return dendrogram object that has the same height for each subcluster
adjustTreeheight <- function(d, scale=1){
  rootDepth <- max_depth(d)
  #step <- scale/rootDepth

  # re-asign labels to dendrogram in the depth first search order
  #   in this case, the label of root node will be 1
  d <- assign_values_to_nodes(d, "label", seq(1, nnodes(d), 1))

  # get all pathes from root to leave nodes
  subtrees <- partition_leaves(d)
  leafNodes <- subtrees[[1]]
  pathRoutes <- function(leafnodes) {
    which(sapply(subtrees, function(x) leafnodes %in% x))
  }
  paths <- lapply(leafNodes, pathRoutes)
  paths <- lapply(paths, function(r){
    names(r) <- seq(1, length(r), 1)
    return(r)
  })

  get_node_level <- function(node, paths){
    for(i in 1:length(paths)){
      if(length(paths[[i]][paths[[i]] == node]) > 0){
        return(names(paths[[i]][paths[[i]] == node]))
      }
    }
  }

  allnodes <- get_nodes_attr(d, "label")
  nodesLevel <- unlist(lapply(allnodes, function(r){
    get_node_level(r, paths)
  }))

  nodesLevel <- as.numeric(nodesLevel) * scale
  levelinfo <- data.frame(nodes=allnodes, nodesLevel=nodesLevel)
  levelinfo$height <- rootDepth - levelinfo$nodesLevel

  d <- assign_values_to_nodes(d, "height", levelinfo$height)
  return(d)
}

#' align multiple dendrograms in the same scale
#' @param treelist list of dendrograms
#' @param upperlevelinfo plot upperlevel annot if true
#' @param plotLeafLabel TRUE/FALSE show leaf label if TRUE
#' @param base scale to adjust plot width
#' @return dendrogram panel aligned at the same scale
alignTreePlot <- function(treelist, plotUpperlevelinfo=TRUE, plotLeafLabel=FALSE, label.size=1, base=10){
  maxdepth <- max(unlist(lapply(treelist, max_depth)))

  # adjust tree height to the same scale
  treePanel <- lapply(treelist, function(t){
    tree.depth <- max_depth(t)
    if(tree.depth < maxdepth){
      diff <- maxdepth - tree.depth
      height <- get_nodes_attr(t, "height")
      height <- height + diff
      t <- assign_values_to_nodes(t, "height", height)
      return(t)
    } else{
      return(t)
    }
  })

  # plot dendrogram
  plotTree <- function(d){
    if(plotLeafLabel){
      plot(d)
    } else{
      plot(d, leaflab="none")
    }
    if(plotUpperlevelinfo){
      upperLevelnodes <- getUpperLevelNode(d, cutoff=0.65)
      text(upperLevelnodes$xy, labels=upperLevelnodes$upperlabel, adj=c(0.5, -0.5), cex=1, col="red")
    }
  }

  # get number of leaves
  leafsizes <- unlist(lapply(treelist, function(r){
    length(get_leaves_attr(r, "label"))
  }))

  leafsizes <- round(leafsizes/base)
  #par(mfrow=c(1, length(treePanel)), cex=label.size)
  par(cex=label.size)
  layout(matrix(rep(seq(1:length(leafsizes)), leafsizes), nrow = 1, byrow = TRUE))
  for(i in 1:length(treePanel)){
    if(plotUpperlevelinfo){
      plotTree(treePanel[[i]])
    } else{
      plotTree(treePanel[[i]], plotUpperlevelinfo=FALSE)
    }
  }
}

#' Calculate leaf statistics in subtrees (below an upperlevel annot)
#' @param dendrogram obj with all the attributes
#' @param heightflag if TRUE, caculate leave depth based on the expression distance
#' @return # of leaves below an upperlevel branch
leavesSubtree <- function(d, heightflag=TRUE){
  # re-asign labels to dendrogram in the depth first search order
  #   in this case, the label of root node will be 1
  d <- assign_values_to_nodes(d, "label", seq(1, nnodes(d), 1))

  # get upperlevel annot nodes information
  upperLevelnodes <- getUpperLevelNode(d, cutoff=0.65)
  nodesLoc <- get_nodes_xy(d)
  labels <- get_nodes_attr(d, "label")
  nodes.info <- cbind(nodesLoc, labels)

  # extract upperlevel nodes info
  nodes.info <- merge(nodes.info, data.frame(upperLevelnodes$xy, upperLevelnodes$upperlabel), by=c(1, 2))
  colnames(nodes.info) <- c("x", "y", "labels", "cellannot")

  # get all pathes from root to leave nodes
  subtrees <- partition_leaves(d)
  leafNodes <- subtrees[[1]]
  pathRoutes <- function(leafnodes) {
    which(sapply(subtrees, function(x) leafnodes %in% x))
  }
  paths <- lapply(leafNodes, pathRoutes)
  paths <- lapply(paths, function(r){
    names(r) <- seq(1, length(r), 1)
  })

  # caculate number of leaves below a subtree
  nleaves.below.subtree <- lapply(1:nrow(nodes.info), function(n){
    node <- nodes.info[n, ]$labels
    nleafnodes <- length(which(unlist(lapply(paths, function(r){length(which(r==node))})) > 0))
    return(data.frame(celltype=nodes.info[n, ]$cellannot, nleaves=nleafnodes))
  })
  nleaves.below.subtree <- dplyr::bind_rows(nleaves.below.subtree)

  # caculate the depth of each leave below a subtree
  if(heightflag){
    root_height <- attr(d, "height")
    labels <- get_nodes_attr(d, "label")
    heights <- get_nodes_attr(d, "height_original")
    leaf.info <- data.frame(label=labels, height=heights)
    leaf.info[is.na(leaf.info)] <- 0

    dleaves.below.subtree <- lapply(1:nrow(nodes.info), function(n){
      node <- nodes.info[n, ]$labels
      nodePaths <- paths[which(unlist(lapply(paths, function(r){length(which(r==node))})) > 0)]
      dplyr::bind_rows(lapply(nodePaths, function(r){
        loc <- length(r)
        leafID <- r[loc]
        height <- round((root_height - leaf.info[leaf.info$label == leafID, ]$height), 2)
        data.frame(celltype=nodes.info[n, ]$cellannot, leaveID=r[loc], depth=height)
      }))
    })
  } else{
    dleaves.below.subtree <- lapply(1:nrow(nodes.info), function(n){
      node <- nodes.info[n, ]$labels
      nodePaths <- paths[which(unlist(lapply(paths, function(r){length(which(r==node))})) > 0)]
      dplyr::bind_rows(lapply(nodePaths, function(r){
        loc <- length(r)
        data.frame(celltype=nodes.info[n, ]$cellannot, leaveID=r[loc], depth=names(r[loc]))
      }))
    })
  }

  dleaves.below.subtree <- dplyr::bind_rows(dleaves.below.subtree)

  return(list(nleaves=nleaves.below.subtree, leaveDepth=dleaves.below.subtree))
}

#' Visualize leaf nodes in the tree based on specific cell groups
#' @param dend dendrogram obj
#' @param cells vector contains cellID
#' @param col node color
#' @param pure.cutoff if the purity of cell composition of a node below the cutoff, not visualize it
#' @return dendrogram obj that mark leaf nodes contain specific cells
#' @export
MarkCells <- function(dend, cells, col="blue", pure.cutoff=0.5){
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

#' fine-tuning for within-species clusters - specific for human, mouse and marmoset
#' Add attribute to monitor within-species composition
#' @param d dendrogram object
#' @return dendrogram with attribute contains within-species cluster percentage
#' @export
addWithinSpeciesComp <- function(d, purityCutoff=0.1, humanAnnot, mouseAnnot, marmoAnnot){
  humanClusterDistr <- table(humanAnnot)
  mouseClusterDistr <- table(mouseAnnot)
  marmoClusterDistr <- table(marmoAnnot)

  cbm <- function(d) {
    if(is.leaf(d)){
      nodeCells <- attr(d, "nodesinfo")
      humanCelldistr <- table(humanAnnot[names(humanAnnot) %in% nodeCells])
      mouseCelldistr <- table(mouseAnnot[names(mouseAnnot) %in% nodeCells])
      marmoCelldistr <- table(marmoAnnot[names(marmoAnnot) %in% nodeCells])

      # caculate percentage
      humanCelldistr <- humanCelldistr / humanClusterDistr[names(humanClusterDistr) %in% names(humanCelldistr)]
      mouseCelldistr <- mouseCelldistr / mouseClusterDistr[names(mouseClusterDistr) %in% names(mouseCelldistr)]
      marmoCelldistr <- marmoCelldistr / marmoClusterDistr[names(marmoClusterDistr) %in% names(marmoCelldistr)]

      # Only get within-species clusters with highest purity
      humanCelldistr <- humanCelldistr[humanCelldistr == max(humanCelldistr)]
      mouseCelldistr <- mouseCelldistr[mouseCelldistr == max(mouseCelldistr)]
      marmoCelldistr <- marmoCelldistr[marmoCelldistr == max(marmoCelldistr)]

      # get within-species clusters that above the purity cutoff
      humanCelldistr <- humanCelldistr[humanCelldistr>purityCutoff]
      mouseCelldistr <- mouseCelldistr[mouseCelldistr>purityCutoff]
      marmoCelldistr <- marmoCelldistr[marmoCelldistr>purityCutoff]

      # add attributes in the nodes
      attr(d, "withinSpeciesComp")$human <- humanCelldistr
      attr(d, "withinSpeciesComp")$marmo <- marmoCelldistr
      attr(d, "withinSpeciesComp")$mouse <- mouseCelldistr
      return(d)
    } else {
      oa <- attributes(d)
      d <- lapply(d,cbm)
      attributes(d) <- oa
      nodeCells <- attr(d, "nodesinfo")
      humanCelldistr <- table(humanAnnot[names(humanAnnot) %in% nodeCells])
      mouseCelldistr <- table(mouseAnnot[names(mouseAnnot) %in% nodeCells])
      marmoCelldistr <- table(marmoAnnot[names(marmoAnnot) %in% nodeCells])
      # caculate percentage
      humanCelldistr <- humanCelldistr / humanClusterDistr[names(humanClusterDistr) %in% names(humanCelldistr)]
      mouseCelldistr <- mouseCelldistr / mouseClusterDistr[names(mouseClusterDistr) %in% names(mouseCelldistr)]
      marmoCelldistr <- marmoCelldistr / marmoClusterDistr[names(marmoClusterDistr) %in% names(marmoCelldistr)]
      # Only get within-species clusters with highest purity
      humanCelldistr <- humanCelldistr[humanCelldistr == max(humanCelldistr)]
      mouseCelldistr <- mouseCelldistr[mouseCelldistr == max(mouseCelldistr)]
      marmoCelldistr <- marmoCelldistr[marmoCelldistr == max(marmoCelldistr)]
      # get within-species clusters that above the purity cutoff
      humanCelldistr <- humanCelldistr[humanCelldistr>purityCutoff]
      mouseCelldistr <- mouseCelldistr[mouseCelldistr>purityCutoff]
      marmoCelldistr <- marmoCelldistr[marmoCelldistr>purityCutoff]
      # add attributes in the nodes
      attr(d, "withinSpeciesComp")$human <- humanCelldistr
      attr(d, "withinSpeciesComp")$marmo <- marmoCelldistr
      attr(d, "withinSpeciesComp")$mouse <- mouseCelldistr
      return(d)
    }
  }
  cbm(d)
}

#' get marker genes that differentiate expressed between two groups - adapted from pagoda2 getDifferentialGenes function
#' @param count gene expression matrix: rows are genes and columns are cells
#' @param ntop get top n marker genes
#' @param groups vector contains cell and group information
#' @return lists of genes that differentiate two cell populations
#' @export
GetDifferentialGenes <- function(count, groups=NULL, upregulated.only=FALSE, z.threshold=3, ntop=15){
  if(is.null(groups)){
    stop("Please provide cell and group information")
  }

  cm <- t(count)
  #taking subset of cells
  cm <- cm[rownames(cm) %in% names(groups),]
  groups <- as.factor(groups[match(rownames(cm), names(groups))])

  lower.lpv.limit <- -100
  # calculate rank per-column (per-gene) average rank matrix
  xr <- pagoda2:::sparse_matrix_column_ranks(cm)
  # calculate rank sums per group
  grs <- pagoda2:::colSumByFac(xr,as.integer(groups))[-1,,drop=F]
  # calculate number of non-zero entries per group
  xr@x <- numeric(length(xr@x))+1
  gnzz <- pagoda2:::colSumByFac(xr,as.integer(groups))[-1,,drop=F]

  group.size <- as.numeric(tapply(groups, groups, length))[1:nrow(gnzz)]; group.size[is.na(group.size)] <- 0
  # add contribution of zero entries to the grs
  gnz <- (group.size-gnzz)
  # rank of a 0 entry for each gene
  zero.ranks <- (nrow(xr)-diff(xr@p)+1)/2 # number of total zero entries per gene
  ustat <- t((t(gnz)*zero.ranks)) + grs - group.size*(group.size+1)/2
  # standardize
  n1n2 <- group.size*(nrow(cm)-group.size)
  # correcting for 0 ties, of which there are plenty
  usigma <- sqrt(n1n2*(nrow(cm)+1)/12)
  usigma <- sqrt((nrow(cm) +1 - (gnz^3 - gnz)/(nrow(cm)*(nrow(cm)-1)))*n1n2/12)
  x <- t((ustat - n1n2/2)/usigma); # standardized U value- z score

  # correct for multiple hypothesis
  x <- matrix(qnorm(pagoda2:::bh.adjust(pnorm(as.numeric(abs(x)), lower.tail = FALSE, log.p = TRUE), log = TRUE),
                    lower.tail = FALSE, log.p = TRUE),ncol=ncol(x))*sign(x)
  rownames(x) <- colnames(cm); colnames(x) <- levels(groups)[1:ncol(x)];

 # add fold change information
  log.gene.av <- log2(Matrix::colMeans(cm));
  group.gene.av <- pagoda2:::colSumByFac(cm,as.integer(groups))[-1,,drop=F] / (group.size+1);
  log2.fold.change <- log2(t(group.gene.av)) - log.gene.av;
  # fraction of cells expressing
  f.expressing <- t(gnzz / group.size);
  max.group <- max.col(log2.fold.change)

  ds <- lapply(1:ncol(x),function(i){
    z <- x[,i];
    vi <- which((if (upregulated.only) z else abs(z)) >= z.threshold);
    r <- data.frame(Z=z[vi],M=log2.fold.change[vi,i],highest=max.group[vi]==i,fe=f.expressing[vi,i], Gene=rownames(x)[vi])
    rownames(r) <- r$Gene;
    r <- r[order(r$Z,decreasing=T),]
    r
 })
   names(ds) <- colnames(x)
   diffgene.list <- list(branch1=ds[[1]][1:ntop, ], branch2=ds[[2]][1:ntop, ])
   return(diffgene.list)
}

#' add marker genes that differentiate sister branches
#' @param d dendrogram object
#' @param count count matrix - rows are genes and columns are cells
#' @param addAUC estimate AUC based on random forest using expression of markers
#'               may take long time to run.
#' @export
AddMarkersToTree <- function(d, count, addAUC=TRUE){
  calculate_marker <- function(d){
    for(i in seq_len(length(d))){
      if(!is.leaf(d[[i]])){
        cells <- c(attr(d[[i]][[1]], "nodesinfo"), attr(d[[i]][[2]], "nodesinfo"))
        groups <- setNames(c(rep("1", length(attr(d[[i]][[1]], "nodesinfo"))),
                             rep("2", length(attr(d[[i]][[2]], "nodesinfo")))), cells)

        #print("calculate markers for each node................")
        diffgene.list <- GetDifferentialGenes(count, groups=groups)
        attr(d[[i]][[1]], "marker") <- diffgene.list[[1]]
        attr(d[[i]][[2]], "marker") <- diffgene.list[[2]]

        markers <- c(as.character(diffgene.list[[1]]$Gene), as.character(diffgene.list[[2]]$Gene))

        # calculate AUC
        if(addAUC){
          #print("run random forest on markers based on 3-fold CV......")
          counts <- t(count)
          counts.bin <- (counts[names(groups), markers, drop=F] > 0)
          #counts.markers <- as.data.frame(counts[names(groups), markers, drop=F])
          #counts.markers$group <- "1"
          #counts.markers[rownames(counts.markers) %in% names(groups[groups==2]), ]$group <- "2"

          ### divide training and testing group randomly (20% for testing)
          #sampleTesting <- sample(1:nrow(counts.markers), floor(nrow(counts.markers)*0.2))
          #trainData <- counts.markers[-sampleTesting, ]
          #testData <- counts.markers[sampleTesting, ]

          #rfModel <- caret::train(trainData[, -1*which(colnames(trainData) == "group")],
          #             as.factor(trainData$group), method="rf",
          #             trControl = caret::trainControl(method = "cv",number = 3))
          #rfPredict <- predict(rfModel, testData[, -1*which(colnames(testData) == "group")])
          #auc <- pROC::auc(pROC::roc(testData$group, as.numeric(rfPredict)))[1]

          counts.bin.Genesums <- Matrix::colSums(counts.bin)
          counts.bin.Groupsums <- Matrix::colSums(counts.bin & (groups == names(table(groups))))
          auc <- mean(apply(counts.bin, 2, function(col) pROC::auc(as.integer(groups), as.integer(col))))
          attr(d[[i]], "auc") <- auc
        }

        d[[i]] <- calculate_marker(d[[i]])
      }
    }
    return(d)
  }
  d <- calculate_marker(d)
  return(d)
}

#' add assoiated genes that differentiate cells from one node from the other cells
#' @param d dendrogram object
#' @param count count matrix - rows are genes and columns are cells from one species
#' @param prefix based factor(species) to extract associated genes
#'
#' @export
AddassGeneToTree <- function(d, count, prefix="human", ntop=30){
  calculate_marker <- function(d){
    for(i in seq_len(length(d))){
      if(!is.leaf(d[[i]])){
        cells <- attr(d[[i]], "nodesinfo")
        cells <- cells[cells %in% colnames(count)]
        cellsOut <- colnames(count)[-1*which(colnames(count) %in% cells)]
        groups <- setNames(c(rep("1", length(cells)), rep("2", length(cellsOut))), c(cells, cellsOut))

        #print("calculate markers for each node................")
        diffgene.list <- GetDifferentialGenes(count, groups=groups, ntop=ntop)
        attr(d[[i]], "assGenes") <- diffgene.list[[1]]
        d[[i]] <- calculate_marker(d[[i]])
      }
    }
    return(d)
  }
  d <- calculate_marker(d)
  return(d)
}

#' Binormial proportion test to identify genes show different expression patterns across factor/species based on
#'   the cells at each node of the tree - (not depend on marker genes)
#' @param d dendrogram object
#' @param count combined raw count matrix - rows are cells, columns are genes
#' @param prefix cell prefix to make the comparison - default: human vs marmoset and human vs mouse
#' @param sampleCells if TRUE, sampled cells based on the cell distribution in group
#' @return dendrogram with proportion test attributes
#' @export
AddProportionTest <- function(d, count, prefix=c("human", "marmoset", "mouse"), sampleCell=TRUE,
                              group=NULL, n.cores=20){
  if(length(prefix) < 3){
    stop("currently only support two group comparison")
  }
  propTestNode <- function(d){
    for(i in seq_len(length(d))){
      if(!is.leaf(d[[i]])){
        cells <- attr(d[[i]], "nodesinfo")
        cellsOut <- rownames(count)[-which(rownames(count) %in% cells)]
        if(sampleCell){
          if(is.null(group)){
            stop("please provide group annotation")
          }
          cells <- sampleCells(cells, group)
          cellsOut <- sampleCells(cellsOut, group)
        }
        group1Cells <- cells[grep(prefix[1], cells)]
        group2Cells <- cells[grep(prefix[2], cells)]
        group3Cells <- cells[grep(prefix[3], cells)]
        group1OutCells <- cellsOut[grep(prefix[1], cellsOut)]
        group2OutCells <- cellsOut[grep(prefix[2], cellsOut)]
        group3OutCells <- cellsOut[grep(prefix[3], cellsOut)]
        cellannot <- setNames(c(rep(prefix[1], length(group1Cells)), rep(prefix[2], length(group2Cells)),
                                rep(prefix[3], length(group3Cells))), c(group1Cells, group2Cells, group3Cells))
        cellannotOut <- setNames(c(rep(prefix[1], length(group1OutCells)), rep(prefix[2], length(group2OutCells)),
                                   rep(prefix[3], length(group3OutCells))), c(group1OutCells, group2OutCells, group3OutCells))

        # gene expression value for each factor(species)
        expByGroup <- conos:::collapseCellsByType(count, cellannot) %>% apply(2, function(r){
          r/sum(r)
        }) %>% t
        #expByGroup <- expByGroup[!is.na(expByGroup[, 1]), ]
        colnames(expByGroup) <- paste(names(table(cellannot)), c("exp"), sep="_")

        # caculate non-zero count per factor(species)
        count.binary <- count
        count.binary@x <- numeric(length(count.binary@x))+1
        groupProp <- conos:::collapseCellsByType(count.binary, cellannot) %>% t %>%
          apply(1, function(r){
            if(min(c(length(group1Cells),length(group2Cells), length(group3Cells))) > 0){
              r/c(length(group1Cells),length(group2Cells), length(group3Cells))
            } else if(length(group2Cells) > 0){
              r/c(length(group1Cells),length(group2Cells))
            } else{
              r/c(length(group1Cells),length(group3Cells))
            }}) %>% t
        colnames(groupProp) <- paste(names(table(cellannot)), c("propIn"), sep="_")

        # caculate non-zero count per factor for cellsOUT groups
        cellsOutCount <- count[rownames(count) %in% cellsOut, ]
        cellsOutCount@x <- numeric(length(cellsOutCount@x))+1
        groupOutProp <- conos:::collapseCellsByType(cellsOutCount, cellannotOut) %>% t %>%
          apply(1, function(r){r/c(length(group1OutCells),
                                   length(group2OutCells), length(group3OutCells))}) %>% t
        colnames(groupOutProp) <- paste(prefix, c("propOut"), sep="_")
        geneFeatures <- cbind(expByGroup, groupProp, groupOutProp) %>% as.data.frame
        # remove NA rows if any
        if(length(which(is.na(geneFeatures))) > 0){
          geneFeatures <- geneFeatures[!is.na(geneFeatures[, 1]), ]
        }

        if(min(c(length(group1Cells), length(group2Cells), length(group3Cells)))>0){
          # filtering criteria (expression of any species larger than 0.05,
          #  expression proportion larger than 0.15 quantile for each species)
          geneFeatures <- geneFeatures[geneFeatures[, 1]>quantile(geneFeatures[, 1])[2] &
                                         geneFeatures[, 2]>quantile(geneFeatures[, 2])[2]
                                       & geneFeatures[, 3]>quantile(geneFeatures[, 3])[2], ]
          geneFeatures <- geneFeatures[geneFeatures[, 4] > 0.3 &
                                         geneFeatures[, 5] > 0.3 &
                                         geneFeatures[, 6] > 0.3,]
          geneFeatures$gene <- rownames(geneFeatures)

          # proportion test
          countSelected <- count[, colnames(count) %in% rownames(geneFeatures)]
          propTest1 <- propTestGene(countSelected, group1Cells, group2Cells, n.cores=n.cores)
          propTest1 <- propTest1[!is.na(propTest1$stat), ]
          propTest2 <- propTestGene(countSelected, group1Cells, group3Cells, n.cores=n.cores)
          propTest2 <- propTest1[!is.na(propTest2$stat), ]

          # combine features
          propTest1 <- merge(propTest1, geneFeatures, by=c("gene"))
          propTest2 <- merge(propTest2, geneFeatures, by=c("gene"))

          # summarize divergence
          div1 <- propTest1[propTest1$p.adjust < 1e-5, ]
          div2 <- propTest2[propTest2$p.adjust < 1e-5, ]
          conv1 <- propTest1[propTest1$p.adjust >= 1e-5, ]
          conv2 <- propTest1[propTest2$p.adjust >= 1e-5, ]

          conservation <- c(length(intersect(conv1$gene, conv2$gene)),
                            (min(nrow(propTest1), nrow(propTest2))-length(intersect(conv1$gene, conv2$gene))-length(intersect(div1$gene, div2$gene))),
                            length(intersect(div1$gene, div2$gene)))
          names(conservation) <- c("conserved in 3", "conserved in 2", "divergent")
          attr(d[[i]], "proptest") <- list(propTest1=propTest1, propTest2=propTest2)
          attr(d[[i]], "conservation") <- conservation
        } else{
          geneFeatures <- geneFeatures[geneFeatures[, 1]>quantile(geneFeatures[, 1])[2] &
                                         geneFeatures[, 2]>quantile(geneFeatures[, 2])[2], ]
          geneFeatures <- geneFeatures[geneFeatures[, 4] > 0.3 &
                                         geneFeatures[, 5] > 0.3,]
          geneFeatures$gene <- rownames(geneFeatures)

          # proportion test
          countSelected <- count[, colnames(count) %in% rownames(geneFeatures)]

          if(length(group3Cells) == 0){
            propTest1 <- propTestGene(countSelected, group1Cells, group2Cells, n.cores=n.cores)
            propTest1 <- propTest1[!is.na(propTest1$stat), ]
          } else{
            propTest1 <- propTestGene(countSelected, group1Cells, group3Cells, n.cores=n.cores)
            propTest1 <- propTest1[!is.na(propTest1$stat), ]
          }

          # combine features
          propTest1 <- merge(propTest1, geneFeatures, by=c("gene"))

          # summarize conservation
          div1 <- propTest1[propTest1$p.adjust < 1e-5, ]
          conv1 <- propTest1[propTest1$p.adjust >= 1e-5, ]

          conservation <- c(0, nrow(conv1), nrow(div1))
          names(conservation) <- c("conserved in 3", "conserved in 2", "divergent")
          attr(d[[i]], "proptest") <- list(propTest1=propTest1)
          attr(d[[i]], "conservation") <- conservation
        }
        print(attr(d[[i]], "size"))
        d[[i]] <- propTestNode(d[[i]])
      }
    }
    return(d)
  }
  d <- propTestNode(d)
  return(d)
}

#' caculate divergent ratio for each node in the tree based on associated genes
#'   based on binormial proportion test
#' @param d dendrogram obj
#' @param count combined raw count matrix
#' @param prefix cell prefix to make the comparison
#' @param sampleCells sample cells in each node based on the cell composition in leaf node
#' @param group factor contains cells and its correspondent groups
#' @return dendrogram obj with divergence ratio added to each node
TreeDivergentRatio <- function(d, count, prefix=c("human", "marmoset", "mouse"), sampleCells=TRUE,
                               group=NULL, n.cores=20){
  if(length(prefix) < 3){
    stop("currently only support two group comparison")
  }
  calculate_divergence <- function(d){
    for(i in seq_len(length(d))){
      if(!is.leaf(d[[i]])){
        cells <- attr(d[[i]], "nodesinfo")
        assGenes <- attr(d[[i]], "assGenes")$Gene
        if(sampleCells){
          if(is.null(group)){
            stop("please provide cell annotation")
          }
          cells <- sampleCells(cells, groups)
        }

        group1Cells <- cells[grep(prefix[1], cells)]
        group2Cells <- cells[grep(prefix[2], cells)]
        group3Cells <- cells[grep(prefix[3], cells)]
        assCount <- count[, assGenes]
        divergentRatio1 <- 0
        divergentRatio2 <- 0

        if(length(group2Cells) > 0){
          propTest1 <- propTestGene(assCount, group1Cells, group2Cells, n.cores=n.cores)
          if(length(which(is.na(propTest1$p.value))) > 0){
            propTest1 <- propTest2[-1*which(is.na(propTest1$p.value)), ]
          }
          divergentRatio1 <- round(nrow(propTest1[propTest1$p.adjust<1e-5, ])/length(assGenes), 2)
        }

        if(length(group3Cells) > 0){
          propTest2 <- propTestGene(assCount, group1Cells, group3Cells, n.cores=n.cores)
          if(length(which(is.na(propTest2$p.value))) > 0){
            propTest2 <- propTest2[-1*which(is.na(propTest2$p.value)), ]
          }
          divergentRatio2 <- round(nrow(propTest2[propTest2$p.adjust<1e-5, ])/length(assGenes), 2)
        }
        attr(d[[i]], "PropTest") <- list(test1=propTest1, test2=propTest2)
        attr(d[[i]], "divergentRatio") <- c(divergentRatio1, divergentRatio2)
        d[[i]] <- calculate_divergence(d[[i]])
      }
    }
    return(d)
  }
  d <- calculate_divergence(d)
  return(d)
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

