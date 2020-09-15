#' add marker genes to gene based on global comparison
#' @param d dendrogram object
#' @param count count matrix - rows are genes and columns are cells
#' @param ntop get ntop markers for each node
#' @param evaltop mean auc based on evaltop genes
#' @export
AddassGeneToTreeGlobal <- function(d, count, ntop=30, addAUC=TRUE,
                             AUCutoff=0.6, evaltop=5){
  # initialize depth
  calculate_marker <- function(d){
    for(i in seq_len(length(d))){
      #if(!is.leaf(d[[i]])){
      if(!is.null(attr(d[[i]], "height"))){
        #print(i)
        cells <- attr(d[[i]], "nodesinfo")
        cells <- cells[cells %in% colnames(count)]
        cellsOut <- colnames(count)[-1*which(colnames(count) %in% cells)]
        groups <- setNames(c(rep("1", length(cells)), rep("2", length(cellsOut))), c(cells, cellsOut))
        
        #print("calculate markers for each node................")
        diffgene.list <- GetDifferentialGenes(count, groups=groups, ntop=ntop, upregulated.only=TRUE)
        markers <- diffgene.list[[1]]
        if(length(which(is.na(markers$Gene)))>0){
          markers <- markers[-which(is.na(markers$Gene)), ]
        }
        
        if(addAUC){
          counts <- Matrix::t(count)
          counts.bin <- (counts[names(groups), markers$Gene, drop=F] > 0)
          counts.bin.Genesums <- Matrix::colSums(counts.bin)
          counts.bin.Groupsums <- Matrix::colSums(counts.bin & (groups == names(table(groups))))
          auc <- apply(counts.bin, 2, function(col) pROC::auc(as.integer(groups), as.integer(col)))
          markers$AUC <- auc
          markers <- markers[order(markers$AUC, decreasing=TRUE), ]
          gene.sel <- markers[markers$AUC>AUCutoff, ]
          
          attr(d[[i]], "auc") <- round(mean(markers[1:evaltop, ]$AUC), 2)
          attr(d[[i]], "passedfilter") <- nrow(gene.sel)
          attr(d[[i]], "passedfilterGene") <- rownames(gene.sel)
        }
        attr(d[[i]], "assGenes") <- markers
        d[[i]] <- calculate_marker(d[[i]])
      }
    }
    return(d)
  }
  d <- calculate_marker(d)
  return(d)
}

#' add tree depth
#' @param d dendrogram object
#' @export
addDepth <- function(d){
  # initialize depth
  attr(d, "depth") <- 0
  attr(d[[1]], "depth") <- 1
  attr(d[[2]], "depth") <- 1
  for(i in seq_len(length(d))){
    calDepth <- function(d){
      if(!is.leaf(d[[i]])){
        attr(d[[i]][[1]], "depth") <- attr(d[[i]], "depth") + 1
        attr(d[[i]][[2]], "depth") <- attr(d[[i]], "depth") + 1
      }
      d[[i]] <- calDepth(d[[i]])
    }
  }
  return(d)
}

#' add AUC to gene lists
#' @param count count matrix - rows are genes and columns are cells
#' @param genes vector of genes to calculate AUC
calAUC <- function(count, groups, genes){
  counts <- Matrix::t(count)
  counts.bin <- (counts[names(groups), genes, drop=F] > 0)
  counts.bin.Genesums <- Matrix::colSums(counts.bin)
  counts.bin.Groupsums <- Matrix::colSums(counts.bin & (groups == names(table(groups))))
  if (requireNamespace("pROC", quietly = TRUE)){
    auc <- suppressWarnings(apply(counts.bin, 2, function(col) pROC::auc(as.integer(groups), as.integer(col), quiet = TRUE)))
  } else{
    message("please install pROC first")
  }
  if(length(which(auc<0.5))>0){
    auc[auc<0.5] <- 1- auc[auc<0.5]
  }
  return(auc)
}

#' add marker genes to gene based on pairwise comparison
#' @export
AddassGeneToTreePairwise <- function(d, count, ntop=30, addAUC=TRUE,
                                     AUCutoff=0.6, evaltop=5){
  
  calculate_marker <- function(d){
    for(i in seq_len(length(d))){
      #print(paste(attr(d[[i]], "depth"), i, sep=" "))
      if(!is.leaf(d[[i]])){
        cells <- c(attr(d[[i]][[1]], "nodesinfo"), attr(d[[i]][[2]], "nodesinfo"))
        groups <- setNames(c(rep("1", length(attr(d[[i]][[1]], "nodesinfo"))),
                             rep("2", length(attr(d[[i]][[2]], "nodesinfo")))), cells)
        
        #print("calculate markers for each node................")
        diffgene.list <- GetDifferentialGenes(count, groups=groups, ntop=ntop, upregulated.only=TRUE)
        diffgene.list <- lapply(diffgene.list, function(r){
          if(length(which(is.na(r$Gene)))>0){
            r[-which(is.na(r$Gene)), ]
          } else{
            r
          }
        })
        
        #markers <- c(as.character(diffgene.list[[1]]$Gene), as.character(diffgene.list[[2]]$Gene))
        if(addAUC){
          diffgene.list[[1]]$AUC <- calAUC(count, groups, diffgene.list[[1]]$Gene)
          diffgene.list[[1]] <- diffgene.list[[1]][order(diffgene.list[[1]]$AUC, decreasing=T), ]
          diffgene.list[[2]]$AUC <- calAUC(count, groups, diffgene.list[[2]]$Gene)
          diffgene.list[[2]] <- diffgene.list[[2]][order(diffgene.list[[2]]$AUC, decreasing=T), ]
          
          attr(d[[i]][[1]], "auc_pairwise") <- round(mean(diffgene.list[[1]][1:min(evaltop, nrow(diffgene.list[[1]])),
                                                                    ]$AUC), 2)
          attr(d[[i]][[2]], "auc_pairwise") <- round(mean(diffgene.list[[2]][1:min(evaltop, nrow(diffgene.list[[2]])), ]$AUC), 2)
        }
        attr(d[[i]][[1]], "assGenes_pairwise") <- diffgene.list[[1]]
        attr(d[[i]][[2]], "assGenes_pairwise") <- diffgene.list[[2]]
        d[[i]] <- calculate_marker(d[[i]])
      }
    }
    return(d)
  }
  # initalization
  cells <- c(attr(d[[1]], "nodesinfo"), attr(d[[2]], "nodesinfo"))
  groups <- setNames(c(rep("1", length(attr(d[[1]], "nodesinfo"))),
                       rep("2", length(attr(d[[2]], "nodesinfo")))), cells)
  diffgene.list <- GetDifferentialGenes(count, groups=groups, ntop=ntop, upregulated.only=TRUE)
  if(addAUC){
    diffgene.list[[1]]$AUC <- calAUC(count, groups, diffgene.list[[1]]$Gene)
    diffgene.list[[1]] <- diffgene.list[[1]][order(diffgene.list[[1]]$AUC, decreasing=T), ]
    diffgene.list[[2]]$AUC <- calAUC(count, groups, diffgene.list[[2]]$Gene)
    diffgene.list[[2]] <- diffgene.list[[2]][order(diffgene.list[[2]]$AUC, decreasing=T), ]
    
    attr(d[[1]], "auc_pairwise") <- round(mean(diffgene.list[[1]][1:evaltop, ]$AUC), 2)
    attr(d[[2]], "auc_pairwise") <- round(mean(diffgene.list[[2]][1:evaltop, ]$AUC), 2)
  }
  attr(d[[1]], "assGenes_pairwsie") <- diffgene.list[[1]]
  attr(d[[2]], "assGenes_pairwise") <- diffgene.list[[2]]
  
  d <- calculate_marker(d)
  return(d)
}

get_max_overlap_group_tree <- function(dend, target.de, top=20){
  calculate_overlap <- function(d){
    for(i in seq_len(length(d))){
      if(!is.leaf(d[[i]])){
        #print(i)
        markers <- attr(d[[i]], "assGenes")
        overlap <- 0
        overlaprate <- lapply(1:length(target.de), function(r){
          celltype <- names(target.de)[r]
          target.de.gene <- target.de[[r]][1:top, 1]
          de.gene <- rownames(markers[1:top, ])
          length(intersect(target.de.gene, de.gene))/top
        }) %>% unlist
        names(overlaprate) <- names(target.de)
        max.overlap <- max(overlaprate)
        max.cl <- names(overlaprate)[which(overlaprate==max.overlap)]
        if(max.overlap==0){
          attr(d[[i]], "overlap") <- paste("NULL", 0, sep=":")
        } else{
          attr(d[[i]], "overlap") <- paste(max.cl[1], max.overlap[1], sep=":")
        }
        d[[i]] <- calculate_overlap(d[[i]])
      }
    }
    return(d)
  }
  d <- calculate_overlap(d)
  return(d)
}

#' prune leaf node without significant marker genes
#' @export
initialPrune <- function(d){
  leaf.labels <- get_leaves_attr(d, "label")
  auc <- get_leaves_attr(d, "auc")
  leaf.info <- data.frame(labels=leaf.labels, auc=auc)
  leaf.info.filtered <- leaf.info[which(is.na(leaf.info$auc)), ]$labels
  d <- prune(d, leaves=as.character(leaf.info.filtered))
  return(d)
}

#' search hierarchical results to get optimal clusters based on marker separation power 
#' for each path from root to leaf, search the optimal node with maximum power to separate branches
#' @param method method for heuristic merging
#'                ave: if average auc of direct child nodes smaller then parent node, then merging
#'                all: if auc of all direct child nodes smaller than parent node, then merging
#'                any: if auc of any direct child node larger than parent node, then keep it
#'                anyFurthur: continue if both of the child nodes has auc larger than futhurCutoff
#' @param opt marker optimization method when do heuristic merging
#' @param hybridMethod vector of length 2, indicate the heuristic criteria - ave, all, any, anyFurthur
#' @export
heuristicMerge <- function(d, method=c("ave", "all", "any", "anyFurthur"), 
                           opt=c("global", "pairwise", "hybrid"), hybridMethod=NULL, furthurCutoff=0.8){
  is.good.split <- function(d, method=method, opt=opt){
    good.split <- FALSE
    if(is.leaf(d)){
      return(good.split)
    }
    if(opt == "global"){
      auc.left.child <- attr(d[[1]], "auc")
      auc.right.child <- attr(d[[2]], "auc")
      auc <- attr(d, "auc")
    } else if(opt == "pairwise"){
      auc.left.child <- attr(d[[1]], "auc_pairwise")
      auc.right.child <- attr(d[[2]], "auc_pairwise")
      auc <- attr(d, "auc_pairwise")
    } 
    
    if(is.na(auc.left.child)){
      auc.left.child <- 0
    }
    if(is.na(auc.right.child)){
      auc.right.child <- 0
    }
    if(length(is.na(auc))==0){
      if(is.null(auc)){
        auc <- 0
      }
    } else{
      if (is.na(auc)){
        auc <- 0
      }
    }
    
    if(method=="ave"){
      if(mean(c(auc.left.child, auc.right.child)) >= auc){
        good.split <- TRUE
      }
    } else if(method=="all"){
      if(auc.left.child >= auc & auc.right.child >= auc){
        good.split <- TRUE
        }
      } else if(method=="any"){
       if(auc.left.child >= auc | auc.right.child >= auc){
         good.split <- TRUE
       }
      } else if(method=="anyFurthur"){
         if(auc.left.child >= auc | auc.right.child >= auc | 
            (auc.left.child > furthurCutoff & auc.right.child > furthurCutoff)){
           good.split <- TRUE
         }
      }
    
    return(good.split)
  }
  
  is.good.split.hybrid <- function(d, method, hybridMethod){
    if(is.null(hybridMethod)){
      m <- c(method, method)
    } else{
      m <- hybridMethod
    }
    good.split <- FALSE
    if(is.leaf(d)){
      return(good.split)
    }
    good.split_global <- is.good.split(d, method=m[1], opt="global")
    good.split_pairwise <- is.good.split(d, method=m[2], opt="pairwise")
    if(good.split_global){
      good.split = TRUE
      return(good.split)
    } else{
      if(good.split_pairwise){
        good.split = TRUE
        return(good.split)
      }
    }
    return(good.split)
  }
  
  is.father.of.subtree.to.merge <- function(d, method, opt, hybridMethod){
    is.father <- FALSE
    if(opt=="global" | opt=="pairwise"){
      if(is.good.split(d, method=method, opt=opt) == FALSE & !is.leaf(d)) is.father <- TRUE
    } else if(opt=="hybrid"){
      if(is.good.split.hybrid(d, method=method, hybridMethod=hybridMethod) == FALSE & !is.leaf(d)) is.father <- TRUE
    } else{
      stop("invalid opt option: provide an option among global, pairwise and hybrid")
    }
    
    return(is.father)
  }
  
  search_subtree <- function(d){
    if(!is.father.of.subtree.to.merge(d, method=method, opt=opt, hybridMethod=hybridMethod)){
      for (i in seq_len(length(d))){
        if(is.leaf(d[[i]])){
          next()
        }
        d[[i]] <- search_subtree(d[[i]])
      }
    } else { # merge tree
      # achieve attributes
      branch.attributes <- attributes(d)
      # get all the leaf labels in this subtree
      leaf.labels <- labels(d)
      # prune subtree
      if(length(leaf.labels) > 1){
        dend <- prune(d, leaf.labels[-1])
        # assign attributes
        attributes(d) <- branch.attributes[-2]
        attr(d, "leaf") <- TRUE
        attr(d, "members") <- 1
        attr(d, "label") <- paste(leaf.labels[1], "m", sep="_")
      }
    }
    return(d)
  }
  new_dend <- search_subtree(d)
  # re-assign space
  new_dend <- ladderize(new_dend, right=FALSE)
  new_dend <- ladderize(fix_members_attr.dendrogram(new_dend), right=FALSE)
}

getleafs <- function (dend, renameClusterlabel = FALSE){
  if (renameClusterlabel) {
    dend <- assign_values_to_nodes(dend, "label", seq(1, 
                                                      nnodes(dend), 1))
    clusters <- get_leaves_attr(dend, "label")
  } else {
    clusters <- get_leaves_attr(dend, "label")
  }
  cells <- get_leaves_attr(dend, "nodesinfo")
  sizes <- get_leaves_attr(dend, "size")
  clusters <- unlist(lapply(1:length(sizes), function(r) {
    rep(clusters[r], sizes[r])
  }))
  names(clusters) <- cells
  return(clusters)
}

#' get overlap rate of reference marker genes
get_max_overlap_group <- function(delist, target.de, top=20, clname="NULL"){
  overlap <- 0
  overlaprate <- lapply(1:length(target.de), function(r){
    celltype <- names(target.de)[r]
    target.de.gene <- target.de[[r]][1:top, 1]
    de.gene <- rownames(delist[1:top, ])
    length(intersect(target.de.gene, de.gene))/top
  }) %>% unlist
  names(overlaprate) <- names(target.de)
  max.overlap <- max(overlaprate)
  max.cl <- names(overlaprate)[which(overlaprate==max.overlap)]
  if(max.overlap==0){
    data.frame(celltype="NULL", overlaprate=0, cluster=clname)
  } else{
    data.frame(celltype=max.cl, overlaprate=max.overlap, cluster=clname)
  }
}

#' get marker gene list from leaf node
#' @return a list contains marker genes for each leaf node 
getMarkers <- function(d, attribute="assGenes"){
  markerlist <- list()
  leafmarkers <- function(d){
    for(i in seq_len(length(d))){
      if(is.leaf(d[[i]])){
        markers <- attr(d[[i]], attribute)
        label <- attr(d[[i]], "label")
        markers$cluster <- label
        markerlist[[length(markerlist)+1]] <<- markers
        
      } else{
        d[[i]] <- leafmarkers(d[[i]])
      }
    }
    return(d)
  } 
  t <- leafmarkers(d)
  names(markerlist) <- dendextend::get_leaves_attr(d, "label")
  return(markerlist)
}

