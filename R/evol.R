#' Compute features of each gene
#' @param count count matrix - rows are cells, columns are genes
#' @param cells cells intested
#' @param geneRanks data.frame of geneRanks
#' @param prefix factor(species) groups
#' @param sampleFac sample cells from each factor(species) so that the cell sizes are the same
#' @return data frame contain several features for each gene
addFeatures <- function(count, cells, group= NULL, sampleCell=TRUE, sampleSize=1000, prefix=c("human", "marmoset", "mouse"), minSize=10,
                        n.cores=30, normalize=FALSE, depthScale=1e3, sizeCutoff=50, sampleFac=FALSE){
  if(length(prefix)!=3){
    stop("currently only support species (factors)")
  }
  cellsOut <- rownames(count)[-which(rownames(count) %in% cells)]
  if(sampleCell){
    if(is.null(group)){
      stop("please provide group annotation")
    }
    cells <- sampleCells(cells, group, sampleSize=sampleSize)
    cellsOut <- sampleCells(cellsOut, group, sampleSize=sampleSize)
  }
  group1Cells <- cells[grep(prefix[1], cells)]
  group2Cells <- cells[grep(prefix[2], cells)]
  group3Cells <- cells[grep(prefix[3], cells)]
  
  if(sampleFac){
    if(min(c(length(group1Cells),length(group2Cells), length(group3Cells))) > 0){
      cellSize <- min(c(length(group1Cells),length(group2Cells), length(group3Cells)))
      if(cellSize > sizeCutoff){
        group1Cells <- sample(group1Cells, cellSize, replace=F)
        group2Cells <- sample(group2Cells, cellSize, replace=F)
        group3Cells <- sample(group3Cells, cellSize, replace=F)
      } else{
        group1Cells <- sample(group1Cells, sizeCutoff, replace=T)
        group2Cells <- sample(group2Cells, sizeCutoff, replace=T)
        group3Cells <- sample(group3Cells, sizeCutoff, replace=T)
      }
    } else{
      sizes <- c(length(group1Cells),length(group2Cells), length(group3Cells))
      cellSize <- min(sizes[sizes>0])
      groups <- list(group1Cells, group2Cells, group3Cells)
      pos <- which(sizes>0)
      if(cellSize > sizeCutoff){
        if(length(group1Cells)==0){
          group2Cells <- sample(group2Cells, cellSize, replace=F)
          group3Cells <- sample(group3Cells, cellSize, replace=F)
        } else if(length(group2Cells)==0){
          group1Cells <- sample(group1Cells, cellSize, replace=F)
          group3Cells <- sample(group3Cells, cellSize, replace=F)
        } else if(length(group3Cells)==0){
          group1Cells <- sample(group1Cells, cellSize, replace=F)
          group2Cells <- sample(group2Cells, cellSize, replace=F)
        }
      } else{
        if(length(group1Cells)==0){
          group2Cells <- sample(group2Cells, cellSize, replace=T)
          group3Cells <- sample(group3Cells, cellSize, replace=T)
        } else if(length(group2Cells)==0){
          group1Cells <- sample(group1Cells, cellSize, replace=T)
          group3Cells <- sample(group3Cells, cellSize, replace=T)
        } else if(length(group3Cells)==0){
          group1Cells <- sample(group1Cells, cellSize, replace=T)
          group2Cells <- sample(group2Cells, cellSize, replace=T)
        }
      }
    }
  }
  
  group1OutCells <- cellsOut[grep(prefix[1], cellsOut)]
  group2OutCells <- cellsOut[grep(prefix[2], cellsOut)]
  group3OutCells <- cellsOut[grep(prefix[3], cellsOut)]
  cellannot <- setNames(c(rep(prefix[1], length(group1Cells)), rep(prefix[2], length(group2Cells)),
                          rep(prefix[3], length(group3Cells))), c(group1Cells, group2Cells, group3Cells))
  cellannotOut <- setNames(c(rep(prefix[1], length(group1OutCells)), rep(prefix[2], length(group2OutCells)),
                             rep(prefix[3], length(group3OutCells))), c(group1OutCells, group2OutCells, group3OutCells))

  # gene expression value for each factor(species)
  if(normalize){
    depth <- Matrix::colSums(count)
    count <- count/as.numeric(depth/depthScale)
    count@x <- as.numeric(log(count@x+1))
    expByGroup <- conos:::collapseCellsByType(count, cellannot) %>% apply(2, function(r){
      r/sum(r)
    }) %>% t
    #expByGroup <- expByGroup[!is.na(expByGroup[, 1]), ]
    colnames(expByGroup) <- paste(names(table(cellannot)), c("exp"), sep="_")
  } else{
    expByGroup <- conos:::collapseCellsByType(count, cellannot) %>% apply(2, function(r){
      r/sum(r)
    }) %>% t
  }
  
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
  return(geneFeatures)
}

#' Binormial Proportion test
prop_test <- function(count_1, count_2){
  # proportion of cells expressed in group1
  p1 <- length(count_1[count_1 > 0])/length(count_1)
  # proportion of cells expressed in group2
  p2 <- length(count_2[count_2 > 0])/length(count_2)

  p <-  prop.test(c(round(length(count_1)*p1, 0), round(length(count_2)*p2, 0)),
                  c(length(count_1), length(count_2)), alternative="two.sided", correct=F)
  return(p)
}

#' proportion test for each individual genes
#' @param count gene count matrix - rows are cells, columns are tree
#' @return prop test with rank and p-values
#' @export
propTestGene <- function(count, group1.cells, group2.cells, n.cores=10, sampleCells=TRUE,
                         sizeCutoff=50){
  genelist <- colnames(count)
  if(sampleCells){
    cellSize <- min(c(length(group1.cells),length(group2.cells)))
    if(cellSize > sizeCutoff){
        group1.cells <- sample(group1.cells, cellSize, replace=F)
        group2.cells <- sample(group2.cells, cellSize, replace=F)
    } else{
        group1.cells <- sample(group1.cells, sizeCutoff, replace=T)
        group2.cells <- sample(group2.cells, sizeCutoff, replace=T)
      }
    } 
  
  gene.test <- parallel::mclapply(1:length(genelist), function(r){
    gene_count <- count[, genelist[r]]
    group1Count <- gene_count[names(gene_count) %in% group1.cells]
    group2Count <- gene_count[names(gene_count) %in% group2.cells]
    p <- prop_test(group1Count, group2Count)
    data.frame(gene=genelist[r], stat=p$statistic, p.value=p$p.value)
  }, mc.cores=n.cores)
  gene.test <- dplyr::bind_rows(gene.test)
  gene.test$p.adjust <- p.adjust(gene.test$p.value, method="BH", n = length(gene.test$p.value))
  ### rank genes
  gene.test <- gene.test[order(gene.test$p.value), ]
  gene.test$rank <- seq(1, nrow(gene.test), 1)
  return(gene.test)
}


#' sample cells based on given cell composition
#' @param cells total cells to be sampled
#' @param groups factor contains cell annotation
#' @return sampled cells
sampleCells <- function(cells, groups, sampleSize=1000){
  cells <- cells[cells %in% names(groups)]
  #sampleSize <- quantile((table(groups)), 0.5)
  cellsannot <- groups[names(groups) %in% cells]
  nodeCellComp <- table(cellsannot)/sum(table(cellsannot))
  celltypeSampleSize <- round(nodeCellComp * sampleSize, 0)

  cellsType <- lapply(names(table(cellsannot)), function(r){
    names(cellsannot[cellsannot == r])
  })
  names(cellsType) <- names(table(cellsannot))

  sampledCells <- lapply(1:length(names(cellsType)), function(r){
    size <- celltypeSampleSize[names(celltypeSampleSize) == names(cellsType)[r]]
    if(size > length(cellsType[[r]])){
      sample(cellsType[[r]], size, replace=TRUE)
    } else{
      sample(cellsType[[r]], size)
    }
  }) %>% unlist
  return(sampledCells)
}


