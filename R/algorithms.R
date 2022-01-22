#' Run recursive leiden clustering algorithm - modified from conos (https://github.com/hms-dbmi/conos)
#'
#' @param graph igraph object
#' @param K number of iterations
#' @param renameCluster if rename clusters
#' @export
rleiden.detection <- function(graph, K=2, renameCluter=TRUE, n.cores=1,
                              min.community.size=10, verbose=FALSE, resolution=1, K.current=1, ...){

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
                          min.community.size=min.community.size, verbose=verbose, n.cores=n.cores)
      },mc.cores=n.cores,mc.allow.recursive = FALSE)
    } else {
      wtl <- lapply(conos:::sn(unique(mem)), function(cluster) {
        cn <- names(mem)[which(mem==cluster)]
        sg <- induced.subgraph(graph,cn)
        rleiden.detection(induced.subgraph(graph,cn), K=K, resolution=resolution, K.current=K.current+1,
                          min.community.size=min.community.size, verbose=verbose, n.cores=n.cores)
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

#' Evaluate expression distances between clusters
#' @param counts raw count matrix - rows are cells and columns are genes
#' @export
cluster.matrix.expression.distances <- function(counts, groups=NULL, useVariablegenes=F, variableGenes=NULL, dist='cor',
                                                use.single.cell.comparisons=FALSE, use.scaled.data=FALSE, min.cluster.size=1, max.n.cells=Inf, n.cells=200) {
  require(abind)
  require(sccore)
  if(is.null(groups)) {
    stop('no groups specified')
  } else {
    groups <- as.factor(groups)
  }

  valid.dists <- c('JS','cor');
  if(!dist %in% valid.dists) stop(paste('only the following distance types are supported:',paste(valid.dists,collapse=', ')))

  cl <- factor(groups[match(rownames(counts), names(groups))],levels=levels(groups))
  tt <- table(cl)
  empty <- tt< min.cluster.size

  # aggregate clusters
  if(any(tt>max.n.cells)) { # need subsampling
    scn <- unlist(tapply(names(cl), cl, function(x) sample(x,min(max.n.cells,length(x)))))
    cl[!(names(cl) %in% scn)] <- NA; tt <- table(cl);
  }

  if(useVariablegenes){
    if(is.null(variableGenes)){
      stop("no variable genesets provided")
    } else{
      counts <- counts[, colnames(counts) %in% variableGenes]
    }
  }

  if(use.single.cell.comparisons) { # sample individual cells and compare

    tcd <- parallel::mclapply(1:n.cells, function(i) {
      scn <- unlist(tapply(names(cl), cl, function(x) sample(x,1)))
      tc <- as.matrix(counts[na.omit(as.character(scn)), ])
      rownames(tc) <- names(scn)[!is.na(scn)]
      tc <- tc[match(levels(cl), rownames(tc)),]
      rownames(tc) <- levels(cl)
      if(dist=='JS') {
        if(use.scaled.data){
          tc <- t(tc)
          tcd <- pagoda2:::jsDist(tc); dimnames(tcd) <- list(colnames(tc),colnames(tc));
        } else{
          tc <- t(tc/pmax(1,rowSums(tc)))
          tcd <- pagoda2:::jsDist(tc); dimnames(tcd) <- list(colnames(tc),colnames(tc));
        }
      } else { # correlation distance
        if(use.scaled.data){
          tcd <- 1-cor(t(tc))
        } else{
          tc <- log10(t(tc/pmax(1,rowSums(tc)))*1e3+1)
          tcd <- 1-cor(tc)
        }
      }
      tcd[empty,] <- tcd[,empty] <- NA;
      tcd
    }, mc.cores=10) %>% abind(along=3) %>% apply(c(1,2),median,na.rm=T)

  } else{
    tc <- sccore:::colSumByFactor(counts, cl);
    tc <- tc[-1,,drop=F]  # omit NA cells
    # correlation distance
    if(dist=='JS') {
      tc <- t(tc/pmax(1,rowSums(tc)))
      tcd <- pagoda2:::jsDist(tc); dimnames(tcd) <- list(colnames(tc),colnames(tc));
    } else { # correlation distance
      if(use.scaled.data){
        tc <- t(tc)
      } else{
        tc <- log10(t(tc/pmax(1,rowSums(tc)))*1e3+1)
      }
      tcd <- 1-cor(tc)
    }
  }
  return(tcd)
}

#'
estimateFoldChanges <- function(cm, markers, annotation) {
  annotation %<>% .[intersect(names(.), rownames(cm))]
  cm %<>% .[names(annotation), markers]
  mean.mtx <- split(1:length(annotation), annotation) %>%
    sapply(function(ids) Matrix::colMeans(cm[ids,]))
  
  if(ncol(mean.mtx) > 2){
    fold.change.mtx <- 1:ncol(mean.mtx) %>%
      sapply(function(i) log2(1e-10 + mean.mtx[,i]) - log2(1e-10 + rowMeans(mean.mtx[,-i]))) %>%
      `colnames<-`(colnames(mean.mtx))
  } else{
    fold.change.mtx <- 1:ncol(mean.mtx) %>%
      sapply(function(i) log2(1e-10 + mean.mtx[,i]) - log2(1e-10 + mean.mtx[,-i])) %>%
      `colnames<-`(colnames(mean.mtx))
  }
  
  return(fold.change.mtx)
}


#' caculate model performance
modelPerformance <- function(pred, pred_prob, class, cat){
  
  confusion <- as.matrix(table(class, pred, deparse.level = 0))
  
  if(cat > 2){
    n <- sum(confusion) # number of instances
    nc <- nrow(confusion) # number of classes
    diag <- diag(confusion) # number of correctly classified instances per class
    accuracy <- sum(diag)/sum(confusion)
    rowsums <- apply(confusion, 1, sum) # number of instances per class
    colsums <- apply(confusion, 2, sum) # number of predictions per class
    p <- rowsums / n # distribution of instances over the actual classes
    q <- colsums / n # distribution of instances over the predicted classes
    precision <- diag / colsums
    recall <- diag / rowsums
    f1 <- 2 * precision * recall / (precision + recall)
    macroPrecision <- mean(precision, na.rm = T)
    macroRecall <- mean(recall, na.rm = T)
    macroF1 <- mean(f1, na.rm = T)
    return(data.frame(accuracy = accuracy, macroPrecision = macroPrecision,
                      macrorecall = macroRecall, macrof1 = macroF1))
  } else {
    accuracy <- (confusion[1,1] + confusion[2,2]) / sum(confusion)
    precision <- confusion[1,1] / (confusion[1,1] + confusion[1,2])
    recall <- confusion[1,1] / (confusion[1,1] + confusion[2,1])
    f1 <- 2 * precision * recall / (precision + recall)
    auc <- pROC::auc(pROC::roc(class, pred_prob))[1]
    return(data.frame(accuracy = accuracy, Precision = precision,
                      recall = recall, f1 = f1, auc = auc))
  }
}



#' train classifier based on genesets
#' @param feature.matrix provided matrix for training, rows are genes, columns are cells
trainClassifier <- function(markers, annotation, feature.matrix, model="rf", subsample=T, ncells=200,
                            n.core=10){
  require(doParallel)
  
  feature.matrix <- feature.matrix[markers, names(annotation)]
  split <- train_test_split(t(feature.matrix), trainSize=0.8)
  train <- split$train
  test <- split$test
  
  train.label <- annotation[match(rownames(train), names(annotation))]
  test.label <- annotation[match(rownames(test), names(annotation))]
  
  if(subsample){
    train.index <- sample(1:nrow(train), ncells, replace = T)
    test.index <- sample(1:nrow(test), ncells, replace = T)
    
    train <- train[train.index, ]
    test <- test[test.index, ]
    train.label <- train.label[train.index]
    test.label <- test.label[test.index]
  }
  
  if(model == "rf"){
    if(n.core>1){
      cl <- makePSOCKcluster(n.core)
      registerDoParallel(cl)
      
      rf <- caret::train(as.matrix(train), as.factor(train.label), method="rf",
                       trControl = caret::trainControl(method = "cv",number = 3))
      
      stopCluster(cl)
    } else{
      rf <- caret::train(as.matrix(train), as.factor(train.label), method="rf",
                         trControl = caret::trainControl(method = "cv",number = 3))
    }
    
    pred <- predict(rf, test)
    pred_prob <- predict(rf, test, "prob")
    
    perf.res <- modelPerformance(pred, pred_prob[, 1], test.label, cat=2)
  }
  
  return(list(model=rf, perf.res=perf.res))
}

#' function to split train-test sets
train_test_split <- function(data, trainSize=0.8){
  train.size <- floor(nrow(data) * trainSize)
  train.index <- sample(nrow(data), train.size, replace=F)
  
  train <- data[train.index, ]
  test <- data[-train.index, ]
  
  return(list(train=train, test=test))
}

#' implement rank-based normalization
#' 
