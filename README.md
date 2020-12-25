# speciesTree
Tree based method for cross-species data integration

speciesTree is initially designed for cross-species data analysis to identify homologous groups and evolutional patterns. 
The approach is used in [BICCN](https://biccn.org/) cross species analysis to identify robust homologous groups in human, mouse and marmoset. 
The method can be extended to the other analysis as well, such as identifying roubust clusters in single and multi-modal datasets.


## Pipeline
![overview](https://github.com/huqiwen0313/speciesTree/blob/master/figs/pipeline.overview.png)

## Installation

### Requirements
* [Conos](https://github.com/hms-dbmi/conos)
* [Pagoda2](https://github.com/hms-dbmi/pagoda2)

```r
install.packages("devtools")
devtools::install_github("huqiwen0313/speciesTree")
library(speciesTree)
```

## Quick start
* build robust homologous clusters across species
```r
library(speciesTree)
library(conos)
library(pagoda2)


# read example files
# gene expression profile
expression.matrix <- readRDS("./example_data/integrated_matrix_all_cells.RDS")
# mete cell clusters
meta_clusters <- readRDS("./example_data/all_cells_cluster_membership.RDS")
# read cell annotations
upperlevelinfo <- readRDS("./example_data/upperlevelinfo.rds")
species_annot <- sapply(strsplit(names(meta_clusters), "_"), `[`, 1)
names(species_annot) <- names(meta_clusters)
# read subsampled cell clusters
subsampledClusters <- readRDS("./example_data/subsampledClusters.rds")

# build tree
d <- cluster.matrix.expression.distances(t(expression.matrix), groups=meta_clusters, dist="cor", 
                                         useVariablegenes=FALSE,  use.scaled.data=TRUE)
dendr <- hclust(as.dist(d), method='ward.D2')
dend <- as.dendrogram(dendr)
tree <- buildSpeciesTree(dend=dend, expMatrix=t(expression.matrix), subsampleClusters=subsampledClusters,
                         cls.groups=meta_clusters, upperlevelannot=upperlevelinfo, species=species_annot,
                         n.cores=10)

# prune tree
tree$dend <- TreeEntropy(tree$dend, entropy.cutoff = 2.9)
dend_pruned <- pruneTreeEntropy(tree$dend, cutoff=2.9)

# get homologous clusters
homo.clusters <- getClusters(dend_pruned, plotTree=TRUE)
```

## Citations




