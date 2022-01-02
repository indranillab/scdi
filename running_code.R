library(SMNN)
library(Seurat)
library(dplyr)
library(Matrix)
library(reshape2)
source("SMNNcorrect.R")
source("unifiedClusterLabelling.R")

# Load the example data data_SMNN
data("data_SMNN")

# Provide the marker genes for cluster matching
markers <- c("Col1a1", "Pdgfra", "Ptprc", "Pecam1")

# Specify the cluster labels for each marker gene
cluster.info <- c(1, 1, 2, 3)


datatype="count"
ident_list=NULL
cluster.use=NULL
cluster.names=NULL
min.exp.thresh=0
min.perc=0.3
min.average.expr=0

features.use<-markers
cluster.labels<-cluster.info

batches<-list(data_SMNN$batch1.mat, data_SMNN$batch2.mat)



