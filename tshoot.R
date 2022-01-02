data1original<-rbind(array(rnorm(500000,-1,1),dim<-c(500,1000)),array(rnorm(500000,-1,1),dim<-c(500,1000)))
data2original<-rbind(array(rnorm(500000,-2,1),dim<-c(500,1000)),array(rnorm(500000,5,1),dim<-c(500,1000)))

data1original<-rbind(array(rnorm(50000,-1,1),dim<-c(50,1000)),array(rnorm(50000,-1,1),dim<-c(50,1000)))
data2original<-rbind(array(rnorm(50000,-1,1),dim<-c(50,1000)),array(rnorm(50000,-1,1),dim<-c(50,1000)))





rownames(data1original)<-1:(dim(data1original)[1])
rownames(data2original)<-1:(dim(data2original)[1])+1000

colnames(data1original)<-1:(dim(data1original)[2])
colnames(data2original)<-1:(dim(data2original)[2])

cols<-c(rep(1,500),rep(2,500),rep(1,500),rep(2,500))
groups<-c(rep(0,1000),rep(1,1000))

library(cccd)

gr<-nng(xxx)
library(igraph)
mmat<-as_adjacency_matrix(gr)



gr<-nng(ddi)
library(igraph)
mmat<-as_adjacency_matrix(gr)



