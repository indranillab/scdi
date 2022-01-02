X1<-array(rnorm(200,-2,1),dim<-c(100,2))
X2<-array(rnorm(200,2,1),dim<-c(100,2))
X<-rbind(X1,X2)




graph_cluster<-function(X,nclust)
{
	ddd<-dist(X)
	mat<-as.matrix(ddd)
	gr1<-mst(mat)
	gr2<-array(as.numeric(gr1),dim=dim(mat))
	mat1<-mat*gr2
	library(igraph)
	gr<-graph_from_adjacency_matrix(mat1,weighted=TRUE)
	i<-1
	while(numb<nclust)
	{
		ind<-which(mat1==max(mat1),arr.ind=TRUE)
		mat1[ind[1,1],ind[1,2]]<-0
		mat1[ind[1,2],ind[1,1]]<-0
		gr<-graph_from_adjacency_matrix(mat1,weighted=TRUE)
		vec<-component(gr)$membership
		numb<-sum(table(gr)$membership>(0.01*dim(X)[1]))
		i<-i+1
		print(i)
	}
	cluster<-components(gr)$membership
	return(cluster)
}






