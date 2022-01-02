
graph_cluster<-function(X,nclust)
{
	library(igraph)
	if("ape"%in%(.packages()))
	{
		detach("package:ape")
	}
	#library(ape)
	ddd<-dist(X)
	mat<-as.matrix(ddd)
	mat1<-graph_from_adjacency_matrix(mat,weighted=TRUE)
	gr1<-minimum.spanning.tree(mat1)
	gr2<-as_adjacency_matrix(gr1)
	mat2<-mat*as.matrix(gr2)
	gr<-graph_from_adjacency_matrix((mat2),weighted=TRUE)
	vec<-components(gr)$membership
	numb<-sum(table(vec)>(0.01*dim(X)[1]))
	i<-1
	while(numb<nclust)
	{
		ind<-which(mat2==max(mat2),arr.ind=TRUE)
		mat2[ind[1,1],ind[1,2]]<-0
		mat2[ind[1,2],ind[1,1]]<-0
		gr<-graph_from_adjacency_matrix(mat2,weighted=TRUE)
		vec<-components(gr)$membership
		numb<-sum(table(vec)>(0.01*dim(X)[1]))
		i<-i+1
		print(i)
	}
	cluster<-components(gr)$membership
	return(cluster)
}


scdi<-function(data1,data2,clus1,clus2)
{
    cl1<-clus1
    cl2<-clus2
    dists<-array(0,dim<-c(length(table(cl1)),length(table(cl2))))
    library(pdist)
    tt1<-unique(cl1)
    tt2<-unique(cl2)
for(i in unique(cl1))
{
	for(j in unique(cl2))
	{
		dat1<-data1[(which(cl1==tt1[i])),]
		dat2<-data2[(which(cl2==tt2[j])),]
		dists[i,j]<-mean(as.matrix(pdist(dat1,dat2)))	
	}
	print(i)
}
nngraph<-array(0,dim<-c(length(table(cl1)),length(table(cl2))))
for(i in 1:(dim(dists)[1]))
{
	x<-dists[i,]
	if(length(x)>2)
	{
	stat<-(mean(x[-which(x==min(x))])-x)/sd(x[-which(x==min(x))])
    if(max(stat)>2)
	{
	 	j<-which(x==min(x))
		nngraph[i,j]<-1
	}else
	{
		j<-which(x==min(x))
		y<-dists[,j]
		if(y[i]==min(y))
		{
			nngraph[i,j]<-1
		}
		
	}
	}else
	{
		j<-which(x==min(x))
		nngraph[i,j]<-1
	}
}
    return(list(nngraph=nngraph,cluster1=cl1,cluster2=cl2))        
}



scdi_integrate<-function(data1,data2,clus1,clus2,ddd1,ddd2)
{

	
for(i in 1:(dim(data1)[1]))
{
	data1[i,]<-(data1[i,]-mean(data1[i,]))/sd(data1[i,])
}
for(i in 1:(dim(data2)[1]))
{
	data2[i,]<-(data2[i,]-mean(data2[i,]))/sd(data2[i,])
}

for(i in 1:(dim(data1)[1]))
{
	data1[i,]<-data1[i,]/sqrt(sum((data1[i,])^2))
}

for(i in 1:(dim(data2)[1]))
{
	data2[i,]<-data2[i,]/sqrt(sum((data2[i,])^2))
}

data<-rbind(data1,data2)



groups<-c(rep(0,(dim(data1)[1])),rep(1,dim(data2)[1]))


library(igraph)
library(maxmatching)	
library(GPLVM)


data[(groups==0),]<-sweep(data[(groups==0),],2,colMeans(data[(groups==0),]))
data[(groups==1),]<-sweep(data[(groups==1),],2,colMeans(data[(groups==1),]))
data3<-data
for(i in 1:(dim(data)[1]))
{
	data3[i,]<-data[i,]/sqrt(sum((data[i,])^2))
}


data1<-data3[(groups==0),]
data2<-data3[(groups==1),]



set.seed(0)
library(tsne)
library(GPLVM)
source("fast_tsne.R")
source("kernels11_random.R")
source("GPLVM10_random.R")
library(NbClust)
library(maxmatching)


xxxnew1<-scdi(data1,data2,clus1,clus2)

xxxnew<-xxxnew1$nngraph
cluster1<-xxxnew1$cluster1
cluster2<-xxxnew1$cluster2+(dim(xxxnew)[1])

set.seed(12345)
set1<-sample(1:(dim(data1)[1]),100)
set2<-sample(1:(dim(data2)[1]),100)
set3<-setdiff(1:(dim(data1)[1]),set1)
set4<-setdiff(1:(dim(data2)[1]),set2)

subdata1<-data1[set1,]
subdata2<-data2[set2,]

dataa<-rbind(subdata1,subdata2)

data<-rbind(subdata1,subdata2)
groups<-c(rep(0,length(set1)),rep(1,length(set2)))

X<-data
if((dim(data)[2])>(dim(data)[1]))
{
X<-t(chol(data%*%t(data)))
}


library(igraph)
colnames(xxxnew)<-(1:(dim(xxxnew)[2]))+(dim(xxxnew)[1])
rownames(xxxnew)<-1:(dim(xxxnew)[1])

A1<-array(0,dim<-c(dim(xxxnew)[1],dim(xxxnew)[1]))
A2<-array(0,dim<-c(dim(xxxnew)[2],dim(xxxnew)[2]))
colnames(A1)<-rownames(xxxnew)
rownames(A1)<-rownames(xxxnew)
colnames(A2)<-colnames(xxxnew)
rownames(A2)<-colnames(xxxnew)
B<-rbind(cbind(t(A1),xxxnew),cbind(t(xxxnew),A2))

comp<-components(graph_from_adjacency_matrix(B))

clustint<-c(cluster1[set1],cluster2[set2])

cluster<-rep(NA,length(clustint))

for(i in 1:(length(comp$membership)))
{
	cluster[clustint==i]<-comp$membership[i]
}

clust<-cluster

library(GPLVM)
source("kernels11_random.R")
source("GPLVM10_random.R")
#xxx<-fit.gplvm1(X[,(order(cluster,groups))],2,iterations=100,Z.normal.prior=TRUE,groups=groups[order(cluster,groups)],cluster=cluster[order(cluster,groups)],sigma1_init=0)
xxx1<-fit.gplvm1(X,4*max(K1,K2),iterations=100,Z.normal.prior=TRUE,Z.init="PCA",groups=groups,sigma1_init=0)


Z1<-xxx1$Z
q<-dim(Z1)[2]
muu<-array(0,dim<-c(length(unique(clust)),q))
sigmaa<-array(0,dim<-c((q*length(unique(clust))),q))
kk<-1
for(k in names(table(clust)))
{
	k<-as.numeric(k)
	if(table(clust)[k]==1)
	{
	    muu[k,]<-Z1[(clust==k),]	
	    sigmaa[(((k-1)*q+1):(k*q)),]<-0
	    next
	}
	muu[kk,]<-colMeans(Z1[(clust==k),])
	if((sum(((clust==k)&(groups==0)))>2)&(sum(((clust==k)&(groups==1)))>2))
	{
	sigmaa[(((k-1)*q+1):(k*q)),]<-(var(Z1[((clust==k)&(groups==0)),])*(sum((clust==k)&(groups==0))-1)+var(Z1[((clust==k)&(groups==1)),])*(sum((clust==k)&(groups==1))-1))/(sum(clust==k)-2)
	}else
	{
		sigmaa[(((kk-1)*q+1):(kk*q)),]<-var(Z1)
	}
	kk<-kk+1
}



Z1<-xxx1$Z
Z2<-Z1
kk<-1
for(k in (names(table(clust))))
{
	k<-as.numeric(k)
if(sum(((groups==0)&(clust==k)))>1)
{
	if(table(clust)[k]==1)
	{
	    m1<-Z1[((groups==0)&(clust==k)),]-Z1[((clust==k)),]
	}else
	{
m1<-colMeans(Z1[((groups==0)&(clust==k)),])-colMeans(Z1[((clust==k)),])
}
newfn<-function(x,vec=m1)
{
    return((x-vec))
}
if(sum(((groups==0)&(clust==k)))>1)
{
Z2[((groups==0)&(clust==k)),]<-t(apply(Z1[((groups==0)&(clust==k)),],1,newfn))
}
}

if(sum(((groups==1)&(clust==k)))>1)
{
	if(table(clust)[k]==1)
	{
	    m2<-Z1[((groups==1)&(clust==k)),]-Z1[((clust==k)),]
	}else
	{
m2<-colMeans(Z1[((groups==1)&(clust==k)),])-colMeans(Z1[((clust==k)),])
}
newfn<-function(x,vec=m2)
{
    return((x-vec))
}
if(sum(((groups==1)&(clust==k)))>1)
{
Z2[((groups==1)&(clust==k)),]<-t(apply(Z1[((groups==1)&(clust==k)),],1,newfn))
}
kk<-kk+1
}
}

Z1<-Z2
for(k in (names(table(clust))))
{
	k<-as.numeric(k)
	v<-var(Z1[(clust==k),])
if(sum(((groups==0)&(clust==k)))>1)
{
	if(table(clust)[k]==1)
	{
	    m1<-Z1[((groups==0)&(clust==k)),]-Z1[((clust==k)),]
	}else
	{
v1<-var(Z1[((groups==0)&(clust==k)),])
m1<-colMeans(Z1[((groups==0)&(clust==k)),])
}
ei<-eigen(v)
ei1<-eigen(v1)
mat<-(ei$vectors[,(ei$values>0)]%*%diag((sqrt(ei$values[ei$values>0])))%*%t(ei$vectors[,(ei$values>0)]))%*%(ei1$vectors[,(ei1$values>0)]%*%diag((1/sqrt(ei1$values[ei1$values>0])))%*%t(ei1$vectors[,(ei1$values>0)]))
newfn<-function(x,sig=mat,vec=m1)
{
    return(vec+mat%*%(x-vec))
}
if(sum(((groups==0)&(clust==k)))>1)
{
Z2[((groups==0)&(clust==k)),]<-t(apply(Z1[((groups==0)&(clust==k)),],1,newfn))
}
}

if(sum(((groups==1)&(clust==k)))>1)
{
	if(table(clust)[k]==1)
	{
	    m2<-Z1[((groups==1)&(clust==k)),]-Z1[((clust==k)),]
	}else
	{
v2<-var(Z1[((groups==1)&(clust==k)),])
m2<-colMeans(Z1[((groups==1)&(clust==k)),])
}
ei<-eigen(v)
ei2<-eigen(v2)
mat<-(ei$vectors[,(ei$values>0)]%*%diag((sqrt(ei$values[ei$values>0])))%*%t(ei$vectors[,(ei$values>0)]))%*%(ei2$vectors[,(ei2$values>0)]%*%diag((1/sqrt(ei2$values[ei2$values>0])))%*%t(ei2$vectors[,(ei2$values>0)]))
newfn<-function(x,sig=mat,vec=m2)
{
    return(vec+mat%*%(x-vec))
}
if(sum(((groups==1)&(clust==k)))>1)
{
Z2[((groups==1)&(clust==k)),]<-t(apply(Z1[((groups==1)&(clust==k)),],1,newfn))
}
}
}
#clust<-graph_cluster(zzz2,K1)
#cluster<-clust
removed<-NULL

ta<-table(clust)
i1<-1
for(i in 1:(length(table(clust))))
{
	if(ta[i]<5)
	{
		removed<-c(removed,(which(clust==i1)))
		dataa<-dataa[-(which(clust==i1)),]
		X<-X[-(which(clust==i1)),]
		Z2<-Z2[-(which(clust==i1)),]
		muu<-muu[-i1,]
		sigmaa<-sigmaa[-(((i1-1)*q+1):(i1*q)),]
		groups<-groups[-(which(clust==i1))]
		clust<-clust[-which(clust==i1)]
		clust[clust>i1]<-clust[clust>i1]-1
		i1<-i1-1
	}
	i1<-i1+1
}



eps<-1

iter<-1
	Z2old<-Z2
	source("GPLVM10_random_prior1.R")
	source("kernels11_random_prior.R")
	xxx1<-fit.gplvm.prior(X,4*max(K1,K2),iterations=100,Z.init=Z2,Z.normal.prior=FALSE,groups=groups,sigma1_init=1,muu=muu,sigmaa=sigmaa,cluster=clust)
	

Z1<-xxx1$Z
q<-dim(Z1)[2]
muu<-array(0,dim<-c(length(unique(clust)),q))
sigmaa<-array(0,dim<-c((q*length(unique(clust))),q))
kk<-1
for(k in names(table(clust)))
{
	k<-as.numeric(k)
	if(table(clust)[kk]==1)
	{
	    muu[k,]<-Z1[(clust==k),]	
	    sigmaa[(((k-1)*q+1):(k*q)),]<-0
	    next
	}
	muu[kk,]<-colMeans(Z1[(clust==k),])
	if((sum(((clust==k)&(groups==0)))>2)&(sum(((clust==k)&(groups==1)))>2))
	{
	sigmaa[(((k-1)*q+1):(k*q)),]<-(var(Z1[((clust==k)&(groups==0)),])*(sum((clust==k)&(groups==0))-1)+var(Z1[((clust==k)&(groups==1)),])*(sum((clust==k)&(groups==1))-1))/(sum(clust==k)-2)
	}else
	{
		sigmaa[(((kk-1)*q+1):(kk*q)),]<-var(Z1)
	}
	kk<-kk+1
}


	
    #clust<-graph_cluster(zzz2,K1)
	#cluster<-clust
Z1<-xxx1$Z
Z2<-Z1
kk<-1
for(k in (names(table(clust))))
{
	k<-as.numeric(k)
if(sum(((groups==0)&(clust==k)))>1)
{
	if(table(clust)[k]==1)
	{
	    m1<-Z1[((groups==0)&(clust==k)),]-Z1[((clust==k)),]
	}else
	{
m1<-colMeans(Z1[((groups==0)&(clust==k)),])-colMeans(Z1[((clust==k)),])
}
newfn<-function(x,vec=m1)
{
    return((x-vec))
}
if(sum(((groups==0)&(clust==k)))>1)
{
Z2[((groups==0)&(clust==k)),]<-t(apply(Z1[((groups==0)&(clust==k)),],1,newfn))
}
}

if(sum(((groups==1)&(clust==k)))>1)
{
	if(table(clust)[k]==1)
	{
	    m2<-Z1[((groups==1)&(clust==k)),]-Z1[((clust==k)),]
	}else
	{
m2<-colMeans(Z1[((groups==1)&(clust==k)),])-colMeans(Z1[((clust==k)),])
}
newfn<-function(x,vec=m2)
{
    return((x-vec))
}
if(sum(((groups==1)&(clust==k)))>1)
{
Z2[((groups==1)&(clust==k)),]<-t(apply(Z1[((groups==1)&(clust==k)),],1,newfn))
}
kk<-kk+1
}
}

Z1<-Z2
for(k in (names(table(clust))))
{
	k<-as.numeric(k)
	v<-var(Z1[(clust==k),])
if(sum(((groups==0)&(clust==k)))>1)
{
	if(table(clust)[k]==1)
	{
	    m1<-Z1[((groups==0)&(clust==k)),]-Z1[((clust==k)),]
	}else
	{
v1<-var(Z1[((groups==0)&(clust==k)),])
m1<-colMeans(Z1[((groups==0)&(clust==k)),])
}
ei<-eigen(v)
ei1<-eigen(v1)
mat<-(ei$vectors[,(ei$values>0)]%*%diag((sqrt(ei$values[ei$values>0])))%*%t(ei$vectors[,(ei$values>0)]))%*%(ei1$vectors[,(ei1$values>0)]%*%diag((1/sqrt(ei1$values[ei1$values>0])))%*%t(ei1$vectors[,(ei1$values>0)]))
newfn<-function(x,sig=mat,vec=m1)
{
    return(vec+mat%*%(x-vec))
}
if(sum(((groups==0)&(clust==k)))>1)
{
Z2[((groups==0)&(clust==k)),]<-t(apply(Z1[((groups==0)&(clust==k)),],1,newfn))
}
}

if(sum(((groups==1)&(clust==k)))>1)
{
	if(table(clust)[k]==1)
	{
	    m2<-Z1[((groups==1)&(clust==k)),]-Z1[((clust==k)),]
	}else
	{
v2<-var(Z1[((groups==1)&(clust==k)),])
m2<-colMeans(Z1[((groups==1)&(clust==k)),])
}
ei<-eigen(v)
ei2<-eigen(v2)
mat<-(ei$vectors[,(ei$values>0)]%*%diag((sqrt(ei$values[ei$values>0])))%*%t(ei$vectors[,(ei$values>0)]))%*%(ei2$vectors[,(ei2$values>0)]%*%diag((1/sqrt(ei2$values[ei2$values>0])))%*%t(ei2$vectors[,(ei2$values>0)]))
newfn<-function(x,sig=mat,vec=m2)
{
    return(vec+mat%*%(x-vec))
}
if(sum(((groups==1)&(clust==k)))>1)
{
Z2[((groups==1)&(clust==k)),]<-t(apply(Z1[((groups==1)&(clust==k)),],1,newfn))
}
}
}

#clust<-graph_cluster(zzz2,K1)
#cluster<-clust

source("kernels11_random.R")
source("GPLVM10_random.R")
#xxx<-fit.gplvm1(X[,(order(cluster,groups))],2,iterations=100,Z.normal.prior=TRUE,groups=groups[order(cluster,groups)],cluster=cluster[order(cluster,groups)],sigma1_init=0)
xxx1<-fit.gplvm1(Z2,2*max(K1,K2),iterations=100,Z.normal.prior=TRUE,Z.init="PCA",groups=groups,sigma1_init=0)
X2<-Z2
Z2<-xxx1$Z




return(list(Z2=Z2,cluster=clust,groups=groups,l=xxx1$l,alpha=xxx1$alpha,sigma=xxx1$sigma,sigma1=xxx1$sigma1,X2=X2,set1=set1,set2=set2,set3=set3,set4=set4,dataa=dataa))


}





