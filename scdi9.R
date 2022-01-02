
graph_cluster<-function(X,nclust)
{
	library(igraph)
	if("ape"%in%(.packages()))
	{
		detach("package:ape")
	}
	library(ape)
	ddd<-dist(X)
	mat<-as.matrix(ddd)
	gr1<-mst(mat)
	gr2<-array(as.numeric(gr1),dim=dim(mat))
	mat1<-mat*gr2
	gr<-graph_from_adjacency_matrix(mat1,weighted=TRUE)
	vec<-components(gr)$membership
	numb<-sum(table(vec)>(0.01*dim(X)[1]))
	i<-1
	while(numb<nclust)
	{
		ind<-which(mat1==max(mat1),arr.ind=TRUE)
		mat1[ind[1,1],ind[1,2]]<-0
		mat1[ind[1,2],ind[1,1]]<-0
		gr<-graph_from_adjacency_matrix(mat1,weighted=TRUE)
		vec<-components(gr)$membership
		numb<-sum(table(vec)>(0.01*dim(X)[1]))
		i<-i+1
		print(i)
	}
	cluster<-components(gr)$membership
	return(cluster)
}


scdi<-function(data1,data2,K1,K2)
{
	x1<-fftRtsne(data1)
	x2<-fftRtsne(data2)
	#opt1<-NbClust(x1,method="kmeans")
	#opt2<-NbClust(x2,method="kmeans")
	library(SC3)
	cluster1<-graph_cluster(x1,K1)
	cluster2<-graph_cluster(x2,K2)
    #cluster1<-as.numeric(as.vector(unlist(colData(xxx1))))
    #cluster2<-as.numeric(as.vector(unlist(colData(xxx2))))
    t1<-table(cluster1)
    t2<-table(cluster2)
    cl1<-graph_cluster(x1,7)
    cl2<-graph_cluster(x2,7)
    dists<-array(0,dim<-c(length(table(cl1)),length(table(cl2))))
    library(pdist)
for(i in 1:length(table(cl1)))
{
	for(j in 1:length(table(cl2)))
	{
		dat1<-data1[(which(cl1==i)),]
		dat2<-data2[(which(cl2==j)),]
		dists[i,j]<-mean(as.matrix(pdist(dat1,dat2)))	
	}
	print(i)
}
nngraph<-array(0,dim<-c(length(table(cl1)),length(table(cl2))))
for(i in 1:(dim(dists)[1]))
{
	x<-dists[i,]
	stat<-(mean(x)-x)/sd(x)
	if(max(stat)>((length(x)-1)/sqrt(length(x))*sqrt(((qt(0.05/length(x),(length(x)-2)))^2/((qt(0.05/length(x),(length(x)-2)))^2+(length(x)-2))))))
	{
		j<-which(x==min(x))
		nngraph[i,j]<-1
	}
}

    return(list(nngraph=nngraph,cluster1=cluster1,cluster2=cluster2))        
}



scdi_integrate<-function(data,groups,K1,K2)
{
	

data1<-data[(groups==0),]
data2<-data[(groups==1),]



library(igraph)
library(maxmatching)	

library(GPLVM)


X<-data
X<-t(chol(data%*%t(data)))

#xxx<-fit.gplvm(X,max(K1,K2),iterations=100,Z.normal.prior=TRUE)
#data<-xxx$Z
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


xxxnew1<-scdi(data1,data2,K1,K2)


cluster1<-xxxnew1$cluster1
cluster2<-xxxnew1$cluster2
xxxnew<-xxxnew1$nngraph

colnames(xxxnew)<-(1:(dim(xxxnew)[2]))
rownames(xxxnew)<-1:(dim(xxxnew)[1])


####### Done till this point 

cluster<-rep(NA,sum(dim(xxxnew)))
clust1<-rep(NA,(dim(xxxnew)[1]))
clust2<-rep(NA,(dim(xxxnew)[2]))
k<-1
dc<-dim(xxxnew)
for(i in 1:(dc[1]))
{
    if(sum(xxxnew[i,])==0)
    {
         clust1[cluster1==i]<-k
         k<-k+1
    }else
    {
    	     j<-which(xxxnew[i,]>0)
    	     clust1[cluster1==i]<-k
    	     clust2[cluster2==j]<-k
    	     k<-k+1
    }
}

for(j in 1:(dc[2]))
{
	if(sum(xxxnew[,j])==0)
	{
		clust2[cluster2==j]<-k
		k<-k+1
	}
}

cluster<-c(clust1,clust2)

names(cluster)<-1:(dim(data)[1])

clust<-cluster

ta<-table(cluster)

for(kk in 1:(max(clust)))
{
	if(ta[kk]<5)
	{
	X<-X[-(which(cluster==kk)),]
	groups<-groups[-which(cluster==kk)]
	cluster<-cluster[-which(cluster==kk)]
	}
}

clust<-cluster


source("kernels11_random.R")
source("GPLVM10_random.R")
#xxx<-fit.gplvm1(X[,(order(cluster,groups))],2,iterations=100,Z.normal.prior=TRUE,groups=groups[order(cluster,groups)],cluster=cluster[order(cluster,groups)],sigma1_init=0)
xxx1<-fit.gplvm1(X,max(K1,K2),iterations=100,Z.normal.prior=TRUE,groups=groups,sigma1_init=0)


Z1<-xxx1$Z
q<-dim(Z1)[2]
muu<-array(0,dim<-c(max(clust),q))
sigmaa<-array(0,dim<-c((q*max(clust)),q))
for(k in 1:(max(clust)))
{
	if(table(clust)[k]==1)
	{
	    muu[k,]<-Z1[(clust==k),]	
	    sigmaa[(((k-1)*q+1):(k*q)),]<-0
	    next
	}
	muu[k,]<-colMeans(Z1[(clust==k),])
	if((sum(((clust==k)&(groups==0)))>2)&(sum(((clust==k)&(groups==1)))>2))
	{
	sigmaa[(((k-1)*q+1):(k*q)),]<-(var(Z1[((clust==k)&(groups==0)),])*(sum((clust==k)&(groups==0))-1)+var(Z1[((clust==k)&(groups==1)),])*(sum((clust==k)&(groups==1))-1))/(sum(clust==k)-2)
	}else
	{
		sigmaa[(((k-1)*q+1):(k*q)),]<-var(Z1)
	}
}



Z1<-xxx1$Z
Z2<-Z1
if(sum(((groups==0)&(clust==k)))>1)
{
for(k in 1:(max(clust)))
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
}
if(sum(((groups==1)&(clust==k)))>1)
{
for(k in 1:(max(clust)))
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
}
}


zzz2<-fftRtsne(Z2)

#clust<-graph_cluster(zzz2,K3)
#cluster<-clust

eps<-1

iter<-1
while(eps>(0.001*(dim(Z2)[1])))
{
	Z2old<-Z2
	source("GPLVM10_random_prior1.R")
	source("kernels11_random_prior.R")
	xxx1<-fit.gplvm.prior(X,max(K1,K2),iterations=100,Z.init=Z2,Z.normal.prior=TRUE,groups=groups,sigma1_init=1,muu=muu,sigmaa=sigmaa,cluster=clust)
	
    #clust<-graph_cluster(zzz2,K3)
	#cluster<-clust
	Z1<-xxx1$Z
q<-dim(Z1)[2]
muu<-array(0,dim<-c(max(cluster),q))
sigmaa<-array(0,dim<-c((q*max(cluster)),q))
for(k in 1:(max(cluster)))
{
	muu[k,]<-colMeans(Z1[(cluster==k),])
	if((sum(((cluster==k)&(groups==0)))>2)&(sum(((cluster==k)&(groups==1)))>2))
	{
	sigmaa[(((k-1)*q+1):(k*q)),]<-(var(Z1[((cluster==k)&(groups==0)),])*(sum((cluster==k)&(groups==0))-1)+var(Z1[((cluster==k)&(groups==1)),])*(sum((cluster==k)&(groups==1))-1))/(sum(cluster==k)-2)
	}else
	{
		sigmaa[(((k-1)*q+1):(k*q)),]<-var(Z1)
	}
}

Z1<-xxx1$Z
Z2<-Z1
for(k in 1:(max(clust)))
{
if(sum(clust==k)>1)
{
m1<-colMeans(Z1[((groups==0)&(clust==k)),])-colMeans(Z1[((clust==k)),])
newfn<-function(x,vec=m1)
{
    return((x-vec))
}
Z2[((groups==0)&(clust==k)),]<-t(apply(Z1[((groups==0)&(clust==k)),],1,newfn))
}
}
for(k in 1:(max(clust)))
{
if(sum(clust==k)>1)
{
m2<-colMeans(Z1[((groups==1)&(clust==k)),])-colMeans(Z1[((clust==k)),])
newfn<-function(x,vec=m2)
{
    return((x-vec))
}
Z2[((groups==1)&(clust==k)),]<-t(apply(Z1[((groups==1)&(clust==k)),],1,newfn))
}
}


zzz2<-fftRtsne(Z2)

#clust<-graph_cluster(zzz2,K3)
#cluster<-clust

eps<-sum(abs(Z2-Z2old)/(dim(Z2)[1]))

print(eps)
#plot(Z2,col=cols)
if(iter>5)
{
	break
}
iter<-iter+1
}



return(Z2)

}





