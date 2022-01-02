scdi<-function(data1,data2)
{
	x1<-fftRtsne(data1)
	x2<-fftRtsne(data2)
	opt1<-NbClust(x1,method="kmeans")
	opt2<-NbClust(x2,method="kmeans")
    cluster1<-opt1$Best.partition
    cluster2<-opt2$Best.partition
    t1<-table(cluster1)
    t2<-table(cluster2)
    num1<-max(as.numeric(names(t1)))
    num2<-max(as.numeric(names(t2)))
    M1<-min(t1)
    M2<-min(t2)
    M<-min(M1,M2)   
    numit<-5
    nngraph<-array(0,dim<-c((dim(data1)[1]),(dim(data2)[1])))
    for(iter in 1:numit)
    {
    	dat1<-array(0,dim<-c((M*num1),(dim(data)[2])))
    	dat2<-array(0,dim<-c((M*num2),(dim(data)[2])))
    	innddexx1<-NULL
    	innddexx2<-NULL
    	for(i in 1:num1)
    	{
    		rand<-((which(cluster1==i))[sample(t1[i],M)])
    		innddexx1<-c(innddexx1,rand)
    	}
    	for(i in 1:num2)
    	{
    		rand<-((which(cluster2==i))[sample(t2[i],M)])
    		innddexx2<-c(innddexx2,rand)
    	}
    	
    		dat1<-data1[innddexx1,]
    		dat2<-data2[innddexx2,]
    		
    	datanew<-rbind(dat1,dat2)
    	groups<-c(rep(0,(dim(dat1)[1])),rep(1,(dim(dat2)[1])))
    	xxx<-fit.gplvm(datanew,2,iterations=100,Z.normal.prior=TRUE,groups=groups,sigma1_init=1)
        xxx1<-xxx$Z[(1:(dim(dat1)[1])),]
        xxx2<-xxx$Z[(((dim(dat1)[1])+1):(dim(dat1)[1]+dim(dat2)[1])),]
        xxx1<-sweep(xxx1,2,colMeans(xxx1))
        xxx2<-sweep(xxx2,2,colMeans(xxx2))
        
        Sigmatrix1<-array(0,dim<-c((2*num1),(2*num1)))
        Sigmatrix2<-array(0,dim<-c((2*num2),(2*num2)))
        for(i in 1:num1)
        {
        	Sigmatrix1[((2*(i-1)+1):(2*i)),((2*(i-1)+1):(2*i))]<-cov(xxx1[(((i-1)*M+1):(i*M)),])
        }
        for(i in 1:num2)
        {
        	Sigmatrix2[((2*(i-1)+1):(2*i)),((2*(i-1)+1):(2*i))]<-cov(xxx2[(((i-1)*M+1):(i*M)),])
        }
        nndist1<-array(0,dim<-c((dim(dat1)[1]),(dim(dat2)[1])))
        nndist2<-array(0,dim<-c((dim(dat2)[1]),(dim(dat1)[1])))
        for(i in 1:(dim(dat1)[1]))
        {
        	for(j in 1:(dim(dat2)[1]))
        	{
        		ind1<-cluster1[i]
        		ind2<-cluster2[j]
        		nndist1[i,j]<-mahalanobis(xxx1[i,],xxx2[j,],Sigmatrix1[((2*(ind1-1)+1):(2*ind1)),((2*(ind1-1)+1):(2*ind1))])
        		nndist2[j,i]<-mahalanobis(xxx1[i,],xxx2[j,],Sigmatrix2[((2*(ind2-1)+1):(2*ind2)),((2*(ind2-1)+1):(2*ind2))])
        	}
        }
        d1<-dim(nndist1)
        d2<-dim(nndist2)
        Z1<-array(0,dim<-c(d1[1],d2[2]))
        Z2<-array(0,dim<-c(d2[1],d1[2]))
        nnd1<-rbind(cbind(Z1,nndist1),cbind(nndist2,Z2))
        nn1<-graph_from_adjacency_matrix(nnd1,weighted=TRUE)
        knngraph1<-maxmatching(nn1,weighted=TRUE)
        gr<-array(0,dim<-c(length(knngraph1$matching),length(knngraph1$matching)))
        for(i in 1:(length(knngraph1$matching)))
        {
        	gr[i,(knngraph1$matching[i])]<-1
        }
        knn1<-gr[(1:(dim(nndist1)[1])),((dim(nndist2)[2]+1):(dim(nndist2)[2]+dim(nndist1)[2]))]
        knn2<-gr[((dim(nndist1)[1]+1):(dim(nndist1)[1]+dim(nndist2)[1])),(1:(dim(nndist2)[2]))]
        for(i in 1:(dim(dat1)[1]))
        {
        	for(j in 1:(dim(dat2)[1]))
        	{
        		nngraph[innddexx1[i],innddexx2[j]]<-nngraph[innddexx1[i],innddexx2[j]]+knn1[i,j]+knn2[j,i]
        	}
        }

    }
    return(nngraph)        
}



scdi_integrate<-function(data,groups)
{
library(igraph)
library(maxmatching)	

data1<-data[(groups==0),]
data2<-data[(groups==1),]



set.seed(0)
library(tsne)
source("fast_tsne.R")
source("kernels11_random.R")
source("GPLVM10_random.R")
library(NbClust)


xxxnew<-scdi(data1,data2)

colnames(xxxnew)<-(1:(dim(xxxnew)[2]))+(dim(xxxnew)[1])
rownames(xxxnew)<-1:(dim(xxxnew)[1])

library(igraph)

set.seed(123)

# generate random bipartite graph.
d<-dim(xxxnew)
Z1<-array(0,dim<-c(d[1],d[1]))
Z2<-array(0,dim<-c(d[2],d[2]))
A<-rbind(cbind(Z1,xxxnew),cbind(t(xxxnew),Z2))

g <- graph.adjacency(A,weighted=TRUE)
# check the type attribute:
V(g)$type

# define color and shape mappings.
col <- c("steelblue", "orange")
shape <- c("circle", "square")

plot(g)




left1<-(rowSums(xxxnew)==0)
left2<-(colSums(xxxnew)==0)

rem1<-(rowSums(xxxnew)>0)
rem2<-(colSums(xxxnew)>0)

nngraph1<-xxxnew[rem1,rem2]

d1<-dim(nngraph1)

Z1<-array(0,dim<-c(d1[1],d1[1]))
Z2<-array(0,dim<-c(d1[2],d1[2]))

nngraph2<-rbind(cbind(Z1,nngraph1),cbind(t(nngraph1),Z2))
colnames(nngraph2)<-c(rownames(nngraph1),colnames(nngraph1))
rownames(nngraph2)<-c(rownames(nngraph1),colnames(nngraph1))

nngraph3<-nngraph2+t(nngraph2)
#nngraph3[nngraph3>0]<-1

#library(statGraph)
#g<-list()
#g[[1]]<-nngraph3
#clust<-graph.cluster(g,3)

library(lsa)
g<-graph_from_adjacency_matrix(nngraph3)
g1<-as.undirected(g)
c1 = cluster_fast_greedy(g1,modularity=TRUE)

cluster<-membership(c1)

dX<-as.matrix(dist(data))


clust<-rep(NA,(dim(data)[1]))

clust[as.numeric(names(cluster))]<-cluster

mat1<-array(0,dim<-c(length(which(left1)),max(cluster)))
mat2<-array(0,dim<-c(length(which(left2)),max(cluster)))

for(j in 1:(dim(mat1)[2]))
{
	mat1[,j]<-rowMeans(dX[which(left1),na.omit(clust==j)])
}

for(j in 1:(dim(mat2)[2]))
{
	mat2[,j]<-rowMeans(dX[which(left2),na.omit(clust==j)])
}
if(dim(mat1)[1]>0)
{
for(i in 1:(dim(mat1)[1]))
{
	clust[as.numeric(names(which(left1)))[i]]<-which(mat1[i,]==min(mat1[i,]))
}
}

if(dim(mat2)[1]>0)
{
for(i in 1:(dim(mat2)[1]))
{
	clust[as.numeric(names(which(left2)))[i]]<-which(mat2[i,]==min(mat2[i,]))
}
}



source("kernels11_random.R")
source("GPLVM10_random.R")
xxx<-fit.gplvm(data,2,iterations=100,Z.normal.prior=TRUE,groups=groups,sigma1_init=1)

Z1<-xxx$Z
q<-dim(Z1)[2]
muu<-array(0,dim<-c(max(clust),q))
sigmaa<-array(0,dim<-c((q*max(clust)),q))
for(k in 1:(max(clust)))
{
	muu[k,]<-colMeans(Z1[(clust==k),])
	if((sum(((clust==k)&(groups==0)))>2)&(sum(((clust==k)&(groups==1)))>2))
	{
	sigmaa[(((k-1)*q+1):(k*q)),]<-(var(Z1[((clust==k)&(groups==0)),])*(sum((clust==k)&(groups==0))-1)+var(Z1[((clust==k)&(groups==1)),])*(sum((clust==k)&(groups==1))-1))/(sum(clust==k)-2)
	}else
	{
		sigmaa[(((k-1)*q+1):(k*q)),]<-var(Z1)
	}
}

Z1<-xxx$Z
Z2<-Z1
for(k in 1:(max(clust)))
{
m1<-colMeans(Z1[((groups==0)&(clust==k)),])-colMeans(Z1[((clust==k)),])
newfn<-function(x,vec=m1)
{
    return((x-vec))
}
Z2[((groups==0)&(clust==k)),]<-t(apply(Z1[((groups==0)&(clust==k)),],1,newfn))
}
for(k in 1:(max(clust)))
{
m2<-colMeans(Z1[((groups==1)&(clust==k)),])-colMeans(Z1[((clust==k)),])
newfn<-function(x,vec=m2)
{
    return((x-vec))
}
Z2[((groups==1)&(clust==k)),]<-t(apply(Z1[((groups==1)&(clust==k)),],1,newfn))
}


eps<-1

while(eps>(0.01*(dim(Z2)[1])))
{
	Z2old<-Z2
	source("GPLVM10_random_prior1.R")
	source("kernels11_random_prior.R")
	xxx1<-fit.gplvm.prior(data,2,iterations=100,Z.init=Z2,Z.normal.prior=TRUE,groups=groups,sigma1_init=1,muu=muu,sigmaa=sigmaa,cluster=clust)
	
	
	cluster<-clust
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
m1<-colMeans(Z1[((groups==0)&(clust==k)),])-colMeans(Z1[((clust==k)),])
newfn<-function(x,vec=m1)
{
    return((x-vec))
}
Z2[((groups==0)&(clust==k)),]<-t(apply(Z1[((groups==0)&(clust==k)),],1,newfn))
}
for(k in 1:(max(clust)))
{
m2<-colMeans(Z1[((groups==1)&(clust==k)),])-colMeans(Z1[((clust==k)),])
newfn<-function(x,vec=m2)
{
    return((x-vec))
}
Z2[((groups==1)&(clust==k)),]<-t(apply(Z1[((groups==1)&(clust==k)),],1,newfn))
}

eps<-sum(abs(Z2-Z2old)/(dim(Z)[1]))

print(eps)
plot(Z2,col=cols)
}

return(Z2)

}