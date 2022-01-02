scdi<-function(data1,data2,K1,K2)
{
	x1<-fftRtsne(data1)
	x2<-fftRtsne(data2)
	#opt1<-NbClust(x1,method="kmeans")
	#opt2<-NbClust(x2,method="kmeans")
    cluster1<-kmeans(x1,K1)$cluster
    cluster2<-kmeans(x2,K2)$cluster
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
    	#xxx<-fit.gplvm1(datanew,2,iterations=100,Z.normal.prior=TRUE,groups=groups,sigma1_init=1)
    	xxx<-fftRtsne(datanew,2,perplexity=10)
        #xxx1<-xxx$Z[(1:(dim(dat1)[1])),]
        #xxx2<-xxx$Z[(((dim(dat1)[1])+1):(dim(dat1)[1]+dim(dat2)[1])),]
        #xxx1<-sweep(xxx1,2,colMeans(xxx1))
        #xxx2<-sweep(xxx2,2,colMeans(xxx2))
        xxx1<-xxx[(1:(dim(dat1)[1])),]
        xxx2<-xxx[(((dim(dat1)[1])+1):(dim(dat1)[1]+dim(dat2)[1])),]
        
        K<-2
        Sigmatrix1<-array(0,dim<-c((K*num1),(K*num1)))
        Sigmatrix2<-array(0,dim<-c((K*num2),(K*num2)))
        for(i in 1:num1)
        {
        	Sigmatrix1[((K*(i-1)+1):(K*i)),((K*(i-1)+1):(K*i))]<-cov(xxx1[(((i-1)*M+1):(i*M)),])
        }
        for(i in 1:num2)
        {
        	Sigmatrix2[((K*(i-1)+1):(K*i)),((K*(i-1)+1):(K*i))]<-cov(xxx2[(((i-1)*M+1):(i*M)),])
        }
        nndist1<-array(0,dim<-c((dim(dat1)[1]),(dim(dat2)[1])))
        nndist2<-array(0,dim<-c((dim(dat2)[1]),(dim(dat1)[1])))
        for(i in 1:(dim(dat1)[1]))
        {
        	for(j in 1:(dim(dat2)[1]))
        	{
        		ind1<-cluster1[i]
        		ind2<-cluster2[j]
        		nndist1[i,j]<-mahalanobis((xxx1[i,]-colMeans(xxx1[(cluster1[innddexx1]==ind1),])),(xxx2[j,]-colMeans(xxx2[(cluster2[innddexx2]==ind2),])),Sigmatrix1[((K*(ind1-1)+1):(K*ind1)),((K*(ind1-1)+1):(K*ind1))])
        		nndist2[j,i]<-mahalanobis((xxx1[i,]-colMeans(xxx1[(cluster1[innddexx1]==ind1),])),(xxx2[j,]-colMeans(xxx2[(cluster2[innddexx2]==ind2),])),Sigmatrix2[((K*(ind2-1)+1):(K*ind2)),((K*(ind2-1)+1):(K*ind2))])
        	}
        }
        d1<-dim(nndist1)
        d2<-dim(nndist2)
        Z1<-array(0,dim<-c(d1[1],d2[2]))
        Z2<-array(0,dim<-c(d2[1],d1[2]))
        nnd1<-rbind(cbind(Z1,(2*max(nndist1)-nndist1)),cbind((2*max(nndist2)-nndist2),Z2))
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
    return(list(nngraph=nngraph,cluster1=cluster1,cluster2=cluster2))        
}



scdi_integrate<-function(data,groups,K1,K2)
{
library(igraph)
library(maxmatching)	

library(GPLVM)


X<-data
X<-t(chol(X%*%t(X)))

xxx<-fit.gplvm(X,max(K1,K2),iterations=100,Z.normal.prior=TRUE)
data<-xxx$Z
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

colnames(xxxnew)<-(1:(dim(xxxnew)[2]))+(dim(xxxnew)[1])
rownames(xxxnew)<-1:(dim(xxxnew)[1])

library(igraph)

set.seed(123)

# generate random bipartite graph.
d<-dim(xxxnew)
Z1<-array(0,dim<-c(d[1],d[1]))
Z2<-array(0,dim<-c(d[2],d[2]))
A<-rbind(cbind(Z1,xxxnew),cbind(t(xxxnew),Z2))
rownames(A)<-1:200
colnames(A)<-1:200

g <- graph.adjacency(A,weighted=TRUE)
# check the type attribute:
V(g)$type

# define color and shape mappings.
col <- c("steelblue", "orange")
shape <- c("circle", "square")

plot(g)


common_member<-array(0,dim<-c(length(table(cluster1)),length(table(cluster2))))
dc<-dim(common_member)
for(i in 1:(dc[1]))
{
	for(j in 1:(dc[2]))
	{
		common_member[i,j]<-sum(xxxnew[(cluster1==i),(cluster2==j)])
	}
}


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
#c1 = cluster_fast_greedy(g1,modularity=TRUE)

#cluster<-membership(c1)

cluster<-rep(NA,sum(dim(xxxnew)))
k<-1
for(i in 1:(dc[1]))
{
	ord<-order(common_member[i,])
	for(j in ord)
	{
		#if((common_member[i,j]/(sum(common_member[i,]))>(1/dc[2]))&(common_member[i,j]/(sum(common_member[,j]))>(1/dc[1])))
		if((common_member[i,j]==max(common_member[i,]))&(common_member[i,j]==max(common_member[,j])))
		{
			cluster[c((cluster1==i),(cluster2==j))]<-k
			k<-k+1
			break
		}
	}
}
clustern<-cluster[(1:(length(cluster1)))]
for(i in 1:(dc[1]))
{
	if(is.na(sum(clustern[cluster1==i])))
	{
		cluster[c((cluster1==i),rep(FALSE,length(cluster2)))]<-k
		k<-k+1
	}
}
clustern<-cluster[-(1:(length(cluster1)))]
for(j in 1:(dc[2]))
{
	if(is.na(sum(clustern[cluster2==j])))
	{
		cluster[c(rep(FALSE,length(cluster1)),(cluster2==j))]<-k
		k<-k+1
	}
}

names(cluster)<-1:(dim(data)[1])



dX<-as.matrix(dist(data))


clust<-rep(NA,(dim(data)[1]))

clust[as.numeric(names(cluster))]<-cluster

mat1<-array(0,dim<-c(length(which(left1)),max(cluster)))
mat2<-array(0,dim<-c(length(which(left2)),max(cluster)))
if(dim(mat1)[1]>0)
{
for(j in 1:(dim(mat1)[2]))
{
	if((dim(mat1)[1])>1)
	{
	mat1[,j]<-rowMeans(dX[which(left1),na.omit(clust==j)])
	}else
	{
	mat1[,j]<-mean(dX[which(left1),na.omit(clust==j)])
	}
}
}

if(dim(mat2)[1]>0)
{
for(j in 1:(dim(mat2)[2]))
{
	if((dim(mat2)[1])>1)
	{
	mat2[,j]<-rowMeans(dX[which(left2),na.omit(clust==j)])
	}else
	{
	mat2[,j]<-mean(dX[which(left2),na.omit(clust==j)])
	}
}
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
#xxx<-fit.gplvm1(X[,(order(cluster,groups))],2,iterations=100,Z.normal.prior=TRUE,groups=groups[order(cluster,groups)],cluster=cluster[order(cluster,groups)],sigma1_init=0)
xxx1<-fit.gplvm1(X,max(K1,K2),iterations=100,Z.normal.prior=TRUE,groups=groups,sigma1_init=0)


Z1<-xxx1$Z
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


#clustold<-clust
#clust<-kmeans(Z2,max(clustold))$cluster
#cluster<-clust

eps<-1

iter<-1
while(eps>(0.001*(dim(Z2)[1])))
{
	Z2old<-Z2
	source("GPLVM10_random_prior1.R")
	source("kernels11_random_prior.R")
	xxx1<-fit.gplvm.prior(X,max(K1,K2),iterations=100,Z.init=Z2,Z.normal.prior=TRUE,groups=groups,sigma1_init=1,muu=muu,sigmaa=sigmaa,cluster=clust)
	
	
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

eps<-sum(abs(Z2-Z2old)/(dim(Z2)[1]))

print(eps)
plot(Z2,col=cols)
if(iter>5)
{
	break
}
iter<-iter+1
}

return(Z2)

}