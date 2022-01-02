

#load("G1_S_bio.Rdata")

#data<-t(cbind(data11[,-1],data21[,-1]))



set.seed(0)
library(mvtnorm)
data1<-rmvnorm(30,rep(-3,300),diag(300))
data2<-rmvnorm(50,rep(3,300),diag(300))
data3<-rmvnorm(70,rep(6,300),diag(300))
data4<-rmvnorm(30,rep(-3,300),diag(300))
data5<-rmvnorm(50,rep(3,300),diag(300))
data6<-rmvnorm(70,rep(6,300),diag(300))
batch1<-rnorm(300,0,10)
batch2<-rnorm(300,0,10)
for(i in 1:(dim(data1)[1]))
{
	data1[i,]<-data1[i,]+batch1
}

for(i in 1:(dim(data2)[1]))
{
	data2[i,]<-data2[i,]+batch1
}

for(i in 1:(dim(data3)[1]))
{
	data3[i,]<-data3[i,]+batch1
}

for(i in 1:(dim(data4)[1]))
{
	data4[i,]<-data4[i,]+batch2
}

for(i in 1:(dim(data5)[1]))
{
	data5[i,]<-data5[i,]+batch2
}

for(i in 1:(dim(data6)[1]))
{
	data6[i,]<-data6[i,]+batch2
}




data<-rbind(data1,data2,data3,data4,data5,data6)



groups<-c(rep(0,150),rep(1,150))


N<-50
cols<-c(rep("red",N),rep("black",N),rep("blue",N),rep("green",N),rep("purple",N),rep("cyan",N))

cols<-c(rep("red",30),rep("black",N),rep("blue",70),rep("green",30),rep("purple",N),rep("cyan",70))

inndd<-sample((dim(data1)[2]),100)
#data<-rbind(data1[sample((dim(data1)[1]),150),inndd],data2[sample((dim(data2)[1]),150),inndd])




#cols<-c(rep("red",150),rep("black",150))

iterations=1000
plot.freq=100
classes=NULL
Z.init="PCA"
num.init.params=100
Z.normal.prior<-TRUE



library(GPLVM)
xxxold<-fit.gplvm(data,2,iterations=100,Z.normal.prior=TRUE)


cluster<-c(rep(1,30),rep(2,N),rep(3,70),rep(1,30),rep(2,N),rep(3,70))

source("kernels11_random.R")
source("GPLVM10_random.R")
xxx<-fit.gplvm(data,2,iterations=100,Z.normal.prior=TRUE,groups=groups,sigma1_init=1)




data1<-data[(groups==0),]
data2<-data[(groups==1),]



set.seed(0)
library(tsne)
source("fast_tsne.R")
source("kernels11_random.R")
source("GPLVM10_random.R")
library(NbClust)
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
    	dat1<-array(0,dim<-c((M*num1),100))
    	dat2<-array(0,dim<-c((M*num2),100))
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
    	
    		dat1<-data1[innddexx1,sample(dim(data1)[2],100)]
    		dat2<-data2[innddexx2,sample(dim(data1)[2],100)]
    		
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
        knngraph1<-array(0,dim<-c((dim(dat1)[1]),(dim(dat2)[1])))
        knngraph2<-array(0,dim<-c((dim(dat2)[1]),(dim(dat1)[1])))
        for(i in 1:(dim(dat1)[1]))
        {
        	for(j in 1:(dim(dat2)[1]))
        	{
        		indd1<-which(nndist1[i,]==min(nndist1[i,]))
        		indd2<-which(nndist2[j,]==min(nndist2[j,]))
        		knngraph1[i,indd1]<-1
        		knngraph2[j,indd2]<-1
        	}
        }
        for(i in 1:(dim(dat1)[1]))
        {
        	for(j in 1:(dim(dat2)[1]))
        	{
        		nngraph[innddexx1[i],innddexx2[j]]<-nngraph[innddexx1[i],innddexx2[j]]+knngraph1[i,j]+knngraph2[i,j]
        	}
        }

    }
    return(nngraph)        
}



xxxnew<-scdi(data1,data2)


library(igraph)

set.seed(123)

# generate random bipartite graph.
g <- graph.adjacency(xxxnew)
# check the type attribute:
V(g)$type

# define color and shape mappings.
col <- c("steelblue", "orange")
shape <- c("circle", "square")

plot(g)












