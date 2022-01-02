
rm(list=ls())
load("pancreas3.Rdata")
rm(list=setdiff(ls(),c("set1","set2")))


load("pan1.Rdata")

xyz1<-which(!((clust1%in%num1)))
xyz2<-which(!((clust2%in%num2)))
wxyz1<-which(((clust1%in%num1)))
wxyz2<-which(((clust2%in%num2)))
wxyz1<-setdiff(wxyz1,set1)
wxyz2<-setdiff(wxyz2,set2)

rm(list=setdiff(ls(),c("xyz1","xyz2","wxyz1","wxyz2","set1","set2")))




load("pancreas4.Rdata")

xyz3<-wxyz1[!match1]
xyz4<-wxyz2[!match2]

inn1<-c(xyz1,xyz3)
inn2<-c(xyz2,xyz4)



l<-Z2$l
alpha<-Z2$alpha
sigma<-Z2$sigma
sigma1<-Z2$sigma1
Z3<-Z2$X2
groups<-Z2$groups
cluster<-Z2$cluster
dataa<-Z2$dataa
X3<-Z2$X2

library(emulator)

Z2Z2<-array(0,dim<-c(length(inn1),dim(X3)[2]))

inndex<-1

dataa1<-dataa
source("add_gplvm.R")

cost<-function(y,d,X1)
{
	s<-0
	for(i in 1:(dim(X1)[1]))
	{
		s<-s+(d[i]-sum((y-X1[i,])^2))^2
	}
	return(s)
}

for(i in 1:(dim(dataa1)[1]))
{
	dataa1[i,]<-dataa1[i,]/sqrt(sum((dataa1[i,])^2))
}

cl1<-Z2$cluster[(Z2$groups==0)]

clustnew11<-NULL
clist<-unique(cl1)

cormat1<-array(0,dim<-c((dim(dataa1)[1]),length(c(inn1))))

for(i in inn1)
{
	a5<-a1[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa1,a5)
	corvec<-cor(a5,t(dataa1))
	disclust<-NULL
	for(j in 1:(length(unique(cl1))))
	{
		disclust[j]<-mean(corvec[cl1==clist[j]])
	}
	kk<-which(disclust==max(disclust))
	corvec1<-corvec[cl1==clist[kk]]
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z4<-X3[(groups==0),]
    Z5<-Z4[(cl1==clist[kk]),]
    if(sum(dim(Z5))>0)
    {
    y<-colMeans(Z5)
    }
    else
    {
    	y<-Z5
    }
    #opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    #Z2Z2[inndex,]<-opt$par
    clustnew11[inndex]<-clist[kk]
    cormat1[,inndex]<-1-corvec
    inndex<-inndex+1	
    print(inndex)
}




Z4Z4<-array(0,dim<-c(length(inn2),dim(X3)[2]))

inndex<-1

dataa3<-dataa
source("add_gplvm.R")

cost<-function(y,d,X1)
{
	s<-0
	for(i in 1:(dim(X1)[1]))
	{
		s<-s+(d[i]-sum((y-X1[i,])^2))^2
	}
	return(s)
}

for(i in 1:(dim(dataa3)[1]))
{
	dataa3[i,]<-dataa3[i,]/sqrt(sum((dataa3[i,])^2))
}

cl1<-Z2$cluster[(Z2$groups==1)]

clustnew12<-NULL
clist<-unique(cl1)

cormat2<-array(0,dim<-c((dim(dataa3)[1]),length(inn2)))
for(i in inn2)
{
	a5<-a2[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa3,a5)
	corvec<-cor(a5,t(dataa3))
	disclust<-NULL
	for(j in 1:(length(unique(cl1))))
	{
		disclust[j]<-mean(corvec[cl1==clist[j]])
	}
	kk<-which(disclust==max(disclust))
	corvec1<-corvec[cl1==clist[kk]]
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z4<-X3[(groups==1),]
    Z5<-Z4[(cl1==clist[kk]),]
    if(sum(dim(Z5))>0)
    {
    y<-colMeans(Z5)
    }
    else
    {
    	y<-Z5
    }
    #opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    #Z4Z4[inndex,]<-opt$par
    clustnew12[inndex]<-clist[kk]
    cormat2[,inndex]<-1-corvec
    inndex<-inndex+1	
    print(inndex)
}


set.seed(12345)

clus1<-Z2$cluster[(Z2$groups==0)]
clus2<-Z2$cluster[(Z2$groups==1)]

clus<-c(clus1,clus2)
Zmat2<-array(0,dim<-c((length(inn1)+length(inn2)),(dim(X3)[2])))
inndex<-1
clustn2<-c(clustnew11,clustnew12)

a<-as.matrix(rbind(a1,a2)[c(inn1,inn2),])
kk<-3
library(mvtnorm)
for(j in 1:(dim(Zmat2)[1]))
{
    Z5<-dataa[(clus==clustn2[inndex]),]
    dis<-NULL
    for(i in 1:(dim(Z5)[1]))
    {
    	    dis[i]<-sum((a[j,]-Z5[i,])^2)
    }
    ss<-sort(dis,index=TRUE)
	mat<-X3[(clus==clustn2[inndex]),][ss$ix[1:kk],]
    Zmat2[inndex,]<-colMeans(mat)+rmvnorm(1,rep(0,(dim(Zmat2)[2])),diag(diag(var(mat))))
    inndex<-inndex+1	
    print(inndex)
}


save.image("pancreas6.Rdata")


