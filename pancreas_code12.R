a1<-read.csv("GSM2230757_human1_umifm_counts.csv")
a2<-read.csv("GSM2230758_human2_umifm_counts.csv")

info1<-a1[,1:3]
info2<-a2[,1:3]

a1<-a1[,-(1:3)]
a2<-a2[,-(1:3)]

set.seed(12345)

set1<-sample(1937,100)
set2<-sample(1724,100)

a3<-a1[set1,]
a4<-a2[set2,]

data1<-as.matrix(a3)
data2<-as.matrix(a4)

colnames(data1)<-colnames(a3)
colnames(data2)<-colnames(a4)


data3<-array(0,dim<-dim(data1))
data4<-array(0,dim<-dim(data2))

for(i in 1:(dim(data1)[1]))
{
	data3[i,]<-as.numeric(data1[i,])
}

for(i in 1:(dim(data2)[1]))
{
	data4[i,]<-as.numeric(data2[i,])
}

for(i in 1:(dim(data3)[1]))
{
	data3[i,]<-data3[i,]/sum(data3[i,])*1000000
}

for(i in 1:(dim(data4)[1]))
{
	data4[i,]<-data4[i,]/sum(data4[i,])*1000000
}

a5<-a1[(sample(1937,1)),]
a5<-as.numeric(a5)
a5<-a5/sum(a5)*1000000


colnames(data3)<-colnames(a3)
colnames(data4)<-colnames(a4)

data1<-data3
data2<-data4

for(i in 1:(dim(data1)[1]))
{
	data1[i,]<-data1[i,]/sqrt(sum((data1[i,])^2))
}
for(i in 1:(dim(data2)[1]))
{
	data2[i,]<-data2[i,]/sqrt(sum((data2[i,])^2))
}


data<-rbind(data1,data2)

groups<-c(rep(0,100),rep(1,100))


save.image("pamcreas1.Rdata")

source("scdi10.R")
K1<-7
K2<-7


Z2<-scdi_integrate(data,groups,K1,K2)

save.image("pancreas2.Rdata")

set3<-setdiff(1:1937,set1)
set4<-setdiff(1:1724,set2)


l<-Z2$l
alpha<-Z2$alpha
sigma<-Z2$sigma
sigma1<-Z2$sigma1
Z3<-Z2$Z2
groups<-Z2$groups
cluster<-Z2$cluster
dataa<-Z2$dataa

library(GPLVM)
X1<-t(chol(dataa[(groups==0),]%*%t(dataa[(groups==0),])))
X2<-t(chol(dataa[(groups==1),]%*%t(dataa[(groups==1),])))
X3<-forwardsolve(X1,diag(length(diag(X1))))
X4<-forwardsolve(X2,diag(length(diag(X2))))

source("gplvm_SE.R")
K<-gplvm.SE(Z3[(groups==0),],l,alpha,sigma)

Kinv<-solve(K)
Kvecmat<-Kinv%*%dataa[(groups==0),]

library(emulator)

Z2Z2<-array(0,dim<-c(length(set3),dim(Z3)[2]))

inndex<-1

dataa1<-dataa[(groups==0),]
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

clustnew1<-NULL

for(i in set3)
{
	a5<-a1[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa1,a5)
	corvec<-cov(a5,t(dataa1))*nrow(data1)
	disclust<-NULL
	for(j in 1:(length(unique(cl1))))
	{
		disclust[j]<-mean(corvec[cl1==j])
	}
	kk<-which(disclust==max(disclust))
	corvec1<-corvec[cl1==kk]
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z4<-Z3[(groups==0),]
    Z5<-Z4[(cl1==kk),]
    y<-colMeans(Z5)
    opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    Z2Z2[inndex,]<-opt$par
    clustnew1[inndex]<-kk
    inndex<-inndex+1	
    print(inndex)
}



dataa2<-a1[set3,]

dataa2<-as.matrix(dataa2)
norm<-sqrt(rowSums(dataa2^2))

for(i in 1:(dim(dataa2)[1]))
{
	dataa2[i,]<-dataa2[i,]/norm[i]
}

Z3Z3<-array(0,dim<-c(length(set1),dim(Z3)[2]))

inndex<-1

clusterold1<-NULL

for(i in set1)
{
	a5<-a1[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa2,a5)
	corvec<-cov(a5,t(dataa2))*nrow(data1)
	disclust<-NULL
	for(j in 1:(length(unique(clustnew1))))
	{
		disclust[j]<-mean(corvec[clustnew1==j])
	}
	kk<-which(disclust==max(disclust))
	corvec1<-corvec[clustnew1==kk]
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z4<-Z2Z2
    Z5<-Z4[(clustnew1==kk),]
    y<-colMeans(Z5)
    opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    Z3Z3[inndex,]<-opt$par
    clusterold1[inndex]<-kk
    inndex<-inndex+1	
    print(inndex)
}



























Z4Z4<-array(0,dim<-c(length(set4),dim(Z3)[2]))

inndex<-1

dataa3<-dataa[(groups==1),]
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

clustnew2<-NULL

for(i in set4)
{
	a5<-a2[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa3,a5)
	corvec<-cor(a5,t(data3))*nrow(data2)
	disclust<-NULL
	for(j in 1:(length(unique(cl1))))
	{
		disclust[j]<-mean(corvec[cl1==j])
	}
	kk<-which(disclust==max(disclust))
	corvec1<-corvec[cl1==kk]
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z4<-Z3[(groups==1),]
    Z5<-Z4[(cl1==kk),]
    y<-colMeans(Z5)
    opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    Z4Z4[inndex,]<-opt$par
    clustnew2[inndex]<-kk
    inndex<-inndex+1	
    print(inndex)
}




dataa4<-a1[set4,]

dataa4<-as.matrix(dataa4)
norm<-sqrt(rowSums(dataa4^2))

for(i in 1:(dim(dataa4)[1]))
{
	dataa4[i,]<-dataa4[i,]/norm[i]
}


Z5Z5<-array(0,dim<-c(length(set2),dim(Z3)[2]))

inndex<-1

clusterold2<-NULL

for(i in set2)
{
	a5<-a2[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa4,a5)
	corvec<-cor(a5,t(dataa4))*nrow(data2)
	disclust<-NULL
	for(j in 1:(length(unique(clustnew2))))
	{
		disclust[j]<-mean(corvec[clustnew2==j])
	}
	kk<-which(disclust==max(disclust))
	corvec1<-corvec[clustnew2==kk]
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z4<-Z4Z4
    Z5<-Z4[(clustnew2==kk),]
    y<-colMeans(Z5)
    opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    Z5Z5[inndex,]<-opt$par
    clusterold2[inndex]<-kk
    inndex<-inndex+1	
    print(inndex)
}








