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
library(tilting)
ei1<-svd(t(dataa[(groups==0),]))
#mat1<-(ei1$u)%*%diag(ei1$d)%*%t(ei1$u)
ei2<-svd(t(dataa[(groups==1),]))
#mat2<-(ei2$u)%*%diag(ei2$d)%*%t(ei2$u)


source("gplvm_SE.R")
K<-gplvm.SE(Z3[(groups==0),],l,alpha,sigma)

Kinv<-solve(K)
Kvecmat<-Kinv%*%dataa[(groups==0),]

library(emulator)

Z2Z2<-array(0,dim<-c(length(set3),dim(Z3)[2]))

inndex<-1

dataa1<-dataa[(groups==0),]
source("add_gplvm.R")
YKZ<-Kinv%*%Z3[(groups==0),]

for(i in set3)
{
	a5<-data3[i,]
	a5<-as.numeric(a5)
	dataanew<-rbind(dataa1,a5)
	dataset1<-t(chol(dataanew%*%t(dataanew)))
	N<-dim(dataset1)
	yvec<-dataset1[N[1],]%*%t(dataset1)
	Znew<-rep(0,(dim(Z3)[2]))
	Z1<-rbind(Z3[(groups==0),],Znew)
	Z2<-sweep(Z1,2,Z1[(N[1]),],"-")
	kvec<-alpha^2*exp(-rowSums((Z2[-N,])^2)/(2*l^2))
	Kinv11<-solve(K-(t(t(kvec))%*%(1/(alpha^2+sigma^2)))%*%t(kvec))
	sum1<-yvec[-N]%*%(Kinv11%*%Z3[(groups==0),])
	eps<-1
	while(eps>0.001)
	{
	Z1<-rbind(Z3[(groups==0),],Znew)
	Z2<-sweep(Z1,2,Z1[(N[1]),],"-")
	kvec<-alpha^2*exp(-rowSums((Z2[-N,])^2)/(2*l^2))
	kvecinv<--Kinv%*%kvec/(dataset1[N])	
	}
	ssx<-add_gplvm(dataset1,Z3[(groups==0),],Kinv,Kvecmat,l,alpha,sigma)
	Z2Z2[inndex,]<-ssx
    inndex<-inndex+1	
    print(inndex)
	
	
}











