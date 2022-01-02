
cluster<-c(rep(1,30),rep(2,N),rep(3,70),rep(1,30),rep(2,N),rep(3,70))

source("kernels11_random.R")
source("GPLVM10_random.R")
xxx<-fit.gplvm(data,2,iterations=100,Z.normal.prior=TRUE,groups=groups,sigma1_init=1)

Z<-xxx$Z
q<-dim(Z)[2]
muu<-array(0,dim<-c(max(cluster),q))
sigmaa<-array(0,dim<-c((q*max(cluster)),q))
for(k in 1:(max(cluster)))
{
	muu[k,]<-colMeans(Z[(cluster==k),])
	sigmaa[(((k-1)*q+1):(k*q)),]<-(var(Z[((cluster==k)&(groups==0)),])*(sum((cluster==k)&(groups==0))-1)+var(Z[((cluster==k)&(groups==1)),])*(sum((cluster==k)&(groups==1))-1))/(sum(cluster==k)-2)
}





X<-data
q<-2
iterations=1000
plot.freq=100
classes=NULL
num.init.params=100
Z.normal.prior=TRUE

sigma1_init=1
sigma2_init=1


source("GPLVM10_random_prior1.R")
source("kernels11_random_prior.R")


xxxnew<-fit.gplvm.prior(X=data,q=2,iterations=100,Z.init=ZZZ,Z.normal.prior=TRUE,groups=groups,sigma1_init=1,muu=muu,sigmaa=sigmaa/1000000,cluster=cluster)



source("GPLVM10_random.R")
source("kernels11_random.R")
xxx2<-fit.gplvm(X=data,q=2,iterations=100,Z.init=ZZZ,Z.normal.prior=TRUE,groups=groups,sigma1_init=1)








