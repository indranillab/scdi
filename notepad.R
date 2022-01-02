

eps<-1
pi<-rep(1/K,K)
ind<-array(1/K,dim<-c(dim(X)[1],K))
mumat<-array(0,dim<-c(K,q))
sigmat<-array(1,dim<-c((K*q),q))

for(i in 1:K)
{
	u<-sample(1:(dim(X)[1]),floor(dim(X)[1]/K))
	mumat[i,]<-colMeans(Z[u,])
	sigmat[((i-1)*q+1):(i*q),]<-var(Z[u,])
}

while(eps>0.001)
{
	mumatold<-mumat
	
	for(i in 1:(dim(Z)[1]))
	{
		s<-0
	    for(k in 1:K)
	    {
		    s<-s+pi[k]*dmvnorm(Z[i,],mumat[k,],sigmat[((k-1)*q+1):(k*q),])
	    }
	    for(k in 1:K)
	    {
		ind[i,k]<-(pi[k]*dmvnorm(Z[i,],mumat[k,],sigmat[((k-1)*q+1):(k*q),]))/s
		}
	}
	pi<-colMeans(ind)
	for(k in 1:K)
	{
        for(k1 in 1:q)
        {
        	mumat[k,q]<-sum(Z[,k1]*ind[,k])/sum(ind[,k])
        }		
        s<-0
        for(i in 1:(dim(Z)[1]))
        {
        	s<-s+ind[i,k]*(Z[i,]-mumat[k,])%*%t(Z[i,]-mumat[k,])
        }
        sigmat[((k-1)*q+1):(k*q),]<-s/(sum(ind[,k]))
        
	}
	
	eps<-sum(abs(mumatold-mumat))
	print(eps)
}





cols<-NULL
for(i in 1:(dim(ind)[1]))
{
	if(ind[i,1]==max(ind[i,]))
	{
		cols[i]<-"black"
	}
	if(ind[i,2]==max(ind[i,]))
	{
		cols[i]<-"blue"
	}
	if(ind[i,3]==max(ind[i,]))
	{
		cols[i]<-"red"
	}
	if(ind[i,4]==max(ind[i,]))
	{
		cols[i]<-"green"
	}
	
}



par<-c(xxx$l,xxx$alpha,xxx$sigma,xxx$sigma1,as.numeric(xxx$Z))
source("GPLVM11.R")






iterations=1000
plot.freq=100
classes=NULL
Z.init="PCA"
num.init.params=100


source("kernels11.R")
library(GPLVM)


source("GPLVM10.R")
xxx<-fit.gplvm(data,2,iterations=100,Z.normal.prior=TRUE,groups=groups)



source("gradient_descent.R")
fn<-gplvm.f
gr<-gplvm.gr
X<-data
q<-2
KK<-3
K<-KK
Z<-xxx$Z
l<-xxx$l
alpha<-xxx$alpha
sigma<-xxx$sigma
sigma1<-xxx$sigma1



eps<-1
pii<-rep(1/K,K)
ind<-array(1/K,dim<-c(dim(X)[1],K))
mumat<-array(0,dim<-c(K,q))
sigmat<-array(1,dim<-c((K*q),q))

vec<-floor(K*runif((dim(X)[1])))+1
for(i in 1:K)
{
	mumat[i,]<-colMeans(Z[(vec==i),])
	sigmat[((i-1)*q+1):(i*q),]<-var(Z[(vec==i),])
}

while(eps>0.001)
{
	mumatold<-mumat
	
	for(i in 1:(dim(Z)[1]))
	{
		s<-0
	    for(k in 1:K)
	    {
		    s<-s+pii[k]*dmvnorm(Z[i,],mumat[k,],sigmat[((k-1)*q+1):(k*q),])
	    }
	    for(k in 1:K)
	    {
		ind[i,k]<-(pii[k]*dmvnorm(Z[i,],mumat[k,],sigmat[((k-1)*q+1):(k*q),]))/s
		}
	}
	pii<-colMeans(ind)
	for(k in 1:K)
	{
        for(k1 in 1:q)
        {
        	mumat[k,q]<-sum(Z[,k1]*ind[,k])/sum(ind[,k])
        }		
        s<-0
        for(i in 1:(dim(Z)[1]))
        {
        	s<-s+ind[i,k]*(Z[i,]-mumat[k,])%*%t(Z[i,]-mumat[k,])
        }
        sigmat[((k-1)*q+1):(k*q),]<-s/(sum(ind[,k]))
        
	}
	
	eps<-sum(abs(mumatold-mumat))
	print(eps)
}



source("gradient_descent.R")
par<-as.numeric(xxx$Z)
fn<-gplvm.f
gr<-gplvm.gr






par






source("GPLVM11.R")
xxx1<-fit.gplvm(data,2,3,iterations=100,Z.init=xxx$Z,groups=groups)







