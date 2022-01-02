#set.seed(0)
library(mvtnorm)
data1<-rmvnorm(20,rep(-3,300),diag(300))
data2<-rmvnorm(50,rep(3,300),diag(300))
data3<-rmvnorm(80,rep(6,300),diag(300))
data4<-rmvnorm(20,rep(-3,300),diag(300))
data5<-rmvnorm(50,rep(3,300),diag(300))
data6<-rmvnorm(80,rep(6,300),diag(300))
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

cols<-c(rep("red",20),rep("black",N),rep("blue",80),rep("green",20),rep("purple",N),rep("cyan",80))




iterations=1000
plot.freq=100
classes=NULL
Z.init="PCA"
num.init.params=100
Z.normal.prior<-TRUE



library(GPLVM)
xxxold<-fit.gplvm(data,2,iterations=100,Z.normal.prior=TRUE)



source("kernels11_random.R")
source("GPLVM10_random.R")
xxx<-fit.gplvm(data,2,iterations=100,Z.normal.prior=TRUE,groups=groups,sigma1_init=1)





source("kernels11_prior.R")
source("GPLVM10_prior.R")
X<-data
Z<-xxx$Z
K<-3
q<-2

eps<-1
pii<-rep(1/K,K)
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






XX<-array(0,dim<-c(dim(data)[1],length(table(groups))))
counter<-1
for(le in (unique(groups)))
{
	XX[,counter]<-as.numeric(groups==le)
	counter<-counter+1
}

data1<-data
for(j in 1:(dim(data)[2]))
{
	y<-data[,j]
	fit<-lm(y~-1+XX)
    data1[,j]<-data[,j]-fit$fitted
}


source("kernels11_prior.R")
source("GPLVM10_prior.R")
q<-2
xxx1<-xxx
for(iter in 1:10)
{
	print(iter)
	plot(xxx1$Z,col=cols)
	dev.new()
	xxx1<-fit.gplvm(data,q,iterations=100,Z.normal.prior=TRUE,groups=groups,mumat=mumat,sigmat=sigmat,pii=pii,Z.init=xxx1$Z)
	Z<-xxx1$Z
	
		for(i in 1:(dim(Z)[1]))
	{
		s<-0
	    for(k in 1:KK)
	    {
		    s<-s+pii[k]*dmvnorm(Z[i,],mumat[k,],sigmat[((k-1)*q+1):(k*q),])
	    }
	    for(k in 1:KK)
	    {
		ind[i,k]<-(pii[k]*dmvnorm(Z[i,],mumat[k,],sigmat[((k-1)*q+1):(k*q),]))/s
		}
	}
	
	s1<-array(0,dim<-c(q,KK))
	s2<-array(0,dim<-c(q,KK))
	
	for(k in 1:KK)
	{
		for(i in 1:(dim(Z)[1]))
		{
			if(groups[i]==1)
			{
				s1[,k]<-s1[,k]+ind[i,k]*Z[i,]
			}else
			{
				s2[,k]<-s2[,k]+ind[i,k]*Z[i,]
			}
		}
	}
	for(k in 1:KK)
	{
		s1[,k]<-s1[,k]/sum(ind[(groups==1),k])
		s2[,k]<-s2[,k]/sum(ind[(groups==0),k])
	}
	batch_effect<-rowMeans(s1-s2)
	
	
	for(i in 1:(dim(Z)[1]))
	{
		if(groups[i]==1)
		{
			Z[i,]<-Z[i,]-batch_effect
		}
	}
	
	
eps<-1

while(eps>0.001)
{
	mumatold<-mumat
	
	for(i in 1:(dim(Z)[1]))
	{
		s<-0
	    for(k in 1:KK)
	    {
		    s<-s+pii[k]*dmvnorm(Z[i,],mumat[k,],sigmat[((k-1)*q+1):(k*q),])
	    }
	    for(k in 1:KK)
	    {
		ind[i,k]<-(pii[k]*dmvnorm(Z[i,],mumat[k,],sigmat[((k-1)*q+1):(k*q),]))/s
		}
	}
	pii<-colMeans(ind)
	for(k in 1:KK)
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


	
	
	
	dis<-as.matrix(dist(xxx1$Z))
	Sigma<-(xxx1$alpha)^2*exp(-dis^2/(2*(xxx1$l)^2))+(xxx1$sigma)^2*diag(dim(xxx1$Z)[1])
	gamma<-array(0,dim<-c(q,(dim(data)[2])))
	mat1<-solve(t(XX)%*%solve(Sigma)%*%(XX))
	mat2<-t(XX)%*%solve(Sigma)
	for(j in 1:(dim(data)[2]))
	{
		gamma[,j]<-mat1%*%mat2%*%data1[,j]
	}
	#data1<-data1-XX%*%gamma
	
	
}











