set.seed(0)

library(mvtnorm)

matt1<-array(c(2,1,1,2),dim<-c(2,2))
matt2<-array(c(10,3,3,10),dim<-c(2,2))


X<-rmvnorm(100,c(0,0),matt1)
Y<-rmvnorm(100,c(10,10),matt2)

data<-rbind(X,Y)



u1<-sort(sample(1:200,size=100))
u2<-setdiff(1:200,u1)


eps<-1
pi<-0.5
ind<-rep(0.5,200)
mu1<-colMeans(data[u1,])
mu2<-colMeans(data[u2,])
sig1<-var(data[u1,])
sig2<-var(data[u2,])

while(eps>0.001)
{
	mu1old<-mu1
	mu2old<-mu2
	for(i in 1:(dim(data)[1]))
	{
		ind[i]<-(pi*dmvnorm(data[i,],mu2,sig2))/(pi*dmvnorm(data[i,],mu2,sig2)+(1-pi)*dmvnorm(data[i,],mu1,sig1))
	}
	pi<-mean(ind)
	mu2<-c((sum(data[,1]*ind)/sum(ind)),(sum(data[,2]*ind)/sum(ind)))
	mu1<-c((sum(data[,1]*(1-ind))/sum(1-ind)),(sum(data[,2]*(1-ind))/sum(1-ind)))
	sig2<-0
	sig1<-0
	for(i in 1:(dim(data)[1]))
	{
		sig2<-sig2+ind[i]*(data[i,]-mu2)%*%t(data[i,]-mu2)
		sig1<-sig1+(1-ind[i])*(data[i,]-mu1)%*%t(data[i,]-mu1)
	}
	sig2<-sig2/sum(ind)
	sig1<-sig1/sum(1-ind)
	eps<-sum(abs(mu1old-mu1))+sum(abs(mu2old-mu2))
	print(eps)
}


cols<-NULL
cols[ind<0.5]<-"red"
cols[ind>=0.5]<-"black"












