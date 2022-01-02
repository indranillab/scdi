

set.seed(12345)
X1<-array(rnorm(10000,0,1),dim<-c(100,100))
X2<-array(rnorm(10000,1,1),dim<-c(100,100))

X<-rbind(X1,X2)


library(GPLVM)

Z<-fit.gplvm(X,2,Z.normal.prior=TRUE,iterations=100)


set.seed(12346)
X3<-array(rnorm(10000,0,1),dim<-c(100,100))
X4<-array(rnorm(10000,1,1),dim<-c(100,100))

alpha<-Z$alpha
l<-Z$l
sigma<-Z$sigma

source("kernels.R")
K<-gplvm.SE(Z$Z,Z$l,Z$alpha,Z$sigma)


matt<-X%*%t(X)

X5<-rbind(X3,X4)

Z1<-array(0,dim<-dim(Z$Z))
for(i in 1:100)
{
	wt<-X1%*%X5[i,]/(dim(X)[2])
	#wt<-(wt-min(wt))/(range(wt)[2]-range(wt)[1])
	Z1[i,]<-(t(Z$Z[1:100,])%*%wt)/sum(wt)
}


Z2<-rbind(Z$Z,Z1[1:100,])


cols<-rep(c(rep(1,100),rep(2,100)),2)






