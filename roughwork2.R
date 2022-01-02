set.seed(12345)
X1<-array(rnorm(10000,0,1),dim<-c(100,100))
X2<-array(rnorm(10000,1,1),dim<-c(100,100))

X<-rbind(X1,X2)


library(GPLVM)

Z<-fit.gplvm(X,2,Z.normal.prior=TRUE,iterations=100)

X3<-array(rnorm(10000,0,1),dim<-c(100,100))
X4<-array(rnorm(10000,1,1),dim<-c(100,100))

X5<-rbind(X3,X4)

classify()


for(i in 1:(dim(X5)[1]))
{
	ind<-classify(X5[i,],X1,X2)
	
	
}









