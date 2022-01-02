X<-array(rnorm(707,0,1),dim<-c(101,7))

X1<-X[1:100,]
x<-X[101,]

d<-NULL
for(i in 1:100)
{
	d[i]<-sum((X1[i,]-x)^2)
}

y<-rep(0,7)
eps<-1


cost<-function(y,d,X1)
{
	s<-0
	for(i in 1:(dim(X1)[1]))
	{
		s<-s+(d[i]-sum((y-X1[i,])^2))^2
	}
	return(s)
}

opt<-optim(y,cost,d=d,X1=X1,method="L-BFGS-B")







