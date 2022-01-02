X<-array(rnorm(202,0,1),dim<-c(101,2))

X1<-X[1:100,]
x<-X[101,]

d<-NULL
for(i in 1:100)
{
	d[i]<-exp(sum((X1[i,]-x)^2))/(1+exp(sum((X1[i,]-x)^2)))
}


y<-rep(0,2)
eps<-1



while(eps>0.000000000001)
{
	yold<-y
	dd<-NULL
	for(i in 1:100)
	{
		dd[i]<-(exp(sum((X1[i,]-y)^2))/(1+exp(sum((X1[i,]-y)^2)))-d[i])*exp(sum((X1[i,]-y)^2))/(1+exp(sum((X1[i,]-y)^2)))^2

	}
	y<-(dd%*%X1)/sum(dd)
	eps<-sum(abs(y-yold))
	print(eps)
	plot(rbind(X1,y),col=c(rep("black",100),"red"))
}





