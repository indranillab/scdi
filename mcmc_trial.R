seqq<-NULL
for(i in 1:1000)
{
	print(i)
	a<-5
	for(j in 1:5000)
	{
		b<-floor(11*runif(1))
		prob<-min(1,dbinom(b,10,1/2)/dbinom(a,10,1/2))
		u<-runif(1)
		if(u<prob)
		{
			a<-b
		}
		
	}
	seqq[i]<-a
}















