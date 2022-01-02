library(mvtnorm)

KK<-3
eps<-1
TT<-2

pii<-rep(1/KK,KK)
thetamat<-array(0,dim<-c(KK,(dim(data)[2])))
betamat<-array(0,dim<-c(TT,(dim(data)[2])))
sigmat<-array(0,dim<-c((KK*(dim(data)[2])),(dim(data)[2])))

xxx<-kmeans(data,KK)
for(i in 1:KK)
{
	thetamat[i,]<-colMeans(data[(xxx$cluster==i),])
	sigmat[((i-1)*dim(data)[2]+1):(i*(dim(data)[2])),]<-(1/(KK*dim(data)[2]))*sum(diag(var(data)))*diag(1,(dim(data)[2]))
}


while(eps>0.000001)
{
	piiold<-pii
	thetamatold<-thetamat
	betamatold<-betamat
	sigmatold<-sigmat
	gamma<-array(0,dim<-c(dim(data)[1],KK))
	for(i in 1:(dim(data)[1]))
	{
		if(i<=150)
		{
			j<-1
		}else
		{
			j<-2
		}
			for(k in 1:KK)
			{
				gamma[i,k]<-pii[k]*dmvnorm(data[i,],(betamat[j,]+thetamat[k,]),sigmat[((k-1)*dim(data)[2]+1):(k*(dim(data)[2])),])
			}
			if(sum(gamma[i,])>0)
			{
			gamma[i,]<-gamma[i,]/sum(gamma[i,])
			}else
			{
				sss<-NULL
				for(k in 1:KK)
				{
					sss[k]<-sum(abs(data[i,]-(betamat[j,]+thetamat[k,])))
				}
				for(k in 1:KK)
				{
					if(sss[k]==min(sss))
					{
						gamma[i,k]<-1
					}else
					{
						gamma[i,k]<-0
					}
				}
			}
	}
	
	for(k in 1:KK)
	{
	    s<-rep(0,(dim(data)[2]))
		frac<-0
		for(j in 1:TT)
		{
			for(i in (((j-1)*150+1):((j)*150)))
			{
				frac<-frac+gamma[i,k]
				s<-s+gamma[i,k]*(data[i,]-betamat[j,])
			}
		}
		thetamat[k,]<-s/frac
	}
	
	for(j in 1:TT)
	{
	    s<-rep(0,(dim(data)[2]))
		frac<-0
		for(k in 1:KK)
		{
			for(i in (((j-1)*150+1):((j)*150)))
			{
				frac<-frac+gamma[i,k]
				s<-s+gamma[i,k]*(data[i,]-thetamat[k,])
			}
		}
		betamat[j,]<-s/frac
	}
	
	for(k in 1:KK)
	{
	    s<-array(0,dim<-c(dim(data)[2],dim(data)[2]))
		frac<-0
		for(j in 1:TT)
		{
			for(i in (((j-1)*150+1):((j)*150)))
			{
				frac<-frac+gamma[i,k]
				s<-s+gamma[i,k]*(data[i,]-betamat[j,]-thetamat[k,])%*%t(data[i,]-betamat[j,]-thetamat[k,])
			}
		}
		sigmat[((k-1)*dim(data)[2]+1):(k*(dim(data)[2])),]<-s/frac
	}
	
    for(k in 1:KK)
    {
    	pii[k]<-sum(gamma[,k])/(dim(data)[1])
    }
    
    eps<-sum(abs(piiold-pii))+sum(abs(thetamatold-thetamat))+sum(abs(betamatold-betamat))+sum(abs(sigmatold-sigmat))
    print(eps)
    print(pii)
    dataadj<-data
    for(i in 1:(dim(data)[1]))
    {
    	if(i<=150)
    	{
    		j<-1
    	}else
    	{
    		j<-2
    	}
    	dataadj[i,]<-data[i,]-betamat[j,]
    }	
}



maxx<-NULL
for(i in 1:(dim(gamma)[1]))
{
	maxx[i]<-which(gamma[i,]==max(gamma[i,]))
}







