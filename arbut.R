
grad_des<-function(par,fn,gr,X,groups,q,KK,mumat,sigmat,pii,l,alpha,sigma,sigma1)
{
	val<-gr(c(l,alpha,sigma,sigma1,par),X,groups,Z.normal.prior=TRUE)[-(1:4)]
	frac<-1
	par1<-par
	for(iter in 1:100)
	{
		print(iter)
		
		mumat<-array(0,dim<-dim(mumat))
        sigmat<-array(1,dim<-dim(sigmat))
        pii<-colMeans(ind)
	    for(k in 1:KK)
	    {
        for(k1 in 1:q)
        {
        	mumat[k,k1]<-sum(Z[,k1]*ind[,k])/sum(ind[,k])
        }		
        s<-0
        for(i in 1:(dim(Z)[1]))
        {
        	s<-s+ind[i,k]*(Z[i,]-mumat[k,])%*%t(Z[i,]-mumat[k,])
        }
        sigmat[((k-1)*q+1):(k*q),]<-s/(sum(ind[,k]))
        
	    }
	    
	    valold<-val
		val<-gr(c(l,alpha,sigma,sigma1,par),X,groups,Z.normal.prior=TRUE)[-(1:4)]
		if(iter<3)
		{
			frac<-0.1/iter
		}else
		{
			frac<-0.1/iter
			#frac<-sum((val-valold)*(par1-par1old))/sum((par1-par1old)^2)
		}
		if(is.nan(frac))
		{
			break
		}
		par1old<-par1
		par1<-par1+frac*gr(c(l,alpha,sigma,sigma1,par),X,groups,Z.normal.prior=TRUE)[-(1:4)]
		
		
		
		Z<-matrix(par1,nrow=(dim(X)[1]))
		plot(Z,col=cols)
		
	}
	
	
	
}






