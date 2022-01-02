library(GPLVM)
source("kernels11.R")
source("GPLVM10.R")
xxx<-fit.gplvm(data,2,iterations=100,Z.normal.prior=TRUE,groups=groups)




dL.dK <- function(X, K, groups) {
  K.chol <- chol(K)
  K.inv <- chol2inv(K.chol)
  temp <- backsolve(K.chol, forwardsolve(t(K.chol), X))
  D <- ncol(X)
  (temp %*% t(temp) - D * K.inv) / 2
}

dK.dZij <- function(Z, K, i, j, l) {
  out <- matrix(0, nrow=nrow(K), ncol=ncol(K))
  out[i, ] <- out[, i] <- (Z[, j] - Z[i, j]) / l^2 * K[i, ]
  return(out)
}


dL.dl <- function(X, Z, l, alpha, sigma, sigma1, K, dL.dK, groups) {
  dK.dl <- dK.dl(Z, K, l)
  out <- sum(dL.dK * dK.dl)
  return(out)
}


dL.dalpha <- function(X, Z, l, alpha, sigma, sigma1, K, dL.dK, groups) {
  dK.dalpha <- dK.dalpha(Z, K, alpha, sigma)
  out <- sum(dL.dK * dK.dalpha)
  return(out)
}


dL.dsigma <- function(X, Z, l, alpha, sigma, sigma1, K, dL.dK, groups) {
  dK.dsigma <- diag(2 * sigma, nrow=nrow(Z))
  out <- sum(dL.dK * dK.dsigma)
  return(out)
}


dL.dsigma1 <- function(X, Z, l, alpha, sigma, sigma1, K, dL.dK, groups) {
	XX<-array(0,dim<-c(dim(Z)[1],length(table(groups))))
	counter<-1
	for(le in (unique(groups)))
	{
		XX[,counter]<-as.numeric(groups==le)
		counter<-counter+1
	}
    dK.dsigma1<-2*sigma1*XX%*%t(XX)

  out <- sum(dL.dK * dK.dsigma1)
  return(out)
}



gplvm.f <- function(par, X,groups,q,KK,mumat,sigmat,pii) {
	l1<-par[1]
	alpha<-par[2]
	sigma<-par[3]
	sigma1<-par[4]
	par<-par[-c(1:4)]
  par1<-par
  Z <- matrix(par1, nrow=nrow(X))
  ind<-array(0,dim<-c(dim(X)[1],KK))
  
  for(i in 1:(dim(ind)[1]))
  {
  	s<-0
  	for(j in 1:(dim(mumat)[1]))
  	{
  		s<-s+pii[j]*dmvnorm(Z[i,],mumat[j,],sigmat[((j-1)*dim(mumat)[2]+1):(j*dim(mumat)[2]),])
  	}
  	if(s>0.00000001)
  	{
  	for(j in 1:(dim(ind)[2]))
  	{
  		ind[i,j]<-pii[j]*dmvnorm(Z[i,],mumat[j,],sigmat[((j-1)*dim(mumat)[2]+1):(j*dim(mumat)[2]),])/s
  	}
  	}else
  	{
  		diss<-NULL
  		for(j in 1:(dim(ind)[2]))
  		{
  			diss[j]<-sum(abs(Z[i,]-mumat[j,]))
  		}
  		for(j in 1:(dim(ind)[2]))
  		{
  			if(diss[j]==min(diss))
  			{
  				ind[i,j]<-1
  			}else
  			{
  				ind[i,j]<-0
  			}
  		}
  	}
  }
  mumat<-array(0,dim<-dim(mumat))
  sigmat<-array(1,dim<-dim(sigmat))
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
	
  gplvm.L(Z, X, l1, alpha, sigma, sigma1,groups,ind,mumat,sigmat,pii)
}


gplvm.L.from.K <- function(Z, X, K, Z.normal.prior=TRUE) {
  X <- as.matrix(X)
  D <- ncol(X)
  N <- nrow(X)
  K.chol <- chol(K)
  K.inv.X <- backsolve(K.chol, forwardsolve(t(K.chol), X))
  K.log.det <- 2 * sum(log(diag(K.chol)))
  Z.prior.term <- 0
  if (Z.normal.prior) {
    Z.prior.term <- - length(Z)/2 * log(2 * 10^2 * pi) - 1 / (2 * 10 ^2) * sum(as.numeric(Z)^2)
  }
  return(-D * N * log(2 * pi) / 2 - D * K.log.det / 2 - 1/2 * sum(K.inv.X * X) + Z.prior.term)
}

gplvm.L <- function(Z, X, l, alpha, sigma,sigma1,groups,ind,mumat,sigmat,pii) {
  K <- gplvm.SE(Z, l, alpha, sigma,sigma1,groups,ind,mumat,sigmat,pii)
  return(gplvm.L.from.K(Z, X, K, Z.normal.prior=TRUE))
}



gplvm.SE <- function(Z, l, alpha, sigma,sigma1,groups,ind,mumat,sigmat,pii) {
  gplvm.SE.dist(dist.matrix=as.matrix(dist(Z)),
                l=l,
                alpha=alpha,
                sigma=sigma,
                sigma1=sigma1,
                groups=groups,
                ind=ind,
                mumat=mumat,
                sigmat=sigmat,
                pii=pii)
}

gplvm.SE.dist <- function(dist.matrix, l, alpha, sigma=0,sigma1,groups,ind,mumat,sigmat,pii) { 
	XX<-array(0,dim<-c(dim(dist.matrix)[1],length(table(groups))))
	counter<-1
	for(le in (unique(groups)))
	{
		XX[,counter]<-as.numeric(groups==le)
		counter<-counter+1
	}


  K <- alpha^2 * exp(-dist.matrix^2/(2*l^2))
  diag(K) <- diag(K) + sigma^2
  K<-K+sigma1^2*XX%*%t(XX)
  return(K)
}




dL.dZ <- function(X, Z, l, alpha, sigma, sigma1, K, dL.dK1, groups,ind,mumat,sigmat,pii) {

  siginvmat<-sigmat
  for(k in 1:(dim(mumat)[1]))
  {
  	siginvmat[((k-1)*dim(mumat)[2]+1):(k*(dim(mumat)[2])),]<-solve(sigmat[((k-1)*dim(mumat)[2]+1):(k*(dim(mumat)[2])),])
  }
  out <- Z
  for (i in 1:nrow(Z)) {
    for (j in 1:ncol(Z)) {
      out[i, j] <- sum(dL.dK1 * dK.dZij(Z, K, i, j, l))
    }
  }
  deriv<-array(0,dim<-dim(out))
    for(k in 1:(dim(mumat)[1]))
    {
    	for(i in 1:(dim(Z)[1]))
    	{
            deriv[i,] <- deriv[i,]-ind[i,k]*siginvmat[((k-1)*dim(mumat)[2]+1):(k*(dim(mumat)[2])),]%*%(Z[i,]-mumat[k,])
    
        }
    }
  return(out+deriv)
  #return(out)
}



gplvm.gr <- function(par, X,groups,q,KK,mumat,sigmat,pii) {
	l1<-par[1]
	alpha<-par[2]
	sigma<-par[3]
	sigma1<-par[4]
	par<-par[-(1:4)]
  Z <- matrix(par, nrow=nrow(X))
  
  K <- gplvm.SE(Z, l1, alpha, sigma,sigma1,groups)
  dL.dK1 <- dL.dK(X, K)
  ind<-array(0,dim<-c(dim(X)[1],KK))
  for(i in 1:(dim(ind)[1]))
  {
  	s<-0
  	for(j in 1:(dim(mumat)[1]))
  	{
  		s<-s+pii[j]*dmvnorm(Z[i,],mumat[j,],sigmat[((j-1)*dim(mumat)[2]+1):(j*dim(mumat)[2]),])
  	}
  	if(s>0.00000001)
  	{
  	for(j in 1:(dim(ind)[2]))
  	{
  		ind[i,j]<-pii[j]*dmvnorm(Z[i,],mumat[j,],sigmat[((j-1)*dim(mumat)[2]+1):(j*dim(mumat)[2]),])/s
  	}
  	}else
  	{
  		diss<-NULL
  		for(j in 1:(dim(ind)[2]))
  		{
  			diss[j]<-sum(abs(Z[i,]-mumat[j,]))
  		}
  		for(j in 1:(dim(ind)[2]))
  		{
  			if(diss[j]==min(diss))
  			{
  				ind[i,j]<-1
  			}else
  			{
  				ind[i,j]<-0
  			}
  		}
  	}
  }
  
  c(dL.dl(X, Z, l1, alpha, sigma, sigma1, K, dL.dK1, groups),
    dL.dalpha(X, Z, l1, alpha, sigma, sigma1, K, dL.dK1, groups),
    dL.dsigma(X, Z, l1, alpha, sigma, sigma1, K, dL.dK1, groups),
    dL.dsigma1(X, Z, l1, alpha, sigma, sigma1, K, dL.dK1, groups),
    as.numeric(dL.dZ(X, Z, l1, alpha, sigma, sigma1, K, dL.dK1, groups,ind,mumat,sigmat,pii)))
}








grad_des<-function(par,fn,gr,X,groups,q,KK,mumat,sigmat,pii,l,alpha,sigma,sigma1)
{
	val<-gr(par,X,groups,q,KK,mumat,sigmat,pii,l,alpha,sigma,sigma1)
	frac<-0.01
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
		val<-gr(par1,X,groups,q,KK,mumat,sigmat,pii,l,alpha,sigma,sigma1)
		if(iter<3)
		{
			frac<-0.01
		}else
		{
			frac<-0.01
			#frac<-sum((val-valold)*(par1-par1old))/sum((par1-par1old)^2)
		}
		if(is.nan(frac))
		{
			break
		}
		par1old<-par1
		par1<-par1+frac*gr(par1,X,groups,q,KK,mumat,sigmat,pii,l,alpha,sigma,sigma1)
		
		
		
		Z<-matrix(par1,nrow=(dim(X)[1]))
		plot(Z,col=cols)
		
	}
	
	
	
}




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


par<-c(l,alpha,sigma,sigma1,as.numeric(xxx$Z))

#par<-c(l,alpha,sigma,sigma1,rnorm(600,0,1))

for(iter in 1:10)
{
	print(iter)
	l1<-l
 out <-optimx(par,gplvm.f,gplvm.gr,X=X,groups=groups,q=q,KK=KK,mumat=mumat,sigmat=sigmat,pii=pii,method="L-BFGS-B", control=list(trace=T,maximize=T,kkt=FALSE,maxit=10,starttests=FALSE))
 
 vec<-as.numeric(out)[-(1:4)]
 Z<-matrix(vec[1:(length(as.numeric(Z)))],nrow=(dim(X)[1]))
 plot(Z,col=cols)
 
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
	pii<-colMeans(ind)
		for(i in 1:(dim(Z)[1]))
	{
		s<-0
	    for(k in 1:K)
	    {
		    s<-s+pii[k]*dmvnorm(Z[i,],mumat[k,],sigmat[((k-1)*q+1):(k*q),])
	    }
	    if(s>0)
	    {
	    for(k in 1:K)
	    {
		ind[i,k]<-(pii[k]*dmvnorm(Z[i,],mumat[k,],sigmat[((k-1)*q+1):(k*q),]))/s
		}
		}else
  	{
  		diss<-NULL
  		for(j in 1:(dim(ind)[2]))
  		{
  			diss[j]<-sum(abs(Z[i,]-mumat[j,]))
  		}
  		for(j in 1:(dim(ind)[2]))
  		{
  			if(diss[j]==min(diss))
  			{
  				ind[i,j]<-1
  			}else
  			{
  				ind[i,j]<-0
  			}
  		}
  	}
	}
	#sigmat<-sigmat/2
	
	
}







