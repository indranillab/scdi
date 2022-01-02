
dL.dK <- function(X, K, groups) {
  K.chol <- chol(K)
  K.inv <- chol2inv(K.chol)
  temp <- backsolve(K.chol, forwardsolve(t(K.chol), X))
  D <- ncol(X)
  (temp %*% t(temp) - D * K.inv) / 2
}


dL.dZ <- function(X, Z, l, alpha, sigma, sigma1, K, dL.dK, groups,ind,mumat,sigmat,pii) {

  siginvmat<-sigmat
  for(k in 1:(dim(mumat)[1]))
  {
  	siginvmat[((k-1)*dim(mumat)[2]+1):(k*(dim(mumat)[2])),]<-solve(sigmat[((k-1)*dim(mumat)[2]+1):(k*(dim(mumat)[2])),])
  }
  out <- Z
  for (i in 1:nrow(Z)) {
    for (j in 1:ncol(Z)) {
      out[i, j] <- sum(dL.dK * dK.dZij(Z, K, i, j, l))
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
}



dL.dmumat<-function(X, Z, l, alpha, sigma, sigma1, K, dL.dK, groups,ind,mumat,sigmat,pii)
{
	deriv<-array(0,dim<-dim(mumat))
	siginvmat<-sigmat
	K<-(dim(sigmat)[1])/(dim(sigmat)[2])
	q<-dim(sigmat)[2]
	for(k in 1:(dim(deriv)[1]))
	{
		for(i in 1:(dim(Z)[1]))
		{
		    deriv[k,]<-deriv[k,]+ind[i,k]*solve(sigmat[((k-1)*q+1):(k*q),])%*%(Z[i,]-mumat[k,])
		}
	}
	return(deriv)
}

dL.dsigmat<-function(X, Z, l, alpha, sigma, sigma1, K, dL.dK, groups,ind,mumat,sigmat,pii)
{
	deriv<-array(0,dim<-dim(sigmat))
	siginvmat<-sigmat
	K<-(dim(sigmat)[1])/(dim(sigmat)[2])
	q<-dim(sigmat)[2]
	for(k in 1:(dim(mumat)[1]))
	{
		S<-0
		sigg<-sigmat[((k-1)*q+1):(k*q),]
		for(i in 1:(dim(Z)[1]))
		{
		    S<-S-ind[i,k]/2*solve(sigg)+ind[i,k]/2*solve(sigg)%*%(Z[i,]-mumat[k,])%*%t(Z[i,]-mumat[k,])%*%solve(sigg)
		}
		deriv[((k-1)*q+1):(k*q),]<-S
	}
	return(deriv)
}

dL.dpii<-function(X, Z, l, alpha, sigma, sigma1, K, dL.dK, groups,ind,mumat,sigmat,pii)
{
    deriv<-rep(0,length(pii))
    for(k in 1:(length(deriv)))
    {
        for(i in 1:(dim(Z)[1]))
        {
        	deriv[k]<-deriv[k]+ind[i,k]/pii[k]-(1-ind[i,k])/(1-pii[k])
        }	
    }
    return(deriv)
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


gplvm.L.from.K <- function(Z, X, K, ind,mumat,sigmat,pii) {
  X <- as.matrix(X)
  D <- ncol(X)
  N <- nrow(X)
  K.chol <- chol(K)
  K.inv.X <- backsolve(K.chol, forwardsolve(t(K.chol), X))
  K.log.det <- 2 * sum(log(diag(K.chol)))
  Z.prior.term <- 0
  for(k in 1:(dim(mumat)[1]))
  {
  	for(i in 1:(dim(Z)[1]))
  	{
  		K1.chol <- chol(sigmat[((k-1)*(dim(mumat)[2])+1):(k*(dim(mumat)[2])),])
        K1.log.det <- 2 * sum(log(diag(K1.chol)))
        Z.prior.term <-Z.prior.term+ind[i,k]*( -log(2*pi)*(dim(mumat)[2])/2-K1.log.det/2 -1/2*t(Z[i,]-mumat[k,])%*%solve(sigmat[((k-1)*(dim(mumat)[2])+1):(k*(dim(mumat)[2])),])%*%(Z[i,]-mumat[k,])+log(pii[k]))
    }
  }
  
  return(as.numeric(-D * N * log(2 * pi) / 2 - D * K.log.det / 2 - 1/2 * sum(K.inv.X * X) + Z.prior.term))
}

gplvm.L <- function(Z, X, l, alpha, sigma,sigma1,groups,ind,mumat,sigmat,pii) {
  K <- gplvm.SE(Z, l, alpha, sigma,sigma1,groups,ind,mumat,sigmat,pii)
  return(gplvm.L.from.K(Z, X, K, ind,mumat,sigmat,pii))
}

gplvm.f <- function(par, X,groups,q,KK) {
  l <- par[1]
  alpha <- par[2]
  sigma <- par[3]
  sigma1 <- par[4]
  par1<-par[-(1:4)]
  Z <- matrix(par1[(1:(dim(X)[1]*q))], nrow=nrow(X))
  par2<-par1[-(1:(dim(X)[1]*q))]
  mumat<-matrix(par2[1:(KK*q)],nrow=KK)
  par3<-par2[-(1:(KK*q))]
  sigmat<-matrix(par3[1:(KK*q^2)],nrow=(KK*q))
  par4<-par3[-(1:(KK*q^2))]
  pii<-par4
  
  for(i in 1:(dim(ind)[1]))
  {
  	s<-0
  	for(j in 1:(dim(mumat)[1]))
  	{
  		s<-s+pii[j]*dmvnorm(Z[i,],mumat[j,],sigmat[((j-1)*dim(mumat)[2]+1):(j*dim(mumat)[2]),])
  	}
  	for(j in 1:(dim(ind)[2]))
  	{
  		ind[i,j]<-pii[j]*dmvnorm(Z[i,],mumat[j,],sigmat[((j-1)*dim(mumat)[2]+1):(j*dim(mumat)[2]),])/s
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
	
  gplvm.L(Z, X, l, alpha, sigma, sigma1,groups,ind,mumat,sigmat,pii)
}




gplvm.gr <- function(par, X,groups,q,KK) {
  l <- par[1]
  alpha <- par[2]
  sigma <- par[3]
  sigma1 <- par[4]
  par1<-par[-(1:4)]
  Z <- matrix(par1[(1:(dim(X)[1]*q))], nrow=nrow(X))
  par2<-par1[-(1:(dim(X)[1]*q))]
  mumat<-matrix(par2[1:(KK*q)],nrow=KK)
  par3<-par2[-(1:(KK*q))]
  sigmat<-matrix(par3[1:(KK*q^2)],nrow=(KK*q))
  par4<-par3[-(1:(KK*q^2))]
  pii<-par4
  
  K <- gplvm.SE(Z, l, alpha, sigma,sigma1,groups)
  dL.dK <- dL.dK(X, K)
  for(i in 1:(dim(ind)[1]))
  {
  	s<-0
  	for(j in 1:(dim(mumat)[1]))
  	{
  		s<-s+pii[j]*dmvnorm(Z[i,],mumat[j,],sigmat[((j-1)*dim(mumat)[2]+1):(j*dim(mumat)[2]),])
  	}
  	for(j in 1:(dim(ind)[2]))
  	{
  		ind[i,j]<-pii[j]*dmvnorm(Z[i,],mumat[j,],sigmat[((j-1)*dim(mumat)[2]+1):(j*dim(mumat)[2]),])/s
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
	
  
  c(dL.dl(X, Z, l, alpha, sigma, sigma1, K, dL.dK,groups),
    dL.dalpha(X, Z, l, alpha, sigma, sigma1, K, dL.dK,groups),
    dL.dsigma(X, Z, l, alpha, sigma, sigma1, K, dL.dK,groups),
    dL.dsigma1(X, Z, l, alpha, sigma, sigma1, K, dL.dK, groups),
    as.numeric(dL.dZ(X, Z, l, alpha, sigma, sigma1, K, dL.dK, groups,ind,mumat,sigmat,pii)),
    as.numeric(dL.dmumat(X, Z, l, alpha, sigma, sigma1, K, dL.dK, groups,ind,mumat,sigmat,pii)),
    as.numeric(dL.dsigmat(X, Z, l, alpha, sigma, sigma1, K, dL.dK, groups,ind,mumat,sigmat,pii)),
    as.numeric(dL.dpii(X, Z, l, alpha, sigma, sigma1, K, dL.dK, groups,ind,mumat,sigmat,pii))
  )
}


gplvm.hp.f <- function(par, Z, X,groups,ind,mumat,sigmat,pii) {
  l <- par[1]
  alpha <- par[2]
  sigma <- par[3]
  sigma1<-par[4]
  gplvm.L(Z, X, l, alpha, sigma,sigma1, groups, ind,mumat,sigmat,pii)
}

gplvm.hp.gr <- function(par, Z, X, groups, Z.normal.prior) {
  l <- par[1]
  alpha <- par[2]
  sigma <- par[3]
  sigma1<-par[4]
  K <- gplvm.SE(Z, l, alpha, sigma,sigma1,groups)
  dL.dK <- dL.dK(X, K)
  c(dL.dl(X, Z, l, alpha, sigma, sigma1, K, dL.dK, groups),
    dL.dalpha(X, Z, l, alpha, sigma, sigma1, K, dL.dK, groups),
    dL.dsigma(X, Z, l, alpha, sigma, sigma1, K, dL.dK, groups),
    dL.dsigma1(X, Z, l, alpha, sigma, sigma1, K, dL.dK, groups)
  )
}


#' Fit an unconstrained GPLVM model
#'
#' @param X
#' @param q
#' @param iterations
#' @param plot.freq
#' @param classes
#' @param Z.init
#' @param num.init.params
#' @param Z.normal.prior
#'
#' @return
#' @export
#'
#' @importFrom optimx optimx
fit.gplvm <- function(X,
                      q,
                      K,
                      iterations=1000,
                      plot.freq=100,
                      classes=NULL,
                      Z.init="PCA",
                      num.init.params=100,
                      groups) {
		Xnew<-array(0,dim<-dim(X))
		for(le in groups)
		{
			X1<-X[(groups==le),]
#	if ('dist' %in% class(X1)) {
#		n = attr(X1,'Size')
#		}else{
#		X1 = as.matrix(X1)
#		X1 = X1 - min(X1)
#		X1 = X1/max(X1)
#		initial_dims = ncol(X1)
#	    X1<- .whiten(as.matrix(X1),n.comp=initial_dims)
#	}
    X1<-t(t(X1)-colMeans(X1))
	Xnew[(groups==le),]<-X1
	}
	#X<-Xnew
	
	
	
KK<-K                      	
                      	
  if (is.null(Z.init)) {
    Z <- matrix(rnorm(nrow(X)*q, sd=0.2), ncol=q)
  } else if (identical(Z.init, "PCA")) {
    X.pca <- prcomp(X)
    Z <- X.pca$x[,1:q]
  } else {
    if (ncol(as.matrix(Z.init)) != q) warning("Mismatch between Z.init and q")
    Z <- as.matrix(Z.init)[,1:q]
  }

  if (plot.freq == 0) {
    plot.freq <- iterations
  }
  
  
  

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
        	mumat[k,k1]<-sum(Z[,k1]*ind[,k])/sum(ind[,k])
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




  if (q > 2) {
    pairs(Z, col=classes)
  } else if (q == 2) {
    plot(Z, col=classes)
  } else {
    plot(as.numeric(Z), rep(0, length(Z)), col=classes)
  }


  par.init <- c(1, 1, 1,1)
  par.init.L <- gplvm.hp.f(par.init, Z=Z, X=X,groups=groups,ind=ind,mumat=mumat,sigmat=sigmat,pii=pii)
  #first optimize the hyperparams for the initial Z
  for (i in seq(num.init.params)) {
    test.par.init <- rgamma(4, 1)
    test.par.init.L <- gplvm.hp.f(test.par.init, Z=Z, X=X,groups=groups, ind=ind,mumat=mumat,sigmat=sigmat,pii=pii)
    if (test.par.init.L > par.init.L) {
      par.init <- test.par.init
      par.init.L <- test.par.init.L
    }
  }


  optout <- optimx(par.init,
                   gplvm.hp.f,
                   Z=Z,
                   X=X,
                   groups=groups,
                   ind=ind,
                   mumat=mumat,
                   sigmat=sigmat,
                   pii=pii,
                   method="L-BFGS-B",
                   control=list(trace=T,
                                maximize=T,
                                kkt=FALSE,
                                maxit=1000,
                                starttests=TRUE)
  )



  init.l <- as.numeric(optout[1])
  init.alpha <- as.numeric(optout[2])
  init.sigma <- as.numeric(optout[3])
  init.sigma1<-as.numeric(optout[4])

  cat("Starting params: l=", init.l, "; alpha=", init.alpha, "; sigma=", init.sigma, fill=TRUE)

  par <- c(init.l, init.alpha, init.sigma, init.sigma1, as.numeric(Z),as.numeric(mumat),as.numeric(sigmat),as.numeric(pii))

  its <- 0
  convcode=1
  while (its < iterations & convcode != 0) {
    out <- optimx(par,
                  gplvm.f,
                  gr=gplvm.gr,
                  X=X,
                  groups=groups,
                  q=q,
                  KK=KK,
                  method="L-BFGS-B",
                  control=list(trace=T,
                               maximize=T,
                               kkt=FALSE,
                               maxit=plot.freq,
                               starttests=FALSE)
    )
    
    
    convcode <- out$convcode
    par <- head(as.numeric(out), length(Z) + 4)
    l <- par[1]
    alpha <- par[2]
    sigma <- par[3]
    sigma1<-par[4]
    par1<-par[-(1:4)]
    Z <- matrix(par1[(1:(dim(X)[1]*q))], ncol=q)
    if (q > 2) {
      pairs(Z,
            col=classes,
            main=paste("l=", signif(l, 3), "; alpha=", signif(alpha, 3), "; sigma=", signif(sigma, 3), sep=""))
    } else if (q == 2) {
      plot(Z,
           col=classes,
           main=paste("l=", signif(l, 3), "; alpha=", signif(alpha, 3), "; sigma=", signif(sigma, 3), sep=""))
    } else {
      plot(as.numeric(Z),
           rep(0, length(Z)),
           col=classes,
           main=paste("l=", signif(l, 3), "; alpha=", signif(alpha, 3), "; sigma=", signif(sigma, 3), sep=""))
    }

    its <- its + plot.freq
  }

  return(list(Z=Z, l=l, alpha=alpha, sigma=sigma, convcode=convcode,sigma1=sigma1))

}


#' Sparse GPLVM
#'
#' @param X
#' @param q
#' @param iterations
#' @param Z.init
#' @param classes
#' @param Z.normal.prior
#'
#' @return
#' @export
sparse.gplvm <- function(X,
                         q,
                         active.set.size,
                         iterations=10,
                         Z.init="PCA",
                         classes=1,
                         Z.normal.prior=T
) {
  if (is.null(Z.init)) {
    Z <- matrix(rnorm(nrow(X)*q, sd=0.2), ncol=q)
  } else if (identical(Z.init, "PCA")) {
    X.pca <- prcomp(X)
    Z <- X.pca$x[,1:q]
  } else {
    if (ncol(as.matrix(Z.init)) != q) warning("Mismatch between Z.init and q")
    Z <- as.matrix(Z.init)[,1:q]
  }

  par <- c(1, 1, 0.1)


  for (it in 1:iterations) {
    # First optimize the hyper parameters
    kernel.function <- function(x1, x2) gplvm.SE.dist(sum((x1-x2)^2), par[1], par[2], 0)
    active.set <- IVM::IVM.regression(Z, active.set.size, par[3], kernel.function)$activeSet

    par.optout <- optimx::optimx(par,
                                 gplvm.hp.f,
                                 gplvm.hp.gr,
                                 Z=Z[active.set, ],
                                 X=X[active.set, ],
                                 Z.normal.prior=Z.normal.prior,
                                 method="L-BFGS-B",
                                 control=list(trace=T,
                                              maximize=T,
                                              kkt=FALSE,
                                              maxit=1000,
                                              starttests=TRUE))
    par <- par.optout[1:3]

    # Now select a new active set and optimize the latent variables
    kernel.function <- function(x1, x2) gplvm.SE.dist(sum((x1-x2)^2), par[1], par[2], 0)
    active.set <- IVM::IVM.regression(Z, active.set.size, par[3], kernel.function)$activeSet

    K.active <- gplvm.SE(Z[active.set, ], par[1], par[2], 0)
    K.chol.active <- chol(K.active)



  }
}
