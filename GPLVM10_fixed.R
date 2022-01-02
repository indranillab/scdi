
dL.dK <- function(X, K, groups) {
  K.chol <- chol(K)
  K.inv <- chol2inv(K.chol)
  temp <- backsolve(K.chol, forwardsolve(t(K.chol), X))
  D <- ncol(X)
  (temp %*% t(temp) - D * K.inv) / 2
}


dL.dZ <- function(X, Z, l, alpha, sigma, sigma1, K, dL.dK, Z.normal.prior, groups,b) {
  out <- Z
    Z1<-Z
  for(i in 1:(dim(Z)[1]))
  {
  	if(groups[i]==1)
  	{
  		Z1[i,]<-Z[i,]+b
  	}else
  	{
  		Z1[i,]<-Z[i,]-b
  	}
  }
  for (i in 1:nrow(Z)) {
    for (j in 1:ncol(Z)) {
      out[i, j] <- sum(dL.dK * dK.dZij(Z1, K, i, j, l,b,groups))
    }
  }
  if (Z.normal.prior) {
    out <- out - Z / 1^2
  }
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

dL.db<-function(X, Z, l, alpha, sigma, sigma1, K, dL.dK, groups,b)
{
    vec<-NULL
    for(i in 1:(length(b)))
    {
    	vec[i]<-sum(dL.dK*dK.dbi(Z,K,i,l,b,groups))
    }
    return(vec)
}
 



dL.dsigma1 <- function(X, Z, l, alpha, sigma, sigma1, K, dL.dK, groups) {
	XX<-array(0,dim<-c(dim(Z)[1],length(table(groups))))
	counter<-1
	for(le in (unique(groups)))
	{
		XX[,counter]<-as.numeric(groups==le)
		counter<-counter+1
	}
    dK.dsigma1<--2*sigma1*(1-XX%*%t(XX))*K

  out <- sum(dL.dK * dK.dsigma1)
  return(out)
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
    Z.prior.term <- - length(Z)/2 * log(2 * 1^2 * pi) - 1 / (2 * 1 ^2) * sum(as.numeric(Z)^2)
  }
  return(-D * N * log(2 * pi) / 2 - D * K.log.det / 2 - 1/2 * sum(K.inv.X * X) + Z.prior.term)
}



gplvm.L <- function(Z, X, l, alpha, sigma,sigma1, Z.normal.prior=TRUE,groups,b) {
  K <- gplvm.SE(Z, l, alpha, sigma,sigma1,groups,b)
  return(gplvm.L.from.K(Z, X, K, Z.normal.prior))
}



gplvm.f <- function(par, X,groups, Z.normal.prior) {
  l <- par[1]
  alpha <- par[2]
  sigma <- par[3]
  sigma1 <- par[4]
  b<-par[(5:(4+q))]
  Z <- matrix(par[-(1:(4+q))], nrow=nrow(X))
  gplvm.L(Z, X, l, alpha, sigma, sigma1, Z.normal.prior,groups,b)
}



gplvm.gr <- function(par, X,groups, Z.normal.prior) {
  l <- par[1]
  alpha <- par[2]
  sigma <- par[3]
  sigma1<-par[4]
  b<-par[(5:(4+q))]
  Z <- matrix(par[-(1:(4+q))], nrow=nrow(X))
  K <- gplvm.SE(Z, l, alpha, sigma,sigma1,groups,b)
  dL.dK <- dL.dK(X, K)
  Z1<-Z
  for(i in 1:(dim(Z)[1]))
  {
  	if(groups[i]==1)
  	{
  		Z1[i,]<-Z[i,]+b
  	}else
  	{
  		Z1[i,]<-Z[i,]-b
  	}
  }
  c(dL.dl(X, Z1, l, alpha, sigma, sigma1, K, dL.dK,groups),
    dL.dalpha(X, Z1, l, alpha, sigma, sigma1, K, dL.dK,groups),
    dL.dsigma(X, Z1, l, alpha, sigma, sigma1, K, dL.dK,groups),
    dL.dsigma1(X, Z1, l, alpha, sigma, sigma1, K, dL.dK, groups),
    dL.db(X, Z, l, alpha, sigma, sigma1, K, dL.dK, groups,b),
    as.numeric(dL.dZ(X, Z, l, alpha, sigma, sigma1, K, dL.dK, Z.normal.prior,groups,b))
  )
}


gplvm.hp.f <- function(par, Z, X,groups, Z.normal.prior) {
  l <- par[1]
  alpha <- par[2]
  sigma <- par[3]
  sigma1<-par[4]
  b<-par[5:(4+dim(Z)[2])]
  gplvm.L(Z, X, l, alpha, sigma,sigma1, Z.normal.prior=Z.normal.prior,groups,b)
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
                      iterations=1000,
                      plot.freq=100,
                      classes=NULL,
                      Z.init="PCA",
                      num.init.params=100,
                      Z.normal.prior,
                      groups,
                      sigma1_init) {
		#Xnew<-array(0,dim<-dim(X))
		#for(le in groups)
		#{
		#	X1<-X[(groups==le),]
#	if ('dist' %in% class(X1)) {
#		n = attr(X1,'Size')
#		}else{
#		X1 = as.matrix(X1)
#		X1 = X1 - min(X1)
#		X1 = X1/max(X1)
#		initial_dims = ncol(X1)
#	    X1<- .whiten(as.matrix(X1),n.comp=initial_dims)
#	}
    #X1<-t(t(X1)-colMeans(X1))
	#Xnew[(groups==le),]<-X1
	#}
	#X<-Xnew
                      	
                      	
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

  if (q > 2) {
    pairs(Z, col=classes)
  } else if (q == 2) {
    plot(Z, col=classes)
  } else {
    plot(as.numeric(Z), rep(0, length(Z)), col=classes)
  }


  par.init <- c(1, 1, 1,1,rep(0,(dim(Z)[2])))
  par.init.L <- gplvm.hp.f(par.init, Z=Z, X=X,groups, Z.normal.prior=Z.normal.prior)
  #first optimize the hyperparams for the initial Z
  for (i in seq(num.init.params)) {
    test.par.init <- c(rgamma(4, 1),rep(0,q))
    test.par.init.L <- gplvm.hp.f(test.par.init, Z=Z, X=X,groups, Z.normal.prior=Z.normal.prior)
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
                   Z.normal.prior=Z.normal.prior,
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
  #init.sigma1<-sigma1_init

  cat("Starting params: l=", init.l, "; alpha=", init.alpha, "; sigma=", init.sigma, fill=TRUE)
  
  init.b<-as.numeric(optout[(5:(4+q))])

  par <- c(init.l, init.alpha, init.sigma, init.sigma1, init.b, as.numeric(Z))

  its <- 0
  convcode=1
  while (its < iterations & convcode != 0) {
    out <- optimx(par,
                  gplvm.f,
                  gplvm.gr,
                  X=X,
                  groups=groups,
                  Z.normal.prior=Z.normal.prior,
                  method="L-BFGS-B",
                  control=list(trace=T,
                               maximize=T,
                               kkt=FALSE,
                               maxit=plot.freq,
                               starttests=FALSE)
    )
    convcode <- out$convcode
    par <- head(as.numeric(out), length(Z) + 4+q)
    l <- par[1]
    alpha <- par[2]
    sigma <- par[3]
    sigma1<-par[4]
    b<-par[(5:(4+q))]
    Z <- matrix(par[-(1:(4+q))], ncol=q)
    if (q > 2) {
      #pairs(Z,
      #      col=classes,
      #      main=paste("l=", signif(l, 3), "; alpha=", signif(alpha, 3), "; sigma=", signif(sigma, 3), sep=""))
    } else if (q == 2) {
      #plot(Z,
      #     col=classes,
      #     main=paste("l=", signif(l, 3), "; alpha=", signif(alpha, 3), "; sigma=", signif(sigma, 3), sep=""))
    } else {
      #plot(as.numeric(Z),
      #     rep(0, length(Z)),
      #     col=classes,
      #     main=paste("l=", signif(l, 3), "; alpha=", signif(alpha, 3), "; sigma=", signif(sigma, 3), sep=""))
    }

    its <- its + plot.freq
  }

  return(list(Z=Z, l=l, alpha=alpha, sigma=sigma, convcode=convcode,sigma1=sigma1,b=b))

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
    kernel.function <- function(x1, x2) gplvm.SE.dist(sum((x1-x2)^2), par[1], par[2], 0,rep(0,(dim(Z)[2])))
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
    kernel.function <- function(x1, x2) gplvm.SE.dist(sum((x1-x2)^2), par[1], par[2], 0,rep(0,(dim(Z)[2])))
    active.set <- IVM::IVM.regression(Z, active.set.size, par[3], kernel.function)$activeSet

    K.active <- gplvm.SE(Z[active.set, ], par[1], par[2], 0)
    K.chol.active <- chol(K.active)



  }
}
