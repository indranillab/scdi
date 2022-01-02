

gplvm.L.from.K <- function(Z1, X, l,alpha,sigma,Kinv,Kvecmat, Z.normal.prior=FALSE) {
  N <- nrow(X)
  m1<-Z1[(dim(Z1)[1]),]
  matt<-sweep(Z1,2,m1,"-")
  mattsum<-rowSums(matt^2)
  kvec<-alpha^2*exp(-mattsum/(2*l^2))
  kvec[N]<-kvec[N]+sigma^2
  X <- as.matrix(X)
  D <- ncol(X)
  #K.chol <- chol(K)
  #K.inv.X <- Kinv%*%X
  #K.log.det <-  sum(log(diag(Kinv)))
  Z.prior.term <- 0
  #if (Z.normal.prior) {
   # Z.prior.term <- - length(Z)/2 * log(2 * 10^2 * pi) - 1 / (2 * 10 ^2) * sum(as.numeric(Z)^2)
  #}
  rr<-sum(t(Kvecmat)%*%kvec[-N])
  rr1<-sum(t(Kvecmat)%*%(kvec[-N]*X[N,-N]))
  return(-D * N * log(2 * pi) / 2 - D * log(kvec[N]-quad.form(Kinv,kvec[-N]))/ 2 - 1/2 * (rr^2/(kvec[N]-quad.form(Kinv,kvec[-N]))-2*rr1/(kvec[N]-quad.form(Kinv,kvec[-N]))+sum((X[N,N])^2)/(kvec[N]-quad.form(Kinv,kvec[-N]))) + Z.prior.term)
}



gplvm.L <- function(Z1, X, l, alpha, sigma,Kinv,Kvecmat,Z.normal.prior=TRUE) {
  return(gplvm.L.from.K(Z1, X, l,alpha,sigma,Kinv,Kvecmat, Z.normal.prior=FALSE))
}




likelihood<-function(Znew,X,Z,Kinv,Kvecmat,l,alpha,sigma)
{
	Z1<-rbind(Z,Znew)
	return(gplvm.L(Z1,X,l,alpha,sigma,Kinv,Kvecmat))
}

likelihood.gr<-function(Znew,X,Z,Kinv,Kvecmat,l,alpha,sigma)
{
	
}



add_gplvm<-function(X,Z,Kinv,Kvecmat,l,alpha,sigma)
{
   Znew<-rep(0,dim(Z)[2])
	optout<-optim(Znew,
	               likelihood,                   
	               X=X,
                   Z=Z,
                   Kinv=Kinv,
                   Kvecmat=Kvecmat,
                   l=l,
                   alpha=alpha,
                   sigma=sigma, 
                   method="L-BFGS-B",
                   control=list(maxit=10))
                   
    return(optout$par)  
	
	
}




