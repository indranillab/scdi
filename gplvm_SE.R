
gplvm.SE.dist <- function(dist.matrix, l, alpha, sigma=0) {

  K <- alpha^2 * exp(-dist.matrix^2/(2*l^2))

  diag(K) <- diag(K) + sigma^2

  return(K)

}





gplvm.SE <- function(Z, l, alpha, sigma=0) {

  gplvm.SE.dist(dist.matrix=as.matrix(dist(Z)),

                l=l,

                alpha=alpha,

                sigma=sigma)

}

