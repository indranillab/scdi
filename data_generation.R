k<-2
M<-100
N<-50


set.seed(0)
mu1<-rnorm(k,0,10)
mu2<-rnorm(k,0,10)
mu3<-rnorm(k,0,10)
batch1<-rnorm(M,0,0)
batch2<-rnorm(M,0,0)
data1<-array(0,dim<-c(N,k))
data2<-array(0,dim<-c(N,k))
data3<-array(0,dim<-c(N,k))
data4<-array(0,dim<-c(N,M))
data5<-array(0,dim<-c(N,M))
data6<-array(0,dim<-c(N,M))
data7<-array(0,dim<-c(N,k))
data8<-array(0,dim<-c(N,k))
data9<-array(0,dim<-c(N,k))
data10<-array(0,dim<-c(N,M))
data11<-array(0,dim<-c(N,M))
data12<-array(0,dim<-c(N,M))


library(mvtnorm)
data1<-rmvnorm(N,mu1,diag(k))
data2<-rmvnorm(N,mu2,diag(k))
data3<-rmvnorm(N,mu3,diag(k))
data7<-rmvnorm(N,mu1,diag(k))
data8<-rmvnorm(N,mu2,diag(k))
data9<-rmvnorm(N,mu3,diag(k))

datalatent<-rbind(data1,data2,data3,data7,data8,data9)

	W1<-array(0,dim<-c(M,2))
	W2<-array(0,dim<-c(M,2))
	W3<-array(0,dim<-c(M,2))
	for(j in 1:M)
	{
		for(l in 1:k)
		{
			W1[j,l]<-rnorm(1,0,0.5)
			W2[j,l]<-rnorm(1,0,0.5)
			W3[j,l]<-rnorm(1,0,0.5)
		}
	}

for(i in 1:N)
{
	data4[i,]<-W1%*%data1[i,]+rnorm(M,0,1)+batch1
	data5[i,]<-W2%*%data2[i,]+rnorm(M,0,1)+batch1
	data6[i,]<-W3%*%data3[i,]+rnorm(M,0,1)+batch1
	
}

for(i in 1:N)
{
	data10[i,]<-W1%*%data7[i,]+rnorm(M,0,1)+batch2
	data11[i,]<-W2%*%data8[i,]+rnorm(M,0,1)+batch2
	data12[i,]<-W3%*%data9[i,]+rnorm(M,0,1)+batch2
	
}

data<-rbind(data4,data5,data6,data10,data11,data12)

library(Rtsne)
xx<-Rtsne(data)


N<-50
cols<-c(rep("red",N),rep("black",N),rep("blue",N),rep("green",N),rep("purple",N),rep("cyan",N))


cols<-c(rep("red",10),rep("black",N),rep("blue",90),rep("green",N),rep("purple",N),rep("cyan",N))



source("tsne_internal.R")


initial_config = NULL
k=2
initial_dims=30
perplexity=30
 max_iter = 1000
 min_cost=0
 epoch_callback=NULL
 whiten=TRUE
 epoch=100 
X<-data


group<-c(rep(0,150),rep(1,150))



pdf("batch_effect_tsne.pdf")
plot(xx$Y,col=cols)
dev.off()



        if ('dist' %in% class(X)) {
		n = attr(X,'Size')
		}else{
		X = as.matrix(X)
		X = X - min(X)
		X = X/max(X)
		initial_dims = min(initial_dims,ncol(X))
		if (whiten){
		X1<-.whiten(as.matrix(X[(1:(dim(X)[1]/2)),]),n.comp=initial_dims)	
		X2<-.whiten(as.matrix(X[((dim(X)[1]/2)+1):(dim(X)[1]),]),n.comp=initial_dims)	
		X<-rbind(X1,X2)
		} 
		n = nrow(X)
	}
	
	matt<-X%*%t(X)

	momentum = .5
	final_momentum = .8
	mom_switch_iter = 250

	epsilon = 500
	min_gain = .01
	initial_P_gain = 4

	eps = 2^(-52) # typical machine precision

	if (!is.null(initial_config) && is.matrix(initial_config)) {
		if (nrow(initial_config) != n | ncol(initial_config) != k){
			stop('initial_config argument does not match necessary configuration for X')
		}
		ydata = initial_config
		initial_P_gain = 1

	} else {
		ydata = matrix(rnorm(k * n),n)
	}

	#P = .x2p(X,perplexity, 1e-5)$P
	#P = .5 * (P + t(P))

	#P[P < eps]<-eps
	#P = P/sum(P)

	#P = P * initial_P_gain
	grads =  matrix(0,nrow(ydata),ncol(ydata))
	incs =  matrix(0,nrow(ydata),ncol(ydata))
	gains = matrix(1,nrow(ydata),ncol(ydata))

	for (iter in 1:max_iter){
		print(iter)
		if (iter %% epoch == 0) { # epoch
			#cost =  sum(apply(P * log((P+eps)/(Q+eps)),1,sum))
			#message("Epoch: Iteration #",iter," error is: ",cost)
			#if (cost < min_cost) break
			if (!is.null(epoch_callback)) epoch_callback(ydata)

		}

		#sum_ydata = apply(ydata^2, 1, sum)
		#num =  1/(1 + sum_ydata +    sweep(-2 * ydata %*% t(ydata),2, -t(sum_ydata)))
		#diag(num)=0
		#Q = num / sum(num)
		#if (any(is.nan(num))) message ('NaN in grad. descent')
		#Q[Q < eps] = eps
		#stiffnesses = 4 * (P-Q) * num
		#for (i in 1:n){
		#	grads[i,] = apply(sweep(-ydata, 2, -ydata[i,]) * stiffnesses[,i],2,sum)
		#}
		one1<-c(rep(1,(dim(X)[1]/2)),rep(0,(dim(X)[1]/2)))
		one2<-c(rep(0,(dim(X)[1]/2)),rep(1,(dim(X)[1]/2)))
		Sigma<-0.5*diag((dim(ydata)[1]))+0.5*one1%*%t(one1)+0.5*one2%*%t(one2)
		
		K<-ydata%*%t(ydata)+diag((dim(ydata)[1]))+Sigma
		K1<-solve(K)
		grads<--(K1%*%matt%*%K1%*%ydata-(dim(X)[2])*K1%*%ydata)

		gains = ((gains + .2) * abs(sign(grads) != sign(incs)) +
				 gains * .8 * abs(sign(grads) == sign(incs)))

		gains[gains < min_gain] = min_gain

		incs = momentum * incs - epsilon * (gains * grads)
		ydata = ydata + incs
		ydata = sweep(ydata,2,apply(ydata,2,mean))
		plot(ydata,col=cols)
		if (iter == mom_switch_iter) momentum = final_momentum

		#if (iter == 100 && is.null(initial_config)) P = P/4

	}
	
	
	
pdf("batch_effect_gplvm.pdf")
plot(ydata,col=cols)
dev.off()



fn<-gplvm.f
gr<-gplvm.gr
hess<-NULL
lower=-Inf
upper=Inf
method<-"L-BFGS-B"
itnmax=NULL
hessian<-FALSE
control=list(trace=T,
                                maximize=T,
                                kkt=FALSE,
                                maxit=10,
                                starttests=TRUE)
optcfg <- optimx.setup(par, fn, gr, hess, lower, upper, method, itnmax, hessian, control=control,X=X,q=q,groups=groups,KK=KK,mumat=mumat,sigmat=sigmat,pii=pii,l=l,alpha=alpha,sigma=sigma,sigma1=sigma1)


optimx.check(par, ufn, ugr, uhess, lower=-Inf, upper=Inf, 

hessian=FALSE, ctrl, have.bounds=FALSE, usenumDeriv=FALSE) 



 optchk <- optimx.check(par, optcfg$ufn, optcfg$ugr, optcfg$uhess, 
            lower, upper, hessian, optcfg$ctrl, have.bounds = optcfg$have.bounds, 
            usenumDeriv = optcfg$usenumDeriv)


