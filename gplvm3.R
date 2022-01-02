tsne <-
function(X, initial_config = NULL, k=2, initial_dims=30, perplexity=30, max_iter = 1000, min_cost=0, epoch_callback=NULL,whiten=TRUE, epoch=100 ){
		Xnew<-array(0,dim<-dim(X))
		for(le in group)
		{
			X1<-X[(group==le),]
	if ('dist' %in% class(X1)) {
		n = attr(X1,'Size')
		}else{
		X1 = as.matrix(X1)
		X1 = X1 - min(X1)
		X1 = X1/max(X1)
		initial_dims = min(initial_dims,ncol(X1))
		if (whiten) X1<-.whiten(as.matrix(X1),n.comp=initial_dims)
		n = nrow(X1)
	}
	Xnew[(group==le),]<-X1
	}
	
	matt<-X%*%t(X)
	XX<-array(0,dim<-c(dim(X)[1],length(table(group))))
	counter<-1
	for(le in (unique(group)))
	{
		XX[,counter]<-as.numeric(group==le)
		counter<-counter+1
	}
	matt1<-XX%*%t(XX)

	momentum = .5
	final_momentum = .8
	mom_switch_iter = 250

	epsilon = 500
	epsilon1<-100
	epsilon2<-100
	epsilon3<-100
	min_gain = .01
	initial_P_gain = 4
	alpha<-1
	beta<-1
	gamma<-1

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
	grads1<-0
	incs1<-0
	gains1<-1
	grads2<-0
	incs2<-0
	gains2<-1
	grads3<-0
	incs3<-0
	gains3<-1

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
		K<-alpha*matt+(1/beta)*diag((dim(ydata)[1]))+gamma*matt1
		K1<-solve(K)
		grads<--(K1%*%matt%*%K1%*%ydata-(dim(X)[2])*K1%*%ydata)
		grads1<-(dim(X)[2])*sum(diag(K1%*%(matt)))/2-sum(diag(K1%*%ydata%*%t(ydata)%*%K1%*%matt))/2
		grads2<--(dim(X)[2])*sum(diag(K1))/(2*beta^2)+sum(diag(K1%*%ydata%*%t(ydata)%*%K1))/(2*beta^2)
		grads3<-(dim(X)[2])*sum(diag(K1%*%(matt1)))/2-sum(diag(K1%*%ydata%*%t(ydata)%*%K1%*%matt1))/2

		gains = ((gains * .2) * abs(sign(grads) != sign(incs)) +
				 gains * .8 * abs(sign(grads) == sign(incs)))
		gains1 = ((gains1 * .2) * abs(sign(grads1) != sign(incs1)) +
				 gains1 * .8 * abs(sign(grads1) == sign(incs1)))
		gains2 = ((gains2 * .2) * abs(sign(grads2) != sign(incs2)) +
				 gains2 * .8 * abs(sign(grads2) == sign(incs2)))
		gains3 = ((gains3 * .2) * abs(sign(grads3) != sign(incs3)) +
				 gains3 * .8 * abs(sign(grads3) == sign(incs3)))

		gains[gains < min_gain] = min_gain

		incs = momentum * incs - epsilon * (gains * grads)
		incs1 = momentum * incs1 - epsilon1 * (gains1 * grads1)
		incs2 = momentum * incs2 - epsilon2 * (gains2 * grads2)
		incs3 = momentum * incs3 - epsilon3 * (gains3 * grads3)
		ydata = ydata + incs
		ydata = sweep(ydata,2,apply(ydata,2,mean))
		alpha<-alpha+incs1
		beta<-beta+incs2
		gamma<-gamma+incs3
		if(alpha<0.01)
		{
			alpha<-0.01
		}
		if(beta<0.01)
		{
			beta<-0.01
		}
		if(gamma<0.01)
		{
			gamma<-0.01
		}
		if(alpha>100)
		{
			alpha<-100
		}
		if(beta>100)
		{
			beta<-100
		}
		if(gamma>100)
		{
			gamma<-100
		}
		plot(ydata,col=cols)
		if (iter == mom_switch_iter) momentum = final_momentum

		#if (iter == 100 && is.null(initial_config)) P = P/4

	}
	
	
	
	lst<-list(ydata=ydata,alpha=alpha,beta=beta)
	return(lst)
}





