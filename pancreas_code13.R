a1<-read.csv("GSM2230757_human1_umifm_counts.csv")
a2<-read.csv("GSM2230758_human2_umifm_counts.csv")

info1<-a1[,1:3]
info2<-a2[,1:3]

a1<-a1[,-(1:3)]
a2<-a2[,-(1:3)]

set.seed(12345)

set1<-sample(1937,100)
set2<-sample(1724,100)

a3<-a1[set1,]
a4<-a2[set2,]

data1<-as.matrix(a3)
data2<-as.matrix(a4)

colnames(data1)<-colnames(a3)
colnames(data2)<-colnames(a4)


data3<-array(0,dim<-dim(data1))
data4<-array(0,dim<-dim(data2))

for(i in 1:(dim(data1)[1]))
{
	data3[i,]<-as.numeric(data1[i,])
}

for(i in 1:(dim(data2)[1]))
{
	data4[i,]<-as.numeric(data2[i,])
}

for(i in 1:(dim(data3)[1]))
{
	data3[i,]<-data3[i,]/sum(data3[i,])*1000000
}

for(i in 1:(dim(data4)[1]))
{
	data4[i,]<-data4[i,]/sum(data4[i,])*1000000
}

a5<-a1[(sample(1937,1)),]
a5<-as.numeric(a5)
a5<-a5/sum(a5)*1000000


colnames(data3)<-colnames(a3)
colnames(data4)<-colnames(a4)

data1<-data3
data2<-data4

for(i in 1:(dim(data1)[1]))
{
	data1[i,]<-data1[i,]/sqrt(sum((data1[i,])^2))
}
for(i in 1:(dim(data2)[1]))
{
	data2[i,]<-data2[i,]/sqrt(sum((data2[i,])^2))
}


data<-rbind(data1,data2)

groups<-c(rep(0,100),rep(1,100))


save.image("pamcreas1.Rdata")

source("scdi18.R")
K1<-9
K2<-9


Z2<-scdi_integrate(data,groups,K1,K2)

save.image("pancreas2.Rdata")

set3<-setdiff(1:1937,set1)
set4<-setdiff(1:1724,set2)


l<-Z2$l
alpha<-Z2$alpha
sigma<-Z2$sigma
sigma1<-Z2$sigma1
Z3<-Z2$Z2
groups<-Z2$groups
cluster<-Z2$cluster
dataa<-Z2$dataa

library(GPLVM)
X1<-t(chol(dataa[(groups==0),]%*%t(dataa[(groups==0),])))
X2<-t(chol(dataa[(groups==1),]%*%t(dataa[(groups==1),])))
X3<-forwardsolve(X1,diag(length(diag(X1))))
X4<-forwardsolve(X2,diag(length(diag(X2))))

source("gplvm_SE.R")
K<-gplvm.SE(Z3[(groups==0),],l,alpha,sigma)

Kinv<-solve(K)
Kvecmat<-Kinv%*%dataa[(groups==0),]

library(emulator)

Z2Z2<-array(0,dim<-c(length(set3),dim(Z3)[2]))

inndex<-1

dataa1<-dataa[(groups==0),]
source("add_gplvm.R")

cost<-function(y,d,X1)
{
	s<-0
	for(i in 1:(dim(X1)[1]))
	{
		s<-s+(d[i]-sum((y-X1[i,])^2))^2
	}
	return(s)
}

for(i in 1:(dim(dataa1)[1]))
{
	dataa1[i,]<-dataa1[i,]/sqrt(sum((dataa1[i,])^2))
}

cl1<-Z2$cluster[(Z2$groups==0)]

clustnew11<-NULL

cormat1<-array(0,dim<-c((dim(dataa1)[1]),length(set3)))
for(i in set3)
{
	a5<-a1[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa1,a5)
	corvec<-cov(a5,t(dataa1))*nrow(data1)
	disclust<-NULL
	for(j in 1:(length(unique(cl1))))
	{
		disclust[j]<-mean(corvec[cl1==j])
	}
	kk<-which(disclust==max(disclust))
	corvec1<-corvec[cl1==kk]
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z4<-Z3[(groups==0),]
    Z5<-Z4[(cl1==kk),]
    y<-colMeans(Z5)
    opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    Z2Z2[inndex,]<-opt$par
    clustnew11[inndex]<-kk
    cormat1[,inndex]<-corvec
    inndex<-inndex+1	
    print(inndex)
}




Z4Z4<-array(0,dim<-c(length(set4),dim(Z3)[2]))

inndex<-1

dataa3<-dataa[(groups==1),]
source("add_gplvm.R")

cost<-function(y,d,X1)
{
	s<-0
	for(i in 1:(dim(X1)[1]))
	{
		s<-s+(d[i]-sum((y-X1[i,])^2))^2
	}
	return(s)
}

for(i in 1:(dim(dataa3)[1]))
{
	dataa3[i,]<-dataa3[i,]/sqrt(sum((dataa3[i,])^2))
}

cl1<-Z2$cluster[(Z2$groups==1)]

clustnew12<-NULL

cormat2<-array(0,dim<-c((dim(dataa3)[1]),length(set4)))
for(i in set4)
{
	a5<-a2[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa3,a5)
	corvec<-cov(a5,t(dataa3))*nrow(data2)
	disclust<-NULL
	for(j in 1:(length(unique(cl1))))
	{
		disclust[j]<-mean(corvec[cl1==j])
	}
	kk<-which(disclust==max(disclust))
	corvec1<-corvec[cl1==kk]
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z4<-Z3[(groups==1),]
    Z5<-Z4[(cl1==kk),]
    y<-colMeans(Z5)
    opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    Z4Z4[inndex,]<-opt$par
    clustnew12[inndex]<-kk
    cormat2[,inndex]<-corvec
    inndex<-inndex+1	
    print(inndex)
}



Z6Z6<-array(0,dim<-c(length(set4),dim(Z3)[2]))

inndex<-1

dataa1<-dataa[(groups==0),]
source("add_gplvm.R")

cost<-function(y,d,X1)
{
	s<-0
	for(i in 1:(dim(X1)[1]))
	{
		s<-s+(d[i]-sum((y-X1[i,])^2))^2
	}
	return(s)
}

for(i in 1:(dim(dataa1)[1]))
{
	dataa1[i,]<-dataa1[i,]/sqrt(sum((dataa1[i,])^2))
}

cl1<-Z2$cluster[(Z2$groups==0)]

clustnew21<-NULL

cormat3<-array(0,dim<-c((dim(dataa1)[1]),length(set4)))
for(i in set4)
{
	a5<-a2[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa1,a5)
	corvec<-cov(a5,t(dataa1))*nrow(data1)
	disclust<-NULL
	for(j in 1:(length(unique(cl1))))
	{
		disclust[j]<-mean(corvec[cl1==j])
	}
	kk<-which(disclust==max(disclust))
	corvec1<-corvec[cl1==kk]
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z4<-Z3[(groups==0),]
    Z5<-Z4[(cl1==kk),]
    y<-colMeans(Z5)
    opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    Z6Z6[inndex,]<-opt$par
    clustnew21[inndex]<-kk
    cormat3[,inndex]<-corvec
    inndex<-inndex+1	
    print(inndex)
}




Z8Z8<-array(0,dim<-c(length(set3),dim(Z3)[2]))

inndex<-1

dataa3<-dataa[(groups==1),]
source("add_gplvm.R")

cost<-function(y,d,X1)
{
	s<-0
	for(i in 1:(dim(X1)[1]))
	{
		s<-s+(d[i]-sum((y-X1[i,])^2))^2
	}
	return(s)
}

for(i in 1:(dim(dataa3)[1]))
{
	dataa3[i,]<-dataa3[i,]/sqrt(sum((dataa3[i,])^2))
}

cl1<-Z2$cluster[(Z2$groups==1)]

clustnew22<-NULL

cormat4<-array(0,dim<-c((dim(dataa3)[1]),length(set3)))
for(i in set3)
{
	a5<-a1[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa3,a5)
	corvec<-cov(a5,t(dataa3))*nrow(data2)
	disclust<-NULL
	for(j in 1:(length(unique(cl1))))
	{
		disclust[j]<-mean(corvec[cl1==j])
	}
	kk<-which(disclust==max(disclust))
	corvec1<-corvec[cl1==kk]
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z4<-Z3[(groups==1),]
    Z5<-Z4[(cl1==kk),]
    y<-colMeans(Z5)
    opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    Z8Z8[inndex,]<-opt$par
    clustnew22[inndex]<-kk
    cormat4[,inndex]<-corvec
    inndex<-inndex+1	
    print(inndex)
}


save.image("pancreas3.Rdata")




density_simulate<-function(ds,N)
{
	weights<-(ds$y[-1]+ds$y[-(length(ds$y))])/2
	samp<-sample(1:length(weights),N,prob=weights,replace=TRUE)
	return((runif(N)*ds$bw+ds$x[samp]))
	
}




set.seed(12345)

clus1<-Z2$cluster[(Z2$groups==0)]
clus2<-Z2$cluster[(Z2$groups==1)]

corrmat<-array(0,dim<-c(length(c(clus1,clus2)),length(c(set3,set4))))

corrmat[(1:(length(clus1))),(1:(length(set3)))]<-cormat1
corrmat[((length(clus1)+1):(length(clus1)+length(clus2))),((length(set3)+1):(length(set3)+length(set4)))]<-cormat2

for(i in 1:(length(table(clus1))))
{
	ind1<-which(clus1==i)
	ind2<-which(clustnew11==i)
	ind3<-which(clustnew12==i)
	for(k1 in ind1)
	{
		vec<-as.numeric(cormat1[k1,ind1])
		ds<-density(vec,from=0,kernel="rectangular")
		rand<-density_simulate(ds,length(ind3))
		corrmat[k1,(length(set3)+ind3[order(cormat3[k1,ind3])])]<-sort(rand)
		rand1<-density_simulate(ds,length(ind2))
		corrmat[k1,ind2[order(cormat1[k1,ind2])]]<-sort(rand1)
	}
	
	
}

for(i in 1:(length(table(clus2))))
{
	ind1<-which(clus2==i)
	ind2<-which(clustnew12==i)
	ind3<-which(clustnew11==i)
	for(k1 in ind1)
	{
		vec<-as.numeric(cormat2[k1,(ind2)])
		ds<-density(vec,from=0,kernel="rectangular")
		rand<-density_simulate(ds,length(ind3))
		corrmat[(length(clus1)+k1),ind3[order(cormat4[k1,ind3])]]<-sort(rand)
		rand1<-density_simulate(ds,length(ind2))
		corrmat[(length(clus1)+k1),ind2[order(cormat1[k1,ind2])]]<-sort(rand1)
	}
	
	
}

clus<-c(clus1,clus2)
Zmat<-array(0,dim<-c(length(c(set3,set4)),(dim(X3)[2])))
inndex<-1
clustn<-c(clustnew11,clustnew12)
for(j in 1:(dim(corrmat)[2]))
{
    corvec<-corrmat[,j]	
	corvec1<-corvec[clus==clustn[inndex]]
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z5<-X3[(clus==clustn[inndex]),]
    y<-colMeans(Z5)
    opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    Zmat[inndex,]<-opt$par
    inndex<-inndex+1	
    print(inndex)
}


save.image("pancreas4.Rdata")






change_sclae<-function(x,m1,m2,sig1,sig2)
{
	ei1<-eigen(sig1)
	ei2<-eigen(sig2)
	L1<-max(ei1$values)
	L2<-max(ei2$values)
	siginv1<-(ei1$vectors[,((ei1$values/L1)>0.00000001)])%*%diag(1/sqrt(ei1$values[((ei1$values/L1)>0.00000001)]))%*%t((ei1$vectors[,((ei1$values/L1)>0.00000001)]))
	sigsq2<-(ei2$vectors[,((ei2$values/L2)>0.00000001)])%*%diag(sqrt(ei2$values[((ei2$values/L2)>0.00000001)]))%*%t((ei2$vectors[,((ei2$values/L2)>0.00000001)]))
	y<-sigsq2%*%siginv1%*%(x+m2-m1)
	#y<-x+m2-m1
	return(y)
}


Z4<-Z3
for(i in 1:(length(table(clus))))
{
	m1<-colMeans(Z3[(clus==i),])
	m2<-colMeans(Zmat[(clustn==i),])
	sig1<-var(Z3[(clus==i),])
	sig2<-var(Zmat[(clustn==i),])
	Z5<-Z3[(clus==i),]
	Z6<-Z5
	for(j in 1:(dim(Z5)[1]))
	{
		Z6[j,]<-change_sclae(Z5[j,],m1,m2,sig1,sig2)
	}
	Z4[(clus==i),]<-Z6
	print(i)
}






dataa2<-a1[set3,]

dataa2<-as.matrix(dataa2)
norm<-sqrt(rowSums(dataa2^2))

for(i in 1:(dim(dataa2)[1]))
{
	dataa2[i,]<-dataa2[i,]/norm[i]
}

Z3Z3<-array(0,dim<-c(length(set1),dim(Z3)[2]))

inndex<-1

clusterold1<-NULL

for(i in set1)
{
	a5<-a1[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa2,a5)
	corvec<-cov(a5,t(dataa2))*nrow(data1)
	disclust<-NULL
	for(j in 1:(length(unique(clustnew11))))
	{
		disclust[j]<-mean(corvec[clustnew11==j])
	}
	kk<-which(disclust==max(disclust))
	corvec1<-corvec[clustnew11==kk]
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z4<-Z2Z2
    Z5<-Z4[(clustnew11==kk),]
    y<-colMeans(Z5)
    opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    Z3Z3[inndex,]<-opt$par
    clusterold1[inndex]<-kk
    inndex<-inndex+1	
    print(inndex)
}






dataa4<-a2[set4,]

dataa4<-as.matrix(dataa4)
norm<-sqrt(rowSums(dataa4^2))

for(i in 1:(dim(dataa4)[1]))
{
	dataa4[i,]<-dataa4[i,]/norm[i]
}


Z5Z5<-array(0,dim<-c(length(set2),dim(Z3)[2]))

inndex<-1

clusterold2<-NULL

for(i in set2)
{
	a5<-a2[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa4,a5)
	corvec<-cor(a5,t(dataa4))*nrow(data2)
	disclust<-NULL
	for(j in 1:(length(unique(clustnew12))))
	{
		disclust[j]<-mean(corvec[clustnew12==j])
	}
	kk<-which(disclust==max(disclust))
	corvec1<-corvec[clustnew12==kk]
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z4<-Z4Z4
    Z5<-Z4[(clustnew12==kk),]
    y<-colMeans(Z5)
    opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    Z5Z5[inndex,]<-opt$par
    clusterold2[inndex]<-kk
    inndex<-inndex+1	
    print(inndex)
}








