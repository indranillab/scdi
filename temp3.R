

set3<-Z2$set3
set4<-Z2$set4

l<-Z2$l
alpha<-Z2$alpha
sigma<-Z2$sigma
sigma1<-Z2$sigma1
Z3<-Z2$Z2
groups<-Z2$groups
cluster<-Z2$cluster
dataa<-Z2$dataa
X3<-Z2$X2
dataa<-Z2$dataa

library(GPLVM)
XX1<-t(chol(data[(groups==0),]%*%t(data[(groups==0),])))
XX2<-t(chol(data[(groups==1),]%*%t(data[(groups==1),])))
XX3<-forwardsolve(XX1,diag(length(diag(XX1))))
XX4<-forwardsolve(XX2,diag(length(diag(XX2))))

source("gplvm_SE.R")
#K<-gplvm.SE(X3[(groups==0),],l,alpha,sigma)

#Kinv<-solve(K)
#Kvecmat<-Kinv%*%data[(groups==0),]

library(emulator)

Z2Z2<-array(0,dim<-c(length(set3),dim(X3)[2]))

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
clist<-unique(cl1)

cormat1<-array(0,dim<-c((dim(dataa1)[1]),length(set3)))
for(i in set3)
{
	a5<-dta1[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa1,a5)
	corvec<-cov(a5,t(dataa1))*nrow(data1)
	disclust<-NULL
	for(j in 1:(length(unique(cl1))))
	{
		disclust[j]<-mean(corvec[cl1==clist[j]])
	}
	kk<-which(disclust==max(disclust))
	corvec1<-corvec[cl1==clist[kk]]
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z4<-X3[(groups==0),]
    Z5<-Z4[(cl1==clist[kk]),]
    if(sum(dim(Z5))>0)
    {
    y<-colMeans(Z5)
    }
    else
    {
    	y<-Z5
    }
    #opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    #Z2Z2[inndex,]<-opt$par
    clustnew11[inndex]<-clist[kk]
    cormat1[,inndex]<-corvec
    inndex<-inndex+1	
    print(inndex)
}




Z4Z4<-array(0,dim<-c(length(set4),dim(X3)[2]))

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
clist<-unique(cl1)

cormat2<-array(0,dim<-c((dim(dataa3)[1]),length(set4)))
for(i in set4)
{
	a5<-dta2[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa3,a5)
	corvec<-cov(a5,t(dataa3))*nrow(data2)
	disclust<-NULL
	for(j in 1:(length(unique(cl1))))
	{
		disclust[j]<-mean(corvec[cl1==clist[j]])
	}
	kk<-which(disclust==max(disclust))
	corvec1<-corvec[cl1==clist[kk]]
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z4<-X3[(groups==1),]
    Z5<-Z4[(cl1==clist[kk]),]
    if(sum(dim(Z5))>0)
    {
    y<-colMeans(Z5)
    }
    else
    {
    	y<-Z5
    }
    #opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    #Z4Z4[inndex,]<-opt$par
    clustnew12[inndex]<-clist[kk]
    cormat2[,inndex]<-corvec
    inndex<-inndex+1	
    print(inndex)
}



Z6Z6<-array(0,dim<-c(length(set4),dim(X3)[2]))

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
clist<-unique(cl1)

cormat3<-array(0,dim<-c((dim(dataa1)[1]),length(set4)))
for(i in set4)
{
	a5<-dta2[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa1,a5)
	corvec<-cov(a5,t(dataa1))*nrow(data1)
	disclust<-NULL
	for(j in 1:(length(unique(cl1))))
	{
		disclust[j]<-mean(corvec[cl1==clist[j]])
	}
	kk<-which(disclust==max(disclust))
	corvec1<-corvec[cl1==clist[kk]]
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z4<-X3[(groups==0),]
    Z5<-Z4[(cl1==clist[kk]),]
    if(sum(dim(Z5))>0)
    {
    y<-colMeans(Z5)
    }
    else
    {
    	y<-Z5
    }
    #opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    #Z6Z6[inndex,]<-opt$par
    clustnew21[inndex]<-clist[kk]
    cormat3[,inndex]<-corvec
    inndex<-inndex+1	
    print(inndex)
}




Z8Z8<-array(0,dim<-c(length(set3),dim(X3)[2]))

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
clist<-unique(cl1)

cormat4<-array(0,dim<-c((dim(dataa3)[1]),length(set3)))
for(i in set3)
{
	a5<-dta1[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa3,a5)
	corvec<-cov(a5,t(dataa3))*nrow(data2)
	disclust<-NULL
	for(j in 1:(length(unique(cl1))))
	{
		disclust[j]<-mean(corvec[cl1==clist[j]])
	}
	kk<-which(disclust==max(disclust))
	corvec1<-corvec[cl1==clist[kk]]
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z4<-X3[(groups==1),]
    Z5<-Z4[(cl1==clist[kk]),]
    if(sum(dim(Z5))>0)
    {
    y<-colMeans(Z5)
    }
    else
    {
    	y<-Z5
    }
    #opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    #Z8Z8[inndex,]<-opt$par
    clustnew22[inndex]<-clist[kk]
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

t1<-as.numeric(names(table(clus1)))
for(i in 1:(length(table(clus1))))
{
	ind1<-which(clus1==t1[i])
	ind2<-which(clustnew11==t1[i])
	ind3<-which(clustnew21==t1[i])
	for(k1 in ind1)
	{
		vec<-as.numeric(cormat1[k1,ind1])
		vec<-vec[vec>0]
		ds<-density(vec,from=0,kernel="rectangular")
		rand<-density_simulate(ds,length(ind3))
		corrmat[k1,(length(set3)+ind3[order(cormat3[k1,ind3])])]<-sort(rand)
		rand1<-density_simulate(ds,length(ind2))
		corrmat[k1,ind2[order(cormat1[k1,ind2])]]<-sort(rand1)
	}	
}

t2<-as.numeric(names(table(clus2)))
for(i in 1:(length(table(clus2))))
{
	ind1<-which(clus2==t2[i])
	ind2<-which(clustnew12==t2[i])
	ind3<-which(clustnew22==t2[i])
	for(k1 in ind1)
	{
		vec<-as.numeric(cormat3[k1,(ind2)])
		vec<-vec[vec>0]
		ds<-density(vec,from=0,kernel="rectangular")
		rand<-density_simulate(ds,length(ind3))
		corrmat[(length(clus1)+k1),(ind3[order(cormat4[k1,ind3])])]<-sort(rand)
		rand1<-density_simulate(ds,length(ind2))
		corrmat[(length(clus1)+k1),(length(set3)+ind2[order(cormat2[(k1),ind2])])]<-sort(rand1)
	}
}

match1<-(clustnew11==clustnew22)
match2<-(clustnew21==clustnew12)
clus<-c(clus1,clus2)
Zmat<-array(0,dim<-c((sum(match1)+sum(match2)),(dim(X3)[2])))
inndex<-1
clustn<-c(clustnew11[match1],clustnew12[match2])
cormat<-corrmat[,c(match1,match2)]
for(j in 1:(dim(cormat)[2]))
{
    corvec<-cormat[,j]	
	corvec1<-corvec[clus==clustn[inndex]]
	M<-min(corvec1[corvec1>0])/1000
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z5<-X3[(clus==clustn[inndex]),]
    if(sum(dim(Z5))>0)
    {
    y<-colMeans(Z5)
    }
    else
    {
    	y<-Z5
    }
    opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    Zmat[inndex,]<-opt$par
    inndex<-inndex+1	
    print(inndex)
}


save.image("pancreas4.Rdata")






set1<-Z2$set1
set2<-Z2$set2

l<-Z2$l
alpha<-Z2$alpha
sigma<-Z2$sigma
sigma1<-Z2$sigma1
Z3<-Z2$Z2
groups<-Z2$groups
cluster<-Z2$cluster
dataa<-Z2$dataa
X3<-Z2$X2

library(GPLVM)
XX1<-t(chol(data[(groups==0),]%*%t(data[(groups==0),])))
XX2<-t(chol(data[(groups==1),]%*%t(data[(groups==1),])))
XX3<-forwardsolve(XX1,diag(length(diag(XX1))))
XX4<-forwardsolve(XX2,diag(length(diag(XX2))))

source("gplvm_SE.R")
K<-gplvm.SE(X3[(groups==0),],l,alpha,sigma)

Kinv<-solve(K)
Kvecmat<-Kinv%*%data[(groups==0),]

library(emulator)

Z2Z2<-array(0,dim<-c(length(set1),dim(X3)[2]))

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
clist<-unique(cl1)

cormat1<-array(0,dim<-c((dim(dataa1)[1]),length(set1)))
for(i in set1)
{
	a5<-dta1[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa1,a5)
	corvec<-cov(a5,t(dataa1))*nrow(data1)
	disclust<-NULL
	for(j in 1:(length(unique(cl1))))
	{
		disclust[j]<-mean(corvec[cl1==clist[j]])
	}
	kk<-which(disclust==max(disclust))
	corvec1<-corvec[cl1==clist[kk]]
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z4<-X3[(groups==0),]
    Z5<-Z4[(cl1==clist[kk]),]
    if(sum(dim(Z5))>0)
    {
    y<-colMeans(Z5)
    }
    else
    {
    	y<-Z5
    }
    #opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    #Z2Z2[inndex,]<-opt$par
    clustnew11[inndex]<-clist[kk]
    cormat1[,inndex]<-corvec
    inndex<-inndex+1	
    print(inndex)
}




Z4Z4<-array(0,dim<-c(length(set2),dim(X3)[2]))

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
clist<-unique(cl1)

cormat2<-array(0,dim<-c((dim(dataa3)[1]),length(set2)))
for(i in set2)
{
	a5<-dta2[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa3,a5)
	corvec<-cov(a5,t(dataa3))*nrow(data2)
	disclust<-NULL
	for(j in 1:(length(unique(cl1))))
	{
		disclust[j]<-mean(corvec[cl1==clist[j]])
	}
	kk<-which(disclust==max(disclust))
	corvec1<-corvec[cl1==clist[kk]]
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z4<-X3[(groups==1),]
    Z5<-Z4[(cl1==clist[kk]),]
    if(sum(dim(Z5))>0)
    {
    y<-colMeans(Z5)
    }
    else
    {
    	y<-Z5
    }
    #opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    #Z4Z4[inndex,]<-opt$par
    clustnew12[inndex]<-clist[kk]
    cormat2[,inndex]<-corvec
    inndex<-inndex+1	
    print(inndex)
}



Z6Z6<-array(0,dim<-c(length(set2),dim(X3)[2]))

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
clist<-unique(cl1)

cormat3<-array(0,dim<-c((dim(dataa1)[1]),length(set2)))
for(i in set2)
{
	a5<-dta2[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa1,a5)
	corvec<-cov(a5,t(dataa1))*nrow(data1)
	disclust<-NULL
	for(j in 1:(length(unique(cl1))))
	{
		disclust[j]<-mean(corvec[cl1==clist[j]])
	}
	kk<-which(disclust==max(disclust))
	corvec1<-corvec[cl1==clist[kk]]
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z4<-X3[(groups==0),]
    Z5<-Z4[(cl1==clist[kk]),]
    if(sum(dim(Z5))>0)
    {
    y<-colMeans(Z5)
    }
    else
    {
    	y<-Z5
    }
    #opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    #Z6Z6[inndex,]<-opt$par
    clustnew21[inndex]<-clist[kk]
    cormat3[,inndex]<-corvec
    inndex<-inndex+1	
    print(inndex)
}




Z8Z8<-array(0,dim<-c(length(set1),dim(X3)[2]))

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
clist<-unique(cl1)

cormat4<-array(0,dim<-c((dim(dataa3)[1]),length(set1)))
for(i in set1)
{
	a5<-dta1[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa3,a5)
	corvec<-cov(a5,t(dataa3))*nrow(data2)
	disclust<-NULL
	for(j in 1:(length(unique(cl1))))
	{
		disclust[j]<-mean(corvec[cl1==clist[j]])
	}
	kk<-which(disclust==max(disclust))
	corvec1<-corvec[cl1==clist[kk]]
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z4<-X3[(groups==1),]
    Z5<-Z4[(cl1==clist[kk]),]
    if(sum(dim(Z5))>0)
    {
    y<-colMeans(Z5)
    }
    else
    {
    	y<-Z5
    }
    #opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    #Z8Z8[inndex,]<-opt$par
    clustnew22[inndex]<-clist[kk]
    cormat4[,inndex]<-corvec
    inndex<-inndex+1	
    print(inndex)
}











density_simulate<-function(ds,N)
{
	weights<-(ds$y[-1]+ds$y[-(length(ds$y))])/2
	samp<-sample(1:length(weights),N,prob=weights,replace=TRUE)
	return((runif(N)*ds$bw+ds$x[samp]))
	
}

set.seed(12345)

clus1<-Z2$cluster[(Z2$groups==0)]
clus2<-Z2$cluster[(Z2$groups==1)]

corrmat<-array(0,dim<-c(length(c(clus1,clus2)),length(c(set1,set2))))

corrmat[(1:(length(clus1))),(1:(length(set1)))]<-cormat1
corrmat[((length(clus1)+1):(length(clus1)+length(clus2))),((length(set1)+1):(length(set1)+length(set2)))]<-cormat2

t1<-as.numeric(names(table(clus1)))
for(i in 1:(length(table(clus1))))
{
	ind1<-which(clus1==t1[i])
	ind2<-which(clustnew11==t1[i])
	ind3<-which(clustnew21==t1[i])
	for(k1 in ind1)
	{
		vec<-as.numeric(cormat1[k1,ind1])
		vec<-vec[vec>0]
		ds<-density(vec,from=0,kernel="rectangular")
		rand<-density_simulate(ds,length(ind3))
		corrmat[k1,(length(set1)+ind3[order(cormat3[k1,ind3])])]<-sort(rand)
		rand1<-density_simulate(ds,length(ind2))
		corrmat[k1,ind2[order(cormat1[k1,ind2])]]<-sort(rand1)
	}	
}

t2<-as.numeric(names(table(clus2)))
for(i in 1:(length(table(clus2))))
{
	ind1<-which(clus2==t2[i])
	ind2<-which(clustnew12==t2[i])
	ind3<-which(clustnew22==t2[i])
	for(k1 in ind1)
	{
		vec<-as.numeric(cormat3[k1,(ind2)])
		vec<-vec[vec>0]
		ds<-density(vec,from=0,kernel="rectangular")
		rand<-density_simulate(ds,length(ind3))
		corrmat[(length(clus1)+k1),(ind3[order(cormat4[k1,ind3])])]<-sort(rand)
		rand1<-density_simulate(ds,length(ind2))
		corrmat[(length(clus1)+k1),(length(set1)+ind2[order(cormat2[(k1),ind2])])]<-sort(rand1)
	}
}

match1<-(clustnew11==clustnew22)
match2<-(clustnew21==clustnew12)
clus<-c(clus1,clus2)
Zmat1<-array(0,dim<-c((sum(match1)+sum(match2)),(dim(X3)[2])))
inndex<-1
clustn1<-c(clustnew11[match1],clustnew12[match2])
cormat<-corrmat[,c(match1,match2)]
for(j in 1:(dim(cormat)[2]))
{
    corvec<-cormat[,j]	
	corvec1<-corvec[clus==clustn1[inndex]]
	M<-min(corvec1[corvec1>0])/1000
    dis<--log(pmax(corvec1,0)/alpha^2)*(2*l^2)
    Z5<-X3[(clus==clustn1[inndex]),]
    if(sum(dim(Z5))>0)
    {
    y<-colMeans(Z5)
    }
    else
    {
    	y<-Z5
    }
    opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    Zmat1[inndex,]<-opt$par
    inndex<-inndex+1	
    print(inndex)
}


save.image("pancreas5.Rdata")



