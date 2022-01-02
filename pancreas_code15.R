a1<-read.csv("GSM2230757_human1_umifm_counts.csv")
a2<-read.csv("GSM2230758_human2_umifm_counts.csv")

info1<-a1[,1:3]
info2<-a2[,1:3]

a1<-a1[,-(1:3)]
a2<-a2[,-(1:3)]

set.seed(12345)

dat1<-as.matrix(a1)
dat2<-as.matrix(a2)

for(i in 1:(dim(dat1)[1]))
{
	dat1[i,]<-(dat1[i,]-mean(dat1[i,]))/sd(dat1[i,])
}

for(i in 1:(dim(dat2)[1]))
{
	dat2[i,]<-(dat2[i,]-mean(dat2[i,]))/sd(dat2[i,])
}

source("fast_tsne.R")
ddd1<-fftRtsne(dat1)
ddd2<-fftRtsne(dat2)

K1<-9
K2<-8
source("scdi19.R")
clust1<-graph_cluster(ddd1,K1)
clust2<-graph_cluster(ddd2,K2)

n1<-0.01*(dim(dat1)[1])
n2<-0.01*(dim(dat2)[1])

ta1<-table(clust1)
ta2<-table(clust2)
num1<-names(ta1)[ta1>n1]
num2<-names(ta2)[ta2>n2]
dta1<-dat1[(clust1%in%num1),]
dta2<-dat2[(clust2%in%num2),]

clus1<-clust1[(clust1%in%num1)]
clus2<-clust2[(clust2%in%num2)]


save.image("pan1.Rdata")

Z2<-scdi_integrate(dta1,dta2,clus1,clus2,ddd1,ddd2)


set1<-Z2$set1
set2<-Z2$set2

a3<-dta1[set1,]
a4<-dta2[set2,]

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

groups<-c(rep(0,(dim(data1)[1])),rep(1,(dim(data2)[1])))


save.image("pamcreas1.Rdata")


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
X3<-Z2$Z2
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
	corvec<-cor(a5,t(dataa1))
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
	corvec<-cor(a5,t(dataa3))
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
	corvec<-cor(a5,t(dataa1))
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
	corvec<-cor(a5,t(dataa3))
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
		vec<-as.numeric(cormat1[k1,ind2])
		vec<-(1+vec)/2
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
		vec<-as.numeric(cormat2[k1,(ind2)])
		vec<-(1+vec)/2
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
    dis<--log(pmax(corvec1,0))
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
X3<-Z2$Z2

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
	corvec<-cor(a5,t(dataa1))
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
	corvec<-cor(a5,t(dataa3))
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
	corvec<-cor(a5,t(dataa1))
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
	corvec<-cor(a5,t(dataa3))
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
		vec<-(1+vec)/2
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
		vec<-as.numeric(cormat2[k1,(ind2)])
		vec<-(1+vec)/2
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
    dis<--log(pmax(corvec1,0))
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


rm(list=ls())
load("pancreas3.Rdata")
rm(list=setdiff(ls(),c("set1","set2")))


load("pan1.Rdata")

xyz1<-which(!((clust1%in%num1)))
xyz2<-which(!((clust2%in%num2)))
wxyz1<-which(((clust1%in%num1)))
wxyz2<-which(((clust2%in%num2)))
wxyz1<-setdiff(wxyz1,set1)
wxyz2<-setdiff(wxyz2,set2)

rm(list=setdiff(ls(),c("xyz1","xyz2","wxyz1","wxyz2","set1","set2")))

load("pancreas4.Rdata")

xyz3<-wxyz1[!match1]
xyz4<-wxyz2[!match2]

inn1<-c(xyz1,xyz3)
inn2<-c(xyz2,xyz4)



l<-Z2$l
alpha<-Z2$alpha
sigma<-Z2$sigma
sigma1<-Z2$sigma1
Z3<-Z2$Z2
groups<-Z2$groups
cluster<-Z2$cluster
dataa<-Z2$dataa
X3<-Z2$Z2

library(emulator)

Z2Z2<-array(0,dim<-c(length(inn1),dim(X3)[2]))

inndex<-1

dataa1<-dataa
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

cormat1<-array(0,dim<-c((dim(dataa1)[1]),length(c(inn1))))

for(i in inn1)
{
	a5<-a1[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa1,a5)
	corvec<-cor(a5,t(dataa1))
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




Z4Z4<-array(0,dim<-c(length(inn2),dim(X3)[2]))

inndex<-1

dataa3<-dataa
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

cormat2<-array(0,dim<-c((dim(dataa3)[1]),length(inn2)))
for(i in inn2)
{
	a5<-a2[i,]
	a5<-as.numeric(a5)
	a5<-a5/sqrt(sum(a5^2))
	dataanew<-rbind(dataa3,a5)
	corvec<-cor(a5,t(dataa3))
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




density_simulate<-function(ds,N)
{
	weights<-(ds$y[-1]+ds$y[-(length(ds$y))])/2
	samp<-sample(1:length(weights),N,prob=weights,replace=TRUE)
	return((runif(N)*ds$bw+ds$x[samp]))
	
}

set.seed(12345)

clus1<-Z2$cluster[(Z2$groups==0)]
clus2<-Z2$cluster[(Z2$groups==1)]

corrmat<-array(0,dim<-c(length(c(clus1,clus2)),length(c(inn1,inn2))))

corrmat[,(1:(length(inn1)))]<-(1+cormat1)/2
corrmat[,((length(inn1)+1):(length(inn1)+length(inn2)))]<-(1+cormat2)/2

t1<-as.numeric(names(table(clus)))
for(i in 1:(length(table(clus))))
{
	ind1<-which(clus==t1[i])
	ind2<-which(clustnew11==t1[i])
	ind3<-which(clustnew12==t1[i])
	for(k1 in ind1)
	{
		vec<-as.numeric(cormat1[k1,ind2])
		if(length(vec)<2)
		{
			next
		}
		ds<-density(vec,from=0,kernel="rectangular")
		rand<-density_simulate(ds,length(ind3))
		corrmat[k1,(length(inn1)+ind3[order(cormat2[k1,ind3])])]<-sort(rand)
		rand1<-density_simulate(ds,length(ind2))
		corrmat[k1,ind2[order(cormat1[k1,ind2])]]<-sort(rand1)
	}	
}


clus<-c(clus1,clus2)
Zmat2<-array(0,dim<-c((length(inn1)+length(inn2)),(dim(X3)[2])))
inndex<-1
clustn2<-c(clustnew11,clustnew12)
for(j in 1:(dim(corrmat)[2]))
{
    corvec<-corrmat[,j]	
	corvec1<-corvec[clus==clustn2[inndex]]
	M<-min(corvec1[corvec1>0])/1000
    dis<--log(pmax(corvec1,0))
    Z5<-X3[(clus==clustn2[inndex]),]
    if(sum(dim(Z5))>0)
    {
    y<-colMeans(Z5)
    }
    else
    {
    	y<-Z5
    }
    opt<-optim(y,cost,d=dis,X1=Z5,method="L-BFGS-B")
    Zmat2[inndex,]<-opt$par
    inndex<-inndex+1	
    print(inndex)
}


save.image("pancreas6.Rdata")








