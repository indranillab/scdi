N<-50
cols<-c(rep("red",N),rep("black",N),rep("blue",N),rep("green",N),rep("purple",N),rep("cyan",N))

group<-c(rep(0,150),rep(1,150))



source("gplvm_polar.R")
source("tsne_internal.R")
xxx<-tsne(data)


ydata1<-xxx$ydata
alpha<-xxx$alpha
beta<-xxx$beta

#sigma<-var(as.numeric(ydata))

XX<-array(0,dim<-c(dim(data)[1],length(unique(group))))
counter<-1
for(le in (unique(group)))
{
	XX[,counter]<-as.numeric(group==le)
	counter<-counter+1
}

ydata2<-ydata1-colMeans(ydata1)
sigma<-NULL

data1<-data-(rep(1,(dim(data)[1])))%*%t(colMeans(data))
sigma1<-NULL
for(j in 1:(dim(data1)[2]))
{
fit<-lm(data1[,j]~-1+XX)
sigma[j]<-sum((as.numeric(fit$residuals))^2)/(dim(XX)[1]-dim(XX)[2])
sigma1[j]<-(sum((as.numeric(data1[,j]))^2)-sum((as.numeric(fit$residuals))^2)-sigma[j])/150
}


tau<-sigma1/sigma


source("gplvm_new.R")
ydata3<-tsne(data,ydata2,sigma,sigma1,group)




source("gplvm3.R")
source("tsne_internal.R")
xxx<-tsne(data)














