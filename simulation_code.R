#set.seed(0)
library(mvtnorm)
data1<-rmvnorm(50,rep(-3,30),diag(30))
data2<-rmvnorm(50,rep(3,30),diag(30))
data3<-rmvnorm(50,rep(6,30),diag(30))
data4<-rmvnorm(50,rep(-3,30),diag(30))
data5<-rmvnorm(50,rep(3,30),diag(30))
data6<-rmvnorm(50,rep(6,30),diag(30))
batch1<-rnorm(30,0,10)
batch2<-rnorm(30,0,10)
for(i in 1:(dim(data1)[1]))
{
	data1[i,]<-data1[i,]+batch1
}

for(i in 1:(dim(data2)[1]))
{
	data2[i,]<-data2[i,]+batch1
}

for(i in 1:(dim(data3)[1]))
{
	data3[i,]<-data3[i,]+batch1
}

for(i in 1:(dim(data4)[1]))
{
	data4[i,]<-data4[i,]+batch2
}

for(i in 1:(dim(data5)[1]))
{
	data5[i,]<-data5[i,]+batch2
}

for(i in 1:(dim(data6)[1]))
{
	data6[i,]<-data6[i,]+batch2
}




data<-rbind(data1,data2,data3,data4,data5,data6)



groups<-c(rep(0,150),rep(1,150))







for(i in 1:(dim(data)[1]))
{
	data[i,]<-data[i,]/sd(data[i,])
}





