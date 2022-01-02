

load("G1_S_bio.Rdata")

data<-t(cbind(data11[,-1],data21[,-1]))
datanew<-array(0,dim<-dim(data))
for(j in 1:(dim(datanew)[2]))
{
	datanew[,j]<-as.numeric(data[,j])
}
data<-datanew
data1<-data[1:96,]
data2<-data[97:192,]


set.seed(0)
library(tsne)
source("fast_tsne.R")
	x1<-fftRtsne(data1)
	x2<-fftRtsne(data2)
	
x<-fftRtsne(as.matrix(data))















