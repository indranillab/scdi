library(scRNAseq)
sce1 <- ZeiselBrainData()
sce2 <- TasicBrainData()
set.seed(19999)
sce1 <- sce1[,sample(ncol(sce1), 200, replace=TRUE)]
sce2 <- sce2[,sample(ncol(sce2), 200, replace=TRUE)]
universe <- intersect(rownames(sce1), rownames(sce2))
sce1 <- sce1[universe,]
sce2 <- sce2[universe,]

data1<-t(assays(sce1)$counts)
data2<-t(assays(sce2)$counts)

for(i in 1:(dim(data1)[1]))
{
	data1[i,]<-data1[i,]/sum(data1[i,])*1000000
}

for(i in 1:(dim(data2)[1]))
{
	data2[i,]<-data2[i,]/sum(data2[i,])*1000000
}

colnames(data1)<-rownames(sce1)
colnames(data2)<-rownames(sce2)

library(tsne)
source("fast_tsne.R")

x1<-fftRtsne(data1)
x2<-fftRtsne(data2)

rho11<-NULL
rho12<-NULL
rho21<-NULL
rho22<-NULL

for(j in 1:(dim(data1)[2]))
{
	rho11[j]<-cor(data1[,j],x1[,1],method="spearman")
}

for(j in 1:(dim(data1)[2]))
{
	rho12[j]<-cor(data1[,j],x1[,2],method="spearman")
}

for(j in 1:(dim(data2)[2]))
{
	rho21[j]<-cor(data2[,j],x1[,1],method="spearman")
}

for(j in 1:(dim(data1)[2]))
{
	rho22[j]<-cor(data2[,j],x1[,2],method="spearman")
}


data31<-data1[,(match(1:100,order(na.omit(abs(rho11)))))]
data32<-data1[,(match(1:100,order(na.omit(abs(rho12)))))]
data41<-data2[,(match(1:100,order(na.omit(abs(rho21)))))]
data42<-data2[,(match(1:100,order(na.omit(abs(rho22)))))]

genenames<-Reduce(union,list(colnames(data31),colnames(data32),colnames(data41),colnames(data42)))

data3<-data1[,match(genenames,colnames(data1))]
data4<-data2[,match(genenames,colnames(data2))]

data1<-data3
data2<-data4

data<-rbind(data1,data2)

groups<-c(rep(0,200),rep(1,200))



