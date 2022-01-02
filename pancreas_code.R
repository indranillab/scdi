a1<-read.csv("GSM2230757_human1_umifm_counts.csv")
a2<-read.csv("GSM2230758_human2_umifm_counts.csv")

info1<-a1[,1:3]
info2<-a2[,1:3]

a1<-a1[,-(1:3)]
a2<-a2[,-(1:3)]

set.seed(12345)

a3<-a1[(sample(1937,200)),]
a4<-a2[(sample(1724,200)),]

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

groups<-c(rep(0,200),rep(1,200))




library(scater)
data7<-data1
data8<-data2

K1<-10
K2<-10
    sce1 <- SingleCellExperiment(assays = list(counts=t(data7)))
    logcounts(sce1)<-log(1+t(data7))
    rowData(sce1)<-paste0("gene",(1:(dim(data7)[2])))
    rowData(sce1)$feature_symbol<-paste0("gene",(1:(dim(data7)[2])))
	
    sce2 <- SingleCellExperiment(assays = list(counts=t(data8)))
    logcounts(sce2)<-log(1+t(data8))
    rowData(sce2)<-paste0("gene",(1:(dim(data8)[2])))
    rowData(sce2)$feature_symbol<-paste0("gene",(1:(dim(data8)[2])))
        
    xxx1<-sc3(sce1,ks=4:7)
    xxx2<-sc3(sce2,ks=4:7)


    cluster1<-as.numeric(as.vector(unlist(colData(xxx1))))
    cluster2<-as.numeric(as.vector(unlist(colData(xxx2))))

#data7<-data7[,(colSums(data7)>0)]
#data8<-data8[,(colSums(data8)>0)]
#library(TSCAN)
#cl1<-exprmclust(t(data7),clusternum=2:9)


library(tsne)
source("fast_tsne.R")

x1<-fftRtsne(data1)

x2<-fftRtsne(data2)




source("scdi8.R")
K1<-5
K2<-5
Z2<-scdi_integrate(data,groups,K1,K2)









