library(splatter)

set.seed(0)
x1<-splatSimulate(params=newSplatParams(),method="groups",batchCells=c(100,100),group.prob=c(0.5,0.5),nGenes=10000,batch.facLoc=1,batch.facScale=1)
x2<-splatSimulate(params=newSplatParams(),method="groups",batchCells=c(100,100),group.prob=c(0.5,0.5),nGenes=10000)


data1<-assays(x1)$counts[,(colData(x1)$Batch=="Batch1")]
data2<-cbind(assays(x1)$counts[,101:200][,(colData(x1)$Group[101:200]=="Group1")],assays(x2)$counts[,101:200][,(colData(x2)$Group[101:200]=="Group2")])



data1<-t(data1)
data2<-t(data2)

data<-rbind(data1,data2)



for(i in 1:(dim(data1)[1]))
{
	data1[i,]<-data1[i,]/sqrt(sum((data1[i,])^2))
}


ta1<-table(colData(x1)$Group[101:200])
ta2<-table(colData(x2)$Group[101:200])
Gr<-c(as.vector(unlist(colData(x1)$Group))[1:100],rep("Group1",ta1[1]),rep("Group2",ta2[2]))

cols<-rep(NA,length(Gr))
cols[Gr=="Group1"]<-1
cols[Gr=="Group2"]<-2


groups1<-c(rep("Batch1",(dim(data1)[1])),rep("Batch2",(dim(data2)[1])))

groups<-rep(NA,length(groups1))

groups[groups1=="Batch1"]<-0
groups[groups1=="Batch2"]<-1

groups1<-c(colData(x1)$Group[1:100],rep("Group1",sum(colData(x1)$Group=="Group1")),rep("Group2",sum(colData(x2)$Group=="Group2")))
cols<-rep(NA,length(groups1))
cols[groups1=="Group1"]<-"red"
cols[groups1=="Group2"]<-"black"
cols[groups1=="Group3"]<-"blue"
cols[groups1=="Group4"]<-"green"
cols[groups1=="Group5"]<-"magenta"


clust<-rep(NA,length(groups1))
clust[groups1=="Group1"]<-0
clust[groups1=="Group2"]<-1
clust[groups1=="Group3"]<-2
clust[groups1=="Group4"]<-3
#clust[groups1=="Group5"]<-4
clust1<-clust[1:500]
clust2<-clust[501:1000]


data<-rbind(data1,data2)


source("scdi12.R")
K1<-2
K2<-2
Z2<-scdi_integrate(data,groups,K1,K2)



cdi<-fftRtsne(Z2$Z2,rand_seed=0)


cols<-rep(NA,length(Gr))
cols[Gr=="Group1"]<-1
cols[Gr=="Group2"]<-2

pdf("scdi_2.pdf")
plot(cdi,col=cols,pch=(1+groups),xlab="Dimension I",ylab="Dimension II",main="SCDI")
dev.off()



library(batchelor)
xsmnn<-mnnCorrect(t(data1),t(data2))


dat<-assays(xsmnn)$corrected

library(tsne)
source("fast_tsne.R")
xdat<-fftRtsne(t(dat))

pdf("smnn_2.pdf")
plot(xdat,col=cols,pch=(1+groups),xlab="Dimension I",ylab="Dimension II",main="SCDI")
dev.off()


library(Seurat)

dat1<-CreateSeuratObject(counts = t(data1), min.cells = 3)
dat2<-CreateSeuratObject(counts = t(data2), min.cells = 3)
lst=list(dat1,dat2)
names(lst)<-c("dat1","dat2")
intanch<-FindIntegrationAnchors(object.list=lst,k.filter=30)

dat<-IntegrateData(intanch)

xdata<-GetAssayData(dat)

source("fast_tsne.R")
xxx<-fftRtsne(t(as.matrix(xdata)),rand_seed=0)

pdf("seurat_2.pdf")
plot(xxx,col=cols,pch=(1+groups),xlab="Dimension I",ylab="Dimension II",main="SCDI")
dev.off()





library(splatter)

set.seed(0)
x1<-splatSimulate(params=newSplatParams(),method="groups",batchCells=c(100,100),group.prob=c(0.33,0.33,0.34),nGenes=10000,batch.facLoc=0.5,batch.facScale=0.5)
x2<-splatSimulate(params=newSplatParams(),method="groups",batchCells=c(100,100),group.prob=c(0.33,0.33,0.34),nGenes=10000)


data1<-assays(x1)$counts[,(colData(x1)$Batch=="Batch1")]
data2<-cbind(assays(x1)$counts[,101:200][,(colData(x1)$Group[101:200]=="Group1")],assays(x2)$counts[,101:200][,(colData(x2)$Group[101:200]=="Group2")],assays(x2)$counts[,101:200][,(colData(x2)$Group[101:200]=="Group3")])



data1<-t(data1)
data2<-t(data2)

data<-rbind(data1,data2)



for(i in 1:(dim(data1)[1]))
{
	data1[i,]<-data1[i,]/sqrt(sum((data1[i,])^2))
}



for(i in 1:(dim(data2)[1]))
{
	data2[i,]<-data2[i,]/sqrt(sum((data2[i,])^2))
}



ta1<-table(colData(x1)$Group[101:200])
ta2<-table(colData(x2)$Group[101:200])
Gr<-c(as.vector(unlist(colData(x1)$Group))[1:100],rep("Group1",ta1[1]),rep("Group2",ta2[2]),rep("Group3",ta2[3]))

cols<-rep(NA,length(Gr))
cols[Gr=="Group1"]<-1
cols[Gr=="Group2"]<-2
cols[Gr=="Group3"]<-3


groups1<-c(rep("Batch1",(dim(data1)[1])),rep("Batch2",(dim(data2)[1])))

groups<-rep(NA,length(groups1))

groups[groups1=="Batch1"]<-0
groups[groups1=="Batch2"]<-1

groups1<-c(colData(x1)$Group[1:100],rep("Group1",sum(colData(x1)$Group=="Group1")),rep("Group2",sum(colData(x2)$Group=="Group2")))
cols<-rep(NA,length(groups1))
cols[groups1=="Group1"]<-"red"
cols[groups1=="Group2"]<-"black"
cols[groups1=="Group3"]<-"blue"
cols[groups1=="Group4"]<-"green"
cols[groups1=="Group5"]<-"magenta"


clust<-rep(NA,length(groups1))
clust[groups1=="Group1"]<-0
clust[groups1=="Group2"]<-1
clust[groups1=="Group3"]<-2
clust[groups1=="Group4"]<-3
#clust[groups1=="Group5"]<-4
clust1<-clust[1:500]
clust2<-clust[501:1000]


data<-rbind(data1,data2)


source("scdi12.R")
K1<-3
K2<-3
Z2<-scdi_integrate(data,groups,K1,K2)



ddi<-fftRtsne(Z2$Z2,rand_seed=0)


cols<-rep(NA,length(Gr))
cols[Gr=="Group1"]<-1
cols[Gr=="Group2"]<-2
cols[Gr=="Group3"]<-3

pdf("scdi_3.pdf")
plot(ddi,col=cols,pch=(1+groups),xlab="Dimension I",ylab="Dimension II",main="SCDI")
dev.off()



library(batchelor)
xsmnn<-mnnCorrect(t(data1),t(data2))


dat<-assays(xsmnn)$corrected

library(tsne)
source("fast_tsne.R")
xdat<-fftRtsne(t(dat))

pdf("smnn_3.pdf")
plot(xdat,col=cols,pch=(1+groups),xlab="Dimension I",ylab="Dimension II",main="SCDI")
dev.off()


library(Seurat)

dat1<-CreateSeuratObject(counts = t(data1), min.cells = 3)
dat2<-CreateSeuratObject(counts = t(data2), min.cells = 3)
lst=list(dat1,dat2)
names(lst)<-c("dat1","dat2")
intanch<-FindIntegrationAnchors(object.list=lst,k.filter=30)

dat<-IntegrateData(intanch)

xdata<-GetAssayData(dat)

source("fast_tsne.R")
xxx<-fftRtsne(t(as.matrix(xdata)),rand_seed=0)

pdf("seurat_2.pdf")
plot(xxx,col=cols,pch=(1+groups),xlab="Dimension I",ylab="Dimension II",main="SCDI")
dev.off()



