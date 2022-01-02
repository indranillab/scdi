rm(list=ls())

load("pancreas6.Rdata")

info3<-info1[set1,]
info4<-info1[set3,]
info5<-info2[set2,]
info6<-info2[set4,]


data1<-dta1
data2<-dta2
for(i in 1:(dim(data1)[1]))
{
	data1[i,]<-(dta1[i,]-mean(dta1[i,]))/sd(dta1[i,])
}
for(i in 1:(dim(data2)[1]))
{
	data2[i,]<-(dta2[i,]-mean(dta2[i,]))/sd(dta2[i,])
}

for(i in 1:(dim(data1)[1]))
{
	data1[i,]<-data1[i,]/sqrt(sum((data1[i,])^2))
}

for(i in 1:(dim(data2)[1]))
{
	data2[i,]<-data2[i,]/sqrt(sum((data2[i,])^2))
}

data<-rbind(data1,data2)



groups<-c(rep(0,(dim(data1)[1])),rep(1,dim(data2)[1]))


library(igraph)
library(maxmatching)	
library(GPLVM)


data[(groups==0),]<-sweep(data[(groups==0),],2,colMeans(data[(groups==0),]))
data[(groups==1),]<-sweep(data[(groups==1),],2,colMeans(data[(groups==1),]))
data3<-data
for(i in 1:(dim(data)[1]))
{
	data3[i,]<-data[i,]/sqrt(sum((data[i,])^2))
}


data1<-data3[(groups==0),]
data2<-data3[(groups==1),]

library(pdist)
#corm1<-pdist(t(dataa),t(data1))
#corm2<-pdist(t(dataa),t(data2))

#corm11<-corm1[(Z2$groups==0),]
#corm22<-corm2[(Z2$groups==1),]

sset1<-NULL
sset2<-NULL

check_equal<-function(x,y)
{
	N<-0
	for(i in 1:(length(x)))
	{
		if(x[i]!=y[i])
		{
			break
		}else
		{
			N<-N+1
		}
	}
	if(N==length(x))
	{
		return(TRUE)
	}else
	{
		return(FALSE)
	}
}


k<-1
for(i in 1:(dim(dataa)[1]))
{
	for(j in 1:(dim(data1)[1]))
	{
		if(check_equal(dataa[i,],data1[j,]))
		{
			sset1[k]<-j
			k<-k+1
		}
	}
    print(i)
}

k<-1
for(i in 1:(dim(dataa)[1]))
{
	for(j in 1:(dim(data2)[1]))
	{
		if(check_equal(dataa[i,],data2[j,]))
		{
			sset2[k]<-j
			k<-k+1
		}
	}
    print(i)
}


save.image("pancreas7.Rdata")

info11<-info1[(clust1%in%num1),]
info21<-info2[(clust2%in%num2),]

info3<-info11[sset1,]
info4<-info11[set3,]
info5<-info21[sset2,]
info6<-info21[set4,]

info7<-info4[match1,]
info8<-info6[match2,]

info9<-info1[inn1,]
info10<-info2[inn2,]

source("fast_tsne.R")
ddd<-fftRtsne(rbind(Zmat,Zmat2),rand_seed=12345)

cols<-c(info4[match1,3],info6[match2,3],info9[,3],info10[,3])

#dddd<-fftRtsne(Zmat)
#cols<-c(info4[match1,3],info6[match2,3])

pdf("pancreas_scdi.pdf")
cols1<-c("red","blue","green","orange","gray","yellow","cyan","black","violet","magenta","purple","brown","orchid","pink")
plot(ddd,col=cols1[cols],xlab="Dimension I",ylab="Dimension II",main="SCDI")
dev.off()

newclust<-graph_cluster(ddd,10)

library(pdfCluster)
scdi_ari<-adj.rand.index(newclust,cols)


ddat1<-rbind(dta1[set3[match1],],dta1[sset1,])
ddat2<-rbind(dta2[set4[match2],],dta2[sset2,])

colnames(ddat1)<-1:(dim(ddat1)[2])
colnames(ddat2)<-1:(dim(ddat2)[2])
rownames(ddat1)<-1:(dim(ddat1)[1])
rownames(ddat2)<-1:(dim(ddat2)[1])


library(Seurat)
ddat11<-CreateSeuratObject(counts = t(ddat1), min.cells = 3)
ddat22<-CreateSeuratObject(counts = t(ddat2), min.cells = 3)

lst=list(ddat11,ddat22)
names(lst)<-c("ddat11","ddat22")
intanch<-FindIntegrationAnchors(object.list=lst,k.filter=30)

dat<-IntegrateData(intanch)

xdata<-GetAssayData(dat)

source("fast_tsne.R")
xxx<-fftRtsne(t(as.matrix(xdata)),rand_seed=12345)

cols<-c(info4[match1,3],info3[,3],info6[match2,3],info5[,3])

pdf("pancreas_seurat.pdf")
cols1<-c("red","blue","green","orange","gray","yellow","cyan","black","violet","magenta","purple","brown","orchid","pink")
plot(xxx,col=cols1[cols],xlab="Dimension I",ylab="Dimension II",main="Seurat")
dev.off()


newclust<-graph_cluster(xxx,8)

library(pdfCluster)
scdi_ari<-adj.rand.index(newclust,cols)




library(batchelor)
xsmnn<-mnnCorrect(t(ddat1),t(ddat2))


dat<-assays(xsmnn)$corrected

save.image("pancreas8.Rdata")

oridata<-fftRtsne(rbind(ddat1,ddat2))

cols<-c(info4[match1,3],info3[,3],info6[match2,3],info5[,3])

rm(list=setdiff(ls(),c("cols","dat")))

library(tsne)
source("fast_tsne.R")
xdat<-fftRtsne(t(dat),rand_seed=12345)

cols1<-c("red","blue","green","orange","gray","yellow","cyan","black","violet","magenta","purple","brown","orchid","pink")


pdf("pancreas_smnn.pdf")
cols1<-c("red","blue","green","orange","gray","yellow","cyan","black","violet","magenta","purple","brown","orchid","pink")
plot(xdat,col=cols1[cols],xlab="Dimension I",ylab="Dimension II",main="SMNN")
dev.off()

rm(list=ls())
load("pancreas8.Rdata")
cols<-c(info4[match1,3],info3[,3],info6[match2,3],info5[,3])

rm(list=setdiff(ls(),c("cols","ddat1","ddat2")))

source("fast_tsne.R")
oridata<-fftRtsne(rbind(ddat1,ddat2),rand_seed=12345)

pdf("pancreas_uncorrected.pdf")
cols1<-c("red","blue","green","orange","gray","yellow","cyan","black","violet","magenta","purple","brown","orchid","pink")
plot(oridata,col=cols1[cols],xlab="Dimension I",ylab="Dimension II",main="Uncorrected")
dev.off()















