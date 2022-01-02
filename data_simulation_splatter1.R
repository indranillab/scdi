library(splatter)

set.seed(0)
x<-splatSimulate(params=newSplatParams(),method="groups",batchCells=c(100,100),group.prob=c(0.2,0.2,0.2,0.2,0.2))

data<-t(assays(x)$counts)


for(i in 1:(dim(data)[1]))
{
	data[i,]<-data[i,]/sum((data[i,]))*1000000
}

dataoriginal<-data


data1<-data[1:100,]
data2<-data[101:200,]

for(j in 1:(dim(data1)[2]))
{
	if(sd(data1[,j])>0)
	{
	#data1[,j]<-data1[,j]/sd(data1[,j])
	}else
	{
	#	data1[,j]<-rep(0,length(data1[,j]))
	}
}
for(j in 1:(dim(data2)[2]))
{
	if(sd(data2[,j])>0)
	{
	#data2[,j]<-data2[,j]/sd(data2[,j])
	}else
	{
	#	data2[,j]<-rep(0,length(data2[,j]))
	}
}


data<-rbind(data1,data2)



for(i in 1:(dim(data)[1]))
{
	data[i,]<-data[i,]/sqrt(sum((data[i,])^2))
}

data1<-data[1:100,]
data2<-data[101:200,]








groups1<-colData(x)$Batch

groups<-rep(NA,length(groups1))

groups[groups1=="Batch1"]<-0
groups[groups1=="Batch2"]<-1

groups1<-colData(x)$Group
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


source("scdi10.R")
K1<-5
K2<-5
Z2<-scdi_integrate(data,groups,K1,K2)


xx<-fftRtsne(Z2$Z2,rand_seed=0)
pdf("5_groups_scdi.pdf")
plot(xx,col=cols,pch=(1+groups),main="SCDI",xlab="Dimension I",ylab="Dimension II",ylim=c(-30,30),xlim=c(-20,40))
legend(25,30,legend=c("Batch1 Group1","Batch1 Group2","Batch1 Group3","Batch1 Group4","Batch1 Group5","Batch2 Group1","Batch2 Group2","Batch2 Group3","Batch2 Group4","Batch2 Group5"),col=c(1,2,3,4,5,1,2,3,4,5),pch=(c(1,1,1,1,1,2,2,2,2,2)))
dev.off()

library(splatter)

set.seed(0)
x<-splatSimulate(params=newSplatParams(),method="groups",batchCells=c(100,100),group.prob=c(0.25,0.25,0.25,0.25))

data<-t(assays(x)$counts)


for(i in 1:(dim(data)[1]))
{
	data[i,]<-data[i,]/sum((data[i,]))*1000000
}

dataoriginal<-data


data1<-data[1:100,]
data2<-data[101:200,]

for(j in 1:(dim(data1)[2]))
{
	if(sd(data1[,j])>0)
	{
	#data1[,j]<-data1[,j]/sd(data1[,j])
	}else
	{
	#	data1[,j]<-rep(0,length(data1[,j]))
	}
}
for(j in 1:(dim(data2)[2]))
{
	if(sd(data2[,j])>0)
	{
	#data2[,j]<-data2[,j]/sd(data2[,j])
	}else
	{
	#	data2[,j]<-rep(0,length(data2[,j]))
	}
}


data<-rbind(data1,data2)



for(i in 1:(dim(data)[1]))
{
	data[i,]<-data[i,]/sqrt(sum((data[i,])^2))
}

data1<-data[1:100,]
data2<-data[101:200,]








groups1<-colData(x)$Batch

groups<-rep(NA,length(groups1))

groups[groups1=="Batch1"]<-0
groups[groups1=="Batch2"]<-1

groups1<-colData(x)$Group
cols<-rep(NA,length(groups1))
cols[groups1=="Group1"]<-"red"
cols[groups1=="Group2"]<-"black"
cols[groups1=="Group3"]<-"blue"
cols[groups1=="Group4"]<-"green"


clust<-rep(NA,length(groups1))
clust[groups1=="Group1"]<-0
clust[groups1=="Group2"]<-1
clust[groups1=="Group3"]<-2
clust[groups1=="Group4"]<-3
#clust[groups1=="Group5"]<-4
clust1<-clust[1:500]
clust2<-clust[501:1000]


source("scdi10.R")
K1<-4
K2<-4
Z2<-scdi_integrate(data,groups,K1,K2)



xx<-fftRtsne(Z2$Z2,rand_seed=0)
pdf("4_groups_scdi.pdf")
plot(xx,col=cols,pch=(1+groups),main="SCDI",xlim=c(-30,50),xlab="Dimension I",ylab="Dimension II")
legend(30,22,legend=c("Batch1 Group1","Batch1 Group2","Batch1 Group3","Batch1 Group4","Batch2 Group1","Batch2 Group2","Batch2 Group3","Batch2 Group3","Batch2 Group4"),col=c(1,2,3,4,1,2,3,4),pch=(c(1,1,1,1,2,2,2,2)))
dev.off()




library(splatter)

set.seed(0)
x<-splatSimulate(params=newSplatParams(),method="groups",batchCells=c(100,100),group.prob=c(0.33,0.33,0.34))

data<-t(assays(x)$counts)


for(i in 1:(dim(data)[1]))
{
	data[i,]<-data[i,]/sum((data[i,]))*1000000
}

dataoriginal<-data


data1<-data[1:100,]
data2<-data[101:200,]

for(j in 1:(dim(data1)[2]))
{
	if(sd(data1[,j])>0)
	{
	#data1[,j]<-data1[,j]/sd(data1[,j])
	}else
	{
	#	data1[,j]<-rep(0,length(data1[,j]))
	}
}
for(j in 1:(dim(data2)[2]))
{
	if(sd(data2[,j])>0)
	{
	#data2[,j]<-data2[,j]/sd(data2[,j])
	}else
	{
	#	data2[,j]<-rep(0,length(data2[,j]))
	}
}


data<-rbind(data1,data2)



for(i in 1:(dim(data)[1]))
{
	data[i,]<-data[i,]/sqrt(sum((data[i,])^2))
}

data1<-data[1:100,]
data2<-data[101:200,]

groups1<-colData(x)$Batch

groups<-rep(NA,length(groups1))

groups[groups1=="Batch1"]<-0
groups[groups1=="Batch2"]<-1

groups1<-colData(x)$Group
cols<-rep(NA,length(groups1))
cols[groups1=="Group1"]<-"red"
cols[groups1=="Group2"]<-"black"
cols[groups1=="Group3"]<-"blue"


clust<-rep(NA,length(groups1))
clust[groups1=="Group1"]<-0
clust[groups1=="Group2"]<-1
clust[groups1=="Group3"]<-2
clust[groups1=="Group4"]<-3
#clust[groups1=="Group5"]<-4
clust1<-clust[1:500]
clust2<-clust[501:1000]


source("scdi10.R")
K1<-3
K2<-3
Z2<-scdi_integrate(data,groups,K1,K2)



xx<-fftRtsne(Z2$Z2,rand_seed=0)
pdf("3_groups_scdi.pdf")
plot(xx,col=cols,pch=(1+groups),main="SCDI",xlab="Dimension 1",ylab="Dimension 2",xlim=c(-20,50))
legend(30,15,legend=c("Batch1 Group1","Batch1 Group2","Batch1 Group3","Batch2 Group1","Batch2 Group2","Batch2 Group3"),col=c(1,2,3,1,2,3),pch=(c(1,1,1,2,2,2)))
dev.off()










library(splatter)

set.seed(0)
x<-splatSimulate(params=newSplatParams(),method="groups",batchCells=c(100,100),group.prob=c(0.5,0.5))

data<-t(assays(x)$counts)


for(i in 1:(dim(data)[1]))
{
	data[i,]<-data[i,]/sum((data[i,]))*1000000
}

dataoriginal<-data


data1<-data[1:100,]
data2<-data[101:200,]

for(j in 1:(dim(data1)[2]))
{
	if(sd(data1[,j])>0)
	{
	#data1[,j]<-data1[,j]/sd(data1[,j])
	}else
	{
	#	data1[,j]<-rep(0,length(data1[,j]))
	}
}
for(j in 1:(dim(data2)[2]))
{
	if(sd(data2[,j])>0)
	{
	#data2[,j]<-data2[,j]/sd(data2[,j])
	}else
	{
	#	data2[,j]<-rep(0,length(data2[,j]))
	}
}


data<-rbind(data1,data2)



for(i in 1:(dim(data)[1]))
{
	data[i,]<-data[i,]/sqrt(sum((data[i,])^2))
}

data1<-data[1:100,]
data2<-data[101:200,]

groups1<-colData(x)$Batch

groups<-rep(NA,length(groups1))

groups[groups1=="Batch1"]<-0
groups[groups1=="Batch2"]<-1

groups1<-colData(x)$Group
cols<-rep(NA,length(groups1))
cols[groups1=="Group1"]<-"red"
cols[groups1=="Group2"]<-"black"


clust<-rep(NA,length(groups1))
clust[groups1=="Group1"]<-0
clust[groups1=="Group2"]<-1
clust[groups1=="Group3"]<-2
clust[groups1=="Group4"]<-3
#clust[groups1=="Group5"]<-4
clust1<-clust[1:500]
clust2<-clust[501:1000]


source("scdi10.R")
K1<-2
K2<-2
Z2<-scdi_integrate(data,groups,K1,K2)



xx<-fftRtsne(Z2$Z2,rand_seed=0)
pdf("2_groups_scdi.pdf")
plot(xx,col=cols,pch=(1+groups),main="SCDI",xlim=c(-30,50),xlab="Dimension 1",ylab="Dimension 2")
legend(30,10,legend=c("Batch1 Group1","Batch1 Group2","Batch2 Group1","Batch2 Group2"),col=c(1,2,1,2),pch=(c(1,1,2,2)))
dev.off()





pdf("2_groups_scdi.pdf")
plot(Z2,col=(1+clust),pch=(1+groups),main="SCDI")
legend(-0.15,0.15,legend=c("Batch1 Group1","Batch1 Group2","Batch2 Group1","Batch2 Group2"),col=c(1,2,1,2),pch=(c(1,1,2,2)))
dev.off()







data<-dataoriginal
data1<-data[1:100,]
data2<-data[101:200,]
source("smnn_5_groups.R")


save.image("smnn_3_groups.Rdata")

source("seurat_workflow.R")



