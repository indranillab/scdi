library(splatter)

set.seed(12345)
x<-splatSimulate(method="groups",batchCells=c(100,100),group.prob=c(0.5,0.5))

data<-t(assays(x)$counts)

for(i in 1:(dim(data)[1]))
{
	data[i,]<-data[i,]/sum(data[i,])*10^6
}

data1<-data[1:100,]
data2<-data[101:200,]

for(j in 1:(dim(data1)[2]))
{
	if(sd(data1[,j])>0)
	{
	data1[,j]<-data1[,j]/sd(data1[,j])
	}else
	{
		data1[,j]<-rep(0,length(data1[,j]))
	}
}
for(j in 1:(dim(data2)[2]))
{
	if(sd(data2[,j])>0)
	{
	data2[,j]<-data2[,j]/sd(data2[,j])
	}else
	{
		data2[,j]<-rep(0,length(data2[,j]))
	}
}

data<-rbind(data1,data2)










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
clust1<-clust[1:100]
clust2<-clust[101:200]



