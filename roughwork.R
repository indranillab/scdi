source("scdi8.R")


cl1<-graph_cluster(x1,7)


cl2<-graph_cluster(x2,7)

dists<-array(0,dim<-c(length(table(cl1)),length(table(cl2))))

library(pdist)
for(i in 1:length(table(cl1)))
{
	for(j in 1:length(table(cl2)))
	{
		dat1<-data1[(which(cl1==i)),]
		dat2<-data2[(which(cl2==j)),]
		dists[i,j]<-mean(as.matrix(pdist(dat1,dat2)))
		
		
	}
	print(i)
}



dat1old<-data1
data2old<-data2
cl2old<-cl2



data2<-data2old[-(which(cl2old==2)),]
cl2<-cl2old[-which(cl2old==2)]








