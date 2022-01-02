data<-iris
data1<-data[,1:4]

data2<-as.matrix(data1)

cl<-kmeans(data2,3)

graph<-array(0,dim<-c(150,150))
for(k in 1:100)
{
	indd<-sample(150,30)
	data3<-data2[indd,]
	ddd<-as.matrix(dist(data3))
	ddd<-ddd+diag(2*max(ddd),(dim(ddd)[1]))
	for(i in 1:(dim(ddd)[1]))
	{
		for(j in 1:(dim(ddd)[2]))
		{
			if((ddd[i,j]==min(ddd[i,]))&(ddd[i,j]==min(ddd[,j])))
			{
				graph[indd[i],indd[j]]<-graph[indd[i],indd[j]]+1
				graph[indd[j],indd[i]]<-graph[indd[j],indd[i]]+1
			}
		}
	}
	print(k)
	
}




