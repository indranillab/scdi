batches<-1+groups
corrected.distance.matrix <- QuantNorm(data,batches,method='row/column', cor_method='pearson', logdat=F,standardize = T,tol=1e-4)

lst<-list(t(data[1:150,]),t(data[151:300,]))
xxx<-BUSgibbs(lst,3)



datanew<-data
for(i in 151:300)
{
	datanew[i,]<-data[i,]-xxx$gamma[,2]
}



ei<-svd(datanew)

