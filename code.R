entropy<-function(X,cols,groups)
{
	en<-0
	for(i in 1:(max(cols)))
	{
		X1<-X[((cols==i)&(groups==0)),]
		X2<-X[((cols==i)&(groups==1)),]
        V1<-cov(X1)
        V2<-cov(X2)
        V<-cov(X[(cols==i),])
	    en<-en+(dim(X1)[1])*log(det(V1))+(dim(X2)[1])*log(det(V2))-(dim(X1)[1]+dim(X2)[1])*log(det(V))
	}
	return(en)
}



ent<-array(0,dim<-c(4,3))
colnames(ent)<-c("SCDD","SMNN","Seurat")
for(i in 2:5)
{
	rm(list=setdiff(ls(),c("i","entropy","ent")))
	filename<-paste0("simulation_",i,".Rdata")
	load(filename)
	if(i==2)
	{
	ent[(i-1),1]<-entropy(cdi,cols,groups)
	}else
	{
	ent[(i-1),1]<-entropy(ddi,cols,groups)
	}
	ent[(i-1),2]<-entropy(xdat,cols,groups)
	ent[(i-1),3]<-entropy(xxx,cols,groups)
	file<-paste0("entropy_",i,".pdf")
	library(ggplot2)
	dtm<-data.frame()
	pdf(file)
	barplot(-ent[(i-1),],col="red",ylim=c(0,600))
	abline(0,0)
	dev.off()
	
	
}




ent<-array(0,dim<-c(4,3))
colnames(ent)<-c("SCDD","SMNN","Seurat")
for(i in 2:5)
{
	rm(list=setdiff(ls(),c("i","entropy","ent")))
	filename<-paste0("simulation_",i,"1.Rdata")
	load(filename)
	if(i==2)
	{
	ent[(i-1),1]<-entropy(cdi,cols,groups)
	}else
	{
	ent[(i-1),1]<-entropy(ddi,cols,groups)
	}
	ent[(i-1),2]<-entropy(xdat,cols,groups)
	ent[(i-1),3]<-entropy(xxx,cols,groups)
	file<-paste0("entropy_",i,"1.pdf")
	library(ggplot2)
	dtm<-data.frame()
	pdf(file)
	barplot(-ent[(i-1),],col="red",ylim=c(0,1300))
	abline(0,0)
	dev.off()
	
	
}















