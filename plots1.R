for(i in 2:5)
{
	rm(list=setdiff(ls(),c("i","entropy","ent")))
	filename<-paste0("simulation_",i,"1.Rdata")
	load(filename)
	library(splatter)
	library(Rtsne)
	matt<-rbind(t(assays(x1)$counts),t(assays(x2)$counts))
	for(i1 in 1:(dim(matt)[1]))
	{
		matt[i1,]<-matt[i1,]/sum(matt[i1,])*1000000
	}
	ccddii<-Rtsne(data)
	pdf(paste0("data_",i,"1.pdf"))
    plot(ccddii$Y,col=cols,pch=(1+groups),xlab="Dimension I",ylab="Dimension II",main="Uncorrected")
    dev.off()
    if(i==2)
    {
    pdf(paste0("scdi_",i,"1.pdf"))
    plot(cdi,col=cols,pch=(1+groups),xlab="Dimension I",ylab="Dimension II",main="SCDI")
    dev.off()
    }else
    {
    pdf(paste0("scdi_",i,"1.pdf"))
    plot(ddi,col=cols,pch=(1+groups),xlab="Dimension I",ylab="Dimension II",main="SCDI")
    dev.off()
    }
    pdf(paste0("smnn_",i,"1.pdf"))
    plot(xdat,col=cols,pch=(1+groups),xlab="Dimension I",ylab="Dimension II",main="SMNN")
    dev.off()
    pdf(paste0("seurat_",i,"1.pdf"))
    plot(xxx,col=cols,pch=(1+groups),xlab="Dimension I",ylab="Dimension II",main="Seurat")
    dev.off()
	
	
}






