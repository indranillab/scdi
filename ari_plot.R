load("pancreas7.Rdata")


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
plot(ddd,col=cols1[cols],xlab="Dimension I",ylab="Dimension II",main="SCDI",pch=c(rep(1,length(match1)),rep(2,length(match2)),rep(1,(dim(info9)[1])),rep(2,dim(info10)[1])))
dev.off()

newclust<-graph_cluster(ddd,10)

library(fossil)
scdi_ari<-rand.index(c(clustn,clustn2),cols)



info9<-info1[inn1,]
info10<-info2[inn2,]

rm(list=setdiff(ls(),c("dta1","dta2","set3","set4","match1","match2","sset1","sset2","info4","info6","info9","info10","graph_cluster","info3","info5","xxx","xdat","oridata","ddd")))

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
plot(xxx,col=cols1[cols],xlab="Dimension I",ylab="Dimension II",main="Seurat",pch=c(rep(1,length(match1)),rep(1,(dim(info3)[1])),rep(2,length(match2)),rep(2,(dim(info5)[1]))))
dev.off()


newclust<-graph_cluster(xxx,8)

library(fossil)
seurat_ari<-rand.index(newclust,cols)




library(batchelor)
xsmnn<-mnnCorrect(t(ddat1),t(ddat2))


dat<-assays(xsmnn)$corrected


oridata<-fftRtsne(rbind(ddat1,ddat2))

cols<-c(info4[match1,3],info3[,3],info6[match2,3],info5[,3])

library(tsne)
source("fast_tsne.R")
xdat<-fftRtsne(t(dat),rand_seed=12345)

cols1<-c("red","blue","green","orange","gray","yellow","cyan","black","violet","magenta","purple","brown","orchid","pink")


pdf("pancreas_smnn.pdf")
cols1<-c("red","blue","green","orange","gray","yellow","cyan","black","violet","magenta","purple","brown","orchid","pink")
plot(xdat,col=cols1[cols],xlab="Dimension I",ylab="Dimension II",main="SMNN",pch=c(rep(1,length(match1)),rep(1,(dim(info3)[1])),rep(2,length(match2)),rep(2,(dim(info5)[1]))))
dev.off()


cols<-c(info4[match1,3],info3[,3],info6[match2,3],info5[,3])

source("fast_tsne.R")
oridata<-fftRtsne(rbind(ddat1,ddat2),rand_seed=12345)

pdf("pancreas_uncorrected.pdf")
cols1<-c("red","blue","green","orange","gray","yellow","cyan","black","violet","magenta","purple","brown","orchid","pink")
plot(oridata,col=cols1[cols],xlab="Dimension I",ylab="Dimension II",main="Uncorrected",pch=c(rep(1,length(match1)),rep(1,(dim(info3)[1])),rep(2,length(match2)),rep(2,(dim(info5)[1]))))
dev.off()


newclust<-graph_cluster(xdat,8)

library(fossil)
smnn_ari<-rand.index(newclust,cols)

cells<-c("Acinar","Stellate","Alpha","Beta","Delta","Ductal","Endfothelial","Epsilon","Gamma","Macrophage","Mast","Quiescent satelite","Schwann","T Cell")
cols1<-c("red","blue","green","orange","gray","yellow","cyan","black","violet","magenta","purple","brown","orchid","pink")

pdf("legend_pancreas.pdf")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =cells, pch=16, pt.cex=2, cex=1, bty='n',
    col = cols1,ncol=4)
mtext("Cell Types", at=0.5, cex=2)
dev.off()

pdf("ari_values.pdf")
barplot(c(scdi_ari,seurat_ari,smnn_ari),names.arg=c("SCDI","Seurat","SMNN"),ylim=c(0,1),xlab="Methods",ylab="Rand Index",main="Rand Index on Pancreas Data",col="magenta")
abline(0,0,lwd=3)
dev.off()


cells<-c("GSM2230757","GSM2230758")
cols1<-c("black","black")

pdf("legend2_pancreas.pdf")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1,pch=c(1,2))
legend("top", legend =cells, pch=c(1,2), pt.cex=1.5, cex=1, bty='n',
    col = cols1,ncol=1)
mtext("Datasets", at=0.5, cex=2)
dev.off()









