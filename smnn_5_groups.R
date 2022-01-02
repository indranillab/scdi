library(batchelor)
xsmnn<-mnnCorrect(t(data1),t(data2))

dat<-assays(xsmnn)$corrected

library(tsne)
source("fast_tsne.R")
xdat<-fftRtsne(t(dat))

pdf("5_groups_smnn.pdf")
plot(xdat,col=cols,pch=(1+groups))
legend(-0.15,0.15,legend=c("Batch1 Group1","Batch1 Group2","Batch1 Group3","Batch1 Group4","Batch1 Group5","Batch2 Group1","Batch2 Group2","Batch2 Group3","Batch2 Group4","Batch2 Group5"),col=c(1,2,3,4,5,1,2,3,4,5),pch=(c(1,1,1,1,1,2,2,2,2,2)))
dev.off()

















