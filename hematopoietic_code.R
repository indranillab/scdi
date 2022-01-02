a1<-read.table("GSE81682_HTSeq_counts.txt")
a3<-a1
a1<-a1[-1,-1]
colnames(a1)<-a3[1,-1]
rownames(a1)<-a3[-1,1]

a2<-read.table("GSE72857_umitab.txt")



