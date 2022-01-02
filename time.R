sink("output.txt")
library(tictoc)
tic("Time taken:")
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
#x<-1:10000
#y<-x^2
mm<-array(pi,dim<-c(15000,15000))
rm1<-mm*mm
rm2<-rm1+rm1
Rprof(NULL)
summaryRprof(tf,memory="both")
toc()
sink()



