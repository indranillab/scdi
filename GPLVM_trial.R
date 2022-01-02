X<-data
q<-2
iterations=1000
plot.freq=100
classes=NULL
Z.init="PCA"
num.init.params=100
Z.normal.prior<-c(0,1)


source("BCGPLVM.R")
source("BCSGPLVM.R")
source("clustersExample.R")
source("GPLVM10.R")
source("kernels.R")
source("LSA-BCSGPLVM.R")
source("missingData.R")
source("noisyComparison.R")
source("structuredGPLVM.R")
source("swissRollExample.R")
source('syntheticData.R')


groups<-c(rep(0,150),rep(1,150))


