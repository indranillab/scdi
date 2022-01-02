sink("output_scdi.txt")
library(tictoc)
tic("Time taken:")
Rprof(tf <- "rprof_scdi.log",memory.profiling=TRUE)
source("time_memory_seurat.R")
Rprof(NULL)
summaryRprof(tf,memory="both")
toc()
sink()



