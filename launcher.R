
# parallel computing
dode<-function(x){
	source("general.pipeline.function.R")
	do.analysis(x)
}

# I had an object with all names
# load(all.inversions.obj)

library(parallel)
cl<-makeCluster(4) 
clusterEvalQ(cl)
# all.inversions is a character vector with the names of inversions
# I recommend test first with only one inversion to check the pipeline works. 
t<-system.time(clusterApplyLB(cl,all.inversions, dode))
write.table(t,"/home/lpantano/time.process")
stopCluster(cl)
