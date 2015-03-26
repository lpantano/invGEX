
# parallel computing
dode<-function(x){
	source("general.pipeline.function.R")
	do.analysis(x)
}
load("all.inversions.obj")

library(parallel)
cl<-makeCluster(4) 
clusterEvalQ(cl)
t<-system.time(clusterApplyLB(cl,all.inversions, dode))
write.table(t,"/home/lpantano/time.process")
stopCluster(cl)
