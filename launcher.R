 library (snowfall)

# initialize cluster
 sfInit (parallel=TRUE , cpus=3)


# parallel computing
dode<-function(x){
	source("general.pipeline.function.R")
	do.analysis(x)
}
load(all.inversions.obj)
system.time(sfLapply(all.inversions, dode))

# stop cluster
sfStop()
