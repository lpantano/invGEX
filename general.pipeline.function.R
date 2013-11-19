"%+%" <- function(x,y) paste(x,y,sep="") 
do.analysis<-function(invname){

	source("parameters.R")
	inversions<-read.table(iMLPA,sep="\t",header=T)
	inversions2<-read.table(dMLPA,sep="\t",header=T)
	inversions<-merge(inversions,inversions2[,c(1,4:ncol(inversions2))],by=1,all=T)

	pos.inv<-read.table(file.pos,header=T)
	pos.inv.aut<-pos.inv
	inversion.run=FALSE
	if (parameter.pipeline=="new"){
		if (!(file.exists(paste(sep="",path_res,invname,"/",type,"/doing")))){
			inversion.run<-TRUE
		}
	}else if (parameter.pipeline=="unfinished"){
		if (!(file.exists(paste(sep="",path_res,invname,"/",type,"/done")))){
			inversion.run<-TRUE
		}
	}else if (parameter.pipeline=="repeat"){
			inversion.run<-TRUE
	}
	if (inversion.run){
		print(invname)
		dir.create(paste(sep="",path_res,invname),showWarnings=F)
		ind<-which(invname==as.character(names(inversions)))
		#info.a<-inversions[,c(1,ind,4,3)]
		info.a<-inversions[,c(1,ind,3,2)]
		names(info.a)[2]<-"condition"

		#dir.create(paste(sep="",path_res,invname),showWarnings = FALSE)
		dir.create(paste(sep="",path_res,invname,"/",type),showWarnings = FALSE)
		write.table("doing",paste(sep="",path_res,invname,"/",type,"/doing"))


		source("info.tr.R")
		info.a<-info.tr(info.a)

		##print(table(info.a[,c(2,3,4)]))

		if (pos.inv.aut[invname,1]=="chrX"){
			print("doing ChrX")
			source("pipeline.deseq.chrX.R")
			pipeline.deseq.chrX(info.a,type,invname)
		}else{
			print("doing autosomic")
			source("pipeline.all.v2.R")
			pipeline.deseq(info.a,invname)
		}
		source("summary.all.R")
		#do.summary(invname)
		write.table("done",paste(sep="",path_res,invname,"/",type,"/done"))
		
	}

	return(invname)
}