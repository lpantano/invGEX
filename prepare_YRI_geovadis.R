

prepare.YRI.geovadis<-function(info.a,obj_exp){
	source("parameters.R")
	table<-read.table(paste(sep="",pathYRI,"genes/counts.",type,".tab"),header=T,sep="\t",row.names=1)
	#table<-table[,1:35]
	#table<-table[1:(nrow(table)-5),]
	#info<-read.table("/home/shareddata/Bioinformatics/iRNASeq/geovadis/YRI/HsInv102_YRI.v2.csv",sep="\t",header=T)
	info<-info.a[info.a[,4]=="YRI",]

	inv.names<-(info[,2])
	names(inv.names)<-(info[,1])
	inv.sex<-info[,3]
	names(inv.sex)<-info[,1]

	table.f<-table[,!is.na(inv.names[names(table)])]
	conditions<-as.vector(inv.names[names(table.f)])
	sex<-as.vector(inv.sex[names(table.f)])
	design<-data.frame(condition=conditions,sex=sex)
	row.names(design)<-names(table.f)
	obj_exp[[1]]<-design
	obj_exp[[2]]<-table.f
	return(obj_exp)
	#continue get.table.filter.R
}