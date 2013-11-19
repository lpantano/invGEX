prepare.CEUTSI.geovadis<-function(info.a,obj_exp){
	source("parameters.R")
	CEU<-read.table(pathCEU,header=T,sep="\t",row.names=1)
	TSI<-read.table(pathTSI,header=T,sep="\t",row.names=1)

	info<-info.a[info.a[,4]=="CEU" | info.a[,4]=="TSI",]

	inv.names<-as.character(info[,2])
	names(inv.names)<-as.character(info[,1])
	inv.pop<-as.character(info[,4])
	names(inv.pop)<-as.character(info[,1])
	inv.sex<-as.character(info[,3])
	names(inv.sex)<-as.character(info[,1])


	all.t<-cbind(CEU,TSI[,2:ncol(TSI)])
	all<-all.t[grep("ENSG",row.names(all.t)),]

	all<-all[,!is.na(inv.names[names(all)])]
	conditions<-as.vector(inv.names[names(all)])
	sex<-as.vector(inv.sex[names(all)])
	pop<-as.vector(inv.pop[names(all)])
	design<-data.frame(condition=conditions,sex=sex,pop=pop)
	design$sex[design$sex=="M"]<-"Male"
	design$sex[design$sex=="F"]<-"Female"
	design$sex<-factor(design$sex,levels=c("Female","Male"))
	row.names(design)<-names(all)
	table.f<-all
	inv.sex<-design$sex
	names(inv.sex)<-row.names(design)
	obj_exp[[1]]<-design
	obj_exp[[2]]<-table.f
	return(obj_exp)
}