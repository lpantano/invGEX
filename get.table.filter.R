
get.table.filter<-function(obj_exp,GEN1,GEN2){

	library(DESeq2)

	design<-obj_exp[[1]]
	table.f<-obj_exp[[2]]

	################
	all.f<-table.f
	all.f<-all.f[grep("ENSG",row.names(all.f)),]
	design2<-design

	#####################################

	std<-row.names(design2[design2$condition==GEN1,])
	het<-row.names(design2[design2$condition==GEN2,])
	go.w.deseq=0
	go.w.deseq.sex=0
	if (length(std)>4 & length(het)>4 ){
		#print(paste("doing:",pop.selected,GEN1,"vs",GEN2))
		go.w.deseq=1
		getmin<-function(x){length(which(x<=1))}
		which.res.std<-apply(all.f[,std],1,getmin)
		which.res.het<-apply(all.f[,het],1,getmin)
		which.res.std<-which.res.std[which.res.std<=length(std)*0.27]
		which.res.het<-which.res.het[which.res.het<=length(het)*0.27]
		genes.f<-unique(c(names(which.res.std),names(which.res.het)))
		all.f<-all.f[genes.f,]
		all.f<-round(all.f,digits=0)

		#####################################
		
		dds <- DESeqDataSetFromMatrix(countData = all.f,
			colData = design2,
			design = ~ condition)
	    dds <- estimateSizeFactors( dds )
		design2<-design2[sizeFactors(dds)<=1.9 & sizeFactors(dds)>=0.5, ]
		design2<-design2[design2$condition==GEN1 | design2$condition==GEN2,]
		design2$condition<-factor(design2$condition,levels=c(GEN1,GEN2))

		if (min(table(design2$condition[design2$sex=="Female"]))<=1){
			#do males
#			print("doing only males")
			design2<-design2[design2$sex=="Male",]
			design2$sex<-factor(design2$sex,levels="Male")
			go.w.deseq.sex<-1
		}else if (min(table(design2$condition[design2$sex=="Male"]))<=1){
			#do female
#			print("doing only females")
			design2<-design2[design2$sex=="Female",]
			design2$sex<-factor(design2$sex,levels="Female")
			go.w.deseq.sex<-1
		}

		ind<-design2[,1]
		names(ind)<-row.names(design2)
		table.fil<-all.f[,!is.na(ind[names(all.f)])]

	}else{
		print("skip this comparison")

	}
	obj_exp[[3]]<-go.w.deseq
	obj_exp[[4]]<-go.w.deseq.sex
	obj_exp[[1]]<-design2
	obj_exp[[2]]<-table.fil

	return(obj_exp)
}
	