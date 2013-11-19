
get.table.filter<-function(obj_exp,sex=F){

	library(DESeq2)

	design<-obj_exp[[1]]
	table.f<-obj_exp[[2]]

	################
	all.f<-table.f
	all.f<-all.f[grep("ENSG",row.names(all.f)),]
	design2<-design

	if (sex==T){
		
		design2<-design2[design2$condition=="STD" | design2$condition=="INV",]
		ind<-design2[,1]
		names(ind)<-row.names(design2)
		all.f<-all.f[,!is.na(ind[names(all.f)])]
		design2$condition<-factor(design2$condition,levels=unique(design2$condition))
	}
	#####################################

	std<-row.names(design2[design2$condition=="STD",])
	het<-row.names(design2[design2$condition=="HET",])
	inv<-row.names(design2[design2$condition=="INV",])
	
	go.w.deseq=0
	go.w.deseq.sex=0
	#print(paste("doing:",pop.selected,GEN1,"vs",GEN2))

	num.gen<-length(unique(design2$condition))

	if (num.gen>=2){
		
		go.w.deseq=1
		
		getmin<-function(x){length(which(x<=0))}
		which.res<-apply(all.f,1,getmin)
		
		which.res<-which.res[which.res<=round(nrow(design2)*0.50)]

		genes.f<-names(which.res)
		all.f<-all.f[genes.f,]

		#####################################
		dds <- DESeqDataSetFromMatrix(countData = all.f,
			colData = design2,
			design = ~ condition)
	    dds <- estimateSizeFactors( dds )
		design2<-design2[sizeFactors(dds)<=1.9 & sizeFactors(dds)>=0.5, ]

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

		if (length(het)<4){
			design2<-design2[design2$condition!="HET",]
		}
		if (length(std)<4){
			design2<-design2[design2$condition!="STD",]
		}
		if (length(inv)<4){
			design2<-design2[design2$condition!="INV",]
		}

		num.gen<-length(unique(design2$condition))

		if (num.gen<2){
			go.w.deseq=0
		}else{
			design2$condition<-as.character(design2$condition)
			design2$condition[design2$condition=="STD"]<-"cSTD"
			design2$condition<-factor(design2$condition,levels=unique(design2$condition))
		}

		ind<-design2[,1]
		names(ind)<-row.names(design2)
		table.fil<-all.f[,!is.na(ind[names(all.f)])]
		obj_exp[[1]]<-design2
		obj_exp[[2]]<-table.fil
	}
	
	obj_exp[[3]]<-go.w.deseq
	obj_exp[[4]]<-go.w.deseq.sex
	

	return(obj_exp)
}
	