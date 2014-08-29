pipeline.deseq<-function(info.a,invname){
	source("parameters.R")
	###GEOVADIS YRI INV/EXON
	pop.selected<-"YRI"

	###GEOVADIS YRI HET/GEN
	source("prepare_YRI_geovadis.R")
	obj_exp<-vector("list")
	obj_exp.raw<-prepare.YRI.geovadis(info.a,obj_exp)
	
	source("get.table.filter.v2.R")
	obj_exp<-get.table.filter(obj_exp.raw)

	outfolder<-paste(sep="",pop.selected,".DESeq2",type,".parametric")
	source("deseq2.v2.R")
	deseq2(obj_exp,outfolder,invname,pop.selected)
	source("deseq2.additive.v2.R")
	deseq2.additive(obj_exp,outfolder,invname,pop.selected)

	#########################
	#########################
	#CEUTSI
	#########################
	###GEOVADIS CEUTSI HET/EXON
	pop.selected<-"CEUTSI"

	###GEOVADIS CEUTSI HET/GENES
	source("prepare_CEUTSI.R")
	obj_exp<-vector("list")
	obj_exp.raw<-prepare.CEUTSI.geovadis(info.a,obj_exp)

	obj_exp<-get.table.filter(obj_exp.raw)

	outfolder<-paste(sep="",pop.selected,".DESeq2",type,".parametric")
	source("deseq2.v2.R")
	deseq2(obj_exp,outfolder,invname,pop.selected)
	source("deseq2.additive.v2.R")
	deseq2.additive(obj_exp,outfolder,invname,pop.selected)


	#########################
}
