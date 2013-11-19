pipeline.deseq<-function(info.a,invname){
	source("parameters.R")
	###GEOVADIS YRI INV/EXON
	pop.selected<-"YRI"

	###GEOVADIS YRI HET/GEN
	source("prepare_YRI_geovadis.R")
	obj_exp<-vector("list")
	obj_exp.raw<-prepare.YRI.geovadis(info.a,obj_exp)
	
	GEN1="STD"
	GEN2="HET"
	source("get.table.filter.v2.R")
	obj_exp<-get.table.filter(obj_exp.raw,GEN1,GEN2)

	outfolder<-paste(sep="",pop.selected,".DESeq2",type,".het.parametric")
	source("deseq2.v2.R")
	deseq2(obj_exp,outfolder,invname,pop.selected)

	# ###GEOVADIS YRI INV/GEN

	GEN1="STD"
	GEN2="INV"
	obj_exp<-get.table.filter(obj_exp.raw,GEN1,GEN2)
	
	outfolder<-paste(sep="",pop.selected,".DESeq2",type,".inv.parametric")
	deseq2(obj_exp,outfolder,invname,pop.selected)

	# ###GEOVADIS YRI INV/GEN

	GEN1="HET"
	GEN2="INV"
	obj_exp<-get.table.filter(obj_exp.raw,GEN1,GEN2)
	
	outfolder<-paste(sep="",pop.selected,".DESeq2",type,".inv.het.parametric")
	deseq2(obj_exp,outfolder,invname,pop.selected)

	#########################
	#########################
	#CEUTSI
	#########################
	###GEOVADIS CEUTSI HET/EXON
	pop.selected<-"CEUTSI"

	###GEOVADIS CEUTSI HET/GENES
	source("prepare_CEUTSI.R")
	obj_exp<-vector("list")
	obj_exp.raw<-prepare.CEU.geovadis(info.a,obj_exp)

	GEN1="STD"
	GEN2="HET"
	obj_exp<-get.table.filter(obj_exp.raw,GEN1,GEN2)

	outfolder<-paste(sep="",pop.selected,".DESeq2",type,".het.parametric")
	deseq2(obj_exp,outfolder,invname,pop.selected)

	###GEOVADIS CEUTSI INV/GENES

	GEN1="STD"
	GEN2="INV"
	obj_exp<-get.table.filter(obj_exp.raw,GEN1,GEN2)

	outfolder<-paste(sep="",pop.selected,".DESeq2",type,".inv.parametric")
	deseq2(obj_exp,outfolder,invname,pop.selected)

	###GEOVADIS CEUTSI INV/GENES

	GEN1="HET"
	GEN2="INV"#GEN is the studied feature
	obj_exp<-get.table.filter(obj_exp.raw,GEN1,GEN2)

	outfolder<-paste(sep="",pop.selected,".DESeq2",type,".inv.het.parametric")
	deseq2(obj_exp,outfolder,invname,pop.selected)

	#########################
}