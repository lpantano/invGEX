source("parameters.R")
library(DESeq2)
library(ggplot2)
library(qvalue)

deseq.random<-function(obj_exp,path,pi0){

	list.qvalues<-vector()
	list.num.sign<-vector()
	dir.create(paste(sep="",path,"/random"))
	design2<-obj_exp[[1]]
	table.fil<-obj_exp[[2]]
	go.w.deseq.sex<-obj_exp[[4]]
	pvalues<-vector()
	pi0.l<-vector("list")
	for (iter in 1:50){
		#print(iter)
		if ((iter/10)%%1==0){print(paste("num pi0.rand below real pi0",min(pvalues),"at",iter))}
		if ( iter>1){
			if (min(pvalues)>=0.15){
				print(paste("stop because not possibility to get lower pval:",min(pvalues)))
				write.table(list.qvalues,paste(sep="",path,"/random.stop"),quote=F,row.names=F)
				#iter<-55
				break
			}
		}
		dse.log<-0
		while(class(dse.log)!="DESeqDataSet"){
		 design.r<-design2
		 design.r$condition<-sample(design2$condition,nrow(design2))
		 write.table(design.r,paste(sep="",path,"/random/",iter,"gen.random"),quote=F,row.names=T)
		 if (go.w.deseq.sex==0){
			dse <- DESeqDataSetFromMatrix(countData = table.fil, colData = design.r,design = ~ sex + condition)
		 }else{
			dse <- DESeqDataSetFromMatrix(countData = table.fil, colData = design.r,design = ~ condition)
		 }
		 dse.log<-try(dse <- DESeq(dse,fitType="parametric",quiet=TRUE),silent=TRUE)
		}

		comp<-resultsNames(dse)
		for (nameres in comp[grep("condition",comp)]){
			#print(nameres)
			res <- results(dse,nameres,independentFiltering=FALSE,cooksCutoff=FALSE)
			qobj <- qvalue(res$pvalue[!is.na(res$pvalue)],fdr.level=0.1)
			#t<-try(qsummary(qobj),silent=TRUE)
			pi0.l[[nameres]]<-c(pi0.l[[nameres]],qobj$pi0)
			pvalues[nameres]<-sum(pi0[nameres]>pi0.l[[nameres]])/51
		}

		if (length(levels(design2$condition))>2){
			res<-results(dse, contrast=c("condition","INV","HET"),
				independentFiltering=FALSE,cooksCutoff=FALSE)
			qobj <- qvalue(res$pvalue[!is.na(res$pvalue)],fdr.level=0.1)
			#t<-try(qsummary(qobj),silent=TRUE)
			nameres<-"condition_INV_vs_HET"
			pi0.l[[nameres]]<-c(pi0.l[[nameres]],qobj$pi0)
			pvalues[nameres]<-sum(pi0[nameres]>pi0.l[[nameres]])/51
		}
		#print(pvalues)

	}
	return(pvalues)
}