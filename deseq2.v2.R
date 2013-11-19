source("parameters.R")
library(DESeq2)
library(ggplot2)
library(qvalue)

send.random<-function(obj_exp,path,pi0){
	source("deseq.random.R")
	pvalues<-deseq.random(obj_exp,path,pi0)
	return(pvalues)
}

save.results<-function(res,path,name,pi0){
	vector<-c(sum(res$padj<=0.1,na.rm=T),sum(res$padj<=0.2,na.rm=T),sum(res$padj<=0.3,na.rm=T),sum(res$padj<=0.4,na.rm=T))
	
	fdr.info<-data.frame(fdr=c(0.1,0.2,0.3,0.4),count=vector)
	write.table(fdr.info,paste(sep="",path,"/",name,".fdr.slides"),quote=F,row.names=F)
	
	qobj <- qvalue(res$pvalue[!is.na(res$pvalue)],fdr.level=0.1)
	qwrite(qobj,filename=paste(sep="",path,"/",name,".qvalue.slides"))
	#t<-try(qsummary(qobj),silent=TRUE)
	pi0[name]<-qobj$pi0
	#num.sign.exp<-sum(qobj$significant)

	pdf(paste(sep="",path,"/",name,".qvalue.pdf"))
		qplot.log<-try(qplot(qobj),silent=TRUE)
	dev.off()

	# Overlaid histograms
	pdf(paste(sep="",path,"/",name,".pvalue.fdr.dist.pdf"))
		hist(res$pvalue,col=rgb(176/255,224/255,230/255,0.5),freq=FALSE,ylim=c(0,2),main="statistic distribution")
		hist(res$padj,col=rgb(102/255,153/255,204/255,0.5),add=T,freq=FALSE)
		legend("top",fill=c(rgb(176/255,224/255,230/255,0.5),rgb(102/255,153/255,204/255,0.5)),legend=c("pvalue","FDR"))
	dev.off()
	return(pi0)
}

deseq2<-function(obj_exp,outfolder,invname,pop){
	design2<-obj_exp[[1]]
	table.fil<-obj_exp[[2]]
	go.w.deseq<-obj_exp[[3]]
	go.w.deseq.sex<-obj_exp[[4]]

	if (go.w.deseq==1){

		#path<-paste(sep="",rootpath,pop,"/",outfolder,"/")
		path<-paste(sep="",path_res,invname,"/",type,"/",outfolder)
		dir.create(paste(sep="",path),showWarnings = FALSE)
		
		#dir.create(paste(sep="",path),showWarnings = FALSE)
		print(path)
		print(nrow(design))
		if (go.w.deseq.sex==0){
			dse <- DESeqDataSetFromMatrix(countData = table.fil, colData = design2,design = ~ sex + condition)
		}else{
			dse <- DESeqDataSetFromMatrix(countData = table.fil, colData = design2,design = ~ condition)
		}
		dse <- DESeq(dse,fitType="parametric",quiet=TRUE)
		save(dse,file=paste(sep="",path,"/dse"))
		rld <- rlogTransformation(dse, blind=TRUE)
		save(rld,file=paste(sep="",path,"/rld"))

		write.table(as.data.frame(table(design2)),paste(sep="",path,"/design.tab"),quote=F,sep="\t")

		comp<-resultsNames(dse)
		res <- results(dse)
		pi0<-vector()
		pvalues<-vector()
		for (nameres in comp[grep("condition",comp)]){
			res <- results(dse,nameres,independentFiltering=FALSE,cooksCutoff=FALSE)

			pi0<-save.results(res,path,nameres,pi0)
			pvalues[nameres]<-1
		}

		if (length(levels(design2$condition))>2){
			res.inv.het<-results(dse, contrast=c("condition","INV","HET"),
				independentFiltering=FALSE,cooksCutoff=FALSE)
			pi0<-save.results(res.inv.het,path,"condition_INV_vs_HET",pi0)
			pvalues["condition_INV_vs_HET"]<-1
		}

		print("Doing 50 permutations")
		if (pi0<1){
			pvalues<-send.random(obj_exp,path,pi0)
		}
		write.table(pvalues,paste(sep="",path,"/pvalues.tab"),quote=F,sep="\t")


	}

}