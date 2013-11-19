deseq2<-function(obj_exp,outfolder,invname,pop){
	source("parameters.R")
	library(DESeq2)
	library(ggplot2)
	library(qvalue)
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
		resultsNames(dse)
		res <- results(dse)
		#head(res[order(res$pvalue),])
		resall<-as.data.frame(counts( dse, normalized=TRUE ))
	
		top10<-row.names((res[order(res$pvalue),])[1:20,])

		#write.table(as.data.frame(res[top10,]),paste(sep="",path,invname,"/top10.tab"),quote=F,sep="\t")
		write.table(as.data.frame(table(design2)),paste(sep="",path,"/design.tab"),quote=F,sep="\t")
		vector<-c(sum(res$padj<=0.1,na.rm=T),sum(res$padj<=0.2,na.rm=T),sum(res$padj<=0.3,na.rm=T),sum(res$padj<=0.4,na.rm=T))
		fdr.info<-data.frame(fdr=c(0.1,0.2,0.3,0.4),count=vector)
		write.table(fdr.info,paste(sep="",path,"/fdr.slides"),quote=F,row.names=F)
		
		qobj <- qvalue(res$pvalue[!is.na(res$pvalue)],fdr.level=0.1)
		qwrite(qobj,filename=paste(sep="",path,"/qvalue.slides"))
		t<-qsummary(qobj)
		pi0<-t$pi0
		num.sign.exp<-sum(qobj$significant)

		pdf(paste(sep="",path,"/qvalue.pdf"))
			qplot.log<-try(qplot(qobj),silent=TRUE)
		dev.off()

		save(dse,file=paste(sep="",path,"/dse"))

	#	pdf(paste(sep="",path,invname,"/pca.pdf"))
	#	vsd <- varianceStabilizingTransformation(dse)
	#	plotPCA(vsd)
	#	dev.off()

		# Overlaid histograms
		pdf(paste(sep="",path,"/pvalue.fdr.dist.pdf"))
			hist(res$pvalue,col=rgb(176/255,224/255,230/255,0.5),freq=FALSE,ylim=c(0,2),main="statistic distribution")
			hist(res$padj,col=rgb(102/255,153/255,204/255,0.5),add=T,freq=FALSE)
			legend("top",fill=c(rgb(176/255,224/255,230/255,0.5),rgb(102/255,153/255,204/255,0.5)),legend=c("pvalue","FDR"))
		dev.off()

		# disp.t<-data.frame(name=row.names(as.data.frame(colData(dse))),disp=as.data.frame(colData(dse))[,3])
		# if (info[1,4]=="YRI"){
		# 	do.graphics.wdisp(top10,paste(sep="",path,invname),resall[top10,],inv.names,inv.sex,as.data.frame(res[top10,]),disp.t)
		# }else{
		# 	do.graphics.wpop(top10,paste(sep="",path,invname),resall[top10,],inv.names,inv.sex,as.data.frame(res[top10,]),disp.t)
		# }

		if (pi0<1){
			print("doing 50 permutations")
			# source("~/scripts/iRNAseq/DE/deseq.random.R")
			# pval.pi0=round(sum(pi0>=list.qvalues)/51,digits=3)
			# write.table(paste(pi0,pval.pi0),paste(sep="",path,invname,"/pi0.pvalue"),quote=F,row.names=F)
			# pdf(paste(sep="",path,invname,"/qvalue.random.dist.pdf"))
			# 	h<-hist(list.qvalues,border="steelblue2",col="steelblue",main="Random pi0 distribution")
			# 	abline(v=pi0,col="red")
			# 	text(pi0,max(h$counts)-2,pos=4,labels=round(pval.pi0,digits=3))
			# dev.off()
			# pval.num.sign=round(sum(num.sign.exp<list.num.sign)/51,digits=3)
			# write.table(paste(num.sign.exp,pval.num.sign),paste(sep="",path,invname,"/num.sign.pvalue"),quote=F,row.names=F)
		}
		#disp.t<-data.frame(name=row.names(as.data.frame(colData(dse))),disp=as.data.frame(colData(dse))[,3])
		#do.graphics(top10,paste(sep="",path,invname),resall[top10,],inv.names,inv.sex,as.data.frame(res[top10,]))

	}

}