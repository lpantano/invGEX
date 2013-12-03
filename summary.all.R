do.summary<-function(invname){
	source("parameters.R")
	#type<-"Genes"
	positions<-read.table("/home/shareddata/Bioinformatics/iRNASeq/positions.all",header=T)
	variations<-read.table("/home/lpantano/iRNAseq/annotation/gene.var.hapmap.20865155.ensembl.sum.tab",header=T,row.names=1)
	gene.pos<-read.table("/home/lpantano/iRNAseq/annotation/gene.pos.all.hg18.bed",sep="\t",header=T,row.names=4)
	gene.sum<-read.table("/home/lpantano/iRNAseq/annotation/gene.summary.all2.tab",sep="\t",header=T,row.names=1)

	table.res<-data.frame()
	#all.inversions<-"HsInv58"
	comp<-vector("list")
	comp[[1]]<-c("INV","cSTD")
	comp[[2]]<-c("INV","HET")
	comp[[3]]<-c("HET","cSTD")
	population<-c("CEUTSI","YRI")
	labs<-c("CEUTSI.inv","CEUTSI.inv.het","CEUTSI.het","YRI.inv","YRI.inv.het","YRI.het")
		print(invname)
		list.g<-vector()
		for (pop in population){
			dse.dir<-paste(sep="",path_res,invname,"/",type,"/",pop,".DESeq2",type,".parametric/dse")
			rld.dir<-paste(sep="",path_res,invname,"/",type,"/",pop,".DESeq2",type,".parametric/rld")
			if (file.exists(dse.dir)){
				load(dse.dir)
				


				for (i in 1:3) {
					gen<-comp[[i]]
					res<-try(results(dse, contrast=c("condition",gen[1],gen[2]),
					independentFiltering=FALSE,cooksCutoff=FALSE),silent=TRUE)
					if(class(res)=="DataFrame"){
						#print(gen)
						res<-res[!is.na(res$padj),]
						res.all<-mcols(dse,use.names=TRUE)
						res.all<-res.all[!is.na(res.all$dispersion),]
						res.all<-res.all[res.all$dispersion<=1.5,]
						#print(sum(res$padj<=0.20,na.rm=T))
						if (sum(res$padj<=0.20,na.rm=T)>0){
							top<-as.data.frame(res[res$padj<=0.20,])
							filtered<-intersect(row.names(top),row.names(res.all))
							
							list.g<-unique(c(list.g,filtered))
						}
					}
				}
			}
			
		}
		

		idx<-0
		fdr.total<-data.frame(genes=list.g,row.names=list.g)
		cf.total<-data.frame(genes=list.g,row.names=list.g)
		pval.total<-data.frame(genes=list.g,row.names=list.g)
		rt.total<-data.frame(genes=list.g,row.names=list.g)
		done<-vector(length=7)
	if (length(list.g)>0){
		for (pop in population){
			idx.pop<-which(population==pop)
			base<-(idx.pop-1)*3
			dse.dir<-paste(sep="",path_res,invname,"/",type,"/",pop,".DESeq2",type,".parametric/dse")
			rld.dir<-paste(sep="",path_res,invname,"/",type,"/",pop,".DESeq2",type,".parametric/rld")
			
			if (file.exists(dse.dir)){
				load(rld.dir)
				load(dse.dir)
				for (i in 1:3) {
					idx<-i+base
					gen<-comp[[i]]
					
					res<-try(results(dse, contrast=c("condition",gen[1],gen[2]),
					independentFiltering=FALSE,cooksCutoff=FALSE),silent=T)
					design<-as.data.frame(colData(dse)[,1:2])
					min.indv<-min(nrow(design[design[,1]==gen[1],]),nrow(design[design[,1]==gen[2],]))
					if(class(res)=="DataFrame" & min.indv>=4){
						print(labs[idx])
						done[idx]<-1
						res<-res[!is.na(res$padj),]
						res.all<-mcols(dse,use.names=TRUE)
						res.all<-res.all[!is.na(res.all$dispersion),]
						res.all<-res.all[res.all$dispersion<=1.5,]
						list.g.t<-intersect(list.g,row.names(res.all))
						
						exp<-as.data.frame(assay(rld[list.g.t,]))
						list.e<-sapply(list.g.t,function(gf){
							#print(gf)
							mean.e<-mean(as.numeric(exp[gf,]))
							q75.e<-quantile(as.numeric(exp[gf,]),.75)
							considered.exp<- mean.e>log2(5) & q75.e>log2(5)
							return(considered.exp)
						})
						names(list.e)<-list.g.t

						list.g.t<-intersect(list.g.t,names(list.e[list.e==TRUE]))
						selected<-as.data.frame(res[list.g.t,])

						if (length(list.g.t)>0){
							fdr<-res[row.names(selected),6]
							fdr<-data.frame(genes=list.g.t,FDR=fdr)
							names(fdr)[2]<-paste(labs[idx],".FDR",sep="")
							fdr.total<-merge(fdr.total,fdr,by=1,all=T)

							pval<-res[row.names(selected),5]
							pval<-data.frame(genes=list.g.t,FDR=pval)
							names(pval)[2]<-paste(labs[idx],".pval",sep="")
							pval.total<-merge(pval.total,pval,by=1,all=T)

							cf<-res[row.names(selected),2]
							cf<-data.frame(genes=list.g.t,FDR=cf)
							names(cf)[2]<-paste(labs[idx],".CF",sep="")
							cf.total<-merge(cf.total,cf,by=1,all=T)

							
							
							
							list.ratio<-sapply(list.g.t,function(gf){
								ind.c1<-row.names(design[design$condition==gen[1],])
								ind.c2<-row.names(design[design$condition==gen[2],])
								gen1.q70<-round(quantile(exp[gf,ind.c1],0.70,na.rm=TRUE),digits=1)
								gen2.q30<-round(quantile(exp[gf,ind.c2],0.30,na.rm=TRUE),digits=1)
								gen1.q30<-round(quantile(exp[gf,ind.c1],0.30,na.rm=TRUE),digits=1)
								gen2.q70<-round(quantile(exp[gf,ind.c2],0.70,na.rm=TRUE),digits=1)
								exp.c1<-as.numeric(exp[gf,ind.c1])
								exp.c2<-as.numeric(exp[gf,ind.c2])
								min.over.1<-max(sum(exp.c1>=gen2.q70),sum(exp.c1>=gen2.q30))
								ratio.1<-min.over.1/length(exp.c1)
								min.over.2<-max(sum(exp.c2>=gen1.q70),sum(exp.c2>=gen1.q30))
								ratio.2<-min.over.2/length(exp.c2)
								ratio<-min(ratio.1,ratio.2)*100
								return(ratio)
							})


							rt<-res[row.names(selected),3]
							rt<-data.frame(genes=list.g.t,FDR=list.ratio)
							names(rt)[2]<-paste(labs[idx],".ratio",sep="")
							rt.total<-merge(rt.total,rt,by=1,all=T)
						}else{
							fdr<-data.frame(gene="g1",FDR=NA)
							names(fdr)[2]<-paste(labs[idx],".FDR",sep="")
							fdr.total<-merge(fdr.total,fdr,by=1,all=T)

							pval<-data.frame(gene="g1",FDR=NA)
							names(pval)[2]<-paste(labs[idx],".pval",sep="")
							pval.total<-merge(pval.total,pval,by=1,all=T)

							cf<-data.frame(gene="g1",FDR=NA)
							names(cf)[2]<-paste(labs[idx],".CF",sep="")
							cf.total<-merge(cf.total,cf,by=1,all=T)

							rt<-data.frame(gene="g1",ratio=NA)
							names(rt)[2]<-paste(labs[idx],".ratio",sep="")
							rt.total<-merge(rt.total,rt,by=1,all=T)
						}
					}else{
						done[idx]<-0
						fdr<-data.frame(gene="g1",FDR=NA)
						names(fdr)[2]<-paste(labs[idx],".FDR",sep="")
						fdr.total<-merge(fdr.total,fdr,by=1,all=T)

						pval<-data.frame(gene="g1",FDR=NA)
						names(pval)[2]<-paste(labs[idx],".pval",sep="")
						pval.total<-merge(pval.total,pval,by=1,all=T)

						cf<-data.frame(gene="g1",FDR=NA)
						names(cf)[2]<-paste(labs[idx],".CF",sep="")
						cf.total<-merge(cf.total,cf,by=1,all=T)

						rt<-data.frame(gene="g1",ratio=NA)
						names(rt)[2]<-paste(labs[idx],".ratio",sep="")
						rt.total<-merge(rt.total,rt,by=1,all=T)
					}
				}
			}else{
						idx<-base+1
						done[idx:(idx+2)]<-0
						fdr<-data.frame(gene="g1",FDR=NA,FDR2=NA,FDR3=NA)
						names(fdr)[2:4]<-paste(labs[idx:(idx+2)],".FDR",sep="")
						fdr.total<-merge(fdr.total,fdr,by=1,all=T)

						pval<-data.frame(gene="g1",FDR=NA,FDR2=NA,FDR3=NA)
						names(pval)[2:4]<-paste(labs[idx:(idx+2)],".pval",sep="")
						pval.total<-merge(pval.total,pval,by=1,all=T)

						cf<-data.frame(gene="g1",FDR=NA,FDR2=NA,FDR3=NA)
						names(cf)[2:4]<-paste(labs[idx:(idx+2)],".CF",sep="")
						cf.total<-merge(cf.total,cf,by=1,all=T)

						rt<-data.frame(gene="g1",FDR=NA,FDR2=NA,FDR3=NA)
						names(rt)[2:4]<-paste(labs[idx:(idx+2)],".ratio",sep="")
						rt.total<-merge(rt.total,rt,by=1,all=T)
			}
		}
	}
	if (nrow(fdr.total)>1){


		row.names(fdr.total)<-fdr.total[,1]
		row.names(pval.total)<-pval.total[,1]
		row.names(cf.total)<-cf.total[,1]
		row.names(rt.total)<-cf.total[,1]

			#cf.total[,1]<-as.character(cf.total[,1])
			table<-apply(cf.total,1,function(x){
				xr<-as.character(x)
				#print(xr)
				for (ni in 2:7){
					#print(fdr.total[xr[1],ni])	
					if (done[ni-1]==0){
						xr[ni]<-"-"
					}else if (is.na(fdr.total[xr[1],ni])){
						
						xr[ni]<-"NaN"
					}
					else if (fdr.total[xr[1],ni]<=0.2){
						xr[ni]<-paste(round(as.numeric(x[ni]),digits=2),"**")
					}else if (pval.total[xr[1],ni]<=0.02){
						xr[ni]<-paste(round(as.numeric(x[ni]),digits=2),"*")
					}else if (pval.total[xr[1],ni]>0.02){
						xr[ni]<-"NO"
					}
				}
				return(xr)

				})




		yridone<-sum(done[5:7])
		eudone<-sum(done[2:4])
		commoneu<-0
		commoneu.r<-0
	
			table<-as.data.frame(t(table))
			names(table)<-names(cf.total)
			idx.ratio<-grep("[::*::]",as.vector(as.matrix(table[,2:7])),invert=T)
			temp.vector<-as.vector(as.matrix(rt.total[,2:7]))
			temp.vector[idx.ratio]<-"-"
			rt.total<-as.data.frame(matrix(temp.vector,nrow=nrow(cf.total)))
			row.names(rt.total)<-row.names(fdr.total)
			names(rt.total)<-names(fdr.total)[2:7]

			idx.ratio<-grep("[::**::]",as.vector(as.matrix(table[,2:7])),invert=T)
			temp.vector<-as.vector(as.matrix(cf.total[,2:7]))
			temp.vector[idx.ratio]<-10000
			temp.vector<-sub("[::**::]","",temp.vector)
			temp.vector<-abs(as.numeric(temp.vector))
			cf.lim<-as.data.frame(matrix(temp.vector,nrow=nrow(cf.total)))
			cf.min<-apply(cf.lim,1,min)
			row.names(cf.lim)<-row.names(cf.total)
			cf.lim$min<-cf.min

			table$POP<-factor("NONE",levels=c("NONE","BOTH","EU","YRI","EU->YRI","YRI->EU"))
			##tendency
			min.eu<-(apply(fdr.total[,2:4],1,function(x){return(which(x<=0.1))}))
			min.eu<-unlist(lapply(min.eu, length))
			min.yri<-(apply(pval.total[,5:7],1,function(x){return(which(x<=0.02))}))
			min.yri<-unlist(lapply(min.yri, length))

			res.f.r<-intersect(names(min.eu[min.eu>0]),names(min.yri[min.yri>0]))

			if (length(res.f.r)>0){
				table[res.f.r,"POP"]<-"EU->YRI"
			}


			min.eu<-(apply(pval.total[,2:4],1,function(x){return(which(x<=0.02))}))
			min.eu<-unlist(lapply(min.eu, length))
			min.yri<-(apply(fdr.total[,5:7],1,function(x){return(which(x<=0.1))}))
			min.yri<-unlist(lapply(min.yri, length))

			res.f.r<-intersect(names(min.eu[min.eu>0]),names(min.yri[min.yri>0]))

			if (length(res.f.r)>0){
				table[res.f.r,"POP"]<-"YRI->EU"
			}


			#both populations
			min.eu<-(apply(fdr.total[,2:4],1,function(x){return(which(x<=0.2))}))
			min.eu<-unlist(lapply(min.eu, length))
			min.yri<-(apply(fdr.total[,5:7],1,function(x){return(which(x<=0.2))}))
			min.yri<-unlist(lapply(min.yri, length))

			res.f.r<-intersect(names(min.eu[min.eu>0]),names(min.yri[min.yri>0]))
			if (length(res.f.r)>0){
				table[res.f.r,"POP"]<-"BOTH"
			}
			#onepopulation
			res.f.r<-setdiff(names(min.eu[min.eu>0]),names(min.yri[min.yri>0]))
			if (length(res.f.r)>0){
				temp<-table$genes[table$POP=="NONE"]
				list.id<-intersect(temp,res.f.r)
				table[list.id,"POP"]<-"EU"
			}
			res.f.r<-setdiff(names(min.yri[min.yri>0]),names(min.eu[min.eu>0]))
			if (length(res.f.r)>0){
				temp<-table$genes[table$POP=="NONE"]
				list.id<-intersect(temp,res.f.r)
				table[list.id,"POP"]<-"YRI"
			}

			table<-table[table$genes!="g1",]		
			table$inv<-invname

			table<-cbind(table,gene.pos[row.names(table),])
			
			posinv<-positions[invname,]
			table$chr<-as.character(table$chr)
			table$start<-as.numeric(table$start)
			table$end<-as.numeric(table$end)
			table$relative<-apply(table[,10:12],1,function(x){
					start<-as.numeric(x[2])
					end<-as.numeric(x[3])
					#print(c(x))
					if (is.na(x[1])){
						dist<-"NaN"
					}else if (as.character(x[1])==as.character(posinv$chr[1])){
						#print(mean(start,end)-posinv$start[1])
						dist<-min(c(abs(mean(start,end)-posinv$start[1]),abs(mean(start,end)-posinv$end[1])))
						dist<-round(dist/1000000,digits=2)
					}else{
						dist<-"TRANS"
					}
					return(dist)
				})
			table<-cbind(table,variations[row.names(table),])
			table<-cbind(table,gene.sum[row.names(table),])
			table$eudone<-eudone
			table$yridone<-yridone
			table$mincf<-cbind(table,cf.lim$min)
			table.res<-rbind(table.res,table)
			table.res<-cbind(table.res,rt.total[row.names(table.res),])
		}


	eunan<-intersect(rownames(table.res)[(table.res$CEUTSI.het.CF)=="-" | (table.res$CEUTSI.het.CF)=="-"],intersect(rownames(table.res)[(table.res$CEUTSI.inv.het.CF)=="-" | (table.res$CEUTSI.inv.het.CF)=="-"],rownames(table.res)[(table.res$CEUTSI.inv.CF)=="-" | (table.res$CEUTSI.inv.CF)=="-"]))
	yrinan<-intersect(rownames(table.res)[(table.res$YRI.het.CF)=="-"  | (table.res$YRI.het.CF)=="-" ],intersect(rownames(table.res)[(table.res$YRI.inv.het.CF)=="-" | (table.res$YRI.inv.het.CF)=="-"],rownames(table.res)[(table.res$YRI.inv.CF)=="-" | (table.res$YRI.inv.CF)=="-"]))

	both<-rownames(table.res[table.res$POP=="EU->YRI" | table.res$POP=="YRI->EU" | table.res$POP=="BOTH" ,])

	euuni<-setdiff(setdiff(rownames(table.res)[table.res$POP=="EU" ],eunan),both)
	yriuni<-setdiff(setdiff(rownames(table.res)[table.res$POP=="YRI" ],yrinan),both)

	table.both<-table.res[unique(c(eunan,yrinan,both)),]
	table.unique<-table.res[unique(c(euuni,yriuni)),]

	write.table(table.both,paste(sep="",path_res,invname,"/",type,"/table.both.tab"),sep="\t",quote=F,row.names=F)
	write.table(table.unique,paste(sep="",path_res,invname,"/",type,"/table.unique.tab"),sep="\t",quote=F,row.names=F)
	write.table(table.res,paste(sep="",path_res,invname,"/",type,"/table.raw.tab"),sep="\t",quote=F,row.names=F)
}