library(DESeq2)
library(HTSFilter)

do.summary.add<-function(invname){
	source("parameters.R")
	#type<-"Genes"

	table.res<-data.frame()
	#all.inversions<-"HsInv58"
	population<-c("CEUTSI","YRI")
	labs<-c("CEUTSI","YRI")
		print(invname)
		list.g<-vector()
		for (pop in population){
      #print(pop)
				dse.dir<-paste(sep="",path_res,invname,"/",type,"/",pop,".DESeq2",type,".parametric/additivedse")
				rld.dir<-paste(sep="",path_res,invname,"/",type,"/",pop,".DESeq2",type,".parametric/addiiverld")
				if (file.exists(dse.dir)){
					load(dse.dir)
					dsef <- try(HTSFilter(dse, s.len=100, plot=FALSE)$filteredData,silent=T)
				 	if (class(dsef)!="try-error"){
				 		res<-try(results(dsef, independentFiltering=FALSE,cooksCutoff=FALSE),silent=TRUE)

				 	}else{
						res<-try(results(dse, independentFiltering=FALSE,cooksCutoff=FALSE),silent=TRUE)
					}
          #print(class(res))
					if(class(res)=="DESeqResults"){
						#print(gen)
						res<-res[!is.na(res$padj),]
						res.all<-as.data.frame(mcols(dse,use.names=TRUE))
						res.all<-res.all[!is.na(res.all$dispersion),]
						res<-res[abs(res$log2FoldChange) > 0.75,]
						#print(sum(res$padj<=0.20,na.rm=T))
						if (sum(res$padj<=0.20,na.rm=T)>0){
              #print(res)
							top<-as.data.frame(res[res$padj<=0.20,])
							filtered<-intersect(row.names(top),row.names(res.all))
							
							list.g<-unique(c(list.g,filtered))
						}
					}
				}
			
			
			
		}
		
    #print(list.g)
		idx<-0
		fdr.total<-data.frame(genes=list.g,row.names=list.g)
		cf.total<-data.frame(genes=list.g,row.names=list.g)
		pval.total<-data.frame(genes=list.g,row.names=list.g)
		rt.total<-data.frame(genes=list.g,row.names=list.g)
		done<-vector(length=3)
	if (length(list.g)>0){
		for (pop in population){
			#print(pop)
			idx.pop<-which(population==pop)
				idx<-idx.pop
				dse.dir<-paste(sep="",path_res,invname,"/",type,"/",pop,".DESeq2",type,".parametric/additivedse")
				rld.dir<-paste(sep="",path_res,invname,"/",type,"/",pop,".DESeq2",type,".parametric/additiverld")
				if (file.exists(dse.dir)){
					load(dse.dir)
					load(rld.dir)
					#print("1 if")
					dsef <- try(HTSFilter(dse, s.len=100, plot=FALSE)$filteredData,silent=T)
				 	if (class(dsef)!="try-error"){
				 		res<-try(results(dsef, independentFiltering=FALSE,cooksCutoff=FALSE),silent=TRUE)

				 	}else{
						res<-try(results(dse, independentFiltering=FALSE,cooksCutoff=FALSE),silent=TRUE)
					}
					design<-as.data.frame(colData(dse)[,1:2])
					min.indv<-min(table(design[,1]))

					if(class(res)=="DESeqResults" & min.indv>=4){
						#print("2 if")
						done[idx]<-1
						res<-res[!is.na(res$padj),]
						res.all<-mcols(dse,use.names=TRUE)
						res.all<-res.all[!is.na(res.all$dispersion),]
						#res.all<-res.all[res.all$dispersion<=2,]
						list.g.t<-intersect(list.g,row.names(res.all))
						
						exp<-as.data.frame(assay(rld[list.g.t,]))
						list.e<-sapply(list.g.t,function(gf){
							#print(gf)
							mean.e<-mean(as.numeric(exp[gf,]))
							q75.e<-quantile(as.numeric(exp[gf,]),.75)
							considered.exp<- mean.e>log2(2) & q75.e>log2(2)
							return(considered.exp)
						})
						names(list.e)<-list.g.t

						list.g.t<-intersect(list.g.t,names(list.e[list.e==TRUE]))
						res<-as.data.frame(res)
                          selected<-res[list.g.t,]
				
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
							

							rt.se<-res[row.names(selected),"lfcSE"]
							list.ratio<-abs(cf[,2])/((abs(cf[,2])+rt.se)-(abs(cf[,2])-rt.se))
							#list.ratio<-abs(cf[,2])-rt.se
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
				}else{
							done[idx:(idx+2)]<-0
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
				for (ni in 2:3){
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




		yridone<-sum(done[3])
		eudone<-sum(done[2])
		commoneu<-0
		commoneu.r<-0
	
			table<-as.data.frame(t(table))
			names(table)<-names(cf.total)
			
			idx.ratio<-grep("[::*::]",as.vector(as.matrix(table[,2:3])),invert=T)
			temp.vector<-as.vector(as.matrix(rt.total[,2:3]))
			temp.vector[idx.ratio]<-"-"

			rt.total<-as.data.frame(matrix(temp.vector,nrow=nrow(cf.total)))
			row.names(rt.total)<-row.names(fdr.total)
			names(rt.total)<-names(fdr.total)[2:3]


			#sc.min<-apply(rt.total,1,min)

			idx.ratio<-grep("[::**::]",as.vector(as.matrix(table[,2:3])),invert=T)
			temp.vector<-as.vector(as.matrix(rt.total))
			temp.vector[idx.ratio]<-10000
			# #temp.vector<-sub("[::**::]","",temp.vector)
			temp.vector<-abs(as.numeric(temp.vector))
			rt.lim<-as.data.frame(matrix(temp.vector,nrow=nrow(rt.total)))
			rt.min<-apply(rt.lim,1,min)
			# row.names(cf.lim)<-row.names(cf.total)
			# cf.lim$min<-cf.min

			table$POP<-factor("NONE",levels=c("NONE","BOTH","EU","YRI","EU->YRI","YRI->EU"))
			##tendency
			min.eu<-(apply(fdr.total[,2,drop=FALSE],1,function(x){return(which(x<=0.1))}))
			min.eu<-unlist(lapply(min.eu, length))
			min.yri<-(apply(pval.total[,3,drop=FALSE],1,function(x){return(which(x<=0.02))}))
			min.yri<-unlist(lapply(min.yri, length))

			res.f.r<-intersect(names(min.eu[min.eu>0]),names(min.yri[min.yri>0]))

			if (length(res.f.r)>0){
				table[res.f.r,"POP"]<-"EU->YRI"
			}


			min.eu<-(apply(pval.total[,2,drop=FALSE],1,function(x){return(which(x<=0.02))}))
			min.eu<-unlist(lapply(min.eu, length))
			min.yri<-(apply(fdr.total[,3,drop=FALSE],1,function(x){return(which(x<=0.1))}))
			min.yri<-unlist(lapply(min.yri, length))

			res.f.r<-intersect(names(min.eu[min.eu>0]),names(min.yri[min.yri>0]))

			if (length(res.f.r)>0){
				table[res.f.r,"POP"]<-"YRI->EU"
			}


			#both populations
			min.eu<-(apply(fdr.total[,2,drop=FALSE],1,function(x){return(which(x<=0.2))}))
			min.eu<-unlist(lapply(min.eu, length))
			min.yri<-(apply(fdr.total[,3,drop=FALSE],1,function(x){return(which(x<=0.2))}))
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

			table<-cbind(table,rt.min)
			table<-table[table$genes!="g1",]		
			table$inv<-invname

			table<-cbind(table,gene.pos[row.names(table),])
			
			posinv<-positions[invname,]
			table$chr<-as.character(table$chr)
			table$start<-as.numeric(table$start)
			table$end<-as.numeric(table$end)
			table$relative<-apply(table[,7:9],1,function(x){
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
			#table$mincf<-cbind(table,cf.lim$min)
			table.res<-rbind(table.res,table)
			#table.res<-cbind(table.res,rt.total[row.names(table.res),])
		}


	eunan<-intersect(rownames(table.res)[(table.res$CEUTSI.het.CF)=="-" | (table.res$CEUTSI.het.CF)=="-"],intersect(rownames(table.res)[(table.res$CEUTSI.inv.het.CF)=="-" | (table.res$CEUTSI.inv.het.CF)=="-"],rownames(table.res)[(table.res$CEUTSI.inv.CF)=="-" | (table.res$CEUTSI.inv.CF)=="-"]))
	yrinan<-intersect(rownames(table.res)[(table.res$YRI.het.CF)=="-"  | (table.res$YRI.het.CF)=="-" ],intersect(rownames(table.res)[(table.res$YRI.inv.het.CF)=="-" | (table.res$YRI.inv.het.CF)=="-"],rownames(table.res)[(table.res$YRI.inv.CF)=="-" | (table.res$YRI.inv.CF)=="-"]))

	both<-rownames(table.res[table.res$POP=="EU->YRI" | table.res$POP=="YRI->EU" | table.res$POP=="BOTH" ,])

	euuni<-setdiff(setdiff(rownames(table.res)[table.res$POP=="EU" ],eunan),both)
	yriuni<-setdiff(setdiff(rownames(table.res)[table.res$POP=="YRI" ],yrinan),both)

	table.both<-table.res[unique(c(eunan,yrinan,both)),]
	table.unique<-table.res[unique(c(euuni,yriuni)),]

	write.table(table.both,paste(sep="",path_res,invname,"/",type,"/table.add.both.tab"),sep="\t",quote=F,row.names=F)
	write.table(table.unique,paste(sep="",path_res,invname,"/",type,"/table.add.unique.tab"),sep="\t",quote=F,row.names=F)
	write.table(table.res,paste(sep="",path_res,invname,"/",type,"/table.add.raw.tab"),sep="\t",quote=F,row.names=F)
}
