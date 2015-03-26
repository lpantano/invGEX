library(DESeq2)
library(HTSFilter)

do.summary<-function(invname){
    source("parameters.R")
    #type<-"Genes"
    positions<-read.table("/home/lpantano/projects/ibb/inputs/positions.all",header=T)
    variations<-read.table("/home/lpantano/projects/ibb/inputs/gene.var.hapmap.20865155.ensembl.sum.tab",header=T,row.names=1)
    gene.pos<-read.table("/home/lpantano/projects/ibb/inputs/gene.pos.all.hg18.bed",sep="\t",header=T,row.names=4)
    gene.sum<-read.table("/home/lpantano/projects/ibb/inputs/gene.summary.all2.tab",sep="\t",header=T,row.names=1)
    
    table.res<-data.frame()
    #all.inversions<-"HsInv58"
    comp<-vector("list")
    comp[[1]]<-c("cSTD","INV")
    comp[[2]]<-c("HET","INV")
    comp[[3]]<-c("cSTD","HET")
    population<-c("CEUTSI","YRI")
    labs<-c("CEUTSI.inv","CEUTSI.inv.het","CEUTSI.het","YRI.inv","YRI.inv.het","YRI.het")
    print(invname)
    list.g<-vector()
    for (pop in population){
        #print(pop)
        for (i in 1:3) {
            gen<-paste0(comp[[i]],collapse="")
            dse.dir<-paste(sep="",path_res,invname,"/",type,"/",pop,".DESeq2",type,".parametric/",gen,"dse")
            rld.dir<-paste(sep="",path_res,invname,"/",type,"/",pop,".DESeq2",type,".parametric/",gen,"rld")
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
                    if (sum(res$padj<=0.20,na.rm=T)>0){
                        #print(res)
                        top<-as.data.frame(res[res$padj<=0.20,])
                        filtered<-intersect(row.names(top),row.names(res.all))
                        
                        list.g<-unique(c(list.g,filtered))
                    }
                }
            }
            
        }
        
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