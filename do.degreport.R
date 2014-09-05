library(DESeq2)
library(HTSFilter)
library(DEGreport)
library(BiocParallel)


do_degreport<-function(tags, dse, res, geno, path_out){
    if (length(tags)>0){
        dir.create(path_out, showWarnings = FALSE)
        design <- colData(dse)
        c <- counts(dse, normalized =TRUE)
        g1<-row.names(design[design$condition == geno[1],])
        g2<-row.names(design[design$condition == geno[2],])
        fc <- res[tags,"log2FoldChange"]
        pvalue <- res$pvalue
        c <- c[row.names(res),]
        createReport(g1, g2, c, tags, pvalue, fc, path = path_out, ncores = 2)
        }
}


do.dge.summary<-function(invname){
    source("parameters.R")
    positions<-read.table("/home/lpantano/projects/ibb/inputs/positions.all",header=T)
    variations<-read.table("/home/lpantano/projects/ibb/inputs/gene.var.hapmap.20865155.ensembl.sum.tab",header=T,row.names=1)
    gene.pos<-read.table("/home/lpantano/projects/ibb/inputs/gene.pos.all.hg18.bed",sep="\t",header=T,row.names=4)
    gene.sum<-read.table("/home/lpantano/projects/ibb/inputs/gene.summary.all2.tab",sep="\t",header=T,row.names=1)
    comp<-vector("list")
    comp[[1]]<-c("cSTD","INV")
    comp[[2]]<-c("HET","INV")
    comp[[3]]<-c("cSTD","HET")
    population<-c("CEUTSI","YRI")
    labs<-c("CEUTSI.inv","CEUTSI.inv.het","CEUTSI.het","YRI.inv","YRI.inv.het","YRI.het")
    print(invname)
    list.g<-vector()
    for (pop in population){
        print(pop)
        for (i in 1:3) {
            gen<-paste0(comp[[i]],collapse="")
            print(gen)
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
                    if(class(res)=="DESeqResults"){
                        res<-res[!is.na(res$padj),]
                        res.all<-as.data.frame(mcols(dse,use.names=TRUE))
                        if (sum(res$padj<=0.10,na.rm=T)>0){
                            top<-as.data.frame(res[res$padj<=0.10,])
                            filtered<-intersect(row.names(top),row.names(res.all))
                            list.g<-unique(filtered)
                            do_degreport(list.g, dse, res, comp[[i]],
                                        paste0(path_res,invname,"/",type,"/",pop,".DESeq2",type,".parametric/",gen,"/"))
                    }
                }
            }
            
        }
        
    }
}
