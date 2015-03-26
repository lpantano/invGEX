type<-"Genes"
pathYRI<-"/home/lpantano/projects/ibb/inputs/YRI"
pathCEU<-paste(sep="","/home/lpantano/projects/ibb/inputs/CEU/counts.",type,".tab")
pathTSI<-paste(sep="","/home/lpantano/projects/ibb/inputs/TSI/counts.",type,".tab")
#rootpath<-"/home/shareddata/Bioinformatics/iRNASeq/test/"
#MLPA<-"/home/lpantano/projects/ibb/inputs/v2_RESULTS_41INVERSIONS_7POPULATIONS_v2_04_02_2014.tab"
MLPA<-"/home/lpantano/projects/ibb/inputs/HsInv45.txt"
#MLPA<-"/home/lpantano/projects/ibb/inputs/HsInv06.txt"
parameter.pipeline<-"repeat"
path_res<-"/home/lpantano/projects/ibb/inputs/res/"
file.pos<-"/home/lpantano/projects/ibb/inputs/positions.all"


positions<-read.table("/home/lpantano/projects/ibb/inputs/positions.all",header=T)
variations<-read.table("/home/lpantano/projects/ibb/inputs/gene.var.hapmap.20865155.ensembl.sum.tab",header=T,row.names=1)
gene.pos<-read.table("/home/lpantano/projects/ibb/inputs/gene.pos.all.hg18.bed",sep="\t",header=T,row.names=4)
gene.sum<-read.table("/home/lpantano/projects/ibb/inputs/gene.summary.all2.tab",sep="\t",header=T,row.names=1)
