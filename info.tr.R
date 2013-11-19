
info.tr<-function(info.a){
 info.a<-info.a[info.a[,2]!="ND",]
 info.a[,2]<-factor(info.a[,2],levels=c("INV","HET","STD"))
 return(info.a)
}