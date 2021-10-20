###SingleR
####cmmm:Macrophage of Sham,3M,7M,14M
####cmtmm:Macrophage of Sham,3M,3T,7M,7T,14M,14T
####
m<-FindAllMarkers(cmmm,min.pct = 0.2,logfc.threshold = 0.2,only.pos = T)
top10 <- m %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
ref<-cmmm
ref$label<-Idents(ref)
mydata<-cmtmm
sceM<-as.SingleCellExperiment(ref[gene,])
#sceM <- scater::logNormCounts(sceM)
a=data.frame()
for (i in 0:c(max(mydata$seurat_clusters%>%as.numeric())-1)) {
  sceG<-as.SingleCellExperiment(subset(mydata,ident=i)[gene,])
  #sceG <- sceG[,colSums(counts(sceG)) > 0]
  #sceG <- scater::logNormCounts(sceG)
  pred.grun <- SingleR(test=sceG, ref=sceM, labels=sceM$label, de.method="wilcox")
  b<-table(pred.grun$labels)%>%as.matrix()%>%as.data.frame()
  b<-t(b)%>%as.data.frame()
  a<-rbind.fill(a,b)
}
a[is.na(a)]<-0 
b<-a
for (i in 1:max(mydata$seurat_clusters%>%as.numeric())) {
  b[i,]<-scale(a[i,]%>%as.numeric())%>%as.numeric()
}
for (i in 1:max(mydata$seurat_clusters%>%as.numeric())) {
  b$type[i]<-names(b[i,]%>%which.max())
}
b$cluster<-0:c(max(mydata$seurat_clusters%>%as.numeric())-1)


for (i in 1:max(mydata$seurat_clusters%>%as.numeric())) {
  Idents(mydata,cell=colnames(subset(mydata,ident=b$cluster[i])))<-b$type[i]
}
