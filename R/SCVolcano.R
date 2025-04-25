#' @title DEGs single conditions
#' @description Performs DEG testing for single clusters as input for the function VolcanoPlot, to plot DEG for cells of a certain metadata group or clusters.
#' @param x SeuratObject
#' @param Condition name of metadata column to test for DEG
#' @param logfc,min.pct Cutoffs for find DEG
#' @param cond1,cond2 Two conditions which are needed to be compared to each other(only if allvs1 is FALSE)
#' @param allvs1 should just marker genes for each cluster be found (see Seurat::FindAllMarkers())
#' @param test which stastitical test to use: wilcoxon Rank sum, DESeq2, MAST, ... (see Seurat::FindAllMarkers())
#' @export
 SCVolcano<-function(x,Condition,cond1,cond2,logfc=0.1,min.pct=0.1,allvs1=T,test="wilcox",only.pos=F){

dataset.markers<-list()

dataset.markers$cluster<-list()

if(allvs1==T){
  ###find marker genes 1 vs rest
  dataset.markers$cluster<-Seurat::FindAllMarkers(object = x,logfc.threshold = logfc,min.pct=min.pct,only.pos = only.pos,test.use=test)

  ###split the whole DEG list into separate lists for each cluster
  dataset.markers$cluster<-base::split(dataset.markers$cluster,dataset.markers$cluster$cluster)

}else{
for (i in levels(Idents(x))){
  print(i)

  ## subset for one cluster
  x2<-subset(x,idents=i)
  ###set the idents to the meta information we want to compare the conditions from
  Idents(object = x2) <- x2[[]][,Condition]



  dataset.markers$cluster[[i]] <-Seurat::FindMarkers(object = x2,ident.1 = cond1,ident.2 = cond2,only.pos = only.pos,logfc.threshold = logfc,min.pct=min.pct,test.use=test)
  # colnames(dataset.markers$cluster[[i]])[6]<-"Condition"

  ###add cluster identity and gene name to the DEG list
  dataset.markers$cluster[[i]]$cluster<-as.factor(i)
  dataset.markers$cluster[[i]]$gene<-rownames(dataset.markers$cluster[[i]])

}


  ##comparison the two conditions, without considering single cluster (bulk DEG testing)
Idents(x)<-x[[]][,Condition]
dataset.markers$cond_total<-Seurat::FindMarkers(x,ident.1 = cond1,ident.2 = cond2,test.use=test,only.pos = only.pos)
dataset.markers$cond_total$cluster<-as.factor(i)
dataset.markers$cond_total$gene<-rownames(dataset.markers$cond_total)
}

return(dataset.markers)
}
