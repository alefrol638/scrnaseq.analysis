#' @title Volcano Plots for single cell clusters
#' @description wrapper for enhancedVolcano, made specifically for DE gene tables in Seurat 4.0 (careful Seurat <4.0 has "avg_logFC"
#' column instead of "avg_log2FC"). Performs DEG testing for single clusters and plots them separately into Volcano plots.
#' @param x SeuratObject
#' @param Condition name of metadata column to test for DEG
#' @param logfc,min.pct Cutoffs for find DEG
#' @param cond1,cond2 Two conditions to plot in Volcano plot
#' @param legend Legend position see ggplot
#' @param titles if used in loop, what titles should the single plots have? If you just look at clusters,this would be levels(x$cluster).
#' @param do.col Should a color code be included?
#' @param lab.size,title.size,axis.font.size Change the fontsizes of the single graph elements
#' @param xlab,ylab Genelabel you want to appear in the plot
#' @param axis.ticks,axis.line Thickness of ticks and lines
#' @param dot.size size of dots representing each gene in plot
#' @param allvs1 should just marker genes for each cluster be found (see Seurat::FindAllMarkers())
#' @export
#'
 SCVolcano<-function(x,Condition,cond1,cond2,logfc=0.1,min.pct=0.1,do.col=F,
lab.size=3,title.size=12,axis.font.size=7,dot.size=1,
FC_cutoff=0.25,p_cutoff=5e-02,legend="none",xlab=NULL,ylab=NULL,axis.ticks=.3,axis.line =.3,allvs1=T){

dataset.markers<-list()

dataset.markers$cluster<-list()

if(allvs1==T){
  dataset.markers$cluster<-Seurat::FindAllMarkers(object = x,logfc.threshold = logfc,min.pct=min.pct,only.pos = T)
  dataset.markers$cluster<-base::split(dataset.markers$cluster,dataset.markers$cluster$cluster)

}else{
for (i in levels(Idents(x))){
  print(i)
  x2<-subset(x,idents=i)
  Idents(object = x2) <- x2[[]][,Condition]



  dataset.markers$cluster[[i]] <-Seurat::FindMarkers(object = x2,ident.1 = cond1,ident.2 = cond2,logfc.threshold = logfc,min.pct=min.pct)
  # colnames(dataset.markers$cluster[[i]])[6]<-"Condition"
  dataset.markers$cluster[[i]]$cluster<-as.factor(i)

}

Idents(x)<-x[[]][,Condition]
dataset.markers$cond_total<-Seurat::FindMarkers(x,ident.1 = cond1,ident.2 = cond2)

rm(x,x2)}


return(dataset.markers)
}
