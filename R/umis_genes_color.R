#' @export
umis_genes_color<-function(dataset,color="percent.mito"){
  fvsu_mito<-as.data.frame(dataset$nCount_RNA)
  fvsu_mito$nFeature<-dataset$nFeature_RNA
  if(color!=F){
    fvsu_mito[,color]<-dataset[[]][,color]
    colnames(fvsu_mito)<-c("nUMIs","nGenes",color)

    ugm<-ggplot2::ggplot(fvsu_mito, ggplot2::aes_string(x="nGenes",y= "nUMIs",color=color))
  }
  else{
    colnames(fvsu_mito)<-c("nUMIs","nGenes")
    ugm<-ggplot2::ggplot(fvsu_mito, ggplot2::aes_string(x="nGenes",y= "nUMIs"))
  }
  ugm<-draw_plot(ugm,"#UMIs per Cell","#Genes per Cell",ggplot2::geom_point())

  # add number of genes and number of UMIs histograms at the fringes
  ggExtra::ggMarginal(ugm, type = "histogram",bins=150)
  }
