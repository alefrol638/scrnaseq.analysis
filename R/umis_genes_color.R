#' @title Scatterplot: #UMIs/cell vs. #genes/cell with color coded metadata
#' @description Creates a Scatterplot, where the x axis represents the number of genes per cell and the y axis the number of UMIs per cell. The color can be set accordingly.
#' usually percentage of mitochondrial genes per cells is set as color, so that low quality cells (dying or damaged cells) can be identified by low number of UMIs,genes and high percentage
#' of mitochondrial genes.
#' @param dataset Seurat Object
#' @param color Name of the metadata column to use for colors in the plot. If FALSE no colorcoding will be performed.
#' @param dot.size The size of the single cells, represented by dots, in the plot
#' @export
umis_genes_color<-function(dataset,color="percent.mito",dot.size=0.5,font.size=17,legend.font=14,legend.key=5){
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
  ugm<-draw_plot(ugm,"#UMIs per Cell","#Genes per Cell",ggplot2::geom_point(size=dot.size)
                 ,x.fontsize = font.size,title.fontsize =font.size, y.fontsize = font.size,legend.fontsize=legend.font,
                 legend.key.size=legend.key)

  # add number of genes and number of UMIs histograms at the fringes
  ggExtra::ggMarginal(ugm, type = "histogram",bins=150)
  }
