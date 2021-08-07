#' @title Gene Biotype distribution plot
#' @description taken from scripts of Jonas Schrepping/DZNE
#' Shows a plot of the distribution of biotypes among the genes in a Seurat object
#' @param seurat_total Seurat object
#' @param mart.use create Mart Object using UseMart()
#' @export
gene_types<-function(seurat_total,mart.use){
  transcript_biotypes<-biomaRt::getBM(attributes=c('external_gene_name', 'transcript_biotype',"chromosome_name"),
                             filters = 'external_gene_name',
                             values = rownames(seurat_total),
                             mart = mart.use)
  # calculate sum of reads and mean # of reads
  expr<-Seurat::GetAssayData(seurat_total,"counts")
  genetypeexpr <- data.frame(SYMBOL=rownames(expr),
                             MEAN=Matrix::rowMeans(expr),
                             SUM=Matrix::rowSums(expr))
  colnames(transcript_biotypes)<-c("SYMBOL","TYPE")
  genetypeexpr<-dplyr::left_join(genetypeexpr,transcript_biotypes)
  # remove outlier with very low expression
  genetypeexpr<-genetypeexpr[genetypeexpr$SUM>1,]

  # filter for protein coding genes only
  return(genetypeexpr)
}
