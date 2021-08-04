#' @title Gene Biotype distribution plot
#' @description taken from scripts of Jonas Schrepping/DZNE
#' Shows a plot of the distribution of biotypes among the genes in a Seurat object
#' @param seurat_total Seurat object
#' @export
gene_types<-function(seurat_total){
  transcript_biotypes<-biomaRt::getBM(attributes=c('external_gene_name', 'transcript_biotype'),
                             filters = 'mgi_symbol',
                             values = rownames(seurat_total),
                             mart = mouse_ensembl)
  # calculate sum of reads and mean # of reads
  expr<-Seurat::GetAssayData(seurat_total,"data")
  genetypeexpr <- data.frame(SYMBOL=rownames(expr),
                             MEAN=rowMeans(expr),
                             SUM=rowSums(expr))
  colnames(transcript_biotypes)<-c("SYMBOL","TYPE")
  genetypeexpr<-dplyr::left_join(genetypeexpr,transcript_biotypes)
  # remove outlier with very low expression
  genetypeexpr<-genetypeexpr[genetypeexpr$SUM>1,]

  # filter for protein coding genes only
  return(genetypeexpr)
}
