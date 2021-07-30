#' taken from scripts of Jonas Schrepping/DZNE
#' @export UMIs_cell
#'
UMIs_cell<-function(seurat_total,cutoff=10000){
  # extract the number of UMIs per cell from seurat object
  data_dge<-as.data.frame(seurat_total$nCount_RNA)
  colnames(data_dge)<-"UMIs"

  # rank the UMis per cell in descending order
  data_dge %>% dplyr::arrange(desc(UMIs))  -> data
  data$index <- 1:nrow(data)

  # plot the UMIs vs the ranked cells, cutoff line inserted at 40 UMIs per cell
  ggplot2::ggplot(data,ggplot2::aes(x=index,y=UMIs)) +
    ggplot2::geom_line() +  ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::ggtitle(paste("Reads per cell distribution (",nrow(data[data$UMIs>cutoff,])," cells)",sep="")) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, size = 12),
          panel.background = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(fill=NA, size=1))+
    ggplot2::geom_hline(yintercept=cutoff, linetype="dashed", color = "red")+
    ggplot2::annotation_logticks()  +
    ggplot2::ylab(label="Total number of sequenced reads")

}
