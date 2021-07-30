#' wrapper for enhancedVolcano, made specifically for DE gene tables in Seurat 4.0 (careful Seurat <4.0 has "avg_logFC"
#' column instead of "avg_log2FC")
#' @export
Volcano_plot<-function(x,Condition,p_cutoff=1e-05,FC_cutoff=1,legend="none",titles=levels(x$cluster),do.col=F){
  # define custom color for conditions
  if(do.col==T)
  {
    geno_col <- x[,Condition]
    names(geno_col)<-geno_col
    levels(geno_col)<-my_cols[1:length(levels(geno_col))]
  }

  plot<-do.call(EnhancedVolcano::EnhancedVolcano,args=c(list(toptable=x,
                                            lab=rownames(x),
                                            x="avg_log2FC",
                                            y="p_val_adj",
                                            labSize=7,
                                            title=titles,
                                            titleLabSize=12,
                                            subtitle=NULL,
                                            caption=NULL,
                                            legendLabSize = 7,
                                            legendIconSize = NULL,
                                            legendLabels = element_blank(),
                                            legendPosition = legend,
                                            axisLabSize = 7,
                                            pCutoff = p_cutoff,
                                            FCcutoff = FC_cutoff),
                                       list(colCustom = geno_col)[do.col]
  ))

}
