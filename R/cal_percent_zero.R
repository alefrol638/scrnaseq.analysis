#' @title Calculate Dropout rate
#' @description calculates percentage of zeros for each object in a list of seurat objects
#' @param x Seurat Object
#' @export
cal_percent_zero <- function(x) {
  expr<-Seurat::GetAssayData(x,slot="data")
  # calculate the percentage of dropouts
  expr[is.na(expr)] <- 0
  percent_zero<-(sum(expr==0)/(nrow(expr)*ncol(expr)))*100
  return(percent_zero)
}
