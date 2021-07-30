#' calculates percentage of zeros for each object in a list of seurat objects
#' @export
cal_percent_zero <- function(x) {
  expr<-Seurat::GetAssayData(x,slot="data")
  # calculate the percentage of dropouts
  expr[is.na(expr)] <- 0
  percent_zero<-(sum(expr==0)/(nrow(expr)*ncol(expr)))*100
  return(percent_zero)
}
