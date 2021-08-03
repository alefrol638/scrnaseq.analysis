#' @title Create ggplot
#' @description wrapper for ggplot, in order to use ggplot for plotting multiple plots using lapply
#' @param col
#' @export
gen_plots<-function (col,data,list,xvalue,yvalue){
  list$col<-ggplot(data,ggplot2::aes_string(x=xvalue,y=yvalue,fill=col,label="Freq"))
  return(list)

}
