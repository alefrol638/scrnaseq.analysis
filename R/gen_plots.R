#' @title Create ggplot
#' @description wrapper for ggplot, in order to use ggplot for plotting multiple plots using lapply
#' @param col Variable in data to use for color coding
#' @param data Dataframe usable with ggplot
#' @param list List, where to save the ggplots into
#' @param xvalue,yvalue see ggplot documentation
#' @export
gen_plots<-function (col,data,list,xvalue,yvalue){
  list$col<-ggplot(data,ggplot2::aes_string(x=xvalue,y=yvalue,fill=col,label="Freq"))
  return(list)

}
