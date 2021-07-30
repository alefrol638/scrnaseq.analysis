#' Publication ready plots using ggplot2 ... any other layer can be added to this function by +, when executing
#' Parameter: y and xlab: axis titles; plot: which plot should be used( basically all geom_ ... possible );
#' legend: if want to remove use "none";do.facet= should facetting be performed, if yes type in facet the formula as string,
#' f.e "~Tissues" make separate plots for each tissue
#' @export
draw_plot<-function(obj,ylab="y",xlab="x",plot=ggplot2::geom_point(),legend="right",do.facet=F,facet=NULL,x.fontsize=15,y.fontsize=15,
                    legend.fontsize=17,
                    facets.fontsize=18,facets.rotation=45,x.rotation=90,title.fontsize=16,legend.key.size=1.5,cluster_colors=NULL){
  obj+plot+
    ggplot2::theme(panel.background=ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank(),axis.line=ggplot2::element_line(colour="black"),
          plot.title = ggplot2::element_text(hjust = 0.5,size=18,lineheight = 1.2),
          legend.title = ggplot2::element_blank(),
          legend.position = legend,
          strip.text =  ggplot2::element_text(size =facets.fontsize,angle=facets.rotation),
          strip.background = ggplot2::element_blank(),
          legend.text = ggplot2::element_text(size = legend.fontsize),
          axis.text.x = ggplot2::element_text(angle = x.rotation, vjust = 0.5, hjust=1,size=x.fontsize),
          axis.text.y = ggplot2::element_text( vjust = 0.5, hjust=1,size=y.fontsize),
          legend.key.size = ggplot2::unit(legend.key.size,"line"),
          axis.title=ggplot2::element_text(size=title.fontsize,face="bold"))+
    ggplot2::labs(y=ylab,x=xlab)+
    {if(do.facet)ggplot2::facet_grid(stats::reformulate(facet))}+
  {if(is.null(cluster_colors))ggplot2::scale_fill_manual(values=cluster_colors)}
  # geom_text(aes(label =Freq),  position = position_stack(vjust = 0.5),size=2)
}
