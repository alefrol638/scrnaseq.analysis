#' @title Publication ready plots using ggplot2
#' @description Any other layer can be added to this function by +, when executing
#' @param obj Data frame, which can be used in ggplot
#' @param ylab,xlab: axis titles; plot: which plot should be used( basically all geom_ ... possible )
#' @param plot Any ggplot (f.e geom_boxplot())
#' @param legend If want to remove use "none"
#' @param do.facet boolean, should facetting be performed? if yes type in facet the formula as string,
#' @param facet character, formula for facetting (subplots based on one variable, f.e "~Tissues" make separate plots for each tissue)
#' @param x.fontsize,y.fontsize,legend.fontsize,facets.fontsize,facets.rotation,x.rotation,title.fontsize,legend.key.size Formatting options sizes, text rotations ...
#' @param cluster_colors named character vector ... ColorIDs named with the cluster ID to use for labelling the clusters, if empty the standard ggplot colors will be used
#' @examples
#' draw_plot(data.frame,plot=geom_point(),legend="none",do.facet=T,facet=".~Genotype",cluster_colors=custom.colors)
#' draw_plot(data.frame,plot=geom_point(),legend="none",do.facet=T,facet="Mouse~Genotype",cluster_colors=custom.colors)
#'
#' @export
draw_plot<-function(obj,ylab="y",xlab="x",plot=ggplot2::geom_point(),legend="right",do.facet=F,facet=NULL,x.fontsize=7,y.fontsize=7,
                    legend.fontsize=7,
                    facets.fontsize=7,facets.rotation=45,x.rotation=90,title.fontsize=8,legend.key.size=2,cluster_colors=NULL){
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
          legend.key.size = ggplot2::unit(legend.key.size,"mm"),
          axis.title=ggplot2::element_text(size=title.fontsize,face="bold")
          )+
    ggplot2::labs(y=ylab,x=xlab)+
    {if(do.facet)ggplot2::facet_grid(stats::reformulate(facet))}+
  {if(!is.null(cluster_colors))ggplot2::scale_fill_manual(values=cluster_colors)}
  # geom_text(aes(label =Freq),  position = position_stack(vjust = 0.5),size=2)
}
