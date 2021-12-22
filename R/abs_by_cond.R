#' @title Compare Cluster abundances across conditions
#' @description This function creates a barplot, with stacks representing the single cell clusters and the x-axis the condition, which you specify by
#' parameter "condition". The y-axis shows the relativ abundance relative to the total cell number (number above bars) in the respective condition.
#' The numbers in the stacks show the cell number in the respective cell cluster
#' @param seuratobject self explanantory
#' @param condition meta data variable you want to use, to compare the cluster abundances
#' @param use_cols named vector with the cluster identities as names and the value a custom color you want to use. If NULL standard color will be used.
#' @param font.total Font size of the total cell number (number above bar plot)
#' @param font.stack Font size of the cluster cell numbers in each stack
#' @param ctrl.group Name of the condition which should appear on the far left side of the plot
#' @param label.stack,label.total Should the number of cells be displayed in each stack (cluster), or above bar (total), respectively?
#' @param legend.key.size,legend.text.size,legend.position Appearance of legend (See ggplot2)
#' @export
abs_by_cond<-function(seuratobject,condition="Genotype",use_cols=NULL,font.total=3,legend.key.size=1,legend.text.size=4,legend.position="right",font.stack=1.5,axis.text.size=7,ctrl.group="WT",label.stack=T,label.total=T,title="Cluster Abundances")
{###extract frequencies of each seurat clusters splitted into conditions as table from Seurat object total
  table_cluster<-table(Idents(seuratobject),seuratobject@meta.data[[condition]])
  table_frame<- as.data.frame(table_cluster)
  colnames(table_frame)<-c("Seurat_cluster",condition,"Freq")
  if(!is.na(ctrl.group)){
  table_frame[[condition]]<-relevel(table_frame[[condition]],ctrl.group)
  }
  ####required for the total cell number labels in the plot later
  table_geno<-table(seuratobject[[condition]])
  table_geno<- as.data.frame(table_geno)
  table_geno$Seurat_cluster<-"act. MG"
  colnames(table_geno)<-c(condition,"Freq","Seurat_cluster")
  if(!is.na(ctrl.group)){
  table_geno[[condition]]<-relevel(table_geno[[condition]],ctrl.group)
  }




  # Columns you want to group by
  grp_cols <- condition

  # Convert character vector to list of symbols
  dots <- lapply(grp_cols, as.symbol)

  ##compare condition

  ###total cell numbers
  props_total_t<-table_geno%>%dplyr::group_by(.dots=dots)%>%dplyr::mutate(total=sum(Freq),prop=sum(Freq/sum(Freq)))

  ### cell percentage for each cluster

  props_total_g<-table_frame%>%dplyr::group_by(.dots=dots)%>%dplyr::mutate(prop=Freq / sum(Freq))
  props_total_g<-as.data.frame(props_total_g)


  plot_list_g<-list()
  plot_list_g<-sapply(colnames(props_total_g),gen_plots,data=props_total_g,list=plot_list_g,xvalue=condition,yvalue="prop")
  # names(plot_list_g)<-colnames(props_total_g)
  ##we need only the plot for the clusters splitted into conditions
  plot_list_g<-plot_list_g[c(1)]

  graph_list_g<-list()
  ###you need to have the same number colors for the clusters otherwise NA will be shown
  # my_cols2<-my_cols2[-c(17,18,19)]

  graph_list_g<-lapply(plot_list_g,scRNAseq.analysis::draw_plot,ylab="rel. Abundance",xlab="",plot=geom_col(),cluster_colors = use_cols)

  ###plot with percentages and absolute cell numbers across conditions and cell clusters
  cluster_abs<-graph_list_g$Seurat_cluster.col+ggplot2::theme(legend.key.size = unit(legend.key.size,"mm"),legend.text = ggplot2::element_text(size=legend.text.size),legend.position = legend.position,
                                                              axis.text.x = ggplot2::element_text(size=axis.text.size),axis.text.y = ggplot2::element_text(size=axis.text.size))+
    ggplot2::ggtitle(title)+
    {if(label.stack==T)ggplot2::geom_text(aes(label =paste(Freq,"cells",sep=" ")),  position = position_stack(vjust = 0.5),size=font.stack)}+
    {if(label.total==T)ggplot2::geom_text(data=props_total_t,aes_string(label ="total",y="prop",x=condition),   vjust = -.25,size=font.total)}
  return(cluster_abs)
}
