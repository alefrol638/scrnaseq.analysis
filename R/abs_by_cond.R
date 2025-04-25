#' @title Compare Cluster abundances across conditions
#' @description This function creates a barplot, with stacks representing the single cell clusters and the x-axis the condition, which you specify by
#' parameter "condition". The y-axis shows the relativ abundance relative to the total cell number (number above bars) in the respective condition.
#' The numbers in the stacks show the cell number in the respective cell cluster
#' @param seuratobject self explanantory
#' @param condition1,condition2 meta data variable you want to use, to compare the cluster abundances, condition 1 will be split on the x axis, condition2 will be used for facetting
#' @param use_cols named vector with the cluster identities as names and the value a custom color you want to use. If NULL standard color will be used.
#' @param font.total Font size of the total cell number (number above bar plot)
#' @param font.stack Font size of the cluster cell numbers in each stack
#' @param ctrl.group1,ctrl.group2 Name of the condition which should appear on the far left side of the plot, refers to condition1 or condition2, respectively
#' @param label.stack,label.total Should the number of cells be displayed in each stack (cluster), or above bar (total), respectively?
#' @param legend.key.size,legend.text.size,legend.position Appearance of legend (See ggplot2)
#' @param abs if True, absolute numbers will be displayed in the y-axis. Otherwise percentages will be displayed
#' @export
#'
#'

abs_by_cond<-function(seuratobject,condition1,condition2=NA,multi_cond=F,use_cols=NULL,font.total=3,legend.key.size=1,legend.text.size=4,legend.position="right",font.stack=1.5,axis.text.size=7,ctrl.group1=NA,ctrl.group2=NA,label.stack=T,label.total=T,title="Cluster Abundances",abs=T)
{###extract frequencies of each seurat clusters splitted into conditions as table from Seurat object total
  if(multi_cond){
  table_cluster<-table(Idents(seuratobject),seuratobject@meta.data[[condition1]],seuratobject@meta.data[[condition2]])
  }else{ table_cluster<-table(Idents(seuratobject),seuratobject@meta.data[[condition1]])
}
  table_frame<- as.data.frame(table_cluster)
  if(multi_cond){
  colnames(table_frame)<-c("Seurat_cluster",condition1,condition2,"Freq")
  }else{colnames(table_frame)<-c("Seurat_cluster",condition1,"Freq")}
  if(!is.na(ctrl.group1)){
  table_frame[[condition1]]<-relevel(table_frame[[condition1]],ctrl.group1)
  }

  if(!is.na(ctrl.group2)){
    table_frame[[condition2]]<-relevel(table_frame[[condition2]],ctrl.group2)
  }
  ####required for the total cell number labels in the plot later
  if(multi_cond){
  table_geno<-table(seuratobject@meta.data[[condition1]],seuratobject@meta.data[[condition2]])
  table_geno<- as.data.frame(table_geno)
  table_geno$Seurat_cluster<-""

  colnames(table_geno)<-c(condition1,condition2,"Freq","Seurat_cluster")
  }else{table_geno<-table(seuratobject@meta.data[[condition1]])
  table_geno<- as.data.frame(table_geno)
  table_geno$Seurat_cluster<-""

  colnames(table_geno)<-c(condition1,"Freq","Seurat_cluster")}

  if(!is.na(ctrl.group1)){
  table_geno[[condition1]]<-relevel(table_geno[[condition1]],ctrl.group1)
  }

  if(!is.na(ctrl.group2)){
    table_geno[[condition2]]<-relevel(table_geno[[condition2]],ctrl.group2)
  }





  # Columns you want to group by
  if(multi_cond){
  grp_cols <- c(condition1,condition2)
}else{grp_cols <- c(condition1)}
  # Convert character vector to list of symbols
  dots <- lapply(grp_cols, as.symbol)

  ##compare condition

  ###total cell numbers
  props_total_t<-table_geno%>%dplyr::group_by(.dots=dots)%>%dplyr::mutate(total=sum(Freq),prop=sum(Freq/sum(Freq)))

  ### cell percentage for each cluster

  props_total_g<-table_frame%>%dplyr::group_by(.dots=dots)%>%dplyr::mutate(prop=Freq / sum(Freq))
  props_total_g<-as.data.frame(props_total_g)


  plot_list_g<-list()
  if(abs!=T){
  plot_list_g<-sapply(colnames(props_total_g),gen_plots,data=props_total_g,list=plot_list_g,xvalue=condition1,yvalue="prop")}
  if(abs==T){
    plot_list_g<-sapply(colnames(props_total_g),gen_plots,data=props_total_g,list=plot_list_g,xvalue=condition1,yvalue="Freq")
  }
  # names(plot_list_g)<-colnames(props_total_g)
  ##we need only the plot for the clusters splitted into conditions
  plot_list_g<-plot_list_g[c(1)]

  graph_list_g<-list()
  ###you need to have the same number colors for the clusters otherwise NA will be shown
  # my_cols2<-my_cols2[-c(17,18,19)]
  graph_list_g<-do.call(lapply,c(list(X=plot_list_g,
                        FUN=draw_plot,
                        ylab="rel. Abundance",xlab="",
                        plot=geom_bar(position="stack", stat="identity"),#.width = 0.2, #stat='bin',
                                     # position=position_dodge(1.5)),
                        cluster_colors = use_cols
                        ),
                   list(do.facet=T,facet=paste("~",condition2,sep=""))[multi_cond]
                        ))
  # graph_list_g<-lapply(plot_list_g,scRNAseq.analysis::draw_plot,ylab="rel. Abundance",xlab="",
  #                      plot=geom_col(),cluster_colors = use_cols,do.facet=T,facet=paste("~",condition2,sep=""))

  ###plot with percentages and absolute cell numbers across conditions and cell clusters
  cluster_abs<-graph_list_g$Seurat_cluster.col+ggplot2::theme(legend.key.size = unit(legend.key.size,"mm"),legend.text = ggplot2::element_text(size=legend.text.size),legend.position = legend.position,
                                                              axis.text.x = ggplot2::element_text(size=axis.text.size),axis.text.y = ggplot2::element_text(size=axis.text.size))+
    ggplot2::ggtitle(title)+ggprism::theme_prism()+
    {if(!is.null(use_cols))scale_fill_manual(values=use_cols)}+
    {if(label.stack==T)ggplot2::geom_text(aes(label =paste(Freq,"cells",sep=" ")),  position = position_stack(vjust = 0.5),size=font.stack)}+
    {if(label.total==T)ggplot2::geom_text(data=props_total_t,aes_string(label ="total",y="prop",x=condition1),   vjust = -.25,size=font.total)}
  return(cluster_abs)
}
