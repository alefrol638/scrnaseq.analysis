#' @title Signifinance Bars for Seurat's Violin Plots
#' @description Takes adjusted p value from the output of FindMarkers and adds p values into the Violin plots (created by VlnPlot) of the respective gene.
#' @param SeuratObject on object of the Seurat class
#' @param gene Gene for which the significance test was performed
#' @param Condition name of metadata column to test for DEG
#' @param DEGs List of Differentially Expressed genes with their respective p values (output from FindMarkers)
#' @param split column name in meta used for split.by in Violin plot
#' @param assay which Seurat assay to use
#' @param slot which slot in the assay to use, possible: "counts", "data" or "scale.data"
#' @param log are the values in Violin plot represented in a logarithmic scale
#' @export
Vln_pvalue<-function(SeuratObject,gene,Condition,DEGs,split="filler",assay="RNA",slot="counts",log=T){

  dataset.markers<-list()
  ###create test object, with positions for p values in plot
  count.data <- GetAssayData(object =SeuratObject[[assay]], slot = slot)

  if(slot=="counts"){
  count.data <- as.matrix(x = count.data + 1)
  }
  if(log==T){
    count.data<-log(count.data,10)
  }
  vln_df = data.frame(expr = count.data[gene,], cluster = Idents(SeuratObject),Condition = SeuratObject[[Condition]])
  names(vln_df)[3]<-"Condition"
  # ggplot(vln_df,aes(x=expr))+geom_histogram(aes(color=genotype),fill="white")+theme_prism()

  ##first calculate p value with rstatix, to get a rstatix test class and thus be able to calculate the position of the p values in the plot
  vln_signif <- vln_df%>%
    rstatix::group_by(cluster) %>%
    rstatix::wilcox_test(expr ~ Condition) %>%
    rstatix::add_xy_position(x = "cluster", dodge = 0.8) # important for positioning!


  #replace the p value calculated by the simple wilcoxon rank sum test with the p values received from FindMarkers(recommended to use MAST test method)
  dataset.markers$selected_genes[[gene]]<-DEGs[grep(gene,DEGs$gene,fixed=T),]

  if(nrow(dataset.markers$selected_genes[[gene]])>0){
    ###add the p values which were found in DEGs
    vln_signif$p.format[vln_signif$cluster%in%dataset.markers$selected_genes[[gene]]$cluster]<-
      dataset.markers$selected_genes[[gene]]$p.format
    vln_signif$p[vln_signif$cluster%in%dataset.markers$selected_genes[[gene]]$cluster]<-
      dataset.markers$selected_genes[[gene]]$p_val
    vln_signif$p.adj[vln_signif$cluster%in%dataset.markers$selected_genes[[gene]]$cluster]<-
      dataset.markers$selected_genes[[gene]]$p_val_adj

    ###replace any genes not found in the DEGs with 1 ( not significant)
    vln_signif$p.format[!vln_signif$cluster%in%dataset.markers$selected_genes[[gene]]$cluster]<-1
    vln_signif$p[!vln_signif$cluster%in%dataset.markers$selected_genes[[gene]]$cluster]<-1
    vln_signif$p.adj[!vln_signif$cluster%in%dataset.markers$selected_genes[[gene]]$cluster]<-1

    vln_signif<-vln_signif  %>% rstatix::add_significance(p.col = "p.format")


    ##because VlnPlot sets the fill aesthetics globally, the variable for fill needs to be present, even if not used.
    vln_signif$split<-split

    dataset.markers$selected_genes[[gene]]<-vln_signif

    return(dataset.markers$selected_genes[[gene]])
  }else{
    ###the gene is not found in the marker gene list from FindMarkers
    return("No DEG found. Please check FC cutoff and min.pct when generating DEG list.")
  }
}
