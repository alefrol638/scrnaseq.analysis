#' @title Standard pipeline for seurat processing
#' @description Wrapper for Seurat tools. Allows to use different normalisation methods by changing one parameter. Furthermore, changes in dimensionality reduction can be adapted.
#'  Magic imputation not developed yet.
#' @param  res Float, resolution for clustering (see FindClusters)
#' @param  dim Sequence of integers, Number of dimensions to use (see ElbowPlot)
#' @param sct Boolean Perform SCTransform?
#' @param regress character vector for metadata columns,which should be regressed for (e.g remove batch effects)
#' @param alg clustering algorithm 1:Louvain 2:multilevel Louvain 3:SLM 4:Leiden(need to be installed)
#' @param low.features TRUE if you are clustering for low number of genes (e.g CITE or Flow cytometry)
#' @param umap "uwot"=R umap "umap-learn"= python package(faster, but needs to be installed)
#' @param red.space Character vector specifying the reduction space to use for FindNeigbors and RunUMAP
#' @param only.var Boolean, should only variable genes returned in scale.data of SCT Assay?
#' @param UMIs set to FALSE, if you have full transcripts instead of 3' 5' captured (f.e SmartSeq data)
#' @param scale should the expression be scaled to unit variance?
#' @param norm.method Which method to use for normalisation, See Seurat::NormalizeData()
#' @param graph graph name for the FindNeighbors() function
#' @export
seurat_flow<-function(x,res=0.7,dim=1:30,sct=T,norm=F,do.pca=T,regress=NULL,alg=1,low.features=F,umap="uwot",red.space="pca",only.var=T,UMIs=T,scale=T,
                      norm.method="LogNormalize",graph="SCT_snn"){

  if(sct==T){
    x<-Seurat::SCTransform(x,do.correct.umi = UMIs,do.scale = scale,vars.to.regress = regress,return.only.var.genes = only.var)
    graph<-"SCT_snn"
  }
  if(norm==T){
    x <- Seurat::NormalizeData(x,normalization.method=norm.method)
    graph<-"RNA_snn"
    if(low.features==F)
    {x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)}

    if(scale==T)
    {
      if(only.var==F | low.features==T){
        x<-Seurat::ScaleData(x,vars.to.regress = regress,
                             features=rownames(x))
      }else{x<-Seurat::ScaleData(x,vars.to.regress = regress)}
    }
  }


  if(do.pca==T){
    if(low.features){
    x<-Seurat::RunPCA(x,features=rownames(x),approx=F)
  }else{x<-Seurat::RunPCA(x,features=VariableFeatures(x))}
  }
  x  <- Seurat::FindNeighbors(x, dims = dim,reduction=red.space)
  x  <- Seurat::FindClusters(x,resolution = res,algorithm=alg,graph.name=graph)
  x  <- Seurat::RunUMAP(x ,dims=dim,umap.method=umap,reduction=red.space)
}
