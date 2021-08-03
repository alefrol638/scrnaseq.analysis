#' Standard pipeline for seurat processing ...
#' @param  res Float, resolution for clustering (see FindClusters)
#' @param  dim Sequence of integers, Number of dimensions to use (see ElbowPlot)
#' @param sct Boolean Perform SCTransform?
#' @param magic Boolean, perform imputation by magic? ... not done yet
#' @param regress character vector for metadata columns,which should be regressed for (e.g remove batch effects)
#' @param alg clustering algorithm 1:Louvain 2:multilevel Louvain 3:SLM 4:Leiden(need to be installed)
#' @param low.features TRUE if you are clustering for low number of genes (e.g CITE or Flow cytometry)
#' @param umap "uwot"=R umap "umap-learn"= python package(faster, but needs to be installed)
#' @param red.space Character vector specifying the reduction space to use for FindNeigbors and RunUMAP
#' @param only.var Boolean, should only variable genes returned in scale.data of SCT Assay?
#'
#' @export
seurat_flow<-function(x,res=0.7,dim=1:30,sct=T,magic=F,norm=F,do.pca=T,regress=NULL,alg=1,low.features=F,umap="uwot",red.space="pca",only.var=T){

  if(magic==T){
    x <- Seurat::NormalizeData(x)
    x<-Seurat::ScaleData(x,vars.to.regress = regress)
    x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
  }
  if(sct==T){
    x<-Seurat::SCTransform(x,do.correct.umi = T,do.scale = T,vars.to.regress = regress,return.only.var.genes = only.var)
  }
  if(norm==T){
    x <- Seurat::NormalizeData(x)
    x<-Seurat::ScaleData(x,vars.to.regress = regress)
    x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
  }
  if(do.pca==T){
    x  <- Seurat::RunPCA(x , features = VariableFeatures(object = x ),approx=low.features)
  }
  x  <- Seurat::FindNeighbors(x, dims = dim,reduction=red.space)
  x  <- Seurat::FindClusters(x,resolution = res,algorithm=alg)
  x  <- Seurat::RunUMAP(x ,dims=dim,umap.method=umap,reduction=red.space)
}

