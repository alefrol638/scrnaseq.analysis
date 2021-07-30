#' correct for ambient genes using SoupX...
#' x:processed SeuratObject; x_raw:SeuratObject containing bad quality cells (before QC processing);
#' cluster:Cluster identity for each cell (standard will take active Idents from SeuratObject x);
#' meta.cols:"columns from metadata to keep (see x[[]])... do not keep nCount_RNA and nFeature_RNA... needs to be recalculated
#' @export
remove_ambient<-function(x,x_raw,cluster=Idents(x),meta.cols=c(1,4:7))
{
  dataset<-x

  # Extract count matrix from raw and QC processed seurat objects
  total_counts_raw<-Seurat::GetAssayData(object = x_raw, slot = 'counts')
  total_counts_norm<-Seurat::GetAssayData(object = dataset, slot = 'counts')

  #both matrices have to have the same number of genes
  total_counts_raw<-total_counts_raw[rownames(total_counts_raw) %in% rownames(total_counts_norm),]


  #Create soupx object
  sc = SoupX::SoupChannel(total_counts_raw,total_counts_norm,)
  #take cluster from processed dataset
  sc =  SoupX::setClusters(sc,cluster)
  # estimate contaimation with ambient genes
  sc =  SoupX::autoEstCont(sc)
  # correct for ambient genes
  out_soup<-  SoupX::adjustCounts(sc)


  #  bugfix do not take all the metadata from uncorrected dataset ... otherwise QC metrics scewed (e.g nCount_RNA, nFeature_RNA)
  dataset_nosoup<-Seurat::CreateSeuratObject(out_soup,meta.data = dataset@meta.data[,meta.cols])
  return(dataset_nosoup)
}

