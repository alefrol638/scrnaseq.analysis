#' @title Ambient Genes Correction
#' @description correct for ambient genes using SoupX (github.com/constantAmateur/SoupX)
#' @param x processed SeuratObject
#' @param x_raw SeuratObject containing bad quality cells (before QC processing)
#' @param cluster Cluster identity for each cell (standard will take active Idents from SeuratObject x)
#' @param meta.cols columns from metadata to keep (see x[[]])... do not keep nCount_RNA and nFeature_RNA... needs to be recalculated
#' @param soupQuantile,tfidfMin Parameters for defining soup
#' @export
remove_ambient<-function(x,x_raw,cluster=Idents(x),meta.cols=c(1,4:7),soupQuantile=0.9,tfidfMin=1)
{

  ###bugfix for error "  Error in quantile.default(soupProf$est, soupQuantile) :
 # missing values and NaN's not allowed if 'na.rm' is FALSE "
  ## This error appeared for me only in 10x datasets
  dataset<-x
  # Extract count matrix from raw and QC processed seurat objects
  total_counts_raw<-Seurat::GetAssayData(object = x_raw, slot = 'counts')
  total_counts_norm<-Seurat::GetAssayData(object = dataset, slot = 'counts')

  #both matrices have to have the same number of genes
  total_counts_raw<-total_counts_raw[rownames(total_counts_raw) %in% rownames(total_counts_norm),]

  cluster=Idents(dataset)


  #Create soupx object
  sc = SoupX::SoupChannel(tod = total_counts_raw,toc=total_counts_norm, keepDroplets = TRUE)
  sc$tod<-total_counts_raw

  toc = sc$toc
  scNoDrops = SoupX::SoupChannel(toc, toc, calcSoupProfile = FALSE)
  # Calculate soup profile
  soupProf = data.frame(row.names = rownames(toc), est = rowSums(toc)/sum(toc), counts = rowSums(toc))
  scNoDrops = SoupX::setSoupProfile(scNoDrops, soupProf)
  sc =  SoupX::setClusters(scNoDrops,cluster)

  head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)

  #plotMarkerDistribution(sc)
  sc =  SoupX::autoEstCont(sc#,soupQuantile=0.7,tfidfMin=0.7
  )

  # correct for ambient genes
  out_soup<-  SoupX::adjustCounts(sc)



  #  bugfix do not take all the metadata from uncorrected dataset ... otherwise QC metrics scewed (e.g nCount_RNA, nFeature_RNA)
  dataset_nosoup<-Seurat::CreateSeuratObject(out_soup,meta.data = dataset@meta.data[,meta.cols])
  return(dataset_nosoup)
}

