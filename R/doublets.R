#' @title  Wrapper for chris-mcginnis-ucsf/DoubletsFinder
#' @description Doublet Finder wrapper for usage on Seurat objects. You need to separately perform QC and clustering on your samples/Conditions before applying this function.
#' @param x Seurat Object
#' @param dim integer sequence ... number of dimensions to use for Doublet finding
#' @param Ncpus integer, Number of cores to use for computation
#' @param pN numeric, metrics for calculation, see Doublet Finder vignette for Details
#' @param exp_doub numeric, expected percentage of doublets based on loaded#Cells/#Beads, formula: doublet rate= (#Cells/#Beads)^K*e^(-#Cells/#Beads)/K!, whereas K is 1, since we do not want to have more than one cell per bead
#' @export
doublets<-function(x,dim=1:30,sct=T,Ncpus=10,pN=0.25,exp_doub=0.014){

  sweep.res.list <- DoubletFinder::paramSweep_v3(x, PCs = dim, sct = sct,num.cores=Ncpus)
  sweep.stats <-  DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <-  DoubletFinder::find.pK(sweep.stats)
  # take the pK for the  maximal BCmetric (here 0.3) for further analysis
  # ggplot(bcmvn_WT,aes(x=pK,y=BCmetric))+geom_line(aes(group=1))

  pK<-as.numeric(levels(bcmvn$pK)[bcmvn[which.max(bcmvn$BCmetric),"pK"]])


  annotations<-Seurat::Idents(x)
  homotypic.prop <-  DoubletFinder::modelHomotypic(annotations)           #  ex: annotations <- seu_kidney@meta.data$ClusteringResults
  # the loading density was 20000 cells on 110000 beads(under poisson the expected doublet rate corresponds to 1.5 % ): doublet rate= (20000/110000)^K*e^(-20000/110000)/K!
  # K equals 1 since we do not want to have more than one cell per bead
  nExp_poi <- round(exp_doub*nrow(x@meta.data))

  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  # x <- doubletFinder_v3(x, PCs = dim, pN = pN, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  # find out the name of the pANN generated above
  # colnames(x@meta.data)
  x <-  DoubletFinder::doubletFinder_v3(x, PCs = dim, pN = pN, pK = 0.03, nExp = nExp_poi.adj, sct = T)
  return(x)
}
