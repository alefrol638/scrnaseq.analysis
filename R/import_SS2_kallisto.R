#' @title Expression Matrix from SS2 kallisto .h5 files
#' @description Function imports single .h5 transcript abundance files and summarises them into a gene expression matrix (values normalised to length and sequencing depth).
#' The ensembl IDs are then transformed into gene names and a Seurat object is created. The values need still to be normalised and log transformed, according to the maintainers of
#' the tximport package (see https://support.bioconductor.org/p/132550/).
#' @param dir Path to folder containing the transcript to genes list "tx2genes.csv" and a subdirectory called kallisto containing each .h5 file in a separate folder (folder name is
#' sample ID, e.g Plate and well number)
#' @param mart.use create Mart object using biomaRt::UseMart()
#' @export
import_SS2_kallisto<-function(dir,mart.use)
{
  samples<-list.files(paste(dir,"kallisto",sep = "/"))

  files<-paste(dir ,"kallisto",samples,"abundance.h5",sep="/")
  names(files)<-samples



  tx2gene <- readr::read_csv(file.path(dir, "tx2genes.csv"))

  txi <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene,countsFromAbundance="lengthScaledTPM")

  ###translate ensembl ID to gene symbol

  rownames(txi$counts)<-unlist(lapply(strsplit(as.character(rownames(txi$counts)), split = ".",fixed = T), `[[`, 1))


  genes_names <- biomaRt::getBM(attributes = c("external_gene_name",
                                               "ensembl_gene_id"), filters = "ensembl_gene_id",
                                values = rownames(txi$counts), mart = mart.use, useCache = F)


  ###keep only genes which were identified
  mtx<-as.data.frame(txi$counts)
  mtx$target_id<-rownames(mtx)
  genes_names$ensembl_gene_id<-make.names(genes_names$ensembl_gene_id, unique=TRUE)
  mtx$target_id<-make.names(mtx$target_id, unique=TRUE)
  mtx<-dplyr::inner_join(mtx,genes_names,by=c("target_id"="ensembl_gene_id"))
  #
  ###duplicated gene names found, thus add identifier (e,g ".1") if gene names duplicated
  rownames(mtx)<-make.names(mtx$external_gene_name, unique=TRUE)
  #remove gene names
  mtx<-mtx[,which(!colnames(mtx)%in%c("target_id","external_gene_name"))]
  ##set NAs to 0
  mtx[is.na(mtx)]<-0

  seurat_total<-Seurat::CreateSeuratObject(mtx,project="SS2_lymphocytes",min.cells = 3)

  return(seurat_total)
}
