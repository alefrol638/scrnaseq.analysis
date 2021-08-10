#' @title import multiple DropSeq-RNA samples using parallelisation
#' @description import expression matrices ... adds metadata from cell names (need to provide the names of the columns) and percentage
#' of zeros for each sample under dataset[[]][,zeros]... it is possible to convert ensmus ids to gene names by setting ensmus=T
#' @param dir directory with the following structure: sample1/umi/matrix.mtx,genes.tsv,barcodes.tsv
#' @param meta_cols character vector, containing column names of the Metadata.
#' if the metadata is contained in the cell names, and separated by _, which information is contained at which position.
#' @param ensmus Boolean, should ensmus IDs be converted to gene symbols?
#' @param mart if ensmus=TRUE, which mart should be used ? (to load mart execute usemaRt())
#' @param genes_names character vector containing gene names ... if you do not want to take the genes names from umi/genes.tsv, set the gene names using this variable
#' @importFrom foreach %dopar%
#' @export read_mtx

read_mtx<-function(dir,meta_cols=c("Genotype","Mouse",NA,"Pool"),ensmus=F,mart,genes_names=NULL){
  samples<-list.files(dir)
  list_seurat <- list()

  list_seurat <- foreach::foreach(i=samples, .packages = c("Matrix","Seurat")) %dopar% {

    mtx <- Matrix::readMM(paste(dir ,i, "umi/matrix.mtx",sep="/"))
    genes <- read.delim(paste(dir,i,"umi/genes.tsv",sep="/"), header = F, stringsAsFactors = F)
    if(ensmus==T&is.null(genes_names)){

      genes_names<-biomaRt::getBM(attributes=c('external_gene_name', 'ensembl_gene_id'),
                            filters = 'ensembl_gene_id',
                            values = genes$V1,
                            mart = mart,
                            useCache = F )


      genes  <- dplyr::left_join(genes,genes_names,by=c("V1"="ensembl_gene_id"))

      }


    bc <- read.delim(paste(dir,i,"umi/barcodes.tsv",sep="/"), header = F, stringsAsFactors = F)
    bc$V1<-paste(bc$V1,i,sep="_")

    if(is.null(genes_names))
    {
      colnames(genes)<-"external_gene_name"
    mtx@Dimnames[[1]] <- as.character(genes$external_gene_name)
    }else
    {
      mtx@Dimnames[[1]] <- as.character(genes_names)
    }
    mtx@Dimnames[[2]] <- as.character(bc$V1)


    seurat <- Seurat::CreateSeuratObject(counts = mtx, project = i, min.cells = 3)
    percent_zero<-cal_percent_zero(seurat)
    meta<-data.frame(zeros=rep(percent_zero,length(colnames(seurat))),library=Idents(seurat),row.names=colnames(seurat),stringsAsFactors = FALSE)
    meta<-tidyr::separate(meta,library,meta_cols,remove=T)
    seurat<-Seurat::AddMetaData(seurat,meta)


  }

  # merge seurat list to one object
  seurat_total<-list_seurat[[1]]
  for(i in 2:length(list_seurat)){
    seurat_total <-merge(seurat_total,list_seurat[[i]])
  }

  return(seurat_total)
}
