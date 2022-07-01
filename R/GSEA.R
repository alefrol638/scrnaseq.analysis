#' @title Gene Set Enrichment Analysis of scRNAseq data
#' @description Performs a GSEA for clusters, Conditions and Conditions between clusters. The GO,KEGG,DO, MSigDb databases are used.
#' For DO and MSigDB mouse genes are transformed to human homologues.
#' OrgDb = org.Mm.eg.db::org.Mm.eg.db need to execute after every restart, otherwise Error: external pointer is not valid.
#' Furthermore the Marts for translation from gene symbol to ensembl ID and finding homologues need to be loaded prior using this function (for usage see script in https://gitlab.dzne.de/frolova/spg15.git/scripts/Seurat_Flow_AF.Rmd)
#' the purpose of this function is to use it in foreach loop, in order to look for enrichment in all available databases in parallel (for usage see script in https://gitlab.dzne.de/frolova/spg15.git/scripts/Seurat_Flow_AF.Rmd).
#' For MSigDb GSEA To use the most up to date databases you might need to download the latest databases from MsigDb (https://www.gsea-msigdb.org/gsea/msigdb/).
#' The Databases used here were downloaded in July 2021.
#' Function adapted from scripts of Jonas Schulte-Schrepping/DZNE.
#' @param object Seurat Object
#' @param genelist A custom gene list to perform GOEA on (f.e from GRN analysis). Needs to contain at least cluster annotation. required columns: gene, cluster (cluster annotation) optional:Condition
#' @param condition name of column in metadata to use
#' @param top Number of genes from FindAllMarkers to use for GSEA
#' @param GeneSets which database to search (possible: GO,DO,KEGG,(from MSigDb:) Hallmark,cannonicalPathways,ImmunoSignatures,(Transcription Factor)Motifs)
#' @param Gontology more information in clusterprofiler vignette
#' @param  pCorrection Choose the p-value adjustment method
#' @param pvalueCutoff,qvalueCutoff set the unadj. or adj. p-value (depending on correction method) and q-value cutoff (FDR corrected)
#' @param min.pct,logfc.threshold parameters for FindAllMarkers (see manual for further information)
#' @param total If TRUE, will perform bulk GSEA, not considering single seurat cluster (useful for total comparison between samples, or just for marker genes of the clusters)
#' @param org Specifying the species: mmu: mouse hsa: human. More information in clusterprofiler vignette
#' @param present_genes=background,# character vector, Background genes used for statistical testing: (universe in enricher() function), gene names need first to be translated into ensemblID using biomaRt
#' @param OrgDb mus musculus database need to be loaded prior to use
#' @param human_ensembl,mouse_ensembl use useMart("ensembl",dataset="...") to load the human and mouse databases
#' @param up should the upregulated or downregulated genes be used?
#' @export
GSEA <- function(object=NULL,
                 genelist=NULL,
                   condition ="Genotype",#name of column in metadata to use
                   top=50, #top genes from FindAllMarkers
                   GeneSets ="GO",# which database to search (possible: GO,DO,KEGG,(from MSigDb:) Hallmark,cannonicalPathways,ImmunoSignatures,(Transcription Factor)Motifs)
                   Gontology = "BP",#more information in clusterprofiler vignette
                   pCorrection = "bonferroni", # choose the p-value adjustment method
                   pvalueCutoff = 0.05, # set the unadj. or adj. p-value cutoff (depending on correction method)
                   qvalueCutoff = 0.05, # set the q-value cutoff (FDR corrected)
                   min.pct=0.1,# param for FindAllMarkers
                   logfc.threshold=0.25,# param for FindAllMarkers
                   total=T,# if T, will perform bulk GSEA, not considering seurat cluster
                   org="mmu",# more information in clusterprofiler vignette
                   present_genes=background,# Background: all genes in dataset (universe in enricher())
                   OrgDb, # mus musculus database need to be loaded prior to use
                   human_ensembl, # # use useMart to load the human database
                   mouse_ensembl,# use useMart to load the mouse database
                   up=T
  ){

    tmp<-object

    markers<-list()

if(length(genelist)==0)
{
    # DE genes for single cluster
    if(total==F){
      for (i in levels(Seurat::Idents(tmp))){
        print(i)
        tmp2<-subset(tmp,idents=i)
        Seurat::Idents(object = tmp2) <- tmp2[[]][,condition]



        markers[[i]] <-Seurat::FindAllMarkers(object = tmp2,only.pos = up,min.pct = min.pct,logfc.threshold = logfc.threshold)
        colnames(markers[[i]])[6]<-"Condition"
        markers[[i]]$cluster<-as.factor(i)

      }
      markers<-do.call(rbind,markers)
      if(up==T){
      markers %>% dplyr::group_by(cluster,Condition) %>% dplyr::slice_max(n = 50, order_by =  avg_log2FC) -> top
      }else{
        markers %>% dplyr::group_by(cluster,Condition) %>% dplyr::slice_min(n = 50, order_by = avg_log2FC) -> top
        }

      entrez <- clusterProfiler::bitr(top$gene, fromType = "SYMBOL", toType="ENTREZID",OrgDb = OrgDb)

      top<-dplyr::left_join(top,entrez,by=c("gene"="SYMBOL"))


    }

    # DE genes bulk
    if(total==T){

      Seurat::Idents(object = tmp) <- tmp[[]][,condition]

      markers <- Seurat::FindAllMarkers(object = tmp, only.pos = T,min.pct = min.pct,logfc.threshold = logfc.threshold)
      if(up==T){
        markers %>% dplyr::group_by(cluster) %>% dplyr::slice_max(n = 50, order_by = avg_log2FC) -> top
      }else{
        markers %>% dplyr::group_by(cluster) %>% dplyr::slice_min(n = 50, order_by = avg_log2FC) -> top
      }

      entrez <- clusterProfiler::bitr(top$gene, fromType = "SYMBOL", toType="ENTREZID",OrgDb = OrgDb)

      top<-dplyr::left_join(top,entrez,by=c("gene"="SYMBOL"))

    }
}else{
  top<-genelist
  entrez <- clusterProfiler::bitr(top$gene, fromType = "SYMBOL", toType="ENTREZID",OrgDb = OrgDb)

  top<-dplyr::left_join(top,entrez,by=c("gene"="SYMBOL"))
}
    results <- list()


    # These databases are for humans ... thus need to find homologues to our mouse genes
    if("Hallmark" %in% GeneSets |
       "DO" %in% GeneSets |
       "cannonicalPathways" %in% GeneSets|
       "ImmunoSignatures" %in% GeneSets |
       "Motifs" %in% GeneSets&org=="mmu"){

      # Get human homologues for mouse genes
      entrez_hsa<- biomaRt::getLDS(attributes = c("entrezgene_id"),
                          filters = "entrezgene_id",
                          values = entrez$ENTREZID,
                          mart = mouse_ensembl,
                          attributesL = c("entrezgene_id"),
                          martL = human_ensembl,
                          uniqueRows=T)

      colnames(entrez_hsa)<-c("entrez_m","entrez_h")
      entrez_hsa$entrez_m<-as.character(entrez_hsa$entrez_m)
      entrez_hsa$entrez_h<-as.character(entrez_hsa$entrez_h)

      top<-dplyr::left_join(top,entrez_hsa,by=c("ENTREZID"="entrez_m"))

    }


    # setting and formatting parameters for GSEA function
    formula_tot<-stats::reformulate(response=ifelse(GeneSets=="GO"|GeneSets=="KEGG","ENTREZID","entrez_h"),termlabels ="cluster")
    formula_cluster<-stats::reformulate(response=ifelse(GeneSets=="GO"|GeneSets=="KEGG","ENTREZID","entrez_h"),termlabels ="cluster+Condition")
    enricher<-ifelse(GeneSets=="GO"|GeneSets=="KEGG"|GeneSets=="DO",paste("clusterProfiler::","enrich",GeneSets,sep=""),
                     "clusterProfiler::enricher")
    if(GeneSets=="DO"){
      enricher<-"DOSE::enrichDO"
    }
    do.MSigDb<-ifelse(enricher=="clusterProfiler::enricher",T,F)
    genes_read<-ifelse(GeneSets=="GO"|GeneSets=="DO",T,F)
    if(GeneSets=="GO"|GeneSets=="KEGG"){
      universe<-present_genes$ENTREZID
    }
    else{
      universe<-present_genes$entrez_h
    }


    # taking different genes, depending on string  in Genesets
    gene_set<-switch(GeneSets,"Hallmark"=hallmark_genes,"cannonicalPathways"=cannonicalPathway_genes,"Motifs"=motifs,"ImmunoSignatures"=immuno_genes)


    # perform GSEA between the levels of defined "Condition"
    # geneClusters=formula.tot only used if total=T
    # geneClusters=formula_cluster used when do.MSigDb
    # readable=T only used when genes_read=T
    #  and so on ...
    results[[GeneSets]]<-do.call(clusterProfiler::compareCluster,c(list(geneClusters=formula_tot)[total],
                                                  list(geneClusters=formula_cluster)[!total],
                                                  list(TERM2GENE=gene_set)[do.MSigDb],
                                                  list(data=top,fun=enricher,universe=universe,
                                                       pAdjustMethod=pCorrection,
                                                       pvalueCutoff=pvalueCutoff,
                                                       qvalueCutoff=qvalueCutoff),
                                                  list(readable=T)[genes_read],
                                                  list(minGSSize=5,
                                                       maxGSSize=500,
                                                       ont=GeneSets)[GeneSets=="DO"],
                                                  list(ont=Gontology,
                                                       OrgDb=OrgDb)[GeneSets=="GO"],
                                                  list(organism=org)[GeneSets=="KEGG"]))



    return(results)
  }
