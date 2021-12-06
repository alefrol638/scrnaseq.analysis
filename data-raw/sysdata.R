####prepare GSEA gene sets ####

###you need to download them from: http://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/7.4/msigdb_v7.4_chip_files_to_download_locally.zip ... the provided files were downloaded in July 2021
###probably you'll need to register there
hallmark_genes <- clusterProfiler::read.gmt(here::here("Metadata/MSigDb/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs","h.all.v7.4.entrez.gmt"))
cannonicalPathway_genes <- clusterProfiler::read.gmt(here::here("Metadata/MSigDb/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs","c2.all.v7.4.entrez.gmt"))
motifs <- clusterProfiler::read.gmt(here::here("Metadata/MSigDb/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs","c3.tft.v7.4.entrez.gmt"))
immuno_genes <- clusterProfiler::read.gmt(here::here("Metadata/MSigDb/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs","c7.all.v7.4.entrez.gmt"))

usethis::use_data(hallmark_genes,cannonicalPathway_genes,motifs,immuno_genes,internal=T)
