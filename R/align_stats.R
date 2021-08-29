#' @title Alignment Statistics
#' @description  finds QC statistics and summarises the information of multiple samples in one plot,function taken from scripts of Jonas Schrepping/DZNE
#' @param path character, Path to DropSeq results folder , ends with .../output/results
#' @param subset character vector, Cells to subset for
#' @param exclude Boolean, if TRUE cells specified in subset will be excluded
#' @param legend.pos the position of legends ("none", "left", "right", "bottom", "top", or two-element numeric vector)
#' @param ymax integer, maximal number of counts shown on y-axis
#' @param title Main title
#' @examples
#' align_stats(path="/data/SeqWell/190613/alignment/2020-07-17/output/results",subset=c("CellID1","CellID2"))
#' @export align_stats
align_stats <- function(path,subset=F,exclude=F,legend.pos="right",ymax=200000000,title="Read Statistics") {
  fastQC_reads <- read.delim(paste(path,"/reports/fastqc_reads_data/multiqc_fastqc.txt",sep=""), stringsAsFactors = T)
  fastQC_barcodes <- read.delim(paste(path,"/reports/fastqc_barcodes_data/multiqc_fastqc.txt",sep=""), stringsAsFactors = T)
  RNAfilt_QC <- read.delim(paste(path,"/reports/RNA_filtering_data/multiqc_cutadapt.txt",sep=""), stringsAsFactors = T)
  starQC <- read.delim(paste(path,"/reports/star_data/multiqc_star.txt",sep=""), stringsAsFactors = T)
  QC <- data.frame(total_reads = fastQC_reads$Total.Sequences,
                   duplicate_barcodes = fastQC_barcodes$Total.Sequences - fastQC_barcodes$Total.Sequences*fastQC_barcodes$total_deduplicated_percentage/100,
                   deduplicated_barcodes = fastQC_barcodes$Total.Sequences*fastQC_barcodes$total_deduplicated_percentage/100,
                   duplicate_reads = fastQC_reads$Total.Sequences - fastQC_reads$Total.Sequences*fastQC_reads$total_deduplicated_percentage/100,
                   deduplicated_reads = fastQC_reads$Total.Sequences*fastQC_reads$total_deduplicated_percentage/100,
                   survived_filtering = starQC$total_reads,
                   aligned_reads = starQC$uniquely_mapped,
                   ID = starQC$Sample)
  data.frame(ID = starQC$Sample,
             total_reads = fastQC_reads$Total.Sequences,
             percent_unique_barcodes = round(fastQC_barcodes$total_deduplicated_percentage,2),
             percent_unique_reads = round(fastQC_reads$total_deduplicated_percentage,2))
  # QC Plot
  QC$ID<-sub("Mouse_Brain_Microglia_","",QC$ID,fixed=T)
  QC$ID<-sub("_S.*$","",QC$ID,fixed=F)
  QC$ID<-sub("_pool_","_",QC$ID,fixed=T)
  if(subset!=F){
    QC<-QC[grep(subset,QC$ID,fixed=T,invert=exclude),]
  }
  ggplot2::ggplot(reshape2::melt(QC), ggplot2::aes(x = ID, y = value, fill = variable)) +
    ggplot2::geom_bar(width=.7,stat = "identity", position = "dodge", colour = "black") +
    ggplot2::ylab("Reads") +
    ggplot2::xlab(NULL) +
    ggplot2::scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE),
                       breaks = c(c(10, 50, 100, 150) * 10^6, -c(10, 50, 100) * 10^6)) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white", colour = "black"),
          legend.position = legend.pos, axis.text.x = element_text(angle = 60, hjust=1,size=5),
          axis.text.y = element_text( vjust = 0.5, hjust=1,size=5),
          axis.title = element_text(size=7),
          legend.title = element_text(size=7),
          legend.key.size = unit(1,"mm"),
          legend.text = element_text(size=5)
    )+ggplot2::ylim(0,ymax)
}

