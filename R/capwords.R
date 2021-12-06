#' @title Capitalise Words
#' @description required for Seurat::AddModuleScore(), because often the gene names from other sources are not consisted in case.
#' @param s Character vector with words to be capitalised
#' @param strict Boolean, TRUE recommended, since it makes sure that only the first letter is capital and other low case.
#' @export
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}
