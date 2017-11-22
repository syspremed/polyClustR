  #' Read a text file of expression values
  #'
  #' @description Reads a gene expression dataset with genes in rows, samples in columns
  #' in .txt format and converts it into an R matrix
  #' 
  #' @param filename A character vector given the file path of the gene expression data.
  #' 
  #' @details Not intended for use outside of a call to \code{polyCluster}.
  #' @return Returns the expression dataset.
#'
read.txt <- function(filename = "NULL") {
  data_set <- read.delim(filename, header=T, quote="", row.names=1, blank.lines.skip=T, comment.char="", as.is=T, check.names=F)
  data_set <- data.matrix(data_set)
  return(data_set)
}

 #' Read a text file of expression values
 #'
 #' @description Reads a gene expression dataset with genes in rows, samples in columns
 #' in .gct format and converts it into an R matrix
 #' 
 #' @param filename A character vector given the file path of the gene expression data.
  #' 
  #' @details Not intended for use outside of a call to \code{polyCluster}.
  #' @return Returns the expression dataset.
  #'
read.gct <- function(filename = "NULL") {
  data_set <- read.delim(filename, header=T, quote="", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T, check.names=F)
  data_set <- data_set[-1]
  data_set <- data.matrix(data_set)
  return(data_set)
}