read.txt <- function(filename = "NULL") {
  #
  # Reads a gene expression dataset with genes in rows, samples in columns
  # in .txt format and converts it into an R matrix
  #
  data_set <- read.delim(filename, header=T, quote="", row.names=1, blank.lines.skip=T, comment.char="", as.is=T, check.names=F)
  data_set <- data.matrix(data_set)
  return(data_set)
}

read.gct <- function(filename = "NULL") {
  #
  # Reads a gene expression dataset with genes in rows, samples in columns
  # in .gct format (which has extra rows and columns without useful data) and converts it into an R matrix
  #
  data_set <- read.delim(filename, header=T, quote="", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T, check.names=F)
  data_set <- data_set[-1]
  data_set <- data.matrix(data_set)
  return(data_set)
}