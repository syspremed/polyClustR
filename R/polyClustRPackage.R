polyCluster <- function(filename, clusterAlg = c("hc", "pm", "km", "nmf"),
                        maxK = 7, reps = 100, phenoFile = NULL, ref = "",
                        nmfData = NULL, interactive = FALSE){
  #
  #'   @param filename     Character string. File name of preprocessed and normalized expression values with genes in rows and samples in columns
  #'   @param clusterAlg   Character vector. Any combination of "hc" (hierarchical clustering), "pm" (partitioning around medoids), "km" (k-means) or "nmf" (nonnegative matrix factorization)
  #'   @param maxK         Integer > 2. The maximum number of clusters of samples to evaluate.
  #'   @param reps         Integer. The number of resampling (for "hc", "pm" and "km") or starting seed ("nmf") iterations.
  #'   @param phenoFile    Character string. File name of known phenotypes, with sample name in first column and phenotype in second column.
  #'   @param ref          Character string. Reference with which to name the output of the analysis.
  #'   @param nmfData      Character string. File path to the result of a previous NMF clustering by this function, in order to save time.
  #'   @param interactive  Logical. If FALSE, networks will be plotted with an automatic layout. If TRUE, a users can customise the layout through tkplot.
  #'   @return Runs the entire polyCluster pipeline and provides output in the current working directory.
  #'   @examples polyCluster('path/to/expression/file', clusterAlg = c('hc', 'km'), ref = 'test_run')
  #
  
  l <- list()

  l$colPal <- c("darkorange1", "mediumvioletred", "seagreen", "powderblue", "rosybrown1")
  l$colRain <- c('red', rainbow(999, start = 0.05, end = 0.7))
  l$clusterAlg <- clusterAlg
  l$maxK <- maxK
  l$phenoFile <- phenoFile
  l$reps <- reps
  l$ref <- ref
  l$interactive <- interactive

  ##
  installed <- installed.packages()[,1]
  if (!'NMF' %in% installed){
    stop('NMF package not detected. Please install NMF package before running this pipeline.')
  }
  library(NMF)

  # Read data
  if (grepl(".txt$", filename)==T){
    l$data <- read.txt(filename)
  }

  else if (grepl(".gct$", filename)==T){
    l$data <- read.gct(filename)
  }

  else if (grepl(".txt$", filename)==F & grepl(".gct$", filename)==F) {
    stop ("Expression file must be in .txt or .gct format")
  }
  ##


  # Create subfolders to store output
  l$wd <- getwd()
  options(warn = -1)

  analysisFolder <- paste(Sys.Date(), ref, "analysis", sep = "_")
  l$analysisTitle <- paste(l$wd, "/Output/", analysisFolder, "/", sep = "")
  dir.create(l$analysisTitle, recursive = TRUE)

  initFolder <- paste(Sys.Date(), ref, "initial", sep = "_")
  hypFolder <- paste(Sys.Date(), ref, "hypergeometric", sep = "_")
  propFolder <- paste(Sys.Date(), ref, "proportion", sep = "_")

  pamFolder <- paste(Sys.Date(), ref, "pam", sep = "_")

  l$initTitle <- paste(l$analysisTitle, initFolder, "/", sep = "")
  l$hypTitle <- paste(l$analysisTitle, hypFolder, "/", sep = "")
  l$propTitle <- paste(l$analysisTitle, propFolder, "/", sep = "")

  l$pamTitle <- paste(l$analysisTitle, pamFolder, "/", sep = "")
  l$knownTitle <- paste0(l$analysisTitle, initFolder, "/", Sys.Date(), "_known_hypergeometric", "/")

  dir.create(l$hypTitle); dir.create(l$propTitle); dir.create(l$initTitle); dir.create(l$knownTitle); dir.create(l$pamTitle)

  
  if ("nmf" %in% clusterAlg){
    nmfFolder <- paste(Sys.Date(), "nmf", sep = "_")
    l$nmfTitle <- paste(l$initTitle, "/", nmfFolder, "/", sep = "")
    dir.create(l$nmfTitle)
    l$nmfBasis <- paste(l$nmfTitle, paste(Sys.Date(), l$ref, 'nmf_basis', sep = '_'), sep = '/')
    dir.create(l$nmfBasis)
  }
  
  options(warn = 0)
  
  ##

  # Create txt file for all class assignments
  l$allClust <- paste(l$analysisTitle, Sys.Date(), "_all_assignments.txt", sep = "")
  cat("samples", colnames(l$data), file = l$allClust, sep = c(rep("\t", ncol(l$data)), "\n"), append = F)
  ##

  # Cluster data
  l <- performClust(l, nmfData)
  save(list = ls(), file = paste0(l$analysisTitle, Sys.Date(), "_session_data.Rdata"))
  ##

  # Reconcile clusters
  l <- testClust(l)
  save(list = ls(), file = paste0(l$analysisTitle, Sys.Date(), "_session_data.Rdata"))
  ##

  # Get PAM centroids
  pamCentroids(l)
  save(list = ls(), file = paste0(l$analysisTitle, Sys.Date(), "_session_data.Rdata"))

  setwd(l$wd)

  writeLines(capture.output(sessionInfo()), paste(l$analysisTitle, Sys.Date(), "_", l$ref, "_session_info.txt", sep = ""))
  writeLines('Success!')
}

