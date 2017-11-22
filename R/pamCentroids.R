  #' PAM centroids
  #'
  #' @description Calculates PAM centroids of community subtypes.
  #'
  #' @param l The output of a call to \code{testClust}.
  #' 
  #' @details Not intended for use outside of a call to \code{polyCluster}.
  #' @return Returns a list of PAM centroids for all communities.
  #'

pamCentroids <- function(l){
  
  allAssign <- data.frame(t(read.delim(l$allClust)))
  geneNames <- rownames(l$data)
  geneID <- rownames(l$data)
  batchLabels <- c(rep(1, ncol(l$data)))
  
  centroids <- lapply(1:ncol(allAssign), function(i, x, label){
    cl <- na.omit(as.numeric(x[,i]))
    dataList <- list(batchlabels = batchLabels, x = l$data[,names(na.omit(x[,i]))], y = cl, genenames = geneNames, geneid = geneID)
    capture.output(dataList <- pamr.knnimpute(dataList))
    capture.output(dataTrain <- pamr.train(dataList))
    capture.output(pamTT <- pamr.listgenes(dataTrain, dataList, threshold = 0) %>%
      set_rownames(extract(.,,1)) %>%
      extract(.,,-1))
    class(pamTT) <- 'numeric'
    pamTTuniq <- pamTT[apply(pamTT, 1, function(x){sum(x > 0) == 1}),]
    write.table(pamTTuniq, paste0(l$pamTitle, Sys.Date(), '_', label[i], '_pam.txt'), sep = '\t', quote = FALSE)
  }, x = allAssign, label = names(allAssign))
}
