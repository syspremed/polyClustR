pamCentroids <- function(l){
  #
  # Calculates PAM centroids of community subtypes.
  #
  
  allAssign <- data.frame(t(read.delim(l$allClust)))
  geneNames <- rownames(l$data)
  geneID <- rownames(l$data)
  batchLabels <- c(rep(1, ncol(l$data)))
  
  centroids <- lapply(1:ncol(allAssign), function(i, x, label){
    cl <- na.omit(as.numeric(x[,i]))
    dataList <- list(batchlabels = batchLabels, x = l$data[,names(na.omit(x[,i]))], y = cl, genenames = geneNames, geneid = geneID)
    dataList <- pamr.knnimpute(dataList)
    dataTrain <- pamr.train(dataList)
    pamTT <- pamr.listgenes(dataTrain, dataList, threshold = 0) %>%
      set_rownames(extract(.,,1)) %>%
      extract(.,,-1)
    class(pamTT) <- 'numeric'
    pamTTuniq <- pamTT[apply(pamTT, 1, function(x){sum(x > 0) == 1}),]
    write.table(pamTTuniq, paste0(l$pamTitle, Sys.Date(), '_', label[i], '_pam.txt'), sep = '\t', quote = FALSE)
  }, x = allAssign, label = names(allAssign))
}
