  #' polyCluster clustering
  #' 
  #' @description Performs clustering via the \code{nmf} and \code{CCP} functions.
  #'  Also calculates gap score and cophenetic corellation coefficient for each consensus matrix,
  #'  and creates pdf plots of these.
  #'
  #' @param l A list of parameters required to carry out the clustering
  #' @param nmfData Character string. File path to the result of a previous NMF clustering by this function, in order to save time.
  #'
  #' @details Not intended for use outside of a call to \code{polyCluster}.
  #' @return Returns an updated list \code{l}.

performClust <- function(l, nmfData = NULL){
  
  setwd(l$initTitle)
  if (!all(l$clusterAlg %in% c("hc", "pm", "km", "nmf"))){
    writeLines(paste(l$clusterAlg[!l$clusterAlg %in% c("hc","pm", "km", "nmf")], "is unsupported. Please choose from \"hc\", \"pm\", \"km\" or \"nmf\"."))
  }
  
  l$consensusMemb <- array(0, c(length(l$clusterAlg), ncol(l$data), l$maxK-1)) # An empty array to store the class membership for all K for all algorithms
  consensusMat <- array(0, c(length(l$clusterAlg), l$maxK-1, ncol(l$data), ncol(l$data))) # An empty array to store the consensus matrices for HC and KM
  
  # Perform the clustering via each algorithm and
  # store the consensus matrices and cophenetic coefficients
  
  if ("hc" %in% l$clusterAlg){
    message("Starting hierarchical clustering...")
    hcccp <- CCP(l$data, maxK = l$maxK, reps = l$reps, clusterAlg = "hc", distance = "euclidean", title=paste(Sys.Date(), "hc", sep="_"))
    l$consensusMemb[which(l$clusterAlg == "hc"),,] <- data.matrix(hcccp[[length(hcccp)]])
    for (i in 2:(length(hcccp)-2)){
      consensusMat[which(l$clusterAlg == "hc"), i-1,,] <- hcccp[[i]]$consensusMatrix
    }
    message("Completed hierarchical clustering.\n")
  }
  
  if ("km" %in% l$clusterAlg){
    message("Starting k-means clustering...")
    kmccp <- CCP(l$data, maxK = l$maxK, reps = l$reps, clusterAlg = "km", distance = "euclidean", title = paste(Sys.Date(), "km", sep="_"))
    l$consensusMemb[which(l$clusterAlg == "km"),,] <- data.matrix(kmccp[[length(kmccp)]])
    for (i in 2:(length(kmccp)-2)){
      consensusMat[which(l$clusterAlg == "km"), i-1,,] <- kmccp[[i]]$consensusMatrix
    }
    message("Completed k-means clustering.\n")
  }
  
  if ("pm" %in% l$clusterAlg){
    message("Starting partitioning around medoids clustering...")
    pmccp <- CCP(l$data, maxK = l$maxK, reps = l$reps, clusterAlg = "pam", distance = "euclidean", title = paste(Sys.Date(), "pm", sep="_"))
    l$consensusMemb[which(l$clusterAlg == "pm"),,] <- data.matrix(pmccp[[length(pmccp)]])
    for (i in 2:(length(pmccp)-2)){
      consensusMat[which(l$clusterAlg == "pm"), i-1,,] <- pmccp[[i]]$consensusMatrix
    }
    message("Completed partitioning around medoids clustering.\n")
  }
  
  if ("nmf" %in% l$clusterAlg){
    message("Starting nonnegative matrix factorization clustering (this may take some time)...")
    
    if(is.null(nmfData)){
      datapos <- nneg(l$data, method = 'pmax')
      nmfccp <- NMF::nmf(datapos, 2:l$maxK, method = 'brunet', nrun = l$reps)
      
      save(nmfccp, file = paste(Sys.Date(), "_nmf_result.rda", sep = ""))
    }
    
    else {
      load(nmfData)
      file.copy(str_replace(nmfData, '[^/\\\\]+\\..+', ''), l$nmfTitle, copy.date = TRUE)
    }
    
    for (i in 1:(length(nmfccp$consensus))){
      consensusMat[which(l$clusterAlg == "nmf"), i,,] <- nmfccp$consensus[[i]]
    }
    
    pdf(paste0(l$nmfTitle, '/', Sys.Date(), l$ref, 'nmf_plots.pdf'), width = 9, height = 6)
    print(plot(nmfccp))
    consensusmap(nmfccp)
    dev.off()
    
    lapply(1:length(nmfccp$fit), function(i){
      write.table(basis(nmfccp$fit[[i]]), paste0(l$nmfBasis, '/', paste(Sys.Date(), 'k', i+1, 'basis.txt', sep = '_')), sep = '\t', quote = FALSE)
      write.table(data.frame(predict(nmfccp$fit[[i]], what = 'features', prob = TRUE)), paste0(l$nmfBasis, '/', paste(Sys.Date(), 'k', i+1, 'metagenes.txt', sep = '_')), sep = '\t', quote = FALSE)
    })
    
    
    l$consensusMemb[which(l$clusterAlg == "nmf"),,] <- sapply(nmfccp$fit, function(x){as.numeric(predict(x, what = 'consensus'))})
    message("Completed nonnegative matrix factorization clustering.\n")
  }
  
  
  ##############################################################################################################
  # Another if() must be added here if another algorithm is required. It must add the class membership to      #
  # l$consensusMemb (samples in rows, k-values in columns) and the consensus matrices to consensusMat.         #
  ##############################################################################################################
  
  ## Cophenetic correlation coefficient for HC, PM and KM
  cophen <- t(apply(consensusMat, c(1, 2), cophcor))
  colnames(cophen) <- l$clusterAlg
  rownames(cophen) <- 2:l$maxK
  
  # Plot cophenetic coefficients
  setwd(l$initTitle)
  pdf(paste(Sys.Date(), "_all_cophenetic.pdf", sep = ""), height = 7, width = 9)
  matplot(rbind(rep(NA, length(l$clusterAlg)), cophen), type = "b", pch = 1, col = l$colPal, lty = 1, xlim = c(2, l$maxK), xlab = "k", ylab = "Cophenetic Correlation Coefficient", main = "Cophenetic Coefficient")
  legend("bottomleft", legend = l$clusterAlg, col=l$colPal, pch=1)
  dev.off()
  
  ## Silhouette width
  sil <- function(clss, cns){
    meanWidth <- mean(silhouette(clss, dmatrix = 1-cns)[,3])
  }
  
  twodClasses <- alply(l$consensusMemb, c(1, 3))
  twodConsensus <- alply(consensusMat, c(1, 2))
  
  twodSw <- mapply(sil, twodClasses, twodConsensus)
  
  silWidth <- sapply(1:length(l$clusterAlg), function(i){
    twodSw[seq(i, length(twodSw), length(l$clusterAlg))]
  })
  
  # Plot silhouette width
  pdf(paste(Sys.Date(), "_all_silhouette.pdf", sep = ""), height = 7, width = 9)
  matplot(rbind(rep(NA, length(l$clusterAlg)), silWidth), type = "b", pch = 1, col = l$colPal, lty = 1, xlim = c(2, l$maxK), xlab = "k", ylab = "Mean Silhouette Width", main = "Silhouette Width")
  legend("bottomleft", legend = l$clusterAlg, col=l$colPal, pch=1)
  dev.off()
  
  
  #####
  ## Gap statistic
  
  gapsArr <- array(0, c((l$maxK-1), 2, length(l$clusterAlg)))
  
  l$clustRes <- list(cophen, l$consensusMemb, gapsArr)
  names(l$clustRes) <- c("cophen", "l$consensusMemb", "gapsArr")
  setwd(l$wd)
  return(l)
}
