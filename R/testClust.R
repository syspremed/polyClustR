  #' Test clustering result.
  #'
  #' @description Compares clusters found by clustering algorithms via the proportion
  #'  of common samples and a hypergeometric test.
  #'
  #' @param l List. The output of a call to performClust
  #' 
  #' @details Not intended for use outside of a call to \code{polyCluster}.
  #' @return Returns an updated list \code{l}.

testClust <- function(l){
  
  sampleClasses <- array(0, c(ncol(l$data), length(l$clusterAlg))) # An array of the best fitting class assignments as found by each algorithm
  rownames(sampleClasses) <- colnames(l$data)
  colnames(sampleClasses) <- l$clusterAlg
  
  # User input for optimal k
  i <- 1
  while (i <= length(l$clusterAlg)){
    
    optK <- inputCond(1, paste("Enter optimal k for ", l$clusterAlg[i], ": ", sep = ""), "x", paste0("x > 1 & x <= ", l$maxK, " & !is.na(x)"))
    
    sampleClasses[,i] <- l$consensusMemb[i,, optK - 1]
    i <- i+1
  }
  
  l$labelledClasses <- data.frame(sampleClasses)
  
  write.table(t(l$labelledClasses), file = l$allClust, quote = FALSE, sep = '\t')
  
  for (i in 1:length(l$clusterAlg)){                                            # Add prefixes to class assignments from each algorithm to
    l$labelledClasses[,i] <- factor(paste(l$clusterAlg[i], l$labelledClasses[,i], sep=""))  # distinguish, e.g., hierarchical cluster 1 from k-means cluster 1
  }
  
  for (i in 1:length(l$clusterAlg)){
    for (j in 1:length(l$clusterAlg)){
      tab <- table(l$labelledClasses[,i], l$labelledClasses[,j])
      if (j == 1){
        clusterRow <- tab
      }
      else {
        clusterRow <- cbind(clusterRow, tab)
      }
    }
    if(i == 1){
      l$clusterTable <- clusterRow
    }
    else{
      l$clusterTable <- rbind(l$clusterTable, clusterRow)
    }
    rm(clusterRow)
  }
  
  # Set colors for plots
  dimension <- sum(sapply(l$labelledClasses, nlevels))
  l$cols <- c(rep(NA, dimension))
  groups <- c(0, cumsum(sapply(l$labelledClasses, nlevels)))
  for (i in 1:length(l$clusterAlg)){
    l$cols[seq(groups[i]+1, groups[i+1])] <- l$colPal[i]
  }
  
  # Get phenotypic data
  
  if (!is.null(l$phenoFile)){
    l$pheno <- as.matrix(read.table(l$phenoFile, header = FALSE, row.names = 1, skip = 1))[, 1]
    if(all(substr(rownames(l$labelledClasses), 1, 1) == 'X') & !all(substr(rownames(l$pheno), 1, 1) == 'X')){
       rownames(l$pheno) <- paste0('X', rownames(l$pheno) 
    }
    l$pheno <- l$pheno[rownames(l$labelledClasses)]
    setwd(l$knownTitle)
    allTabs <- sapply(l$labelledClasses, function(x){table(x, l$pheno)}, simplify = FALSE)
    allKnownHyper <- sapply(1:length(allTabs), function(i){hyperClust3(allTabs[[i]], paste(Sys.Date(), names(allTabs[i]), 'hyper_known', sep = '_'), l = l)}, simplify = FALSE) %>%
      do.call(rbind, .)
    
    kSubtypes <- colnames(allKnownHyper)[apply(allKnownHyper, 1, function(x){
      minP <- min(x)
      if (sum(x == minP) == 1){
        return(which.min(x))
      }
      else {
        return(NA)
      }
    })]
    
    pBreaks <- cut(0:1, c(0, 0.05, 1), include.lowest = TRUE, right = FALSE)
    
    pVals <- apply(allKnownHyper, 1, function(x){
      minP <- min(x)
      if (sum(x == minP) == 1 & minP <= 0.05){
        return(as.character(pBreaks[1]))
      }
      else if (sum(x == minP) == 1 & minP > 0.05){
        return(as.character(pBreaks[2]))
      }
      else {
        return(NA)
      }
    })
    
    l$knownSubsP <- data.frame(kSubtypes, pVals)
    
    setwd(l$wd)
  }
  
  else {
    l$pheno <- NULL
  }
  
  l$numTable <- table(unlist(l$labelledClasses))
  
  # Perform the reconciliations
  
  setwd(l$propTitle)
  l$distArr <- commonProp(l)$distArr
  setwd(l$wd)
  
  setwd(l$hypTitle)
  l$adjustP <- hyperClust2(l)$adjust.p
  setwd(l$wd)
  
  nodeEdge(l, c("hyper", "minprop"))
  
  setwd(l$wd)
  return(l)
}
