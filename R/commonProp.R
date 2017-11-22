  #' Calculating the proportion of maximum intersection (PMI)
  #'  
  #' @description Finds the relative proportion of common samples between
  #'  clusters from different algorithms, i.e. (A ∩ B)/min\{|A|, |B|\}.
  #'  Then calculates the euclidean and cosine distance metrics for the propotions.
  #'  Creates a pdf heatmap to show distances.
  #'
  #' @param l The output of a call to \code{performClust}.
  #' @details Not intended for use outside of a call to \code{polyCluster}.
  #' @return Returns an updated list \code{l}.
  #'

commonProp <- function(l){
  
  # Divide the number of shared samples by the number of samples in the smallest cluster
  clusterMin <- l$clusterTable
  clusterMin[lower.tri(clusterMin)] <- combn(l$numTable, 2, min)
  clusterMin <- as.matrix(forceSymmetric(clusterMin, "L"))
  clusterProp <- l$clusterTable/clusterMin
  diag(clusterProp) <- NA
  write.table(clusterProp, file = paste(Sys.Date(), "shared_proportion.txt", sep = "_"), quote = F, sep = "\t", na = "1")
  
  # Heatmap of proportion
  pdf(paste(l$propTitle, Sys.Date(), "_min_cluster_prop.pdf", sep=""), height=13, width=13)
  image((clusterProp), xaxt = "n", yaxt = "n", col = rev(l$colRain))
  mtext(rownames(clusterProp), side = 1, at=seq(0, 1, (1/(nrow(clusterProp)-1))), col = l$cols)
  mtext(rownames(clusterProp), side = 2, at=seq(0, 1, (1/(nrow(clusterProp)-1))), col = l$cols)
  heatmap.2(clusterProp, trace = "none", density.info = "density", denscol = "black", scale = "none", ColSideColors = l$cols, RowSideColors = l$cols, col = rev(l$colRain))
  dev.off()
  
  # Create and heatmap distance matrices of the proportion of common samples using cosine and euclidean distance
  diag(clusterProp) <- 1
  dissimCos <- cosine(clusterProp)
  dissimCos <- 1 - as.matrix(dissimCos)
  diag(dissimCos) <- NA
  
  # Plots
  pdf(paste(l$propTitle, Sys.Date(), "_cluster_cosine.pdf", sep=""), height=13, width=13)
  image(dissimCos, xaxt = "n", yaxt = "n", col = l$colRain)
  mtext(rownames(clusterProp), side = 1, at = seq(0, 1, (1/(nrow(clusterProp)-1))), col = l$cols)
  mtext(rownames(clusterProp), side = 2, at = seq(0, 1, (1/(nrow(clusterProp)-1))), col = l$cols)
  heatmap.2(dissimCos, trace = "none", density.info = "density", denscol = "black", scale = "none", ColSideColors = l$cols, RowSideColors = l$cols, col = l$colRain)
  dev.off()
  write.table(dissimCos, file = paste(Sys.Date(), "cosine_similarity.txt", sep = "_"), quote = F, sep = "\t", na = "1")
  
  diag(clusterProp) <- 1
  dissimEuc <- dist(clusterProp, method = "euclidean")
  dissimEuc <- as.matrix(dissimEuc)
  diag(dissimEuc) <- NA
  
  pdf(paste(l$propTitle, Sys.Date(), "_cluster_euclidean.pdf", sep=""), height=13, width=13)
  image(dissimEuc, xaxt = "n", yaxt = "n", col = l$colRain)
  mtext(rownames(clusterProp), side = 1, at = seq(0, 1, (1/(nrow(clusterProp)-1))), col = l$cols, cex.axis = 0.5)
  mtext(rownames(clusterProp), side = 2, at = seq(0, 1, (1/(nrow(clusterProp)-1))), col = l$cols, cex.axis = 0.5)
  heatmap.2(dissimEuc, trace = "none", density.info = "density", denscol = "black", scale = "none", ColSideColors = l$cols, RowSideColors = l$cols, col = l$colRain)
  dev.off()
  write.table(dissimEuc, file = paste(Sys.Date(), "euclidean_distance.txt", sep = "_"), quote = F, sep = "\t", na = "0")
  
  l$distArr <- 1-clusterProp
  return(l)
}