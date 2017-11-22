  #' Hypergeometric test of cluster overlap
  #'
  #' @description Compares all clusters found by different algorithms via a
  #'  hypergeometric test and creates a pdf heatmap to show
  #'  the resulting p values.
  #'  where n is the total number of clusters found across the algorithms.
  #'  
  #' @param l The output of a call to \code{performClust}.
  #' 
  #' @details Not intended for use outside of a call to \code{polyCluster}.
  #' @return Returns an updated list \code{l}.
  #'
  
hyperClust2 <- function(l){
  
  # Do a hypergeometric test between each pair of clusters
  raw.p <- matrix(NA, nrow(l$clusterTable), ncol(l$clusterTable))
  for (i in 1:nrow(l$clusterTable)){
    for(j in 1:ncol(l$clusterTable)){
      raw.p[i, j] <- phyper((l$clusterTable[i, j])-1, (sum(l$clusterTable[,j]))/length(l$clusterAlg), nrow(l$labelledClasses)-((sum(l$clusterTable[,j]))/length(l$clusterAlg)), (sum(l$clusterTable[i,]))/length(l$clusterAlg), lower.tail=FALSE)
    }
  }
  
  # Adjust p-values
  l$adjust.p <- matrix(p.adjust(c(raw.p), "fdr"), nrow(l$clusterTable), ncol(l$clusterTable))
  rownames(l$adjust.p) <- rownames(raw.p) <- rownames(l$clusterTable)
  colnames(l$adjust.p) <- colnames(raw.p) <- colnames(l$clusterTable)
  
  write.table(l$adjust.p, file = paste(Sys.Date(), "cluster", "adjustp.txt", sep = "_"), quote = F, sep = "\t")
  
  adjust.p2 <- l$adjust.p
  diag(adjust.p2) <- NA
  
  # Plots of p-values
  
  pdf(paste(Sys.Date(), "cluster_hypergeometric.pdf", sep="_"), height=13, width=13)
  image(adjust.p2, xaxt = "n", yaxt = "n", col = l$colRain, breaks = c(0, seq(0.05, 1, length.out = length(l$colRain))))
  mtext(rownames(adjust.p2), side = 1, at = seq(0, 1, (1/(nrow(adjust.p2)-1))), col = l$cols, cex.axis = 0.5)
  mtext(rownames(adjust.p2), side = 2, at = seq(0, 1, (1/(nrow(adjust.p2)-1))), col = l$cols, cex.axis = 0.5)
  heatmap.2(adjust.p2, trace = "none", density.info = "density", denscol = "black", scale = "none", ColSideColors = l$cols[1:ncol(adjust.p2)], RowSideColors = l$cols[1:nrow(adjust.p2)], col = l$colRain, breaks = c(0, seq(0.05, 1, length.out = length(l$colRain))))
  dev.off()
  
  return(l)
}