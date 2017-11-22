hyperClust3 <- function(clusterTable, outfile, l){
  #
  # Compares clusters discovered from reconciliation methods with
  # known phenotypes via a hypergeometric test and creates a
  # pdf heatmap to show the resulting p values.
  # Returns a matrix of the resulting adjusted p-values.
  #
  # clusterTable  Numeric matrix. A matrix giving the number of samples shared between all clusters and all phenotypes
  # outfile       Character string. Prepended to the p value and heatmap filenames
  #
  
  # Do a hypergeometric test between each pair of clusters
  raw.p <- matrix(NA, nrow(clusterTable), ncol(clusterTable))
  for (i in 1:nrow(clusterTable)){
    for(j in 1:ncol(clusterTable)){
      raw.p[i, j] <- phyper((clusterTable[i, j])-1, (sum(clusterTable[,j])), nrow(l$labelledClasses)-(sum(clusterTable[,j])), (sum(clusterTable[i,])), lower.tail=FALSE)
    }
  }
  
  # Adjust p-values
  adjust.p <- matrix(p.adjust(c(raw.p), "fdr"), nrow(clusterTable), ncol(clusterTable))
  rownames(adjust.p) <- rownames(raw.p) <- rownames(clusterTable)
  colnames(adjust.p) <- colnames(raw.p) <- colnames(clusterTable)
  
  write.table(adjust.p, file = paste(outfile, "adjustp.txt", sep = "_"), quote = F, sep = "\t")
  
  # Plots of p-values
  plotfile <- unlist(strsplit(outfile,".txt"))[1]
  
  pdf(paste(plotfile, ".pdf", sep=""), height=13, width=13)
  image(adjust.p, xaxt = "n", yaxt = "n", col = l$colRain, breaks = c(0, seq(0.05, 1, length.out = length(l$colRain))))
  mtext(rownames(adjust.p), side = 1, at = seq(0, 1, (1/(nrow(adjust.p)-1))), col = "black", cex.axis = 0.5)
  mtext(colnames(adjust.p), side = 2, at = seq(0, 1, (1/(ncol(adjust.p)-1))), col = "black", cex.axis = 0.5)
  mtext("known", side = 1, outer = T)
  mtext("nmf", side = 2, outer = T)
  
  tryCatch(heatmap.2(adjust.p, trace = "none", density.info = "density",
                     denscol = "black", scale = "none", col = l$colRain,
                     breaks = c(0, seq(0.05, 1, length.out = length(l$colRain)))), error = function(e) e)
  
  dev.off()
  
  return(adjust.p)
}
