nodeEdge <- function(l, recons){
  #
  # Creates a pdf node-edge graph from each reconciliation method
  # and identifies which meta-clusters each sample belongs to. Then uses
  # hyperClust3 to perform a hypergeometric test betwen meta-clusters
  # and known phenotypes.
  #
  # weightMeasures   List. Each element of the list should be a numeric matrix giving the different measures
  #                  by which the thickness of edges should be calculated, e.g. p-values
  # recons           Character vector. The names of the reconcilliation methods to be used (for naming/annotation purposes only)
  # labelledClasses  Character matrix. Algorithms in columns and samples in rows, showing the class membershi[]
  #                  of each sample in each algorithm. Clusters from each algorithm should be easily distinguishable
  #                  ('labelled'), so that e.g. heirarchical cluster 1 and k-means cluster 1 have different names
  # clusterSize      Numeric vector. The number of samples in each cluster.
  # titles           Character vector. The folders each reconciliation should be saved to
  # phenoF           Character vector. The known phenotypes of the samples, with the name of each element being the sample name
  #
  
  weightMeasures <- list(l$adjustP, l$distArr)
  titles <- c(l$hypTitle, l$propTitle)
  
  desat <- function(cols, sat=0.5) {
    X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
    hsv(X[1,], X[2,], X[3,])
  }
  
  rescale <- function(x, newmin, newmax){
    scaled <- (newmax - newmin)/(max(x) - min(x))*(x - min(x)) + newmin
    return(scaled)
  }
  
  for(i in 1:length(recons)){
    
    setwd(titles[i])
    wasError <- TRUE
    
    while(wasError == TRUE){
      
      classPairsM <- data.frame(t(combn(colnames(weightMeasures[[i]]), 2)), weightMeasures[[i]][lower.tri(weightMeasures[[i]])])
      names(classPairsM) <- c("class1", "class2", "weight")
      
      # Remove edges that fall outside cutoff
      classPairsCut <- classPairsM[classPairsM[,3] < 1,]
      classPairsLog <- classPairsCut
      classPairsLog[classPairsLog == 0] <- .Machine$double.eps
      classPairsLog[,3] <- -log10(classPairsLog[,3])
      names(classPairsLog) <- c("class1", "class2", "measure")
      edgeWeight <- classPairsLog[,3]
      
      
      # Plot node edge graph #####
      pars <- par(no.readonly = T)
      namcl <- cbind(colnames(l$adjustP), l$cols)
      NEgraph <- graph.data.frame(classPairsCut[,-3], vertices = namcl, directed = F)
      V(NEgraph)$shape <- "circle"
      V(NEgraph)$color <- V(NEgraph)$frame.color <- namcl[,2]
      V(NEgraph)$size <- rescale(l$numTable, 20, 40)[V(NEgraph)$name]
      V(NEgraph)$label.cex <- (V(NEgraph)$size)/25
      V(NEgraph)$label.family <- "ArialMT"
      V(NEgraph)$label.color <- "black"
      E(NEgraph)$width <- edgeWeight
      E(NEgraph)$width[E(NEgraph)$width >= -log10(0.05)] <- -log10(0.05)
      E(NEgraph)$width <- rescale(E(NEgraph)$width, 5, 15)
      E(NEgraph)$weight <- edgeWeight
      
      if(!is.null(l$phenoFile)){
        subCols <- brewer.pal(length(unique(na.omit(l$pheno))), 'Set2')
        
        names(subCols) <- unique(na.omit(l$pheno))
        
        V(NEgraph)$color <- V(NEgraph)$frame.color <- sapply(1:nrow(l$knownSubsP), function(i){
          if (is.na(l$knownSubsP$pVals[i])){
            V(NEgraph)$color[i] <- 'gray65'
          }
          else if (l$knownSubsP$pVals[i] == '[0,0.05)') {
            V(NEgraph)$color[i] <- subCols[as.character(l$knownSubsP$kSubtypes)[i]]
          }
          else if (l$knownSubsP$pVals[i] == '[0.05,1]') {
            desat(subCols[as.character(l$knownSubsP$kSubtypes)[i]], 0.3)
          }
        })
        
        V(NEgraph)$label.color <- namcl[,2]
      }
      
      NEgraph$layout <- layout_with_fr(NEgraph)
      edgeBreaks <- c(seq(0, -log10(0.05), length.out = 9), -log10(.Machine$double.xmin))
      edgeBins <- cut(edgeWeight, edgeBreaks)
      allEdgeCol <- rev(rainbow(nlevels(edgeBins), start = 0.05, end = 0.6))
      edgeCol <- allEdgeCol[edgeBins]
      E(NEgraph)$color <- edgeCol
      pdf(paste(Sys.Date(), recons[i], "node_edge.pdf", sep = "_"))
      commun <- label.propagation.community(NEgraph)
      plot(NEgraph, mark.groups = commun, mark.shape = 1, mark.border = NA, mark.expand = 30, mark.col = alpha(rep('gray80', length(commun)), 0.7))
      par(fig=c(0, 1, 0, 1/3.1))
      par(new=T)
      par(plt=c(.139,.897,.3,.35))
      key <- matrix(seq(1:length(allEdgeCol)), ncol = 1)
      image(key, col = allEdgeCol, xlim = c(0,1), axes = F, yaxt = "n")
      legend <- as.character(signif((10^-edgeBreaks)[round(seq(1, length(edgeBreaks), length.out = 5), 1)], 1))
      legend[length(legend)] <- '< 0.05'
      mtext(legend, side = 1, at = seq(0, 1, 1/(length(legend)-1)))
      par(pars)
      
      if(!is.null(l$phenofile)){
        legend('topleft', na.omit(names(subCols)), col = subCols, pch = 19, bty = 'n', title = expression(p <= 0.05))
        legend('bottomleft', na.omit(names(subCols)), col = desat(subCols, 0.3), pch = 19, bty = 'n', title = 'p > 0.05')
        legend('topright', paste('Modularity:', signif(modularity(commun), 3)), bty = 'n')
      }
      #####
      dev.off()
      
      if(l$interactive == TRUE){
        tkID <- tkplot(NEgraph)
        inputCond(0, 'Please type 1 to close the window once you have finished arranging the plot: ', 'x', 'x == 1')
        intCoords <- tkplot.getcoords(tkID)
        tkplot.close(tkID)
        dev.off()
        
        NEgraph$layout <- intCoords
        pdf(paste(Sys.Date(), recons[i], "interactive_node_edge.pdf", sep = "_"))
        plot(NEgraph, mark.groups = commun, mark.shape = 1, mark.border = NA, mark.expand = 30, mark.col = alpha(rep('gray80', length(commun)), 0.7))
        par(fig=c(0, 1, 0, 1/3.1))
        par(new=T)
        par(plt=c(.139,.897,.3,.35))
        key <- matrix(seq(1:length(allEdgeCol)), ncol = 1)
        image(key, col = allEdgeCol, xlim = c(0,1), axes = F, yaxt = "n")
        legend <- as.character(signif((10^-edgeBreaks)[round(seq(1, length(edgeBreaks), length.out = 5), 1)], 1))
        legend[length(legend)] <- expression(""<= 0.05)
        mtext(legend, side = 1, at = seq(0, 1, 1/(length(legend)-1)))
        par(pars)
        legend('topleft', na.omit(names(subCols)), col = subCols, pch = 19, bty = 'n', title = expression(p <= 0.05))
        legend('bottomleft', na.omit(names(subCols)), col = desat(subCols, 0.3), pch = 19, bty = 'n', title = 'p > 0.05')
        legend('topright', paste('Modularity:', signif(modularity(commun), 3)), bty = 'n')
        dev.off()
      }
      
      ## Finding samples in node-edge metaclusters
      
      names <- names(l$numTable)
      
      clustList <- communities(commun)
      wasError <- FALSE
      
      names(clustList) <- LETTERS[1:length(clustList)]
      
      NEClust <- data.frame(algClust = names, NEclust = NA)
      for (j in 1:nrow(NEClust)){
        NEClust[j, 2] <- names(clustList[which(lapply(clustList, function(x) names[j] %in% x) == T)])
      }
      write.table(NEClust, paste(Sys.Date(), recons[i], "NE_metaclusters.txt", sep = "_"), sep = "\t", quote = F)
      
      
      NEmap <- data.frame(lapply(l$labelledClasses, function(x){
        y <- NEClust$NEclust[match(x, NEClust$algClust)]
        return(as.factor(y))
      })) %>%
        set_rownames(rownames(l$labelledClasses))
      
      NEall <- apply(NEmap, 1, function(x){
        if (anyDuplicated(x) == 0){
          z <- NA
        }
        else {
          w <- table(x)
          y <- which(w == max(w))
          if (length(y) != 1){
            z <- NA
          }
          else {
            z <- names(which.max(table(x)))
          }
        }
        return(z)
      })
      
      # Get phenotypic data
      l$pheno1 <- NULL
      if (!is.null(l$pheno)){
        l$pheno1 <- as.character(l$pheno[names(NEall)])
        knownSub <- cbind(l$pheno1, NEall)
        
        # Hypergeometric test between known and discovered subtypes
        l[[paste0('knownTable_', recons[i])]] <- table(knownSub[,1], knownSub[,2])
        knownHyper <- hyperClust3(l[[paste0('knownTable_', recons[i])]], paste(Sys.Date(), recons[i], "known", sep = "_"), l)
        
      }
      
      pdf(paste(Sys.Date(), "_", recons[i], "_NE_silhouette.pdf", sep = ""), height = 7, width = 9)
      sk <- silhouette(na.omit(as.numeric(as.factor(NEall))), dist(t(l$data[,!is.na(NEall)])))
      plot(sk, main = "Silhouette width")
      dev.off()
      
      cat(recons[i], as.character(NEall), file = l$allClust, sep = c(rep("\t", length(NEall)), "\n"), append = T)
    }
    setwd(l$wd)
  }
  return(l)
}