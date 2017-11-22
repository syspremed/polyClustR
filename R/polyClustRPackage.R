polyCluster <- function(filename, clusterAlg = c("hc", "pm", "km", "nmf"),
                        maxK = 7, reps = 100, phenoFile = NULL, ref = "",
                        nmfData = NULL, interactive = FALSE){
  #
  #   filename     Character string. File name of preprocessed and normalized expression values with genes in rows and samples in columns
  #   clusterAlg   Character vector. Any combination of "hc" (hierarchical clustering), "pm" (partitioning around medoids), "km" (k-means) or "nmf" (nonnegative matrix factorization)
  #   maxK         Integer > 2. The maximum number of clusters of samples to evaluate.
  #   reps         Integer. The number of resampling (for "hc", "pm" and "km") or starting seed ("nmf") iterations.
  #   phenoFile    Character string. File name of known phenotypes, with sample name in first column and phenotype in second column.
  #   ref          Character string. Reference with which to name the output of the analysis.
  #   nmfData      Character string. File path to the result of a previous NMF clustering by this function, in order to save time.
  #   interactive  Logical. If FALSE, networks will be plotted with an automatic layout. If TRUE, a users can customise the layout through tkplot.
  #


CCP <- function (d = NULL, maxK = maxK, reps = reps, pItem = 0.8, pFeature = 1,
                 clusterAlg = "hc", title = "untitled_consensus_cluster",
                 innerLinkage = "average", finalLinkage = "average", distance = "euclidean",
                 ml = NULL, tmyPal = NULL, seed = NULL, plot = "pdf", writeTable = TRUE,
                 weightsItem = NULL, weightsFeature = NULL, verbose = F, corUse = "everything") {
  #
  # Takes an expression matrix and finds clusters via hierarchical clustering, k-means or partitioning around medoids.
  # Adapted from ConsensusClusterPlus of the ConsensusClusterPlus package.
  # Returns a list including a lot of information for each k, most importantly the consensus matrices and class assignments
  #
  # From the ConsensusClusterPlus help page:
  # d              Numeric matrix. Data to be clustered, where columns=items/samples and rows are features
  # maxK           Integer. Maximum number of clusters to evaluate
  # reps           Integer. Number of subsamples so consensus can be evaluated
  # pItem          Numerical value. Proportion of items (columns) to sample in each subsampling
  # pFeature       Numerical value. Proportion of features (rows) to sample in each subsampling
  # clusterAlg     Character string. Cluster algorithm: 'hc' heirarchical (hclust), 'pam' for paritioning around medoids, 'km' for k-means
  # title          Character string. Name for output directory. Directory is created only if plot is not NULL or writeTable is TRUE. This title can be an abosulte or relative path.
  # innerLinkage   heirarchical linkage method for subsampling
  # finalLinkage   heirarchical linkage method for consensus matrix
  # distance       Character string. 'pearson': (1 - Pearson correlation), 'spearman' (1 - Spearman correlation), 'euclidean', 'binary', 'maximum', 'canberra', 'minkowski" or custom distance function.
  # ml             Optional. Prior result. If supplied then only do graphics and tables
  # tmyPal         Optional. Character vector. Colors for consensus matrix
  # seed           Optional. Numerical. Sets random seed for reproducible results
  # plot           Character string. NULL - print to screen, 'pdf', 'png', 'pngBMP' for bitmap png, helpful for large datasets
  # writeTable     Logical. TRUE - write ouput and log to csv
  # weightsItem    Optional. Numerical vector. Weights to be used for sampling items
  # weightsFeature Optional. Numerical vector. Weights to be used for sampling features
  # verbose        Logical. If TRUE, print messages to the screen to indicate progress. This is useful for large datasets
  # corUse         Optional. Character string. Specifies how to handle missing data in correlation distances 'everything','pairwise.complete.obs', 'complete.obs'
  #

  ## CCP Nested Functions

  ccRun <- function (d = d, maxK = NULL, repCount = NULL, diss = inherits(d, "dist"),
                     pItem = NULL, pFeature = NULL, innerLinkage = NULL,
                     distance = NULL, clusterAlg = NULL, weightsItem = NULL, weightsFeature = NULL,
                     verbose = NULL, corUse = NULL){
    m = vector(mode = "list", repCount)
    ml = vector(mode = "list", maxK)
    n <- ifelse(diss, ncol(as.matrix(d)), ncol(d))
    mCount = mConsist = matrix(c(0), ncol = n, nrow = n)
    ml[[1]] = c(0)
    if (is.null(distance))
      distance <- "euclidean"
    acceptable.distance <- c("euclidean", "maximum", "manhattan",
                             "canberra", "binary", "minkowski", "pearson", "spearman")
    main.dist.obj <- NULL
    if (diss) {
      main.dist.obj <- d
      if ((!is.null(pFeature)) && (pFeature < 1)) {
        message("user-supplied data is a distance matrix; ignoring user-specified pFeature parameter\n")
        pFeature <- 1
      }
      if (!is.null(weightsFeature)) {
        message("user-supplied data is a distance matrix; ignoring user-specified weightsFeature parameter\n")
        weightsFeature <- NULL
      }
    }
    else {
      if ((clusterAlg != "km") && (is.null(pFeature) || ((pFeature ==
                                                            1) && is.null(weightsFeature)))) {
        if (inherits(distance, "character")) {
          if (!distance %in% acceptable.distance & (class(try(get(distance),
                                                              silent = T)) != "function"))
            stop("unsupported distance.")
          if (distance == "pearson" | distance == "spearman") {
            main.dist.obj <- as.dist(1 - cor(d, method = distance,
                                             use = corUse))
          }
          else if (class(try(get(distance), silent = T)) ==
                     "function") {
            main.dist.obj <- get(distance)(t(d))
          }
          else {
            main.dist.obj <- dist(t(d), method = distance)
          }
          attr(main.dist.obj, "method") <- distance
        }
        else stop("unsupported distance specified.")
      }
      else {
      }
    }
    for (i in 1:repCount) {
      if (verbose) {
        message(paste("random subsample", i))
      }
      sample_x = sampleCols(d, pItem, pFeature, weightsItem,
                            weightsFeature)
      this_dist = NA
      if (!is.null(main.dist.obj)) {
        boot.cols <- sample_x$subcols
        this_dist <- as.matrix(main.dist.obj)[boot.cols,
                                              boot.cols]
        if (clusterAlg != "km") {
          this_dist <- as.dist(this_dist)
          attr(this_dist, "method") <- attr(main.dist.obj,
                                            "method")
        }
      }
      else {
        if (clusterAlg != "km") {
          if (!distance %in% acceptable.distance & (class(try(get(distance),
                                                              silent = T)) != "function"))
            stop("unsupported distance.")
          if ((class(try(get(distance), silent = T)) ==
                 "function")) {
            this_dist <- get(distance)(t(sample_x$submat))
          }
          else {
            if (distance == "pearson" | distance == "spearman") {
              this_dist <- as.dist(1 - cor(sample_x$submat,
                                           use = corUse, method = distance))
            }
            else {
              this_dist <- dist(t(sample_x$submat), method = distance)
            }
          }
          attr(this_dist, "method") <- distance
        }
        else {
          if (is.null(pFeature) || ((pFeature == 1) &&
                                      is.null(weightsFeature))) {
            this_dist <- d[, sample_x$subcols]
          }
          else {
            if (is.na(sample_x$submat)) {
              stop("error submat is NA")
            }
            this_dist <- sample_x$submat
          }
        }
      }
      this_cluster = NA
      if (clusterAlg == "hc") {
        this_cluster = hclust(this_dist, method = innerLinkage)
      }
      mCount <- connectivityMatrix(rep(1, length(sample_x[[3]])),
                                   mCount, sample_x[[3]])
      for (k in 2:maxK) {
        if (verbose) {
          message(paste("  k =", k))
        }
        if (i == 1) {
          ml[[k]] = mConsist
        }
        this_assignment = NA
        if (clusterAlg == "hc") {
          this_assignment = cutree(this_cluster, k)
        }
        else if (clusterAlg == "kmdist") {
          this_assignment = kmeans(this_dist, k, iter.max = 10,
                                   nstart = 1, algorithm = c("Hartigan-Wong"))$cluster
        }
        else if (clusterAlg == "km") {
          this_assignment <- kmeans(t(this_dist), k, iter.max = 10,
                                    nstart = 1, algorithm = c("Hartigan-Wong"))$cluster
        }
        else if (clusterAlg == "pam") {
          this_assignment <- pam(x = this_dist, k, diss = TRUE,
                                 metric = distance, cluster.only = TRUE)
        }
        else {
          this_assignment <- get(clusterAlg)(this_dist,
                                             k)
        }
        ml[[k]] <- connectivityMatrix(this_assignment, ml[[k]],
                                      sample_x[[3]])
      }
    }
    res = vector(mode = "list", maxK)
    for (k in 2:maxK) {
      tmp = triangle(ml[[k]], mode = 3)
      tmpCount = triangle(mCount, mode = 3)
      res[[k]] = tmp/tmpCount
      res[[k]][which(tmpCount == 0)] = 0
    }
    return(res)
  }

  sampleCols <- function (d, pSamp = NULL, pRow = NULL, weightsItem = NULL, weightsFeature = NULL){
    space <- ifelse(inherits(d, "dist"), ncol(as.matrix(d)),
                    ncol(d))
    sampleN <- floor(space * pSamp)
    sampCols <- sort(sample(space, sampleN, replace = FALSE,
                            prob = weightsItem))
    this_sample <- sampRows <- NA
    if (inherits(d, "matrix")) {
      if ((!is.null(pRow)) && ((pRow < 1) || (!is.null(weightsFeature)))) {
        space = nrow(d)
        sampleN = floor(space * pRow)
        sampRows = sort(sample(space, sampleN, replace = FALSE,
                               prob = weightsFeature))
        this_sample <- d[sampRows, sampCols]
        dimnames(this_sample) <- NULL
      }
      else {
      }
    }
    return(list(submat = this_sample, subrows = sampRows, subcols = sampCols))
  }

  connectivityMatrix <- function (clusterAssignments, m, sampleKey){
    names(clusterAssignments) <- sampleKey
    cls <- lapply(unique(clusterAssignments), function(i) as.numeric(names(clusterAssignments[clusterAssignments %in% i])))
    for (i in 1:length(cls)) {
      nelts <- 1:ncol(m)
      cl <- as.numeric(nelts %in% cls[[i]])
      updt <- outer(cl, cl)
      m <- m + updt
    }
    return(m)
  }

  triangle <- function (m, mode = 1){
    n = dim(m)[1]
    nm = matrix(0, ncol = n, nrow = n)
    fm = m
    nm[upper.tri(nm)] = m[upper.tri(m)]
    fm = t(nm) + nm
    diag(fm) = diag(m)
    nm = fm
    nm[upper.tri(nm)] = NA
    diag(nm) = NA
    vm = m[lower.tri(nm)]
    if (mode == 1) {
      return(vm)
    }
    else if (mode == 3) {
      return(fm)
    }
    else if (mode == 2) {
      return(nm)
    }
  }

  myPal <- function (n = 10){
    seq = rev(seq(0, 255, by = 255/(n)))
    palRGB = cbind(seq, seq, 255)
    rgb(palRGB, maxColorValue = 255)
  }

  setClusterColors <- function (past_ct, ct, colorU, colorList){
    newColors = c()
    if (length(colorList) == 0) {
      newColors = colorU[ct]
      colori = 2
    }
    else {
      newColors = rep(NULL, length(ct))
      colori = colorList[[2]]
      mo = table(past_ct, ct)
      m = mo/apply(mo, 1, sum)
      for (tci in 1:ncol(m)) {
        maxC = max(m[, tci])
        pci = which(m[, tci] == maxC)
        if (sum(m[, tci] == maxC) == 1 & max(m[pci, ]) ==
              maxC & sum(m[pci, ] == maxC) == 1) {
          newColors[which(ct == tci)] = unique(colorList[[1]][which(past_ct ==
                                                                      pci)])
        }
        else {
          colori = colori + 1
          newColors[which(ct == tci)] = colorU[colori]
        }
      }
    }
    return(list(newColors, colori, unique(newColors)))
  }

  CDF <- function (ml, breaks = 100){
    plot(c(0), xlim = c(0, 1), ylim = c(0, 1), col = "white",
         bg = "white", xlab = "consensus index", ylab = "CDF",
         main = "consensus CDF", las = 2)
    k = length(ml)
    this_colors = rainbow(k - 1)
    areaK = c()
    for (i in 2:length(ml)) {
      v = triangle(ml[[i]], mode = 1)
      h = hist(v, plot = FALSE, breaks = seq(0, 1, by = 1/breaks))
      h$counts = cumsum(h$counts)/sum(h$counts)
      thisArea = 0
      for (bi in 1:(length(h$breaks) - 1)) {
        thisArea = thisArea + h$counts[bi] * (h$breaks[bi +
                                                         1] - h$breaks[bi])
        bi = bi + 1
      }
      areaK = c(areaK, thisArea)
      lines(h$mids, h$counts, col = this_colors[i - 1], lwd = 2,
            type = "l")
    }
    legend(0.8, 0.5, legend = paste(rep("", k - 1), seq(2, k,
                                                        by = 1), sep = ""), fill = this_colors)
    deltaK = areaK[1]
    for (i in 2:(length(areaK))) {
      deltaK = c(deltaK, (areaK[i] - areaK[i - 1])/areaK[i -
                                                           1])
    }
    plot(1 + (1:length(deltaK)), y = deltaK, xlab = "k", ylab = "relative change in area under CDF curve",
         main = "Delta area", type = "b")
  }

  clusterTrackingPlot <- function (m){
    plot(NULL, xlim = c(-0.1, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "samples", ylab = "k", main = "tracking plot")
    for (i in 1:nrow(m)) {
      rect(xleft = seq(0, 1 - 1/ncol(m), by = 1/ncol(m)), ybottom = rep(1 -
                                                                          i/nrow(m), ncol(m)), xright = seq(1/ncol(m), 1, by = 1/ncol(m)),
           ytop = rep(1 - (i - 1)/nrow(m), ncol(m)), col = m[i,
                                                             ], border = NA)
    }
    xl = seq(0, 1 - 1/ncol(m), by = 1/ncol(m))
    segments(xl, rep(-0.1, ncol(m)), xl, rep(0, ncol(m)), col = "black")
    ypos = seq(1, 0, by = -1/nrow(m)) - 1/(2 * nrow(m))
    text(x = -0.1, y = ypos[-length(ypos)], labels = seq(2, nrow(m) +
                                                           1, by = 1))
  }

  CDFdk <- function (ml, breaks = 100){
    k = length(ml)
    areaK = c()
    for (i in 2:length(ml)) {
      v = triangle(ml[[i]], mode = 1)
      h = hist(v, plot = FALSE, breaks = seq(0, 1, by = 1/breaks))
      h$counts = cumsum(h$counts)/sum(h$counts)
      thisArea = 0
      for (bi in 1:(length(h$breaks) - 1)) {
        thisArea = thisArea + h$counts[bi] * (h$breaks[bi + 1] - h$breaks[bi])
        bi = bi + 1
      }
      areaK = c(areaK, thisArea)
    }
    deltaK = areaK[1]
    for (i in 2:(length(areaK))) {
      deltaK = c(deltaK, (areaK[i] - areaK[i - 1])/areaK[i - 1])
    }
    return(deltaK)
  }


  ## CCP Main Body

  if (is.null(seed) == TRUE) {
    seed = timeSeed = as.numeric(Sys.time())
  }

  set.seed(seed)
  if (is.null(ml) == TRUE) {
    if (!class(d) %in% c("dist", "matrix", "ExpressionSet")) {
      stop("d must be a matrix, distance object or ExpressionSet (eset object)")
    }
    if (inherits(d, "dist")) {
      if (is.null(attr(d, "method"))) {
        attr(d, "method") <- distance <- "unknown - user-specified"
      }
      if (is.null(distance) || (distance != attr(d, "method"))) {
        distance <- attr(d, "method")
      }
      if ((!is.null(pFeature)) && (pFeature < 1)) {
        message("Cannot use the pFeatures parameter when specifying a distance matrix as the data object\n")
        pFeature <- 1
      }
      if (!is.null(weightsFeature)) {
        message("Cannot use the weightsFeature parameter when specifying a distance matrix as the data object\n")
        weightsFeature <- NULL
      }
      if (clusterAlg == "km") {
        message("Note: k-means will cluster the distance matrix you provided.  This is similar to kmdist option when suppling a data matrix")
      }
    }
    else {
      if (is.null(distance)) {
        distance <- "pearson"
      }
    }
    if ((clusterAlg == "km") && inherits(distance, "character") && (distance != "euclidean")) {
      message("Note: The km (kmeans) option only supports a euclidean distance metric when supplying a data matrix.  If you want to cluster a distance matrix using k-means use the 'kmdist' option, or use a different algorithm such as 'hc' or 'pam'.  Changing distance to euclidean")
      distance <- "euclidean"
    }
    if (inherits(d, "ExpressionSet")) {
      d <- exprs(d)
    }

    ml <- ccRun(d = d, maxK = maxK, repCount = reps, diss = inherits(d, "dist"),
                pItem = pItem, pFeature = pFeature, innerLinkage = innerLinkage,
                clusterAlg = clusterAlg, weightsFeature = weightsFeature,
                weightsItem = weightsItem, distance = distance, verbose = verbose,
                corUse = corUse)
  }
  res = list()
  if ((is.null(plot) == FALSE | writeTable) & !file.exists(paste(title, sep = ""))) {
    dir.create(paste(title, sep = ""))
  }
  log <- matrix(ncol = 2, byrow = T, c("title", title, "maxK",
                                       maxK, "input matrix rows", ifelse(inherits(d, "matrix"), nrow(d), "dist-mat"), "input matrix columns", ifelse(inherits(d, "matrix"), ncol(d), ncol(as.matrix(d))), "number of bootstraps",
                                       reps, "item subsampling proportion", pItem, "feature subsampling proportion",
                                       ifelse(is.null(pFeature), 1, pFeature), "cluster algorithm",
                                       clusterAlg, "inner linkage type", innerLinkage, "final linkage type",
                                       finalLinkage, "correlation method", distance, "plot",
                                       if (is.null(plot)) NA else plot, "seed", if (is.null(seed)) NA else seed))
  colnames(log) = c("argument", "value")
  if (writeTable) {
    write.csv(file = paste(title, "/", title, ".log.csv",
                           sep = ""), log, row.names = F)
  }
  if (is.null(plot)) {
  }
  else if (plot == "pngBMP") {
    bitmap(paste(title, "/", "consensus%03d.png", sep = ""))
  }
  else if (plot == "png") {
    png(paste(title, "/", "consensus%03d.png", sep = ""))
  }
  else if (plot == "pdf") {
    pdf(onefile = TRUE, paste(title, "/", "consensus.pdf", sep = ""))
  }
  else if (plot == "ps") {
    postscript(onefile = TRUE, paste(title, "/", "consensus.ps", sep = ""))
  }
  colorList = list()
  colorM = rbind()
  thisPal <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
               "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6",
               "#6A3D9A", "#FFFF99", "#B15928", "#bd18ea", "#2ef4ca",
               "#f4cced", "#f4cc03", "#05188a", "#e5a25a", "#06f106",
               "#85848f", "#000000", "#076f25", "#93cd7f", "#4d0776",
               "#ffffff")
  colBreaks = NA
  if (is.null(tmyPal) == TRUE) {
    colBreaks = 10
    tmyPal = myPal(colBreaks)
  }
  else {
    colBreaks = length(tmyPal)
  }
  sc = cbind(seq(0, 1, by = 1/(colBreaks)))
  rownames(sc) = sc[, 1]
  sc = cbind(sc, sc)
  heatmap(sc, Colv = NA, Rowv = NA, symm = FALSE, scale = "none",
          col = tmyPal, na.rm = TRUE, labRow = rownames(sc), labCol = F,
          main = "consensus matrix legend")
  for (tk in 2:maxK) {
    if (verbose) {
      message(paste("consensus ", tk))
    }
    fm = ml[[tk]]
    hc = hclust(as.dist(1 - fm), method = finalLinkage)
    ct = cutree(hc, tk)
    names(ct) = colnames(d)
    if (class(d) == "dist") {
      names(ct) = colnames(as.matrix(d))
    }
    c = fm
    colorList = setClusterColors(res[[tk - 1]][[3]], ct,
                                 thisPal, colorList)
    pc = c
    pc = pc[hc$order, ]
    if (!is.null(plot) && plot == "pngBMP") {
      pc = pc[, hc$order]
      pc = rbind(pc, 0)
      oc = colorList[[1]][hc$order]
      heatmap(pc, Colv = NA, Rowv = NA, symm = FALSE, scale = "none",
              col = tmyPal, na.rm = TRUE, labRow = F, labCol = F,
              mar = c(5, 5), main = paste("consensus matrix k=", tk, sep = ""), ColSideCol = oc)
    }
    else {
      pc = rbind(pc, 0)
      heatmap(pc, Colv = as.dendrogram(hc), Rowv = NA,
              symm = FALSE, scale = "none", col = tmyPal, na.rm = TRUE,
              labRow = F, labCol = F, mar = c(5, 5), main = paste("consensus matrix k=", tk, sep = ""), ColSideCol = colorList[[1]])
    }
    legend("topright", legend = unique(ct), fill = unique(colorList[[1]]), horiz = FALSE)
    res[[tk]] = list(consensusMatrix = c, consensusTree = hc, consensusClass = ct, ml = ml[[tk]], clrs = colorList)
    colorM = rbind(colorM, colorList[[1]])
  }


  deltaK <- CDFdk(ml)

  CDF(ml)
  clusterTrackingPlot(colorM[, res[[length(res)]]$consensusTree$order])
  if (is.null(plot) == FALSE) {
    dev.off()
  }
  res[[1]] = colorM
  res[[length(res)+1]] = deltaK

  classes <- data.frame(rep(NA, ncol(d)), check.names=F)
  for (i in 2:(length(res)-1)){
    classes[,i-1] <- res[[i]][[3]]
  }
  colnames(classes) <- c(seq(2, maxK))

  res[[length(res)+1]] = classes

  if (writeTable) {
    for (i in 2:(length(res)-2)) {
      write.csv(file = paste(title, "/", title, ".k=", i, ".consensusMatrix.csv", sep = ""), res[[i]][[1]])
      write.table(file = paste(title, "/", title, ".k=", i, ".consensusClass.csv", sep = ""), res[[i]][[3]], col.names = F, sep = ",")
      write.table(file = paste(title, "/", title, "deltaCDF.csv", sep = ""), res[[length(res)]], col.names = F, sep = ",")
    }
  }
  return(res)
}

commonProp <- function(l){
  #
  # Finds the relative proportion of common samples between
  # clusters from different algorithms, i.e. (A ∩ B)/min{|A|, |B|}.
  # Then calculates the euclidean and cosine distance metrics for the propotions.
  # Creates a pdf heatmap to show distances.
  # Returns an array of dimension n x n x 2 containing the distances between
  # all combinations of the n clusters, in the two distance metrics.
  #
  # labelledClasses  Character matrix. Algorithms in columns and samples in rows, showing the class membership
  #                  of each sample in each algorithm. Clusters from each algorithm should be easily distinguishable
  #                  ('labelled'), so that e.g. heirarchical cluster 1 and k-means cluster 1 have different names
  # clusterTable     Numeric matrix. An n x n matrix giving the number of samples shared between two clusters
  #

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

hyperClust2 <- function(l){
  #
  # Compares all clusters found by different algorithms via a
  # hypergeometric test and creates a pdf heatmap to show
  # the resulting p values.
  # Returns an n x n matrix of the resulting adjusted p-values,
  # where n is the total number of clusters found across the algorithms.
  #
  # clusterTable  Numeric matrix. An n x n matrix giving the number of samples shared between two clusters
  # outfile       Character string. Prepended to the p value and heatmap filenames
  #

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

  tryCatch(heatmap.2(adjust.p, trace = "none", density.info = "density", denscol = "black", scale = "none", col = l$colRain, breaks = c(0, seq(0.05, 1, length.out = length(l$colRain)))), error = function(e) e)

#   correctInput <- FALSE
#   wasError <- FALSE
#   browser()
#   while (correctInput == FALSE){
#     wasError <- tryCatch(heatmap.2(adjust.p, trace = "none", density.info = "density", denscol = "black", scale = "none", col = l$colRain), error = function(e) {print("The cutoff chosen only creates one metacluster. Please try again with a lower cutoff."); wasError <- TRUE; return(wasError)})
#     if (wasError == FALSE){
#       correctInput <- TRUE
#     }
#   }
  dev.off()

  return(adjust.p)
}

inputCond <- function(init = 0, prompt = "Please enter a value for x: ", variable = "x", condition = "x < 3"){
  #
  # Takes user input and tests it against a given condition.
  #
  # init      Numeric. The initial value of the variable
  # prompt    Character string. Printed to prompt the user for input
  # variable  Character string. The name of the variable being tested in the condition
  # condition Character string. The logical test of the variable. If true, the user input is returned.
  #
  x <- variable
  assign(x, init)
  while (x == init){

    x <- as.numeric(readline(paste(prompt)))
    if (eval(parse(text = condition)) == F){
      cat(paste("Please choose another value such that ", condition, sep = ""))
      x <- init
    }
    else {
      return(x)
    }
  }
}

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
      #classPairsCut <- classPairsM[classPairsM[,3] <= 0.2,]
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
        #subCols <- c("#EB1E28FF", "#7CCB2CFF", "#5C1219FF", "#3C4BA7FF", "#B834A1FF")
        #names(subCols) <- c('Enterocyte', 'Goblet-like', 'Inflammatory', 'Stem-like' ,'TA')
        #names(subCols) <- c('basal', 'erbb2', 'lumA', 'lumb', 'norm')
        names(subCols) <- unique(na.omit(l$pheno))

        #         V(NEgraph)$color <- V(NEgraph)$frame.color <- sapply(1:nrow(l$knownSubsP), function(i){
        #           desat(subCols[as.character(l$knownSubsP$kSubtypes)[i]], 1-rescale(l$knownSubsP$pVals, 0, 0.9)[i])
        #           })
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
      #edgeBreaks <- cumsum(seq(0, 2, by = 0.1)) %>%
      #extract(., c(1, which(max(E(NEgraph)$weight) > .)+1))
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

      #      NElist <- list()

      names <- names(l$numTable)

      clustList <- communities(commun)
      wasError <- FALSE

      names(clustList) <- LETTERS[1:length(clustList)]

      NEClust <- data.frame(algClust = names, NEclust = NA)
      for (j in 1:nrow(NEClust)){
        NEClust[j, 2] <- names(clustList[which(lapply(clustList, function(x) names[j] %in% x) == T)])
      }
      write.table(NEClust, paste(Sys.Date(), recons[i], "NE_metaclusters.txt", sep = "_"), sep = "\t", quote = F)

#       old <- FALSE ####
#       if (old == TRUE){
#           # Count how many times a sample appears in a metacluster ###
#         NEarr <- matrix(0, nrow(l$labelledClasses), length(unique(NEClust[,2])))
#         rownames(NEarr) <- rownames(l$labelledClasses)
#         colnames(NEarr) <- LETTERS[1:ncol(NEarr)]
#         for (h in 1:nrow(NEarr)){
#           for (j in 1:ncol(NEarr)){
#             for (m in 1:nrow(NEClust)){
#               if ((NEClust[m, 2] == LETTERS[j]) && (as.character(NEClust[m, 1]) %in% as.matrix(l$labelledClasses[h,]))){
#                 NEarr[h, j] <- NEarr[h, j] + 1
#               }
#             }
#           }
#         }
#         write.table(NEarr, paste(Sys.Date(), recons[i], "NE_samplecount.txt", sep = "_"), sep = "\t", quote = F)
#
#         # Find which samples are common to clusters in a metacluster
#         NEind <- which((apply(NEarr, 1, max) != 1) & apply(NEarr, 1, function(x){sum(x == max(x))}) == 1)
#         NEsub <- NEarr[NEind,]
#         NEsub <- as.matrix(NEsub)
#         colnames(NEsub) <- colnames(NEarr)
#         #NEmemb <- rep(NA, length.out = nrow(NEsub))
#         NEmemb <- colnames(NEsub)[apply(NEsub, 1, function(x){which(x == max(x))})]
#         names(NEmemb) <- rownames(NEsub)
# #        #       NEexp <- data_set[,NEind]##
#         #       for(j in 1:ncol(NEsub)){
#         #         for(m in 1:nrow(NEsub)){
#         #           if(NEsub[m, j] == max(NEsub[m,])){
#         # #             colnames(NEexp)[m] <- paste(colnames(NEsub)[j], colnames(NEexp)[m])
#         #             NEmemb[rownames(NEsub)[m]] <- colnames(NEsub)[j]
#         #           }
#         #         }
#         #       }
#
#
#         # Find outlying metaclusters (those which do not include each algorithm at least once)
#         # and which samples are common to all clusters within them
#         #NEoutclust <- which((apply(NEarr, 2, max) != 1) & (apply(NEarr, 2, max) != length(l$clusterAlg)))
#         NEuniqalg <- lapply(clustList, function(x){length(unique(str_extract(x, '[A-Z]+')))})
#         NEoutclust <- names(NEuniqalg)[NEuniqalg != 1 & NEuniqalg != length(l$clusterAlg)]
#
#         if (length(NEoutclust) != 0){
#           NEoutsub <- NEarr[,NEoutclust]
#           NEoutsub <- as.matrix(NEoutsub)
#           NEoutlier <- which((apply(NEoutsub, 1, max) != 0) & (apply(NEoutsub, 1, max) != 1))
#
#           if (length(NEoutlier) != 0){
#             NEoutsub <- as.matrix(NEoutsub[NEoutlier,])
#             colnames(NEoutsub) <- colnames(NEarr)[NEoutclust]
#             NEoutmemb <- colnames(NEoutsub)[apply(NEoutsub, 1, function(x){which(x == max(x))})]
#             names(NEoutmemb) <- rownames(NEoutsub)
# #            #           NEoutexp <- data_set[,NEoutlier]##
#             #           for(j in 1:ncol(NEoutsub)){
#             #             for(m in 1:nrow(NEoutsub)){
#             #               if(NEoutsub[m, j] == max(NEoutsub[m,])){
#             # #                 colnames(NEoutexp)[m] <- paste(colnames(NEoutsub)[j], colnames(NEoutexp)[m])
#             #                 NEoutmemb[rownames(NEoutsub)[m]] <- colnames(NEoutsub)[j]
#             #               }
#             #             }
#             #           }
#           }
#
#           # Metacluster assignments of all samples
#           NEall <- rep(NA, nrow(NEarr))
#           names(NEall) <- rownames(NEarr)
#           NEall[c(names(NEmemb), names(NEoutmemb))] <- c(NEmemb, NEoutmemb)
#         }
#
#         #      else {
#         NEinclust <- which(apply(NEarr, 2, max) != 1)
#         NEinsub <- as.matrix(NEarr[,NEinclust])
#         NEinsamp <- which(apply(NEinsub, 1, max) == length(l$clusterAlg))
#         NEinsub <- as.matrix(NEinsub[NEinsamp,])
#         colnames(NEinsub) <- colnames(NEarr)[NEinclust]
#         NEall <- rep(NA, length.out = nrow(NEarr))
#         names(NEall) <- rownames(NEarr)
#         NEall[rownames(NEinsub)] <- colnames(NEinsub)[apply(NEinsub, 1, which.max)]
#         #      }###
#         # cat(recons[i], NEall, file = allClust, sep = c(rep("\t", nrow(l$labelledClasses)), "\n"), append = T)
#
#         #       NEexp <- NEexp[,order(colnames(NEexp))]
#         #       write.table(NEexp, file = paste(Sys.Date(), recons[i], cut, "node_edge_expression.txt", sep = "_"), sep = "\t", quote = F)
#       }

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

#        # Dealing with single clusters#####
        #         if (length(single != 0)){
        #
        #           allSingleSamples <- vector()
        #           singleMat <- matrix(nrow = nrow(l$labelledClasses), ncol = 1)
        #           rownames(singleMat) <- rownames(l$labelledClasses)
        #           singleClust <- LETTERS[(ncol(NEarr)+1-length(single)):(ncol(NEarr))]
        #           browser()
        #           for (p in 1:length(single)){
        #             singleSamples <- names(which(rowSums(l$labelledClasses == singleNames[p]) == 1))
        #             singleMat[singleSamples,] <- singleClust[p]
        #
        #             allSingleSamples <- c(allSingleSamples, singleSamples)
        #           }
        #
        #           allSingleSamples <- allSingleSamples[!(allSingleSamples %in% allSingleSamples[duplicated(allSingleSamples)])]
        #           singleSub <- cbind(l$pheno1, singleMat)
        #           l$singleTable <- table(singleSub[,1], singleSub[,2])
        #
        #           hyperClust3(t(l$singleTable), paste(Sys.Date(), "singleclusters", recons[i], cut, "known", sep = "_"), l)
        #         }####
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

pamCentroids <- function(l){
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

# performClust <- function(data, maxK, reps, l$clusterAlg = c("hc", "km", "nmf"), nmfData = NULL){
performClust <- function(l, nmfData = NULL){
  #
  # Performs clustering via the nmf and CCP functions.
  # Also calculates gap score and cophenetic corellation coefficient for each consensus matrix,
  # and creates pdf plots of these.
  # Returns a list of: a matrix of the cophenetic coefficients (algorithms in columns, k-values in rows); an
  # an array of all cluster memberships (algorithms x samples x k-value); and a matrix of gap scores.
  #
  # data        Numeric matrix. Expression values with genes in rows and samples in columns
  # maxK        Integer. The maximum number of clusters to evaluate
  # reps        Integer. The number of resamplings (for hc, km and pm) or repetitions (for nmf) for consensus clustering
  # l$clusterAlg  Character vector. The clustering algorithms to be applied to the data. Any combination of "hc" (heirarchical),
  #             "km" (k-means), "pm" (partitioning around medoids), "nmf" (nonnegative matrix factorisation)
  #

  setwd(l$initTitle)
  if (!all(l$clusterAlg %in% c("hc", "pm", "km", "nmf"))){
    writeLines(paste(l$clusterAlg[!l$clusterAlg %in% c("hc","pm", "km", "nmf")], "is unsupported. Please choose from \"hc\", \"pm\", \"km\" or \"nmf\"."))
  }

  l$consensusMemb <- array(0, c(length(l$clusterAlg), ncol(l$data), l$maxK-1)) # An empty array to store the class membership for all K for all algorithms
  consensusMat <- array(0, c(length(l$clusterAlg), l$maxK-1, ncol(l$data), ncol(l$data))) # An empty array to store the consensus matrices for HC and KM

  # Perform the clustering via each algorithm and
  # store the consensus matrices and cophenetic coefficients

  if ("hc" %in% l$clusterAlg){
    writeLines("Starting hierarchical clustering...")
    hcccp <- CCP(l$data, maxK = l$maxK, reps = l$reps, clusterAlg = "hc", distance = "pearson", title=paste(Sys.Date(), "hc", sep="_"))
    l$consensusMemb[which(l$clusterAlg == "hc"),,] <- data.matrix(hcccp[[length(hcccp)]])
    for (i in 2:(length(hcccp)-2)){
      consensusMat[which(l$clusterAlg == "hc"), i-1,,] <- hcccp[[i]]$consensusMatrix
    }
    writeLines("Completed hierarchical clustering.")
  }

  if ("km" %in% l$clusterAlg){
    writeLines("Starting k-means clustering...")
    kmccp <- CCP(l$data, maxK = l$maxK, reps = l$reps, clusterAlg = "km", distance = "euclidean", title = paste(Sys.Date(), "km", sep="_"))
    l$consensusMemb[which(l$clusterAlg == "km"),,] <- data.matrix(kmccp[[length(kmccp)]])
    for (i in 2:(length(kmccp)-2)){
      consensusMat[which(l$clusterAlg == "km"), i-1,,] <- kmccp[[i]]$consensusMatrix
    }
    writeLines("Completed k-means clustering.")
  }

  if ("pm" %in% l$clusterAlg){
    writeLines("Starting partitioning around medoids clustering...")
    pmccp <- CCP(l$data, maxK = l$maxK, reps = l$reps, clusterAlg = "pam", distance = "euclidean", title = paste(Sys.Date(), "pm", sep="_"))
    l$consensusMemb[which(l$clusterAlg == "pm"),,] <- data.matrix(pmccp[[length(pmccp)]])
    for (i in 2:(length(pmccp)-2)){
      consensusMat[which(l$clusterAlg == "pm"), i-1,,] <- pmccp[[i]]$consensusMatrix
    }
    writeLines("Completed partitioning around medoids clustering.")
  }

  if ("nmf" %in% l$clusterAlg){
    writeLines("Starting nonnegative matrix factorization clustering (this may take some time)...")

    if(is.null(nmfData)){
#       nmfccp <- NMF::nmf(l$data, k.init = 2, k.final = l$maxK, num.clusterings = l$reps, maxniter = 1000, error.function = "euclidean", doc.string = paste(Sys.Date(), "nmf", sep="_"))
#       save(nmfccp, file = paste(Sys.Date(), "_nmf_result.rda", sep = ""))
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
    writeLines("Completed nonnegative matrix factorization clustering.")
  }


  ##############################################################################################################
  # Another if() must be added here if another algorithm is required. It must add the class membership to      #
  # l$consensusMemb (samples in rows, k-values in columns) and the consensus matrices to consensusMat.         #
  ##############################################################################################################

  ## Cophenetic correlation coefficient for HC, PM and KM
  cophen <- t(apply(consensusMat, c(1, 2), cophcor))
  colnames(cophen) <- l$clusterAlg
  rownames(cophen) <- 2:l$maxK



#   for (i in 1:dim(consensusMat)[1]){ ####
#     kVec <- vector(mode = "numeric", length = l$maxK - 1)
#     cophCoef <- vector(mode = "numeric", length = l$maxK - 1)
#
#     for (j in 2:l$maxK){
#       # Distance matrix
#       distMat <- 1 - consensusMat[i, j - 1,,]
#       distMat <- as.dist(distMat)
#       distHC <- hclust(distMat, method="average")
#
#       # Calculate coefficient
#       distCoph <- cophenetic(distHC)
#       kVec[j - 1] <- j
#       cophCoef[j - 1] <- cor(distMat, distCoph)
#       cophCoef[j - 1] <- signif(cophCoef[j - 1], digits = 4)
#       cophen[,i] <- cophCoef
#     }
#   }

#   if ("nmf" %in% l$clusterAlg){
#     cophen[,ncol(cophen)] <- nmfccp$measures$cophenetic
#   } ####

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

#   for (i in 1:length(l$clusterAlg)){ ####
#     clusFun <- function(data, k){
#       # The clustering function for use in clusGap simply returns the already calculated assignments
#       memb <- list(cluster = l$consensusMemb[i,,(k-1)])
#       names(memb$cluster) <- rownames(data)
#       return(memb)
#     }
#     options(warn = -1)
#     gaps <- clusGap(t(data), FUNcluster = clusFun, K.max = maxK, B = 50, verbose = T) # Calculate gap scores
#     options(warn = 0)
#     gapsTab <- gaps$Tab[-1,]
#     gapsTab <- gapsTab[,3:4]
#     gapsArr[,,i] <- gapsTab
#   }
#
#   CImat <- data.frame(K = rep(2:maxK, length(l$clusterAlg)), Gap = as.vector(gapsArr[,1,]), Err = as.vector(gapsArr[,2,]))
#
#   # Plot gap scores
#   pdf(paste(Sys.Date(), "_all_gap.pdf", sep = ""), height = 7, width = 9)
#   plotCI(CImat[,1], CImat[,2], uiw = CImat[,3], liw = CImat[,3], col = c(rep(l$colPal, each = (maxK-1))), lty = "11", lwd = 2, gap = 0.2, xlim = c(2, maxK), xlab = "k", ylab = "Gap Score", main = "Gap Score")
#   matplot(rbind(rep(NA, length(l$clusterAlg)), gapsArr[,1,]), type = "c", col = l$colPal, lty = 1, add = T)
#   legend("topleft", legend = l$clusterAlg, col=l$colPal, pch=1)
#   dev.off()
#####

  l$clustRes <- list(cophen, l$consensusMemb, gapsArr)
  names(l$clustRes) <- c("cophen", "l$consensusMemb", "gapsArr")
  setwd(l$wd)
  return(l)
}

#protoType <- function(l){
  #
  # Compares clusters found by different algorithms via the protovectors method
  # described in "Using Cluster Ensemble and Validation to Identify
  # Subtypes of Pervasive Developmental Disorders", Shen et. al. (2007)
  # Then uses hyperClust3 to do a hypergeometric test between prototype subtypes
  # and known phenotypes.
  # Returns a list of all the samples present in each protovector.
  #
  # samplesClasses  Numeric matrix. With algorithms in columns and samples in rows, this gives the class assignments
  #                 of each sample according to each algorithm. Classes should be named by numbers only
  # phenoF          Character vector. The known phenotypes of the samples, with the name of each element being the sample name.
  #

  # Find the prototypes represented in the sample
  sampleClasses <- sapply(l$labelledClasses, str_extract_all, pattern='[0-9]{1,2}')
  sampleClasses <- sapply(data.frame(sampleClasses), as.numeric)
  rownames(sampleClasses) <- rownames(l$labelledClasses)
  protoVec <- data.frame(unique(sampleClasses))
  agglist <- list()
  for(i in 1:ncol(sampleClasses)){
    agglist[[i]] <- as.numeric(sampleClasses[,i])
  }

  agg <- aggregate(sampleClasses, by = agglist, FUN = sum)
  agg$freq <- agg[,ncol(agg)]/agg[,ncol(agg)-length(l$clusterAlg)] # Find how many samples are represented by each prototype

  for (i in length(l$clusterAlg):1){
    protoVec <- protoVec[order(protoVec[,i]),]
    agg <- agg[order(agg[,i]),]
  }

  # Order prototypes by frequency and name as letters
  protoFreq <- data.frame(protoVec)
  protoFreq$freq <- agg$freq
  protoVec <- protoVec[order(protoFreq$freq, decreasing=T),]
  protoFreq <- protoFreq[order(protoFreq$freq, decreasing=T),]

  names <- LETTERS
  for (i in 1:(floor(nrow(protoVec)/26))){
    names <- c(names, paste(LETTERS, i, sep = ""))
  }

  rownames(protoVec) <- names[1:nrow(protoVec)]
  rownames(protoFreq) <- names[1:nrow(protoVec)]

  # Plot prototype frequency
  pdf(paste(l$protTitle, Sys.Date(), "_", "prototype_frequency.pdf", sep=""), height=13, width=13)
  plot(protoFreq$freq, xlab = "Protovector", ylab = "Sample Frequency", col = "transparent", col.axis = "transparent")
  axis(1, labels = rownames(protoFreq), at = 1:nrow(protoFreq))
  axis(2, labels = 1:max(protoFreq$freq), at = 1:(max(protoFreq$freq)))
  for (i in 1:nrow(protoFreq)){
    abline(v = i, col = 8)
  }
  for (i in 1:max(protoFreq$freq)){
    abline(h = i, col = 8)
  }
  lines(protoFreq$freq, lwd = 1.5)
  points(protoFreq$freq)
  dev.off()

  # Find which samples are represented by which prototypes
  protoSamples <- array(NA, c(nrow(sampleClasses), nrow(protoVec)))
  for (i in 1:nrow(protoVec)){
    for (j in 1:nrow(sampleClasses)){
      if ((sum(protoVec[i,] == sampleClasses[j,]) == length(l$clusterAlg))){
        if (is.na(protoSamples[1, i]) == F){
          protoSamples[min(which(is.na(protoSamples[,i]))), i] <- rownames(sampleClasses)[j]
        }
        else {
          protoSamples[1, i] <- rownames(sampleClasses)[j]
        }
      }
    }
  }

  # Convert array with columns of differing numbers of NAs for each prototype into a neat list of samples in each prototype
  protoSamplesL <- list()
  for (i in 1:ncol(protoSamples)){
    protoSamplesL[[i]] <- na.omit(protoSamples[,i])
    attr(protoSamplesL[[i]], "na.action") <- NULL
  }
  names(protoSamplesL) <- rownames(protoVec)

  # Enter optimal number of prototype vectors
  cond <- paste("x > 1 & x <= ", nrow(protoVec), " & !is.na(x)", sep = "")
  protoOpti <- inputCond(0, "Enter cutoff number of prototype vectors: ", "x", cond)

  # Use k-modes clustering to pair the small prototypes to a larger prototype
  kmod <- kmodes(protoVec, protoVec[1:protoOpti,], 100)
  kmod$cluster <- rownames(protoVec)[kmod$cluster]
  names(kmod$cluster) <- rownames(protoVec)
  protoClose <- kmod$cluster

  write.table(protoClose, file = paste(Sys.Date(), "prototype_pairs.txt", sep = "_"), quote = F, sep = "\t")

  # Append those samples from small prototypes to their closest larger prototype
  protoApp <- protoSamplesL
  for (i in (protoOpti+1):length(protoClose)){
    protoApp[[protoClose[i]]] <- append(protoApp[[protoClose[i]]], protoApp[[names(protoClose)[i]]])
  }
  protoApp <- protoApp[1:protoOpti]

  # Vector of prototype assignments for all samples
  protoTab <- vector("character", nrow(l$labelledClasses))
  names(protoTab) <- rownames(l$labelledClasses)
  for (i in 1:length(protoTab)){
    for (j in 1:length(protoApp)){
      if (names(protoTab)[i] %in% protoApp[[j]]){
        protoTab[i] <- names(protoApp)[j]
      }
    }
  }

  # Plot silhouette width
  pdf(paste(Sys.Date(), "_", "prototypes_silhouette.pdf", sep = ""), height = 7, width = 9)
  sk <- silhouette(na.omit(as.numeric(as.factor(protoTab))), dist(t(l$data[,!is.na(protoTab)])))
  plot(sk, main = "Silhouette width")
  dev.off()

  cat("prototypes", protoTab, file = l$allClust, sep = c(rep("\t", length(protoTab)), "\n"), append = T)
  write.table(sapply(protoSamplesL, '[', seq(max(sapply(protoSamplesL, length)))), file = paste(Sys.Date(), protoOpti, "prototype_samples.txt", sep = "_"), quote = F, sep = "\t", na = "")
  write.table(sapply(protoApp, '[', seq(max(sapply(protoApp, length)))), file = paste(Sys.Date(), protoOpti, "combined_prototype_samples.txt", sep = "_"), quote = F, sep = "\t", na = "")
  write.table(protoFreq, file = paste(Sys.Date(), protoOpti, "prototype_frequency.txt", sep = "_"), quote = F, sep = "\t")
  write.table(protoClose, file = paste(Sys.Date(), protoOpti, "prototype_pairs.txt", sep = "_"), quote = F, sep = "\t")

  # Hypergeometric test between known and found subtypes
  if (!is.null(l$phenoFile)){
    hyperClust3(table(l$pheno[names(protoTab)], protoTab), paste(Sys.Date(), "_", protoOpti, "_prototypes_known", sep = ""), l)
  }

  return(protoSamplesL)
}

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

testClust <- function(l){
  #
  # Compares clusters found by clustering algorithms via: the proportion of common samples; a hypergeometric test;
  # and the protovectors method described in "Using Cluster Ensemble and Validation to Identify
  # Subtypes of Pervasive Developmental Disorders", Shen et. al. (2007)
  #
  # clustRes  List. The output of a call to performClust
  #

    sampleClasses <- array(0, c(ncol(l$data), length(l$clusterAlg))) # An array of the best fitting class assignments as found by each algorithm
    rownames(sampleClasses) <- colnames(l$data)
    colnames(sampleClasses) <- l$clusterAlg

    # User input for optimal k
    i <- 1
    while (i <= length(l$clusterAlg)){

      optK <- inputCond(1, paste("Enter optimal k for ", l$clusterAlg[i], ": ", sep = ""), "x", paste0("x > 1 & x <= ", l$maxK, " & !is.na(x)"))

      sampleClasses[,i] <- l$consensusMemb[i,, optK - 1]
      #cat(l$clusterAlg[i], l$consensusMemb[i,, optK - 1], file = l$allClust, sep = c(rep("\t", ncol(data_set)), "\n"), append = T)
      i <- i+1
    }

    l$labelledClasses <- data.frame(sampleClasses)

    write.table(t(l$labelledClasses), file = l$allClust, quote = FALSE, sep = '\t')

    for (i in 1:length(l$clusterAlg)){                                            # Add prefixes to class assignments from each algorithm to
      l$labelledClasses[,i] <- factor(paste(l$clusterAlg[i], l$labelledClasses[,i], sep=""))  # distinguish, e.g., hierarchical cluster 1 from k-means cluster 1
    }
  #
  #   # Create a table showing the number of samples shared between pairs of clusters from two algorithms
  #   for (i in 1:length(clusterAlg)){
  #     for (j in 1:(length(clusterAlg))){
  #       table <- table(labelledClasses[,i], labelledClasses[,j])
  #       if (j == 1){
  #         clusterRow <- table
  #       }
  #       else {
  #         clusterRow <- cbind(clusterRow, table)
  #       }
  #     }
  #     if(i == 1){
  #       clusterTable <- clusterRow
  #     }
  #     else{
  #       clusterTable <- rbind(clusterTable, clusterRow)
  #     }
  #     rm(clusterRow)
  #   }
  #

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
    l$pheno <- l$pheno[rownames(l$labelledClasses)]
    setwd(l$knownTitle)
    allTabs <- sapply(l$labelledClasses, function(x){table(x, l$pheno)}, simplify = FALSE)
    #for (i in 1:length(allTabs)){hyperClust3(allTabs[[i]], paste(Sys.Date(), names(allTabs[i]), 'hyper_known', sep = '_'), l = l)}
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

  #setwd(l$protTitle)
  #protoApp <- protoType(l)

  setwd(l$wd)
  return(l)
}


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

  #sds <- apply(as.matrix(l$data), 1, sd, na.rm = TRUE)
  #l$data <- l$data[sds >= quantile(sds, probs = 0.9, na.rm = TRUE),]

  # Create subfolders to store output
  l$wd <- getwd()
  options(warn = -1)

  analysisFolder <- paste(Sys.Date(), ref, "analysis", sep = "_")
  l$analysisTitle <- paste(l$wd, "/Output/", analysisFolder, "/", sep = "")
  dir.create(l$analysisTitle, recursive = TRUE)

  initFolder <- paste(Sys.Date(), ref, "initial", sep = "_")
  hypFolder <- paste(Sys.Date(), ref, "hypergeometric", sep = "_")
  propFolder <- paste(Sys.Date(), ref, "proportion", sep = "_")
  #protFolder <- paste(Sys.Date(), ref, "prototypes", sep = "_")
  pamFolder <- paste(Sys.Date(), ref, "pam", sep = "_")

  l$initTitle <- paste(l$analysisTitle, initFolder, "/", sep = "")
  l$hypTitle <- paste(l$analysisTitle, hypFolder, "/", sep = "")
  l$propTitle <- paste(l$analysisTitle, propFolder, "/", sep = "")
  #l$protTitle <- paste(l$analysisTitle, protFolder, "/", sep = "")
  l$pamTitle <- paste(l$analysisTitle, pamFolder, "/", sep = "")
  l$knownTitle <- paste0(l$analysisTitle, initFolder, "/", Sys.Date(), "_known_hypergeometric", "/")

  dir.create(l$hypTitle); dir.create(l$propTitle); dir.create(l$initTitle); dir.create(l$knownTitle); dir.create(l$pamTitle)
  #dir.create(l$protTitle)
  
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

