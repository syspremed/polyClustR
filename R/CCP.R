  #' Consensus clustering
  #'
  #' @description Takes an expression matrix and finds clusters via hierarchical clustering, k-means or partitioning around medoids. 
  #'  Adapted from ConsensusClusterPlus of the ConsensusClusterPlus package, see \code{?ConsensusClusterPlus}.
  #'  
  #' @param d              Numeric matrix. Data to be clustered, where columns=items/samples and rows are features
  #' @param maxK           Integer. Maximum number of clusters to evaluate
  #' @param reps           Integer. Number of subsamples so consensus can be evaluated
  #' @param pItem          Numerical value. Proportion of items (columns) to sample in each subsampling
  #' @param pFeature       Numerical value. Proportion of features (rows) to sample in each subsampling
  #' @param clusterAlg     Character string. Cluster algorithm: 'hc' heirarchical (hclust), 'pam' for paritioning around medoids, 'km' for k-means
  #' @param title          Character string. Name for output directory. Directory is created only if plot is not NULL or writeTable is TRUE. This title can be an abosulte or relative path.
  #' @param innerLinkage   heirarchical linkage method for subsampling
  #' @param finalLinkage   heirarchical linkage method for consensus matrix
  #' @param distance       Character string. 'pearson': (1 - Pearson correlation), 'spearman' (1 - Spearman correlation), 'euclidean', 'binary', 'maximum', 'canberra', 'minkowski" or custom distance function.
  #' @param ml             Optional. Prior result. If supplied then only do graphics and tables
  #' @param tmyPal         Optional. Character vector. Colors for consensus matrix
  #' @param seed           Optional. Numerical. Sets random seed for reproducible results
  #' @param plot           Character string. NULL - print to screen, 'pdf', 'png', 'pngBMP' for bitmap png, helpful for large datasets
  #' @param writeTable     Logical. TRUE - write ouput and log to csv
  #' @param weightsItem    Optional. Numerical vector. Weights to be used for sampling items
  #' @param weightsFeature Optional. Numerical vector. Weights to be used for sampling features
  #' @param verbose        Logical. If TRUE, print messages to the screen to indicate progress. This is useful for large datasets
  #' @param corUse         Optional. Character string. Specifies how to handle missing data in correlation distances 'everything','pairwise.complete.obs', 'complete.obs'
  #' 
  #' @details Not intended for use outside of a call to \code{polyCluster}.
  #' @return Returns a list including a lot of information for each k, most importantly the consensus matrices and class assignments
  #'
 
 CCP <- function (d = NULL, maxK = maxK, reps = reps, pItem = 0.8, pFeature = 1,
                  clusterAlg = "hc", title = "untitled_consensus_cluster",
                  innerLinkage = "average", finalLinkage = "average", distance = "euclidean",
                  ml = NULL, tmyPal = NULL, seed = NULL, plot = "pdf", writeTable = TRUE,
                  weightsItem = NULL, weightsFeature = NULL, verbose = F, corUse = "everything") {
  
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