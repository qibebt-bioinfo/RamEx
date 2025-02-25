#' Perform Louvain Clustering on Spectral Data
#'
#' This function performs Louvain clustering on the spectral data stored in a Ramanome
#' object. It first reduces the dimensionality of the data using PCA, then calculates
#' the k-nearest neighbors for each point, and finally applies the Louvain algorithm
#' to detect communities at different resolutions.
#'
#' @param object A Ramanome object containing the spectral data.
#' @param resolutions A vector of resolutions at which to perform the Louvain clustering.
#' @param npc The number of principal components to retain in the PCA. Defaults to 10.
#' @param threshold The minimum number of samples required to consider a cluster valid.
#'                  Defaults to 0.001.
#' @param k The number of nearest neighbors to consider. Defaults to 30.
#' @param n_tree The number of trees to use in the Approximate Nearest Neighbor search.
#'               Defaults to 50.
#' @return A data frame containing the cluster assignments for each resolution.
#' @export Phenotype.Analysis.Louvaincluster
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph cluster_louvain
#' @importFrom igraph membership
#' @importFrom Matrix sparseMatrix
#' @importFrom foreach %dopar%
#' @import foreach foreach
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores
#' @importFrom RcppAnnoy AnnoyEuclidean
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[data_cleaned$index_good,]
#' data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
#' #options(mc.cores = 2)
#' #clusters_Louvaincluster <- Phenotype.Analysis.Louvaincluster(object = data_cleaned, resolutions = c(0.8))

Phenotype.Analysis.Louvaincluster <- function(object, resolutions,npc=10, threshold=0.001, k=30, n_tree=50) {

  data.red <- prcomp_irlba(get.nearest.dataset(object), n = npc, center = TRUE, scale. = TRUE)

  matrix <- data.red$x

  cat('Creating the communities...\n')

  # Calculate the k-nearest neighbor via Annoy
  n <- nrow(matrix)
  ann_index <- new(AnnoyEuclidean, ncol(matrix))

  for (i in 1:n) {
    ann_index$addItem(i - 1, as.numeric(matrix[i, ]))
  }

  ann_index$build(n_tree)

  nearest_neighbors <- t(sapply(0:(n-1), function(i) ann_index$getNNsByItem(i, k)))

  # parallel calculating
  num_cores <- detectCores() - 2
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  clusterExport(cl, list("matrix", "nearest_neighbors", "n"), envir = environment())

  similarities <- foreach(i = 1:n, .combine = 'rbind', .packages = c("Matrix")) %dopar% {
    rows <- rep(i, k)
    cols <- nearest_neighbors[i, ] + 1

    numerator <- rowSums(matrix[rows, ] * matrix[cols, ])
    denominator <- sqrt(rowSums(matrix[rows, ]^2)) * sqrt(rowSums(matrix[cols, ]^2))
    sim_values <- numerator / denominator  # Cosine correlation
    sim_values[sim_values < 0] <- 0
    data.frame(rows = rows, cols = cols, sim = sim_values)
  }

  stopCluster(cl)

  rows <- unlist(similarities$rows)
  cols <- unlist(similarities$cols)
  sims <- unlist(similarities$sim)

  sparse_matrix <- sparseMatrix(i = rows, j = cols, x = sims)

  g <- graph_from_adjacency_matrix(sparse_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)

  min.sample <- floor(threshold * n)
  results <- data.frame(matrix(ncol = length(resolutions), nrow = n))

  for (j in seq_along(resolutions)) {
    resolution <- resolutions[j]
    cat('Processing resolution:', resolution, '\n')

    louvain_clusters <- cluster_louvain(g, resolution = resolution)  #这里好像可以改成leiden聚类，但是我没试过
    clusters <- membership(louvain_clusters)

    leaved_clusters <- which(table(clusters) > min.sample)
    clusters_index <- clusters %in% leaved_clusters
    if (length(leaved_clusters) == 1) {
      warning('There is only one cluster, please set a larger resolution')
    }

    cat('Merging the isolate samples...\n')
    iso_clusters <- which(table(clusters) <= min.sample)
    if (length(iso_clusters) >= 0.1 * n) {
      stop('Too many clusters for the dataset, please set a smaller resolution!')
    }

    if (length(iso_clusters) >= 1) {
      cat(length(iso_clusters), 'clusters contain less than', min.sample, 'samples, So they were merged into the nearest community.\n')
      for (cc in iso_clusters) {
        isolate_index <- clusters == cc
        similarities_temp <- sparse_matrix[isolate_index, clusters_index, drop = FALSE]
        closest_node <- if (length(which(isolate_index)) == 1) {
          which.max(similarities_temp)
        } else {
          which.max(colMeans(as.matrix(similarities_temp)))
        }
        clusters[isolate_index] <- clusters[clusters_index][closest_node]
      }
    }

    results[, j] <- clusters
  }

  colnames(results) <- paste0("Resolution_", resolutions)

  object@meta.data <- cbind(object@meta.data, results)

  return(results)
}



#' k-Means Clustering Analysis
#'
#' A centroid-based clustering algorithm that partitions
#' data into a predefined number of clusters by assigning
#' sample to the nearest center
#'
#' @param object A Ramanome object containing the dataset.
#' @return The datafram which contains Kmeans result.
#' @export Phenotype.Analysis.Kmeans
#' @importFrom stats kmeans
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[data_cleaned$index_good,]
#' data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
#' clusters_kmneans <- Phenotype.Analysis.Kmeans(data_cleaned)

Phenotype.Analysis.Kmeans <- function(object, k) {
  data_x <- get.nearest.dataset(object)
  cl <- kmeans(data_x, k)
  return(list(clusters = cl$cluster, centers = cl$centers))
}


#' Perform Hierarchical Clustering Analysis (HCA)
#'
#' This function performs hierarchical clustering on the dataset retrieved from a
#' Ramanome object using the average linkage method. It calculates the Euclidean distance
#' matrix and then applies the hierarchical clustering algorithm to create a dendrogram.
#' The function also plots the resulting dendrogram.
#'
#' @param object A Ramanome object containing the dataset.
#' @return The hierarchical clustering fit object.
#' @export Phenotype.Analysis.Hca
#' @importFrom stats hclust
#' @importFrom vegan vegdist
#' @importFrom graphics plot
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[data_cleaned$index_good,]
#' data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
#' clusters_hca <- Phenotype.Analysis.Hca(data_cleaned)

Phenotype.Analysis.Hca <- function(object) {
  dataset <- object@datasets$normalized.data
  distance.matrix <- vegdist(dataset, method = "euclidean")
  fit.average <- hclust(distance.matrix, method="average")
  plot(
    fit.average,
    hang=-1,
    cex=.8,
    main="Vaerage Linkage Clustering"
  )
  return(fit.average)
}
