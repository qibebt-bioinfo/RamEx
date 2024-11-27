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
#' @export louvain_clustering
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph cluster_louvain
#' @importFrom Matrix sparseMatrix
#' @importFrom foreach %dopar%
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores
louvain_clustering <- function(object, resolutions,npc=10, threshold=0.001, k=30, n_tree=50) {

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

  # 构建稀疏相似度矩阵，这里是不是可以写数据量<10000的时候做非稀疏的、数据量大的时候做稀疏的矩阵
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


#' Perform Hierarchical Clustering Analysis (HCA)
#'
#' This function performs hierarchical clustering on the dataset retrieved from a
#' Ramanome object using the average linkage method. It calculates the Euclidean distance
#' matrix and then applies the hierarchical clustering algorithm to create a dendrogram.
#' The function also plots the resulting dendrogram.
#'
#' @param object A Ramanome object containing the dataset.
#' @return The hierarchical clustering fit object.
#' @export hca_clustering
#' @importFrom stats hclust
#' @importFrom vegan vegdist
#' @importFrom graphics plot
hca_clustering <- function(object) {
  dataset <- get.nearest.dataset(object)
  distance.matrix <- vegdist(dataset, method = "euclidean")
  fit.average <- hclust(distance.matrix, method="average")
  plot(
    fit.average,
    hang=-1,
    cex=.8,
    main="Vaerage Linkage Clustering"
  )

}


#' Perform Non-negative Matrix Factorization (NMF)
#'
#' This function applies Non-negative Matrix Factorization to the dataset retrieved from
#' a Ramanome object. NMF is a technique used to decompose a non-negative matrix into
#' two non-negative matrices, often used for data dimensionality reduction and feature
#' extraction.
#'
#' @param object A Ramanome object containing the dataset.
#' @return A list containing the results of the NMF decomposition, including the basis
#' matrices and the reconstructed matrix.
#' @export nmf_matrix
#' @importFrom NMF nmf
nmf_matrix <- function(object) {
  dataset <- get.nearest.dataset(object)
  res <- nmf(dataset, 2, method="ns", seed=123456)
  return(res)
}
