#' Louvain Clustering
#'
#' A graph-based community detection algorithm that optimizes modularity to identify clusters by iteratively merging nodes to maximize within-cluster connections.
#'
#' @param object A Ramanome object.
#' @param resolutions A vector of resolutions at which to perform the Louvain clustering.
#' @param n_pc The number of principal components to retain in the PCA. Defaults to 10.
#' @param threshold The minimum number of samples required to consider a cluster valid. Defaults to 0.001.
#' @param k Number of Nearest Neighbors, determines how many nearest neighbors each data point will consider when building the k-NN graph. A larger k will result in a denser graph, where each node is connected to more neighbors. This can capture broader local structures but may blur fine-grained distinctions. A smaller k will produce a sparser graph, emphasizing local structures but potentially missing broader patterns.
#' @param n_tree The number of trees to use in the Approximate Nearest Neighbor search. Defaults to 50.
#'                A higher n_tree improves the accuracy of the nearest neighbor search but increases computation time and memory usage. A lower n_tree speeds up the search but may return less accurate neighbors.
#' @param seed The seed for the random number generator.
#' @return A data frame containing the cluster assignments for each resolution.
#' @export Phenotype.Analysis.Louvaincluster
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph cluster_louvain
#' @importFrom igraph membership
#' @importFrom Matrix sparseMatrix
#' @importFrom parallel makeCluster
#' @importFrom parallel detectCores
#' @importFrom RcppAnnoy AnnoyEuclidean
#' @importFrom data.table rbindlist
#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' #options(mc.cores = 2)
#' #clusters_Louvaincluster <- Phenotype.Analysis.Louvaincluster(object = data_processed, resolutions = c(0.8))

Phenotype.Analysis.Louvaincluster <- function(object, resolutions,n_pc=10, threshold=0.001, k=30, n_tree=50, n_cores=1,seed=42) {
  set.seed(seed)
  if (!is.null(object@reductions$PCA)) {
    if(ncol(object@reductions$PCA) < n_pc) {
      data.red <- prcomp_irlba(get.nearest.dataset(object), n = n_pc, center = TRUE, scale. = TRUE)$x
    } else {
      data.red <- object@reductions$PCA[,1:n_pc]
    }
  } else {
    data.red <- prcomp_irlba(get.nearest.dataset(object), n = n_pc, center = TRUE, scale. = TRUE)$x
  }
  
  matrix <- as.matrix(data.red)
  
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
  num_cores <- n_cores
  cl <- makeCluster(num_cores)
  
  clusterExport(cl, list("matrix", "nearest_neighbors", "n", 'k'), envir = environment())
  
  similarities_list <- parLapply(cl, 1:n, function(i) {
    rows <- rep(i, k)
    cols <- nearest_neighbors[i, ] + 1
    
    numerator <- rowSums(matrix[rows, ] * matrix[cols, ])
    denominator <- sqrt(rowSums(matrix[rows, ]^2)) * sqrt(rowSums(matrix[cols, ]^2))
    sim_values <- numerator / denominator  # Cosine correlation
    sim_values[sim_values < 0] <- 0
    data.frame(rows = rows, cols = cols, sim = sim_values)
  })
  
  stopCluster(cl)
  
  similarities <- rbindlist(similarities_list)
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
    clusters <- as.factor(membership(louvain_clusters))
    
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
#' A centroid-based partitioning algorithm that assigns data points to clusters by minimizing the sum of squared distances between points and their cluster centers
#'
#' @param object A Ramanome object.
#' @param k the number of clusters.
#' @param n_pc The number of principal components to retain in the PCA. Defaults to 10.

#' @return A list containing:
#' \describe{ 
#'   \item{clusters}{The cluster assignments for each data point}
#'   \item{centers}{The coordinates of the cluster centers}
#' }
#' 
#' @export Phenotype.Analysis.Kmeans
#' @importFrom stats kmeans
#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' clusters_kmneans <- Phenotype.Analysis.Kmeans(data_processed,5)

Phenotype.Analysis.Kmeans <- function(object, k, n_pc=10) {
  if (!is.null(object@reductions$PCA)) {
    if(ncol(object@reductions$PCA) < n_pc) {
      data.red <- prcomp_irlba(get.nearest.dataset(object), n = n_pc, center = TRUE, scale. = TRUE)$x
    } else {
      data.red <- object@reductions$PCA[,1:n_pc]
    }
  } else {
    data.red <- prcomp_irlba(get.nearest.dataset(object), n = n_pc, center = TRUE, scale. = TRUE)$x
  }
  cl <- kmeans(as.matrix(data.red), k)
  return(list(clusters = cl$cluster, centers = cl$centers))
}


#' Hierarchical Clustering Analysis (HCA)
#'
#' An iterative clustering method that builds a nested tree of clusters by progressively merging similar samples or dividing dissimilar ones, visualized as a dendrogram showing relationships at multiple scales.
#' The function also plots the resulting dendrogram.
#'
#' @param object A Ramanome object.
#' @param show Whether to plot the dendrogram.
#' @return The hierarchical clustering fit object.
#' @export Phenotype.Analysis.Hca
#' @importFrom stats hclust
#' @importFrom vegan vegdist
#' @importFrom graphics plot
#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' clusters_hca <- Phenotype.Analysis.Hca(data_processed)

Phenotype.Analysis.Hca <- function(object, show = FALSE) {
  dataset <- get.nearest.dataset(object)
  distance.matrix <- vegdist(dataset, method = "euclidean")
  fit.average <- hclust(distance.matrix, method="average")
  if(show){
    plot(
    fit.average,
    hang=-1,
    cex=.8,
    main="Vaerage Linkage Clustering"
  )}
  return(fit.average)
}


#' Gaussian Mixture Model (GMM)
#'
#' A probabilistic model that assumes data is generated from a mixture of Gaussian distributions. 
#' Models are estimated by EM algorithm initialized by hierarchical model-based agglomerative clustering. The optimal model is then selected according to BIC.
#'
#' @param object The Ramanome object
#' @param n_pc The number of principal components to use
#' 
#' @return A vector containing the cluster assignments for each data point
#' @export Phenotype.Analysis.Gmm
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' cluster_gmm <- Phenotype.Analysis.Gmm(data_processed)
Phenotype.Analysis.Gmm <- function(object, n_pc = 20, seed=42) {
  set.seed(seed)
  if (!is.null(object@reductions$PCA)) {
    if(ncol(object@reductions$PCA) < n_pc) {
      data.red <- prcomp_irlba(get.nearest.dataset(object), n = n_pc, center = TRUE, scale. = TRUE)$x
    } else {
      data.red <- object@reductions$PCA[,1:n_pc]
    }
  } else {
    data.red <- prcomp_irlba(get.nearest.dataset(object), n = n_pc, center = TRUE, scale. = TRUE)$x
  }
  
  gmm_model <- Mclust(data.red, verbose = FALSE)

  return(gmm_model$classification)
}

