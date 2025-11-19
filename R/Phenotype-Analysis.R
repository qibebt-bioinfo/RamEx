#' Louvain Clustering
#'
#' Build a shared nearest neighbor (SNN) graph from the latest spectral matrix
#' stored in a Ramanome object and run Louvain community detection across a
#' vector of resolution values.
#'
#' @param object A Ramanome object.
#' @param resolutions Numeric vector of resolution values to evaluate.
#' @param n_pc Number of principal components used to build the neighbor graph.
#' @param threshold Minimum cluster size. Values < 1 are treated as the fraction
#' of cells and values >= 1 are treated as absolute counts.
#' @param k Number of nearest neighbors per cell when constructing the k-NN
#' graph. Larger values encourage coarser clusters.
#' @param n_tree Number of trees used by the Annoy index.
#' @param prune.SNN Pruning threshold applied when computing the SNN Jaccard
#' similarity matrix.
#' @param seed Seed used to make the Louvain results reproducible.
#' @param verbose Whether to print progress messages.
#'
#' @return A data frame whose columns contain the Louvain cluster labels for
#' each requested resolution.
#' @export Phenotype.Analysis.Louvain
#' @importFrom igraph graph_from_adjacency_matrix cluster_louvain membership
#' @importFrom RcppAnnoy AnnoyEuclidean
#' @importFrom irlba prcomp_irlba
#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' clusters <- Phenotype.Analysis.Louvain(
#'   object = data_processed,
#'   resolutions = c(0.5, 0.8, 1.2)
#' )
Phenotype.Analysis.Louvain <- function(
    object,
    resolutions,
    n_pc = 10,
    threshold = 0.001,
    k = 30,
    n_tree = 50,
    prune.SNN = 1 / 15,
    seed = 42,
    verbose = TRUE
) {
  stopifnot(inherits(object, "Ramanome"))
  if (missing(resolutions) || length(resolutions) == 0) {
    stop("Please provide at least one resolution value.")
  }
  if (!is.numeric(resolutions) || any(!is.finite(resolutions))) {
    stop("`resolutions` must be a numeric vector without missing values.")
  }

  embedding <- .ramex_prepare_embedding(object, n_pc = n_pc, verbose = verbose)
  n_cells <- nrow(embedding)
  if (n_cells < 2) {
    stop("At least two spectra are required to perform Louvain clustering.")
  }

  k <- as.integer(k)
  if (k < 1) {
    stop("`k` must be at least 1.")
  }
  if (k >= n_cells) {
    if (verbose) {
      message("Reducing k to the number of spectra minus one.")
    }
    k <- n_cells - 1L
  }

  neighbor_idx <- .ramex_annoy_knn(embedding, k = k, n_trees = n_tree)
  snn_matrix <- .ramex_build_snn_graph(
    neighbor_idx = neighbor_idx,
    prune = prune.SNN,
    cell_names = rownames(embedding)
  )
  graph <- igraph::graph_from_adjacency_matrix(
    snn_matrix,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  cluster_list <- lapply(
    resolutions,
    FUN = function(res) {
      ids <- .ramex_run_louvain(graph = graph, resolution = res, seed = seed)
      ids <- .ramex_filter_clusters(ids, threshold = threshold, n_cells = n_cells)
      factor(ids, exclude = NULL)
    }
  )
  names(cluster_list) <- vapply(resolutions, .ramex_resolution_label, character(1))
  clusters <- as.data.frame(
    cluster_list,
    row.names = rownames(embedding),
    optional = FALSE,
    check.names = FALSE
  )
  return(clusters)
}

# Prepare a PCA embedding for Louvain clustering --------------------------------
.ramex_prepare_embedding <- function(object, n_pc, verbose) {
  dataset <- get.nearest.dataset(object)
  dataset <- as.matrix(dataset)
  if (nrow(dataset) == 0) {
    stop("No spectra available for clustering.")
  }
  if (is.null(rownames(dataset))) {
    rownames(dataset) <- paste0("sample_", seq_len(nrow(dataset)))
  }

  existing_pca <- object@reductions$PCA
  if (!is.null(existing_pca) && ncol(existing_pca) >= n_pc) {
    embedding <- as.matrix(existing_pca[, seq_len(n_pc), drop = FALSE])
  } else {
    max_comp <- min(n_pc, ncol(dataset), max(1, nrow(dataset) - 1))
    if (verbose) {
      message("Computing PCA embedding with ", max_comp, " components.")
    }
    pca_fit <- irlba::prcomp_irlba(dataset, n = max_comp, center = TRUE, scale. = TRUE)
    embedding <- as.matrix(pca_fit$x[, seq_len(max_comp), drop = FALSE])
  }
  if (is.null(rownames(embedding))) {
    rownames(embedding) <- rownames(dataset)
  }
  storage.mode(embedding) <- "double"
  return(embedding)
}

# Build a k-NN graph with Annoy -------------------------------------------------
.ramex_annoy_knn <- function(embedding, k, n_trees) {
  annoy <- RcppAnnoy::AnnoyEuclidean$new(ncol(embedding))
  for (i in seq_len(nrow(embedding))) {
    annoy$addItem(i - 1L, embedding[i, ])
  }
  annoy$build(n_trees)
  idx <- matrix(NA_integer_, nrow = nrow(embedding), ncol = k)
  for (i in seq_len(nrow(embedding))) {
    neighbors <- annoy$getNNsByItem(
      item = i - 1L,
      n = k + 1L
    )
    neighbors <- neighbors[neighbors != (i - 1L)]
    if (!length(neighbors)) {
      neighbors <- rep(i - 1L, k)
    } else if (length(neighbors) < k) {
      neighbors <- c(neighbors, rep(tail(neighbors, 1), k - length(neighbors)))
    }
    idx[i, ] <- neighbors[seq_len(k)] + 1L
  }
  return(idx)
}

# Construct the SNN matrix from neighbor indices --------------------------------
.ramex_build_snn_graph <- function(neighbor_idx, prune, cell_names) {
  snn <- ComputeSNN(nn_ranked = neighbor_idx, prune = prune)
  rownames(snn) <- cell_names
  colnames(snn) <- cell_names
  return(snn)
}

# Run Louvain with (optional) resolution parameter ------------------------------
.ramex_run_louvain <- function(graph, resolution, seed) {
  args <- list(graph = graph, weights = igraph::E(graph)$weight)
  louvain_formals <- names(formals(igraph::cluster_louvain))
  if ("resolution" %in% louvain_formals) {
    args$resolution <- resolution
  } else if ("resolution_parameter" %in% louvain_formals) {
    args$resolution_parameter <- resolution
  } else if (!isTRUE(all.equal(resolution, 1))) {
    warning(
      "The installed igraph version does not expose a Louvain resolution ",
      "parameter; results will match the default resolution."
    )
  }
  set.seed(seed + as.integer(round(resolution * 1000)))
  clusters <- do.call(igraph::cluster_louvain, args)
  ids <- igraph::membership(clusters)
  return(ids)
}

# Remove clusters smaller than the requested threshold --------------------------
.ramex_filter_clusters <- function(ids, threshold, n_cells) {
  labels <- as.character(ids)
  min_size <- if (threshold < 1) ceiling(threshold * n_cells) else ceiling(threshold)
  min_size <- max(1, min(n_cells, min_size))
  cluster_sizes <- table(labels)
  small <- names(cluster_sizes[cluster_sizes < min_size])
  if (length(small)) {
    labels[labels %in% small] <- NA_character_
  }
  return(labels)
}

# Format column names for each resolution ---------------------------------------
.ramex_resolution_label <- function(resolution) {
  paste0("res.", format(resolution, trim = TRUE, scientific = FALSE))
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
