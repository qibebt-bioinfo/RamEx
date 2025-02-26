#' Calculate Mahalanobis distance for every observation
#'
#' This function calculates the Mahalanobis distance for every observation in the given data matrix.
#' @importFrom proxy as.matrix
#' @param data The input data matrix.
#'
#' @return A vector containing the Mahalanobis distance for each observation.


mahalDistForEvery <- function(data) {
  data <- as.matrix(data)
  inv_cov <- solve(cov(data))
  x_mu <- colMeans(data)
  I <- dim(data)[1]
  distance_matrix <- matrix(rep(NA, I), nrow = I)
  for (i in 1:I) {
    x_minus_y <- data[i,] - x_mu
    mdist <- t(x_minus_y) %*% inv_cov %*% x_minus_y
    distance_matrix[i] <- mdist
  }
  return(sqrt(distance_matrix))
}

#' Detect Outliers Using Mahalanobis Distance and Chi-squared Test
#'
#' This function detects outliers in a dataset using the Mahalanobis distance and chi-squared test.
#' Outliers are defined as observations with Mahalanobis distances that have chi-squared test p-values below a certain threshold.
#'
#' @param data A numeric matrix or data frame.
#' @param p Threshold probability (default is 0.1).
#' @return A vector containing the indices of the detected outliers.

outliers_maha_chisquare <- function(data, p = 0.1) {
  del <- which(colMeans(data) < 0.02 | colMeans(data) == 1)
  if (length(del) != 0) { data <- data[, -del] }
  outliers <- NULL
  for (i in 1:floor(ncol(data) / 2)) {
    maha_matrix <- mahalDistForEvery(data[, (2 * i - 1):(2 * i)])
    chi_p <- stats::dchisq(maha_matrix, df = 2)
    outliers <- c(outliers, which(chi_p < p))
  }
  if (ncol(data) %% 2 != 0) {
    maha_matrix <- mahalDistForEvery(data[, (ncol(data) - 1):(ncol(data))])
    chi_p <- stats::dchisq(maha_matrix, df = 2)
    outliers <- c(outliers, which(chi_p < p))
  }
  return(unique(outliers))
}


#' Calculate Pearson's and Spearman's Correlation Coefficients for a Data Frame
#'
#' This function calculates Pearson's and Spearman's correlation coefficients for a given data frame.
#' The function uses the `rcorr` function from the Hmisc package to perform the correlation analysis.
#' @importFrom Hmisc rcorr
#' @importFrom proxy as.matrix
#' @param df A data frame.
#' @return A list containing the correlation matrix, p-values, and the number of observations used for each correlation.

rcorr_df <- function(df) {
  df_mat <- as.matrix(df)
  return(rcorr(df_mat))
}


#' Vectorize Distance Matrix
#'
#' This function vectorizes a distance matrix, optionally grouping the samples by a metadata category.
#'
#' @param dm A distance matrix.
#' @param group A vector specifying the sample categories. (default is NULL)
#' @param duplicate A logical value indicating whether to include duplicate pairs. (default is TRUE)
#' @importFrom reshape2 melt
#' @return A data frame containing the vectorized distance matrix.

vectorize_dm <- function(dm, group = NULL, duplicate = TRUE) {
  if (ncol(dm) != nrow(dm))
    stop('The distance matrix is not squared')
  dm <- data.matrix(dm)
  if (!is.null(group)) {
    if (length(unique(group)) < 2)
      stop('At least two levels for a given sample category are required.')
    if (length(group) != nrow(dm))
      stop('The number of rows in the metadata and distance matrix is not equal')
    if (is.factor(group) & nlevels(group) > length(group) * 0.9)
      stop('The number of levels in a certain category cannot exceed 90% of the total number of samples')

    colnames(dm) <- rownames(dm) <- paste(rownames(dm), group, sep = "____")
    if (duplicate) {
      melt_dm <- subset(melt(dm))
      melt_dm$value <- as.numeric(melt_dm$value)
    } else {
      dm[lower.tri(dm)] <- NA
      diag(dm) <- NA
      melt_dm <- subset(melt(dm), !is.na(value))
      melt_dm$value <- as.numeric(melt_dm$value)
    }

    Row_Info <- data.frame(do.call(rbind, strsplit(as.character(melt_dm[, 1]), "____", fixed = TRUE)))
    Col_Info <- data.frame(do.call(rbind, strsplit(as.character(melt_dm[, 2]), "____", fixed = TRUE)))
    VS <- paste0(Row_Info[, 2], "_VS._", Col_Info[, 2])
    dm_value <- data.frame(VS, Row_Info, Col_Info, d = melt_dm$value)

    colnames(dm_value) <- c("GroupPair", "Sample_1", "Group_1", "Sample_2", "Group_2", "value")
    if (is.factor(group)) {
      DistType <- as.factor(dm_value$Group_1 == dm_value$Group_2)
      DistType <- factor(
        DistType,
        levels = levels(DistType),
        labels = c("AllBetween", "AllWithin")
      )
      dm_value <- data.frame(DistType, dm_value)
      dm_value$value <- as.numeric(dm_value$value)
    }
  } else {
    if (duplicate) {
      dm_value <- subset(melt(dm))
      dm_value$value <- as.numeric(dm_value$value)
    } else {
      dm[lower.tri(dm)] <- NA
      diag(dm) <- NA
      dm_value <- subset(melt(dm), !is.na(value))
      dm_value$value <- as.numeric(dm_value$value)
    }
    colnames(dm_value) <- c("Sample_1", "Sample_2", "value")
  }

  return(dm_value)
}


#' Vectorize Correlation Coefficients and p-values from a Correlation Matrix
#'
#' This function vectorizes the correlation coefficients and p-values from a correlation matrix,
#' optionally grouping the samples by a metadata category.
#'
#' @param dm A list containing correlation coefficients (r) and p-values (P).
#' @param group A vector specifying the sample categories. (default is NULL)
#' @param duplicate A logical value indicating whether to include duplicate pairs. (default is TRUE)
#' @return A data frame containing the vectorized correlation coefficients and p-values.

vectorize_dm_rcorr <- function(dm, group = NULL, duplicate = TRUE) {
  dm_r <- data.matrix(dm$r)
  dm_pvalue <- data.matrix(dm$P)
  dm_r_value <- vectorize_dm(dm_r, group = group, duplicate = FALSE)
  dm_p_value <- vectorize_dm(dm_pvalue, group = group, duplicate = FALSE)
  dm_value <- cbind(dm_r_value, dm_p_value[, "value"])
  colnames(dm_value)[ncol(dm_value)] <- "p_value"
  return(dm_value)
}


#' Plot Global Chord Diagram with Statistical Significance
#'
#' This function plots a global chord diagram based on a correlation matrix, highlighting edges with statistical significance.
#'
#' @param Corr_mat The correlation matrix for plotting the chord diagram.
#' @param a_degree A dataframe containing the degree information.
#' @param a_SumCorr A dataframe containing the sum correlation information.
#' @param a_MeanCorr A dataframe containing the mean correlation information.
#' @param main_title The main title for the plot.
#' @param Pos_Edge Logical indicating whether to include positive edges.
#' @param Neg_Edge Logical indicating whether to include negative edges.
#' @param Threshold The threshold value for considering edges as statistically significant.
#' @param bands_ann A dataframe containing the band annotation information.
#' @return No explicit return value. The function generates the plot directly.
#' @importFrom stringr str_replace_all
#' @import circlize
#' @import circlize

Plot_global_chordDiagram_rpvalue <- function(
    Corr_mat,
    a_degree = NULL,
    a_SumCorr = NULL,
    a_MeanCorr = NULL,
    main_title = NULL,
    Pos_Edge = FALSE,
    Neg_Edge = FALSE,
    Threshold = 0.6,
    bands_ann = NULL
) {
  w_num <- sort(unique(as.numeric(str_replace_all(rownames(Corr_mat), "[A-Z]", ""))))
  if (length(unique(w_num)) < length(w_num))
    stop("Row/column names of input correlation matrix should be unique. Please check!")

  bands_ord <- data.frame(Band = rownames(Corr_mat), ord = order(w_num))
  v_Corr_mat <- vectorize_dm(Corr_mat, group = bands_ord$ord, duplicate = FALSE)

  if (Pos_Edge & Neg_Edge) {
    v_Corr_mat <- subset(v_Corr_mat, value > Threshold | value < (-Threshold))
  } else if (Pos_Edge) {
    v_Corr_mat <- subset(v_Corr_mat, value > Threshold & value < 1)
  } else if (Neg_Edge) {
    v_Corr_mat <- subset(v_Corr_mat, value > -1 & value < (-Threshold))
  } else stop('Please check the pos/neg edges parameters')

  v_Corr_mat$Group_1 <- as.numeric(str_replace_all(v_Corr_mat$Sample_1, "[A-Z]", ""))
  v_Corr_mat$Group_2 <- as.numeric(str_replace_all(v_Corr_mat$Sample_2, "[A-Z]", ""))

  par(mar = c(1, 1, 1, 1))
  circos.clear()
  circos.par("start.degree" = 90)
  raw_wn <- as.numeric(str_replace_all(colnames(Corr_mat), "[A-Z]", ""))
  circos.initialize("a", xlim = c(raw_wn[1] - 1, raw_wn[length(raw_wn)] + 0.1 * (raw_wn[length(raw_wn)] - raw_wn[1]))) # 'a' just means there is one sector

  if (Pos_Edge & Neg_Edge) {
    col_line <- "grey20"
  } else if (Pos_Edge) {
    col_line <- "#FF6666"
  } else if (Neg_Edge) {
    col_line <- "#66CCFF"
  } else stop('Please check the pos/neg edges parameters')

  col_mat <- colorRamp2(
    c(-1, -Threshold, -Threshold + 1e-10, Threshold - 1e-10, Threshold, 1),
    c("#0033FF", "#66CCFF", "#FFFFFF00", "#FFFFFF00", "#FF6666", "#FF3300")
  )

  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = 0.3,
    bg.border = NA,
    panel.fun = function(x, y) {
      xlim <- get.cell.meta.data("xlim")

      if (!is.null(a_degree)) {
        circos.trackLines(
          a_degree$factor,
          a_degree$x,
          a_degree$y,
          type = "h",
          area = TRUE,
          col = col_line,
          border = col_line
        )
        circos.text(
          xlim[1] - 10,
          0.3,
          expression('Degree'),
          facing = "downward",
          adj = c(1, 1),
          col = col_line,
          cex = 1
        )
      }

      if (!is.null(a_SumCorr)) {
        circos.trackLines(
          a_SumCorr$factor,
          a_SumCorr$x,
          a_SumCorr$y,
          type = "h",
          area = TRUE,
          col = col_line,
          border = col_line
        )
        circos.text(
          xlim[1] - 10,
          0.3,
          expression('SumCorr'),
          facing = "downward",
          adj = c(1, 1),
          col = col_line,
          cex = 1
        )
      }

      if (!is.null(a_MeanCorr)) {
        circos.trackLines(
          a_MeanCorr$factor,
          a_MeanCorr$x,
          a_MeanCorr$y,
          type = "h",
          area = TRUE,
          col = col_line,
          border = col_line
        )
        circos.text(
          xlim[1] - 10,
          0.3,
          expression('MeanCorr'),
          facing = "downward",
          adj = c(1, 1),
          col = col_line,
          cex = 1
        )
      }

      circos.lines(
        c(raw_wn[1], raw_wn[length(raw_wn)]),
        c(0, 0), col = "#CCCCCC"
      )

      circos.text(
        xlim[1] - 10,
        0.1,
        expression('Shift (cm'^-1*')'),
        facing = "downward",
        adj = c(1, 1),
        col = "grey60",
        cex = 0.8
      )

      breaks <- round(c(raw_wn[1], c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1) * raw_wn[length(raw_wn)]), 0)
      circos.axis(
        h = 0,
        major.at = breaks,
        labels = breaks,
        col = "grey50",
        direction = "inside",
        labels.facing = "inside",
        major.tick.length = 0.02,
        labels.col = "grey50",
        labels.cex = 0.7
      )

      if (nrow(v_Corr_mat) > 1) {
        apply(v_Corr_mat, 1, function(i) {
          circos.link("a", as.numeric(i[3]), "a", as.numeric(i[5]), col = col_mat(as.numeric(i[6])))
        })
      }
    }
  )

  circos.clear()
}


#' Plot local chord diagram with p-value thresholds
#'
#' This function plots a local chord diagram based on a correlation matrix and p-value thresholds.
#'
#' @param x A correlation matrix.
#' @param Pos_Edge Logical value to indicate whether to include positive edges. Default is FALSE.
#' @param Neg_Edge Logical value to indicate whether to include negative edges. Default is FALSE.
#' @param Threshold The threshold value for correlation. Default is 0.6.
#' @param bands_ann Annotation for the bands.
#'
#' @return A chord diagram plot showing correlations between variables. The plot includes:
#'   \item{Correlation links}{Colored edges showing correlations between variables, with colors indicating correlation strength}
#'   \item{Band annotations}{Labels and groupings for the variables}
#'   \item{Sector highlights}{Visual grouping of related variables}
#'   The function has no explicit return value as it generates the plot directly.
#'
#' @import RColorBrewer
#' @import circlize
#' @importFrom graphics par
#' @importFrom graphics plot

Plot_local_chordDiagram_rpvalue <- function(
    x,
    Pos_Edge = FALSE,
    Neg_Edge = FALSE,
    Threshold = 0.6,
    bands_ann = NULL
) {
  # Compute correlation matrix and p-values
  Corr_mat <- rcorr_df(x)[[1]] # r
  Corr_mat_pvalue <- rcorr_df(x)[[3]] # p.value

  # Assign bands labels to rows and columns of the correlation matrix
  Bands_label <- bands_ann$Wave_num
  colnames(Corr_mat) <- rownames(Corr_mat) <- Bands_label
  colnames(Corr_mat_pvalue) <- rownames(Corr_mat_pvalue) <- Bands_label

  # Set non-significant correlations to a value greater than 1
  Corr_mat[which(seq_along(Corr_mat) %in% which(Corr_mat_pvalue < 0.05) == FALSE)] <- 1.1

  # Define color palette
  colours <- c(
    brewer.pal(8, "Dark2"),
    brewer.pal(9, "Set1"),
    brewer.pal(9, "Pastel1"),
    brewer.pal(12, "Paired"),
    brewer.pal(8, "Set2"),
    brewer.pal(12, "Set3"),
    brewer.pal(8, "Accent"),
    brewer.pal(8, "Pastel2")
  )

  # Determine color scale depending on the edge type
  if (Pos_Edge & Neg_Edge) {
    col_mat <- colorRamp2(
      c(-1, 0, 1),
      c("Blue4", "White", "Red2")
    )
    Corr_mat[which(Corr_mat < Threshold & Corr_mat > -(Threshold))] <- 1.1 # r cutoff
  } else if (Pos_Edge) {
    col_mat <- colorRamp2(
      c(-1, -1e-10 + Threshold, Threshold, 1),
      c("#FFFFFF00", "#FFFFFF00", "Red", "DarkRed")
    )
    Corr_mat[which(Corr_mat < Threshold)] <- 1.1 # r cutoff
  } else if (Neg_Edge) {
    col_mat <- colorRamp2(
      c(-1, -Threshold, -Threshold + 1e-10, 1),
      c("DarkBlue", "Blue", "#FFFFFF00", "#FFFFFF00")
    )
    Corr_mat[which(Corr_mat > (-Threshold))] <- 1.1 # r cutoff
  } else {
    stop('Please check the pos/neg edges parameters')
  }

  col_link <- col_mat(Corr_mat)
  col_link[which(Corr_mat == 1.1)] <- "#FFFFFF00"

  # Create grid of bands annotations
  bands_ann <- bands_ann[order(bands_ann$Group),]
  bands_ann$Group <- factor(bands_ann$Group)
  grid.col <- colours[unclass = bands_ann$Group]

  # Set plot parameters
  op <- par(mar = c(0, 0, 0, 0))

  # Compute circle gap widths
  n_bands_Group <- as.numeric(lapply(levels(bands_ann$Group), function(x) { nrow(bands_ann[which(bands_ann$Group == x),]) })) - 1
  x_vec <- NULL
  for (i in n_bands_Group) {
    tem <- rep(1, i)
    x_vec <- c(x_vec, tem, 8)
  }

  # Generate chord diagram plot
  circos.clear()
  circos.par(gap.degree = x_vec) # circle gap widths

  chordDiagram(
    Corr_mat,
    order = as.character(bands_ann$Wave_num),
    symmetric = TRUE,
    link.visible = TRUE,
    keep.diagonal = TRUE,
    transparency = 0,
    col = col_link,
    grid.col = grid.col,
    grid.border = NULL,
    annotationTrack = "grid",
    preAllocateTracks = list(list(track.height = 0.02))
  )

  # Add band labels
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.index = 2,
    panel.fun = function(x, y) {
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      sector.index <- get.cell.meta.data("sector.index")
      circos.text(
        mean(xlim),
        mean(ylim),
        sector.index,
        col = "White",
        font = 2,
        cex = 0.4,
        facing = "bending.inside",
        niceFacing = TRUE
      )
    },
    bg.border = NA
  )

  # Highlight sectors based on band groups
  invisible(lapply(levels(bands_ann$Group), function(x) {
    b <- as.character(subset(bands_ann, Group == x)$Wave_num)
    highlight.sector(
      sector.index = b,
      track.index = 1,
      col = "grey80",
      lwd = 0.001,
      text = x,
      text.col = "grey30",
      text.vjust = -0.5,
      facing = "bending.inside",
      niceFacing = TRUE,
      cex = 0.6
    )
  }))

  # Clear circos and restore plot parameters
  circos.clear()
  par(op)
}

#' Perform IRCA analysis on global scale
#'
#' This function performs IRCA analysis on the global scale using the input dataset and group information.
#' It generates a chord diagram that represents the significant correlations between variables.
#'
#' @param dataset The input dataset for analysis
#' @param group The group information for the dataset
#'
#' @return A dataframe containing the significant correlations between variables
#' @importFrom Hmisc rcorr

Intraramanome.Analysis.Irca.Global.draw <- function(dataset, group) {
  outliers <- outliers_maha_chisquare(dataset)
  if (length(outliers) != 0) {
    dataset <- dataset[-outliers,]
  }
  cor_list <- rcorr(dataset, type = 'pearson')
  corr_matrix <- cor_list[[1]]
  corr_matrix[cor_list[[3]] > 0.05] <- 0
  rownames(corr_matrix) <- colnames(corr_matrix)
  Plot_global_chordDiagram_rpvalue(
    corr_matrix,
    Neg_Edge = TRUE,
    Threshold = 0.6
  )
  corr_matrix[!upper.tri(corr_matrix, diag = TRUE)] <- 0
  locs <- which(corr_matrix < -0.6, arr.ind = TRUE)
  waves <- round(as.numeric(colnames(corr_matrix)), 0)
  interests <- data.frame(
    wave1 = waves[locs[, 1]],
    wave2 = waves[locs[, 2]],
    corr = corr_matrix[locs],
    group = group
  )
  return(interests)
}

Intraramanome.Analysis.Irca.Global.cal <- function(dataset, group, threshold=0.6) {
  outliers <- outliers_maha_chisquare(dataset)
  if (length(outliers) != 0) {
    dataset <- dataset[-outliers,]
  }
  cor_list <- rcorr(dataset, type = 'pearson')
  corr_matrix <- cor_list[[1]]
  corr_matrix[cor_list[[3]] > 0.05] <- 0
  rownames(corr_matrix) <- colnames(corr_matrix)
  corr_matrix[!upper.tri(corr_matrix, diag = TRUE)] <- 0
  locs <- which(corr_matrix < -0.6, arr.ind = TRUE)
  waves <- round(as.numeric(colnames(corr_matrix)), 0)
  interests <- data.frame(
    wave1 = waves[locs[, 1]],
    wave2 = waves[locs[, 2]],
    corr = corr_matrix[locs],
    group = group
  )
  return(interests)
}

#' Perform IRCA analysis on local scale
#'
#' This function performs IRCA analysis on the local scale using the input dataset and bands annotation information.
#' It generates a chord diagram that represents the significant correlations between variables within each band.
#'
#' @param dataset The input dataset for analysis
#' @param bands_ann The bands annotation information
#' @return A chord diagram plot showing the significant correlations between variables within each band
#' @examples
#' # Create sample dataset
#' n <- 100  # samples
#' p <- 10   # variables
#' dataset <- matrix(rnorm(n*p), nrow=n)
#' colnames(dataset) <- paste0("W", 1:p)
#'
#' # Create bands annotation
#' bands_ann <- data.frame(
#'   Wave_num = paste0("W", 1:p),
#'   Group = rep(c("A", "B"), each=p/2)
#' )
#'

Intraramanome.Analysis.Irca.Local.draw <- function(dataset, bands_ann) {
  outliers <- outliers_maha_chisquare(dataset)
  if (length(outliers) != 0) {
    dataset <- dataset[-outliers,]
  }
  Plot_local_chordDiagram_rpvalue(
    dataset,
    Neg_Edge = TRUE,
    Threshold = 0.6,
    bands_ann = bands_ann
  )
}

#' Perform global IRCA analysis
#'
#' This function performs global IRCA analysis on the given Raman spectroscopy data object.
#' It generates negative IRCA plots for each group and saves them as jpeg files.
#'
#' @param object The Raman spectroscopy data object.
#'
#' @return A matrix containing the global IRCA interests.
#'
#' @export Intraramanome.Analysis.Irca.Global
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocessing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocessing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocessing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocessing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- RamEx_data[qc_icod$quality,]
#' IRCA.interests <- Intraramanome.Analysis.Irca.Global(data_cleaned)

Intraramanome.Analysis.Irca.Global <- function(object) {
  dataset <- get.nearest.dataset(object)
  waves <- object@wavenumber
  IRCA.interests <- lapply(unique(object@meta.data$group), function(x) {
    temp_data <- dataset[object@meta.data$group == x,]
    jpeg(
      filename = paste0('IRCA_global_negative_', x, '.jpeg'),
      width = 800,
      height = 800,
      quality = 100,
      res = 200
    )
    interests <- Intraramanome.Analysis.Irca.Global.draw(temp_data, x)
    dev.off()
    return(interests)
  })
  IRCA.interests <- do.call(rbind, rlist::list.map(IRCA.interests, .))
  return(IRCA.interests)
}

#' Perform local IRCA analysis
#'
#' This function performs local IRCA analysis on the given Raman spectroscopy data object.
#' It generates negative IRCA plots for each group and saves them as png files.
#'
#' @importFrom Cairo Cairo
#' @param object The Raman spectroscopy data object.
#' @param bands_ann The band annotation data.
#' @return No return value, called for side effects. The function generates and saves PNG files
#'         containing negative IRCA plots for each group in the dataset.
#' @export Intraramanome.Analysis.Irca.Local
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocessing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocessing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocessing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocessing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- RamEx_data[qc_icod$quality,]
#' bands_ann <- data.frame(rbind(cbind(c(742,850,872,971,997,1098,1293,1328,1426,1576),'Nucleic acid'),cbind(c(824,883,1005,1033,1051,1237,1559,1651),'Protein'),cbind(c(1076,1119,1370,2834,2866,2912),'Lipids')))
#' colnames(bands_ann) <- c('Wave_num', 'Group')
#' Intraramanome.Analysis.Irca.Local(data_cleaned, bands_ann = bands_ann)

Intraramanome.Analysis.Irca.Local <- function(object, bands_ann) {
  dataset <- get.nearest.dataset(object)
  waves <- round(object@wavenumber, 0)
  locs <- unlist(lapply(as.numeric(bands_ann$Wave_num), function(x)which.min(abs(object@wavenumber - x))))
  bands_ann$Wave_num <- waves[locs]
  dataset <- dataset[, waves %in% bands_ann$Wave_num]
  lapply(unique(object@meta.data$group), function(x) {
    temp_data <- dataset[object@meta.data$group == x,]
    Cairo(
      file = paste0('IRCA_local_negative_', x, '.png'),
      units = "in",
      dpi = 300,
      width = 3,
      height = 3,
      type = 'png'
    )
    Intraramanome.Analysis.Irca.Local.draw(temp_data, bands_ann)
    grDevices::dev.off()
  })
}


#' 2D Correlation Spectroscopy Analysis (2D-COS)
#' Captures both synchronous (simultaneous changes) and
#' asynchronous (sequential changes) relationships, providing
#' detailed insights into spectral changes
#'
#' @param object A Ramanome object containing the dataset.
#' @return synchronous and asynchronous correlation spectra result
#' @export Intraramanome.Analysis.2Dcos
#' @import corr2D
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocessing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocessing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocessing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocessing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- RamEx_data[qc_icod$quality,]
#' data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
#' #reslut_2dos <- Intraramanome.Analysis.2Dcos(data_cleaned)
Intraramanome.Analysis.2Dcos <- function(object) {
  data <- object@datasets$normalized.data
  twod<-corr2d(data)
  plot_corr2d(twod, Legend = FALSE)
  return(twod)
}
