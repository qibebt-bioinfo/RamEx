#' Plot Heatmap of Markers
#'
#' @param object A Ramanome object
#' @param markers A list of markers (wave numbers)
#' @param group group information to show the heatmap and marker intensities
#' @param scale Whether to scale the data (default: "row")
#' @param show_legend Whether to show legend (default: TRUE)
#' @param show_rownames Whether to show row names (Wave Number) (default: FALSE)
#' @param show_colnames Whether to show column names (Group levels) (default: TRUE)
#' @param wave_gap_threshold Threshold for wave number gap to split wave number into groups (default: 50)
#' @return A heatmap plot
#' @export
Plot.Heatmap.Markers <- function(object, markers, group=object$group,
                            scale = "row", show_legend = TRUE,
                            show_rownames = FALSE, show_colnames = TRUE,
                            wave_gap_threshold = 20) {
    data_after <- Feature.Reduction.Intensity(object, as.list(markers))
    markers_intensity <- data.frame(
        group = data_after@meta.data$group,
        data_after@interested.bands
    )
    markers_intensity <- aggregate(
        markers_intensity[,-1],
        list(markers_intensity$group),
        mean
    )
    row.names(markers_intensity) <- markers_intensity$Group.1
    
    waves <- as.numeric(gsub('X', '', colnames(markers_intensity[,-1])))
    wave_order <- order(waves)
    waves_sorted <- waves[wave_order]
    
    wave_diff <- diff(waves_sorted)
    split_points <- which(wave_diff > wave_gap_threshold)

    groups <- cut(1:length(waves_sorted), c(0, split_points, length(waves_sorted)))
    
    group_ranges <- tapply(waves_sorted, groups, function(x) {
        paste0(round(min(x), 0), "~", round(max(x), 0))
    })
    
    bands_ann <- data.frame(
        bands = factor(groups, levels = levels(groups), labels = group_ranges)[order(wave_order)]
    )
    row.names(bands_ann) <- colnames(markers_intensity[,-1])[wave_order]
    
    n_groups <- length(unique(bands_ann$bands))
    annotation_colors <- RamEx.colors[1:n_groups]
    names(annotation_colors) <- unique(bands_ann$bands)
    annotation_colors <- list(bands = annotation_colors)
    
    colors <- colorRampPalette(c("#F7F650FF", '#51C56AFF', "#2D708EFF"))(100)
    p <- pheatmap(
        t(markers_intensity[,-1][,wave_order]),
        color = colors,
        cluster_cols = FALSE,
        cluster_rows = TRUE,
        legend = show_legend,
        scale = scale,
        trace = "none",
        angle_col = "45",
        show_rownames = show_rownames,
        show_colnames = show_colnames,
        clustering_distance_rows = "correlation",
        treeheight_row = 0,
        density.info = "none",
        annotation_row = bands_ann,
        annotation_colors = annotation_colors
    )
    return(p)
}

#' Mark the Ramanome Markers on Average Spectrum
#' Mark the Ramanome Markers on Average Spectrum
#'
#' @param object A Ramanome object
#' @param markers_single A list of single markers (wave numbers)
#' @param markers_paired A list of paired markers (wave numbers) or NULL/NA
#' @param gap Vertical gap between spectra (default: 0.6)
#' @param cols Color palette for markers (default: c('#F18656','#397FC7'))
#' @return A ggplot object
#' @export
Plot.Markers.Spectrum <- function(object, markers, gap = 0.6, cols = c('#F18656','#397FC7')) {
    mean_spectrum <- data.frame(
        wave = object@wavenumber,
        mean = colMeans(get.nearest.dataset(object))
    )
    
    markers_single <- markers$correlations_singular$wave
    markers_paired <- markers$correlations_paired

    mean_spectrum$markers_single <- mean_spectrum$mean
    mean_spectrum$markers_single[!mean_spectrum$wave %in% markers_single] <- NA
    
    if (!is.null(markers_paired) && !all(is.na(markers_paired))) {
        markers_paired <- as.vector(as.matrix(markers_paired))
        
        mean_spectrum$markers_paired <- mean_spectrum$mean
        mean_spectrum$markers_paired[!mean_spectrum$wave %in% markers_paired] <- NA
        plot_data <- reshape2::melt(mean_spectrum, id.vars = c("wave", "mean"))
        plot_data$variable <- factor(plot_data$variable, 
                                   levels = c("markers_single", "markers_paired"))

        plot_data$mean <- plot_data$mean + gap * (2 - as.numeric(plot_data$variable))
        plot_data$value <- plot_data$value + gap * (2 - as.numeric(plot_data$variable))
    } else {
        plot_data <- data.frame(
            wave = mean_spectrum$wave,
            mean = mean_spectrum$mean,
            value = mean_spectrum$markers_single,
            variable = "markers_single"
        )
    }
    
    p <- ggplot(plot_data, aes(x = wave, y = mean, group = variable)) +
        geom_line(linewidth = 0.8, color = "#999999") +
        geom_point(aes(wave, value, color = variable, fill = variable),
                  alpha = 0.5, size = 1, na.rm = TRUE) +  # 添加 na.rm = TRUE
        scale_fill_manual(values = cols[1:length(unique(plot_data$variable))]) +
        scale_color_manual(values = cols[1:length(unique(plot_data$variable))]) +
        scale_x_continuous(breaks = seq(500, 3000, by = 500)) +
        theme_classic()
    
    return(p)
}

#' Plot Violin and Box Plots
#'
#' @param data A data frame or matrix of data
#' @param group A vector or data frame of group information
#' @param cols Optional color palette for groups
#' @param violin_width Width of violin plot (default: 0.75)
#' @param box_width Width of box plot (default: 0.12)
#' @return A ggplot object
#' @export
Plot.ViolinBox <- function(data, group, cols = NULL, violin_width = 0.75, box_width = 0.12) {
    data <- assign_variable(data, group)
    
    if(is.null(cols)){
        cols <- RamEx.colors
    }
    
    if(length(cols) < length(unique(data$variable))){
        warning("The number of colors is not equal to the number of groups")
        cols <- cols[c(1:length(unique(data$variable) -1 )) %% length(cols) +1]
    }
    
    p <- ggplot(data, aes(x = group, y = value, fill = variable)) +
        geom_violin(trim = FALSE, alpha = 0.5,
                   position = position_dodge(width = violin_width)) +
        geom_boxplot(width = box_width,
                    outlier.shape = NA,
                    position = position_dodge(width = violin_width)) +
        scale_fill_manual(values = cols) +
        theme_classic()
    
    return(p)
}

#' @noRd 
assign_variable <- function(data, group){
    if(is.null(dim(group)) || (is.null(ncol(group)) && length(group) == 1)){
        if(is.null(dim(data)) || (is.null(ncol(data)) && length(data) == 1)){
            data <- data.frame(value = data, group = group, variable = group)
        } else {  
            data_df <- data.frame(data, group = group)
            n_vars <- ncol(data_df) - 1  
            data <- data.frame(
                value = unlist(data_df[, -ncol(data_df)]),
                group = rep(data_df$group, n_vars),
                variable = rep(names(data_df)[-ncol(data_df)], each = nrow(data_df))
            )
        }
    } else if(!is.null(ncol(group)) && ncol(group) == 2){
        if(is.null(dim(data)) || (is.null(ncol(data)) && length(data) == 1)){
            data <- data.frame(value = data, group = group[,1], variable = group[,2])
        } else {
            warning("data and group both have more than one column, only the first column in data will be used")
            data <- data.frame(value = data[,1], group = group[,1], variable = group[,2])
        }
    } else {
        warning("group must have at most two columns, only the first two columns will be used")
        data <- data.frame(value = data, group = group[,1], variable = group[,2])
    }
    return(data)
}

#' Plot Scatter Plot from Reductions
#'
#' @param object A Ramanome object
#' @param reduction Name of the reduction to plot (default: "UMAP")
#' @param cols Optional color palette for groups
#' @param point_size Size of points (default: 1)
#' @param point_alpha Transparency of points (default: 0.5)
#' @return A ggplot object
#' @export
Plot.reductions <- function(object, reduction = "UMAP", dims=c(1,2), group = object$group, cols=RamEx.colors,
                        point_size = 1, point_alpha = 0.5) {
  plot_data <- object@reductions[[reduction]][,dims]
  plot_data$group <- group
  names <- colnames(plot_data)[1:2]
  
  p <- ggplot(plot_data, aes(
    x = get(names[1]),
    y = get(names[2]),
    color = group,
    fill = group
  )) +
    geom_point(size = point_size, alpha = point_alpha) +
    theme_classic()
  
  return(p)
}


#' Plot histogram of group distributions
#'
#' @param group_1 A vector of group assignments
#' @param group_2 A vector of group assignments
#' @param cols Optional color palette for groups
#' @param width Width of the alluvium (default: 0.9)
#' @param alpha Transparency of the alluvium (default: 0.5)
#' @param curve_type Type of curve for the alluvium (default: "sine")
#' @return A ggplot object
#' @export
Plot.Cluster.Distribution <- function(group_1, group_2, cols = RamEx.colors, 
                                    width = 0.9, alpha = 0.5, 
                                    curve_type = "sine") {
    # 计算比例
    plot_data <- data.frame(
        table(group_1, group_2) / 
        rowSums(table(group_1, group_2))
    )
    colnames(plot_data) <- c('Group_1', 'Group_2', 'Prop')
    
    if (is.null(cols)) {
        cols <- RamEx.colors
    }
    
    p <- ggplot(plot_data, aes(
        x = Group_1,
        y = Prop,
        fill = Group_2,
        stratum = Group_2,
        alluvium = Group_2
    )) +
        geom_stratum(width = width, color = 'white') +
        geom_alluvium(
            alpha = alpha,
            width = width,
            curve_type = curve_type
        ) +
        scale_fill_manual(values = cols) +
        scale_color_manual(values = cols) +
        theme_classic() +
        labs(
            x = "Time",
            y = "Proportion",
            fill = "Cluster"
        )
    
    return(p)
}



#' export RamEx.colors
RamEx.colors <- c('#64B5F6FF','#FF7080FF','#6BD76BFF','#8888FFFF', 
                  '#F18656','#397FC7','#999999', '#B41BBC',
                  '#1790a3', '#DF6182', '#f4c61f', '#B07836',
                  '#2b6584', '#A1CC44', '#F88E1B', '#1bf8c8')