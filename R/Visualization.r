#' Plot Heatmap of Markers
#'
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
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
                            scale = "row", show_legend = TRUE, cluster_cols = FALSE, cluster_rows = TRUE,
                            show_rownames = FALSE, show_colnames = TRUE, show = TRUE, save = FALSE,
                            wave_gap_threshold = 20, width = 10, height = 10) {
    data_after <- Feature.Reduction.Intensity(object, as.list(markers))
    markers_intensity <- data.frame(
        group = group,
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
    if(save){cat('Saving heatmap to the current working directory: ', getwd(), '\n')}
    if(length(split_points) < 3){
        p <- pheatmap(
        t(markers_intensity[,-1]),
        cluster_cols = FALSE,
        cluster_rows = TRUE,
        legend = show_legend,
        scale = scale,
        clustering_distance_rows = "correlation",
        silent = TRUE)
        wave_order <- p$tree_row$order
        waves_sorted <- waves[wave_order]

        group_ranges <- merge_wave_numbers(waves_sorted, wave_gap_threshold)
        bands_ann <- data.frame(
            bands = factor(group_ranges$assignments, levels = unique(group_ranges$assignments), labels = group_ranges$intervals)
        )
    } else {
        groups <- cut(1:length(waves_sorted), c(0, split_points, length(waves_sorted)))
        group_ranges <- tapply(waves_sorted, groups, function(x) {
            paste0(round(min(x), 0), "~", round(max(x), 0))
        })
        bands_ann <- data.frame(
        bands = factor(groups, levels = levels(groups), labels = group_ranges)[order(wave_order)]
    )
    }
    
    row.names(bands_ann) <- colnames(markers_intensity[,-1])[wave_order]
    
    n_groups <- length(unique(bands_ann$bands))
    annotation_colors <- RamEx.color(n_groups)
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
        silent = !show,
        filename = ifelse(save, 'Raman_Markers_heatmap.png', NA),
        width = width,
        height = height
    )
    return(p)
}


#' Mark the Ramanome Markers on Average Spectrum
#' Mark the Ramanome Markers on Average Spectrum
#'
#' @param object A Ramanome object
#' @param markers list object obtained from 'Raman markers' module.
#' @param gap Vertical gap between spectra (default: 0.6)
#' @param cols Color palette for markers (default: c('#F18656','#397FC7'))
#' @return A ggplot object
#' @export
Plot.Markers.Spectrum <- function(object, markers, gap = 0.6, cols = c('#F18656','#397FC7')) {
    mean_spectrum <- data.frame(
        wave = object@wavenumber,
        mean = colMeans(get.nearest.dataset(object))
    )
    
    markers_singular <- markers$markers_singular$wave
    markers_paired <- markers$markers_paired

    mean_spectrum$markers_singular <- mean_spectrum$mean
    mean_spectrum$markers_singular[!mean_spectrum$wave %in% markers_singular] <- NA
    
    if (!is.null(markers_paired) && !all(is.na(markers_paired))) {
        markers_paired <- as.vector(as.matrix(markers_paired))
        
        mean_spectrum$markers_paired <- mean_spectrum$mean
        mean_spectrum$markers_paired[!mean_spectrum$wave %in% markers_paired] <- NA
        plot_data <- reshape2::melt(mean_spectrum, id.vars = c("wave", "mean"))
        plot_data$variable <- factor(plot_data$variable, 
                                   levels = c("markers_singular", "markers_paired"))

        plot_data$mean <- plot_data$mean + gap * (2 - as.numeric(plot_data$variable))
        plot_data$value <- plot_data$value + gap * (2 - as.numeric(plot_data$variable))
    } else {
        plot_data <- data.frame(
            wave = mean_spectrum$wave,
            mean = mean_spectrum$mean,
            value = mean_spectrum$markers_singular,
            variable = "markers_singular"
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
        cols <- RamEx.color(length(unique(data$variable)))      
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
#' @param reduction Name of the reduction to plot (default: "UMAP"). This must be one of 'PCA', 'UMAP', 'PCoA', 'tSNE'
#' @param dims Dimensions to plot (default: c(1,2))
#' @param color Group information or intensity values to fill the points
#' @param cols Optional color palette for drawing.
#' @param point_size Size of points (default: 1)
#' @param point_alpha Transparency of points (default: 0.5)
#' @param quantile_range Range of quantiles for color scaling when given a continuous values as 'color'(default: c(0.05, 0.95))
#' @return A ggplot object
#' @export
Plot.reductions <- function(object, reduction = "UMAP", dims = c(1,2), color = object$group, 
                        cols = NULL, point_size = 1, point_alpha = 0.5,
                        quantile_range = c(0.05, 0.95)) {
    plot_data <- object@reductions[[reduction]][,dims]
    plot_data$group <- color
    names <- colnames(plot_data)[1:2]
    
    if (is.numeric(color) && length(unique(color)) > 10) {
        if (is.null(cols)) {
            cols <- colorRampPalette(c("#F7F650FF", '#51C56AFF', "#2D708EFF"))(100)
        }
        
        p <- ggplot(plot_data, aes(
            x = get(names[1]),
            y = get(names[2]),
            color = color
        )) +
            geom_point(size = point_size, alpha = point_alpha) +
            scale_color_gradientn(
                colors = cols,
                limits = quantile(color, quantile_range),
                oob = scales::squish
            ) +
            labs(x = names[1], y = names[2],color = 'Group') +
            theme_classic()
    } else {
        n_groups <- length(unique(color))
        if (is.null(cols)) {
            cols <- RamEx.color(n_groups)
        }
        
        p <- ggplot(plot_data, aes(
            x = get(names[1]),
            y = get(names[2]),
            color = color
        )) +
            geom_point(size = point_size, alpha = point_alpha) +
            scale_color_manual(values = cols) +
            labs(x = names[1], y = names[2],color = 'Group') +
            theme_classic()
    }
    
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
Plot.Distribution <- function(group_1, group_2, cols = NULL, width = 0.9) {
    plot_data <- data.frame(
        table(group_1, group_2) / 
        rowSums(table(group_1, group_2))
    )
    colnames(plot_data) <- c('Group_1', 'Group_2', 'Prop')
    
    if (is.null(cols)) {
        cols <- RamEx.color(length(unique(plot_data$Group_2)))
    }
    
    p <- ggplot(plot_data, aes(
        x = Group_1,
        y = Prop,
        fill = Group_2
    )) +
        geom_bar(stat = "identity", position = "stack", width = width) +
        scale_fill_manual(values = cols) +
        theme_classic() +
        labs(
            x = "Group_1",
            y = "Proportion",
            fill = "Group_2"
        )
    
    return(p)
}

#' Plot scatter plot
#'
#' @param x A vector of x values
#' @param y A vector of y values
#' @param cols Optional color palette for points
#' @param point_size Size of points (default: 1)
#' @param point_alpha Transparency of points (default: 0.5)
#' @return A ggplot object
#' @export
Plot.scatter <- function(x, y, cols = NULL, point_size = 1, point_alpha = 0.5) {
    plot_data <- data.frame( x, y)
    p <- ggplot(plot_data, aes(x = x, y = y, color = cols, fill = cols)) +
        geom_point(size = point_size, alpha = point_alpha) +
        scale_color_continuous(low = '#bebef1', high = '#0449a9') +
        scale_fill_continuous(low = '#bebef1', high = '#0449a9') +
        theme_classic()
    return(p)
}

#' Get RamEx colors
#'
#' @param n Number of colors to return
#' @return A vector of colors
#' @export
RamEx.color <- function(n) {
    colors <- c('#64B5F6FF','#FF7080FF','#6BD76BFF','#8888FFFF', 
                  '#F18656','#397FC7','#999999', '#B41BBC',
                  '#1790a3', '#DF6182', '#f4c61f', '#B07836',
                  '#2b6584', '#A1CC44', '#F88E1B', '#1bf8c8')
    if (n <= 16) {
        return(colors[1:n])
    } else {
        return(colorRampPalette(c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", 
                                "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", 
                                "#bcbd22", "#17becf"))(n))
    }
}


#' Merge wave numbers into groups based on the gap threshold
#' @noRd 
merge_wave_numbers <- function(wave_numbers, gap_threshold = 20) {
    groups <- list()  
    current_group <- c(wave_numbers[1])  

    for (i in 2:length(wave_numbers)) {  
        if (abs(wave_numbers[i] - tail(current_group, 1)) <= gap_threshold) {  
            current_group <- c(current_group, wave_numbers[i])  
        } else {  
            groups <- append(groups, list(current_group))  
            current_group <- c(wave_numbers[i])  
        }  
    }  
    groups <- append(groups, list(current_group))  

    merged_groups <- merge_groups(groups, gap_threshold)  

    group_intervals <- sapply(merged_groups, function(group) {  
      if (length(group) > 1) {  
        paste(round(min(group), 0), round(max(group), 0), sep = "~")  
      } else {  
        as.character(group)  
      }  
    })  

    group_assignments <- rep(NA, length(wave_numbers))  
    group_number <- 1    

    for (i in seq_along(merged_groups)) {  
      for (value in merged_groups[[i]]) {  
        group_assignments[which(wave_numbers == value)[1]] <- group_number  
      }  
      group_number <- group_number + 1  
    }  
    return(list(assignments = group_assignments, intervals = group_intervals))
}

#' Merge groups based on the gap threshold
#' @noRd 
merge_groups <- function(groups, gap_threshold) {  
  merged <- list(groups[[1]])  
  
  for (i in 2:length(groups)) {  
    current_group <- groups[[i]]  
    merged_last <- tail(merged, 1)[[1]]  
    
    if (any(abs(outer(merged_last, current_group, `-`)) <= gap_threshold)) {  
      merged[[length(merged)]] <- c(merged_last, current_group)  
    } else {  
      merged <- append(merged, list(current_group))  
    }  
  }  
  
  further_merge_needed <- FALSE  
  for (i in 2:length(merged)) {  
    prev_group <- merged[[i - 1]]  
    curr_group <- merged[[i]]  
    if (any(abs(outer(prev_group, curr_group, `-`)) <= gap_threshold)) {  
      further_merge_needed <- TRUE  
      break  
    }  
  }  
  
  if (further_merge_needed) {  
    return(merge_groups(merged, gap_threshold))  
  } else {  
    return(merged)  
  }  
}  