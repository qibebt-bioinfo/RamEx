#' Ramanome class
#'
#' This class represents a Ramanome object, which contains Raman data and associated metadata.
#'
#' @slot wavenumber Numeric vector indicating the wavenumbers of the Raman spectra
#' @slot datasets A list of Raman datasets
#' @slot meta.data A data frame containing metadata associated with the Raman datasets
#' @slot reductions A list of data reduction results
#' @slot interested.bands A list of interested Raman bands
#'
#' @return A new Ramanome object
#' @examples
#' # Create sample data
#' wavenumbers <- seq(500, 3500, by = 10)
#' spectra <- matrix(rnorm(100 * length(wavenumbers)), nrow = 100)
#' metadata <- data.frame(
#'   group = factor(rep(c("Control", "Treatment"), each = 50)),
#'   filenames = paste0("sample_", 1:100, ".txt")
#' )
#' 
#' # Create a Ramanome object
#' raman_obj <- new("Ramanome",
#'   datasets = list(raw = spectra),
#'   wavenumber = wavenumbers,
#'   meta.data = metadata
#' )
#'
#' @importFrom methods setClass
#' @importFrom methods setMethod
#' @export

Ramanome <- setClass(
  Class = 'Ramanome',
  slots = c(
    wavenumber = "numeric",
    datasets = 'list',
    meta.data = 'data.frame',
    reductions = 'list',
    interested.bands = 'list'
  )
)

#' Access slots of Ramanome object using $ operator
#'
#' This method allows accessing slots of a Ramanome object using the $ operator.
#' It can access datasets, meta.data columns, reductions, and interested.bands.
#'
#' @param x A Ramanome object   
#' @param name The name of the slot to access
#' @return The requested slot content
#' @export
setMethod("$", "Ramanome", function(x, name) {
  if (name %in% names(x@datasets)) {
    return(x@datasets[[name]])
  } else if (name %in% colnames(x@meta.data)) {
    return(x@meta.data[[name]])
  } else if (name %in% names(x@reductions)) {
    return(x@reductions[[name]])
  } else if (name %in% names(x@interested.bands)) {
    return(x@interested.bands[[name]])
  } else {
    stop(paste("No slot named", name, "found in Ramanome object"))
  }
})

#' Get available slot names for Ramanome object
#'
#' This method returns all available slot names that can be accessed using the $ operator.
#' It enables tab-completion in RStudio and other IDEs.
#'
#' @param x A Ramanome object
#' @return A character vector of available slot names
#' @export
setMethod("names", "Ramanome", function(x) {
  c(
    names(x@datasets),
    colnames(x@meta.data),
    names(x@reductions),
    names(x@interested.bands)
  )
})

#' Show method for Ramanome objects
#' This method will give brief information about the Ramanome object
#' @param object The Ramanome object
#' @return No return value, prints object information to console
#' @examples
#' # Create sample Ramanome object
#' wavenumbers <- seq(500, 3500, by = 10)
#' spectra <- matrix(rnorm(100 * length(wavenumbers)), nrow = 100)
#' metadata <- data.frame(
#'   group = factor(rep(c("Control", "Treatment"), each = 50))
#' )
#' raman_obj <- new("Ramanome",
#'   datasets = list(raw = spectra),
#'   wavenumber = wavenumbers,
#'   meta.data = metadata
#' )
#' 
#' # Display object information
#' show(raman_obj)
#' @export
setMethod(
  f = "show",
  signature = "Ramanome",
  definition = function(object) {
    num.spec <- nrow(get.nearest.dataset(object))
    num.wave <- length(object@wavenumber)
    num.group <- length(unique(object@meta.data$group))
    cat(class(x = object), "Object", "\n")
    cat(
      num.wave,
      'Raman features across',
      num.spec,
      'samples within',
      num.group,
      ifelse(test = num.group == 1, yes = 'group', no = 'groups'),
      "\n"
    )
    cat('\n')
  }
)

#' Subset method for Ramanome objects
#'
#' @param x The Ramanome object
#' @param i The indices for subsetting
#' @param j Index or name of the slot to subset
#' @param ... Additional arguments
#' @param drop Boolean indicating whether to drop unused slots (default: TRUE)
#' @return A new Ramanome object containing the subset of data
#' @examples
#' # Create sample Ramanome object
#' wavenumbers <- seq(500, 3500, by = 10)
#' spectra <- matrix(rnorm(100 * length(wavenumbers)), nrow = 100)
#' metadata <- data.frame(
#'   group = factor(rep(c("Control", "Treatment"), each = 50))
#' )
#' raman_obj <- new("Ramanome",
#'   datasets = list(raw = spectra),
#'   wavenumber = wavenumbers,
#'   meta.data = metadata
#' )
#' 
#' @export
setMethod("[", "Ramanome", function(x, i, j, ..., drop = TRUE) {
  new.wavenumber = x@wavenumber
  new.datasets = lapply(x@datasets,'[',i,)
  new.meta.data = x@meta.data[i,]
  new.reductions = lapply(x@reductions, function(reduction) {
    if (is.matrix(reduction)) {
      reduction[i, , drop = FALSE]
    } else {
      reduction
    }
  })
  new.interested.bands = lapply(x@interested.bands, function(bands) {
    if (is.matrix(bands)) {
      bands[i, , drop = FALSE]
    } else {
      bands
    }
  })
  new_obj <- new("Ramanome", wavenumber = new.wavenumber, datasets = new.datasets, meta.data=new.meta.data,
                 reductions=new.reductions, interested.bands=new.interested.bands )
  return(new_obj)
})

#' Length method for Ramanome objects
#'
#' @param x The Ramanome object
#' @return An integer indicating the number of samples in the Ramanome object
#' @export
setMethod("length", "Ramanome", function(x) {
  nrow(get.nearest.dataset(x))
})



#' Combine Ramanome objects by row binding
#'
#' @param x The first Ramanome object
#' @param y The second Ramanome object
#' @return A new Ramanome object containing the combined data from both input objects
#' @export

setMethod("rbind2", signature(x = "Ramanome", y = "Ramanome"),
          function(x, y) {
            Ramanome(wavenumber = x@wavenumber ,
                      datasets = list(data = rbind(get.nearest.dataset(x), get.nearest.dataset(y))), 
                      meta.data = rbind(x@meta.data, y@meta.data))
          })

#' Plot mean spectra for Ramanome objects
#'
#' Mean spectra of each group in a given Ramanome
#' The latest spectral matrix in datasets slot and the group information in meta.data slot are used for plotting
#'
#' @param x A Ramanome object
#' @param y Not used (required for plot generic)
#' @param ... Additional arguments passed to mean.spec()
#' @return A plot showing mean spectra and the diversity of each group
#' @export
setMethod("plot", "Ramanome", function(x, y, ...) {
  # 获取光谱数据和分组信息
  spectra <- get.nearest.dataset(x)
  groups <- x@meta.data$group
  
  # 使用 mean.spec 绘制平均光谱
  mean.spec(spectra, groups, ...)
})