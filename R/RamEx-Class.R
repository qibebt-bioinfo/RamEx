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
#'   meta.data = metadata,
#'   reductions = list(),
#'   interested.bands = list()
#' )
#'
#' @importFrom hyperSpec nrow
#' @importFrom hyperSpec print
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
  new.reductions = lapply(x@reductions,'[',i,)
  new.interested.bands = lapply(x@interested.bands,'[',i,)
  new_obj <- new("Ramanome", wavenumber = new.wavenumber, datasets = new.datasets, meta.data=new.meta.data,
                 reductions=new.reductions, interested.bands=new.interested.bands )
  return(new_obj)
})

#' Length method for Ramanome objects
#'
#' @param x The Ramanome object
#' @return An integer indicating the number of samples in the Ramanome object
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
#' # Get number of samples
#' n_samples <- length(raman_obj)
#' print(n_samples)
#' @export

setMethod("length", "Ramanome", function(x) {
  nrow(get.nearest.dataset(x))
})

#' #' Plot method for Ramanome objects
#' #'
#' #' @param x The Ramanome object
#' #' @param y Character specifying the plot type
#' #' @importFrom hyperSpec plot
#' #' @importFrom hyperSpec print
#' #' @importFrom methods setMethod
#'
#'
#' setMethod("plot",
#'           signature(x = "Ramanome", y = "character"),
#'           function(x, y, ...) {
#'             tmp <- hyperSpec::plot(x, y, ...)
#'             if (is(tmp, "trellis"))
#'               hyperSpec::print(tmp)
#'             invisible(tmp)
#'           }
#' )

#' Combine Ramanome objects by row binding
#'
#' @param x The first Ramanome object
#' @param y The second Ramanome object
#' @return A new Ramanome object containing the combined data from both input objects
#' @examples
#' # Create two sample Ramanome objects
#' wavenumbers <- seq(500, 3500, by = 10)
#' spectra1 <- matrix(rnorm(50 * length(wavenumbers)), nrow = 50)
#' spectra2 <- matrix(rnorm(50 * length(wavenumbers)), nrow = 50)
#' metadata1 <- data.frame(
#'   group = factor(rep("Control", 50))
#' )
#' metadata2 <- data.frame(
#'   group = factor(rep("Treatment", 50))
#' )
#' 
#' raman_obj1 <- new("Ramanome",
#'   datasets = list(raw = spectra1),
#'   wavenumber = wavenumbers,
#'   meta.data = metadata1
#' )
#' 
#' raman_obj2 <- new("Ramanome",
#'   datasets = list(raw = spectra2),
#'   wavenumber = wavenumbers,
#'   meta.data = metadata2
#' )
#' 
#' # Combine the two objects
#' combined_obj <- rbind2(raman_obj1, raman_obj2)
#' @export

setMethod("rbind2", signature(x = "Ramanome", y = "Ramanome"),
          function(x, y) {
            Ramanome(wavenumber = x@wavenumber ,datasets = list(data = rbind(get.nearest.dataset(x),
                                                                             get.nearest.dataset(y))), meta.data = rbind(x@meta.data, y@meta.data))
          })
