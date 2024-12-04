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
#' @export


setMethod("rbind2", signature(x = "Ramanome", y = "Ramanome"),
          function(x, y) {
            Ramanome(wavenumber = x@wavenumber ,datasets = list(data = rbind(get.nearest.dataset(x),
                                                                             get.nearest.dataset(y))), meta.data = rbind(x@meta.data, y@meta.data))
          })


