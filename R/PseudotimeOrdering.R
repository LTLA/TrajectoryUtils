#' The PseudotimeOrdering class
#'
#' The PseudotimeOrdering class defines a two-dimensional object where rows represent cells and columns represent paths through a trajectory (i.e., \dQuote{lineages}).
#' It is expected to contain a numeric matrix of pseudotime orderings for each cell (row) in each path (column).
#' If a cell is on a path, it should have a valid pseudotime for the corresponding column; otherwise its entry should be set to \code{NA}.
#' Cells may lie on multiple paths if those paths span shared regions of the trajectory.
#'
#' @section Constructor:
#' \code{PseudotimeOrdering(pathStats, cellData=NULL, pathData=NULL, metadata=list()} will construct a PseudotimeOrdering object given:
#' \itemize{
#' \item \code{pathStats}, a (usually numeric) matrix-like object of pseudotime orderings as described above.
#' Alternatively, a list of such matrices can be supplied if multiple statistics are associated with each cell/path combination.
#' By convention, the first matrix in such a list should contain the pseudotime orderings.
#' \item \code{cellData}, a \linkS4class{DataFrame} of cell-level metadata.
#' This should have number of rows equal to the number of cells.
#' \item \code{pathData}, a \linkS4class{DataFrame} of path-level metadata.
#' This should have number of rows equal to the number of paths.
#' \item \code{metadata}, a list of any additional metadata to be stored in the object.
#' }
#' 
#' @section Getting/setting path statistics:
#' In the following code chunks, \code{x} is a PseudotimeOrdering object.
#' \describe{
#' \item{\code{pathStat(x, i=1L, withDimnames=TRUE)}:}{Returns a (usually numeric) matrix-like object containing some path statistics.
#' The default of \code{i=1L} will extract the first matrix of path statistics - by convention, this should contain the pseudotime orderings.
#' \code{i} may also be a string if the path statistics in \code{x} are named.
#' The dimnames of the output matrix are guaranteed to be the same as \code{dimnames(x)} if \code{withDimnames=TRUE}.}
#' \item{\code{pathStat(x, i=1L, withDimnames=TRUE) <- value}:}{Replaces the path statistics at \code{i} in \code{x} with the matrix \code{value}.
#' This should have the same dimensions as \code{x}, and if \code{withDimnames=TRUE}, it should also have the same dimnames.}
#' \item{\code{pathStats(x, withDimnames=TRUE)}:}{Returns a list of matrices containing path statistics.
#' The dimnames of each matrix are guaranteed to be the same as \code{dimnames(x)} if \code{withDimnames=TRUE}.}
#' \item{\code{pathStats(x, withDimnames=TRUE) <- value}:}{Replaces the path statistics in \code{x} with those in the list \code{value}.
#' Each entry of \code{value} should have the same dimensions as \code{x}.
#' The dimnames of each matrix should also be the same as \code{dimnames(x)} if \code{withDimnames=TRUE}.}
#' \item{\code{pathStatNames(x)}:}{Returns a character vector containing the names for each matrix of path statistics.}
#' \item{\code{pathStatNames(x) <- value}:}{Replaces the names of the path statistics with those in the character vector \code{value}.}
#' }
#'
#' @section Getting/setting path metadata:
#' In the following code chunks, \code{x} is a PseudotimeOrdering object.
#' \describe{
#' \item{\code{npaths(x)}:}{Returns an integer scalar containing the number of paths in \code{x}.
#' This is the same as \code{ncol(x)}.}
#' \item{\code{pathnames(x)}:}{Returns a character vector containing the names of paths in \code{x} (or \code{NULL}, if no names are available).
#' This is the same as \code{colnames(x)}.}
#' \item{\code{pathnames(x) <- value}:}{Replaces the path names in \code{x} with those in the character vector \code{value} (or \code{NULL}, to unname the paths).
#' This is the same as \code{colnames(x) <- value}.}
#' \item{\code{pathData(x, use.names=TRUE)}:}{Returns a DataFrame containing the path-level metadata of \code{x}.
#' This has the same number of rows as the number of columns in \code{x}.
#' Row names are guaranteed to be equal to \code{pathnames(x)} if \code{use.names=TRUE}.}
#' \item{\code{pathData(x) <- value}:}{Replaces the path-level metadata of \code{x} with a DataFrame \code{value} containing the same number of rows.}
#' }
#'
#' @section Getting/setting cell metadata:
#' In the following code chunks, \code{x} is a PseudotimeOrdering object.
#' \describe{
#' \item{\code{ncells(x)}:}{Returns an integer scalar containing the number of cells in \code{x}.
#' This is the same as \code{nrow(x)}.}
#' \item{\code{cellnames(x)}:}{Returns a character vector containing the names of cells in \code{x} (or \code{NULL}, if no names are available).
#' This is the same as \code{rownames(x)}.}
#' \item{\code{cellnames(x) <- value}:}{Replaces the cell names in \code{x} with those in the character vector \code{value} (or \code{NULL}, to unname the cells).
#' This is the same as \code{rownames(x) <- value}.}
#' \item{\code{cellData(x, use.names=TRUE)}:}{Returns a DataFrame containing the cell-level metadata of \code{x}.
#' This has the same number of rows as \code{x}.
#' Row names are guaranteed to be equal to \code{cellnames(x)} if \code{use.names=TRUE}.}
#' \item{\code{cellData(x) <- value}:}{Replaces the cell-level metadata of \code{x} with a DataFrame \code{value} containing the same number of rows.}
#' }
#'
#' @section Further operations:
#' In the following code chunks, \code{x} is a PseudotimeOrdering object.
#'
#' \code{x$name} and \code{x$name <- value} will get and set, respectively, the \code{name}d field of the \code{cellData}.
#' This is primarily provided for convenience.
#'
#' Subsetting operations (e.g., \code{x[i, j]}) and combining operations (\code{rbind(x, ...)}, \code{cbind(x, ...)}) will return the expected PseudotimeOrdering object.
#'
#' \code{metadata(x)} and \code{metadata(x) <- value} will get and set, respectively, the metadata of \code{x}.
#'
#' @section Comments on advanced usage:
#' The PseudotimeOrdering class is actually just a reskin of the widely-used \linkS4class{SummarizedExperiment} class.
#' We re-use the same underlying data structure and simply rename \code{row} and \code{col} to cell and path, respectively.
#' This means that any method that operates on a SummarizedExperiment can also - in theory - be applied to PseudotimeOrdering.
#'
#' We chose to do this reskinning to provide a clear conceptual break between the two classes.
#' The PseudotimeOrdering's dimensions do not follow the SummarizedExperiment's conventional \dQuote{samples as columns} philosophy, as each row instead represents a cell/sample.
#' Similarly, it is hard to argue that the paths are really interpretable as \dQuote{features} in any meaningful sense.
#' By reskinning, we hide the SummarizedExperiment implementation from the end-user and avoid any confusion with the interpretation of PseudotimeOrdering's dimensions.
#'
#' Of course, we could just transpose the inputs to get them to fit into a SummarizedExperiment.
#' However, the use of rows as cells is convenient as we often have many cells but few paths; it is easier to inspect the pseudotime ordering matrix with this orientation.
#' It also allows us to store the PseudotimeOrdering as a column in the \code{\link{colData}} of a SummarizedExperiment.
#' In this manner, datasets can be easily annotated with pseudotime orderings from trajectory reconstruction methods.
#'
#' @author Aaron Lun
#' @examples
#' # Make up a matrix of pseudotime orderings.
#' ncells <- 200
#' npaths <- 5
#' orderings <- matrix(rnorm(1000), ncells, npaths)
#'
#' # Default constructor:
#' (pto <- PseudotimeOrdering(orderings))
#' (pto <- PseudotimeOrdering(list(ordering=orderings)))
#' 
#' # Adding some per-cell metadata:
#' pto$cluster <- sample(LETTERS, ncells, replace=TRUE)
#' table(pto$cluster)
#'
#' # Adding some per-path metadata:
#' pathData(pto)$description <- c("EMT", "differentiatoin", "activation", "other", "?")
#' pathData(pto)
#'
#' # Subsetting and combining works fine:
#' rbind(pto, pto)
#' cbind(pto, pto)
#' pto[1:10,]
#' pto[,1:2]
#'
#' @aliases
#' PseudotimeOrdering
#' PseudotimeOrdering-class
#' show,PseudotimeOrdering-method
#'
#' ncells
#' npaths
#' pathStatNames
#' pathStatNames<-
#' pathStatNames<-,PseudotimeOrdering-method
#'
#' pathStat
#' pathStat<-
#' pathStat<-,PseudotimeOrdering-method
#'
#' pathStats
#' pathStats<-
#' pathStats<-,PseudotimeOrdering-method
#'
#' pathnames
#' pathnames<-
#' pathnames<-,PseudotimeOrdering-method
#'
#' pathData
#' pathData<-
#' pathData<-,PseudotimeOrdering-method
#' 
#' cellnames
#' cellnames<-
#' cellnames<-,PseudotimeOrdering-method
#'
#' cellData
#' cellData<-
#' cellData<-,PseudotimeOrdering-method
#' 
#' $,PseudotimeOrdering-method
#' $<-,PseudotimeOrdering-method
#'
#' @docType class
#' @name PseudotimeOrdering
NULL

#' @export
setClass("PseudotimeOrdering", "SummarizedExperiment")

#' @export
PseudotimeOrdering <- function(pathStats, cellData=NULL, pathData=NULL, metadata=list()){ 
    if (is.null(pathData)) pathData <- substitute()
    se <- SummarizedExperiment(pathStats, rowData=cellData, colData=pathData, metadata=metadata)
    as(se, "PseudotimeOrdering")
}

#' @export
setMethod("show", "PseudotimeOrdering", function(object) {
    cat("class:", class(object), "\n")
    cat("dim:", dim(object), "\n")

    expt <- names(metadata(object))
    if (is.null(expt)) {
        expt <- character(length(metadata(object)))
    }
    coolcat("metadata(%d): %s\n", expt)

    nms <- pathStatNames(object)
    if (is.null(nms)) {
        nms <- character(length(pathStats(object, withDimnames = FALSE)))
    }
    coolcat("pathStats(%d): %s\n", nms)

    cellnames <- cellnames(object)
    if (!is.null(cellnames)) {
        coolcat("cellnames(%d): %s\n", cellnames)
    } else {
        cat("cellnames: NULL\n")
    }
    coolcat("cellData names(%d): %s\n", names(cellData(object, use.names = FALSE)))

    pathnames <- pathnames(object)
    if (!is.null(pathnames)) {
        coolcat("pathnames(%d): %s\n", pathnames)
    } else {
        cat("pathnames: NULL\n")
    }
    coolcat("pathData names(%d): %s\n", names(pathData(object)))
})

#########################################

#' @export
ncells <- function(x) nrow(x)

#' @export
cellData <- function(x, use.names=TRUE) rowData(x, use.names=use.names)

#' @export
cellnames <- function(x) rownames(x)

#' @export
setGeneric("cellData<-", function(x, value) standardGeneric("cellData<-"))

#' @export
setReplaceMethod("cellData", "PseudotimeOrdering", function(x, value) {
    rowData(x) <- value
    x
})

#' @export
setGeneric("cellnames<-", function(x, value) standardGeneric("cellnames<-"))

#' @export
setReplaceMethod("cellnames", "PseudotimeOrdering", function(x, value) {
    rownames(x) <- value
    x
})

#########################################

#' @export
npaths <- function(x) ncol(x)

#' @export
pathData <- function(x, use.names=TRUE) colData(x, use.names=use.names)

#' @export
pathnames <- function(x) colnames(x)

#' @export
setGeneric("pathData<-", function(x, value) standardGeneric("pathData<-"))

#' @export
setReplaceMethod("pathData", "PseudotimeOrdering", function(x, value) {
    colData(x) <- value
    x
})

#' @export
setGeneric("pathnames<-", function(x, value) standardGeneric("pathnames<-"))

#' @export
setReplaceMethod("pathnames", "PseudotimeOrdering", function(x, value) {
    colnames(x) <- value
    x
})

#########################################

#' @export
pathStat <- function(x, i=1L, withDimnames=TRUE) assay(x, i, withDimnames=withDimnames)

#' @export
pathStats <- function(x, withDimnames=TRUE) assays(x, withDimnames=withDimnames)

#' @export
pathStatNames <- function(x) assayNames(x)

#' @export
setGeneric("pathStat<-", function(x, ..., value) standardGeneric("pathStat<-"))

#' @export
setReplaceMethod("pathStat", "PseudotimeOrdering", function(x, i, withDimnames=TRUE, ..., value) {
    assay(x, i, withDimnames=withDimnames, ...) <- value
    x
})

#' @export
setGeneric("pathStats<-", function(x, ..., value) standardGeneric("pathStats<-"))

#' @export
setReplaceMethod("pathStats", "PseudotimeOrdering", function(x, withDimnames=TRUE, ..., value) {
    assays(x, withDimnames=withDimnames, ...) <- value
    x
})

#' @export
setGeneric("pathStatNames<-", function(x, ..., value) standardGeneric("pathStatNames<-"))

#' @export
setReplaceMethod("pathStatNames", "PseudotimeOrdering", function(x, value) {
    assayNames(x) <- value
    x
})

#########################################

#' @export
setMethod("$", "PseudotimeOrdering", function(x, name) cellData(x)[[name]])

#' @export
setReplaceMethod("$", "PseudotimeOrdering", function(x, name, value) {
    cellData(x)[[name]] <- value
    x
})

