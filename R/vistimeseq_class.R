###############################################################################
#' @title Checks validity of \code{vistimeseq} class.
#'
#' @param object a data object
#'
#' @return \code{TRUE} if valid and a character vector with errors.
#'
check_vistimeseq <- function(object){
  errors <- character()
  ncol_raw <- ncol(object@raw.data)
  nrow_raw <- nrow(object@raw.data)
  nrow_sample_data <- nrow(object@sample.data)
  nrow_feature_data <- nrow(object@feature.data)
  if (nrow_sample_data > 0) {
    if(ncol_raw != nrow_sample_data) {
      msg <- paste0("Number of rows in sample.data is ", nrow_sample_data,
                    ". Should be equal to number of columns in raw.data, ",
                    ncol_raw)
      errors <- c(errors, msg)
    }
    if (! all(rownames(object@sample.data) == colnames(object@raw.data))) {
      msg <- "Inconsistent row names in sample.data and column names in raw.data"
      errors <- c(errors, msg)
    }
    if (all(!is.null(object@data),
            !all(rownames(object@sample.data) == colnames(object@data)))) {
      msg <- "Inconsistent row names in sample.data and column names in data"
      errors <- c(errors, msg)
    }
  }
  if (nrow_feature_data > 0) {
    if (nrow_raw != nrow_feature_data) {
      msg <- paste0("Number of rows in feature.data is ", nrow_feature_data,
                    ". Should be equal to number of rows in raw.data, ",
                    nrow_raw
      )
      errors <- c(errors, msg)
    }
    if (!all(rownames(object@feature.data) == rownames(object@raw.data))) {
      msg <- "Inconsistent row names in feature.data and row names in raw.data"
      errors <- c(errors, msg)
    }
    if (all(!is.null(object@data),
            !all(rownames(object@feature.data) == rownames(object@data)))) {
      msg <- "Inconsistent row names in feature.data and row names in data"
      errors <- c(errors, msg)
    }
  }
  if (!is.null(object@data.collapsed)) {
    if (!all(rownames(object@sample.data.collapsed) ==
             colnames(object@data.collapsed))) {
      msg <- paste("Inconsistent row names in sample.data.collapsed",
                   "and column names in data.collapsed")
      errors <- c(errors, msg)
    }
  }
  length_group <- length(object@group)
  length_replicate <- length(object@replicate)
  length_time <- length(object@time)
  if (all(length_group > 0, length_group != ncol_raw)) {
    msg <- paste0("Group is length ", length_group,
                  ". Should be equal to number of rows in raw.data, ", ncol_raw)
    errors <- c(errors, msg)
  }
  if (all(length_replicate > 0, length_replicate != ncol_raw)) {
    msg <- paste0("Replicate is length ", length_replicate,
                  ". Should be equal to number of rows in raw.data, ", ncol_raw)
    errors <- c(errors, msg)
  }
  if (all(length_time > 0, length_time != ncol_raw)) {
    msg <- paste0("Time is length ", length_time,
                  ". Should be equal to number of rows in raw.data, ", ncol_raw)
    errors <- c(errors, msg)
  }
  if (length(errors) == 0) TRUE else errors
}


###############################################################################
#' @title The Time Course Data Class
#'
#' @description The \code{vistimeseq} object is the main object in the time course
#' experiment analysis. It stores all relevant information associated with
#' the dataset, including the raw data, group, replicate and time associated
#' with each sample (column of the data). The object includes also slots for
#' results from some class-specific methods.
#'
#' Each vistimeseq object has a number of key slots listed below.
#'
#' @slot project.name A character string indicating the project name.
#' @slot raw.data Raw data with columns corresponding to samples
#' (observations) and rows to features.
#' @slot norm.data Normalized data.
#' @slot sample.data A \code{data.frame} object were rows are samples
#' (observations) and columns are sample attributes (e.g. group/condition,
#' replicate, time)
#' @slot sample.data.collapsed A \code{data.frame} with rows corresponding
#' to samples aggregated over replicates. The columns indicate group membership
#' and time.
#' @slot feature.data A \code{data.frame} object were rows are features and
#' columns are feature names under different conventions (e.g. Enterez IDs,
#' gene symbols) or feature attributes (e.g. chromosome, location in the
#' genome, biotype, gc content etc.)
#' @slot group A character vector of length \code{ncol(raw.data)} indicating
#' the group membership of the sample
#' @slot replicate Vector of length \code{ncol(raw.data)} indicating
#' the replicate from which the sample came from
#' @slot time Vector of length \code{ncol(raw.data)} indicating the time
#' at which the sample was collected
#' @slot timecourse.data List of data in a time-course format. Each element
#' of the list is a \code{data.frame} where the data is organized in with
#' first three columns indicating feature, group, replicate, and the remaining
#' ones providing the data at each available time point.
#' \code{vistimeseq} methods will typically generate elements named:
#' 'tc', 'tc_with_lags', 'tc_collapsed' and 'tc_collapsed_with_lags'.
#' @slot dim.red List of stored dimmensional reductions; named by technique
#' @slot cluster.features A list storing results of timecourse feature clustering.
#' @slot diff.expr A \code{data.frame} storing results of differential
#' expression analysis
#'
#' @name vistimeseq
#' @rdname vistimeseq
#' @aliases vistimeseq-class
#' @exportClass vistimeseq
#' @useDynLib vistimeseq
#'
vistimeseq <- setClass(
  "vistimeseq",
   slots = c(
    project.name = "character",
    raw.data = "ANY",
    data = "ANY",
    sample.names = "character",
    feature.names = "character",
    sample.data = "data.frame",
    feature.data = "data.frame",
    data.collapsed = "ANY",
    sample.data.collapsed = "data.frame",
    group = "ANY",
    replicate = "ANY",
    time = "numeric",
    timecourse.data = "list",
    dim.red = "list",
    cluster.features = "list",
    diff.expr = "data.frame"
  ),
  prototype = list(
    project.name = character(),
    raw.data = NULL,
    data = NULL,
    sample.names = character(),
    feature.names = character(),
    sample.data = data.frame(),
    feature.data = data.frame(),
    data.collapsed = NULL,
    sample.data.collapsed = data.frame(),
    group = character(),
    replicate = character(),
    time = numeric(),
    timecourse.data = list(),
    dim.red = list(),
    cluster.features = list(),
    diff.expr = data.frame()
  ),
  validity = check_vistimeseq
)

###############################################################################
#' show method for vistimeseq
#'
#' @param object A vistimeseq object
#' @name show
#' @aliases show,vistimeseq-method
#' @docType methods
#' @rdname show-methods
#'
setMethod(
  f = "show",
  signature = "vistimeseq",
  definition = function(object) {
    cat(
      "An object of class \"",
      class(object),
      "\". \nProject name: ",
      object@project.name,
      " \nRaw data: ",
      dim(object@raw.data)[1],
      " genes across ",
      dim(object@raw.data)[2],
      " samples.\n",
      sep = ""
    )
    invisible(x = NULL)
  }
)

