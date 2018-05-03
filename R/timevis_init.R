###############################################################################
#' @title Initialize and setup the timevis object
#'
#' @description Initializes the timevis object and populates
#' the time, replicate, and group slots.
#'
#' @param project Project name (string)
#' @param raw.data Raw input data with columns corresponding to samples
#' (observations) and rows to features
#' @param sample.data Optional. A \code{data.frame} object were rows are
#' samples (observations) and columns are sample attributes
#' (e.g. group/condition, replicate, time)
#' @param feature.data Optional. A \code{data.frame} object were rows are
#' features and columns are feature names under different conventions (e.g.
#' Enterez IDs,  gene symbols) or feature attributes (e.g. chromosome,
#' location in the genome, biotype, gc content etc.)
#' @param time Vector of length \code{ncol(raw.data)} indicating the time
#' at which the sample was collected. If not given, \code{time_column} MUST
#' be specified.
#' @param time_column A character string equal to one of the column names
#' of \code{sample.data}. This is an alternative, way of specifying the time
#' corresponding to samples, if \code{time} argument is not provided.
#' @param replicate Vector of length \code{ncol(raw.data)} indicating
#' the replicate from which the sample came from. If not given,
#' \code{replicate_column} will be used.
#' @param replicate_column A character string equal to one of the column names
#' of \code{sample.data}. This is an alternative, way of specifying the
#' replicate corresponding to samples, if \code{replicate} argument is not
#' provided. If both \code{replicate} and \code{replicate_column} are not set,
#' the function assumes all samples come from the same replicate, and assigns
#' a replicate name 'R1' to all samples.
#' @param group A character vector of length \code{ncol(raw.data)} indicating
#' the group membership of the sample. If not given, \code{group_column}
#' will be used.
#' @param group_column A character string equal to one of the column names
#' of \code{sample.data}. This is an alternative, way of specifying the
#' replicate corresponding to samples, if \code{replicate} argument is not
#' provided. If both \code{group} and \code{group_column} are not set, the
#' function assumes all samples come from the same group, and assigns a group
#' name 'G1' to all samples.
#'
#' @return Returns a timevis object with the raw data stored in object@@raw.data.
#' object@@group, object@@replicate, and object@@timealso initialized.
#'
#' @export
#'
#' @examples
#' endoderm_raw <- read.table(
#'   file = system.file('extdata', 'raw_endoderm_small.txt',
#'     package = 'timevis'),
#'   as.is = TRUE
#' )
#' time <- as.numeric(gsub("(.*)\\_D", "", colnames(endoderm_raw)))
#' replicate <- gsub("\\_(.*)", "", colnames(endoderm_raw))
#' group <- substring(replicate, 1, 1)
#' endoderm_small <- timevis(
#'   raw.data = endoderm_raw,
#'   time = time,
#'   replicate = replicate,
#'   group = group)
#' endoderm_small
#'
timevis_init <- function(
  raw.data,
  project = "timevis time course project",
  sample.data = NULL,
  feature.data = NULL,
  time = NULL,
  time_column = NULL,
  replicate = NULL,
  replicate_column = FALSE,
  group = NULL,
  group_column = NULL
) {
  nSamples <- ncol(raw.data)
  nFeatures <- nrow(raw.data)
  object <- new(
    Class = "timevis",
    project.name = project,
    raw.data = raw.data
  )
  object@sample.data <- data.frame(sample = colnames(raw.data),
                                   row.names = colnames(raw.data),
                                   stringsAsFactors = FALSE)
  if (!is.null(sample.data)) {
    object <- add_sample_data(object = object, sampledata = sample.data)
  }
  object@feature.data <- data.frame(feature = rownames(raw.data),
                                    row.names = rownames(raw.data),
                                    stringsAsFactors = FALSE)
  if (!is.null(feature.data)) {
    object <- add_feature_data(object = object, featuredata = feature.data)
  }
  if(!is.null(time)) {
    if (!is.numeric(time))
      stop("The argument time, if specified must be a numeric vector.")
    if (length(time) != nSamples)
      stop(paste0("Length of time is ", length(time),
                  ". Should be equal to the number",
                  " of columns in raw.data, ", nSamples))
    object@time <- time
  } else {
    if (is.null(time_column))
      stop("Either time or time_column must be specified")
    if (! time_column %in% colnames(sample.data))
      stop("No time_column, ", time_column, " in data frame sample.data")
    time <- suppressMessages(as.numeric(sample.data[, time_column]))
    if (all(is.na(time)))
      stop("Invalid time input, must be numeric.")
    object@time <- time
  }
  if (!is.null(replicate)) {
    if (length(replicate) != nSamples)
      stop(paste0("Length of replicate is ", length(replicate),
                  ". Should be equal to the number",
                  " of columns in raw.data, ", nSamples))
    object@replicate <- replicate
  } else {
    if (!is.null(replicate_column)) {
      if (! replicate_column %in% colnames(sample.data))
        stop("No replicate_column, ", replicate_column,
             " in data frame sample.data")
      object@replicate <- sample.data[, replicate_column]
    } else {
      object@replicate <- rep("R1", nSamples)
    }
  }
  if (!is.null(group)) {
    if (length(group) != nSamples)
      stop(paste0("Length of group is ", length(group),
                  ". Should be equal to the number",
                  " of columns in raw.data, ", nSamples))
    object@group <- group
  } else {
    if (!is.null(group_column)) {
      if (! group_column %in% colnames(sample.data))
        stop("No group_column, ", group_column, " in data frame sample.data")
      object@group <- sample.data[, group_column]
    } else {
      object@group <- rep("G1", nSamples)
    }
  }
  return(object)
}
