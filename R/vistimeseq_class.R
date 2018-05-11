###############################################################################
#' @title Checks validity of \code{vistimeseq} class.
#'
#' @description function for checking \code{vistimeseq} vaidity
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
#' @slot cluster.features A list storing results from timecourse (feature)
#' clustering.
#' @slot diff.expr A list storing results of differential expression analysis.
#'
#' @name vistimeseq
#' @rdname vistimeseq
#' @aliases vistimeseq-class
#' @exportClass vistimeseq
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
    diff.expr = "list"
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
    diff.expr = list()
  ),
  validity = check_vistimeseq
)


###############################################################################
#' @title vistimeseq object and constructors
#'
#' @description Initializes the vistimeseq object and populates
#' the time, replicate, and group slots.
#'
#' @param raw.data Raw input data with columns corresponding to samples
#' (observations) and rows to features
#' @param project Project name (string)
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
#' @return Returns a vistimeseq object with the raw data stored in
#' object@@raw.data, object@@sample.data object@@group, object@@replicate,
#' and object@@time are also initialized.
#'
#' @export
#'
#' @examples
#' raw <- matrix(runif(3000), ncol = 30)
#' time <- rep(rep(1:5, each = 3), 2)
#' replicate <- rep(1:3, 10)
#' group <- rep(1:2, each = 15)
#' test_vistimeseq <- vistimeseq(
#' raw.data = raw,
#' time = time,
#' replicate = replicate,
#' group = group)
#' test_vistimeseq
#'
vistimeseq <- function(
  raw.data,
  project = "'vistimeseq' time course project",
  sample.data = NULL,
  feature.data = NULL,
  time = NULL,
  time_column = NULL,
  replicate = NULL,
  replicate_column = NULL,
  group = NULL,
  group_column = NULL
) {
  nSamples <- ncol(raw.data)
  nFeatures <- nrow(raw.data)

  if (all(is.null(time), is.null(time_column))){
    stop("Either time or time_column must be specified")
  }
  if(all(!is.null(time), !is.numeric(time))){
    stop("if specified, \"time\" must be a numeric vector.")
  }
  if(all(!is.null(time), (length(time) != nSamples))){
    stop("Length of time is ", length(time),". Should be equal to the number",
         " of columns in raw.data, ", nSamples)
  }
  if(!is.null(time_column)){
    if (!time_column %in% colnames(sample.data)){
      stop("No time_column, ", time_column, " in data frame sample.data")
    }
    time <- suppressMessages(as.numeric(sample.data[, time_column]))
    if (all(is.na(time)))
      stop("Invalid time input, must be numeric.")
  }
  if (all(!is.null(replicate), (length(replicate) != nSamples))) {
    stop("Length of replicate is ", length(replicate),
         ". Should be equal to the number of columns in raw.data, ", nSamples)
  }
  if (!is.null(replicate_column)) {
    if (! replicate_column %in% colnames(sample.data)) {
      stop("No replicate_column, ", replicate_column,
           " in data frame sample.data")
    }
    replicate <- sample.data[, replicate_column]
  }
  if (all(!is.null(group), (length(group) != nSamples))) {
    stop("Length of group is ", length(group),
         ". Should be equal to the number of columns in raw.data, ", nSamples)
  }
  if (!is.null(group_column)) {
    if (!group_column %in% colnames(sample.data)){
      stop("No group_column, ", group_column, " in data frame sample.data")
    }
    group <- sample.data[, group_column]
  }
  if(is.null(replicate)){
    replicate <- rep("R1", nSamples)
  }
  if(is.null(group)){
    group <- rep("G1", nSamples)
  }
  if(is.null(colnames(raw.data))) {
    colnames(raw.data) <- paste0("S", seq_len(nSamples))
  }
  if(is.null(rownames(raw.data))) {
    rownames(raw.data) <- paste0("F", seq_len(nFeatures))
  }
  sample_data <- data.frame(
    sample = colnames(raw.data),
    row.names = colnames(raw.data),
    stringsAsFactors = FALSE)

  feature_data <- data.frame(
    feature = rownames(raw.data),
    row.names = rownames(raw.data),
    stringsAsFactors = FALSE)

  if(is.null(time)) {
    time <- numeric()
  }
  object <- new(
    Class = "vistimeseq",
    project.name = project,
    raw.data = raw.data,
    data = raw.data,
    sample.names = colnames(raw.data),
    feature.names = rownames(raw.data),
    sample.data  = sample_data,
    feature.data = feature_data,
    time = time,
    replicate = replicate,
    group = group
  )
  if (!is.null(sample.data)) {
    object <- add_sample_data(object = object, sampledata = sample.data)
  }
  if (!is.null(feature.data)) {
    object <- add_feature_data(object = object, featuredata = feature.data)
  }
  return(object)
}


#' @title vistimeseq object and constructor from ExpressionSet
#'
#' @description Initializes the vistimeseq object from ExpressionSet
#' and populates the time, replicate, and group slots.
#'
#' @param eset ExpressionSet object
#' @param project Project name (string)
#' @param time_column A character string equal to one of the column names
#' of \code{sample.data}. This is an alternative, way of specifying the time
#' corresponding to samples, if \code{time} argument is not provided.
#' @param replicate_column A character string equal to one of the column names
#' of \code{sample.data}. This is an alternative, way of specifying the
#' replicate corresponding to samples, if \code{replicate} argument is not
#' provided. If both \code{replicate} and \code{replicate_column} are not set,
#' the function assumes all samples come from the same replicate, and assigns
#' a replicate name 'R1' to all samples.
#' @param group_column A character string equal to one of the column names
#' of \code{sample.data}. This is an alternative, way of specifying the
#' replicate corresponding to samples, if \code{replicate} argument is not
#' provided. If both \code{group} and \code{group_column} are not set, the
#' function assumes all samples come from the same group, and assigns a group
#' name 'G1' to all samples.
#'
#' @return Returns a vistimeseq object with the raw data stored in
#' object@@raw.data, object@@sample.data object@@group, object@@replicate,
#' and object@@time are also initialized.
#'
#' @importFrom Biobase fData pData exprs
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr rename
#' @export
#'
#' @examples
#' library(Biobase)
#' cop1_eset <- readRDS(
#' file = system.file("extdata", "NGS1471_esetCounts_fltr.rds",
#'   package = "vistimeseq", mustWork = TRUE)
#' )
#' # convert timpoint data T0, T2.5, T4, ..., T13 to numeric values
#' pData(cop1_eset)$time <- gsub("T", "", pData(cop1_eset)$time)
#' pData(cop1_eset)$time <- as.numeric(pData(cop1_eset)$time)
#' cop1_vistimeseq <- vistimeseqFromExpressionSet(
#'   cop1_eset, time = "time", group = "genotype", replicate = "individual")
#' cop1_vistimeseq
#'
vistimeseqFromExpressionSet <- function(
  eset,
  time_column,
  replicate_column = NULL,
  group_column = NULL,
  project = "'vistimeseq' time course project"
) {
  gene_data <- fData(eset) %>%
    rownames_to_column("feature")
  smp_data <- pData(eset) %>%
    rownames_to_column("sample")

  if(!is.null(replicate_column)) {
    if(replicate_column %in% colnames(smp_data)){
      smp_data[["replicate"]] <- smp_data[[replicate_column]]
    } else {
      stop("\"replicate_column\" not found.")
    }
  } else {
    smp_data[["replicate"]] <- rep("R1", nrow(smp_data))
  }

  if(!is.null(group_column)) {
    if(group_column %in% colnames(smp_data)){
      smp_data[["group"]] <- smp_data[[group_column]]
    } else {
      stop("\"group_column\" not found.")
    }
  } else {
    smp_data[["group"]] <- rep("G1", nrow(smp_data))
  }

  if(!is.numeric(smp_data[[time_column]])) {
    stop("time column must be numeric.")
  }
  smp_data[["time"]] <- smp_data[[time_column]]

  cnts <- exprs(eset)[, smp_data$sample]

  object <- vistimeseq(
    project = project,
    raw.data = cnts,
    feature.data = gene_data,
    sample.data = smp_data,
    time_column = "time",
    replicate_column = "replicate",
    group_column = "group"
  )
  return(object)
}
