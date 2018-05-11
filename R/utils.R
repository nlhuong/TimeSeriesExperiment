#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @return None
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' @title Variance stabilization.
#'
#' @description This function performs variance stabilization
#' on assay data.
#'
#' @param X an assay data matrix or data.frame where columns correspond
#' @param method Method for variance stabilization (VST).
#' Currently, supports "none" (no VST), "log1p" (log plus one), "asinh"
#' (inverse hyperbolic sine) or "deseq"
#' (\code{\link[DESeq2]{varianceStabilizingTransformation}} function from
#' \code{DESeq2} package). Default is "log1p".
#'
#' @return Returns a varianced stabilized data matrix.
#'
#' @importFrom DESeq2 varianceStabilizingTransformation
#' @export
#' @examples
#' X <- sapply(exp(rlnorm(10)), function(m) rnbinom(20, size = 1, mu = m))
#' head(X)
#' Y <- variance_stabilization(X, method = "asinh")
#' head(Y)
#'
variance_stabilization <- function(X, method = "log1p") {
  X <- as.matrix(X)
  if(method == "none") {
    Y <- X
  } else if (method == "log1p"){
    Y <- log1p(X)
  } else if (method == "asinh"){
    Y <- asinh(X)
  } else if (method == "deseq"){
    Y <- varianceStabilizingTransformation(X)
  } else {
    stop("Unsupported variance stabilization method.")
  }
  return(Y)
}


#' @title Melt data matrix.
#'
#' @description Convert data matrix to a long format
#'
#' @param X a data matrix.
#'
#' @return a data matrix in a long format.
#'
#' @importFrom tibble rownames_to_column as_data_frame
#' @importFrom tidyr gather
#'
#' @export
#' @examples
#' Z <- matrix(rnorm(100), 20)
#' Z.m <- melt_matrix(Z)
#' head(Z.m)
#'
melt_matrix <- function(X) {
  y <- X %>%
    as.data.frame(
      row.names = rownames(X),
      col.names = colnames(X),
      stringsAsFactors = FALSE
    ) %>%
    rownames_to_column("sample") %>%
    as_data_frame() %>%
    gather(key = "feature", value = "value", -sample) %>%
    as.data.frame(stringsAsFactors = FALSE)
  return(y)
}


#' @title Normalize Data
#'
#' @description Normalize data the assay data.
#'
#' @param object A \code{vistimeseq} object or a data matrix/data frame.
#' @param sample.norm.method Method for sample normalization.
#' Currently supports only scaling to a common factor,
#' "scale_common_factor" which with  \code{column.scale.factor} = 1e+06
#' is equivalent to CPM normalization.
#' @param column.scale.factor Sets the scale factor for sample-level
#' normalization
#'
#' @return Returns \code{vistimeseq} object after normalization.
#' Normalized data is stored \code{data} slot.
#'
#' @importFrom methods slot<- validObject
#' @export
#' @examples
#' endoderm_small
#' endoderm_small <- normalize_data(endoderm_small)
#' head(get_data(endoderm_small))
#'
normalize_data <- function(object,
                           sample.norm.method = "scale_common_factor",
                           column.scale.factor = 1e+06) {
  if (!class(object) %in% c("vistimeseq", "data.frame", "matrix")){
    stop("The argument 'object' must be either in either 'data.frame',",
         " 'matrix', or 'vistimeseq' class.")
  }
  if (class(object) == "vistimeseq") {
    if (!validObject(object))
      stop("Invalid vistimeseq object.")
    if (is.null(get_data(object, raw = TRUE))) {
      stop("Raw data for \"vistimeseq\" object has not been set")
    }
    raw.data <- normalized.data <- get_data(object, raw = TRUE)
  } else {
    raw.data <- normalized.data <- object
  }
  raw.data <- as.matrix(raw.data)
  if (sample.norm.method == "scale_common_factor") {
    normalized.data <- column.scale.factor *
      sweep(raw.data, 2, colSums(raw.data), "/")
  }
  slot(object, name = "data", check = TRUE) <- as.data.frame(normalized.data)
  return(object)
}


#' @title Collapse data over replicates.
#'
#' @description This function aggregates the data over replicates,
#' i.e. returns collapse data for each group and at each time point.
#'
#' @param object A \code{vistimeseq} object
#' @param FUN the aggreagate function. Default is mean

#' @return Returns \code{vistimeseq} object after collapsing
#' over replicates. Collapsed data is stored \code{sample.data.collapsed}
#' and \code{data.collapsed} slots.
#'
#' @importFrom dplyr mutate select
#' @importFrom  methods slot<-
#' @importFrom methods validObject
#' @importFrom stats aggregate
#' @export
#' @examples
#' endoderm_small
#' endoderm_small <- normalize_data(endoderm_small)
#' endoderm_small <- collapse_replicates(endoderm_small)
#' head(collapsed_data(endoderm_small))
#'
collapse_replicates <- function(object, FUN = mean) {
  group <- time <- NULL
  if (!validObject(object)){
    stop("Invalid vistimeseq object.")
  }
  dat <- get_data(object)
  sample.data.collapsed <-
    expand.grid(group = unique(get_group(object)),
                time = unique(get_time(object)),
                stringsAsFactors = FALSE) %>%
    mutate(sample = paste0(group, "_", time)) %>%
    select(sample, group, time)
  rownames(sample.data.collapsed) <- sample.data.collapsed$sample
  slot(object, name = "sample.data.collapsed", check = TRUE) <-
    sample.data.collapsed

  if (nrow(sample.data.collapsed) == n_samples(object)) {
    warning("Only single replicate per group found. ",
            "Collapsed data same as data.")
    colnames(dat) <- paste0(get_group(object), "_", get_time(object))
    dat <- dat[, sample.data.collapsed$sample]
    slot(object, name = "data.collapsed", check = TRUE) <- dat
    return(object)
  }
  data.collapsed <- aggregate(
    t(dat), list(get_group(object), get_time(object)), FUN)
  rownames(data.collapsed) <- paste0(data.collapsed$Group.1, "_",
                                     data.collapsed$Group.2)
  data.collapsed <- data.collapsed[, seq(3, ncol(data.collapsed))]
  data.collapsed <- t(data.collapsed)[, sample.data.collapsed$sample]
  slot(object, name = "data.collapsed", check = TRUE) <-
    as.data.frame(data.collapsed)
  return(object)
}


#' @title Data to time-course
#' @description Function that splits data to time series for each
#' replicate included.
#'
#' @param X a data matrix or data.frame where columns correspond to
#' samples and rows to features.
#' @param time a vector of length equal to ncol(X) indicating the time
#' variable corresponding to the sample (data column).
#' @param replicate a vector of length equal to ncol(X) indicating
#' the replicate variable corresponding to the sample (data column).
#' @param group a vector of length equal to ncol(X) indicating the group
#' membership corresponding to the sample (data column).
#' @return a data.frame with ordered time series for each replicate.
#'
#' @importFrom dplyr left_join select
#' @importFrom tidyr spread
#' @export
#' @examples
#' X <- matrix(rnorm(1000), ncol = 50)
#' group <- rep(c("A", "B"), each = 25)
#' replicate <- rep(paste0("rep", 1:5), each = 5)
#' time <- rep(1:5, 10)
#' tc <- data_to_tc(X, time, replicate, group)
#' head(tc)
#'
data_to_tc <- function(X, time, replicate = NULL, group = NULL){
  value <- feature <- NULL
  if (is.null(group)) group <- rep("G1", ncol(X))
  if (is.null(replicate)) replicate <- rep("R1", ncol(X))
  if (is.null(colnames(X))) colnames(X) <- seq_len(ncol(X))
  time.names <- as.character(sort(unique(time)))
  DF <- data.frame(
    sample = colnames(X),
    group, replicate, time,
    stringsAsFactors = FALSE)
  tc <- suppressMessages(
    melt_matrix(t(X)) %>%
      left_join(DF) %>%
      select(group, replicate, time, value, feature) %>%
      spread(key = time, value = value) %>%
      arrange(feature, group, replicate) %>%
      select("feature", "group", "replicate", time.names)
  )
  return(tc)
}


#' @title Convert data to time-course.
#'
#' @description This function converts the wide data matrix
#' to time-course long \code{data.frame} format where each
#' row gives data values over time (at each time point) for each
#' feature, group, and replicate.
#'
#' @param object A \code{vistimeseq} object
#' @param feature.trans.method Method for feature normalization. Default "none".
#' Currently supports only "none" (no transformation),
#' "scale_feat_sum" (scaling by feature sum),
#' or "var_stab" (variance stabilization). Default is "var_stab".
#' @param var.stabilize.method Method for variance stabilization (VST).
#' Currently, supports "none" (no VST), "log1p" (log plus one), "asinh"
#' (inverse hyperbolic sine) or "deseq"
#' (\code{\link[DESeq2]{varianceStabilizingTransformation}} function from
#' \code{DESeq2} package). Default is "log1p".
#'
#' @return Returns \code{vistimeseq} object after conversion to
#' time-course format. Converted data is stored in
#' \code{timecourse.data} slot.
#'
#' @importFrom  methods slot<-
#' @importFrom methods validObject
#' @export
#' @examples
#' endoderm_small
#' endoderm_small <- normalize_data(endoderm_small)
#' endoderm_small <- collapse_replicates(endoderm_small)     # Optional
#' endoderm_small <- convert_to_timecourse(endoderm_small)
#' head(time_course(endoderm_small))
#'
convert_to_timecourse <- function(
  object, feature.trans.method = "var_stab", var.stabilize.method = "log1p") {
  if (!validObject(object))
    stop("Invalid vistimeseq object.")
  dat <- get_data(object)
  dat_collapsed <- collapsed_data(object)

  if (feature.trans.method == "scale_feat_sum") {
    dat <- sweep(dat, 1, rowSums(dat), "/")
    if (!is.null(dat_collapsed)) {
      dat_collapsed <- sweep(dat_collapsed, 1, rowSums(dat_collapsed), "/")
    }
  } else if (feature.trans.method == "var_stab") {
    dat <- variance_stabilization(dat, var.stabilize.method)
    if (!is.null(dat_collapsed)) {
      dat_collapsed <- variance_stabilization(
        dat_collapsed, var.stabilize.method)
    }
  } else if (feature.trans.method != "none"){
    stop("Unsupported \"feature.trans.method\" chosen.")
  }
  timecourse.data <- list()
  timecourse.data[["tc"]] <- data_to_tc(
    dat,
    time = get_time(object),
    replicate = get_replicate(object),
    group = get_group(object)
  )
  if (!is.null(dat_collapsed)) {
    tc_cllps <- data_to_tc(
      dat_collapsed,
      time = collapsed_sample_data(object)$time,
      replicate = rep("Collapsed", ncol(dat_collapsed)),
      group = collapsed_sample_data(object)$group
    )
    timecourse.data[["tc_collapsed"]] <- tc_cllps
  }
  slot(object, name = "timecourse.data", check = TRUE) <- timecourse.data
  return(object)
}


#' @title Add differences to the time-course data.
#'
#' @description This function add lag difference features
#' to the time-course data.
#'
#' @param timecourse a data matrix or data.frame where each
#' column corresponds to consecutive time point.
#' @param lambda Weights for each lag difference, for time-course data.
#' Length of \code{lambda} specifies number of lags to include. By default
#' lag of order one and two are included with coefficients 0.5 and 0.25
#' respectively.
#'
#' @return a data matrix with added difference lags.
#'
add_lags_to_tc <- function(timecourse, lambda) {
  if(is.null(lambda)) {
    stop("Need to specify weights for lags.")
  }
  nT <- ncol(timecourse)
  timeNames <- colnames(timecourse)
  timecourse <- as.matrix(timecourse)
  lags <- lapply(seq_along(lambda), function(i) {
    ilag <- lambda[i] * t(diff(t(timecourse), lag = i))
    colnames(ilag) <- paste0("Lag_",
                             timeNames[seq((i+1), nT)], "_",
                             timeNames[seq(1,(nT-i))])
    return(ilag)
  })
  lags <- do.call("cbind", lags)
  res <- as.data.frame(cbind(timecourse, lags))
  return(res)
}


#' @title Add differences to the time-course data.
#'
#' @description This function concatenates lags to time-course
#' data stored in elements of the \code{timecourse.data} slot.
#'
#' @param object A \code{vistimeseq} object
#' @param lambda Weights for each lag difference, for time-course data.
#' Length of \code{lambda} specifies number of lags to include. By default
#' lag of order one and two are included with coefficients 0.5 and 0.25
#' respectively.
#'
#' @return Returns \code{vistimeseq} object with lags added to elements
#' in \code{timecourse.data} slot.
#'
#' @importFrom dplyr select contains
#' @importFrom  methods slot<-
#' @importFrom methods validObject
#' @export
#' @examples
#' endoderm_small
#' endoderm_small <- normalize_data(endoderm_small)
#' endoderm_small <- convert_to_timecourse(endoderm_small)
#' endoderm_small <- add_lags(endoderm_small)
#' head(time_course(endoderm_small, collapsed = FALSE))
#' head(time_course(endoderm_small, collapsed = TRUE))
#'
#' endoderm_small <- collapse_replicates(endoderm_small)
#' endoderm_small <- convert_to_timecourse(endoderm_small)
#' endoderm_small <- add_lags(endoderm_small)
#' head(time_course(endoderm_small, collapsed = TRUE))
#'
add_lags <- function(object, lambda = c(0.5, 0.25)) {
  group <- feature <- NULL
  if (!validObject(object)){
    stop("Invalid vistimeseq object.")
  }
  if(is.null(time_course(object))) {
    stop("No data in 'timecourse.data' slot. Use 'convert_to_timecourse()', ",
         "function to convert the data first.")
  }
  timecourse.data <- list()
  timecourse.data[["tc"]] <- time_course(object)
  timecourse.data[["tc_collapsed"]] <- time_course(object, collapsed = TRUE)

  for(tc_name in names(timecourse.data)) {
    tc <- timecourse.data[[tc_name]] %>%
      select(-contains("Lag_"))
    tc_with_lags <- add_lags_to_tc(
      tc %>% select(-feature, -group, -replicate), lambda = lambda)
    tc_with_lags <- cbind(tc %>% select(feature, group, replicate),
                          tc_with_lags)
    timecourse.data[[tc_name]] <- tc_with_lags
  }
  slot(object, name = "timecourse.data", check = TRUE) <- timecourse.data
  return(object)
}





