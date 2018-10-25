#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @return None
#' @export
#' @importFrom magrittr %>%
#' @examples
#' matrix(sample(30), 10) %>% head
#' @usage lhs \%>\% rhs
NULL


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
#' Y <- varianceStabilization(X, method = "asinh")
#' head(Y)
#'
varianceStabilization <- function(X, method = "asinh") {
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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
#' Z.m <- meltMatrix(Z)
#' head(Z.m)
#'
meltMatrix <- function(X) {
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @title Normalize Data
#'
#' @description Normalize data the assay data.
#'
#' @param object A \code{TimeSeriesExperiment} object or a data matrix/data frame.
#' @param sample.norm.method Method for sample normalization.
#' Currently supports only scaling to a common factor,
#' "scale_common_factor" which with  \code{column.scale.factor} = 1e+06
#' is equivalent to CPM normalization.
#' @param column.scale.factor Sets the scale factor for sample-level
#' normalization
#'
#' @return Returns \code{TimeSeriesExperiment} object after normalization.
#' Normalized data is stored \code{data} slot.
#'
#' @importFrom methods slot<- validObject
#' @importFrom methods is
#' @importFrom SummarizedExperiment assays Assays
#' @export
#' @examples
#' data("endoderm_small")
#' names(assays(endoderm_small))
#' endoderm_small <- normalizeData(endoderm_small)
#' assays(endoderm_small)$norm[1:10, 1:6]
#'
normalizeData <- function(object, sample.norm.method = "scale_common_factor",
                          column.scale.factor = 1e+06) 
{
    if (!any(is(object, "TimeSeriesExperiment"),
             is(object, "data.frame"),
             is(object, "matrix"))){
        stop("The argument 'object' must be either in either 'data.frame',",
             " 'matrix', or 'TimeSeriesExperiment' class.")
    }
    message("Normalizing data...")
    if (is(object, "TimeSeriesExperiment")) {
        if (!validObject(object))
            stop("Invalid TimeSeriesExperiment object.")
        if (is.null(assays(object)$raw)) {
            stop("Raw data for 'TimeSeriesExperiment' object has not been set")
        }
        curr.assays <- assays(object)
        raw.data <- normalized.data <- curr.assays$raw
    } else {
        raw.data <- normalized.data <- object
    }
    normalized.data <- as.matrix(raw.data)
    if (sample.norm.method == "scale_common_factor") {
        normalized.data <- column.scale.factor *
            sweep(raw.data, 2, colSums(raw.data), "/")
    }
    if (is(object, "TimeSeriesExperiment") ) {
        curr.assays$norm <- normalized.data
        slot(object, name = "assays", check = TRUE) <- Assays(curr.assays)
    }
    return(object)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @title Collapse data over replicates.
#'
#' @description This function aggregates the data over replicates,
#' i.e. returns collapse data for each group and at each time point.
#'
#' @param object A \code{TimeSeriesExperiment} object
#' @param FUN the aggreagate function. Default is mean

#' @return Returns \code{TimeSeriesExperiment} object after collapsing
#' over replicates. Collapsed data is stored \code{sample.data.collapsed}
#' and \code{data.collapsed} slots.
#'
#' @importFrom dplyr mutate select
#' @importFrom  methods slot<-
#' @importFrom methods validObject
#' @importFrom stats aggregate
#' @importFrom SummarizedExperiment assays
#' @export
#' @examples
#' data("endoderm_small")
#' endoderm_small <- collapseReplicates(endoderm_small)
#' assayCollapsed(endoderm_small)[1:10, 1:6]
#'
collapseReplicates <- function(object, FUN = mean) {
    if (!is(object, "TimeSeriesExperiment")) 
      stop("Input must be a 'TimeSeriesExperiment' object.")
    if (!validObject(object)) 
      stop("Invalid TimeSeriesExperiment object.")
    if(!"norm" %in% names(assays(object))) {
        object <- normalizeData(object)
    }
    message("Aggregating across replicates...")
    group <- timepoint <- NULL
    sample.data.collapsed <-
        expand.grid(
          group = unique(groups(object)),
          timepoint = unique(timepoints(object)),
          stringsAsFactors = FALSE) %>%
        mutate(sample = paste0(group, "_", timepoint)) %>%
        select(sample, group, timepoint)
    slot(object, name = "colDataCollapsed", check = TRUE) <-
        DataFrame(sample.data.collapsed)
    curr.assays <- assays(object)
    norm.assay <- curr.assays$norm
    if (nrow(sample.data.collapsed) == ncol(object)) {
        warning("Only single replicate per group found => ",
                "collapsed data same as data. Renaming samples.")
        colnames(norm.assay) <- paste0(groups(object), "_", timepoints(object))
        collapsed.assay <- norm.assay[, sample.data.collapsed$sample]
    } else {
        collapsed.assay <- aggregate(
          t(norm.assay), 
          list(groups(object), timepoints(object)),
          FUN, na.rm = TRUE)
        rownames(collapsed.assay) <- paste0(
          collapsed.assay$Group.1, "_", collapsed.assay$Group.2)
        collapsed.assay <- collapsed.assay[, seq(3, ncol(collapsed.assay))]
        collapsed.assay <- t(collapsed.assay)[, sample.data.collapsed$sample]
    }
    slot(object, name = "assayCollapsed", check = TRUE) <- collapsed.assay
    return(object)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @title Data to time-series
#' @description Function that splits data to time series for each
#' replicate included.
#'
#' @param X a data matrix or data.frame where columns correspond to
#' samples and rows to features.
#' @param timepoint a vector of length equal to ncol(X) indicating the time
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
#' tc <- dataToTimeSeries(X, time, replicate, group)
#' head(tc)
#'
dataToTimeSeries <- function(X, timepoint, group = NULL, replicate = NULL){
    value <- feature <- NULL
    if (is.null(group)) group <- rep("G1", ncol(X))
    if (is.null(replicate)) replicate <- rep("R1", ncol(X))
    if (is.null(colnames(X))) colnames(X) <- seq_len(ncol(X))
    DF <- data.frame(
        sample = colnames(X),
        group, replicate, timepoint,
        stringsAsFactors = FALSE)
    time.names <- as.character(sort(unique(timepoint)))
    
    ts <- suppressMessages(
        meltMatrix(t(X)) %>%
            left_join(DF) %>%
            select(group, replicate, timepoint, value, feature) %>%
            spread(key = timepoint, value = value) %>%
            arrange(feature, group, replicate) %>%
            select("feature", "group", "replicate", time.names)
    )
    return(ts)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @title Convert data to time-series
#'
#' @description This function converts the wide data matrix
#' to time-course long \code{data.frame} format where each
#' row gives data values over time (at each time point) for each
#' feature, group, and replicate.
#'
#' @param object A \code{TimeSeriesExperiment} object
#' @param feature.trans.method Method for feature normalization. 
#' Default "none". Currently supports only "none" (no transformation),
#' "scale_feat_sum" (scaling by feature sum),
#' or "var_stab" (variance stabilization). Default is "var_stab".
#' @param var.stabilize.method Method for variance stabilization (VST).
#' Currently, supports "none" (no VST), "log1p" (log plus one), "asinh"
#' (inverse hyperbolic sine) or "deseq"
#' (\code{\link[DESeq2]{varianceStabilizingTransformation}} function from
#' \code{DESeq2} package). Default is "log1p".
#'
#' @return Returns \code{TimeSeriesExperiment} object after conversion to
#' time-course format. Converted data is stored in
#' \code{timecourse.data} slot.
#'
#' @importFrom  methods slot<-
#' @importFrom methods validObject
#' @importFrom SummarizedExperiment assays 
#' @importFrom S4Vectors DataFrame
#' @export
#' @examples
#' data("endoderm_small") 
#' endoderm_small <- makeTimeSeries(endoderm_small)
#' names(timeSeries(endoderm_small))
#' head(timeSeries(endoderm_small)[[1]])
#'
makeTimeSeries <- function(
    object, feature.trans.method = "var_stab", var.stabilize.method = "asinh")
{   
    if (!is(object, "TimeSeriesExperiment")) 
        stop("Input must be a 'TimeSeriesExperiment' object.")
    if (!validObject(object)) 
        stop("Invalid TimeSeriesExperiment object.")
    if(!"norm" %in% names(assays(object))) {
        message("Normalizing data...")
        object <- normalizeData(object)
    }
    message("Converting to timeseries format...")
    norm.assay <- assays(object)$norm
    collapsed.assay <- assayCollapsed(object)

    if (feature.trans.method == "scale_feat_sum") {
        trans.assay <- sweep(norm.assay, 1, rowSums(norm.assay), "/")
        if (!is.null(collapsed.assay)) {
            collapsed.trans.assay <- sweep(
                collapsed.assay, 1, rowSums(collapsed.assay), "/")
        }
    } else if (feature.trans.method == "var_stab") {
        trans.assay <- varianceStabilization(
            norm.assay, var.stabilize.method)
        if (!is.null(collapsed.assay)) {
            collapsed.trans.assay <- varianceStabilization(
                collapsed.assay, var.stabilize.method)
        }
    } else if (feature.trans.method != "none"){
        stop("Unsupported 'feature.trans.method' chosen.")
    }
    timeseries.data <- list()
    timeseries.data[["ts"]] <- dataToTimeSeries(
      norm.assay,
      timepoint = timepoints(object),
      group = groups(object),
      replicate = replicates(object)
    )
    timeseries.data[["ts_trans"]] <- dataToTimeSeries(
        trans.assay,
        timepoint = timepoints(object),
        group = groups(object),
        replicate = replicates(object)
    )
    if (dim(collapsed.assay)[[1]]) {
        timeseries.data[["ts_collapsed"]] <- dataToTimeSeries(
            collapsed.trans.assay,
            timepoint = colDataCollapsed(object)$timepoint,
            replicate = rep("Collapsed", ncol(collapsed.assay)),
            group = colDataCollapsed(object)$group
        )
    }
    slot(object, name = "timeSeries", check = TRUE) <- timeseries.data
    return(object)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @title Add differences to the time-course data.
#'
#' @description This function add lag difference features
#' to the time-course data.
#'
#' @param timeseries a data matrix or data.frame where each
#' column corresponds to consecutive time point.
#' @param lambda Weights for each lag difference, for time-course data.
#' Length of \code{lambda} specifies number of lags to include. By default
#' lag of order one and two are included with coefficients 0.5 and 0.25
#' respectively.
#'
#' @return a data matrix with added difference lags.
#'
.addLagsToTimeSeries <- function(timeseries, lambda) {
    if(is.null(lambda)) stop("Need to specify weights for lags.")
    nT <- ncol(timeseries)
    time_names <- colnames(timeseries)
    timeseries <- as.matrix(timeseries)
    lags <- lapply(seq_along(lambda), function(i) {
        ilag <- lambda[i] * t(diff(t(timeseries), lag = i))
        colnames(ilag) <- paste0(
          "Lag_", time_names[seq((i+1), nT)], "_", time_names[seq(1,(nT-i))])
        return(ilag)
    })
    lags <- do.call("cbind", lags)
    res <- as.data.frame(cbind(timeseries, lags))
    return(res)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @title Add differences to the time-course data.
#'
#' @description This function concatenates lags to time-course
#' data stored in elements of the \code{timeSeries} slot.
#'
#' @param object A \code{TimeSeriesExperiment} object
#' @param lambda Weights for each lag difference, for time-course data.
#' Length of \code{lambda} specifies number of lags to include. By default
#' lag of order one and two are included with coefficients 0.5 and 0.25
#' respectively.
#'
#' @return Returns \code{TimeSeriesExperiment} object with lags added 
#' to elements in \code{timeSeries} slot.
#'
#' @importFrom dplyr select contains
#' @importFrom  methods slot<-
#' @importFrom methods validObject
#' @importFrom SummarizedExperiment assays 
#' @export
#' @examples
#' data("endoderm_small")
#' endoderm_small <- collapseReplicates(endoderm_small)
#' endoderm_small <- makeTimeSeries(endoderm_small)
#' endoderm_small <- addLags(endoderm_small)
#' head(timeSeries(endoderm_small, "ts"))
#' head(timeSeries(endoderm_small, "ts_collapsed"))
#'
addLags <- function(object, lambda = c(0.5, 0.25)) {
    if (!is(object, "TimeSeriesExperiment")) 
        stop("Input must be a 'TimeSeriesExperiment' object.")
    if (!validObject(object)) 
        stop("Invalid TimeSeriesExperiment object.")
    group <- feature <- NULL
    if(length(timeSeries(object)) == 0) {
        stop("No data in 'timeSeries' slot. Use 'makeTimeSeries()' function ",
             "to convert the data first.")
    }
    message("Adding lags with coefficients: ", 
            paste0(lambda, collapse = " "), "...")
    timeseries <- timeSeries(object, NULL)
    for(ts_name in names(timeseries)) {
        ts <- timeseries[[ts_name]] %>%
            select(-contains("Lag_"))
        ts_with_lags <- .addLagsToTimeSeries(
            ts %>% select(-feature, -group, -replicate), lambda = lambda)
        ts_with_lags <- cbind(
          ts %>% select(feature, group, replicate), ts_with_lags)
        timeseries[[ts_name]] <- ts_with_lags
    }
    slot(object, name = "timeSeries", check = TRUE) <- timeseries
    return(object)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @title Filter features
#'
#' @description Filters \code{TimeSeriesExperiment} object to keep only chosen 
#' features. All relevant slots are updates. 
#' @details The slots for collapsed data and time series formatted data are
#' filtered accordingly, but \code{dimensionReduction}, 
#' \code{clusterAssignment} and \code{differentialExpression} are reset to 
#' \code{NULL} as different set of features would output in different results.
#'
#' @param object TimeSeriesExperiment object
#' @param features features (genes) to keep
#'
#' @return "TimeSeriesExperiment" object
#' @importFrom  dplyr filter
#' @importFrom methods validObject new
#' 
#' @export
#' @examples
#' data("endoderm_small")
#' features <- 1:100
#' endoderm_small <- filterFeatures(endoderm_small, features)
#' endoderm_small
#'
filterFeatures <- function (object, features) {
    feature <- NULL
    if (!is(object, "TimeSeriesExperiment")) 
      stop("Input must be a 'TimeSeriesExperiment' object.")
    if (!validObject(object))
      stop("Invalid 'TimeSeriesExperiment' object.")
    
    if (all(is.numeric(features), 
            !all(features %in% seq_along(rownames(object))))){
      stop("Some 'features' not found in rownames(object).")
    }
    if(is.numeric(features)){
      features <- rownames(object)[features]
    }
    if (!all(features %in% rownames(object))){
      stop("Some 'features' not found in rownames(object).")
    }
    
    object <- object[features, ]
    fltr_collapsed_data <- NULL
    if(nrow(assayCollapsed(object)) != 0) {
      fltr_collapsed_data <- assayCollapsed(object)[features, ]
      slot(object, "assayCollapsed", check = TRUE) <- fltr_collapsed_data
    }
    
    fltr_timecourse <- list()
    for(ts_name in names(timeSeries(object))) {
        fltr_timecourse[[ts_name]] <- timeSeries(object, ts_name) %>%
            filter(feature %in% features)
    }
    slot(object, "timeSeries", check = TRUE) <- fltr_timecourse
    
    if(length(dimensionReduction(object)) > 0 ) {
      message("Dimensionality reduction results are reset due to feature ",
              "filtering and need to be recomputed.")
    }
    if(length(clusterAssignment(object)) > 0) {
      message("Feature clustering results are reset due to feature ",
              "filtering and need to be recomputed.")
    }
    if(length(differentialExpression(object)) > 0) {
      message("Differential expression results are reset due to feature ",
              "filtering and need to be recomputed.")
    }
    slot(object, name = "dimensionReduction", check = TRUE) <- list()
    slot(object, name = "clusterAssignment", check = TRUE) <- list()
    slot(object, name = "differentialExpression", check = TRUE) <- list()
  
    validObject(object)
    return(object)
}



