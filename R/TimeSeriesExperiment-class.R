### =========================================================================
### TimeSeriesExperiment objects
### -------------------------------------------------------------------------
#' @title Time Series Experiment Class
#'
#' @details The \code{TimeSeriesExperiment} class is an object 
#' (the main data container) used for in the time series/course experiment 
#' analysis. It stores all relevant information associated with the 
#' dataset, including the raw data, group, replicate and time associated 
#' with each sample (column of the data). The object includes also slots 
#' for results from some class-specific methods.
#'
#' @description \code{TimeSeriesExperiment} is an extension of 
#' \code{SummarizedExperiment} class with the following new slots 
#' (in addition to SummarizedExperiment slots).
#'
#' @slot timepoint A vector indicating the time-point of each sample collection.
#' @slot group A vector indicating the group membership of each sample.
#' @slot replicate A vector indicating the replicate id of each sample.
#' @slot assayCollapsed A matrix with assay data aggregated over replicates.
#' @slot colDataCollapsed A \link{DataFrame} where rows correspond
#' to samples aggregated over replicates and columns indicate group 
#' membership and time-point.
#' @slot timeSeries A list of time-course formatted data. Each element
#' of the list is a \link{DataFrame} with the first three columns 
#' indicating feature, group, replicate, and the remaining ones storing
#' the assay data at consecutive time points. \code{TimeSeriesExperiment} 
#' methods will typically generate elements of the list named: 'ts',
#' 'ts_collapsed'.
#' @slot dimensionReduction A list of results from applying 
#' dimmensionality reduction methods; elements named by technique used.
#' @slot clusterAssignment A list of results from clustering of the time-series 
#' features (rows) containing elements: 'settings', 'hclust', 'clust_map'
#' and 'clust_centroids'.
#' @slot differentialExpression A list of results from differential expression
#' analysis. Either from point-wise ('timepoint_de') or trajectory
#' ('trajectory_de') differential expression analysis.
#' @name TimeSeriesExperiment-class
#' @rdname TimeSeriesExperiment
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @aliases TimeSeriesExperiment-class
#' @exportClass TimeSeriesExperiment
#'
.TimeSeriesExperiment <- setClass(
    "TimeSeriesExperiment",
    slots = representation(
      timepoint = "numeric",
      group = "ANY",
      replicate = "ANY",
      assayCollapsed = "matrix",
      colDataCollapsed = "DataFrame",
      timeSeries = "list",
      dimensionReduction = "list",
      clusterAssignment = "list",
      differentialExpression = "list"
    ),
    contains="SummarizedExperiment"
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.TimeSeriesExperiment <- function(object)
{
  msg <- NULL
  nSamples <- ncol(object)
  for (sl in c("timepoint", "group", "replicate")) {
    n <- length(slot(object, sl))
    if (n != nSamples) {
      msg <- c(msg, paste0(
          "Length of ", sl, ", ", n, ", is not be equal to the number ",
          "of samples in the TimeSeriesExperiment object, ", nSamples))
    }
  }
  if(!is(slot(object, "colDataCollapsed"), "DataFrame")) {
    msg <- c(msg, paste0("'colDataCollapsed' must be a 'DataFrame'"))
  }
  if(!is(slot(object, "assayCollapsed"), "matrix")) {
    msg <- c(msg, paste0("'assayCollapsed' must be a 'matrix'"))
  }
  if(ncol(slot(object, "assayCollapsed")) != 
     nrow(slot(object, "colDataCollapsed"))) {
      msg <- c(
        msg, paste0("ncol() of 'assayCollapsed' not equal to nrow",
                    "of 'colDataCollapsed' slot."))
  }
  if (length(msg)) msg else TRUE
}

#' @importFrom S4Vectors setValidity2
#' @export
setValidity2("TimeSeriesExperiment", .valid.TimeSeriesExperiment)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @title Imputing missing samples
#' @description Estimates assay values for missing timepoints
#' @details Sometimes datasets include a missing sample. Here we impute
#' the expression values by taking the mean expression of two samples
#' from the same group and replicate (same individual) at two times 
#' surrounding the time of the missing sample
#' @param raw.data Raw input data with columns corresponding to samples
#' (observations) and rows to features
#' @param sample.data Optional. A \code{data.frame} object were rows are
#' samples (observations) and columns are sample attributes
#' (e.g. group/condition, replicate, timepoint)
#' @return list of modified raw.data and sample.data
#' @importFrom dplyr filter mutate
#' @importFrom S4Vectors DataFrame
#' @importFrom tibble add_row
.imputeData <- function(raw.data, sample.data) 
{
    group <- timepoint <- category <- NULL
    if(!any(c("group", "replicate", "timepoint") %in% colnames(sample.data))) {
        stop("Missing columns 'group', 'replicate', or ",
             "'timepoint' in 'sample.data'")
    }

    time_levels <- sort(unique(sample.data$timepoint))
    sample.names <- rownames(sample.data)
    sample.data <- sample.data %>%
        as.data.frame() %>%
        mutate(category = paste0(group, "_", replicate))

    for(id in unique(sample.data$category)) {
        smp_gr_timepoint <- filter(sample.data, category == id)
        gr <- strsplit(id, "_")[[1]][1]
        rep <- strsplit(id, "_")[[1]][2]
        missing_timepoint <- !(time_levels %in% smp_gr_timepoint$timepoint)
        if(any(missing_timepoint)) {
            message("Missing sample for group: ", gr, " replicate ", rep,
                    " at timepoint ", time_levels[missing_timepoint], ".")
            if (any(missing_timepoint[-1] & 
                    missing_timepoint[-length(missing_timepoint)])) {
                stop("Two or more consecutive samples are missing for the same",
                     "group-replicate. Please fix the experimental design.")
            }
            message("Imputing counts for the missing sample...")
            for (i in which(missing_timepoint)) {
                t_surround <- time_levels[
                    c(max(1, i -1), min(i+1, length(time_levels)))]
                surrounding_samples <- smp_gr_timepoint %>% 
                  filter(timepoint %in% t_surround)
                surrounding_samples <- surrounding_samples[["sample"]]
                if (length(surrounding_samples) == 1) {
                    missing_sample <- raw.data[, surrounding_samples]
                } else {
                    missing_sample <- rowMeans(raw.data[, surrounding_samples])
                }
                raw.data <- cbind(raw.data, missing_sample)
                missing_sample_name <- paste0(
                  "Imputed_G", gr, "_R", rep, "_T", time_levels[i])
                colnames(raw.data)[ncol(raw.data)] <- missing_sample_name
                sample.data <- add_row(
                    sample.data, sample = missing_sample_name,
                    group = gr, replicate = rep, timepoint = time_levels[i]
                )
                sample.names <- c(sample.names, missing_sample_name)
            }
        }
    }
    rownames(sample.data) <- sample.names
    return(list(raw.data = raw.data, sample.data = DataFrame(sample.data)))
}


#' @importFrom S4Vectors DataFrame SimpleList
#' @importMethodsFrom SummarizedExperiment colData colData<-
#' @importMethodsFrom SummarizedExperiment rowData rowData<-
#' @importMethodsFrom  SummarizedExperiment assays assays<- 
.processSE <- function(se, timepoint, group, replicate) {
      nSamples <- ncol(se)
      nFeatures <- nrow(se)
      
      if(is.null(colnames(se))) colnames(se) <- paste0("S", seq_len(nSamples))
      if(is.null(rownames(se))) rownames(se) <- paste0("F", seq_len(nFeatures))
      
      if(all(!is.null(timepoint), timepoint %in% colnames(colData(se)))){
        timepoint <- suppressMessages(as.numeric(colData(se)[[timepoint]]))
      }
      if (all(!is.null(replicate), replicate %in% colnames(colData(se)))) {
        replicate <- colData(se)[[replicate]]
      }
      if (all(!is.null(group), group %in% colnames(colData(se)))){
        group <- colData(se)[[group]]
      }
      if(any(is.null(timepoint), !is.numeric(timepoint), 
             all(is.na(timepoint)), length(timepoint) != nSamples)) {
        stop("Wrong timepoint input.")
      }
      if(is.null(replicate)) replicate <- rep("R1", nSamples)
      if(is.null(group)) group <- rep("G1", nSamples)
      
      timepoint <- as.numeric(timepoint)
      rowData(se)$feature   <- rownames(se)    
      colData(se)$sample    <- colnames(se)
      colData(se)$timepoint <- timepoint
      colData(se)$group     <- group
      colData(se)$replicate <- replicate
      
      if(length(assays(se)) > 1) {
        message("Keeping only the first supplied assay.")
      }
      imputed <- .imputeData(assays(se)[[1]], colData(se))
      se <- SummarizedExperiment(
        (list(raw = as.matrix(imputed$raw.data))),
        rowData = rowData(se),
        colData  = imputed$sample.data)
      return(se)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @title TimeSeriesExperiment constructor 
#' @description Constructor for 'TimeSeriesExperiment' object which stores
#' same data as \code{SummarizedExperiment} but also slots for time, group, 
#' replicate  and other time-series formated data useful for applyting 
#' data-analysis.
#' @details \code{TimeSeriesExperiment} constructor initializes the 
#' \code{TimeSeriesExperiment} TimeSeriesExperiment object and populates
#' the time, replicate, and group slots.
#' @name TimeSeriesExperiment-constructor
#' @rdname TimeSeriesExperiment
#' @param timepoint a vector indicating timepoint at which each sample 
#' was collected or a character string equal to one of the column names 
#' of a supplied colData.
#' @param group a vector indicating a group membership for each sample 
#' or a character string equal to one of the column names of a supplied 
#' colData. If not specified, the group is set to 'G1' for each sample.
#' @param replicate a vector indicating a replicate id of each sample 
#' or a character string equal to one of the column names of a supplied 
#' colData. If not specified, the replicate is set to 'R1' for each sample.
#' @param ... For \link[SummarizedExperiment:RangedSummarizedExperiment-class]{SummarizedExperiment:: SummarizedExperiment}, 
#' S4 methods \code{list} and \code{matrix}, arguments identical to those 
#' of the \code{SimpleList} method.
#' @return Returns an initialized TimeSeriesExperiment object.
#' @importFrom SummarizedExperiment SummarizedExperiment colData
#' @importFrom S4Vectors DataFrame
#' @export
#' @examples
#' raw <- matrix(runif(3000), ncol = 30)
#' timepoint <- rep(rep(1:5, each = 3), 2)
#' replicate <- rep(1:3, 10)
#' group <- rep(1:2, each = 15)
#' test_TimeSeriesExperiment <- TimeSeriesExperiment(
#'     assays = list(raw),
#'     timepoint = timepoint,
#'     replicate = replicate,
#'     group = group)
#' test_TimeSeriesExperiment
#' 
TimeSeriesExperiment <- function(
  ..., 
  timepoint=numeric(0), 
  group=character(0),
  replicate=character(0))
{
    se <- SummarizedExperiment(...)
    se <- .processSE(se, timepoint, group, replicate)
    object <- .TimeSeriesExperiment(
      se, 
      timepoint=colData(se)[["timepoint"]], 
      group=colData(se)[["group"]], 
      replicate=colData(se)[["replicate"]],  
      assayCollapsed = matrix(0,0,0),
      colDataCollapsed = DataFrame(),
      timeSeries = list(),
      dimensionReduction = list(),
      clusterAssignment = list(),
      differentialExpression = list())
    return(object)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' 
#' @title TimeSeriesExperiment constructor from SummarizedExperiment
#' @description \code{TimeSeriesExperiment} constructor initializes the 
#'  object from 
#' \link[SummarizedExperiment:RangedSummarizedExperiment-class]{SummarizedExperiment::SummarizedExperiment}
#' and populates the time, replicate, and group slots.
#' @details \code{TimeSeriesExperiment} is an extension of 
#' \code{SummarizedExperiment} class.
#' @param se  \link[SummarizedExperiment:RangedSummarizedExperiment-class]{SummarizedExperiment::SummarizedExperiment}
#' object
#' @param timepoint a vector indicating timepoint at which each sample 
#' was collected or a character string equal to one of the column names 
#' of a supplied \code{SummarizedExperiment}.
#' @param group a vector indicating a group membership for each sample 
#' or a character string equal to one of the column names of a supplied 
#' \code{SummarizedExperiment} If not specified, the group is set to 'G1' for
#' each sample.
#' @param replicate a vector indicating a replicate id of each sample 
#' or a character string equal to one of the column names of a supplied 
#' \code{SummarizedExperiment} If not specified, the replicate is set to 'R1' 
#' for each sample.
#'
#' @return Returns an initialized \code{TimeSeriesExperiment} object.
#' @importMethodsFrom  SummarizedExperiment colData
#' @importFrom S4Vectors DataFrame
#' @export
#'
#' @examples
#' raw <- matrix(runif(3000), ncol = 30)
#' pheno.data <- data.frame(
#'    time = rep(rep(1:5, each = 3), 2),
#'    replicate = rep(1:3, 10),
#'    group = rep(1:2, each = 15))
#' feature.data <- data.frame(
#'    feature = paste0("F", 1:100)
#' )
#' test_sumexp <- SummarizedExperiment::SummarizedExperiment(
#'    assays = list(counts = raw),
#'    rowData = feature.data, colData = pheno.data)
#' test_TimeSeriesExperiment <- 
#' makeTimeSeriesExperimentFromSummarizedExperiment(
#'    test_sumexp, timepoint = "time", group = "group",
#'    replicate = "replicate")
#' test_TimeSeriesExperiment
#' 
makeTimeSeriesExperimentFromSummarizedExperiment <- function(
  se, timepoint = NULL, group = NULL, replicate = NULL) 
{
    if(!is(se, "SummarizedExperiment")) 
        stop("'se' must be a SummarizedExperiment object")
    se <- .processSE(se, timepoint, group, replicate)
    object <- .TimeSeriesExperiment(
        se, 
        timepoint=colData(se)[[timepoint]], 
        group=colData(se)[[group]], 
        replicate=colData(se)[[replicate]],  
        assayCollapsed = matrix(0, 0, 0),
        colDataCollapsed = DataFrame(),
        timeSeries = list(),
        dimensionReduction = list(),
        clusterAssignment = list(),
        differentialExpression = list())
    return(object)
}
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#'
#' @title TimeSeriesExperiment constructor from ExpressionSet
#' @description \code{TimeSeriesExperiment} constructor initializes the 
#' \code{TimeSeriesExperiment} object from 
#' \code{ExpressionSet}
#' and populates the time, replicate, and group slots.
#'
#' @details \code{TimeSeriesExperiment} is an extension of 
#' \code{SummarizedExperiment} class.
#' 
#' @param eset ExpressionSet object
#' @param timepoint a vector indicating timepoint at which each sample 
#' was collected or a character string equal to one of the column names 
#' of a supplied \code{ExpressionSet}.
#' @param group a vector indicating a group membership for each sample 
#' or a character string equal to one of the column names of a supplied 
#' \code{ExpressionSet} If not specified, the group is set to 'G1' for
#' each sample.
#' @param replicate a vector indicating a replicate id of each sample 
#' or a character string equal to one of the column names of a supplied 
#' \code{ExpressionSet} If not specified, the replicate is set to 'R1' 
#' for each sample.
#'
#' @return Returns an initialized TimeSeriesExperiment object.
#' @importFrom SummarizedExperiment makeSummarizedExperimentFromExpressionSet 
#' @export
#'
#' @examples
#' raw <- matrix(runif(3000), ncol = 30)
#' pheno.data <- data.frame(
#'    time = rep(rep(1:5, each = 3), 2),
#'    replicate = rep(1:3, 10),
#'    group = rep(1:2, each = 15))
#' feature.data <- data.frame(
#'    feature = paste0("F", 1:100)
#' )
#' \dontrun{
#'     library(Biobase)
#'     test_eset <- ExpressionSet(
#'          raw, phenoData = AnnotatedDataFrame(pheno.data),
#'          featureData = AnnotatedDataFrame(feature.data))
#'     test_TimeSeriesExperiment <- makeTimeSeriesExperimentFromExpressionSet(
#'        test_eset, timepoint = "time", group = "group",
#'        replicate = "replicate")
#'     test_TimeSeriesExperiment
#' }
#'
#' 
makeTimeSeriesExperimentFromExpressionSet <- function(
    eset, timepoint = NULL, group = NULL, replicate = NULL) 
{
    se <- makeSummarizedExperimentFromExpressionSet(eset)
    object <- makeTimeSeriesExperimentFromSummarizedExperiment(
        se, timepoint, group, replicate)
    return(object)
}

