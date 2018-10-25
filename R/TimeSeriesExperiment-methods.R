### =========================================================================
### Getter and setter methods for TimeSeriesExperiment
### -------------------------------------------------------------------------

### Helpers - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

.resetResults <- function(object){
    slot(object, name = "assays", check = TRUE) <- 
      Assays(SimpleList(list(raw = object@assays[["raw"]])))
    slot(object, name = "colDataCollapsed", check = TRUE) <- DataFrame()
    slot(object, name = "timeSeries", check = TRUE) <- list()
    slot(object, name = "dimensionReduction", check = TRUE) <- list()
    slot(object, name = "clusterAssignment", check = TRUE) <- list()
    slot(object, name = "differentialExpression", check = TRUE) <- list()
    return(object)
}


.replaceNames <- function(query, orig, replace) {
    if(length(orig) != length(replace)) {
        stop("orig and replace must be of the same length")
    }
    if (!all(query %in% orig)) {
        stop("values in query not in orig")
    }
    map <- replace
    names(map) <- orig
    return(map[query])
}


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

#' @importFrom methods slot slot<-
#' @export
setReplaceMethod(
  "colnames", "TimeSeriesExperiment", function(x, value) {
    if(length(value) != ncol(x))
      stop("Wrong length of substitute vector for colnames.")
    newobject <- x
    colnames(newobject) <- value
    # Additionally modify the dimensionRediuction slot
    dimreds <- slot(newobject, "dimensionReduction")
    if(is.null(names(dimreds))) return(newobject)
    for (name_mat in names(dimreds)[grep("_sample", names(dimreds))]) {
        base::rownames(dimreds[[name_mat]]) <- value
    }
    slot(newobject, "dimensionReduction") <- dimreds
    validObject(newobject)
    return(newobject)
  }
)

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

#' @importFrom methods slot slot<-
#' @export
setReplaceMethod(
    "rownames", "TimeSeriesExperiment", function(x, value) {
    if(length(value) != nrow(x))
        stop("Wrong length of substitute vector for rownames.")
    old_rownames <- rownames(x)
    newobject <- x
    rownames(newobject) <- value
    
    ts <- slot(newobject, "timeSeries")
    for(i in seq_along(slot(newobject, "timeSeries"))){
        ts[[i]]$feature <- .replaceNames(
            ts[[i]]$feature, old_rownames, value)
    }
    slot(newobject, "timeSeries") <- ts
    
    cl <- slot(newobject , "clusterAssignment")
    if("cluster_map" %in% names(cl)) {
        cl[["cluster_map"]]$feature <- .replaceNames(
            cl[["cluster_map"]]$feature, old_rownames, value)
    }
    if("hclust" %in% names(cl)) {
        cl[["hclust"]]$labels <- .replaceNames(
            cl[["hclust"]]$labels, old_rownames, value)
    }
    slot(newobject, "clusterAssignment") <- cl   

    dimreds <- slot(newobject, "dimensionReduction")
    if(is.null(names(dimreds))) return(newobject)
    for (name_mat in names(dimreds)[grep("_feature", names(dimreds))]) {
        rownames(dimreds[[name_mat]]) <- value
    }
    slot(newobject, "dimensionReduction") <- dimreds
    

    de_res <- slot(newobject, "differentialExpression")
    if ("timepoint_de" %in% names(de_res)){
        for(tmp in names(de_res$timepoint_de)){
            de_res$timepoint_de[[tmp]]$feature <-
                .replaceNames(de_res$timepoint_de[[tmp]]$feature,
                              orig = old_rownames, replace = value)
        }
    }
    if ("trajectory_de" %in% names(de_res)){
        tmp_de$trajectory_de$feature <-
            .replaceNames(tmp_de$trajectory_de$feature,
                          orig = old_rownames,
                          replace = value)
    }
    slot(newobject, "differentialExpression") <- de_res
    
    validObject(newobject)
    return(newobject)
  }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

#' @importMethodsFrom SummarizedExperiment "rowData<-"
#' @export
setReplaceMethod("rowData", "TimeSeriesExperiment", function(x, ..., value) {
    if(nrow(value) != nrow(x))
        stop("nrow(value) does not match the input object dimensions.")
    if(!"feature" %in% base::colnames(value)) {
        stop("rowData data must contain columns 'feature'.")
    }
    if(!all(value$feature == rownames(x))) {
      rownames(x) <- value$feature
    }
    newobject <- callNextMethod()
    validObject(newobject)
    return(newobject)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @importMethodsFrom SummarizedExperiment "colData<-"
#' @export
setReplaceMethod(
  "colData", "TimeSeriesExperiment",
  function(x, ..., value) {
        cols <- c("sample", "timepoint", "group", "replicate")
        if(!all(cols %in% base::colnames(value))) {
          stop("colData() must contain columns: 'sample', 'timepoint', 
               'group', ", "and 'replicate'.")
        }
        old_colData <- colData(x)
        print(value)
        newobject <- callNextMethod()
      
        mgs <- NULL
        if(!all(value$group == old_colData$group)){
          msg <- c(msg, "new group assignment")
          groups(newobject) <- value$group
        }
        if(!all(value$replicate == old_colData$replicate)){
          msg <- c(msg, "new replicate assignment")
          replicates(newobject) <- value$replicate
        }
        if(!all(value$timepoint == old_colData$timepoint)){
          msg <- c(msg, "new timepoint assignment")
          timepoints(newobject) <- as.numeric(value$timepoint)
        }
        if(length(msg)) {
          message(msg)
          newobject <- .resetResults(newobject)
        }
        validObject(newobject)
        return(newobject)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @title Timepoint information
#' @docType methods
#' @rdname timepoints
#' @param object a \code{TimeSeriesExperiment} object.
#' @param value a numeric vector with new time information.
#' @param ... argiments to other functions.
#' 
#' @return a numeric vector
#' @export
#' @examples
#' data("endoderm_small")
#' head(timepoints(endoderm_small))
#' timepoints(endoderm_small) <- sample(1:ncol(endoderm_small))
#' head(timepoints(endoderm_small))
#'
setGeneric("timepoints", function(object, ...) standardGeneric("timepoints"))

#' @rdname timepoints
#' @export
setMethod("timepoints", "TimeSeriesExperiment", function(object) {
    out <- slot(object, name = "timepoint")
    names(out) <- colnames(object)
    out
})

.setTime <- function(object, value) {
    tmp <- slot(object, "timepoint")
    if(length(value) != length(tmp)) {
        stop("[setting timepoints] vector length doesn't match object", 
             " dimensions")
    }
    if(!all(is.numeric(value))) stop("'timepoint' must be numeric")
  
    slot(object, name = "timepoint", check = TRUE) <- value
    colData(object)$timepoint <- value
    message("new timepoint assignment; all prior results reset to NULL.")
    object <- .resetResults(object)
    return(object)
}

#' @rdname timepoints
#' @export
setGeneric("timepoints<-", function(object, ..., value) {
    standardGeneric("timepoints<-")})

#' @rdname timepoints
#' @exportMethod "timepoints<-"
setReplaceMethod(f = "timepoints", signature = "TimeSeriesExperiment",
                 definition = .setTime)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @title Group information 
#' @docType methods
#' @rdname groups
#'
#' @param object a \code{TimeSeriesExperiment} object.
#' @param value a character vector with new group membership
#' @param ... argiments to other functions.
#'
#' @return a character vector
#' @export
#' @examples
#' data("endoderm_small")
#' head(groups(endoderm_small))
#' groups(endoderm_small) <- sample(c("A", "B"), 
#'     ncol(endoderm_small), replace = TRUE)
#' head(groups(endoderm_small))
#'
setGeneric("groups", function(object, ...) standardGeneric("groups"))

#' @rdname groups
#' @export
setMethod("groups", "TimeSeriesExperiment", function(object) {
    out <- object@group
    names(out) <- colnames(object)
    out
  })

.setGroup <- function(object, value) {
    if(length(value) != length(object@group)) {
        stop("[setting groups] vector length doesn't match object dimensions")
    }
    slot(object, name = "group", check = TRUE) <- value
    message("new group assignment; all prior results reset to NULL.")
    object <- .resetResults(object)
    return(object)
}

#' @rdname groups
#' @export
setGeneric("groups<-", function(object, ..., value) {
  standardGeneric("groups<-")})

#' @rdname groups
#' @exportMethod "groups<-"
setReplaceMethod(f = "groups", signature = "TimeSeriesExperiment",
                 definition = .setGroup
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @title Replicate information
#' @docType methods
#' @rdname replicates
#'
#' @param object a \code{TimeSeriesExperiment} object.
#' @param value a character vector with new replicate ids.
#' @param ... argiments to other functions.
#'
#' @return a character vector
#' @export
#' @examples
#' data("endoderm_small")
#' head(replicates(endoderm_small))
#' replicates(endoderm_small) <- sample(c("R1", "R2", "R3"), 
#'     ncol(endoderm_small), replace = TRUE)
#' head(replicates(endoderm_small))
setGeneric("replicates", function(object, ...) standardGeneric("replicates"))

#' @rdname replicates
#' @export
setMethod("replicates", "TimeSeriesExperiment", function(object) {
  out <- object@replicate
  names(out) <- colnames(object)
  out
})

.setReplicate <- function(object, value) {
    if(length(value) != length(object@replicate)) {
        stop("[setting replicates] vector length doesn't match object ", 
             "dimensions")
    }
  slot(object, name = "replicate", check = TRUE) <- value
  message("new group assignment; all prior results reset to NULL.")
  object <- .resetResults(object)
  return(object)
}

#' @rdname replicates
#' @export
setGeneric("replicates<-", function(object, ..., value)
  standardGeneric("replicates<-"))

#' @rdname replicates
#' @exportMethod "replicates<-"
setReplaceMethod(f = "replicates", signature = "TimeSeriesExperiment",
                 definition = .setReplicate
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @title Collapsed data.
#' @description Assay and colData collapsed over replicates. The values 
#' can be computed with set with \link{collapseReplicates} function.
#' @docType methods
#' @rdname collapsed-data
#' 
#' @param object a \code{TimeSeriesExperiment} object.
#' @param value a numerical matrix
#' @param ... argiments to other functions.
#' 
#' @return a \link{DataFrame}
#' @export
#' @examples
#' data("endoderm_small")
#' endoderm_small <- collapseReplicates(endoderm_small)
#' head(assayCollapsed(endoderm_small))
#' 
setGeneric("assayCollapsed", function(object, ...) {
  standardGeneric("assayCollapsed")
})

#' @rdname collapsed-data
#' @export
setMethod(f = "assayCollapsed", signature = "TimeSeriesExperiment",
          function(object) { return(object@assayCollapsed) }
)

#' @rdname collapsed-data
#' @export
setGeneric("assayCollapsed<-", function(object, ..., value) {
  standardGeneric("assayCollapsed<-")})

.setAssayCollapsed <- function(object, value) {
  if(ncol(object@assayCollapsed) != nrow(object@colDataCollapsed)) {
    stop("dimensions of collapsed assay data and collapsed colData does ",
         "not match.")
  }
  slot(object, name = "assayCollapsed", check = TRUE) <- value
  return(object)
}

#' @rdname collapsed-data
#' @exportMethod "assayCollapsed<-"
setReplaceMethod(f = "assayCollapsed", signature = "TimeSeriesExperiment",
                 definition = .setAssayCollapsed)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @docType methods
#' @rdname collapsed-data
#' @export
#' @examples
#' data("endoderm_small")
#' endoderm_small <- collapseReplicates(endoderm_small)
#' head(colDataCollapsed(endoderm_small))
#' 
setGeneric("colDataCollapsed", function(object, ...) {
  standardGeneric("colDataCollapsed")
})

#' @rdname collapsed-data
#' @export
setMethod(f = "colDataCollapsed", signature = "TimeSeriesExperiment",
          function(object) { return(object@colDataCollapsed) }
)

#' @rdname collapsed-data
#' @export
setGeneric("colDataCollapsed<-", function(object, ..., value) {
  standardGeneric("colDataCollapsed<-")})

.setColDataCollapsed <- function(object, value) {
  if(ncol(object@assayCollapsed) != nrow(object@colDataCollapsed)) {
    stop("dimensions of collapsed assay data and collapsed colData does ",
         "not match.")
  }
  slot(object, name = "colDataCollapsed", check = TRUE) <- value
  return(object)
}

#' @rdname collapsed-data
#' @exportMethod "colDataCollapsed<-"
setReplaceMethod(f = "colDataCollapsed", signature = "TimeSeriesExperiment",
                 definition = .setColDataCollapsed)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @title Time series formatted data.
#'
#' @description Getter and setter methods for \code{timeSeries} slot of
#' a \code{TimeSeriesExperiment} object.
#' @details \code{timeSeries} slot is a list with 'ts' and (optionally) 
#' 'ts_collapsed' storing data formatted as time-series/time-courses.
#' @docType methods
#' @rdname timeSeries
#' @param object a \code{TimeSeriesExperiment} object.
#' @param name a character string, one of 'ts', 'ts_with_lags', 'ts_collapsed'
#' and 'ts_collapsed_with_lags'. If NULL, all elements are returned.
#' @param value replacement list
#' @param ... argiments to other functions.
#' 
#' @return a \code{data.frame}
#' @export
#' @examples
#' data("endoderm_small")
#' endoderm_small <- makeTimeSeries(endoderm_small)
#' head(timeSeries(endoderm_small))
#' head(timeSeries(endoderm_small, name = 'ts'))
#' 
setGeneric("timeSeries", function(object, ...) standardGeneric("timeSeries"))

#' @rdname timeSeries
#' @export
setMethod("timeSeries", "TimeSeriesExperiment", function(object, name = NULL) {
    if(is.null(name)) return(object@timeSeries)
    tsNames <- c('ts', 'ts_trans', 'ts_collapsed')
    if(!name %in% tsNames) {
        stop("'", name, "' not in 'timeSeries' slot")
    }
    else return(object@timeSeries[[name]])
})

#' @rdname timeSeries
#' @export
setGeneric("timeSeries<-", function(object, ..., value) {
    standardGeneric("timeSeries<-")})

.setTimeSeries <- function(object, value) {
    if (!all(is(value, "list"), names(value) %in% c("ts", "ts_collapsed"))) {
      stop("dimensions of collapsed assay data and collapsed colData does ",
           "not match.")
    }
    slot(object, "timeSeries", check = TRUE) <- value
    return(object)
}

#' @rdname timeSeries
#' @exportMethod "timeSeries<-"
setReplaceMethod(f = "timeSeries", signature = "TimeSeriesExperiment",
                 definition = .setTimeSeries)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @title Dimension reduction results
#' @description Getter methods for \code{dimensionReduction} slot
#' of a \code{TimeSeriesExperiment} object. The slot is a list of 
#' \code{data.frames}: 'pca_sample', 'pca_feature' and 'pca_eigs' storing 
#' results from a PCA projection.
#' @docType methods
#' @rdname dimensionReduction
#' @param object a \code{TimeSeriesExperiment} object.
#' @param name one of elements of 'dimensionReduction' slot: 'pca_sample', 
#' 'pca_feature' and 'pca_eigs' for returning the entire list. If NULL, 
#' all elements are returned.
#' @param ... argiments to other functions.
#' 
#' @return a \code{data.frame} or a list of \code{data.frame}s
#' @export
#' @examples
#' data("endoderm_small")
#' endoderm_small <- runPCA(endoderm_small)
#' head(dimensionReduction(endoderm_small, "pca_sample")[, 1:3])
#'
setGeneric("dimensionReduction", function(object, ...) {
  standardGeneric("dimensionReduction")})


#' @rdname dimensionReduction
#' @export
setMethod(
    "dimensionReduction", "TimeSeriesExperiment", 
    function(object, name = NULL) {
        if(is.null(name)) return(object@dimensionReduction)
        if(!name %in% c('pca_sample', 'pca_feature', 'pca_eigs')){
            stop("'name': ", name, " not in 'dimensionReduction' slot")
        }
        else return(object@dimensionReduction[[name]])
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @title Cluster analysis results
#' 
#' @description Getter methods for \code{clusterAssignment} slot
#' of a \code{TimeSeriesExperiment} object. The slot is a list with 
#' with elements named: 'settings', 'hclust', 'cluster_map',
#' 'clust_centroids' storing results from running \link{clusterTimeSeries}
#' function.
#'
#' @docType methods
#' @rdname clusterAssignment
#'
#' @param object a \code{TimeSeriesExperiment} object.
#' @param name one of elements of 'clusterAssignment' slot: 'settings', 
#' 'hclust', 'cluster_map', 'clust_centroids'. If NULL, all elements are 
#' returned.
#' @param ... argiments to other functions.
#'
#' @return a \code{data.frame}
#' @export
#' @examples
#'
#' data("endoderm_small")
#' endoderm_small <- clusterTimeSeries(endoderm_small)
#' clusterAssignment(endoderm_small, name = 'settings')
#' head(clusterAssignment(endoderm_small, name = 'final_cluster_map'))
#' head(clusterAssignment(endoderm_small, name = 'clust_centroids'))
#'
setGeneric("clusterAssignment", function(object, ...){
  standardGeneric("clusterAssignment")})

#' @rdname clusterAssignment
#' @export
setMethod(
  "clusterAssignment", "TimeSeriesExperiment", function(object, name = NULL) {
      if(is.null(name)) return(object@clusterAssignment)
      clust_res <- names(object@clusterAssignment)
      if(!name %in% clust_res){
          stop(name, " not in 'clusterAssignment' slot")
      }
      else return(object@clusterAssignment[[name]])
})

#' @rdname clusterAssignment
#' @export
setGeneric("clusterMap", function(object, ...){
  standardGeneric("clusterMap")})

#' @rdname clusterAssignment
#' @export
setMethod("clusterMap", "TimeSeriesExperiment", function(object) {
    if(!"final_cluster_map" %in% names(object@clusterAssignment)) {
        NULL
    }
    else return(object@clusterAssignment[["final_cluster_map"]])
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @title Differential Expression for \code{TimeSeriesExperiment}
#'
#' @description Getter method for \code{differentialExpression} slot 
#' of a \code{TimeSeriesExperiment}. The slot is a list with differential 
#' expression results possibly containing elements named 'timepoint_de', and 
#' 'trajectory_de' computed with \link{timepointDE} and \link{trajectoryDE}
#' functions.
#'
#' @docType methods
#' @rdname differentialExpression
#'
#' @param object a \code{TimeSeriesExperiment} object.
#' @param name one of elements of 'differentialExpression' slot: 
#' 'timepoint_de', 'trajectory_de'. If NULL, all elements are returned.
#' @param ... argiments to other functions.
#'
#' @return a \code{data.frame} or a list of \code{data.frame}s
#' @export
#' @examples
#' data("endoderm_small")
#' endoderm_small <- trajectoryDE(endoderm_small)
#' head(differentialExpression(endoderm_small, "trajectory_de"))
#'
#' @export
setGeneric("differentialExpression", function(object, ...) {
    standardGeneric("differentialExpression")})

#' @rdname differentialExpression
#' @export
setMethod(
    "differentialExpression", "TimeSeriesExperiment", 
    function(object, name = NULL) {
        if(is.null(name)) return(object@differentialExpression)
        else return(object@differentialExpression[[name]])
    }
)


### ===========================================================================
### Show method for TimeSeriesExperiment
### ---------------------------------------------------------------------------
#' @title show method for \code{TimeSeriesExperiment}
#' @docType methods
#' @rdname show-methods
#' @aliases show,TimeSeriesExperiment-method
#' @param object A TimeSeriesExperiment object
#' @return nothing, just prints to console
#' @importMethodsFrom SummarizedExperiment show
#' @export
setMethod("show", "TimeSeriesExperiment", function(object) {
    callNextMethod()
    scat <- function(fmt, vals=character(), exdent=2, ...)
    {
      vals <- ifelse(nzchar(vals), vals, "''")
      lbls <- paste(S4Vectors:::selectSome(vals), collapse=" ")
      txt <- sprintf(fmt, length(vals), lbls)
      cat(strwrap(txt, exdent=exdent, ...), sep="\n")
    }
    cat("========== \n")
    ## timepoints
    dlen <- length(timepoints(object))
    if (dlen) scat("timepoints(%d): %s\n", timepoints(object))
    else scat("timepoint: NULL\n")
    
    ## groups
    dlen <- length(groups(object))
    if (dlen) scat("groups(%d): %s\n", groups(object))
    else scat("groups: NULL\n")
    
    ## replicates
    dlen <- length(replicates(object))
    if (dlen) scat("replicates(%d): %s\n", replicates(object))
    else scat("replicates: NULL\n")
    
    cat("----- \n")
    # collapsed data
    dlen <- dim(assayCollapsed(object))
    if (dlen[[1]]) cat("assayCollapsed dim:", dlen[1], dlen[2], "\n")
    else scat("assayCollapsed: NULL\n")
    
    # collapsed colData
    dlen <- dim(colDataCollapsed(object))
    if (dlen[[1]]) scat("colDataCollapsed names(%d): %s\n", 
                        colnames(colDataCollapsed(object)))
    else scat("colDataCollapsed: NULL\n")
    
    # time-series
    dlen <- length(timeSeries(object))
    if (dlen) scat("timeSeries(%d): %s\n", names(timeSeries(object)))
    else scat("timeSeries: NULL\n")
    
    # dimensionality reduction results
    dlen <- length(dimensionReduction(object))
    if (dlen) scat("dimensionReduction(%d): %s\n", 
                   names(dimensionReduction(object)))
    else scat("dimensionReduction: NULL\n")
    
    # clusterin results for features (rows)
    dlen <- length(clusterAssignment(object))
    if (dlen) scat("clusterAssignment(%d): %s\n", 
                   names(clusterAssignment(object)))
    else scat("clusterRows: NULL\n")
    
    # differential expression results
    dlen <- length(differentialExpression(object))
    if (dlen) scat("differentialExpression(%d): %s\n", 
                   names(differentialExpression(object)))
    else scat("differentialExpression: NULL\n")
})













