
reset_results <- function(object){
    object@sample.data.collapsed <- data.frame()
    object@data.collapsed <- NULL
    object@timecourse.data = list()
    object@dim.red = list()
    object@cluster.features = list()
    object@diff.expr = list()
    return(object)
}


replace_names <- function(query, orig, replace) {
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

################################################################################

#' show method for \code{vistimeseq}
#'
#' @docType methods
#' @name show
#' @rdname show-methods
#' @aliases show,vistimeseq-method
#'
#' @param object A vistimeseq object
#'
#' @return nothing, just prints to console
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

##############################################################################

#' Number of samples and features
#'
#' Returns number of samples or features in the data.
#'
#' @docType methods
#' @name n_samples
#' @rdname counters
#'
#' @param object a \code{vistimeseq} object.
#'
#' @return an integer
#' @export
#' @examples
#'
#' endoderm_small
#' n_samples(endoderm_small)
#'
setGeneric("n_samples", function(object) standardGeneric("n_samples"))

#' @rdname counters
#' @aliases n_samples
setMethod(f = "n_samples", signature = "vistimeseq",
          function(object) { return(length(object@sample.names)) }

)

#' @docType methods
#' @name n_features
#' @rdname counters
#'
#' @export
#' @examples
#' endoderm_small
#' n_features(endoderm_small)
#'
setGeneric("n_features", function(object) standardGeneric("n_features"))

#' @rdname counters
#' @aliases n_features
setMethod(f = "n_features", signature = "vistimeseq",
          function(object) { return(length(object@feature.names)) }

)

##############################################################################

#' Accessors for the 'sample.names' slot of a \code{vistimeseq} object.
#'
#' 'sample.names' slots holds sample names
#'
#' @docType methods
#' @name sample_names
#' @rdname sample_names
#'
#' @param object a \code{vistimeseq} object.
#' @param value a character vector with new sample names
#'
#' @return a character vector
#' @export
#' @examples
#' endoderm_small
#' head(sample_names(endoderm_small))
#'
setGeneric("sample_names", function(object) standardGeneric("sample_names"))

#' @rdname sample_names
#' @aliases sample_names
setMethod(f = "sample_names",signature = "vistimeseq",
          function(object) {return(object@sample.names)}
)

set_sample_names <- function(object, value) {
    if(length(value) != length(object@sample.names)) {
        stop("substitute sample names length does not agree with object ",
                 "dimensions")
    }
    object@sample.names <- value
    colnames(object@data) <- colnames(object@raw.data) <- value
    rownames(object@sample.data) <- object@sample.data$sample <- value
    dim_red_names <- names(object@dim.red)
    if(is.null(dim_red_names)) return(object)
    for (name_mat in dim_red_names[grep("_sample", dim_red_names)]) {
        mat <- object@dim.red[[name_mat]]
        rownames(mat) <- value
        object@dim.red[[name_mat]] <- mat
    }
    validObject(object)
    return(object)
}

#' @rdname sample_names
setGeneric("sample_names<-", function(object, value)
    standardGeneric("sample_names<-"))

#' @rdname sample_names
#' @exportMethod "sample_names<-"
setReplaceMethod(f = "sample_names", signature = "vistimeseq",
                 definition = set_sample_names
)


#' Accessors for the 'sample.data' slot of a \code{vistimeseq} object.
#'
#' 'sample.data' slots holds information on individual samples including
#' corresponding sample name in column 'sample' as well as time, group
#' and replicate.
#'
#' @docType methods
#' @name sample_data
#' @rdname sample_data
#'
#' @param object a \code{vistimeseq} object
#' @param value a \code{data.frame} with new sample data
#'
#' @return a \code{data.frame}
#' @examples
#'
#' endoderm_small
#' head(sample_data(endoderm_small))
#'
#' @export
#'
setGeneric("sample_data", function(object) standardGeneric("sample_data"))

#' @rdname sample_data
#' @aliases sample_data
setMethod(f = "sample_data", signature = "vistimeseq",
          function(object) {return(object@sample.data)}
)

#' @rdname sample_data
setGeneric("sample_data<-", function(object, value)
    standardGeneric("sample_data<-"))

set_sample_data <- function(object, value) {
    if(!all(c("sample", "group", "time", "replicate") %in% colnames(value))) {
        stop("sample data must contain columns ", 
             "\"sample\", \"group\", \"time\", and \"replicate\".")
    }
    if(!all(value$sample == object@sample.names)) {
        object <- set_sample_names(object, value$sample)
    }
    object@sample.data <- value

    reset <- FALSE
    if(!all(value$group == object@sample.data$group)){
        message("new group assignment, resetting all affected results.")
        object@group <-value$group
        reset <- TRUE
    }
    if(!all(value$replicate == object@sample.data$replicate)){
        message("new replicate assignment, resetting all affected results.")
        object@replicate <-value$replicate
        reset <- TRUE
    }
    if(!all(value$time == object@sample.data$time)){
        message("new time assignment, resetting all affected results.")
        object@time <-value$time
        reset <- TRUE

    }
    if(reset) {
        object <- reset_results(object)
    }
    validObject(object)
    return(object)
}

#' @rdname sample_data
#' @exportMethod "sample_data<-"
setReplaceMethod(f = "sample_data", signature = "vistimeseq",
                 definition = set_sample_data
)

################################################################################

#' Accessors for the 'feature.names' slot of a \code{vistimeseq} object.
#'
#' 'feature.names' slots holds feature names
#'
#' @docType methods
#' @name feature_names
#' @rdname feature_names
#'
#' @param object a \code{vistimeseq} object.
#' @param value a character vectore with new feature names
#'
#' @return a character vector
#' @export
#' @examples
#' endoderm_small
#' head(feature_names(endoderm_small))
#'
setGeneric("feature_names", function(object)
    standardGeneric("feature_names"))

#' @rdname feature_names
#' @aliases feature_names
setMethod(f = "feature_names", signature = "vistimeseq",
          function(object) {return(object@feature.names)}
)

#' @rdname feature_names
#' @export
setGeneric("feature_names<-", function(object, value)
    standardGeneric("feature_names<-"))

set_feature_names <- function(object, value) {
    if(length(value) != length(object@feature.names)) {
        stop("substitute feature names length does not agree with object",
                 " dimensions")
    }
    object@feature.names <- value
    rownames(object@data) <- rownames(object@raw.data) <- value
    rownames(object@feature.data) <- object@feature.data$feature <- value

    if(!is.null(object@data.collapsed)) {
        rownames(object@data.collapsed) <- value
    }
    if(!is.null(object@timecourse.data)) {
        for(i in seq_along(object@timecourse.data)){
            object@timecourse.data[[i]]$feature <-
                replace_names(object@timecourse.data[[i]]$feature,
                              orig = object@feature.names,
                              replace = value)
        }
    }
    if(!is.null(object@cluster.features)){
        object@cluster.features$cluster_map$feature <-
            replace_names(object@cluster.features$cluster_map$feature,
                          orig = object@feature.names,
                          replace = value)
        object@cluster.features$hclust$labels <-
            replace_names(object@cluster.features$hclust$labels,
                          orig = object@feature.names,
                          replace = value)
    }

    dim_red_names <- names(object@dim.red)
    for (name_mat in dim_red_names[grep("_feature")]) {
        mat <- object@dim.red[[name_mat]]
        rownames(mat) <- value
        object@dim.red[[name_mat]] <- mat
    }

    if(!is.null(object@diff.expr$timepoint_de)){
        for(tmp in names(object@diff.expr$timepoint_de)){
            object@diff.expr$timepoint_de[[tmp]]$feature <-
                replace_names(object@diff.expr$timepoint_de[[tmp]]$feature,
                              orig = object@feature.names,
                              replace = value)
        }
    }

    if(!is.null(object@diff.expr$trajectory_de)){
        object@diff.expr$trajectory_de$feature <-
            replace_names(object@diff.expr$trajectory_de$feature,
                          orig = object@feature.names,
                          replace = value)
    }
    validObject(object)
    return(object)
}

#' @rdname feature_names
#' @exportMethod "feature_names<-"
setReplaceMethod(f = "feature_names", signature = "vistimeseq",
                 definition = set_feature_names
)


#' Accessors for the 'feature.data' slot of a \code{vistimeseq} object.
#'
#' 'feature.data' slots holds information on individual features including
#' corresponding feature name in column 'feature'.
#'
#' @docType methods
#' @name feature_data
#' @rdname feature_data
#'
#' @param object a \code{vistimeseq} object.
#' @param value a \code{data.frame} with new feature data
#'
#' @return a \code{data.frame}
#' @export
#' @examples
#' endoderm_small
#' head(feature_data(endoderm_small))
#'
setGeneric("feature_data", function(object) standardGeneric("feature_data"))

#' @rdname feature_data
#' @aliases feature_data
setMethod(f = "feature_data", signature = "vistimeseq",
          function(object) {return(object@feature.data)}
)

#' @rdname feature_data
#' @export
setGeneric("feature_data<-", function(object, value)
    standardGeneric("feature_data<-"))

set_feature_data <- function(object, value) {
    if(!"feature" %in% colnames(value)) {
        stop("feature data must contain columns \"feature\".")
    }
    if(!all(value$feature == object@feature.names)) {
        object <- set_feature_names(object, value$feature)
    }
    object@feature.data <- value
    return(object)
}

#' @rdname feature_data
#' @exportMethod "feature_data<-"
setReplaceMethod(f = "feature_data", signature = "vistimeseq",
                 definition = set_feature_data
)

################################################################################

#' Accessors for the 'group' slot of a \code{vistimeseq} object.
#'
#' 'grop' slots holds feature names
#'
#' @docType methods
#' @name get_group
#' @rdname group_data
#'
#' @param object a \code{vistimeseq} object.
#' @param value a character vectore with new group membership
#'
#' @return a character vector
#' @export
#' @examples
#' endoderm_small
#' head(get_group(endoderm_small))
#'
setGeneric("get_group", function(object) standardGeneric("get_group"))

#' @rdname group_data
#' @aliases get_group
setMethod(f = "get_group", signature = "vistimeseq",
          function(object) {return(object@group)}
)


set_group <- function(object, value) {
    if(length(value) != length(object@group)) {
        stop("substitute groups length does not agree with object",
             " dimensions")
    }
    object@group <- value
    message("new group assignment, resetting all affected results.")
    object <- reset_results(object)
    return(object)
}

#' @rdname group_data
#' @export
setGeneric("set_group<-", function(object, value) {
    standardGeneric("set_group<-")})

#' @rdname group_data
#' @exportMethod "set_group<-"
setReplaceMethod(f = "set_group", signature = "vistimeseq",
                 definition = set_group
)

################################################################################

#' @rdname replicate_data
#' @export
setGeneric("set_replicate<-", function(object, value)
    standardGeneric("set_replicate<-"))

#' Accessors for the 'replicate' slot of a \code{vistimeseq} object.
#'
#' 'grop' slots holds feature names
#'
#' @docType methods
#' @name get_replicate
#' @rdname replicate_data
#'
#' @param object a \code{vistimeseq} object.
#' @param value a character vector with new replicate ids.
#'
#' @return a character vector
#' @export
#' @examples
#' endoderm_small
#' head(get_replicate(endoderm_small))
#'
setGeneric("get_replicate", function(object) standardGeneric("get_replicate"))

#' @rdname replicate_data
#' @aliases get_replicate
setMethod(f = "get_replicate", signature = "vistimeseq",
          function(object) {return(object@replicate)}
)


set_replicate <- function(object, value) {
    if(length(value) != length(object@replicate)) {
        stop("substitute replicates length does not agree with object ",
             "dimensions")
    }
    object@replicate <- value
    message("new replicate assignment, resetting all affected results.")
    object <- reset_results(object)
    return(object)
}

#' @rdname replicate_data
#' @exportMethod "set_replicate<-"
setReplaceMethod(f = "set_replicate", signature = "vistimeseq",
                 definition = set_replicate
)

################################################################################

#' Accessors for the 'time' slot of a \code{vistimeseq} object.
#'
#' 'grop' slots holds feature names
#'
#' @docType methods
#' @name get_time
#' @rdname time_data
#'
#' @param object a \code{vistimeseq} object.
#' @param value a numeric vector with new time information.
#'
#' @return a numeric vector
#' @export
#' @examples
#' endoderm_small
#' head(get_time(endoderm_small))
#'
setGeneric("get_time", function(object) standardGeneric("get_time"))

#' @rdname time_data
#' @aliases get_time
setMethod(f = "get_time", signature = "vistimeseq",
          function(object) {return(object@time)}
)


set_time <- function(object, value) {
    if(length(value) != length(object@time)) {
        stop("substitute times length does not agree with object",
                 " dimensions")
    }
    if(!all(is.numeric(value))) {
        stop("time must be numeric values.")
    }
    object@time <- value
    message("new time assignment, resetting all affected results.")
    object <- reset_results(object)
    return(object)
}

#' @rdname time_data
#' @export
setGeneric("set_time<-", function(object, value) standardGeneric("set_time<-"))


#' @rdname time_data
#' @exportMethod "set_time<-"
setReplaceMethod(f = "set_time", signature = "vistimeseq",
                 definition = set_time
)

################################################################################

#' @title Accessors for the 'raw.data' and 'data' slots of a \code{vistimeseq}
#' object.
#'
#' @description  The 'raw.data' of 'data' slots hold data where columns are
#' samples taken at specific timepoints and rows are features. If normalization
#' was performed 'data' holds normalized data
#'
#' @docType methods
#' @name get_data
#' @rdname vistimeseq_data
#'
#' @param object a \code{vistimeseq} object.
#' @param raw whether raw data should be returned
#' @param value a \code{data.frame}
#'
#' @return a \code{data.frame}
#'
#' @export
#' @examples
#' endoderm_small
#' head(get_data(endoderm_small))
#'
setGeneric("get_data", 
           function(object, raw = FALSE) standardGeneric("get_data"))

#' @rdname vistimeseq_data
#' @aliases get_data
setMethod(f = "get_data", signature = "vistimeseq",
          definition = function(object, raw = FALSE) {
              if(raw) {return(object@raw.data)}
              return(object@data)
          })

set_vistimeseq_data <- function(object, value, raw) {
    if(raw){
        object@raw.data <- value
    }
    object@data <- value
    object <- reset_results(object)
    feature_names(object) <- rownames(value)
    sample_names(object) <- colnames(value)

    validObject(object)
    return(object)
}

#' @rdname vistimeseq_data
#' @export
setGeneric("set_raw_data<-", function(object, value)
    standardGeneric("set_raw_data<-"))

#' @rdname vistimeseq_data
#' @exportMethod "set_raw_data<-"
setReplaceMethod(f = "set_raw_data", signature = "vistimeseq",
                 definition = function(object, value) {
                    set_vistimeseq_data(object, value, TRUE)
                 })

#' @rdname vistimeseq_data
#' @export
setGeneric("set_data<-", function(object, value) standardGeneric("set_data<-"))

#' @rdname vistimeseq_data
#' @exportMethod "set_data<-"
setReplaceMethod(f = "set_data", signature = "vistimeseq",
                 definition = function(object, value) {
                     set_vistimeseq_data(object, value, FALSE)
               })

################################################################################

#' @title Accessors for the 'data.collapsed' slot of a \code{vistimeseq} object.
#'
#' @description 'data.collapsed' slots holds dataset collapsed by replicates
#'
#' @docType methods
#' @name collapsed_data
#' @rdname collapsed_data
#'
#' @param object a \code{vistimeseq} object.
#'
#' @return a \code{data.frame}
#' @export
#' @examples
#' endoderm_small
#' endoderm_small <- collapse_replicates(endoderm_small)
#' head(collapsed_data(endoderm_small))
#'
setGeneric("collapsed_data", function(object) standardGeneric("collapsed_data"))


#' @rdname collapsed_data
#' @aliases collapsed_data
setMethod(f = "collapsed_data", signature = "vistimeseq",
          function(object) {return(object@data.collapsed)}
)

#' @title Accessors for the 'data.collapsed' slot of a \code{vistimeseq} object.
#'
#' @description 'sample.data.collapsed' slots holds sample data for data
#' collapsed by replicates
#'
#' @docType methods
#' @name collapsed_sample_data
#' @rdname sample_collapsed_data
#'
#' @param object a \code{vistimeseq} object.
#'
#' @return a \code{data.frame}
#' @export
#' @examples
#' endoderm_small
#' endoderm_small <- collapse_replicates(endoderm_small)
#' head(collapsed_sample_data(endoderm_small))
#'
setGeneric("collapsed_sample_data", function(object)
    standardGeneric("collapsed_sample_data"))

#' @rdname sample_collapsed_data
#' @aliases collapsed_sample_data
setMethod(f = "collapsed_sample_data", signature = "vistimeseq",
          function(object) {return(object@sample.data.collapsed)}
)

################################################################################
#' Accessors for the 'timecourse.data' slot of a \code{vistimeseq} object.
#'
#' 'timecourse.data' slots is a list with 'tc' and (optionally) 'tc_collapsed'
#' storing timce-course format of the data.
#'
#' @docType methods
#' @name time_course
#' @rdname time_course
#'
#' @param object a \code{vistimeseq} object.
#' @param collapsed whether collapsed time-course should be returned
#'
#' @return a \code{data.frame}
#'
#' @export
#' @examples
#' endoderm_small
#' endoderm_small <- convert_to_timecourse(endoderm_small)
#' head(time_course(endoderm_small))
#' head(time_course(endoderm_small, collapsed = TRUE))
#'
setGeneric("time_course", function(object, collapsed = FALSE)
    standardGeneric("time_course"))

#' @rdname time_course
#' @aliases time_course
setMethod(f = "time_course", signature = "vistimeseq",
          function(object, collapsed = FALSE) {
              if(collapsed) {
                  return(object@timecourse.data$tc_collapsed)
              }
              return(object@timecourse.data$tc)}
)

################################################################################

#' Accessors for the 'dim.red' slot of a \code{vistimeseq} object.
#'
#' 'dim.red' slots is a list with 'pca_sample', 'pca_feature' and
#' 'pca_eigs' results from PCA projection.
#'
#' @docType methods
#' @name get_dim_reduced
#' @rdname dim_reduce
#'
#' @param object a \code{vistimeseq} object.
#' @param type one of elements of 'dim.red' slot: 'pca_sample', 'pca_feature'
#' and 'pca_eigs' or 'all' for returning the entire list.
#'
#' @return a \code{data.frame} or a list of \code{data.frame}s
#' @export
#' @examples
#' endoderm_small
#' endoderm_small <- run_pca(endoderm_small)
#' head(get_dim_reduced(endoderm_small, "pca_sample")[, 1:3])
#'
setGeneric("get_dim_reduced", function(object, type = "all")
    standardGeneric("get_dim_reduced"))


#' @rdname dim_reduce
#' @aliases get_dim_reduced
setMethod(f = "get_dim_reduced", signature = "vistimeseq",
          function(object, type = "all") {
              if(!type %in% c('all', 'pca_sample', 'pca_feature', 'pca_eigs')){
                  stop("type must be one of 'all', 'pca_sample', 'pca_feature'",
                           "and 'pca_eigs'.")
              }
              if(type == "all") {
                  return(object@dim.red)
              }
              if(type == "pca_sample") {
                  return(object@dim.red$pca_sample)
              }
              if(type == "pca_feature") {
                  return(object@dim.red$pca_feature)
              }
              if(type == "pca_eigs") {
                  return(object@dim.red$pca_eigs)
              }
          }
)

################################################################################

#' Accessors for the "cluster_map" element of 'cluster.features' slot of
#' a \code{vistimeseq} object.
#'
#' 'cluster.features' slots is a list with '"hclust", "cluster_map"
#' where first is an hclust object and second is a \code{data.frame} with 
#' cluster assignement, both computed with \code{cluster_timecourse_features()}
#' function.
#'
#' @docType methods
#' @name get_cluster_map
#' @rdname feat_clust
#'
#' @param object a \code{vistimeseq} object.
#'
#' @return a \code{data.frame}
#' @export
#' @examples
#'
#' endoderm_small
#' endoderm_small <- cluster_timecourse_features(endoderm_small)
#' head(get_cluster_map(endoderm_small))
#'
setGeneric("get_cluster_map", function(object)
    standardGeneric("get_cluster_map"))

#' @rdname feat_clust
#' @aliases feat_clust
setMethod(f = "get_cluster_map", signature = "vistimeseq",
          function(object) { return(object@cluster.features$cluster_map) }

)

#' Accessors for the "hclust" element of 'cluster.features' slot of
#' a \code{vistimeseq} object.
#'
#' 'cluster.features' slots is a list with '"hclust", "cluster_map"
#' where first is an hclust object and second is a \code{data.frame} with
#' cluster assignement, both computed with \code{cluster_timecourse_features()}
#' function.
#'
#' @docType methods
#' @name get_cluster_hclust
#' @rdname feat_clust
#'
#' @return an hclust object
#' @export
#' @examples
#'endoderm_small
#'endoderm_small <- cluster_timecourse_features(endoderm_small)
#'plot(get_cluster_hclust(endoderm_small), labels = FALSE, xlab = "genes",
#'       sub = "")
#'
setGeneric("get_cluster_hclust", function(object)
    standardGeneric("get_cluster_hclust"))

#' @rdname feat_clust
#' @aliases get_cluster_hclust
setMethod(f = "get_cluster_hclust", signature = "vistimeseq",
          function(object) { return(object@cluster.features$hclust) }

)


################################################################################

#' Accessors for the 'diff.expr' slot of a \code{vistimeseq} object.
#'
#' 'diff.expr' slots is a list with differential expression results
#' possibly containing elements named 'timepoint_de', and 'trajectory_de'.
#'
#' @docType methods
#' @name get_diff_expr
#' @rdname diff_expr
#'
#' @param object a \code{vistimeseq} object.
#' @param type one of elements of 'dim.red' slot: 'timepoint_de',
#' 'trajectory_de' or 'all' for returning the entire list.
#'
#' @return a \code{data.frame} or a list of \code{data.frame}s
#' @export
#' @examples
#' endoderm_small <- trajectory_de(endoderm_small)
#' head(get_diff_expr(endoderm_small, "trajectory_de"))
#'
setGeneric("get_diff_expr", function(object, type = "all")
    standardGeneric("get_diff_expr"))

#' @rdname diff_expr
#' @aliases get_diff_expr
setMethod(f = "get_diff_expr", signature = "vistimeseq",
          function(object, type = "all") {
              if(!type %in% c('all', 'timepoint_de', 'trajectory_de')){
                  stop("type must be one of 'all', 'timepoint_de', or
                           'trajectory_de'.")
              }
              if(type == "all") {
                  return(object@diff.expr)
              }
              if(type == "timepoint_de") {
                  return(object@diff.expr$timepoint_de)
              }
              if(type == "trajectory_de") {
                  return(object@diff.expr$trajectory_de)
              }
          }
)

################################################################################

#' @title Project name
#'
#' @description Returns project name.
#'
#' @docType methods
#' @name project_name
#' @rdname proj_name
#'
#' @param object a \code{vistimeseq} object.
#' @param value a character string
#'
#' @return a character string
#' @export
#' @examples
#' endoderm_small
#' project_name(endoderm_small)
#'
setGeneric("project_name", function(object) standardGeneric("project_name"))

#' @rdname proj_name
#' @aliases project_name
setMethod(f = "project_name", signature = "vistimeseq",
          function(object) { return(object@project.name) }

)

#' @rdname proj_name
#' @export
setGeneric("project_name<-", function(object, value)
    standardGeneric("project_name<-"))

#' @rdname proj_name
#' @exportMethod "set_group<-"
setReplaceMethod(f = "project_name", signature = "vistimeseq",
                 definition = function(object, value) {
                     return(object@project.name <- value) }
)

