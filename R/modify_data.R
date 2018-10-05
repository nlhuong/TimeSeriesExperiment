#' @title Add Featuredata
#'
#' @description Adds additional data for samples to \code{vistimeseq} object.
#' Can be any piece of information associated with a sample (e.g. subject
#' clinical or nonclinical data). The additional data can be used in
#' plotting functions.
#'
#' @param object vistimeseq object
#' @param sampledata Data frame where the row names are sample names
#' and the columns are additional sample attributes
#' @param col.name Column names in sampledata to add to the object
#'
#' @return \code{vistimeseq} object where the additional sample data has been
#' added as columns in object@@sample.data
#'
#' @importFrom methods validObject
#' @export
#' @examples
#' endoderm_small
#' idx <- sample(1:26, length(sample_names(endoderm_small)), replace = TRUE)
#' random_letters <- data.frame(
#'   letters = LETTERS[idx],
#'   row.names = sample_names(endoderm_small)
#'   )
#' endoderm_small <- add_sample_data(
#'   object = endoderm_small,
#'   sampledata = random_letters
#'   )
#' head(sample_data(endoderm_small))
#'
add_sample_data <- function (
    object, sampledata, col.name = colnames(sampledata)) {
    if (!validObject(object)){
        stop("Invalid 'vistimeseq' object.")
    }
    if(nrow(sampledata) != length(sample_names(object))) {
        stop("new sample data dimension does not agree ",
             "with object@sample.data.")
    }
    smp <- sample_data(object)
    intersect_cols <- intersect(col.name, colnames(smp))
    smp[, intersect_cols] <- sampledata[, intersect_cols, drop = FALSE]
    smp2 <- sampledata[, setdiff(col.name, colnames(smp)), drop = FALSE]
    smp <- data.frame(smp, smp2)
    sample_data(object) <- smp
    return(object)
}


#' @title Add Featuredata
#'
#' @description Adds additional data for features to \code{vistimeseq}  object.
#' Can be any piece of information associated with a feature (e.g. new
#' featureIDs, prevalence, total counts). The additional data can be used in
#' plotting functions.
#'
#' @param object vistimeseq object
#' @param featuredata Data frame where the row names are feature names
#' and the columns are additional feature attributes
#' @param col.name Column names in featuredata to add to the object
#'
#' @return Seurat object where the additional feature data has been added
#' as columns in object@@feature.data
#' @importFrom methods validObject
#' @export
#' @examples
#' idx <- sample(1:26, n_features(endoderm_small), replace = TRUE)
#' random_letters <- data.frame(
#'    letters = LETTERS[idx],
#'    row.names = feature_names(endoderm_small)
#' )
#' endoderm_small <- add_feature_data(
#'     object = endoderm_small,
#'     featuredata = random_letters
#' )
#' head(feature_data(endoderm_small))
#'
add_feature_data <- function(
    object, featuredata, col.name = colnames(featuredata)) {
    if (!validObject(object)){
        stop("Invalid vistimeseq object.")
    }
    if(nrow(featuredata) != length(feature_names(object))) {
        stop("new feature data dimension does not agree",
             " with object@feature.data")
    }
    fdata <- feature_data(object)
    intersect_cols <- intersect(col.name, colnames(fdata))
    fdata[, intersect_cols] <- featuredata[, intersect_cols, drop = FALSE]
    fdata2 <- featuredata[, setdiff(col.name, colnames(fdata)), drop = FALSE]
    fdata <- data.frame(fdata, fdata2)
    feature_data(object) <- fdata
    return(object)
}


#' @title Filter features
#'
#' @description Filters "vistimeseq" object to keep only chosen features.
#' All relevant slots are updates.
#'
#' @param object vistimeseq object
#' @param features features (genes) to keep
#'
#' @return "vistimeseq" object
#' @importFrom  dplyr filter
#' @importFrom methods validObject new
#' @export
#' @examples
#' features <- 1:100
#' endoderm_small
#' endoderm_small <- filter_features(endoderm_small, features)
#' endoderm_small
#'
filter_features <- function (object, features) {
    feature <- NULL
    if (!validObject(object)){
        stop("Invalid 'vistimeseq' object.")
    }
    if (all(is.numeric(features),
                    !all(features %in% seq_len(n_features(object))))){
        stop("Some \"features\" not included in \"vistimeseq\" object.")
    }
    if(is.numeric(features)){
        features <- feature_names(object)[features]
    }
    if (!all(features %in% feature_names(object))){
        stop("Some \"features\" not included in \"vistimeseq\" object.")
    }

    fltr_feature_data <- feature_data(object) %>% filter(feature %in% features)
    rownames(fltr_feature_data) <- fltr_feature_data$feature
    fltr_feature_data <- fltr_feature_data[features, , drop = FALSE]
    fltr_raw_data <- get_data(object, raw = TRUE)[features, , drop = FALSE]
    fltr_data <- get_data(object, raw = FALSE)[features, , drop = FALSE]

    fltr_collapsed_data <- NULL
    fltr_collapsed_sample_data <- data.frame()
    if(!is.null( collapsed_data(object) )) {
        fltr_collapsed_data <- collapsed_data(object)[features, ]
        fltr_collapsed_sample_data <- collapsed_sample_data(object)
    }

    fltr_timecourse <- list()
    if(!is.null(time_course(object))){
        fltr_timecourse[["tc"]] <- time_course(object) %>%
            filter(feature %in% features)
    }
    if(!is.null(time_course(object, collapsed = TRUE))){
        fltr_timecourse[["tc_collapsed"]] <-
            time_course(object, collapsed = TRUE) %>%
                filter(feature %in% features)
    }

    if(length(get_dim_reduced(object)) > 0 ) {
        message("Dimensionality reduction results are reset due to feature ",
                        "filtering. Need to be recomputed.")
    }
    if(length(get_cluster_map(object)) > 0) {
        message("Feature clustering results are reset due to feature ",
                        "filtering. Need to be recomputed.")
    }
    if(length(get_diff_expr(object)) > 0) {
        message("Differential expression results are reset due to feature ",
                        "filtering. Need to be recomputed.")
    }

    fltr_object <- new(
        Class = "vistimeseq",
        project.name = project_name(object),
        raw.data = fltr_raw_data,
        data = fltr_data,
        sample.names = sample_names(object),
        feature.names = features,
        sample.data = sample_data(object),
        feature.data = fltr_feature_data,
        data.collapsed = fltr_collapsed_data,
        sample.data.collapsed = fltr_collapsed_sample_data,
        group = get_group(object),
        replicate = get_replicate(object),
        time = get_time(object),
        timecourse.data = fltr_timecourse,
        dim.red = list(),
        cluster.features = list(),
        diff.expr = list()
    )
    validObject(object)
    return(fltr_object)
}
