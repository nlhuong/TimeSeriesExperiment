#' @title Add Featuredata
#'
#' @description Adds additional data for samples to the \code{vistimeseq} object.
#' Can be any piece of information associated with a sample (e.g. subject
#' clinical or nonclinical data). The additional data can be used in
#' plotting functions.
#'
#' @param object vistimeseq object
#' @param sampledata Data frame where the row names are sample names
#' and the columns are additional sample attributes
#' @param col.name Column names in sampledata to add to the object
#'
#' @return Seurat object where the additional sample data has been added as
#' columns in object@@sample.data
#'
#' @export
#'
#' @examples
#' idx <- sample(1:26, nrow(endoderm_small@raw.data), replace = TRUE)
#' random_letters <- data.frame(
#'   letters = LETTERS[idx],
#'   row.names = colnames(endoderm_small@raw.data)
#' )
#' endoderm_small <- add_sample_data(
#'    object = endoderm_small,
#'    sampledata = random_letters
#' )
#' head(endoderm_small@sample.data)
#'
add_sample_data <- function (object, sampledata,
                            col.name = colnames(sampledata)) {
  if (!validObject(object))
    stop("Invalid vistimeseq object.")
  cols.add <- intersect(col.name, colnames(sampledata))
  object@sample.data[, cols.add] <-
    sampledata[rownames(x = object@sample.data), cols.add]
  return(object)
}


#' @title Add Featuredata
#'
#' @description Adds additional data for features to the vistimeseq object.
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
#'
#' @examples
#' idx <- sample(1:26, nrow(endoderm_small@raw.data), replace = TRUE)
#' random_letters <- data.frame(
#'   letters = LETTERS[idx],
#'   row.names = rownames(endoderm_small@raw.data)
#'  )
#' endoderm_small <- add_feature_data(
#'    object = endoderm_small,
#'    featuredata = random_letters
#' )
#' head(endoderm_small@feature.data)
#'
add_feature_data <- function (object, featuredata,
                            col.name = colnames(featuredata)) {
  if (!validObject(object))
    stop("Invalid vistimeseq object.")
  cols.add <- intersect(col.name, colnames(featuredata))
  object@feature.data[, cols.add] <-
    featuredata[rownames(x = object@feature.data), cols.add]
  return(object)
}
