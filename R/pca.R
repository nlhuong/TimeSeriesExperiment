#' @title Run PCA scores.
#'
#' @description Compute the projection of either samples or time series
#' features onto a PCA space. In case, of PCA for features, the PCA can be
#' computed for individual sample group as indicated or for samples from all
#' groups. In either, case the data is first collapsed over replicates, so
#' that each gene is represented as a vector of a single time course.
#'
#' @param object A \code{vistimeseq} object
#' @param collapse.replicates Whether PCA should be computed on the data
#' with replicates aggregated.
#' @param group.selected An optional character string indicating a particular
#' group of samples PCA should be applied to. By default all groups are
#' included.
#' @param var.stabilize.method Method for variance stabilization (VST).
#' Currently, supports "none" (no VST), "log1p" (log plus one), "asinh"
#' (inverse hyperbolic sine) or "deseq"
#' (\code{\link[DESeq2]{varianceStabilizingTransformation}} function from
#' \code{DESeq2} package). Default is "log1p".
#'
#' @return Returns \code{vistimeseq} object with PCA results in \code{dim.red}
#' slot, a lists containing matrices of coordinates 'pca_sample', and
#' 'pca_features', as well as a vector 'pca_eigs'.
#'
#' @importFrom DESeq2 varianceStabilizingTransformation
#' @importFrom dplyr select summarize_all group_by
#' @importFrom tibble column_to_rownames
#' @export
#' @examples
#' endoderm_small
#' endoderm_small <- run_pca(endoderm_small)
#' head(endoderm_small@dim.red$pca_sample[, 1:5])
#' head(endoderm_small@dim.red$pca_eigs)
#'
run_pca <- function(object,
                    collapse.replicates = FALSE,
                    group.selected = NULL,
                    var.stabilize.method = "log1p") {
  if (!validObject(object)){
    stop("Invalid vistimeseq object.")
  }
  if (all(!is.null(group.selected), !group.selected %in% object@group)){
    stop("\"group.selected\", ", group.selected, ", is not in the data.")
  }
  if(is.null(object@data)) {
    message("First converting raw-data to CPMs.")
    object <- normalize_data(object)
  }
  if(all(collapse.replicates, is.null(object@data.collapsed))){
    message("First, collapsing over replicates.")
    object <- collapse_replicates(object)
  }

  if(!collapse.replicates){
    X <- as.matrix(object@data)
    if (!is.null(group.selected)) {
      X <- X[, object@group == group.selected]
    }
  } else {
    X <- as.matrix(object@data.collapsed)
    if (!is.null(group.selected)) {
      X <- X[, object@sample.data.collapsed$group == group.selected]
    }
  }
  X <- t(variance_stabilization(X, var.stabilize.method))
  pca.res <- prcomp(X)
  dim.red <- list()
  eigs <- pca.res$sdev^2
  names(eigs) <- paste0("eig_", 1:length(eigs))
  dim.red[["pca_eigs"]] <- eigs
  dim.red[["pca_sample"]] <- pca.res$x
  dim.red[["pca_feature"]] <- pca.res$rotation
  # dim.red[[paste0("pca_", type)]] <-
  #   list("pca.scores" = pca.scores, "pca.eigs" = pca.eigs)
  slot(object, name = "dim.red", check = TRUE) <- dim.red
  return(object)
}
