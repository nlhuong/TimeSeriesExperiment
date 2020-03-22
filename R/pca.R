#' @title Run PCA scores.
#'
#' @description Compute the projection of either samples or time series
#' features onto a PCA space. In case, of PCA for features, the PCA can be
#' computed for individual sample group as indicated or for samples from all
#' groups. In either, case the data is first collapsed over replicates, so
#' that each gene is represented as a vector of a single time course.
#'
#' @param object A \code{TimeSeriesExperiment} object
#' @param collapse.replicates Whether PCA should be computed on the data
#' with replicates aggregated.
#' @param groups.selected An optional character string indicating a particular
#' group of samples PCA should be applied to. By default set to NULL and all 
#' groups are included.
#' @param var.stabilize.method Method for variance stabilization (VST).
#' Currently, supports "none" (no VST), "log1p" (log plus one), "asinh"
#' (inverse hyperbolic sine) or "deseq"
#' (\code{\link[DESeq2]{varianceStabilizingTransformation}} function from
#' \code{DESeq2} package). Default is "log1p".
#'
#' @return Returns \code{TimeSeriesExperiment} object with PCA results 
#' in \code{dim.red} slot, a lists containing matrices of coordinates 
#' 'pca_sample', and 'pca_features', as well as a vector 'pca_eigs'.
#'
#' @importFrom DESeq2 varianceStabilizingTransformation
#' @importFrom dplyr select summarize_all group_by
#' @importFrom stats prcomp
#' @importFrom  methods slot<-
#' @importFrom methods validObject
#' @importFrom SummarizedExperiment assays
#' @export
#' @examples
#' data("endoderm_small")
#' endoderm_small <- runPCA(endoderm_small)
#' head(dimensionReduction(endoderm_small, "pca_sample")[, 1:5])
#' head(dimensionReduction(endoderm_small, "pca_eigs"))
#'
runPCA <- function(object, collapse.replicates = FALSE, 
                   groups.selected = NULL, var.stabilize.method = "log1p") 
{
    if (!is(object, "TimeSeriesExperiment")) 
        stop("Input must be a 'TimeSeriesExperiment' object.")
    if (!validObject(object)) 
        stop("Invalid TimeSeriesExperiment object.")
    if (all(!is.null(groups.selected),
            !any(groups.selected %in% groups(object))))
        stop("One of the groups: ", groups.selected, ", is not in the data.")
    
    if(!"norm" %in% names(assays(object))) {
        object <- normalizeData(object)
    }
    if(all(collapse.replicates, nrow(assayCollapsed(object)) == 0)){
        object <- collapseReplicates(object)
    }

    if(!collapse.replicates){
        X <- as.matrix(assays(object)$norm)
        if (!is.null(groups.selected)) {
            X <- X[, groups(object) %in% groups.selected]
        }
    } else {
        X <- as.matrix(assayCollapsed(object))
        if (!is.null(groups.selected)) {
            X <- X[, colDataCollapsed(object)$group %in% groups.selected]
        }
    }
    X <- t(varianceStabilization(X, var.stabilize.method))
    pca.res <- prcomp(X)
    dim.red <- dimensionReduction(object)
    dim.red$pca_eigs <- pca.res$sdev^2
    names(dim.red$pca_eigs) <- paste0("eig_", seq_along(pca.res$sdev))
    dim.red$pca_sample <- pca.res$x
    dim.red$pca_feature <- pca.res$rotation
    slot(object, name = "dimensionReduction", check = TRUE) <- dim.red
    return(object)
}
