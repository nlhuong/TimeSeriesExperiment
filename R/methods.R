#' @title Run PCA scores.
#'
#' @description Compute the projection of either samples or time series
#' features onto a PCA space. In case, of PCA for features, the PCA can be
#' computed for individual sample group as indicated or for samples from all
#' groups. In either, case the data is first collapsed over replicates, so
#' that each gene is represented as a vector of a single time course.
#'
#' @param object A \code{timevis} object
#' @param type Either "sample" or "feature" indicating whether to compute PCA
#' for observations or for time series features.
#' @param group.selected An optional character string indicating which group
#' to subset the data to. By default all groups are considered. Only considered
#' for \code{type} = "feature" case.
#' @param var.stabilize.method Method for variance stabilization (VST)
#' if computing PCA for \code{type} = "sample". Currently, supports log plus
#' one normalization, "log", or DESeq2::varianceStabilizingTransformation,
#' "deseq". By default "lognorm" VST is performed.
#' @param pseudocount A numerical indicating the pseudocount in log
#' transformation
#' @param log.base A numerical indicating the base for log
#' transformation
#'
#' @return Returns \code{timevis} object with PCA results in the \code{dim.red}
#' slot. The PCA results are lists of 'pca.scores' and 'pca.eigs', and will
#' enter elements of \code{dim.red} slot named "pca_samples" or "pca_features".
#'
#' @importFrom DESeq2 varianceStabilizingTransformation
#' @importFrom dplyr select summarize_all group_by
#' @importFrom tibble column_to_rownames
#' @export
#' @examples
#' endoderm_small
#' endoderm_small <- normalize_data(endoderm_small)
#' endoderm_small <- run_pca(endoderm_small, type = "sample")
#' head(endoderm_small@dim.red$pca_sample[["pca.scores"]][, 1:5])
#' head(endoderm_small@dim.red$pca_sample[["pca.eigs"]])
#'
run_pca <- function(object, type, group.selected = NULL,
                    var.stabilize.method = "log",
                    pseudocount = 1, log.base = exp(1)) {
  if (!validObject(object))
    stop("Invalid timevis object.")
  if(is.null(object@data)) {
    stop("Data has not been normalized. Normalize data first.")
  }
  if (!type %in% c("sample", "feature")) {
    stop("The 'type' argument must be 'sample' or 'feature'.")
  }

  if (type == "sample") {
    X <- t(variance_stabilization(as.matrix(object@data),
                                  method = var.stabilize.method,
                                  pseudocount = pseudocount,
                                  log.base = log.base))
  } else if (type == "feature") {
    if (! "tc_collapsed" %in% names(object@timecourse.data)) {
      stop("No 'tc_collapsed' entry in object@timecourse.data.",
           " Use 'collapse_replicates()' and 'convert_to_timecourse()'",
           " first to collapse the data over replicates and convert to",
           " time-course format.")
    }
    X <- object@timecourse.data$tc_collapsed %>%
      select(-replicate)
    if (all(!is.null(group.selected), !group.selected %in% X$group)){
      stop("The 'group.selected', ", group.selected,
           ", is not inlcuded in the data.")
    }
    if (is.null(group.selected)) {
      X <- X %>% group_by(feature) %>%
        select(-group) %>%
        summarize_all(mean) %>%
        as.data.frame() %>%
        column_to_rownames("feature")
    } else {
      X <- X %>%
        filter(group == group.selected) %>%
        select(-group) %>%
        column_to_rownames("feature")
    }
  }
  pca.res <- prcomp(X)
  pca.eigs <- pca.res$sdev^2
  pca.scores <- pca.res$x
  dim.red <- list()
  dim.red[[paste0("pca_", type)]] <-
    list("pca.scores" = pca.scores, "pca.eigs" = pca.eigs)
  slot(object, name = "dim.red", check = TRUE) <- dim.red
  return(object)
}


#' @title Compute cluster assignment.
#'
#' @description This function computes the cluster assignment for hierarchical
#' clustering results using static or dynamics methods.
#'
#' @param hclst an object of class hclust representing the clustering of
#' features of X.
#' @param dynamic whether dynamic branch cutting should be done for cluster
#' assignment.
#' @param h a numeric scalar for the height of the tree if static cutting is
#' done for cluster assignment.
#' @param ... other parameters for dynamicTreeCut::cutreeDynamic().
#'
#' @return a matrix of cluster assignments.
#' @importFrom dynamicTreeCut cutreeDynamicTree
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join arrange
#' @export
assign_cluster <- function(hclst, dynamic = FALSE, h = 0.5, ...) {
  leaf_ix <- hclst$order
  leaf_order <- hclst$labels[leaf_ix]
  if(!dynamic) {
    cluster <- data.frame("cluster" = cutree(hclst, h = h))
  } else{
    cluster <- data.frame("cluster" = cutreeDynamicTree(hclst, ...))
    rownames(cluster) <- hclst$labels
  }
  cluster <- cluster %>%
    rownames_to_column(var = "feature") %>%
    left_join(data_frame(feature = leaf_order, leaf_ix = leaf_ix)) %>%
    arrange(leaf_ix)
  cluster$cluster <- paste0("C", cluster$cluster)
  return(cluster)
}


#' @title Cluster time series data.
#'
#' @description Find the cluster assignment for time series data.
#' If 'subset' is specified the clsutering is performed on the subset
#' of the data, and the rest of the rows are assigned based on the distances
#' to the centroids of computed clusters.
#'
#' @param X a data matrix or data.frame where rows are time series.
#' @param subset indices or rownames of rows to keep.
#' @param alpha_coeff the coefficients for the added diff/lag features, by default set
#' as in 'add_diff' function.
#' @param method the distance metric for the dissimilarity used for clustering.
#' @param ... other parameters for 'assign_cluster' function
#'
#' @return a list with the hclust object and the clust.assignment data.frame.
#'
#' @importFrom proxy dist
#' @importFrom dplyr select left_join group_by summarise_all
#' @export
cluster_data <- function(X, subset = 1:nrow(X), alpha_coeff = NULL,
                         method = "euclidean",  ...) {
  X.subset <- data.frame(X[subset, ], check.names = FALSE)
  # Add differences (lags)
  if(!is.null(alpha_coeff)) {
    X.subset <- add_diff(X.subset, alpha = alpha_coeff)
  }

  # Perform hierarchical clustering
  D <- dist(X.subset, method = method)
  hclust <- hclust(D)

  # Find clusters
  clst.mapping <- assign_cluster(hclust, ...)

  if(nrow(clst.mapping) < nrow(X)) {
    # Since, 0 in dynamicCutTree means the object was unassigned,
    # we remove these features.
    clst.mapping <- clst.mapping[clst.mapping$cluster != "C0", ]

    # Compute cluster centroids
    clst.centroids <- clst.mapping %>%
      select(feature, cluster) %>%
      left_join(X.subset %>% rownames_to_column("feature")) %>%
      select(-feature) %>%
      group_by(cluster) %>%
      summarise_all(mean) %>%
      data.frame(check.names = FALSE)

    # Assign the rest of the genes to the closest cluster centroid
    D <- proxy::dist(X, clst.centroids %>%  column_to_rownames("cluster"))
    clst.mapping <- apply(D, 1, function(x) {colnames(D)[which.min(x)]})
    clst.mapping <- data.frame(feature = names(clst.mapping),
                               cluster = clst.mapping,
                               stringsAsFactors = FALSE)
  }
  return(list(hclust = hclust, clst.assignment = clst.mapping))
}



#' @title Static differential expression testing.
#'
#' @description This is a wrapper around 'limma' + 'voom' functions for testing
#' differential expression. We are performing tests for each single timepoint.
#'
#' @param X a data matrix or data.frame where rows are time series.
#' @param sample_data a data.frame of sample data.
#' @param alpha a scalar for level of significance.
#' @param group_var a character scalar (one of the column names in
#' 'sample_data') indicating the group sample variable for which DE is tested.
#' @param test_coeff a character scalar for test coefficient.
#'
#' @return a list with the hclust object and the clust.assignment data.frame.
#'
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom limma voom lmFit eBayes topTable
#' @export
static_diff_expr <- function(X, sample_data, gene_data = NULL, alpha = 0.05,
                             group_var = NULL, test_coeff = NULL) {
  if(is.null(group_var))
    group_var <- colnames(sample_data)[1]

  sample_data$category <- sample_data[[group_var]]

  if(is.null(gene_data))
    gene_data <- data.frame(gene = rownames(X))

  limma.res <- lapply(unique(sample_data$time), function(t) {
    smp.t <- sample_data[sample_data$time == t, ]
    x.t <- X[, smp.t$sample]

    # Construct a DGEList object
    dge <- DGEList(counts = x.t,
                   genes = gene_data,
                   samples = smp.t)
    # Compute size factor (lib sizes)
    dge <- calcNormFactors(dge)

    # Construct a design model matrix
    dsgn <- model.matrix(~ category, dge$samples)

    # Voom normalization
    v <- voom(dge, dsgn, plot = FALSE)

    # Fit the model
    fit <- lmFit(v, dsgn)

    # Compute test statistics
    fit <- eBayes(fit)

    # Statistics table
    t.res.limma <- topTable(fit, adjust="BH", number = Inf,
                            p.value = 1, coef = test_coeff)
    return(t.res.limma)
  })
  names(limma.res) <- unique(sample_data$time)
  limma.signif <- sapply(limma.res, function(tab) {
    tab <- tab[, c("gene", "adj.P.Val")]
    tab <- tab[order(tab$adj.P.Val), ]
    tab[tab$adj.P.Val >= alpha, "gene"] <- NA
    return(tab$gene)
  })
  colnames(limma.signif) <- names(limma.res)
  return(list(limma.res.lst = limma.res, limma.signif.df  = limma.signif))
}


#' @title Differences in expression profiles.
#'
#' @description This is function computes differences in gene
#' expression profiles between two selected groups.
#'
#' @param X a data matrix or data.frame where rows are time series.
#' @param group a vector of length equal to ncol(X) indicating
#' the group membership of the sample.
#' @param replicate a vector of length equal to ncol(X) indicating
#' the sample replicate (e.g. individual).
#' @param time a vector of length equal to ncol(X) indicating the time
#' variable corresponding to the sample (data column).
#' @param reference_name a character scalar indicating the reference group.
#' @param test_group_name a character scalar indicating the test group.
#' @param alpha_coeff the weights for each lag difference, by default
#' 0.5^(1:(ncol(X)-1))
#'
#' @return a data.frame with columns 'gene', 'diff.mean', 'diff.sd',
#' 'ref.diff.mean', 'ref.diff.sd', 'test.diff.mean', 'test.diff.sd', indicating
#' respectively the mean and sd of the differences between time series
#' gene expression profiles of replicates between two groups, within the
#' reference, and within the test group.
#'
#' @importFrom dplyr filter
#' @export
time_series_diff <- function(X, time, replicate,
                             group = NULL,
                             reference_name = NULL,
                             test_group_name = NULL,
                             alpha_coeff = c(0.5, 0.25)) {
  if(is.null(group))
    group <- rep("G1", ncol(X))
  if(is.null(reference_name))
    reference_name <- unique(group)[1]
  if(is.null(test_group_name))
    test_group_name <- unique(group)[2] # This will be NA if only one group

  if(ncol(X) != length(group))
    stop("ncol(X) must match length(group)")
  if(ncol(X) != length(replicate))
    stop("ncol(X) must match length(replicate)")
  if(ncol(X) != length(time))
    stop("ncol(X) must match length(time)")
  if(!reference_name %in% group)
    stop("reference_name must be a member of group.")

  time <- as.numeric(time)
  X <- X[, order(replicate, time)]
  map.df <- data.frame(group, replicate)
  map.df <- map.df[!duplicated(map.df), ]

  Y <- lapply(map.df$replicate, function(lab) {
    y <- X[, replicate == lab]
    colnames(y) <- time[replicate == lab]
    y <- add_diff(y, alpha = alpha_coeff)
    df <- data.frame(y, check.names = FALSE)
    return(df)
  })
  names(Y) <- map.df$replicate

  # Distance between replicates within a reference group
  ref.lst <- Y[names(Y)[map.df$group == reference_name]]
  pairs.df <- expand.grid(1:length(ref.lst), 1:length(ref.lst)) %>%
    filter(Var1 < Var2)
  D.ref <- mapply(ref.lst[pairs.df$Var1], ref.lst[pairs.df$Var2],
                  FUN = function(x, y) {sqrt(rowSums((x - y)^2))})

  diff.df <- data.frame(gene = rownames(D.ref),
                        ref.diff.mean = rowMeans(D.ref),
                        ref.diff.sd = apply(D.ref, 1, sd),
                        stringsAsFactors = FALSE)
  if(test_group_name %in% group) {
    test.lst <- Y[names(Y)[map.df$group == test_group_name]]
    # Distance between replicates within a test group
    pairs.df <- expand.grid(1:length(test.lst), 1:length(test.lst)) %>%
      filter(Var1 < Var2)
    D.test <- mapply(test.lst[pairs.df$Var1], test.lst[pairs.df$Var2],
                     FUN = function(x, y) {sqrt(rowSums((x - y)^2))})

    # Distance between each replicate in a reference and replicate row in a test group
    pairs.df <- expand.grid(1:length(ref.lst), 1:length(test.lst))
    D <- mapply(ref.lst[pairs.df$Var1], test.lst[pairs.df$Var2],
                FUN = function(x, y) {sqrt(rowSums((x - y)^2))})


    diff.df <- data.frame(diff.df,
                          test.diff.mean = rowMeans(D.test),
                          test.diff.sd = apply(D.test, 1, sd),
                          diff.mean = rowMeans(D),
                          diff.sd = apply(D, 1, sd),
                          stringsAsFactors = FALSE)
  }

  return(diff.df)
}


#' @title Adonis for time series data.
#'
#' @description A function computes distances between time series data
#' (with added lag differences) and uses vegan::adonis() function
#' to test for group differences.
#'
#' @param X a data matrix or data.frame where rows are time series.
#' @param time a vector of length equal to ncol(X) indicating the time
#' variable corresponding to the sample (data column).
#' @param replicate a vector of length equal to ncol(X) indicating
#' the sample replicate (e.g. individual).
#' @param group a vector of length equal to ncol(X) indicating
#' the group membership of the sample.
#' @param alpha_coeff weights for each lag difference, by default NULL.
#' @param dist_method the name of any method used in vegdist to calculate
#' pairwise distances, "euclidean" by defaults.
#' @param p_adj_method a correction method. Can be abbreviated.
#' @param ... other options to adonis function.
#'
#' @return a data.frame with adonis results for all features.
#'
#' @importFrom vegan adonis
#' @importFrom dplyr filter rename arrange mutate
#' @importFrom tibble rownames_to_column
#' @export
ts_adonis <- function(X, time, replicate = NULL, group = NULL,
                      alpha_coeff = NULL, dist_method = "euclidean",
                      p_adj_method = "BH",
                      verbose = TRUE, ...) {

  time.ord <- sort(unique(time))
  feature.names <- rownames(X)

  X.ts <- data_to_ts(X, time = time,
                     replicate = replicate,
                     group = group)
  if(!is.null(alpha_coeff)) {
    X.diff <- add_diff(X.ts[, as.character(time.ord)],
                       alpha_coeff = alpha_coeff)
    X.ts <- cbind(X.ts[, 1:2], X.diff)
  }

  adonis.res <- lapply(1:length(feature.names), function(i) {

    if((i %% 500) == 0 && verbose) message("Feature ", i)

    ix.ts <- X.ts %>% filter(feature == feature.names[i])

    iadonis <- suppressMessages(
      adonis(ix.ts[, -c(1,2)] ~ group, data = ix.ts, method = dist_method,
             ...)
    )
    ires <- iadonis$aov.tab[1, ]
    rownames(ires) <- feature.names[i]
    return(ires)
  })

  adonis.res <- do.call("rbind", adonis.res) %>%
    data.frame(check.names = FALSE) %>%
    rownames_to_column("feature") %>%
    rename("pval" = "Pr(>F)") %>%
    arrange(-R2) %>%
    mutate(p.adj = p.adjust(pval),
           method = p_adj_method)
  return(adonis.res)
}
