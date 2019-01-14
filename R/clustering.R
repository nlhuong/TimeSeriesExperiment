#' @title Cluster assignment using static branch cutting.
#'
#' @description This function computes the cluster assignment for hierarchical
#' clustering results using static branch cutting.
#'
#' @param hclst an object of class \code{hclust} representing the clustering 
#' of features of X.
#' @param h a the fraction of the max tree height at which to cut and assign
#' labels.
#' @param k an integer scalar or vector with the desired number of groups
#'
#' @return a data.frame containing cluster assignments.
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join arrange data_frame
#' @importFrom stats cutree
#' @export
assignClusterStatic <- function(hclst, h = NULL, k = NULL) {
    if(!is.null(h)) {
        h <- h*max(hclst$height)
    }
    clst <- cutree(hclst, k = k, h = h)
    cluster <- suppressMessages(
        data.frame("cluster" = clst) %>%
            rownames_to_column(var = "feature") %>%
            left_join(data_frame("feature" = hclst$labels, 
                                 "leaf_order" = hclst$order)))
            # arrange(leaf_ix) 
    return(cluster)
}


#' @title Cluster assignment using dynamic branch cutting.
#'
#' @description This function computes the cluster assignment for hierarchical
#' clustering results using dynamic branch cutting.
#'
#' @param hclst an object of class \code{hclust} representing the clustering of
#' features of X.
#' @param max_height a fraction of the total tree height used to compute
#' \code{maxTreeHeight} argument for \code{\link[dynamicTreeCut]{cutreeDynamic}}
#' function from \code{dynamicTreeCut} package.
#' @param ... other parameters for \code{\link[dynamicTreeCut]{cutreeDynamic}}.
#'
#' @return a data.frame containing cluster assignments.
#' @importFrom dynamicTreeCut cutreeDynamicTree
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join arrange data_frame
#'
assignClusterDynamic <- function(hclst, max_height = 0.9, ...){
    max_height <- max_height*max(hclst$height)
    clst <- cutreeDynamicTree(hclst, maxTreeHeight = max_height, ...)
    cluster <- suppressMessages(
        data.frame("cluster" = clst, row.names = hclst$labels) %>%
            rownames_to_column(var = "feature") %>%
            left_join(
                data_frame("feature" = hclst$labels, "leaf_order" = hclst$order)
            ))
    return(cluster)
}


#' @title Cluster time series data.
#'
#' @description Perform cluster assignment using hierarchical clustering
#' and static or dynamic branch cutting.
#' If 'subset' is specified the clsutering is performed on the subset
#' of the data, and the rest of the rows are assigned based on the distances
#' to the centroids of computed clusters.
#'
#' @param X a data matrix or data.frame where rows are time series.
#' @param dist the distance metric for the dissimilarity used for clustering.
#' @param dynamic whether dynamic branch cutting should be done for cluster
#' assignment.
#' @param hclust_params parameters for \code{\link[stats]{hclust}} function.
#' @param static_cut_params parameters for \code{\link{assignClusterStatic}}.
#' @param dynamic_cut_params parameters for
#' \code{\link{assignClusterDynamic}}.
#'
#' @return a list with the \code{hclust} object, as well as \code{clust_map}
#' and \code{clust_centroids} data.frames.
#'
#' @importFrom dplyr select left_join group_by summarise_all
#' @importFrom tibble rownames_to_column
#' @importFrom stats hclust
#' @export
#' @examples
#' data("endoderm_small")
#' endoderm_small <- normalizeData(endoderm_small)
#' X <- SummarizedExperiment::assays(endoderm_small)$norm
#' clust_res <- clusterData(X)
#' head(clust_res$clust_centroids)
#' head(clust_res$clust_map)
#'
clusterData <- function(X, dist = "euclidean", dynamic = FALSE,
                        hclust_params = list(), 
                        static_cut_params = list(h = 0.5),
                        dynamic_cut_params = list(max_height = 0.9)) 
{
  
    cluster <- feature <- NULL
    X <- as.data.frame(X)
    dynparams <- paste(names(dynamic_cut_params), "=", dynamic_cut_params,
                       collapse = ", ")
    statparams <- paste(names(static_cut_params), "=", static_cut_params,
                       collapse = ", ")
    settings <- ifelse(
      dynamic, 
      paste0("Dynamic cluster assignment with params: ", dynparams),
      paste0("Static cluster assignment with params: ", statparams)
    )
    
    # Perform hierarchical clustering
    D <- dist(X, method = dist)
    hclust_params[["d"]] <- D
    hclst <- do.call(hclust, hclust_params)

    # Find clusters
    static_cut_params[["hclst"]] <- dynamic_cut_params[["hclst"]] <- hclst
    if(!dynamic) {
        clst.mapping <- do.call(assignClusterStatic, static_cut_params)
    } else {
        clst.mapping <- do.call(assignClusterDynamic, dynamic_cut_params)
        # Label 0 in dynamicCutTree means the object was unassigned,
        # so we remove these features.
        clst.mapping <- clst.mapping[clst.mapping$cluster != "C0", ]
    }
    
    # Compute cluster centroids
    clst.centroids <- suppressMessages(
        clst.mapping %>%
            select(feature, cluster) %>%
            left_join(X %>% rownames_to_column("feature")) %>%
            select(-feature) %>%
            group_by(cluster) %>%
            summarise_all(mean) %>%
            data.frame(check.names = FALSE)
    )
    
    return(list(settings = settings, hclust = hclst, 
                clust_map = clst.mapping, clust_centroids = clst.centroids))
}


#' @title Cluster time series features.
#'
#' @description Find the cluster assignment for timecourse features.
#' Clustering computed for top "n.top.feat" features most variable over time
#' in each of the selected "groups" using time-series expression 
#' (collpased over replicates). The cluster assignment of the remaining
#' genes is based on the distance to the closest cluster centroid previously
#' obtained. Hierarchical clustering is performed and both static and dynamic
#' branch cutting algorithm are available for assigning cluster membership.
#'
#' @param object A \code{TimeSeriesExperiment} object
#' @param n.top.feat A number of top most variable time-course features to use
#' for clustering.
#' @param groups.selected One or multiple groups from \code{object@group} 
#' to take into account when aggregating time-course features.
#' @param lambda Weights for each lag difference, for time-course data.
#' Length of \code{lambda} specifies number of lags to include. By default
#' lag of order one and two are included with coefficients 0.5 and 0.25
#' respectively.
#' @param clust.params A list contating arguments for hierarchical clustering.
#' For details see \code{\link{clusterData}}.
#'
#' @return a \code{TimeSeriesExperiment} object with cluster assignment stored
#' in \code{cluster.map} slot.
#'
#' @importFrom stats sd
#' @importFrom dplyr select 
#' @importFrom dplyr left_join group_by summarise_all starts_with top_n
#' @importFrom tibble column_to_rownames
#' @importFrom proxy dist
#' @importFrom methods slot<- validObject
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment rowData
#' @export
#' @examples
#' data("endoderm_small")
#' endoderm_small <- clusterTimeSeries(endoderm_small)
#' head(clusterMap(endoderm_small))
#'
clusterTimeSeries <- function(object, n.top.feat = 1000, 
                              groups.selected = "all", lambda = c(0.5, 0.25), 
                              clust.params = list())
{
    if (!is(object, "TimeSeriesExperiment")) 
        stop("Input must be a 'TimeSeriesExperiment' object.")
    if (!validObject(object)) 
        stop("Invalid TimeSeriesExperiment object.")
    
    clust_params_default <- list(
        dist = "euclidean",
        dynamic = FALSE,
        hclust_params = list(),
        static_cut_params = list(h = 0.5),
        dynamic_cut_params = list()
    )
    clust_params_update  <- clust_params_default
    for (ele in names(clust.params)) {
        if (ele %in% names(clust_params_update)){
            clust_params_update[[ele]] <- clust.params[[ele]]
        }
    }
    
    group <- cluster <- feature <- NULL
    n.top.feat <- min(n.top.feat, length(rownames(object)))
    if (any(groups.selected == "all")){
        groups.selected <- unique(groups(object))
    }
    
    if (is.null(timeSeries(object, name = "ts_collapsed"))) {
        object <- collapseReplicates(object)
        object <- makeTimeSeries(object)
    }
    
    ts_collapsed <- timeSeries(object, "ts_collapsed") %>%
        select(-replicate) %>%
        filter(group %in% groups.selected)  # filter to only the chosen groups
    # Find top "n.top.feat" most variable features
    feat_tmps <- ts_collapsed %>%
        select(-starts_with("Lag"))
    feat_tmps$sd <- apply(feat_tmps %>% select(-feature, -group),
                          1, sd, na.rm = TRUE)

    feat_sd <- feat_tmps %>%
        group_by(group) %>%
        top_n(n = n.top.feat, wt = sd) 
    top_features <- unique(feat_sd[["feature"]])

    # Aggregate timecourses across all "groups" and recompute lags
    if(length(groups.selected) > 1) {
        message("Averaging timecourses over all 'groups' selected ",
                "and recomputing lags with coefficients: ",
                 paste0(lambda, collapse = " "))
        ts_collapsed <- feat_tmps %>%
            select(-group, -sd) %>%
            group_by(feature) %>%
            summarise_all(mean)
        ts_with_lags <- .addLagsToTimeSeries(
            ts_collapsed %>% select(-feature), lambda = lambda)
        ts_collapsed <- cbind(feature = ts_collapsed$feature, ts_with_lags)
    } else {
        if(!any(grepl("Lag_", colnames(timeSeries(object, "ts_collapsed"))))){
            object <- addLags(object, lambda = lambda)
        }
        ts_collapsed <- select(feat_tmps, -group)
    }
    # Cluster a subset of features
    ts_subset <- ts_collapsed %>%
        filter(feature %in% top_features) %>%
        column_to_rownames("feature")
    clust_params_update[["X"]] <- ts_subset
    res_cluster_subset <- do.call(clusterData, clust_params_update)
    cluster_hclust <- res_cluster_subset$hclust
    clust_map <- res_cluster_subset$clust_map %>%
        select(feature, cluster)
    clust_centroids <-  res_cluster_subset$clust_centroids %>%
        column_to_rownames("cluster")

    # Assign the rest of the genes to the closest cluster centroid
    ts_remain <- ts_collapsed %>%
        filter(!feature %in% clust_map$feature) %>%
        column_to_rownames("feature")
    dist_to_nearest_clust <- proxy::dist(ts_remain, clust_centroids)
    clst.remain <- apply(dist_to_nearest_clust, 1, function(x) {
        colnames(dist_to_nearest_clust)[which.min(x)]
    })
    clust_map_remain <- data.frame(
        feature = names(clst.remain),
        cluster = clst.remain)
    clust_map$used_for_hclust <- rep(TRUE, nrow(clust_map))
    clust_map_remain$used_for_hclust <- rep(FALSE, nrow(clust_map_remain))
    final_cluster_map = rbind(clust_map, clust_map_remain)
   
    freq_df <- final_cluster_map %>%
      group_by(cluster) %>%
      summarise(freq = n()) %>%
      arrange(desc(freq)) %>%
      mutate(cluster_name = paste0("C", seq_len(nrow((.)))))

    final_cluster_map <- final_cluster_map %>%
      mutate(cluster = factor(
        cluster, levels = freq_df$cluster, labels = freq_df$cluster_name))
    
    rowData(object) <- suppressMessages(
      DataFrame(as.data.frame(rowData(object)) %>% 
                  left_join(final_cluster_map))
    )

    res_cluster_subset$clust_map <- res_cluster_subset$clust_map %>%
      mutate(cluster = factor(
        cluster, levels = freq_df$cluster, labels = freq_df$cluster_name))
    
    res_cluster_subset$final_cluster_map <- final_cluster_map
    slot(object, name = "clusterAssignment", check = TRUE) <- 
      res_cluster_subset
    return(object)
}


