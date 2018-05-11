#' @title Cluster assignment using static brnach cutting.
#'
#' @description This function computes the cluster assignment for hierarchical
#' clustering results using static branch cutting.
#'
#' @param hclst an object of class \code{hclust} representing the clustering of
#' features of X.
#' @param h a the fraction of the max tree height at which to cut and assign
#' labels.
#' @param k an integer scalar or vector with the desired number of groups
#'
#' @return a data.frame containing cluster assignments.
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join arrange data_frame
#' @importFrom stats cutree
#'
assign_cluster_static <- function(hclst, h = NULL,  k = NULL) {
  if(!is.null(h)) {
    h <- h*max(hclst$height)
  }
  clst <- cutree(hclst, k = k, h = h)
  cluster <- suppressMessages(
    data.frame("cluster" = clst) %>%
      rownames_to_column(var = "feature") %>%
      left_join(
        data_frame("feature" = hclst$labels, "leaf_order" = hclst$order)
      ) %>% #arrange(leaf_ix) %>%
      mutate(cluster =  paste0("C", cluster))
  )
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
assign_cluster_dynamic <- function(hclst, max_height = 0.9, ...){
  max_height <- max_height*max(hclst$height)
  clst <- cutreeDynamicTree(hclst, maxTreeHeight = max_height, ...)
  cluster <- suppressMessages(
    data.frame("cluster" = clst, row.names = hclst$labels) %>%
      rownames_to_column(var = "feature") %>%
      left_join(
        data_frame("feature" = hclst$labels, "leaf_order" = hclst$order)
      ) %>% #arrange(leaf_ix) %>%
      mutate(cluster =  paste0("C", cluster))
  )
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
#' @param static_cut_params parameters for \code{\link{assign_cluster_static}}.
#' @param dynamic_cut_params parameters for
#' \code{\link{assign_cluster_dynamic}}.
#'
#' @return a list with the \code{hclust} object, as well as \code{clust_map}
#' and \code{clust_centroids} data.frames.
#'
#' @importFrom dplyr select left_join group_by summarise_all
#' @importFrom tibble rownames_to_column
#' @importFrom stats hclust
#' @export
#' @examples
#' X <- get_data(endoderm_small)
#' clust_res <- cluster_data(X)
#' head(clust_res$clust_centroids)
#' head(clust_res$clust_map)
#'
cluster_data <- function(
  X, dist = "euclidean",
  dynamic = FALSE,
  hclust_params = list(),
  static_cut_params = list(h = 0.5),
  dynamic_cut_params = list(max_height = 0.9)) {

  cluster <- feature <- NULL

  # Perform hierarchical clustering
  D <- dist(X, method = dist)
  hclust_params[["d"]] <- D
  hclst <- do.call(hclust, hclust_params)

  # Find clusters
  static_cut_params[["hclst"]] <- dynamic_cut_params[["hclst"]] <- hclst
  if(!dynamic) {
    clst.mapping <- do.call(assign_cluster_static, static_cut_params)
  } else {
    clst.mapping <- do.call(assign_cluster_dynamic, dynamic_cut_params)
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
  return(list(hclust = hclst, clust_map = clst.mapping,
              clust_centroids = clst.centroids))
}


#' @title Cluster timecourse features.
#'
#' @description Find the cluster assignment for timecourse features.
#' Clustering computed on top "n_top_feat" features most variable over time
#' in each of the selected "groups". The cluster assignment of the remaining
#' genes is based on the distance to the closest cluster centroid previously
#' obtained. Hierarchical clustering is performed and both static and dynamic
#' branch cutting algorithm are available for assigning cluster membership.
#'
#' @param object A \code{vistimeseq} object
#' @param n_top_feat A number of top most variable time-course features to use
#' for clustering.
#' @param groups One or multiple groups from \code{object@group} to take
#' into account when aggregating time-course features.
#' @param lambda Weights for each lag difference, for time-course data.
#' Length of \code{lambda} specifies number of lags to include. By default
#' lag of order one and two are included with coefficients 0.5 and 0.25
#' respectively.
#' @param clust_params A list contating arguments for hierarchical clustering.
#' For details see \code{\link{cluster_data}}.
#'
#' @return a \code{vistimeseq} object with cluster assignment stored
#' in \code{cluster.map} slot.
#'
#' @importFrom stats sd
#' @importFrom dplyr select left_join group_by summarise_all contains top_n
#' @importFrom tibble column_to_rownames
#' @importFrom proxy dist
#' @importFrom  methods slot<-
#' @importFrom methods validObject
#' @export
#' @examples
#' endoderm_small
#' endoderm_small <- cluster_timecourse_features(endoderm_small)
#' head(get_cluster_map(endoderm_small))
#'
cluster_timecourse_features <- function(
  object, n_top_feat = 1000, groups = "all", lambda = c(0.5, 0.25),
  clust_params = list()){

  group <- cluster <- feature <- NULL
  clust_params_default <- list(
    dist = "euclidean",
    dynamic = FALSE,
    hclust_params = list(),
    static_cut_params = list(h = 0.5),
    dynamic_cut_params = list()
  )
  clust_params_update  <- clust_params_default
  for (ele in names(clust_params)) {
    if (ele %in% names(clust_params_update)){
      clust_params_update[[ele]] <- clust_params[[ele]]
    }
  }

  if (!validObject(object)){
    stop("Invalid 'vistimeseq' object.")
  }
  if (is.null(time_course(object, collapsed = TRUE))) {
    message("Aggregating across replicates.")
    object <- collapse_replicates(object)
    message("Converting to timecorse format.")
    object <- convert_to_timecourse(object)
  }
  if(length(grep("Lags_", colnames(time_course(object)))) ){
    message("Adding lags with coefficients: ",
            paste0(lambda, collapse = " "))
    object <- add_lags(object, lambda = lambda)
  }

  tc_data <- time_course(object, collapsed = TRUE)
  # Find top "n_top_feat" most variable features
  n_top_feat <- min(n_top_feat, length(feature_names(object)))
  timepoints <- tc_data %>%
    select(-(feature:replicate), -contains("Lag")) %>%
    colnames()
  features_sd <- tc_data %>%
    select(-replicate, -contains("Lag"))
  features_sd[["sd"]] <- apply(features_sd[, timepoints], 1, sd)
  features_sd <- features_sd %>%
    group_by(group) %>%
    top_n(n = n_top_feat, wt = sd)
  top_features <- unique(features_sd[["feature"]])

  if (any(groups == "all")){
    groups <- unique(get_group(object))
  }
  tc_data <- tc_data %>%
    select(-replicate) %>%
    filter(group %in% groups)

  # Aggregate timecourses across all "groups" and recompute lags
  if(length(groups) > 1) {
    message("Averaging timecourses over all \"groups\" selected ",
            "and recomputing lags with coefficients: ",
            paste0(lambda, collapse = " "))
    tc_data <- tc_data %>%
      select(-contains("Lag_"), -group) %>%
      group_by(feature) %>%
      summarise_all(mean)
    tc_lags <- add_lags_to_tc(tc_data %>% select(-feature), lambda = lambda)
    tc_data <- cbind(feature = tc_data$feature, tc_lags)
  } else {
    tc_data <- select(tc_data, -group)
  }

  # Cluster a subset of features
  tc_subset <- tc_data %>%
    filter(feature %in% top_features) %>%
    column_to_rownames("feature")
  clust_params_update[["X"]] <- tc_subset
  res_cluster_subset <- do.call(cluster_data, clust_params_update)
  cluster_hclust <- res_cluster_subset$hclust
  clust_map <- res_cluster_subset$clust_map %>%
    select(feature, cluster)
  clust_centroids <-  res_cluster_subset$clust_centroids %>%
    column_to_rownames("cluster")

  # Assign the rest of the genes to the closest cluster centroid
  tc_remain <- tc_data %>%
    filter(!feature %in% clust_map$feature) %>%
    column_to_rownames("feature")
  dist_to_nearest_clust <- proxy::dist(tc_remain, clust_centroids)
  clst.remain <- apply(dist_to_nearest_clust, 1, function(x) {
    colnames(dist_to_nearest_clust)[which.min(x)]
  })
  clust_map_remain <- data.frame(
    feature = names(clst.remain),
    cluster = clst.remain,
    stringsAsFactors = FALSE)

  clust_res <- list(hclust = cluster_hclust,
    cluster_map = rbind(clust_map, clust_map_remain))
  slot(object, name = "cluster.features", check = TRUE) <- clust_res
  return(object)
}


