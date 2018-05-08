#' @title Plot a heatmap
#'
#' @description Generates a heatmap of expression data.
#'
#' @param object A \code{vistimeseq} object
#' @param num.feat Number of top most features to use.
#' @param scale Whether to scale the data (by features) before plotting.
#' @param feat_desc One of the column names from \code{object@feature.data}
#' to describe the features.
#' @param sample_desc ne of the column names from \code{object@sample.data}
#' to describe the samples.
#' @param ... other parameters ggplot.
#'
#' @return Returns a \code{ggplot2} objet.
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom viridis viridis
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar
#' @export
#' @examples
#' endoderm_small
#' endoderm_small <- normalize_data(endoderm_small)
#' plot_heatmap(endoderm_small)
#'
plot_heatmap <- function(object, num.feat = 200, scale = TRUE,
                         feat_desc = "feature", sample_desc = "sample"){
  if (!validObject(object)) {
    stop("Invalid vistimeseq object.")
  }
  if (!feat_desc %in% colnames(object@feature.data)){
    stop("\"feat_desc\" not in object@feature.data")
  }
  if (!sample_desc %in% colnames(object@sample.data)){
    stop("\"sample_desc\" not in object@sample.data")
  }
  if(is.null(object@data)){
    cnts <- object@data
  } else {
    cnts <- object@raw.data
  }
  rownames(cnts) <- object@feature.data[, feat_desc]
  colnames(cnts) <- object@sample.data[, sample_desc]
  top_feat <- apply(cnts, 1, sd)
  top_feat <- names(top_feat[order(-top_feat)])[1:num.feat]
  cols_ordered <- order(object@group, object@replicate, object@time)
  Y <- cnts[top_feat, cols_ordered]
  if (scale) {
    Y <- t(scale(t(Y), center = TRUE, scale = TRUE))
  }

  smpdf <- object@sample.data %>%
    select(group, replicate, time) %>%
    arrange(group, replicate, time) %>%
    mutate(time = as.numeric(time))

  n_group <- length(unique(object@group))
  group_cols <- colorRampPalette(
    colors = brewer.pal(9, name = "Set1")[1:min(9, n_group)])(n_group)
  names(group_cols) <- unique(object@group)

  n_replicates <- length(unique(object@replicate))
  rep_cols <- colorRampPalette(
    colors = brewer.pal(8, name = "Set3")[1:min(8, n_replicates)])(n_replicates)
  names(rep_cols) <- unique(object@replicate)

  time_cols <- colorRamp2(
    breaks = seq(min(object@time), max(object@time), length.out = 10),
    colors =  viridis(10))

  ha1 <- ComplexHeatmap::HeatmapAnnotation(smpdf,
    col = list(group = group_cols, replicate = rep_cols, time = time_cols))

  ComplexHeatmap::Heatmap(
    Y, name = "Z-score", cluster_columns = FALSE,
    top_annotation = ha1, row_names_gp = gpar(fontsize = 8))
}


#' @title Plot a standard PCA
#'
#' @description Generates a standard PCA plot of observations in the dataset.
#'
#' @param object A \code{vistimeseq} object
#' @param axis An integer vector indicating principal components to use for
#' plotting, by default 1:2.
#' @param col.var A character string indicating a column
#' from object@@sample.data which should be used for coloring
#' the points. By default NULL.
#' @param ... other parameters ggplot.
#'
#' @return Returns a \code{ggplot2} objet.
#'
#' @importFrom viridis viridis
#' @importFrom dplyr left_join
#' @importFrom tibble rownames_to_column column_to_rownames
#' @export
#' @examples
#' endoderm_small
#' endoderm_small <- run_pca(endoderm_small)
#' plot_sample_pca(endoderm_small, col.var = "group")
#'
plot_sample_pca <- function(object, axis = 1:2, col.var = NULL, ...) {
  if(! "pca_sample" %in% names(object@dim.red)) {
    stop("No 'pca_sample' in object@dim.red. Run PCA for samples first.")
  }
  axis <- if(is.numeric(axis)) paste0("PC", axis) else axis
  pca.sample <- object@dim.red$pca_sample
  pca.eigs <- object@dim.red$pca_eigs

  pca.scores <- pca.sample[, axis] %>%
    as.data.frame() %>%
    rownames_to_column("sample")

  if(all(pca.scores$sample %in% object@sample.data$sample)) {
    pca.scores <- suppressMessages(
      pca.scores %>%
        left_join(object@sample.data) %>%
        column_to_rownames("sample")
    )
  } else if (all(pca.scores$sample %in% object@sample.data.collapsed$sample)) {
    pca.scores <- suppressMessages(
      pca.scores %>%
        left_join(object@sample.data.collapsed) %>%
        column_to_rownames("sample")
    )
  } else {
    stop("Sample names in sample data and PCA coordinates disagree.")
  }

  axis_label <- paste0(
    colnames(pca.scores)[1:2], " [",
    signif(pca.eigs[1:2]/sum(pca.eigs)*100, 3), "%]"
  )
  plt <- ggplot(
    data = pca.scores,
    aes(x = pca.scores[[1]], y = pca.scores[[2]])
    ) +
    geom_point(aes_string(fill = col.var), color = "grey80", pch = 21, ...) +
    geom_hline(aes(yintercept =0), size=.2) +
    geom_vline(aes(xintercept = 0), size=.2) +
    xlab(axis_label[1]) +
    ylab(axis_label[2]) +
    coord_fixed(sqrt(pca.eigs[2]/pca.eigs[1]))

  if(all(!is.null(pca.scores), is.numeric(pca.scores[, col.var]))){
    plt <- plt + viridis::scale_fill_viridis()
  }
  return(plt)
}


#' @title Overlay (time) series over PCA grid
#'
#' @param object A \code{vistimeseq} object.
#' @param axis An integer vector indicating principal components to use for
#' plotting, by default 1:2.
#' @param m a number of tiles in a grid in the horizontal direction.
#' @param n a number of tiles in a grid in the vertical direction.
#' @param group.highlight An optional character string indicating the group
#' subset for which the time-course trends should be plotted. By default all
#' time-course trends are plotted for all groups.
#' @param linecol a vector indicating the color of the gene profile trend line,
#' different for each group.
#' @param ... other parameters for the line plots.
#'
#' @importFrom Hmisc subplot
#' @importFrom proxy dist
#' @importFrom viridis viridis
#' @export
#' @examples
#' endoderm_small
#' endoderm_small <- run_pca(endoderm_small)
#' plot_ts_pca(endoderm_small)
#'
plot_ts_pca <- function(object, axis = 1:2, m = 20, n = 20,
                        group.highlight = NULL, linecol = NULL,
                        ...) {
  if (!validObject(object)){
    stop("Invalid vistimeseq object.")
  }
  if(! "pca_feature" %in% names(object@dim.red)) {
    stop("No 'pca_feature' in object@dim.red. Run PCA for features first.")
  }
  if (! "tc_collapsed" %in% names(object@timecourse.data)) {
    message("First, aggregating replicates and converting",
            "to time-course format.")
    object <- collapse_replicates(object, FUN = mean)
    object <- convert_to_timecourse(object)
  }
  # Prepare scores data
  tc <- object@timecourse.data$tc_collapsed
  tc <- tc[, !grepl("Lag", colnames(tc))]
  tmp <- as.numeric(colnames(tc %>% select(-feature, -group, -replicate)))
  axis <- if(is.numeric(axis)) paste0("PC", axis) else axis
  pca.feature <- object@dim.red$pca_feature
  pca.eigs <- object@dim.red$pca_eigs
  pca.loadings <- suppressMessages(
    pca.feature[, axis] %>%
    as.data.frame() %>%
    rownames_to_column("feature") %>%
    left_join(object@feature.data) %>%
    column_to_rownames("feature")
  )
  colnames(pca.loadings)[1:2] <-
    paste0(colnames(pca.loadings)[1:2], " [",
           signif(pca.eigs[1:2]/sum(pca.eigs)*100, 3), "%]")

  # Create a grid over score values
  mins <- apply(pca.loadings[, 1:2], 2, min)
  maxes <- apply(pca.loadings[, 1:2], 2, max)
  x <- seq(mins[1], maxes[1], length.out = m)
  y <- seq(mins[2], maxes[2], length.out = n)
  dx <- x[2] - x[1]; dy <- y[2] - y[1]
  grid <- expand.grid(x, y)

  # Find gene closest to the grid center
  xD <- proxy::dist(grid[, 1], pca.loadings[, 1])
  yD <- proxy::dist(grid[, 2], pca.loadings[, 2])
  D <- proxy::dist(grid, pca.loadings[, 1:2])

  min_dists <- apply(D, 1, min)
  min_dists_ix <- apply(D, 1, which.min)
  x_min_dists <- sapply(1:nrow(xD), function(i) xD[i, min_dists_ix[i]])
  y_min_dists <- sapply(1:nrow(yD), function(i) yD[i, min_dists_ix[i]])
  min_dists_ix[x_min_dists > dx/2 | y_min_dists > dy/2] <- NA

  # Plot all points corresponding to each feature
  par(mar=par()$mar * c(2, 1.2, 1.2, 1.2), xpd=TRUE, #c(par()$mar[1], 0, 0, 0)
      cex = 0.7, cex.main = 2, cex.axis = 1.5, cex.lab = 1.5)
  plot(pca.loadings[, 1:2], type = "p", pch = 16,
       xlim = c(mins[1] - dx/2, maxes[1] + dx/2),
       ylim = c(mins[2] - dy/2, maxes[2] + dy/2),
       asp=max(sqrt(pca.eigs[2]/pca.eigs[1]), 0.5), ...)

  groups.unique <-
    if (is.null(group.highlight)) unique(object@group) else group.highlight

  if(is.null(linecol)) {
    linecol <- viridis(length(groups.unique))
    names(linecol) <- groups.unique
  }
  ylimits <- c(min(tc %>% select(-feature, -group, -replicate)),
               max(tc %>% select(-feature, -group, -replicate)))

  # Plot all the time-course profiles
  for(i in seq_along(min_dists_ix)) {
    igene <- rownames(pca.feature)[min_dists_ix[i]]
    if(is.na(min_dists_ix[i])) next
    for (gr in groups.unique) {
      gTC <- tc %>%
        filter(feature == igene, group == gr) %>%
        select(-feature, -group, -replicate) %>%
        as.numeric()
      Hmisc::subplot(
        plot(tmp, gTC, type = "l", lwd = 2,
             col =  linecol[gr], frame = F, axes = F,
             xlab = "", ylab = "", ylim = ylimits),
             x = c(grid[i, 1] - dx/2, grid[i, 1] + dx/2),
             y = c(grid[i, 2] - dy/2, grid[i, 2] + dy/2)
      )
    }
  }
  # Add a legend
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
      mar = c(1, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", inset = c(0, 0), horiz = TRUE,
         groups.unique, col = linecol, xpd = TRUE,
         lty = c(1,1), lwd = c(3,3))
}

#' @title Plot (time) series over clusters.
#'
#' @description Plots timecourse (aggregated over replicates) feature data
#' faceted by computed clusters and experimental group.
#'
#' @param object A \code{vistimeseq} object.
#' @param alpha transparency of trajectory lines.
#' @param ncol number of columns in the facte plot.
#'
#' @return ggplot object
#' @importFrom dplyr left_join mutate select group_by summarize_all arrange
#' @importFrom tidyr gather
#'
#' @export
#' @examples
#' endoderm_small <- cluster_timecourse_features(endoderm_small)
#' plot_ts_clusters(endoderm_small)
#'
plot_ts_clusters <- function(object, alpha = 0.5, ncol = 4) {
  if (!validObject(object)){
    stop("Invalid vistimeseq object.")
  }
  if(! "cluster_map" %in% names(object@cluster.features)) {
    stop("No 'cluster_map' in object@cluster.features. Perform
         clustering with 'cluster_timecourse_features()' first.")
  }
  if (! "tc_collapsed" %in% names(object@timecourse.data)) {
    stop("No 'tc_collapsed' in object@timecourse.data, but 'cluster.features'
         slot is non-emty. Check if data 'vistimeseq' object is corrupted.")
  }
  cluster_map <- object@cluster.features$cluster_map
  freq_df <- cluster_map %>%
    group_by(cluster) %>%
    summarise(freq = n()) %>%
    arrange(desc(freq))

  cat_levels <- paste0("[", freq_df$cluster, ": ", freq_df$freq, "]")
  cat_levels <- paste(rep(cat_levels, each = length(unique(object@group))),
                      rep(unique(object@group), length(cat_levels)))
  tc_data <- suppressMessages(
    object@timecourse.data$tc_collapsed %>%
      select(-replicate, -contains("Lag_")) %>%
      left_join(cluster_map) %>%
      gather(
        key = "time", value  = "value",
        -feature, -group, -cluster
      ) %>%
      left_join(freq_df) %>%
      mutate(
        time = as.numeric(time),
        value = as.numeric(value),
        category = paste0("[", cluster, ": ", freq, "] ", group)
      ) %>%
      arrange(cluster, group, feature, time) %>%
      mutate(
        category = factor(category, levels = cat_levels))
  )
  # Compute cluster mean expression profile for each group
  tc_cluster_mean <- suppressMessages(
    tc_data %>%
      select(-feature) %>%
      group_by(cluster, group, category, time) %>%
      summarize_all(mean) %>%
      left_join(freq_df) %>%
      arrange(cluster, group)
  )

  plt <- ggplot(
    tc_data,
    aes(y = value , x = time, color = group)
    ) +
    geom_line(
      aes(group = feature), alpha = alpha
    ) +
    geom_point() +
    geom_line(
      data = tc_cluster_mean, lwd = 1.5, color = "grey50",
      aes(group = group)
    ) +
    facet_wrap(~category, scales = "free", ncol = ncol)
  return(plt)
}



#' @title Plot selected time series.
#'
#' @param X a data matrix or data.frame
#' @param time a vector of length equal to ncol(X) indicating the time
#' variable corresponding to the sample (data column).
#' @param keep_rows slected features to plot.
#' @param group a vector of length equal to ncol(X) indicating
#' the group membership of the sample.
#' @param replicate a vector of length equal to ncol(X) indicating
#' the replicate (e.g. individual).
#' @param feature_data a data.frame contating the information on features.
#' @param ncol an integer indicating the number of columns for facetting
#'
#' @return list of ggplot objects
#'
#' @importFrom dplyr filter group_by mutate left_join select summarise_all
#' @importFrom tibble column_to_rownames
#' @export
plot_time_series <- function(X, time, keep_rows = rownames(X),
                             group = NULL, replicate = NULL,
                             feature_data = NULL, ncol = 5,
                             smooth = TRUE){
  if(is.null(group))
    group <- rep("G1", ncol(X))
  if(is.null(replicate))
    replicate <- rep("R1", ncol(X))
  if(ncol(X) != length(group))
    stop("ncol(X) must match length(group)")
  if(ncol(X) != length(replicate))
    stop("ncol(X) must match length(replicate)")
  if(ncol(X) != length(time))
    stop("ncol(X) must match length(time)")
  if(!is.null(keep_rows) & !all(keep_rows %in% rownames(X)))
    stop("keep_rows must be a subset of rownames(X)")

  sample.data <- data.frame(sample = colnames(X),
                            group, replicate, time,
                            stringsAsFactors = FALSE)
  feature.data <- data.frame(row.names = rownames(X),
                             feature = rownames(X))

  if (!is.null(feature_data))
    feature.data <- data.frame(feature.data, feature_data,
                               stringsAsFactors = FALSE)
  if(!"symbol" %in% colnames(feature.data))
    feature.data$symbol <- rownames(feature.data)

  feature.data <- feature.data[keep_rows, ]

  X.long <- suppressMessages(
    melt_matrix(t(X[keep_rows, ])) %>%
      left_join(sample.data) %>%
      left_join(feature.data) %>%
      mutate(
        symbol = factor(symbol, levels = feature.data$symbol)
      )
  )

  plt <- ggplot(X.long, aes(x = time, y = value, color = group)) +
      geom_point(size = 1) +
      facet_wrap(~ symbol, scales = "free", ncol = ncol)
  if(length(unique(replicate)) > 1) {
    plt <- plt + geom_line(aes(group = replicate), lty = 3, alpha = 0.7)
  }
  if(smooth) {
    plt <- plt + geom_smooth(aes(x = as.numeric(time)), lwd = 1.5)
  } else {
    X.long.mean <- X.long %>%
      mutate(time = as.numeric(time)) %>%
      group_by(symbol, group, time) %>%
      select(symbol, time, group, value) %>%
      summarise_all(mean)
    plt <- plt +
      geom_line(data = X.long.mean, lwd = 1.5,
                aes(x = time, y = value, color = group))
  }
  return(plt)
}


#' @title Plot Gene Ontology limma::goana enrichment results.
#'
#' @param DF a data matrix or data.frame results from limma::goana
#' with gene list. Must contain columns Term, DE, and P.DE.
#'
#' @return a ggplot object
#'
#' @importFrom dplyr filter group_by mutate left_join select summarise_all
#' @importFrom tibble column_to_rownames
#' @export
plot_goana <- function(DF, nGSmax = 15) {
  DF <- DF %>%  arrange(-DE, P.DE) %>%
    mutate(Term = paste0(Term, " (", DE, "/", N , ")"),
           Term = factor(Term, level = Term))

  plt <- ggplot(DF[1:min(nGSmax, nrow(DF)), ],
                aes(y = Term, x = -log10(P.DE),
                    size = N, color = DE/N)) +
    geom_point() +  viridis::scale_color_viridis()
  return(plt)
}


#' #' @title Venn diagram for features significant at each timpoints
#' #'
#' #' @param df a data.frame of features by timepoints listing significant
#' #' features for each timepoint, NA entries stands for non-significant,
#' #' features do not have to be in the same order.
#' #' @param size plot size, in centimeters.
#' #' @param transparency transparency for the color(s) specified with zcolor
#' #' @param cexil a character expansion for the intersection labels
#' #' @param cexsn a character expansion for the set names
#' #' @param zcolor a vector of colors for the custom zones,
#' #' or predefined colors if "style"
#' #' @param ... other parameters for the line plots.
#' #'
#' #' @importFrom venn venn
#' #' @importFrom viridis viridis
#' #' @export
#' #' @examples
#' #'
#' timepoint_venn <- function(df, size = 50, transparency = 0.5,
#'                            cexil = .9, cexsn = 1,
#'                            zcolor = viridis(length(venn.lst)), ...) {
#'   venn.lst <- lapply(colnames(df), function(i) {
#'     df[!is.na(df[, i]), i]
#'   })
#'   names(venn.lst) <- colnames(df)
#'   venn(venn.lst, size = size, transparency = transparency,
#'        cexil = cexil, cexsn = cexsn, zcolor = viridis(length(venn.lst)), ...)
#' }


