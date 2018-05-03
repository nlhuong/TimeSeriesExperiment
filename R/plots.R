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
#' endoderm_small <- normalize_data(endoderm_small)
#' endoderm_small <- run_pca(endoderm_small, type = "sample")
#' plot_sample_pca(endoderm_small)
#'
plot_sample_pca <- function(object, axis = 1:2, col.var = NULL, ...) {
  if (!validObject(object))
    stop("Invalid vistimeseq object.")
  if(! "pca_sample" %in% names(object@dim.red)) {
    stop("No 'pca_sample' in object@dim.red. Run PCA for samples first.")
  }
  axis <- if(is.numeric(axis)) paste0("PC", axis) else axis
  pca.sample <- object@dim.red$pca_sample
  pca.eigs <- pca.sample$pca.eigs
  pca.scores <- suppressMessages(
    pca.sample$pca.scores[, axis] %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    left_join(object@sample.data) %>%
    column_to_rownames("sample")
  )
  axis_name <- paste0(
    colnames(pca.scores)[1:2], " [",
    signif(pca.eigs[1:2]/sum(pca.eigs)*100, 3), "%]"
  )

  plt <- ggplot(
    data = pca.scores,
    aes(x = pca.scores[[1]], y = pca.scores[[2]])
    ) +
    geom_point(aes_string(color = col.var), ...) +
    geom_hline(aes(yintercept =0), size=.2) +
    geom_vline(aes(xintercept = 0), size=.2) +
    xlab(axis_name[1]) +
    ylab(axis_name[2]) +
    coord_fixed(pca.eigs[2]/pca.eigs[1])

  if(is.numeric(pca.scores[[col.var]]))
    plt <- plt + viridis::scale_color_viridis()
  return(plt)
}


#' @title Overlay (time) series over PCA grid
#'
#' @param object A \code{vistimeseq} object.
#' @param axis An integer vector indicating principal components to use for
#' plotting, by default 1:2.
#' @param m a number of tiles in a grid in the horizontal direction.
#' @param n a number of tiles in a grid in the vertical direction.
#' @param group.selected An optional character string indicating the group
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
plot_ts_pca <- function(object, axis = 1:2, m = 20, n = 20,
                        group.selected = NULL, linecol = NULL,
                        ...) {
  if (!validObject(object))
    stop("Invalid vistimeseq object.")
  if(! "pca_feature" %in% names(object@dim.red)) {
    stop("No 'pca_feature' in object@dim.red. Run PCA for features first.")
  }
  if (! "tc_collapsed" %in% names(object@timecourse.data)) {
    stop("No 'tc_collapsed' entry in object@timecourse.data.",
         " Use 'collapse_replicates()' and 'convert_to_timecourse()'",
         " first to collapse the data over replicates and convert to",
         " time-course format.")
  }
  # Prepare scores data
  tc <- object@timecourse.data$tc_collapsed
  tc <- tc[, !grepl("Lag", colnames(tc))]
  tmp <- as.numeric(colnames(tc %>% select(-feature, -group, -replicate)))
  axis <- if(is.numeric(axis)) paste0("PC", axis) else axis
  pca.feature <- object@dim.red$pca_feature
  pca.eigs <- pca.feature$pca.eigs
  pca.scores <- suppressMessages(
    pca.feature$pca.scores[, axis] %>%
    as.data.frame() %>%
    rownames_to_column("feature") %>%
    left_join(object@feature.data) %>%
    column_to_rownames("feature")
  )
  colnames(pca.scores)[1:2] <-
    paste0(colnames(pca.scores)[1:2], " [",
           signif(pca.eigs[1:2]/sum(pca.eigs)*100, 3), "%]")

  # Create a grid over score values
  mins <- apply(pca.scores[, 1:2], 2, min)
  maxes <- apply(pca.scores[, 1:2], 2, max)
  x <- seq(mins[1], maxes[1], length.out = m)
  y <- seq(mins[2], maxes[2], length.out = n)
  dx <- x[2] - x[1]; dy <- y[2] - y[1]
  grid <- expand.grid(x, y)

  # Find gene closest to the grid center
  xD <- proxy::dist(grid[, 1], pca.scores[, 1])
  yD <- proxy::dist(grid[, 2], pca.scores[, 2])
  D <- proxy::dist(grid, pca.scores[, 1:2])

  min_dists <- apply(D, 1, min)
  min_dists_ix <- apply(D, 1, which.min)
  x_min_dists <- sapply(1:nrow(xD), function(i) xD[i, min_dists_ix[i]])
  y_min_dists <- sapply(1:nrow(yD), function(i) yD[i, min_dists_ix[i]])
  min_dists_ix[x_min_dists > dx/2 | y_min_dists > dy/2] <- NA

  # Plot all points corresponding to each feature
  par(mar=par()$mar + c(0.5*par()$mar[1], 0, 0, 0), xpd=TRUE)
  plot(pca.scores[, 1:2], type = "p", pch = 16,
       xlim = c(mins[1] - dx/2, maxes[1] + dx/2),
       ylim = c(mins[2] - dy/2, maxes[2] + dy/2), ...)

  groups.unique <-
    if (is.null(group.selected)) unique(object@group) else group.selected

  if(is.null(linecol)) {
    linecol <- viridis(length(groups.unique))
    names(linecol) <- groups.unique
  }
  ylimits <- c(min(tc %>% select(-feature, -group, -replicate)),
               max(tc %>% select(-feature, -group, -replicate)))

  # Plot all the time-course profiles
  for(i in seq_along(min_dists_ix)) {
    igene <- rownames(X)[min_dists_ix[i]]
    if(is.na(min_dists_ix[i])) next
    for (gr in groups.unique) {
      gTC <- tc %>%
        filter(feature == igene, group == gr) %>%
        select(-feature, -group, -replicate) %>%
        as.numeric()
      Hmisc::subplot(plot(tmp, gTC, type = "l", lwd = 2,
                          col =  linecol[gr], frame = F, axes = F,
                          xlab = "", ylab = "", ylim = ylimits),
                     x = c(grid[i, 1] - dx/2, grid[i, 1] + dx/2),
                     y = c(grid[i, 2] - dy/2, grid[i, 2] + dy/2))
    }
  }

  # Add a legend
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(1, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", inset = c(0, 0), horiz = TRUE,
         groups.unique, col = linecol, xpd = TRUE,
         lty = c(1,1), lwd = c(3,3))
}

#' @title Plot (time) series over clusters.
#'
#' @param X a data matrix or data.frame
#' @param cluster a data.frame containing 'feature' and 'cluster' columns.
#' @param time a vector of length equal to ncol(X) indicating the time
#' variable corresponding to the sample (data column).
#' @param group a vector of length equal to ncol(X) indicating
#' the group membership of the sample.
#' @param freq.df a data.frame contating the 'cluster' and 'freq' columns
#' which indicates the count of genes in each cluster.
#'
#' @return list of ggplot objects
#' @importFrom dplyr filter group_by left_join mutate select summarize_all
#' @export
plot_ts_clusters <- function(X, cluster, time, group = NULL,
                             freq_df = NULL, plot = TRUE,
                             facet = TRUE, textsize = 5) {
  if(is.null(group))
    group <- rep("G1", ncol(X))
  if(is.null(time))
    time <- 1:ncol(X)
  if(nrow(X) != nrow(cluster))
    stop("nrow(X) must match nrow(cluster)")
  if(ncol(X) != length(group))
    stop("ncol(X) must match length(group)")
  if(ncol(X) != length(time))
    stop("ncol(X) must match length(time)")

  mid.t <- sort(unique(time))
  mid.t <- mid.t[floor(length(mid.t)/2)]

  sample.data <- data.frame(sample = colnames(X), group, time,
                            stringsAsFactors = FALSE)
  data.long <- suppressMessages(
    join_sources(X, sample.data, cluster) %>%
      mutate(ts = paste0(group, "_", feature)) %>%
      data.frame(stringsAsFactors = FALSE)
  )
  # Compute cluster mean expression profile for each group
  clst.mean <- data.long %>%
    select(cluster, time, value, group) %>%
    group_by(cluster, group, time) %>%
    summarize_all(mean)

  if(!is.null(freq_df))
    clst.mean <- suppressMessages(clst.mean %>% left_join(freq_df))
  clst.mean <- data.frame(clst.mean, stringsAsFactors = FALSE)

  plt.lst <- lapply(unique(clst.mean$cluster), function(clst) {
    df <- data.long %>% filter(cluster == clst)
    df.mean <- clst.mean %>% filter(cluster == clst)
    df.text <- df.mean %>%
      filter(time == df.mean$time[1]) %>%
      mutate(max_y = max(df$value))
    plt <- ggplot(df, aes(x = time, y = value, color = group)) +
      geom_point() +
      geom_line(aes(group = ts), alpha = 0.2) +
      geom_line(data = df.mean, lwd = 2,
                aes(x = time, y = value, group = group)) +
      theme(legend.position="none")
    if("freq" %in% colnames(df.text)){
      plt <- plt + geom_text(data = df.text, size = textsize,
                             aes(x = mid.t, y = 0.9*max_y,
                                 label = paste0(cluster, ": ", freq)))
    }
    if(facet)
      plt <- plt + facet_wrap(~group, ncol = 4)
    return(plt)
  })
  if(plot) {
    do.call("grid.arrange", c(plt.lst, ncol=2))
  }
  names(plt.lst) <- unique(clst.mean$cluster)
  return(plt.lst)
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


