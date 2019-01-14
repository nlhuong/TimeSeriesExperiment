#' @title Plot a heatmap
#'
#' @description Generates a heatmap of expression data.
#'
#' @param object A \code{TimeSeriesExperiment} object
#' @param num.feat Number of top most features to use.
#' @param scale Whether to scale the data (by features) before plotting.
#' @param feat_desc One of the column names from \code{feature_data(object)}
#' to describe the features.
#' @param sample_desc ne of the column names from \code{sample_data(object)}
#' to describe the samples.
#' @param ... other parameters ggplot.
#'
#' @return Returns a \code{ggplot2} objet.

#' @importFrom SummarizedExperiment assays rowData
#' @importFrom methods validObject
#' @importFrom viridis viridis
#' @importFrom stats sd
#' @export
#' @examples
#' data("endoderm_small")
#' endoderm_small <- normalizeData(endoderm_small)
#' \dontrun{
#'     plotHeatmap(endoderm_small)
#' }
plotHeatmap <- function(object, num.feat = 200, scale = TRUE, 
                        feat_desc = "feature", sample_desc = "sample")
{   
  
    if (!is(object, "TimeSeriesExperiment")) 
        stop("Input must be a 'TimeSeriesExperiment' object.")
    if (!validObject(object)) 
        stop("Invalid TimeSeriesExperiment object.")
    if (!feat_desc %in% colnames(rowData(object)))
        stop("'feat_desc' not in colnames(rowData(object))")
    if (!sample_desc %in% colnames(colData(object)))
        stop("'sample_desc' not in colnames(colData(object))")
  
    pkgs_needed <- c("ComplexHeatmap", "circlize", "grid", "grDevices",
                     "RColorBrewer")
    pkgs_missing <- setdiff(pkgs_needed, installed.packages()) 
    if (length(pkgs_missing) > 0) {
      stop("Packages:", paste(pkgs_missing, collapse = ", "), 
           " needed for this function to work. ",
           "Please install them.", call. = FALSE)
    }
    
    group <- timepoint <- NULL
    if(!"norm" %in% names(assays(object))) {
        cnts <- assays(object)$raw
    } else {
        cnts <- assays(object)$norm
    }
    rownames(cnts) <- rowData(object)[, feat_desc]
    colnames(cnts) <- colData(object)[, sample_desc]
    top_feat <- apply(cnts, 1, sd)
    top_feat <- names(top_feat)[order(-top_feat)[seq_len(num.feat)]]
    cols_ordered <- order(groups(object), replicates(object), 
                          timepoints(object))
    Y <- cnts[top_feat, cols_ordered]
    if (scale) {
        Y <- t(scale(t(Y), center = TRUE, scale = TRUE))
    }
    smpdf <- colData(object) %>%
        as.data.frame() %>%
        select(group, replicate, timepoint) %>%
        arrange(group, replicate, timepoint) %>%
        mutate(timepoint = as.numeric(timepoint))

    n_group <- length(unique(groups(object)))
    cols <- RColorBrewer::brewer.pal(9, name = "Set1")[
      seq_len(min(9, n_group))]
    group_cols <- grDevices::colorRampPalette(colors = cols)(n_group)
    names(group_cols) <- unique(groups(object))

    n_replicates <- length(unique(replicates(object)))
    cols <- RColorBrewer::brewer.pal(8, name = "Set3")[
      seq_len(min(8, n_replicates))]
    rep_cols <- grDevices::colorRampPalette(colors = cols)(n_replicates)
    names(rep_cols) <- unique(replicates(object))

    time_cols <- circlize::colorRamp2(
        breaks = seq(min(timepoints(object)), max(timepoints(object)),
                     length.out = 10),
        colors = viridis(10))

    ha1 <- ComplexHeatmap::HeatmapAnnotation(
        df = smpdf, 
        col = list(group = group_cols, replicate = rep_cols, 
                   timepoint = time_cols))

    ComplexHeatmap::Heatmap(
        Y, name = "Z-score", cluster_columns = FALSE,
        top_annotation = ha1, row_names_gp = grid::gpar(fontsize = 8))
}


#' @title Plot a standard PCA
#'
#' @description Generates a standard PCA plot of observations in the dataset.
#'
#' @param object A \code{TimeSeriesExperiment} object
#' @param axis An integer vector indicating principal components to use for
#' plotting, by default 1:2.
#' @param col.var A character string indicating a column
#' from sample_data(object) which should be used for coloring
#' the points. By default NULL.
#' @param ... other parameters ggplot.
#'
#' @return Returns a \code{ggplot2} objet.
#'
#' @importFrom ggplot2 ggplot aes aes_string geom_point geom_hline
#' @importFrom ggplot2 geom_vline xlab ylab coord_fixed
#' @importFrom viridis scale_fill_viridis
#' @importFrom dplyr left_join
#' @importFrom tibble rownames_to_column column_to_rownames
#' @export
#' @examples
#' data("endoderm_small")
#' endoderm_small <- runPCA(endoderm_small)
#' plotSamplePCA(endoderm_small, col.var = "group")
#'
plotSamplePCA <- function(object, axis = c(1, 2), col.var = NULL, ...) {
    if (!is(object, "TimeSeriesExperiment")) 
        stop("Input must be a 'TimeSeriesExperiment' object.")
    if (!validObject(object)) 
        stop("Invalid TimeSeriesExperiment object.")
    if(is.null(dimensionReduction(object, "pca_sample"))) 
        stop("No 'pca_sample' available. Run PCA for samples first.")
    
    axis <- if(is.numeric(axis)) paste0("PC", axis) else axis
    pca.sample <- dimensionReduction(object, "pca_sample")
    pca.eigs <- dimensionReduction(object, "pca_eigs")

    pca.scores <- pca.sample[, axis] %>%
        as.data.frame() %>%
        rownames_to_column("sample")

    if(all(pca.scores$sample %in% colnames(object))) {
        pca.scores <- suppressMessages(
            pca.scores %>%
                left_join(as.data.frame(colData(object)))%>%
                column_to_rownames("sample")
        )
    } else if (all(pca.scores$sample %in% colDataCollapsed(object)$sample)) {
        pca.scores <- suppressMessages(
            pca.scores %>%
                left_join(as.data.frame(colDataCollapsed(object)))%>%
                column_to_rownames("sample")
        )
    } else {
        stop("Sample names in sample data and PCA coordinates disagree.")
    }

    axis_label <- paste0(
        colnames(pca.scores)[c(1, 2)], " [",
        signif(pca.eigs[c(1, 2)]/sum(pca.eigs)*100, 3), "%]"
    )
    plt <- ggplot(
        data = pca.scores,
        aes(x = pca.scores[[1]], y = pca.scores[[2]])
        ) +
        geom_point(
          aes_string(fill = col.var), color = "grey80", pch = 21, ...) +
        geom_hline(aes(yintercept =0), size=.2) +
        geom_vline(aes(xintercept = 0), size=.2) +
        xlab(axis_label[1]) +
        ylab(axis_label[2]) +
        coord_fixed(1)  # ratio must reflect variances of new PCs from prcomp

    if(all(!is.null(pca.scores), is.numeric(pca.scores[, col.var]))){
        plt <- plt + scale_fill_viridis()
    }
    return(plt)
}


#' @title Overlay (time) series over PCA grid
#' @description PCA plot for data features, with time-series levels overlayed
#' on top.
#'
#' @param object A \code{TimeSeriesExperiment} object.
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
#' @importFrom graphics plot legend par
#' @importFrom methods validObject
#' @importFrom viridis viridis
#' @importFrom SummarizedExperiment rowData
#'
#' @return None
#' @export
#' @examples
#' data("endoderm_small")
#' endoderm_small <- runPCA(endoderm_small)
#' plotTimeSeriesPCA(endoderm_small)
#'
plotTimeSeriesPCA <- function(object, axis = c(1, 2), m = 20, n = 20, 
                              group.highlight = NULL, linecol = NULL, ...) 
{
    feature <- group <- NULL
    if (!is(object, "TimeSeriesExperiment")) 
        stop("Input must be a 'TimeSeriesExperiment' object.")
    if (!validObject(object)) 
        stop("Invalid TimeSeriesExperiment object.")
    if(is.null(dimensionReduction(object, "pca_sample"))) 
        stop("No 'pca_sample' available. Run PCA for samples first.")
    if (is.null(timeSeries(object, "ts_collapsed"))) {
        object <- collapseReplicates(object)
        object <- makeTimeSeries(object)
    }
    
    pkgs_needed <- c("Hmisc", "proxy")
    pkgs_missing <- setdiff(pkgs_needed, installed.packages()) 
    if (length(pkgs_missing) > 0) {
      stop("Packages:", pkgs_missing, " needed for this function to work. ",
           "Please install them.", call. = FALSE)
    }
    
    # Prepare scores data
    ts <- timeSeries(object, "ts_collapsed")
    ts <- ts[, !grepl("Lag_", colnames(ts))]
    tmp <- as.numeric(colnames(ts %>% select(-feature, -group, -replicate)))
    
    axis <- if(is.numeric(axis)) paste0("PC", axis) else axis
    pca.feature <- dimensionReduction(object, "pca_feature")
    pca.eigs <- dimensionReduction(object, "pca_eigs")
    pca.loadings <- suppressMessages(
        pca.feature[, axis] %>%
        as.data.frame() %>%
        rownames_to_column("feature") %>%
        left_join(as.data.frame(rowData(object)))%>%
        column_to_rownames("feature")
    )
    colnames(pca.loadings)[c(1, 2)] <-
        paste0(colnames(pca.loadings)[c(1, 2)], " [",
                     signif(pca.eigs[c(1, 2)]/sum(pca.eigs)*100, 3), "%]")

    # Create a grid over score values
    mins <- apply(pca.loadings[, c(1, 2)], 2, min)
    maxes <- apply(pca.loadings[, c(1, 2)], 2, max)
    x <- seq(mins[1], maxes[1], length.out = m)
    y <- seq(mins[2], maxes[2], length.out = n)
    dx <- x[2] - x[1]; dy <- y[2] - y[1]
    grid <- expand.grid(x, y)

    # Find gene closest to the grid center
    xD <- proxy::dist(grid[, 1], pca.loadings[, 1])
    yD <- proxy::dist(grid[, 2], pca.loadings[, 2])
    D <- proxy::dist(grid, pca.loadings[, c(1, 2)])

    min_dists <- apply(D, 1, min)
    min_dists_ix <- apply(D, 1, which.min)
    x_min_dists <- vapply(seq_len(nrow(xD)),
                          function(i) xD[i, min_dists_ix[i]], numeric(1))
    y_min_dists <- vapply(seq_len(nrow(yD)),
                          function(i) yD[i, min_dists_ix[i]], numeric(1))
    min_dists_ix[x_min_dists > dx/2 | y_min_dists > dy/2] <- NA

    # Plot all points corresponding to each feature
    par(mar=par()$mar * c(2, 1.2, 1.2, 1.2), xpd = TRUE,
        cex = 0.7, cex.main = 2, cex.axis = 1.5, cex.lab = 1.5)
    plot(pca.loadings[, c(1, 2)], type = "p", pch = 16,
             xlim = c(mins[1] - dx/2, maxes[1] + dx/2),
             ylim = c(mins[2] - dy/2, maxes[2] + dy/2),
             asp=1,  # ratio must reflect variances of new PC from prcomp
             ...)

    if (is.null(group.highlight)){
        groups.unique <- unique(groups(object))
    } else {
        groups.unique <- group.highlight
    }
    if(is.null(linecol)) {
        linecol <- viridis::viridis(length(groups.unique))
        names(linecol) <- groups.unique
    }
    ylimits <- c(min(ts %>% select(-feature, -group, -replicate)),
                 max(ts %>% select(-feature, -group, -replicate)))

    # Plot all the time-course profiles
    for(i in seq_along(min_dists_ix)) {
        igene <- rownames(pca.feature)[min_dists_ix[i]]
        if(is.na(min_dists_ix[i])) next
        for (gr in groups.unique) {
            gTC <- ts %>%
                filter(feature == igene, group == gr) %>%
                select(-feature, -group, -replicate) %>%
                as.numeric()
            Hmisc::subplot(
                plot(tmp, gTC, type = "l", lwd = 2,
                     col =  linecol[gr], frame = FALSE, axes = FALSE,
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
            names(linecol), col = linecol, xpd = TRUE,
            lty = c(1,1), lwd = c(3,3))
}

#' @title Plot (time) series over clusters.
#'
#' @description Plots timecourse (aggregated over replicates) feature data
#' faceted by computed clusters and experimental group.
#'
#' @param object A \code{TimeSeriesExperiment} object.
#' @param features A vector of names of selected features to plot.
#' @param transparency transparency of trajectory lines.
#' @param ncol number of columns in the factet plot.
#' @param scales character scalar indecating facet scales, by default "free".
#'
#' @return ggplot object
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_smooth facet_wrap
#' @importFrom dplyr filter left_join mutate select group_by
#' @importFrom dplyr summarize_all arrange desc n contains
#' @importFrom tidyr gather
#' @importFrom methods validObject
#'
#' @export
#' @examples
#' data("endoderm_small")
#' endoderm_small <- clusterTimeSeries(endoderm_small)
#' plotTimeSeriesClusters(endoderm_small)
#'
plotTimeSeriesClusters <- function(object, features = NULL,
                                   transparency = 0.5, ncol = 4,
                                   scales = "free") 
{
    feature <- cluster <- freq <- group <- timepoint <- 
      value <- category <- used_for_hclust <- NULL
    if (!is(object, "TimeSeriesExperiment")) 
        stop("Input must be a 'TimeSeriesExperiment' object.")
    if (!validObject(object))
        stop("Invalid TimeSeriesExperiment object.")
    if(!"final_cluster_map" %in% names(clusterAssignment(object))){
        stop("No clustering results not available. Perform ",
              "clustering with 'clusterTimeSeries()' first.")
    }
    if (is.null(timeSeries(object, "ts_collapsed"))) {
        stop("'ts_collapsed' not in 'timeSeries' slot, but ",
             "'cluster.features' slot is non-emty. Check if ", 
             "'TimeSeriesExperiment' object is corrupted.")
    }
    if (is.null(features)){
        features <- rownames(object)
    }
    if (!all(features %in% rownames(object))) {
        stop("One or more feature in 'features' not found in ",
             "rownames(object).")
    }
    cluster_map <- clusterMap(object) %>%
        filter(feature %in% features) %>%
        select(-used_for_hclust)

    freq_df <- cluster_map %>%
        group_by(cluster) %>%
        summarise(freq = n()) %>%
        arrange(desc(freq))

    groups.labs <- unique(groups(object))
    cat_levels <- paste0("[", freq_df$cluster, ": ", freq_df$freq, "]")
    cat_levels <- paste(rep(cat_levels, each = length(groups.labs)),
                        rep(groups.labs, length(cat_levels)))
    ts_data <- suppressMessages(
        timeSeries(object, "ts_collapsed") %>%
            filter(feature %in% features) %>%
            select(-replicate, -contains("Lag_")) %>%
            left_join(cluster_map) %>%
            gather(
                key = "timepoint", value  = "value",
                -feature, -group, -cluster
            ) %>%
            left_join(freq_df) %>%
            mutate(
                timepoint = as.numeric(timepoint),
                value = as.numeric(value),
                category = paste0("[", cluster, ": ", freq, "] ", group)
            ) %>%
            arrange(cluster, group, feature, timepoint) %>%
            mutate(
                category = factor(category, levels = cat_levels))
    )
    # Compute cluster mean expression profile for each group
    ts_cluster_mean <- suppressMessages(
        ts_data %>%
            select(-feature) %>%
            group_by(cluster, group, category, timepoint) %>%
            summarize_all(mean) %>%
            left_join(freq_df) %>%
            arrange(cluster, group)
    )
    plt <- ggplot(ts_data, aes(y = value , x = timepoint, color = group)) + 
        geom_line(aes(group = feature), alpha = transparency) +
        geom_point() +
        geom_line(
            data = ts_cluster_mean, lwd = 1.5, color = "grey50",
            aes(group = group)
        ) +
        facet_wrap(~category, scales = scales, ncol = ncol)
    return(plt)
}


#' @title Plot selected time series.
#' @description Plotting expression over time for selected genes curve
#' and colors correspond to distinct groups.
#'
#' @param object A \code{TimeSeriesExperiment} object.
#' @param features A vector of names of selected features to plot.
#' @param trans A boolean indicating whether (TRUE) transformed, variance 
#' stabilized, assay values should be printed or (FALSE) just normalized by
#' sample depth.
#' @param smooth If TRUE a smoothed line is plotted for each gene
#' and each group, else a piecewise linear average (over replicates) curve
#' is plotted.
#' @param ncol An integer indicating the number of columns for facetting.
#' Default is 5.
#' @param scales character scalar indecating facet scales, by default "free".
#'
#' @return list of ggplot objects
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_smooth facet_wrap
#' @importFrom dplyr filter select mutate left_join group_by summarise_all 
#' @importFrom dplyr starts_with
#' @importFrom tidyr gather
#' @importFrom SummarizedExperiment rowData
#' @export
#' @examples
#' data("endoderm_small")
#' feat_to_plot <- rownames(endoderm_small)[1:10]
#' plotTimeSeries(endoderm_small, features = feat_to_plot, smooth = FALSE)
#'
plotTimeSeries <- function(object, features = rownames(object), 
                           trans = FALSE, smooth = TRUE, ncol = 5,
                           scales = "free")
{
    feature <- symbol <- timepoint <- value <- group <- category <- NULL
    if (!is(object, "TimeSeriesExperiment")) 
        stop("Input must be a 'TimeSeriesExperiment' object.")
    if (!validObject(object))
        stop("Invalid TimeSeriesExperiment object.")
    if(!all(features %in% rownames(object)))
        stop("'features' must be a subset of rownames(object)")
    
    if (!"ts" %in% names(timeSeries(object))) {
        object <- makeTimeSeries(object)
    }
    feature_data <- rowData(object) %>% 
        as.data.frame() %>%
        filter(feature %in% features) %>%
        arrange(factor(feature, levels = features))
    
    if(!"symbol" %in% colnames(feature_data)){
        feature_data$symbol <- feature_data$feature
    }
    
    if(trans) {
        ts_data <- timeSeries(object, "ts_trans")
    } else {
        ts_data <- timeSeries(object, "ts")
    }
    ts_data <- suppressMessages(
        ts_data %>%
            filter(feature %in% features) %>%
            select(-starts_with("Lag_")) %>%
            gather(key = "timepoint", value = "value", -(feature:replicate)) %>%
            left_join(feature_data %>% select(feature, symbol)) %>%
            mutate(
                symbol = factor(symbol, levels = feature_data$symbol),
                timepoint = as.numeric(timepoint),
                category = paste0(group, "_", replicate)) 
    )
    plt <- ggplot(
        ts_data,
        aes(x = timepoint, y = value, color = group)) +
        geom_point(size = 1) +
        facet_wrap(~ symbol, scales = scales, ncol = ncol)

    if(length(unique(ts_data$replicate)) > 1) {
        plt <- plt + geom_line(aes(group = category), lty = 3, alpha = 0.7)
    }
    if(smooth) {
        plt <- plt + geom_smooth(aes(x = timepoint ), lwd = 1.5)
    } else {
        ts_data_mean <- suppressMessages(
            ts_data %>%
                select(-replicate) %>%
                group_by(feature, group, timepoint ) %>%
                summarise(value = mean(value)) %>%
                left_join(feature_data)
        )
        plt <- plt + geom_line(data = ts_data_mean, lwd = 1.5)
    }
    return(plt)
}


#' @title Plot enrichment results.
#' @description Plotting top most enriched terms found with DE methods
#' and tested for overrepresentation in GO/KEGG db using
#' goana/kegga from limma package.
#'
#' @param enrich a data matrix or data.frame with enrichment result -
#' outputs from \code{\link{pathwayEnrichment}} function or
#' \code{\link[limma]{goana}}, \code{\link[limma:goana]{limma::kegga()}}.
#' Must contain columns Term, DE, and P.DE.
#' @param n_max max number of terms to show
#'
#' @return a ggplot object
#'
#' @importFrom dplyr arrange mutate
#' @importFrom viridis scale_color_viridis
#' @importFrom ggplot2 ggplot aes geom_point
#' @export
#' @examples
#' data("endoderm_small")
#' selected_genes <- c('114299', '2825', '3855', '221400', '7941',
#'                     '6164', '1292', '6161', '6144', '23521')
#' enrich_res <- pathwayEnrichment(
#'   object = endoderm_small, clustered = FALSE,
#'   features = selected_genes,
#'   species = "Hs", ontology = "BP", fltr_DE = 0,
#'   fltr_N = Inf, fltr_P.DE = 0.05)
#' plotEnrichment(enrich = enrich_res, n_max = 15)
#'
plotEnrichment <- function(enrich, n_max = 15) {
    DE <- N <- P.DE <- Term <- NULL
    enrich <- enrich %>%
        arrange(-DE, P.DE) %>%
        mutate(
            Term = paste0(Term, " (", DE, "/", N , ")"),
            Term = factor(Term, levels = Term)
        )
    plt <- ggplot(
        enrich[seq(1, min(n_max, nrow(enrich))), ],
        aes(y = Term, x = -log10(P.DE), size = N, color = DE/N)
        ) +
        geom_point() +
        scale_color_viridis()
    return(plt)
}
