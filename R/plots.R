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
#' @importFrom grDevices colorRampPalette
#' @importFrom methods validObject
#' @importFrom stats sd
#' @export
#' @examples
#' endoderm_small
#' endoderm_small <- normalize_data(endoderm_small)
#' \dontrun{
#' plot_heatmap(endoderm_small)
#' }
plot_heatmap <- function(
    object, num.feat = 200, scale = TRUE, feat_desc = "feature",
    sample_desc = "sample"){
    group <- time <- NULL
    if (!validObject(object)) {
        stop("Invalid vistimeseq object.")
    }
    if (!feat_desc %in% colnames(feature_data(object))){
        stop("\"feat_desc\" not in object@feature.data")
    }
    if (!sample_desc %in% colnames(sample_data(object))){
        stop("\"sample_desc\" not in object@sample.data")
    }
    cnts <- get_data(object)
    rownames(cnts) <- feature_data(object)[, feat_desc]
    colnames(cnts) <- sample_data(object)[, sample_desc]
    top_feat <- apply(cnts, 1, sd)
    top_feat <- names(top_feat[order(-top_feat)])[seq_len(num.feat)]
    cols_ordered <- order(get_group(object), get_replicate(object),
                                                get_time(object))
    Y <- cnts[top_feat, cols_ordered]
    if (scale) {
        Y <- t(scale(t(Y), center = TRUE, scale = TRUE))
    }

    smpdf <- sample_data(object) %>%
        select(group, replicate, time) %>%
        arrange(group, replicate, time) %>%
        mutate(time = as.numeric(time))

    n_group <- length(unique(get_group(object)))
    cols <- brewer.pal(9, name = "Set1")[seq_len(min(9, n_group))]
    group_cols <- colorRampPalette(colors = cols)(n_group)
    names(group_cols) <- unique(get_group(object))

    n_replicates <- length(unique(get_replicate(object)))
    rep_cols <- colorRampPalette(colors = brewer.pal(8, name = "Set3")[
        seq_len(min(8, n_replicates))])(n_replicates)
    names(rep_cols) <- unique(get_replicate(object))

    time_cols <- colorRamp2(
        breaks = seq(min(get_time(object)), max(get_time(object)), 
                     length.out = 10),
        colors =    viridis(10))

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
#' @importFrom ggplot2 ggplot aes aes_string geom_point geom_hline
#' @importFrom ggplot2 geom_vline xlab ylab coord_fixed
#' @importFrom viridis viridis
#' @importFrom dplyr left_join
#' @importFrom tibble rownames_to_column column_to_rownames
#' @export
#' @examples
#' endoderm_small
#' endoderm_small <- run_pca(endoderm_small)
#' plot_sample_pca(endoderm_small, col.var = "group")
#'
plot_sample_pca <- function(object, axis = c(1, 2), col.var = NULL, ...) {
    if(is.null(get_dim_reduced(object, "pca_sample"))) {
        stop("No 'pca_sample' in object@dim.red. Run PCA for samples first.")
    }
    axis <- if(is.numeric(axis)) paste0("PC", axis) else axis
    pca.sample <- get_dim_reduced(object, "pca_sample")
    pca.eigs <- get_dim_reduced(object, "pca_eigs")

    pca.scores <- pca.sample[, axis] %>%
        as.data.frame() %>%
        rownames_to_column("sample")

    if(all(pca.scores$sample %in% sample_names(object))) {
        pca.scores <- suppressMessages(
            pca.scores %>%
                left_join(sample_data(object)) %>%
                column_to_rownames("sample")
        )
    } else if (all(pca.scores$sample %in% 
                   collapsed_sample_data(object)$sample)) {
        pca.scores <- suppressMessages(
            pca.scores %>%
                left_join(collapsed_sample_data(object)) %>%
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
        plt <- plt + viridis::scale_fill_viridis()
    }
    return(plt)
}


#' @title Overlay (time) series over PCA grid
#' @description PCA plot for data features, with time-series levels overlayed
#' on top.
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
#' @importFrom graphics plot legend par
#' @importFrom proxy dist
#' @importFrom viridis viridis
#' @importFrom methods validObject
#'
#' @return None
#' @export
#' @examples
#' endoderm_small
#' endoderm_small <- run_pca(endoderm_small)
#' plot_ts_pca(endoderm_small)
#'
plot_ts_pca <- function(
    object, axis = c(1, 2), m = 20, n = 20, group.highlight = NULL,
    linecol = NULL, ...) {
    feature <- group <- NULL
    if (!validObject(object)){
        stop("Invalid vistimeseq object.")
    }
    if(is.null(get_dim_reduced(object, type = "pca_feature"))) {
        stop("No 'pca_feature' in object@dim.red. Run PCA for features first.")
    }
    if (is.null(time_course(object, collapsed = TRUE))) {
        object <- collapse_replicates(object, FUN = mean)
        object <- convert_to_timecourse(object)
        message("Aggregated data over replicates and converted",
                        "to time-course format.")
    }
    # Prepare scores data
    tc <- time_course(object, collapsed = TRUE)
    tc <- tc[, !grepl("Lag_", colnames(tc))]
    tmp <- as.numeric(colnames(tc %>% select(-feature, -group, -replicate)))
    axis <- if(is.numeric(axis)) paste0("PC", axis) else axis
    pca.feature <- get_dim_reduced(object, type = "pca_feature")
    pca.eigs <- get_dim_reduced(object, type = "pca_eigs")
    pca.loadings <- suppressMessages(
        pca.feature[, axis] %>%
        as.data.frame() %>%
        rownames_to_column("feature") %>%
        left_join(feature_data(object)) %>%
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
        groups.unique <- unique(get_group(object))
    } else {
        groups.unique <- group.highlight
    }
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
                 groups.unique, col = linecol, xpd = TRUE,
                 lty = c(1,1), lwd = c(3,3))
}

#' @title Plot (time) series over clusters.
#'
#' @description Plots timecourse (aggregated over replicates) feature data
#' faceted by computed clusters and experimental group.
#'
#' @param object A \code{vistimeseq} object.
#' @param features A vector of names of selected features to plot.
#' @param transparency transparency of trajectory lines.
#' @param ncol number of columns in the facte plot.
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
#' endoderm_small <- cluster_timecourse_features(endoderm_small)
#' plot_ts_clusters(endoderm_small)
#'
plot_ts_clusters <- function(
    object, features = NULL, transparency = 0.5, ncol = 4) {
    feature <- cluster <- freq <- group <- time <- value <- category <- NULL
    if (!validObject(object)){
        stop("Invalid vistimeseq object.")
    }
    if(is.null(get_cluster_map(object))) {
        stop("No 'cluster_map' in object@cluster.features. Perform
                 clustering with 'cluster_timecourse_features()' first.")
    }
    if (is.null(time_course(object, collapsed = TRUE))) {
        stop("No 'tc_collapsed' in object@timecourse.data, but ",
             "'cluster.features' slot is non-emty. Check if data 'vistimeseq' ",
             "object is corrupted.")
    }
    if (is.null(features)){
        features <- feature_names(object)
    }
    if (!all(features %in% feature_names(object))) {
        stop("One or more feature in \"features\" not found in ",
                 "'object@feature.names'")
    }
    cluster_map <- get_cluster_map(object) %>%
        filter(feature %in% features)

    freq_df <- cluster_map %>%
        group_by(cluster) %>%
        summarise(freq = n()) %>%
        arrange(desc(freq))

    groups <- unique(get_group(object))
    cat_levels <- paste0("[", freq_df$cluster, ": ", freq_df$freq, "]")
    cat_levels <- paste(rep(cat_levels, each = length(groups)),
                        rep(groups, length(cat_levels)))
    tc_data <- suppressMessages(
        time_course(object, collapsed = TRUE) %>%
            filter(feature %in% features) %>%
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
            aes(group = feature), alpha = transparency
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
#' @description Plotting expression over time for selected genes curve
#' and colors correspond to distinct groups.
#'
#' @param object A \code{vistimeseq} object.
#' @param features A vector of names of selected features to plot.
#' @param smooth If TRUE a smoothed line is plotted for each gene
#' and each group, else a piecewise linear average (over replicates) curve
#' is plotted.
#' @param ncol An integer indicating the number of columns for facetting.
#' Default is 5.
#'
#' @return list of ggplot objects
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_smooth facet_wrap
#' @importFrom dplyr filter select mutate left_join group_by summarise_all
#' @importFrom tidyr gather
#' @export
#' @examples
#' feat_to_plot <- feature_names(endoderm_small)[1:10]
#' plot_time_series(endoderm_small, features = feat_to_plot)
#'
plot_time_series <- function(
    object, features = feature_names(object), smooth = TRUE, ncol = 5){
    feature <- symbol <- time <- value <- group <- category <- NULL
    if(!all(features %in% feature_names(object))){
        stop("\"features\" must be a subset of object@feature.names")
    }
    feature_data <- feature_data(object) %>% filter(feature %in% features)

    if(!"symbol" %in% colnames(feature_data)){
        feature_data$symbol <- feature_data$feature
    }

    tc_data <- data_to_tc(
        X = get_data(object)[features, ], time = get_time(object),
        group = get_group(object), replicate = get_replicate(object))

    tc_data <- suppressMessages(
        tc_data %>%
            gather(key = "time", value = "value", -(feature:replicate)) %>%
            left_join(feature_data %>% select(feature, symbol)) %>%
            mutate(
                symbol = factor(symbol, levels = feature_data$symbol),
                time = as.numeric(time))
    )
    tc_data <- tc_data %>%
        mutate(category = paste0(group, "_", replicate))
    plt <- ggplot(
        tc_data,
        aes(x = time, y = value, color = group)) +
        geom_point(size = 1) +
        facet_wrap(~ symbol, scales = "free", ncol = ncol)

    if(length(unique(tc_data$replicate)) > 1) {
        plt <- plt + geom_line(aes(group = category), lty = 3, alpha = 0.7)
    }
    if(smooth) {
        plt <- plt + geom_smooth(aes(x = time), lwd = 1.5)
    } else {
        tc_data_mean <- suppressMessages(
            tc_data %>%
                select(-replicate) %>%
                group_by(feature, group, time) %>%
                summarise(value = mean(value)) %>%
                left_join(feature_data)
        )
        plt <- plt + geom_line(data = tc_data_mean, lwd = 1.5)
    }
    return(plt)
}


#' @title Plot enrichment results.
#' @description Plotting top most enriched terms found with DE methods
#' and tested for overrepresentation in GO/KEGG db using
#' goana/kegga from limma package.
#'
#' @param enrich a data matrix or data.frame with enrichment result -
#' outputs from \code{\link{pathway_enrichment}} function or
#' \code{\link[limma]{goana}}, \code{\link[limma]{kegga}}.
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
#' selected_genes <- c('114299', '2825', '3855', '221400', '7941',
#'                     '6164', '1292', '6161', '6144', '23521')
#' enrich_res <- pathway_enrichment(
#'   object = endoderm_small, clustered = FALSE,
#'   features = selected_genes,
#'   species = "Hs", ontology = "BP", fltr_DE = 0,
#'   fltr_N = Inf, fltr_P.DE = 0.05)
#' plot_enrichment(enrich = enrich_res, n_max = 15)
#'
plot_enrichment <- function(enrich, n_max = 15) {
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

