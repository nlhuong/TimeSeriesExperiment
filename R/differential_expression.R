#' @title Differential timepoint expression testing.
#'
#' @description This is a wrapper around \code{\link[lmFit]{limma}}
#' and \code{\link[limma]{voom}} functions from \code{limma} package
#' for testing differential expression at specified timepoints.
#'
#' @param object A \code{vistimeseq} object.
#' @param timepoints Vector of timepoints to test at.
#' @param min_gene_sum A scalar for filtering sparse genes before DE testing.
#' Default is 1.
#' @param alpha A scalar for level of significance. Default is 0.05.
#'
#' @return a \code{vistimeseq} object with timepoint differential expression
#' testing results stored in 'timepoint_de' element in \code{diff.expr} slot.
#'
#' @importFrom dplyr filter rename
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom stats model.matrix
#' @importFrom  methods slot<-
#' @importFrom methods validObject
#' @export
#' @examples
#' data("endoderm_small")
#' endoderm_small <- timepoint_de(endoderm_small, timepoint = 1.0)
#' head(get_diff_expr(endoderm_small, "timepoint_de")$`1`)
#'
timepoint_de <- function(
    object, timepoints = "all", min_gene_sum = 1, alpha = 0.05) {
    feature <- group <- time <- NULL
    if (!validObject(object)){
        stop("Invalid 'vistimeseq' object.")
    }
    if (any(timepoints == "all")) {
        timepoints <- sort(unique(get_time(object)))
    }
    if(!all(timepoints %in% unique(get_time(object)))) {
      stop("One or more entries of \"timepoints\" not found in 'object@time'.")
    }
    # we rename "group" to "condition" beacuse in 'model.matrix()' misuses
    # 'group' in design argument
    sample_data <- sample_data(object) %>%
        rename(condition = group)
    feature_data <- feature_data(object)

    limma.res <- lapply(timepoints, function(t) {
        message("testing timepoint: ", t)
        smp.t <- sample_data %>% filter(time == t)
        x.t <- get_data(object, raw = TRUE)[, smp.t$sample]
        # filter out very sparse genes:
        x.t <- x.t[rowSums(x.t) >= min_gene_sum, ]  
        gdata.t <- feature_data %>% filter(feature %in% rownames(x.t))

        # Construct a DGEList object
        dge <- DGEList(counts = x.t, genes = gdata.t, samples = smp.t)

        # Compute size factor (lib sizes)
        dge <- calcNormFactors(dge)

        # Construct a design model matrix
        dsgn <- model.matrix(~ condition, dge$samples)

        # Voom normalization
        v <- voom(dge, dsgn, plot = FALSE)

        # Fit the model
        fit <- lmFit(v, dsgn)

        # Compute test statistics
        fit <- eBayes(fit)

        # Statistics table
        t.res.limma <- suppressMessages(
            topTable(fit, adjust.method="BH", number = Inf, p.value = alpha)
        )
        return(t.res.limma)
    })
    names(limma.res) <- timepoints

    diff_exp <- get_diff_expr(object, "all")
    diff_exp[["timepoint_de"]] <- limma.res
    slot(object, name = "diff.expr", check = TRUE) <- diff_exp
    return(object)
}


#' @title Differential trajectory testing.
#'
#' @description Performs differential trajectory testing for timecourse
#' data using \code{\link[vegan]{adonis}} method.
#'
#' @param object A \code{vistimeseq} object.
#' @param dist_method the name of any method used in vegdist to calculate
#' pairwise distances, "euclidean" by defaults.
#' @param p_adj_method a correction method. See details in
#' \code{\link[stats]{p.adjust}}. Default is "BH".
#' @param lambda Weights for each lag difference, for time-course data.
#' Length of \code{lambda} specifies number of lags to include. Default
#' is c(0.5, 0.25) for lag 1 and 2. Used only if 'timecourse.data' slot not
#' initialized.
#' @param verbose whether code comments should be printed. Default is TRUE.
#' @param ... other options to \code{\link[vegan]{adonis}} function from
#' \code{vegan}.
#'
#' @return a data.frame with adonis results for all features.
#'
#' @importFrom vegan adonis
#' @importFrom dplyr filter rename arrange mutate bind_rows
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom stats p.adjust
#' @importFrom  methods slot<-
#' @importFrom methods validObject
#' @export
#' @examples
#' data("endoderm_small")
#' endoderm_small <- normalize_data(endoderm_small)
#' endoderm_small <- trajectory_de(endoderm_small)
#' head(get_diff_expr(endoderm_small, "trajectory_de"))
#'
trajectory_de <- function(
    object, dist_method = "euclidean", p_adj_method = "BH", 
    lambda = c(0.5, 0.25), verbose = TRUE, ...) {
    feature <- Df <- pval <- R2 <- NULL
    if (!validObject(object)){
        stop("Invalid 'vistimeseq' object.")
    }
    if (is.null(time_course(object))) {
        object <- convert_to_timecourse(object)
        message("Converted to timecourse format.")
        object <- add_lags(object, lambda = lambda)
        message("Added lags with coefficients: ",
                paste0(lambda, collapse = " "))
    }

    message("Testing differential feature trajectories.")
    tc <- time_course(object)
    feature_names <- unique(tc$feature)
    adonis.res <- lapply(seq_along(feature_names), function(i) {
        feat <- feature_names[i]
        if((i %% 500) == 0 && verbose) {
            message("Feature: ", i)
        }
        itc <- tc %>%
            filter(feature == feat)
        iadonis <- suppressMessages(
            adonis(
              formula = itc %>% select(-(feature:replicate)) ~ group,
              data = itc %>% select(feature:replicate),
              method = dist_method, ...)
        )
        ires <- iadonis$aov.tab[1, ]
        ires$feature <- feature_names[i]
        return(ires)
    })

    adonis.res.df <- bind_rows(adonis.res) %>%
        rename("pval" = "Pr(>F)") %>%
        select(feature, Df:pval) %>%
        arrange(-R2) %>%
        mutate(p.adj = p.adjust(pval, method = p_adj_method))

    diff_exp <- get_diff_expr(object, type = "all")
    diff_exp[["trajectory_de"]] <- adonis.res.df
    slot(object, name = "diff.expr", check = TRUE) <- diff_exp
    return(object)
}

