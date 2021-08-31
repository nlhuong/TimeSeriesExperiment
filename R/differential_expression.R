#' @title Differential timepoint expression testing.
#'
#' @description This is a wrapper around \code{\link[limma]{lmFit}}
#' and \code{\link[limma]{voom}} functions from \code{limma} package
#' for testing differential expression at specified timepoints.
#'
#' @param object A \code{TimeSeriesExperiment} object.
#' @param timepoints Vector of timepoints to test at.
#' @param min_gene_sum A scalar for filtering sparse genes before DE testing.
#' Default is 1.
#' @param alpha A scalar for level of significance. Default is 0.05.
#'
#' @return a \code{TimeSeriesExperiment} object with timepoint differential 
#' expression testing results stored in 'timepoint_de' element in 
#' \code{diff.expr} slot.
#'
#' @importFrom dplyr filter rename
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom stats model.matrix
#' @importFrom  methods slot<-
#' @importFrom methods validObject
#' @importFrom SummarizedExperiment assays rowData
#' @export
#' @examples
#' data("endoderm_small")
#' endoderm_small <- timepointDE(endoderm_small, timepoint = 1.0)
#' head(differentialExpression(endoderm_small, "timepoint_de")$`1`)
#'
timepointDE <- function(object, timepoints = "all", 
                         min_gene_sum = 1, alpha = 0.05) 
{
    if (!is(object, "TimeSeriesExperiment")) 
        stop("Input must be a 'TimeSeriesExperiment' object.")
    if (!validObject(object)) 
        stop("Invalid TimeSeriesExperiment object.")
    feature <- group <- time <- NULL
    if (any(timepoints == "all")) {
        timepoints <- sort(unique(timepoints(object)))
    }
    if(!all(timepoints %in% unique(timepoints(object)))) {
        stop("One or more entries of 'timepoints' not found in ",
             "timepoints(object).")
    }
    # we rename "group" to "condition" beacuse in 'model.matrix()' misuses
    # 'group' in design argument
    sample_data <- colData(object) %>%
        as.data.frame() %>%
        rename(condition = group)
    feature_data <- rowData(object) %>%
        as.data.frame()

    limma.res <- lapply(timepoints, function(t) {
        message("testing timepoint: ", t)
        smp.t <- sample_data %>% filter(time == t)
        x.t <- assays(object)$raw[, smp.t$sample]
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
        v <- voom(counts = dge, design = dsgn, plot = FALSE)

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

    diff_exp <- differentialExpression(object)
    diff_exp$timepoint_de <- limma.res
    slot(object, name = "differentialExpression", check = TRUE) <- diff_exp
    return(object)
}


#' @title Differential trajectory testing.
#'
#' @description Performs differential trajectory testing for timecourse
#' data using \code{\link[vegan]{adonis}} method.
#'
#' @param object A \code{TimeSeriesExperiment} object.
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
#' @importFrom SummarizedExperiment rowData
#' @export
#' @examples
#' data("endoderm_small")
#' endoderm_small <- makeTimeSeries(endoderm_small)
#' \dontrun{
#'    endoderm_small <- trajectoryDE(endoderm_small)
#'    head(differentialExpression(endoderm_small, "trajectory_de"))
#' }
trajectoryDE <- function(object, dist_method = "euclidean", 
                         p_adj_method = "BH", lambda = c(0.5, 0.25), 
                         verbose = TRUE, ...) 
{
    if (!is(object, "TimeSeriesExperiment")) 
        stop("Input must be a 'TimeSeriesExperiment' object.")
    feature <- Df <- pval <- R2 <- NULL
    if (!"ts_trans" %in% names(timeSeries(object))) {
        object <- makeTimeSeries(object)
    }
    if (!any(grepl("Lag_", colnames(timeSeries(object)[["ts_trans"]])))) {
        object <- addLags(object, lambda = lambda)
    }
    message("Testing differential feature trajectories...")
    ts_with_lags <- timeSeries(object, "ts_trans")
    feature_names <- unique(ts_with_lags$feature)
    adonis.res <- lapply(seq_along(feature_names), function(i) {
        feat <- feature_names[i]
        if((i %% 500) == 0 && verbose) message("Feature: ", i)
        itc <- ts_with_lags %>% filter(feature == feat)
        iadonis <- suppressMessages(
            vegan::adonis(
              formula = itc %>% dplyr::select(-(`feature`:`replicate`)) ~ group,
              data = itc %>% select(`feature`:`replicate`),
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
    
    adonis.res.df <- suppressMessages(
      adonis.res.df %>%
          left_join(as.data.frame(rowData(object)))
    )
    diff_exp <- differentialExpression(object)
    diff_exp[["trajectory_de"]] <- adonis.res.df
    slot(object, name = "differentialExpression", check = TRUE) <- diff_exp
    return(object)
}

