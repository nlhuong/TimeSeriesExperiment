
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
#' @return a list with the \code{hclust} object and the clust.assignment data.frame.
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
