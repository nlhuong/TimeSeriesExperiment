#' @title Pathway enrichment testing.
#'
#' @description This is a wrapper around \code{\link[limma]{goana}} function
#' from \code{limma} package for testing over-representation of gene ontology
#'(GO) terms or KEGG pathways in sets of genes.
#'
#' @param object A \code{vistimeseq} object.
#' @param features A vector of names of selected features to plot.
#' @param species A character string  specifying the species.
#' See \code{\link[limma]{goana}} for details.
#' @param clustered Whether \code{features} should be grouped based on
#' cluster assignment (stored in \code{object@cluster.features$cluster_map})
#' and tested as separate sets.
#' @param kegg Whether KEGG pathways should be used instead of gene ontology
#' (GO)
#' @param ontology character vector of ontologies to be included in output.
#' Elements should be one or more of "BP", "CC" or "MF".
#' @param fltr_DE A scalar fraction of the number of genes in tested set
#' to use as a threshold for filtering genes based on "DE" column
#' (the number of genes in the DE set). Default is 0.1, i.e. at least
#' 0.1 genes in the set must be present in the pathway.
#' @param fltr_N A number of genes used as a threshold to filter out all
#' pathway terms of size greater than the threshold. Default is 500.
#' @param fltr_P.DE A p-value threshold to filter out terms in the enrichment
#' results. Default is 0.05.
#' @param ... other parameters for \code{\link[limma]{goana}} or
#' \code{\link[limma]{kegga}}.
#'
#' @return a \code{data.frame} or list of \code{data.frame}s with enrichment
#'results.
#'
#' @importFrom dplyr left_join filter arrange group_by summarise n
#' @importFrom limma goana kegga
#' @importFrom utils installed.packages
#' @export
#' @examples
#' \dontrun{
#' endoderm_small
#' endoderm_small <- normalize_data(endoderm_small)
#' endoderm_small <- trajectory_de(endoderm_small)
#' genes_with_de_trajectory <- get_diff_expr(endoderm_small, "trajectory_de") %>%
#'   filter(pval <= max(0.05, min(pval)), R2 > 0.7) %>%
#'   arrange(-R2)
#' res <- pathway_enrichment(
#'   object = endoderm_small, clustered = FALSE,
#'   features = genes_with_de_trajectory$feature,
#'   species = "Hs", fltr_DE = 0, fltr_N = Inf, fltr_P.DE = 0.05)
#' head(res)
#'}
#'
pathway_enrichment <- function(object, features, species,
                               clustered = TRUE, kegg = FALSE,
                               ontology = c("BP", "CC", "MF"),
                               fltr_DE = 0.1, fltr_N = 500, fltr_P.DE = 0.05,
                               ...){
  feature <- cluster <- DE <- N <- P.DE <- NULL
  if(all(clustered, is.null(get_cluster_map(object)))) {
    stop("No 'cluster_map' in object@cluster.features. Perform
         clustering with 'cluster_timecourse_features()' first.")
  }
  if(clustered) {
    cluster_map <- get_cluster_map(object)
  } else {
    cluster_map <- data.frame(feature = features, cluster = "C1")
  }
  if(!kegg) {
    if(nchar(species) != 2) { stop("wrong \"species\" name.") }
    species_db_pkg <- paste0("org.", species, ".eg.db")
    if(!(species_db_pkg %in% installed.packages())){
      stop(species_db_pkg, "must be installed use the following commands",
           "source(\"https://bioconductor.org/biocLite.R\")
            biocLite(", species_db_pkg, ")")

    }
  }
  feature_df <- suppressMessages(
    data.frame(feature = features, stringsAsFactors = FALSE) %>%
      left_join(cluster_map %>% mutate(feature = as.character(feature))) %>%
      arrange(cluster)
  )
  freq_clust <- feature_df %>%
    group_by(cluster) %>%
    summarise(freq = n())

  res <- vector("list", length(unique(feature_df$cluster)))
  names(res) <- unique(feature_df$cluster)
  for (clst in names(res)) {
    if(length(res) > 1){
      message("Testing enrichment for cluster: ", clst)
    }
    clst_df <- feature_df %>%
      filter(feature_df$cluster == clst)
    if(!kegg){
      clust_res <- goana(
        de = clst_df$feature,
        universe = feature_names(object),
        species = species, ...)
    } else{
      kegg_species <- ifelse(nchar(species) == 3, species, NULL)
      clust_res <- kegga(
        de = clst_df$feature,
        universe = feature_names(object),
        species = species,
        species.KEGG = kegg_species, ...)
    }
    n_DE <- freq_clust %>% filter(cluster == clst)
    n_DE <- fltr_DE * n_DE[["freq"]]
    res[[clst]] <- clust_res %>%
      filter(Ont %in% ontology, DE > n_DE, N <= fltr_N, P.DE <= fltr_P.DE)
  }
  if(length(res) == 1) { res <- res[[1]] }
  return(res)
}



