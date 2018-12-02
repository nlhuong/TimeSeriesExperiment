#' @title Pathway enrichment testing.
#'
#' @description This is a wrapper around \code{\link[limma]{goana}} function
#' from \code{limma} package for testing over-representation of gene ontology
#'(GO) terms or KEGG pathways in sets of genes.
#'
#' @param object A \code{TimeSeriesExperiment} object.
#' @param features A vector of ENTREZID for enrichment testing.
#' @param species A character string    specifying the species.
#' See \code{\link[limma]{goana}} for details.
#' @param feature_column the feature column in 'feature.data' slot holding
#' ENTREZIDs, by deafult the 'feature' column.
#' @param universe A vector of genes in the universe. By default
#' all the genes in the 'raw.data' slot.
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
#' \code{\link[limma:goana]{limma::kegga()}}.
#'
#' @return a \code{data.frame} or list of \code{data.frame}s with enrichment
#'results.
#'
#' @importFrom dplyr left_join filter arrange group_by summarise n
#' @importFrom limma goana kegga
#' @importFrom utils installed.packages
#' @importFrom SummarizedExperiment rowData
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
#' head(enrich_res)
#'
pathwayEnrichment <- function(object, features, species, 
                              feature_column = "feature", universe = NULL, 
                              clustered = TRUE, kegg = FALSE,
                              ontology = c("BP", "CC", "MF"), fltr_DE = 0.1, 
                              fltr_N = 500, fltr_P.DE = 0.05, ...)
{
    if (!is(object, "TimeSeriesExperiment")) 
        stop("Input must be a 'TimeSeriesExperiment' object.")
    if (!validObject(object))
        stop("Invalid TimeSeriesExperiment object.")
    
    feature <- cluster <- Ont <- DE <- N <- P.DE <- NULL
    if(is.null(universe)) {
        universe <- rownames(object)
        if (!all(features %in% universe)) {
            stop("Some selected 'features' are not in the 'universe'.")
        }
    }
    if(all(clustered, is.null(clusterMap(object)))) {
        stop("No 'cluster_map' in object@cluster.features. Perform, ",
             "clustering with 'cluster_timecourse_features()' first.")
    }
    feature_df <- as.data.frame(rowData(object))
    feature_df$feature <- feature_df[[feature_column]]
    if(!clustered) feature_df$cluster <- "OneCluster"
  
    feature_df <- feature_df %>%
        filter(feature %in% features) %>%
        arrange(cluster)

    freq_clust <- feature_df %>%
        group_by(cluster) %>%
        summarise(freq = n())

    if(!kegg) {
        if(nchar(species) != 2) { stop("wrong \"species\" name.") }
        species_db_pkg <- paste0("org.", species, ".eg.db")
        if (!requireNamespace(species_db_pkg, quietly = TRUE)) {
          stop("Package ", species_db_pkg, " needed for this function to work.",
               "Please install it using: \n", 
               "if (!requireNamespace(\"BiocManager\", quietly = TRUE))\n",
               "    install.packages(\"BiocManager\")\n",
               "BiocManager::install(\"", species_db_pkg, "\")", call. = FALSE)
        }
    }
    res <- vector("list", length(unique(feature_df$cluster)))
    names(res) <- unique(feature_df$cluster)
    for (clst in names(res)) {
        if(length(res) > 1){
            message("Testing enrichment for cluster: ", clst)
        }
        clst_df <- feature_df %>%
            filter(feature_df$cluster == clst)
        if(nrow(clst_df) == 0) {
            res[[clst]] <- NULL
            next
        }
        if(kegg){
            kegg_species <- ifelse(nchar(species) == 3, species, NULL)
            clust_res <- kegga(
                de = clst_df$feature,
                universe = universe,
                species = species,
                species.KEGG = kegg_species, ...)
        } else{
            clust_res <- goana(
                de = clst_df$feature,
                universe = universe,
                species = species, ...)
        }
        n_DE <- freq_clust %>% filter(cluster == clst)
        n_DE <- fltr_DE * n_DE[["freq"]]
        res[[clst]] <- clust_res %>%
            filter(Ont %in% ontology, DE > n_DE, N <= fltr_N, 
                   P.DE <= fltr_P.DE)
    }
    if(length(res) == 1) { res <- res[[1]] }
    return(res)
}
