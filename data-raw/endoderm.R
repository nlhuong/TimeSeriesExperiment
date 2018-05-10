# This script specified how "vistimeseq" data object: "endoderm_small" was built
# using the data from a study by Blake et al. (2017)
# (https://www.biorxiv.org/content/early/2017/05/09/135442) with accession
# number: GSE98411, which we downloaded from Gene Ontology Onmibus (GEO)
# repository. "endoderm_samll" is used as dataset for running examples.

library(dplyr)
library(tidyr)
library(tibble)
library(biomaRt)

url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE98411&format=file&file=GSE98411%5FRNA%5Fcounts%2Etxt%2Egz"
tmp <- tempfile()
download.file(url,tmp)

cnts <- read.table(
  gzfile(tmp),
  header=TRUE,
  stringsAsFactors=FALSE
)
dim(cnts) # [1] 30030    64

smpDF <- data.frame(
  sample = colnames(cnts)
  )  %>%
  mutate(
    time = gsub("(.*)\\_", "", sample),
    time = substr(time, 2, 2),
    replicate = gsub("\\_(.*)", "", sample),
    group = substr(replicate, 1, 1)
  ) %>%
  column_to_rownames("sample")
head(smpDF)

# Filter out the genes that are very low expression levels.
# Keep genes with minimum mean count of 'minCntsMeans' in at least in one group
minCntsMeans <- 10
groupCntsMeans <- data.frame(row.names = rownames(cnts))
for(g in unique(smpDF$group)) {
  groupCntsMeans[, g] = apply(cnts[ , which(smpDF$group == g)], 1, mean)
}
groupCntsMax <- apply(as.vector(groupCntsMeans), 1, max)
isexpr <- (groupCntsMax > minCntsMeans)
cnts <- cnts[isexpr, ]
dim(cnts) #[1] 14694    64

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(cnts)
gene_df <- getBM(
  filters= "ensembl_gene_id",
  attributes= c("ensembl_gene_id", "entrezgene", "description"),
  values=genes, mart= mart
)
gene_df <- gene_df[!duplicated(gene_df$entrezgene), ]
gene_df <- gene_df[!is.na(gene_df$entrezgene), ]
dim(gene_df)

cnts <- cnts[gene_df$ensembl_gene_id, ]
rownames(cnts) <- gene_df$entrezgene

endoderm <- vistimeseq(
  project = "[Subset] Endoderm differntiation in human and chimpanzee",
  raw.data = cnts,
  sample.data = smpDF,
  time_column = "time",
  replicate_column = "replicate",
  group_column = "group"
)
#save(endoderm, file = "data/endoderm.rda")

# Pick the top 250 most variable genes and save as smaller size data
cpm <- apply(cnts, 2, function(x) 1e6*x/sum(x))
top250 <- names(sort(apply(log(cpm + 1), 1, sd), decreasing = TRUE))[1:250]
endoderm_small <- filter_features(endoderm, top250)
endoderm_small
save(list = "endoderm_small", file = "data/endoderm_small.rda")
