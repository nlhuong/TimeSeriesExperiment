library(tidyverse)

url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE114762&",
              "format=file&file=GSE114762_raw_counts.csv.gz")
url2 <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE114762&",
               "format=file&file=GSE114762_gene_data.csv.gz")
url3 <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE114762&",
               "format=file&file=GSE114762_sample_info.csv.gz")

raw.data <- read_csv(url) %>% 
  column_to_rownames("X1") %>%
  as.matrix()
feature.data <- read_csv(url2) %>% 
  as.data.frame() %>%
  column_to_rownames("X1")
sample.data <- read_csv(url3) %>% 
  as.data.frame() %>%
  column_to_rownames("X1")

cop1_eset <- ExpressionSet(
  as.matrix(raw.data), 
  phenoData = AnnotatedDataFrame(sample.data),
  featureData = AnnotatedDataFrame(feature.data))

