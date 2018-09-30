raw_data_dir <- "~/Downloads/Cop1_data"
raw_data_file <- file.path(raw_data_dir, "GSE114762.tar")
data_url <- paste0("http://www.ncbi.nlm.nih.gov/geo/download/", 
                   "?acc=GSE114762&format=file&token=wfonuieepfwvjit")

dir.create(path = raw_data_dir)
download.file(url = data_url, destfile = raw_data_file)
untar(tarfile=raw_data_file, exdir=raw_data_dir)

smp_url <- paste0("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
                  "GSE114762&targ=gsm&form=text&view=brief&token=wfonuieepfwvjit")
smp_data_file <- file.path(raw_data_dir, "GSE114762_sample_data.txt")
download.file(url=smp_url, destfile=smp_data_file)

sampleDat <- getGEO(filename=smp_data_file)
# Just select the sample titles as these contain the meta data we need
sampleDat <- Meta(sampleDat)

#ID's
gsm <- sampleDat$geo_accession

# Get the sample meta data
meta <- sampleDat$characteristics_ch1
