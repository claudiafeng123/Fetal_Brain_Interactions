

## ---- Nowakowski2017_SetVariables

#study <- "Nowakowski2017"

#folderPath <- paste0("data/FET/", study, "/")
#rawCounts.URL <-"https://cells.ucsc.edu/cortex-dev/exprMatrix.tsv.gz"
#rawCounts.downloadSuffix <- ".gz"
#rawCounts.matrixPath <- paste0(folderPath, study, "_rawCounts.tsv")
#metadata.URL <- "https://cells.ucsc.edu/cortex-dev/meta.tsv"
#metadata.matrixPath <- paste0(folderPath, "/", study, "_metadata.tsv")
#saveTo = paste0(folderPath, "allData.RData")

#ndims = 40
#min.cells = 3
#min.features = 200

## ---- Nowakowski2017_LoadData

loadNowakowski2017 <- function(study = "Nowakowski2017", folderPath = paste0("data/FET/", study),
                               rawCounts.URL = NULL, rawCounts.downloadSuffix = NULL, 
                               metadata.URL = NULL, 
                               rawCounts.matrixPath, metadata.matrixPath,
                               min.cells = 3, min.features = 200,
                               saveTo = paste0(folderPath, "allData.RData"), 
                               DOWNLOAD = !(file.exists(folderPath)), SAVE = (!(file.exists(saveTo)) || DOWNLOAD)){
  
  ## study: c('Nowakowski2017', 'Polioudakis2019')
  ## folderPath: folder where both matrices will be saved to (by default, data/Nowakowski2017/)
  ## rawCounts.URL: URL for the raw counts
  ## rawCounts.downloadSuffix: .gz, .zip (data is often saved in different formats)
  ## rawCounts.matrixPath: path to saved rawCount matrix
  ## metadata.URL: URL for metadata
  ## metadata.matrixPath: path to saved metadata matrix
  ## saveTo: save metadata and rawcounts as R object
  
  if ((DOWNLOAD == TRUE)){
    
    #raw counts matrix
    download.file(rawCounts.URL, destfile = paste0(rawCounts.matrixPath, rawCounts.downloadSuffix))
    gunzip(paste0(rawCounts.matrixPath, rawCounts.downloadSuffix))
    #metadata
    download.file(metadata.URL, destfile = metadata.matrixPath)
    
  }
  
  if (!(file.exists(saveTo)) | (DOWNLOAD == TRUE)){
    
    rawCountsMatrix <- fread(rawCounts.matrixPath)
    geneList <- as.vector(t(rawCountsMatrix[,1]))
    rawCountsMatrix <- as.matrix(rawCountsMatrix[, -1])
    row.names(rawCountsMatrix) <- geneList
    
    metadata <- as.data.frame(fread(metadata.matrixPath))
    
    #redo id's because rna magnet was unhappy
    newIDs <- gsub(unlist(lapply(lapply(strsplit(metadata$Cell, "[.]"), unique), paste, collapse = "-")), pattern = "_", replacement = "-")
    newIDs <- paste(newIDs, metadata$WGCNAcluster, "2017", sep = "_")
    row.names(metadata) <- newIDs
    colnames(rawCountsMatrix) <- newIDs[match(colnames(rawCountsMatrix), metadata$Cell)]
    
    ##
    #log transform
    tpmMatrix = log(rawCountsMatrix + 1)
    
    ##create seurat object
    Nowakowski2017 <- CreateSeuratObject(counts = tpmMatrix, project = study, 
                                       min.cells = min.cells, min.features = min.features,
                                       meta.data = metadata, names.field = 2)
    
    if (SAVE == TRUE){saveRDS(Nowakowski2017, file = saveTo)}
    
  } else {
    Nowakowski2017 <- readRDS(saveTo)
  }
  
  return(Nowakowski2017)
  
}
