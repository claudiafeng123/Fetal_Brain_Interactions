## ---- Polioudakis2019_SetVariables

#study <- "Polioudakis2019"

#folderPath <- paste0("data/FET/", study)
#data.URL <- "http://solo.bmap.ucla.edu/shiny/webapp/session/4f81dc6713be6627658f457bbb100974/download/downloadData?w="
#data.downloadSuffix <- ".zip"
#rawCounts.matrixPath <- paste0(folderPath, "/raw_counts_mat.rdata")
#metadata.matrixPath <- paste0(folderPath, "/cell_metadata.csv")

#saveTo = "data/FET/Polioudakis2019/allData.RData"
#saveTo = paste0(folderPath, "/allData.RData")

#ndims = 40
#min.cells = 3
#min.features = 200

## ---- Polioudakis2019_LoadData

loadPolioudakis2019 <- function(study = "Polioudakis2019", folderPath = paste0("data/FET/", study),
                               data.URL = NULL, data.downloadSuffix = ".zip", 
                               rawCounts.matrixPath, metadata.matrixPath,
                               saveTo = paste0(folderPath, "/allData.RData"), 
                               DOWNLOAD = !(file.exists(folderPath)), SAVE = (!(file.exists(saveTo)) || DOWNLOAD),
                               min.cells = 3, min.features = 200){
  
  ## study: c('Nowakowski2017', 'Polioudakis2019')
  ## folderPath: folder where both matrices will be saved to (by default, data/Nowakowski2017/)
  ## data.URL: URL for the data (it's in one zip folder)
  ## rawCounts.downloadSuffix: .gz, .zip (data is often saved in different formats)
  ## rawCounts.matrixPath: path to saved rawCount matrix
  ## metadata.URL: URL for metadata
  ## metadata.matrixPath: path to saved metadata matrix
  ## saveTo: save metadata and rawcounts as R object
  
  if ((DOWNLOAD == TRUE)){
    
    #raw counts matrix
    download.file(data.URL, destfile = paste0(folderPath, data.downloadSuffix))
    unzip(paste0(folderPath, data.downloadSuffix), exdir = folderPath)
    
  }
  
  if (!(file.exists(saveTo)) | (DOWNLOAD == TRUE)){
    
    #load metadata
    metadata <- fread(metadata.matrixPath)
    
    #load counts matrix
    load(rawCounts.matrixPath)
    rawCountsMatrix <- as.matrix(raw_counts_mat)
    
    #metadata missing data, get rid of extra entries from counts matrix
    rawCountsMatrix <- rawCountsMatrix[, which(colnames(rawCountsMatrix) %in% metadata$Cell)]
    metadata <- metadata[which(metadata$Cell %in% colnames(rawCountsMatrix)),]
    
    #adjust cell names
    metadata$Cell_ID <- metadata$Cell
    metadata$Cell <- paste(metadata$Cell_ID, metadata$Cluster, metadata$Donor, sep = "_")
    colnames(rawCountsMatrix) <- metadata$Cell[match(colnames(rawCountsMatrix), metadata$Cell_ID)]#checked that this doesn't mess  up order of celsl
    row.names(metadata) <- metadata$Cell
    #metadata <- metadata[which(metadata$Cell %in% colnames(rawCountsMatrix)),]
    
    ## QC + Normalize
    data <- CreateSeuratObject(counts = rawCountsMatrix, project = study, 
                               min.cells = min.cells, min.features = min.features, 
                               meta.data = metadata, names.delim = "_", names.field = 2)
    data[["MT_percent"]] <- PercentageFeatureSet(data, pattern = "^MT-")
    
    data <- SplitObject(data, split.by = "Donor")
    
    for (i in 1:length(data)){
      data[[i]] <- NormalizeData(data[[i]], verbose = FALSE)
      data[[i]] <- FindVariableFeatures(data[[i]], selection.method = "vst", nfeatures = dim(data[[i]][["RNA"]]@counts)[1])
    }
    
    data <- FindIntegrationAnchors(data, anchor.features = dim(data[[i]][["RNA"]]@counts)[1])
    data <- IntegrateData(data, dims = 1:ndims)
    
    DefaultAssay(data) <- "integrated"
    
    Polioudakis2019 <- data
    
    if (SAVE == TRUE){saveRDS(Polioudakis2019, file = saveTo)}
    
  } else {
    Polioudakis2019 <- readRDS(saveTo)
  }
  
  return(Polioudakis2019)
  
}
