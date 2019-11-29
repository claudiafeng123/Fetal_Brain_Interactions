
## ---- IXN_GeneCellStats-CODE
CheckCPDBQuality <- function(fetusData.path = "data/FET/Polioudakis2019/allData.RData", 
                             organoids = sample.ids, organoidData.folder = "data/ORG/seurat_objects/", 
                             CPDB.resultsFolder = paste0(results.folder, "cpdb_input/"),
                             writeTo = NULL){
  
  rtn <- data.frame(nGenes = numeric(0),
                    nCells = integer(0),
                    nLRPairs = integer(0))
  
  for (i in 1:length(organoids)){
    sample <- organoids[i]
    if (sample == "fetus"){
      data <- readRDS(fetusData.path)
    } else {
      data <- readRDS(paste0(organoidData.folder, sample, ".RData"))
    }
    nGenes <- mean(data$nFeature_RNA)
    nCells <- dim(data@meta.data)[1]
    nLRPairs <- dim(fread(paste0(CPDB.resultsFolder, sample, "/means.txt")))[1]
    rtn[i, ] <- c(nGenes, nCells, nLRPairs)
  }
  
  row.names(rtn) <- sample.ids
  
  if (!(is.null(writeTo))){fwrite(rtn, writeTo, sep = "\t", row.names = TRUE)}
  return(rtn)
}

