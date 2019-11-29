## ---- SetVariables

#SeuratObjectList <- SeuratList.Filtered

## ---- FET_IntegrateData

IntegrateDatasets <- function(SeuratObjectList, 
                              ndims = 30,
                              saveTo = "data/FET/allData.RData",
                              REDO = !(file.exists(saveTo)), SAVE = !(file.exists(saveTo))){
  
  if (REDO == TRUE | !(file.exists(saveTo))){
    
    for (i in 1:length(SeuratObjectList)){
      study <- names(SeuratObjectList)[i]
      SeuratObjectList[[i]] <- AddMetaData(SeuratObjectList[[i]], metadata = rep(study, dim(SeuratObjectList[[i]]@meta.data)[1]), col.name = "source")
      SeuratObjectList[[i]] <- AddMetaData(SeuratObjectList[[i]], metadata = Idents(SeuratObjectList[[i]]), col.name = "cell_type")
    }
    
    anchors <- FindIntegrationAnchors(object.list = SeuratObjectList, dims = 1:ndims)
    allData <- IntegrateData(anchorset = anchors, dims = 1:ndims)
    
    DefaultAssay(allData) <- "integrated"
    
    allData <- ScaleData(allData, verbose = FALSE)
    allData <- RunPCA(allData, npcs = ndims, verbose = FALSE)
    allData <- RunUMAP(allData, reduction = "pca", dims = 1:ndims)
    
  } else {
    allData <- readRDS(saveTo)
  }
  
  if (SAVE == TRUE){ saveRDS(allData, saveTo) }
  return(allData)
  
}



