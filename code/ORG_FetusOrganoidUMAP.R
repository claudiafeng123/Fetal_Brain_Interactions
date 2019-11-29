## ---- SetVariables

#organoidData.paths <- paste0("data/ORG/seurat_objects/", list.files("data/ORG/seurat_objects/"))
#saveTo = paste0(savedData.folder, "allData.RDS")
#.plotBy = c('cell_type', 'source')

#savePDF.path = paste0(savePDF.folder, "UMAP_fetus-organoid.pdf")
#savePDF.width = 9
#savePDF.height = 6

## ---- ORG_FetusOrganoidUMAP-code

FetusOrganoidUMAP <- function(organoidData.paths = NULL,
                              fetus.data = NULL, .plotBy = c('cell_type', 'source'),
                              ndims = 30,
                              saveTo = paste0(savedData.folder, "allData.RDS"), REDO = FALSE, SAVE = TRUE,
                              savePDF.path = paste0(savePDF.folder, "UMAP_fetus-organoid.pdf"), savePDF.width = 10, savePDF.height = 6){
  
  if ((REDO == TRUE) | !(file.exists(saveTo))){
    organoid.data <- mapply(organoidData.paths, FUN = readRDS)
    names(organoid.data) <- organoid.list
    
    allData <- organoid.data
    allData[["Polioudakis2019"]] <- fetus.data
    
    #integrate data
    allData <- FindIntegrationAnchors(allData, dims = 1:ndims)
    allData <- IntegrateData(allData, dims = 1:ndims)
    DefaultAssay(allData) <- "integrated"
    
    cellType <- allData@meta.data$predicted.id
    cellType[which(is.na(allData@meta.data$predicted.id))] <- allData@meta.data$Cluster[which(is.na(allData@meta.data$predicted.id))]
    dataSource <- allData@meta.data$ID
    dataSource[which(is.na(allData@meta.data$ID))] <- "Polioudakis2019"
    allData <- AddMetaData(allData, cellType, "cell_type")
    allData <- AddMetaData(allData, dataSource, "source")
    
    #for visualizing data later
    allData <- ScaleData(allData, verbose = FALSE)
    allData <- RunPCA(allData, npcs = ndims, verbose = FALSE)
    allData <- RunUMAP(allData, reduction = "pca", dims = 1:ndims)
    
    if (SAVE == TRUE){saveRDS(allData, saveTo)}
    
  } else {
    allData <- readRDS(saveTo)
  }
  
  if (!is.null(savePDF.path)){pdf(savePDF.path, width = savePDF.width, height = savePDF.height)}
  for (feature in .plotBy){
    print(DimPlot(allData, reduction = 'umap', group.by = feature))
  }
  while(!(names(dev.cur()) %in% c('RStudioGD', 'null device'))){ dev.off() }
  
  return(allData)
    
}

