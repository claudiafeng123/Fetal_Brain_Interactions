#organoid.id <- organoid.list[1]
#organoid.exprMatrix.path <- organoid.filePaths[1]
#organoid.id <- organoid.list[1]

## ---- ORG_LabelCellTypes

LabelCellTypes <- function(organoid.id, organoid.exprMatrix.path,
                           fetus.data,
                           ndims, batch.sep = "_",
                           REDO = FALSE, saveTo = NULL){
  
  if ((REDO == TRUE) | !(file.exists(saveTo))){
    
    study <- unlist(lapply(strsplit(organoid.id, split = "_"), "[[", 1))
    cell.line <- unlist(lapply(strsplit(unlist(lapply(strsplit(organoid.id, split = "_"), "[[", 2)), split = "-"), "[[", 1))
    time.point <- unlist(lapply(strsplit(unlist(lapply(strsplit(organoid.id, split = "_"), "[[", 2)), split = "-"), "[[", 2))
    batch <- unlist(lapply(strsplit(unlist(lapply(strsplit(organoid.id, split = "_"), "[[", 2)), split = "-"), "[[", 3))
    
    organoid.exprMatrix <- fread(organoid.exprMatrix.path)
    
    #convert 
    geneList <- as.vector(t(organoid.exprMatrix[,1]))
    organoid.exprMatrix <-  as.matrix(organoid.exprMatrix[, -1])
    row.names(organoid.exprMatrix) <- geneList
    cellIDs <- colnames(organoid.exprMatrix)
    
    organoid.data <- CreateSeuratObject(counts = organoid.exprMatrix, project = study)
    
    #add metadata
    ncells <- dim(organoid.exprMatrix)[2]
    organoid.data <- AddMetaData(organoid.data, metadata = rep(study, ncells), col.name = "source")
    organoid.data <- AddMetaData(organoid.data, metadata = rep(cell.line, ncells), col.name = "cell_line")
    organoid.data <- AddMetaData(organoid.data, metadata = rep(time.point, ncells), col.name = "time_point")
    organoid.data <- AddMetaData(organoid.data, metadata = rep(batch, ncells), col.name = "batch")
    organoid.data <- AddMetaData(organoid.data, metadata = rep(organoid.id, ncells), col.name = "ID")
    sub.batch <- unlist(lapply(strsplit(colnames(organoid.exprMatrix), batch.sep), "[[", 2))
    organoid.data <- AddMetaData(organoid.data, metadata = sub.batch, col.name = "sub_batch")
    
    organoid.data <- NormalizeData(organoid.data, selection.method = "vst")
    
    anchors <- FindTransferAnchors(reference = fetus.data, query = organoid.data, dims = 1:ndims)
    integrated.data <- TransferData(anchorset = anchors, refdata = Idents(fetus.data), dims = 1:ndims)
    organoid.data <- AddMetaData(organoid.data, metadata = integrated.data)
    
    organoid.data <- ScaleData(organoid.data)
    organoid.data <- FindVariableFeatures(organoid.data)
    organoid.data <- RunPCA(organoid.data, npcs = ndims)
    organoid.data <- RunUMAP(organoid.data, reduction = "pca", dims = 1:ndims)
    
    if (is.null(saveTo) == FALSE){ saveRDS(organoid.data, saveTo)}
    
  } else {
    organoid.data <- readRDS(saveTo)
  }
  
  return(organoid.data)
  
}

