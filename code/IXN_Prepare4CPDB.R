## ---- SetVariables

#fetusData.path = "data/FET/Polioudakis2019/allData.RData"
#organoidData.folder = "data/ORG/seurat_objects/"

#writeTo = CPDB.folderPath
#meta.fileName <- "_meta.txt"
#counts.fileName <- "_counts.txt"
#organoids = sample.ids

## ---- IXN_RunCPDB
#for each cohort of data, print a counts matrix and metadata file

Prepare4CPDB <- function(fetusData.path = "data/FET/Polioudakis2019/allData.RData", 
                         organoidData.folder = "data/ORG/seurat_objects/", organoids = sample.ids,
                         meta.fileName = "_meta.txt", counts.fileName = "_counts.txt",
                         writeTo = paste0(results.folder, "cpdb_input/"),
                         REDO = FALSE ) {
  
  #fetus
  
  if ((!(file.exists(paste0(writeTo, "fetus", counts.fileName)))) | (REDO == TRUE)){
    
    fetus.data <- readRDS(fetusData.path)
    fetus.cellCounts <- GetAssayData(fetus.data, slot = "data")
    cellNames <- colnames(fetus.cellCounts)
    
    #need to convert genes to ENSEMBL IDs
    geneList <- row.names(fetus.cellCounts)
    geneList <- genes(EnsDb.Hsapiens.v86, filter = list(GeneNameFilter(geneList),
                                                        GeneIdFilter("ENSG", "startsWith")), return.type = "data.frame", columns = c("gene_id"))
    
    ## change row and column names
    fetus.cellCounts <- fetus.cellCounts[which(row.names(fetus.cellCounts) %in% geneList$gene_name),]
    fetus.cellCounts <- data.frame(Gene = row.names(fetus.cellCounts), fetus.cellCounts)
    fetus.cellCounts$Gene <- geneList$gene_id[match(fetus.cellCounts$Gene, geneList$gene_name)]
    colnames(fetus.cellCounts) <- c('Gene', cellNames)
    
    meta <- data.frame(row.names(fetus.data@meta.data), fetus.data@meta.data[, "orig.ident"])
    names(meta) <- c('Cell', 'cell_type')
    
    fwrite(meta, paste0(writeTo, "fetus", meta.fileName), sep = "\t")
    fwrite(fetus.cellCounts, paste0(writeTo, "fetus", counts.fileName), sep = "\t")
    
  }
  
  #organoids
  if (is.null(organoids) == TRUE){ 
    organoidFilePaths <- paste0(organoidData.folder, list.files(organoidData.folder))
    } else {organoidFilePaths <- paste0(organoidData.folder, organoids, ".RData")}
  
  for (organoidFile in organoidFilePaths){
    
    organoid.id <- gsub(unlist(lapply(strsplit(organoidFile, "/"), "[[", 4)), pattern = ".RData", replacement = "")
    
    if ((!(file.exists(paste0(writeTo, organoid.id, counts.fileName)))) | (REDO == TRUE)){
      
      organoid.data <- readRDS(organoidFile)
      organoid.cellCounts <- GetAssayData(organoid.data, slot = "data")
      cellNames <- colnames(organoid.cellCounts)
      
      #convert genes to ENSEMBL IDs
      geneList <- row.names(organoid.cellCounts)
      geneList <- genes(EnsDb.Hsapiens.v86, filter = list(GeneNameFilter(geneList),
                                                          GeneIdFilter("ENSG", "startsWith")), return.type = "data.frame", columns = c("gene_id"))
      
      organoid.cellCounts <- organoid.cellCounts[which(row.names(organoid.cellCounts) %in% geneList$gene_name),]
      organoid.cellCounts <- data.frame(Gene = row.names(organoid.cellCounts), organoid.cellCounts)
      organoid.cellCounts$Gene <- geneList$gene_id[match(organoid.cellCounts$Gene, geneList$gene_name)]
      colnames(organoid.cellCounts) <- c('Gene', cellNames)
      
      meta <- data.frame(row.names(organoid.data@meta.data), organoid.data@meta.data[, "predicted.id"])
      names(meta) <- c('Cell', 'cell_type')
      
      fwrite(meta, paste0(writeTo, organoid.id, meta.fileName), sep = "\t")
      fwrite(organoid.cellCounts, paste0(writeTo, organoid.id, counts.fileName), sep = "\t")
    }
    
  }
  
}



