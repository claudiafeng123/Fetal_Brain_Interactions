## ---- SetVariables

#SeuratObjectList <- readRDS("data/FET/allData_renamed.RData")
#predictionScoreRatio.cutoff <- 0.8

#saveTo <- "allData_FilterCells.RData"

## ---- FET_FilterCells

FilterCells <- function(SeuratObjectList, main.study,
                        PLOT = TRUE, REDO = FALSE, 
                        save.fileName = "allData_FilterCells.RData", SAVE = (REDO | !(file.exists(save.fileName))),
                        savePDF = PLOT, savePDF.folder = fig.folder, savePDF.width = 13, savePDF.height = 5){
  
  rtn <- list()
  
  for (i in 1:length(SeuratObjectList)){
    
    study <- names(SeuratObjectList)[i]
    data <- SeuratObjectList[[i]]
    
    ##FILTER
    
    if (PLOT == TRUE){
      if (study == main.study){category = "orig.ident"} else {category = "predicted.id"}
      if (savePDF == TRUE){pdf(paste0(savePDF.folder, "QC_", study, "_FilterCells.pdf"), width = savePDF.width, height = savePDF.height)}
      if ("MT_percent" %in% colnames(data@meta.data)){
        print(VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "MT_percent"), group.by = category,
                      ncol = 3, pt.size = 0.02, cols = ColorChooser()))
      } else {print(VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA"), group.by = category,
                            ncol = 3, pt.size = 0.02, cols = ColorChooser()))}
      
      while(!names(dev.cur()) %in% c('RStudioGD', 'null device')){dev.off()}
    }
    
    eval(parse(text = paste0("rtn[['", study, "']] <- data")))
    
  }
  
  if (SAVE == TRUE){ saveRDS(rtn, paste0(saveRDS.folder, "/", save.fileName)) }
  return(rtn)
  
}

#data <- subset(data, subset = (nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5))

#violin plot
#nFeature
#nCount
#pct_MT
#get rid of badly mapped cells
#return objects