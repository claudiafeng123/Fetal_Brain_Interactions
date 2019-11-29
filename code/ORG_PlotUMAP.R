#organoid.id <- organoid.list[1]

## ---- ORG_PlotUMAP-CODE

PlotUMAP <- function(organoid.id, .plotBy = c("predicted.id", "sub_batch"),
                     savePDF = TRUE, 
                     savePDF.path, savePDF.width = 8, savePDF.height = 6){
  
  expr.data <- readRDS(paste0(savedData.folder, "seurat_objects/", organoid.id, ".RData"))
  
  if (savePDF == TRUE){pdf(savePDF.path, width = savePDF.width, height = savePDF.height)}
  for (feature in .plotBy){
    print(DimPlot(expr.data, reduction = 'umap', group.by = feature))
  }
  while(!(names(dev.cur()) %in% c('RStudioGD', 'null device'))){dev.off()}
  
}



