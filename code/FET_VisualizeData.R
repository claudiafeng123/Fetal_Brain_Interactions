

## ---- SetVar

#.plotBy = c("cell_type", "source")
#reduction = "umap"
#saveTo.folder = fig.folder

## ---- FET_VisualizeData

VisualizeData <- function(allData, .plotBy = c("cell_type", "source"),
                          reduction = "umap",
                          ndims = 30, 
                          saveFigPDF = TRUE, 
                          fig.width = 8, fig.height = 6,
                          saveTo.folder = fig.folder){
  
  for (feature in .plotBy){
    if (saveFigPDF == TRUE){pdf(paste0(fig.folder, "UMAP_", gsub(feature, pattern = "_", replacement = "-"), ".pdf"), width = fig.width, height = fig.height)}
    print(DimPlot(allData, reduction =  reduction, group.by = feature))
    while(!(names(dev.cur()) %in% c('RStudioGD', 'null device'))){dev.off()}
  }
  
}


