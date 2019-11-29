## ---- SetVariables

#mappedObject <- readRDS("data/ORG/seurat_objects/Velasco2019_PGP1-3mo-Batch2.RData")


#mappedObjects.paths <- paste0("data/ORG/seurat_objects/", list.files("data/ORG/seurat_objects/"))
#saveMetadata.path = paste0(savedData.folder, "organoid-metadata.RData")

#.plotBy = c("cell_line", "time_point")

## ---- ORG_EvaluateMappedLabels

EvaluateMappedLabels <- function(mappedObjects.paths,
                                  savePDF = TRUE, savePDF.path = paste0(savePDF.folder, "LAB_prediction-score-histograms.pdf"), savePDF.width = 8, savePDF.height = 6,
                                 .plotBy = c("cell_line", "time_point", 'predicted.id'),
                                 saveMetadata.path = paste0(savedData.folder, "organoid-metadata.RData"), saveMetadata = !(file.exists(saveMetadata.path))){
  
  #load metadata
  if (saveMetadata == TRUE | !(file.exists(saveMetadata.path))){
    metadata <- foreach (x = mappedObjects.paths, .combine = 'rbind') %do% {
      mappedObject <- readRDS(x)
      return(mappedObject@meta.data)
    }
    saveRDS(metadata, saveMetadata.path)
  } else {
    metadata <- readRDS(saveMetadata.path)
  }
  
  if (savePDF == TRUE){pdf(savePDF.path, width = savePDF.width, height = savePDF.height)}
  hist(metadata$prediction.score.max, main= "Distribution of Prediction Scores", xlab = "Prediction Scores", ylab = "Frequency", breaks = 40)
  for (column.name in .plotBy){
    eval(parse(text = paste0("print(ggplot(metadata, aes(x=prediction.score.max, color = ", column.name, ")) + ",
                             "xlab('Prediction Score') + ylab('Density') + ",
      "geom_density(alpha = 0.3) + theme_classic())")))
  }
  while(!(names(dev.cur()) %in% c('null device', 'RStudioGD'))){dev.off()}
  
  return(metadata)
  
}






