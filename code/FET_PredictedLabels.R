## ---- SetVariables

#main.study <- "Polioudakis2019"
#save.fileName = "allData_integrated.RDS"
#secondBest.cutoff = 0.8

#savePDF.width = 8; savePDF.height = 6

#cellTypeList <- unique(Polioudakis2019$Cluster)

## ---- FET_PredictLabels

PredictLabels <- function(studies, main.study,
                          secondBest.cutoff = 0.8,
                          ndims = 30,
                          REDO = FALSE, save.fileName = "allData_PredictedLabels.RData", SAVE = TRUE,
                          savePDF = TRUE, savePDF.width = 8, savePDF.height = 6){
  
  cellTypeList <- as.character(sort(unique(Idents(get(main.study)))))
  rtn <- list()
  eval(parse(text = paste0("rtn[['", main.study, "']] <- ", main.study)))
  
  for (study in studies[-(which(studies == main.study))]){
    
    if ((REDO == TRUE) || (!file.exists(paste0(saveRDS.folder, study, "/", save.fileName)))){
      eval(parse(text = paste0("anchors <- FindTransferAnchors(reference = ", main.study, ", query = ", study, ", dims = 1:", ndims, ")")))
      eval(parse(text = paste0("predictions <- TransferData(anchorset = anchors, refdata = Idents(", main.study, "), dims = 1:", ndims, ")")))
      eval(parse(text = paste0(study, "<- AddMetaData(object = ", study, ", metadata = predictions)")))
      eval(parse(text = paste0("Idents(", study, ") <- 'predicted.id'")))
      
      eval(parse(text = paste0(
        "predictions <- ", study, "@meta.data[, grep(colnames(", study, "@meta.data), pattern = 'predict')]"
      )))
      
      #calculates the ratio of the best match with the second best match
      second.best <- apply(predictions[, -1], 1, FUN = function(vec){
        return(sort(vec, decreasing = TRUE)[3])
      })
      second.best.frac <- second.best/predictions$prediction.score.max
      
      #poor.maps.frac is the fraction of cells of each type whose prediction score if more than some cutoff
      poor.maps.frac <- t(mapply(seq(0, by = 0.01, 1), FUN = function(predictedScore.cutoff){
        (table(c(predictions$predicted.id[which(predictions$prediction.score.max < predictedScore.cutoff)], 
                 cellTypeList))-1)/table(predictions$predicted.id)
      })); row.names(poor.maps.frac) <- seq(0, by = 0.01, 1)
      
      eval(parse(text = paste0(study, " <- AddMetaData(object = ", study, ", metadata = second.best.frac, col.name = 'prediction.ratio')")))
      if (SAVE == TRUE){ saveRDS(get(study), paste0(saveRDS.folder, study, "/", save.fileName)) }
      
    } else { 
      eval(parse(text = paste0(study, " <- readRDS('", paste0(saveRDS.folder, study, "/", save.fileName), "')"))) 
      eval(parse(text = paste0(
        "predictions <- ", study, "@meta.data[, grep(colnames(", study, "@meta.data), pattern = 'predict')]"
      )))
      poor.maps.frac <- t(mapply(seq(0, by = 0.01, 1), FUN = function(predictedScore.cutoff){
        (table(c(predictions$predicted.id[which(predictions$prediction.score.max < predictedScore.cutoff)], 
                 cellTypeList))-1)/table(predictions$predicted.id)
      })); row.names(poor.maps.frac) <- seq(0, by = 0.01, 1)
      eval(parse(text = paste0("second.best.frac <- ", study, "$prediction.ratio")))
    }
    
    ##PLOT stuff
    if (savePDF == TRUE){pdf(paste0(fig.folder, "QC_", study, "_PredictLabels.pdf"), width = savePDF.width, height = savePDF.height)}
    
    #poor.maps.frac
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
    plot(NULL, xlim = c(0, 1), ylim = c(0,1), 
         main = "Fraction of Cells Mapped", xlab = "Prediction Score", ylab = "Proportion of Cells")
    mapply(ind = 1:dim(poor.maps.frac)[2], FUN = function(ind){
      lines(x = seq(0, by = 0.01, 1), y = poor.maps.frac[, ind], col = ColorChooser(colnames(poor.maps.frac)[ind]), lwd = 2)
    })
    legend("topright", inset=c(-0.25,0), legend = cellTypeList, fill = ColorChooser(cellTypeList))
    
    #second.best
    par(mar=c(5.1, 4.1, 4.1, 2.1), xpd = FALSE)
    hist(second.best.frac, breaks = 40,
         main = "Prediction Score Fitness", xlab = "Second Best/Best Prediction Score")
    lines(x = rep(secondBest.cutoff, 2), y = c(0, 10000), col = "red", lty = 2, lwd = 3)
    
    while(!names(dev.cur()) %in% c('RStudioGD', 'null device')){dev.off()}
    
    eval(parse(text = paste0("rtn[['", study, "']] <- ", study)))
    
  }
  
  if (SAVE == TRUE){ saveRDS(rtn, paste0(saveRDS.folder, "/", save.fileName)) }
  return(rtn)
  
}




