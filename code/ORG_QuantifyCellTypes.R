#fetus.metadata = fetus.data@meta.data
#organoid.metadata <- readRDS(paste0(savedData.folder, "organoid-metadata.RData"))

## ---- ORG_QuantifyCellTypes-CODE

QuantifyCellTypes <- function(organoid.metadata, 
                              fetus.metadata, cellTypeList,
                              saveRDS.path = paste0("data/ORG/", "cell-counts.RData"),
                              savePDF = TRUE, 
                              savePDF.path = paste0(savePDF.folder, "CC_fraction-variance.pdf"), savePDF.width = 8, savePDF.height = 6){
  
  fetal.cellBreakdown <- table(fetus.metadata$orig.ident)/dim(fetus.metadata)[1]
  
  organoid.metadata$toSplit <- paste(organoid.metadata$ID, organoid.metadata$sub_batch, sep = ":")
  organoid.metadata <- split(organoid.metadata, f = organoid.metadata$toSplit)
  organoid.cellBreakdown <- t(mapply(organoid.metadata, FUN = function(x){
    return((table(c(cellTypeList, x$predicted.id))-1))
  }))
  
  organoid.cellFrac <- apply(organoid.cellBreakdown, 2, "/", apply(organoid.cellBreakdown, 1, sum))
  organoid.cellFrac <- apply(organoid.cellFrac, 1, "/", fetal.cellBreakdown)
  
  #barplot to see the variance in proportion over different samples
  claudia <- data.frame(log(organoid.cellFrac))
  claudia$cell_type <- row.names(claudia)
  claudia <- data.table::melt(claudia, id.vars = "cell_type")
  claudia <- claudia[which(is.finite(claudia$value)), ]
  if (savePDF == TRUE){pdf(savePDF.path, width = savePDF.width, height = savePDF.height)}
  par(mar=c(4.1, 5.1, 4.1, 8.1), mgp = c(4, 1, 0), xpd=TRUE)
  boxplot(value~cell_type, data = claudia,
          main = 'Log-Change in Organoid Cell Composition', xlab = '', ylab = 'Cell Type',
          col = ColorChooser(cellTypeList),
          horizontal = TRUE, las = 2)
  legend("topright", inset=c(-0.25,0), legend = cellTypeList, fill = ColorChooser(cellTypeList))
  heatmap(cor(organoid.cellBreakdown), main = "Cell Count Correlation")
  heatmap(cor(apply(organoid.cellBreakdown, 2, "/", apply(organoid.cellBreakdown, 1, sum))),
          main = "Cell Proportion Correlation")
  while(!(names(dev.cur()) %in% c('RStudioGD', 'null device'))){dev.off()}
  
  if (!(is.null(saveRDS.path))){saveRDS(organoid.cellBreakdown, saveRDS.path)}
  
  return(organoid.cellBreakdown)
  
}

