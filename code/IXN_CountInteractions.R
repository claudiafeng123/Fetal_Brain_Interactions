#threshold <- 0.05


## ---- IXN_CountInteractions-CODE


CountInteractions <- function(organoid.id,
                              threshold = 0.05,
                              dataFolder = CountInteractions.dataFolder,
                              writeTo = NULL,
                              savePDF.path = paste0(figure.folder, "CPDB_ixn-counts_", organoid.id, ".pdf"), savePDF.width = 7, savePDF.height = 6){
  
  cpdb.result <- fread(paste0(CountInteractions.dataFolder, organoid.id, ".csv"))
  cpdb.result$significant <- (cpdb.result$p_value < threshold)
  cpdb.result <- cpdb.result[significant == TRUE,]
  
  ixns <- matrix(0, nrow = length(cellTypes), ncol = length(cellTypes))
  row.names(ixns) <- colnames(ixns) <- cellTypes
  
  for (a in 1:length(cellTypes)){
    cell.a <- cellTypes[a]
    for (b in 1:length(cellTypes)){
      cell.b <- cellTypes[b]
      ixns[a,b] <- length(which(paste(cpdb.result$cell_type_a, cpdb.result$cell_type_b, sep = ":") == paste(cell.a, cell.b, sep = ":")))
    }
  }
  
  pdf(savePDF.path, width = savePDF.width, height = savePDF.height)
  heatmap(ixns)
  dev.off()
  
  if (!(is.null(writeTo))){ write.csv(ixns, writeTo) }
  return(ixns)
  
}
