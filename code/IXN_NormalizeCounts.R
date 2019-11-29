## ---- SetVariables

#interactionTable.path = paste0(results.folder, "LR_interactions_annotated.csv") 

## ---- RSC_NormalizeCounts-CODE

NormalizeCounts <- function(interactionTable.path = paste0(results.folder, "LR_interactions_annotated.csv"),
                            writeTo = paste0(results.folder, "LR_interactions_normalized.csv")){
  
  ixnTable <- fread(interactionTable.path)
  
  counts <- select(ixnTable, colnames(ixnTable)[which(!(colnames(ixnTable) %in% c(LRMetadataColumns,
                                                                                  'gene_a.anno', 'gene_b.anno')))])
  
  totals <- as.vector(counts[1,])
  counts <- counts[-1,]
  counts <- mapply('/', counts, totals)
  
  meta <- select(ixnTable, c(LRMetadataColumns,
                             'gene_a.anno', 'gene_b.anno'))
  
  rtn <- cbind(meta[-1,], counts)
  if (is.null(writeTo) == FALSE){ write.csv(rtn, writeTo, row.names = FALSE) }
  return(rtn)
  
}

