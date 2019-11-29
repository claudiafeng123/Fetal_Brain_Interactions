## ---- SetVariables

#interactionTable.path = paste0(results.folder, "annotated_spreadsheets/LR_interactions_annotated.csv") 
#writeTo = paste0(results.folder, "annotated_spreadsheets/LR_interactions_subtype-scores.csv")

## ---- RSC_SubtypeSpecificityScores-CODE


SubtypeSpecificityScores <- function(interactionTable.path = paste0(results.folder, "LR_interactions_annotated.csv"),
                                   ps.count = 1,
                                   writeTo = paste0(results.folder, "LR_interactions_subtype-scores.csv")){
  ixnTable <- fread(interactionTable.path)
  
  meta <- select(ixnTable, c(LRMetadataColumns,
                             'gene_a.anno', 'gene_b.anno'))
  counts <- as.matrix(select(ixnTable, colnames(ixnTable)[which(!(colnames(ixnTable) %in% c(LRMetadataColumns,
                                                                                            'gene_a.anno', 'gene_b.anno')))]))
  totals <- counts[1,]
  row.names(counts) <- meta$id_cp_interaction
  
  #non-cell-specific interactions
  non.counts <- mapply(1:dim(counts)[2], FUN = function(i){
    counts[, 1] - counts[, i]
  })
  colnames(non.counts) <- colnames(counts); row.names(non.counts) <- row.names(counts)
  
  rtn.counts <- counts/(non.counts+1)
  rtn.counts[1,] <- totals
  rtn.counts[which(is.nan(rtn.counts))] <- 0
  
  
  rtn <- cbind(meta, rtn.counts)
  if (is.null(writeTo) == FALSE){ write.csv(rtn, writeTo, row.names = FALSE) }
  return(rtn)
}



