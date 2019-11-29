## ---- SetVariables

#organoid.id <- sample.ids[1]
#thresh <- 0.05

## ---- IXN_ReformatCPDB-CODE

ReformatCPDB <- function(organoid.id, CPDB.resultsPath = "data/IXN/cpdb_results/",
                         thresh = 0.05,
                         saveTo = paste0(ReformatCPDB.writeTo, organoid.id, ".csv")){
  
  
  p.values <- fread(paste0(CPDB.resultsPath, organoid.id, "/pvalues.txt"))
  means <- fread(paste0(CPDB.resultsPath, organoid.id, "/means.txt"))
  #deconvoluted <- fread(paste0(CPDB.resultsPath, organoid.id, "/deconvoluted.txt"))
  
  p.values.reformatted <- melt(p.values, id = grep(colnames(p.values), pattern = "[|]", value = TRUE, invert = TRUE))
  p.values.reformatted$`cell_a` <- unlist(lapply(strsplit(as.character(p.values.reformatted$variable), split = "[|]"), "[[", 1))
  p.values.reformatted$`cell_b` <- unlist(lapply(strsplit(as.character(p.values.reformatted$variable), split = "[|]"), "[[", 2))
  p.values.reformatted$`p_value` <- p.values.reformatted$value
  p.values.reformatted$variable <- p.values.reformatted$value <- NULL
  
  means.reformatted <- melt(means, id = grep(colnames(means), pattern = "[|]", value = TRUE, invert = TRUE))
  means.reformatted$`cell_a` <- unlist(lapply(strsplit(as.character(means.reformatted$variable), split = "[|]"), "[[", 1))
  means.reformatted$`cell_b` <- unlist(lapply(strsplit(as.character(means.reformatted$variable), split = "[|]"), "[[", 2))
  means.reformatted$`mean` <- means.reformatted$value
  means.reformatted$variable <- means.reformatted$value <- NULL
  
  cellphoneDB.result <- merge(p.values.reformatted, means.reformatted, all = TRUE)
  cellphoneDB.result$significant <- (cellphoneDB.result$p_value < thresh)
  
  #convert ensembl id to gene name
  cellphoneDB.result$geneid_a <- unlist(lapply(strsplit(cellphoneDB.result$interacting_pair, split = "_"), "[[", 1))
  cellphoneDB.result$geneid_b <- unlist(lapply(strsplit(cellphoneDB.result$interacting_pair, split = "_"), "[[", 2))
  cellphoneDB.result$gene_a_summary <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", cellphoneDB.result$geneid_a)
  cellphoneDB.result$gene_b_summary <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", cellphoneDB.result$geneid_b)
  
  #rearrange order of columns
  
  rtn <- cbind(cellphoneDB.result[, c('cell_a', 'cell_b')],
               cellphoneDB.result[, c('geneid_a', 'geneid_b')],
               cellphoneDB.result[, c('mean', 'p_value', 'significant')],
               cellphoneDB.result[, c('id_cp_interaction', 'gene_a', 'gene_b', 'partner_a', 'partner_b')],
               cellphoneDB.result[, c('secreted', 'receptor_a', 'receptor_b', 'annotation_strategy', 'is_integrin')],
               cellphoneDB.result[, c('gene_a_summary', 'gene_b_summary')])
  names(rtn) <- c('cell_type_a', 'cell_type_b',
                  'gene_a', 'gene_b',
                  'mean', 'p_value', 'significant',
                  'cp_interaction_id', 'ensembl_id_a', 'ensembl_id_b', 'partner_a', 'partner_b',
                  'secreted', 'receptor_a', 'receptor_b', 'annotation_strategy', 'is_integrin',
                  'gene_a_summary', 'gene_b_summary')
  
  write.csv(rtn, saveTo, row.names = FALSE)
  
  return(rtn)
  
}



