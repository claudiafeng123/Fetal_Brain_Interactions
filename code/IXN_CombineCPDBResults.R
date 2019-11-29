## ---- SetVariables

#CPDB.resultsFolder <- "data/IXN/cpdb_results/"
#organoid.ids = sample.ids
#organoid <- organoid.ids[2]
#threshold = 0.05

## ---- IXN_CombineCPDBResults-CODE

CombineCPDBResults <- function(organoid.ids = sample.ids, threshold = 0.05,
                               CPDB.resultsFolder = "data/IXN/cpdb_results/",
                               fetus.weight = 1,
                               writeTo = "data/IXN/interaction_counts/interactions-by-LR-pair.csv"){
  interactionPairList <- c()
  for (cty in cellTypes){
    interactionPairList <- sort(unique(c(interactionPairList, paste(cty, cellTypes, sep = "|"))))
  }
  rtn <- matrix(0, nrow = 0, ncol = length(interactionPairList))
  colnames(rtn) <- interactionPairList
  totals <- rep(0, length(interactionPairList)); names(totals) <- interactionPairList
  
  LRList <- data.frame(id_cp_interaction = character(0), 
                       interacting_pair = character(0),
                       partner_a  = character(0),
                       partner_b  = character(0),
                       gene_a = character(0),
                       gene_b = character(0),
                       secreted  = logical(0),
                       receptor_a = logical(0),
                       receptor_b = logical(0),
                       annotation_strategy = character(0),
                       is_integrin = logical(0))
  columnNames <- colnames(LRList)
  
  #count total number of interactions per ligand by cell
  for (organoid in organoid.ids){
    
    p.values <- fread(paste0(CPDB.resultsFolder, organoid, "/pvalues.txt"))
    means <- fread(paste0(CPDB.resultsFolder, organoid, "/means.txt"))
    meta <- eval(parse(text = paste0(
      "p.values[," ,
      "c(", paste0(grep(colnames(p.values), pattern = "[|]", invert = TRUE), collapse = ", "), ")",
      "]"
    )))
    LRList <- base::merge(LRList, meta, all = TRUE)
    
    sig.ixns <- mapply(i=grep(colnames(p.values), pattern = "[|]", value = TRUE), FUN = function(i, thresh = threshold){
      if (organoid == "fetus"){ rtn <- rep(fetus.weight, dim(p.values)[1]) } else {rtn <- rep(1, dim(p.values)[1])}
      eval(parse(text = paste0("rtn[which(p.values[,'", i, "'] > thresh)] <- 0")))
      return(rtn)
    })
    if (organoid == "fetus"){ totals[which(names(totals) %in% colnames(sig.ixns))] <- totals[which(names(totals) %in% colnames(sig.ixns))] + 1 + (fetus.weight - 1) } else {
      totals[which(names(totals) %in% colnames(sig.ixns))] <- totals[which(names(totals) %in% colnames(sig.ixns))] + 1
    }
    
    sig.ixns <- base::merge(LRList, cbind(meta, sig.ixns), by = columnNames, all = TRUE)
    for (ixn.pair in interactionPairList){
      if (!(ixn.pair %in% colnames(sig.ixns))){ eval(parse(text = paste0("sig.ixns$`", ixn.pair, "` <- 0")))}
    }; sig.ixns <- sig.ixns[, c(columnNames, interactionPairList)]
    sig.ixns <- replace(sig.ixns,is.na(sig.ixns),0)
    if (dim(rtn)[1] == 0){ rtn <- sig.ixns } else {
      rtn <- base::merge(LRList, rtn, by = c(columnNames), all = TRUE)
      rtn[is.na(rtn)] <- 0
      
      numbers.1 <- sig.ixns[, -(1:length(columnNames))]
      numbers.2 <- rtn[, -(1:length(columnNames))]
      
      counts <- numbers.1 + numbers.2
      
      rtn <- cbind(LRList, counts)
      
    }
    
  }
  
  rtn <- rbind(c(rep(NA, length(columnNames)), totals), rtn)
  counts <- rtn[, -(1:length(columnNames))]
  total <- apply(counts, 1, sum)
  
  #calculate total interactions per cell
  total.by.cell <- matrix(0, ncol = length(cellTypes), nrow = dim(rtn)[1]); colnames(total.by.cell) <- cellTypes
  for (i in 1:length(cellTypes)){
    cty <- cellTypes[i]
    mat <- eval(parse(text = paste0(
      "counts[," ,
      "c(", paste0(grep(colnames(counts), pattern = cty), collapse = ", "), ")",
      "]"
    )))
    total.by.cell[,i] <- apply(mat, 1, sum)
  }
  
  #convert ensembl id to gene name
  LRList$geneid_a <- unlist(lapply(strsplit(LRList$interacting_pair, split = "_"), "[[", 1))
  LRList$geneid_b <- unlist(lapply(strsplit(LRList$interacting_pair, split = "_"), "[[", 2))
  LRList$gene_a_summary <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", LRList$geneid_a)
  LRList$gene_b_summary <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", LRList$geneid_b)
  
  
  LRList <- rbind(c('', ''), LRList)
  rtn <- cbind(LRList, total, total.by.cell, counts)
  if (!(is.null(writeTo))){ write.csv(rtn, writeTo, row.names = FALSE)}
  
  return(rtn)
}






