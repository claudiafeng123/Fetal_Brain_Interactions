## ---- SetVariables

#dataFolder = InteractionsByCohort.dataFolder
#writeTo = InteractionsByCohort.writeTo
#organoid.id <- sample.ids[5]
#threshold = 0.05
#i <- 5
#dataList <- LRCounts

## ---- IXN_InteractionsByCohort-CODE

InteractionsByCohort <- function(sample.ids,
                                 threshold = 0.05,
                                 dataFolder = InteractionsByCohort.dataFolder,
                                 writeTo = InteractionsByCohort.writeTo){
  
  rtn.tracker <- 0
  for (i in 1:length(sample.ids)){
    
    organoid.id = sample.ids[i]
    
    p.values <- fread(paste0(dataFolder, organoid.id, "/pvalues.txt"))
    means <- fread(paste0(dataFolder, organoid.id, "/means.txt"))
    
    meta <- eval(parse(text = paste0(
      "p.values[," ,
      "c(", paste0(grep(colnames(p.values), pattern = "[|]", invert = TRUE), collapse = ", "), ")",
      "]"
    )))
    
    sig.means <- mapply(i=grep(colnames(p.values), pattern = "[|]", value = TRUE), FUN = function(i, thresh = threshold){
      rtn <- rep(1, dim(p.values)[1])
      eval(parse(text = paste0("rtn[which(p.values[,'", i, "'] > thresh)] <- 0")))
      return(rtn)
    })
    
    names(sig.means) <- grep(colnames(p.values), pattern = "[|]", value = TRUE)
    total_ixns_pLR <- apply(sig.means, 1, sum)
    names(total_ixns_pLR) <- p.values$id_cp_interaction
    sig.ixns <- cbind(meta, total = total_ixns_pLR)
    
    names(sig.ixns)[dim(sig.ixns)[2]] <- organoid.id
    if (rtn.tracker == 0){rtn <- sig.ixns} else {rtn <- base::merge(rtn, sig.ixns, all = TRUE)}
    rtn.tracker <- rtn.tracker + 1
    
  }
  rtn[is.na(rtn)] <- 0
  #rtn <- as.data.frame(apply(rtn, 2, FUN = function(x){x[which(is.na(x))] <- 0; return(x)}))
  
  lr.counts <- eval(parse(text = paste0(
    "rtn[," ,
    "c(", paste0(which(names(rtn) %in% sample.ids), collapse = ", "), ")",
    "]"
  )))
  rtn$average <- apply(lr.counts, 1, sum)/length(sample.ids)
  rtn$var <- apply(lr.counts, 1, var)
  
  write.csv(rtn, writeTo, row.names = FALSE)
  return(rtn)
  
}


