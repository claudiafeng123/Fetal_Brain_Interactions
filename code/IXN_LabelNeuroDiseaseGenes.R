
## ---- SetVariables

interactionTable.path <- "data/IXN/interaction_counts/interactions-by-LR-pair.csv"

disease <- openTargetDiseaseList[1]
disease <- openTargetDiseaseList[2]
score_threshold = 0.05

## ---- LabelNeuroDiseaseGenes

LabelNeuroDiseaseGenes <- function(interactionTable.path = "data/IXN/interaction_counts/interactions-by-LR-pair.csv",
                                   score_threshold = 0.05, 
                                   writeTo = paste0("data/RSC/LR_interactions_annotated.csv")){
  
  ixnTable <- fread(interactionTable.path)
  counts <- select(ixnTable, colnames(ixnTable)[which(!(colnames(ixnTable) %in% LRMetadataColumns))])
  
  meta <- select(ixnTable, LRMetadataColumns)
  meta <- meta[-1,]
  
  gene_a.anno_ALL <- gene_b.anno_ALL <- rep("", dim(meta)[1])
  
  for (disease in openTargetDiseaseList){
    
    riskGenesTable <- fread(paste0(resource.folder, "open_targets/targets_associated_with_", gsub(disease, pattern = " ", replacement = "_"), ".csv"))
    riskGenes <- as.character(riskGenesTable$target.gene_info.symbol[riskGenesTable$association_score.datatypes.genetic_association > score_threshold])
    riskGenes.ensID <- as.character(riskGenesTable$target.id[riskGenesTable$association_score.datatypes.genetic_association > score_threshold])
    
    riskGeneScores <- as.character(round(riskGenesTable$association_score.datatypes.genetic_association,2)[riskGenesTable$association_score.datatypes.genetic_association > score_threshold])
    
    diseaseAnnotation <- data.frame(gene = riskGenes,
                                    ensembl.id = riskGenes.ensID,
                                    risk = paste0(disease, ": ", riskGeneScores), stringsAsFactors = FALSE)
    
    gene.a_anno <- mapply(meta$gene_a, FUN = function(x){
      if (x %in% diseaseAnnotation$ensembl.id){ 
        ind <-  which(x == diseaseAnnotation$ensembl.id)
        return(paste0(diseaseAnnotation$risk[ind], "; "))} else {return("")}
    })
    
    gene.b_anno <- mapply(meta$gene_b, FUN = function(x){
      if (x %in% diseaseAnnotation$ensembl.id){ 
        ind <-  which(x == diseaseAnnotation$ensembl.id)
        return(paste0(diseaseAnnotation$risk[ind], "; "))} else {return("")}
    })
    
    gene_a.anno_ALL <- paste0(gene_a.anno_ALL, gene.a_anno); gene_b.anno_ALL <- paste0(gene_b.anno_ALL, gene.b_anno)
    
  }
  
  ## add brain size
  for (disease in otherDiseases){
    if (disease == "Brain size"){
      riskGenesTable <- read_excel(paste0(resource.folder, list.files(path = resource.folder, pattern = gsub(disease, pattern = " ", replacement = "_"))), skip = 1)
      riskGenesTable <- riskGenesTable[!(duplicated(riskGenesTable$GENE)),]
      riskGenes <- riskGenesTable$GENE
      diseaseAnnotation <- data.frame(gene = riskGenesTable$SYMBOL,
                                      ensembl.id = riskGenesTable$GENE,
                                      risk = disease)
    }
    
    gene.a_anno <- mapply(meta$gene_a, FUN = function(x){
      if (x %in% diseaseAnnotation$ensembl.id){ 
        ind <-  which(x == diseaseAnnotation$ensembl.id)
        return(paste0(diseaseAnnotation$risk[ind], "; "))} else {return("")}
    })
    
    gene.b_anno <- mapply(meta$gene_b, FUN = function(x){
      if (x %in% diseaseAnnotation$ensembl.id){ 
        ind <-  which(x == diseaseAnnotation$ensembl.id)
        return(paste0(diseaseAnnotation$risk[ind], "; "))} else {return("")}
    })
    
    gene_a.anno_ALL <- paste0(gene_a.anno_ALL, gene.a_anno); gene_b.anno_ALL <- paste0(gene_b.anno_ALL, gene.b_anno)
    
  }
  
  
  
  meta <- data.frame(meta, gene_a.anno = gene_a.anno_ALL, gene_b.anno = gene_b.anno_ALL, stringsAsFactors = FALSE)
  meta <- rbind(rep('', length(LRMetadataColumns)), meta)
  
  rtn <- cbind(meta, counts)
  if (is.null(writeTo) == FALSE){ write.csv(rtn, writeTo, row.names = FALSE)}
  return(rtn)
  
}

