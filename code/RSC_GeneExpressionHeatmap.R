
## ---- RSC_GeneExpressionHeatmap-CODE
#checks that the selected genes are actually expressed in the fetus
#it's okay if it isn't, but nice to look at which cells are expressing which of the chosen genes

GeneExpressionHeatmap <- function(fetus.data, geneList,
                                  savePDF = TRUE, savePDF.path = paste0(figure.folder, "panel-gene-expression.pdf"), savePDF.width = 6, savePDF.height = 5){
  cluster.averages <- AverageExpression(fetus.data, return.seurat = TRUE)
  rtn <- DoHeatmap(cluster.averages, features = geneList$gene, size = 3, draw.lines =  FALSE)
  if (savePDF == TRUE){pdf(savePDF.path, width = savePDF.width, height = savePDF.height)}
  rtn
  while(!(names(dev.cur()) %in% c('null device', 'RStudioGD'))){dev.off()}
  return(rtn)
}

