## this bit of the pipeline preprocesses the data from the Nowakowski 2017 and Polioudakis 2019 studies
## steps can be outlined as follows:

## 1.) 

## ---- SetWD
setwd("~/rotation1/")

## ---- LoadLibraries

library(data.table)
library(Seurat)
library(ggplot2)


## ---- SetVariables

fetus.source <- "Polioudakis2019"
geneList <- fread("data/RSC/panel-genes.txt")
figure.folder <- "figures/RSC/"

fetus.data <- readRDS(paste0("data/FET/", fetus.source, "/allData.RData"))



#GeneExpressionHeatmap
.GeneExpressionHeatmap = TRUE
GeneExpressionHeatmap.savePDFPath = paste0(figure.folder, "panel-gene-expression.pdf")

### ACTUALLY START DOING THINGS

## ---- RSC_GeneExpressionHeatmap

if (.GeneExpressionHeatmap == TRUE){
  source('code/RSC_GeneExpressionHeatmap.R')
  GeneExpressionHeatmap(geneList = geneList, fetus.data = fetus.data,
                        savePDF.path = GeneExpressionHeatmap.savePDFPath)
}

rm(list = ls(pattern = "GeneExpressionHeatmap"))

#need to create dot plot in terminal

