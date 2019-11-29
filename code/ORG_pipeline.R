## this bit of the pipeline integrates organoid data with reference fetal data
## steps can be outlined as follows:

## 1.) 

## ---- SetWD
setwd("~/rotation1/")


## ---- SetVariables

ndims = 30

#fetus.data should be a Seurat Object for the reference sc data
fetus.data <- readRDS("data/FET/Polioudakis2019/allData.RData")
cellTypeList <- sort(unique(as.character(Idents(fetus.data))))

organoid.filePaths <- paste0("data/ORG/raw_data/", list.files("data/ORG/raw_data/"))
organoid.list <- gsub(unlist(lapply(strsplit(organoid.filePaths, split = "/"), "[[", 4)), pattern = ".txt", replacement = "")

#folder where all figures are saved
savePDF.folder <- "figures/ORG/"
savedData.folder <- "data/ORG/"


####
## step-specific variables

#.DownloadData
#already done
#you'll have to go into the code to redo any sort of thing you want changec
.DownloadData <- FALSE

#.LabelCellTypes
.LabelCellTypes <- FALSE
batch.seps <- list("Giandomenico2019" = "-", "Velasco2019" = "_")

#.EvaluateMappedLabels
.EvaluateMappedLabels <- FALSE
organoidData.paths <- paste0("data/ORG/seurat_objects/", list.files("data/ORG/seurat_objects/"))
features = c("cell_line", "time_point", "predicted.id")


#.PlotUMAP
.PlotUMAP <- FALSE

#.FetusOrganoidUMAP
.FetusOrganoidUMAP <- TRUE
organoidSeurat.paths <- organoidData.paths

#.QuantifyCellTypes
.QuantifyCellTypes <- FALSE
cell.counts.SaveTo <- paste0("data/ORG/", "cell-counts.RData")


## ---- LoadLibraries

library(data.table)
library(R.utils)
library(stringr)
library(Seurat)
library(scales)
library(dplyr)
library(ggplot2)
library(doParallel)

source("code/GEN_ColorChooser.R")

### ACTUALLY START DOING THINGS

## ---- ORG_DownloaData

if (.DownloadData == TRUE){
  source("code/ORG_DownloadData.R")
}

rm(.DownloadData)

## ---- ORG_LabelCellTypes

if (.LabelCellTypes == TRUE){
  source("code/ORG_LabelCellTypes.R")
  for (i in 1:length(organoid.list)){
    study = unlist(lapply(strsplit(organoid.list[i], split = "_"), "[[", 1))
    LabelCellTypes(organoid.id = organoid.list[i], organoid.exprMatrix.path = organoid.filePaths[i],
                   fetus.data = fetus.data,
                   ndims = ndims, batch.sep = eval(parse(text = (paste0("batch.seps$", study)))),
                   saveTo = paste0(savedData.folder, "seurat_objects/", organoid.list[i], '.RData'),
                   REDO = TRUE)
  }
  rm(i)
  rm(LabelCellTypes)
}

rm(.LabelCellTypes)
rm(batch.seps)

## ---- ORG_EvaluateMappedLabels

if (.EvaluateMappedLabels == TRUE){
  source("code/ORG_EvaluateMappedLabels.R")
  organoid.metadata <- EvaluateMappedLabels(mappedObjects.paths = organoidData.paths, saveMetadata = FALSE)
  rm(EvaluateMappedLabels)
}

rm(.EvaluateMappedLabels)
rm(features); rm(organoidData.paths)

## ---- ORG_PlotUMAP

if (.PlotUMAP == TRUE){
  source('code/ORG_PlotUMAP.R')
  for (i in 1:length(organoid.list)){
    PlotUMAP(organoid.id = organoid.list[i],
             savePDF.path = paste0(savePDF.folder, "UMAP_", organoid.list[i], ".pdf"))
  }
  rm(i)
  rm(PlotUMAP)
}

rm(.PlotUMAP)

## ---- ORG_FetusOrganoidUMAP

if (.FetusOrganoidUMAP == TRUE){
  source('code/ORG_FetusOrganoidUMAP.R')
  FetusOrganoidUMAP(organoidData.paths = organoidSeurat.paths, fetus.data = fetus.data,
                    REDO = FALSE)
  rm(FetusOrganoidUMAP)
}

rm(organoidSeurat.paths)


## ---- ORG_QuantifyCellTypes

if (.QuantifyCellTypes == TRUE){
  source('code/ORG_QuantifyCellTypes.R')
  if (!exists("organoid.metadata")){ organoid.metadata <- readRDS(paste0(savedData.folder, "organoid-metadata.RData")) }
  cell.counts <- QuantifyCellTypes(organoid.metadata = organoid.metadata,
                                        fetus.metadata = fetus.data@meta.data, cellTypeList = cellTypeList,
                                   saveRDS.path = cell.counts.SaveTo)
  rm(QuantifyCellTypes)
}

rm(cell.counts.SaveTo)
rm(.QuantifyCellTypes)

## ---- ORG_Rider
rm(savedData.folder); rm(savePDF.folder)
rm(organoid.list); rm(organoid.filePaths)
