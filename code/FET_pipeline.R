## this bit of the pipeline preprocesses the data from the Nowakowski 2017 and Polioudakis 2019 studies
## steps can be outlined as follows:

## 1.) 

## ---- SetWD
setwd("~/rotation1/")

## ---- SetVariables

#these variables will not be deleted after each step
studies <- c('Nowakowski2017', 'Polioudakis2019')
main.study <- "Polioudakis2019"

ndims = 30
min.cells = 3
min.features = 200

saveRDS.folder <- "data/FET/"
fig.folder <- "figures/FET/"

#these  variables will be deleted after they are no longer needed
##Cell names should be in the form:
## CELLID_TYPE_DONOR

#.LoadData
.LoadData <- FALSE
#Polioudakis2019
Polioudakis2019.folderPath <- paste0(saveRDS.folder, "Polioudakis2019") #this should NOT have a backslash after it
Polioudakis2019.data.URL <- "http://solo.bmap.ucla.edu/shiny/webapp/session/638ff98242f1a20b5cb4143847f7152c/download/downloadData?w="
Polioudakis2019.data.downloadSuffix <- ".zip"
Polioudakis2019.rawCounts.matrixPath <- paste0(Polioudakis2019.folderPath, "/raw_counts_mat.rdata")
Polioudakis2019.metadata.matrixPath <- paste0(Polioudakis2019.folderPath, "/cell_metadata.csv")
#Nowakowski2017
Nowakowski2017.folderPath <- paste0(saveRDS.folder, "Nowakowski2017/")
Nowakowski2017.rawCounts.URL <-"https://cells.ucsc.edu/cortex-dev/exprMatrix.tsv.gz"
Nowakowski2017.rawCounts.downloadSuffix <- ".gz"
Nowakowski2017.rawCounts.matrixPath <- paste0(Nowakowski2017.folderPath, "Nowakowski2017_rawCounts.tsv")
Nowakowski2017.metadata.URL <- "https://cells.ucsc.edu/cortex-dev/meta.tsv"
Nowakowski2017.metadata.matrixPath <- paste0(Nowakowski2017.folderPath, "Nowakowski2017_metadata.tsv")

#.PredictLabels
.PredictLabels <- FALSE
cutoff <- 0.8

#.FilterCells
.FilterCells <- FALSE
FilterCells.DataPath <- paste0(saveRDS.folder, "allData_PredictedLabels.RData")
predictionScoreRatio.cutoff <- cutoff

#.IntegrateDatasets
.IntegrateDatasets <- FALSE
IntegrateCells.DataPath <- paste0(saveRDS.folder, "allData_FilterCells.RData")
IntegrateDatasets.saveTo <- paste0(saveRDS.folder, "allData.RDS")

#.VisualizeData
.VisualizeData <- TRUE
features = c("cell_type", "source")#these must match the column name in metadata


## ---- LoadLibraries

library(Matrix)
library(dplyr)
library(Seurat)
library(data.table)
library(R.utils)
library(ggplot2)
library(cowplot)
library(scales)

source("code/GEN_ColorChooser.R")


### ACTUALLY START DOING THINGS

## ---- LoadData
#returns two matrices for each study: metadata and raw counts
#also returns something named SeuratObject which is the normalized and qc'd data that should be ready to be integrated

if (.LoadData == TRUE){
  #Polioudakis2019
  source("code/FET_LoadPolioudakis2019.R")
  Polioudakis2019 <- loadPolioudakis2019(study = "Polioudakis2019", folderPath = Polioudakis2019.folderPath,
                                        data.URL = Polioudakis2019.data.URL, data.downloadSuffix = Polioudakis2019.data.downloadSuffix,
                                        rawCounts.matrixPath = Polioudakis2019.rawCounts.matrixPath, metadata.matrixPath = Polioudakis2019.metadata.matrixPath,
                                        saveTo = paste0(Polioudakis2019.folderPath, "/allData.RData"))
  #Nowakowski2017
  source("code/FET_LoadNowakowski2017.R")
  Nowakowski2017 <- loadNowakowski2017(study = "Nowakowski2017", folderPath = Nowakowski2017.folderPath,
                                        rawCounts.URL = Nowakowski2017.rawCounts.URL, rawCounts.downloadSuffix = Nowakowski2017.rawCounts.downloadSuffix, 
                                        metadata.URL = Nowakowski2017.metadata.URL,
                                        rawCounts.matrixPath = Nowakowski2017.rawCounts.matrixPath, metadata.matrixPath = Nowakowski2017.metadata.matrixPath,
                                        saveTo = paste0(Nowakowski2017.folderPath, "allData.RData"))
  rm("loadPolioudakis2019"); rm("loadNowakowski2017")
}

rm(list = ls(pattern = "Polioudakis2019[.]")); rm(list = ls(pattern = "Nowakowski2017[.]"))


## ---- PredictLabels
#uses Seurat's TransferData() function to predict labels based on a reference dataset

if (.PredictLabels == TRUE){
  source("code/FET_PredictedLabels.R")
  for (study in studies){eval(parse(text = paste0(
    "if (!(exists('", study, "'))){", study, " <- readRDS('", saveRDS.folder, study, "/allData.RData')}"
  )))}
  SeuratList <- PredictLabels(studies = studies, main.study = main.study,
                              secondBest.cutoff = cutoff,
                              REDO = FALSE, SAVE = FALSE)
  rm(PredictLabels)
}

rm(.PredictLabels); rm(predictionScoreRatio.cutoff)

## ---- FilterCells

if (.FilterCells == TRUE){
  source("code/FET_FilterCells.R")
  if (!(exists("SeuratList"))) { SeuratList <- readRDS(FilterCells.DataPath) }
  SeuratList.Filtered <- FilterCells(SeuratList, main.study = main.study)
  rm(FilterCells)
}

rm(.FilterCells); rm(FilterCells.DataPath)

## ---- IntegrateDatasets
#uses Seurat to integrate the datasets available
#also runs pca + UMAP to avoid creating too many RData objects

if (.IntegrateDatasets == TRUE){
  source("code/FET_IntegrateDatasets.R")
  if (!(exists("SeuratList.Filtered"))) { SeuratList.Filtered <- readRDS(IntegrateCells.DataPath) }
  allData <- IntegrateDatasets(SeuratObjectList = SeuratList.Filtered, 
                               saveTo = IntegrateDatasets.saveTo, REDO = TRUE, SAVE = TRUE)
  rm(IntegrateDatasets)
}

rm(.IntegrateDatasets); rm(IntegrateDatasets.DataPath); rm(IntegrateDatasets.saveTo)

## ---- DataVisualization

if (.VisualizeData == TRUE){
  source("code/FET_VisualizeData.R")
  if (!exists("allData")){ allData = readRDS(allData.path) }
  VisualizeData(allData = allData, .plotBy = features)
  rm(VisualizeData)
}


## ---- FET_Rider

rm(Polioudakis2019); rm(Nowakowski2017)
rm(SeuratList); rm(SeuratList.Filtered)

rm(studies); rm(main.study)
rm(ndims); rm(min.cells); rm(min.features)
rm(saveRDS.folder); rm(fig.folder)

