## this bit of the pipeline works out cellular interaction networks
#we'll probably start out with cellphone DB

## ---- SetWD
setwd("~/rotation1/")


## ---- LoadLibraries

library(data.table)
library(Seurat)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(biomaRt)
library(dplyr)
library(VennDiagram)
library(xlsx)
library(readxl)

source("code/GEN_ColorChooser.R")


## ---- SetVariables

#Polioudakis2019 <- readRDS("data/FET/Polioudakis2019/allData.RData")$SeuratOutput

results.folder <- "data/IXN/"
figure.folder <- "figures/IXN/"
fetus.source <- "Polioudakis2019" #set as "" for all datasets
CPDBResultsFolder <- paste0(results.folder, "cpdb_results/")
resource.folder <- "data/IXN/resources/"
OpenTargets.folder <- paste0(resource.folder, "open_targets/")

sample.ids <- c('fetus', 
                'Giandomenico2019_H9-75d-Batch1',
                'Velasco2019_HUES66-3mo-Batch1', 'Velasco2019_PGP1-3mo-Batch1', 'Velasco2019_PGP1-3mo-Batch2',
                'Velasco2019_11a-6mo-Batch1', 'Velasco2019_GM-6mo-Batch1', 'Velasco2019_PGP1-6mo-Batch1', 'Velasco2019_PGP1-6mo-Batch3')

cellTypes = c("End", "ExDp1", "ExDp2", "ExM", "ExM-U", "ExN",
              "InCGE", "InMGE", "IP", 
              "Mic", 
              "OPC", "oRG", 
              "Per", "PgG2M", "PgS", 
              "vRG")

LRMetadataColumns <- c('id_cp_interaction', 'interacting_pair',
                       'partner_a', 'partner_b',
                       'gene_a', 'gene_b',
                       'secreted', 'is_integrin',
                       'receptor_a', 'receptor_b',
                       'annotation_strategy',
                       'geneid_a', 'geneid_b', 'gene_a_summary', 'gene_b_summary')

openTargetDiseaseList <- gsub(gsub(gsub(list.files(OpenTargets.folder, pattern = "targets_associated"), pattern = "targets_associated_with_", replacement = ""),
                                   pattern = ".csv", replacement = ""),
                              pattern = "_", replacement = " ")
otherDiseases <- gsub(gsub(gsub(list.files(resource.folder, pattern = "targets_associated_with_"),
                                pattern = "targets_associated_with_", replacement = ""),
                           pattern = ".xlsx", replacement = ""),
                      pattern = "_", replacement = " ")



## for executing stuff

#Prepare4CPDB: converts Seurat Objects into count + metadata matrices for CPDB
.Prepare4CPDB = FALSE
Prepare4CPDB.fetusDataPath <- paste0("data/FET/", fetus.source, "/allData.RData")
Prepare4CPDB.organoidDataFolder <- "data/ORG/seurat_objects/"
Prepare4CPDB.writeTo <- paste0(results.folder, "cpdb_input/")
Prepare4CPDB.REDO <- FALSE

#CheckCPDBQuality: produces the number of 
.CheckCPDBQuality <- FALSE
CheckCPDBQuality.fetusDataPath <- Prepare4CPDB.fetusDataPath
CheckCPDBQuality.organoidFolder <- Prepare4CPDB.organoidDataFolder
CheckCPDBQuality.CPDBResultsFolder <- CPDBResultsFolder
CheckCPDBQuality.writeTo <- paste0(results.folder, "cpdb_counts.csv")


#ReformatCPDB: reformats CPDB output to something, as its default is a bit annoying
.ReformatCPDB = FALSE
ReformatCPDB.threshold <- 0.05
ReformatCPDB.resultsFolder <- CPDBResultsFolder
ReformatCPDB.writeTo <- paste0(results.folder, "reformatted_results/")

#CountInteractions: counts the number of interactions each cell type has with others
.CountInteractions = FALSE
CountInteractions.dataFolder <- ReformatCPDB.writeTo
CountInteractions.threshold <- ReformatCPDB.threshold
CountInteractions.writeTo <- paste0(results.folder, "interaction_counts/cell-interactions-byOrganoid.xlsx")


#InteractionsByCohort: list number of interactions per LR pair by cohort
.InteractionsByCohort <- TRUE
InteractionsByCohort.threshold <- 0.05
InteractionsByCohort.dataFolder <- paste0(results.folder, "cpdb_results/")
InteractionsByCohort.writeTo <- paste0(results.folder, "interaction_counts/cell-interactions-byCohort.csv")


#CombineCPDBResults: combine results from the 
.CombineCPDBResults = TRUE
CombineCPDBResults.CPDBResultsFolder <- CPDBResultsFolder
CombineCPDBResults.threshold <- 0.05
CombineCPDBResults.fw = 8 #1
if (CombineCPDBResults.fw == 1){CombineCPDBResults.writeTo <- paste0(results.folder, "interaction_counts/interactions-by-LR-pair.csv")} else {CombineCPDBResults.writeTo <- paste0(results.folder, "interaction_counts/interactions-by-LR-pair_weighted.csv")}

#.LabelNeuroDiseaseGenes
.LabelNeuroDiseaseGenes <- TRUE
LabelNeuroDiseaseGenes.interactionTablePath <- CombineCPDBResults.writeTo
LabelNeuroDiseaseGenes.threshold <- 0.05
LabelNeuroDiseaseGenes.writeTo = paste0(results.folder, "annotated_spreadsheets/LR_weighted-interactions_annotated.csv")

#.NormalizeCounts
.NormalizeCounts <- TRUE
NormalizeCounts.interactionTablePath <- LabelNeuroDiseaseGenes.writeTo
NormalizeCounts.writeTo = paste0(results.folder, "annotated_spreadsheets/LR_weighted-interactions_normalized.csv")

#.SubtypeSpecificityScores
.SubtypeSpecificityScores <- TRUE
SubtypeSpecificityScores.interactionTablePath <- LabelNeuroDiseaseGenes.writeTo
SubtypeSpecificityScores.psCount <- 1
SubtypeSpecificityScores.writeTo = paste0(results.folder, "annotated_spreadsheets/LR_weighted-interactions_subtype-scores.csv")



### ACTUALLY START DOING THINGS

## ---- IXN_RunCPDB
#CellPhone DB is actually run in Unix
#run "code/IXN_RunCPDB.sh" after

if (.Prepare4CPDB == TRUE){
  source("code/IXN_Prepare4CPDB.R")
  Prepare4CPDB(fetusData.path = Prepare4CPDB.fetusDataPath,
               organoidData.folder = Prepare4CPDB.organoidDataFolder,
               organoids = sample.ids,
               writeTo = Prepare4CPDB.writeTo, REDO = Prepare4CPDB.REDO)
}

rm(list = ls(pattern = "Prepare4CPDB"))


## ---- IXN_GeneCellStats

if (.CheckCPDBQuality == TRUE){
  source("code/IXN_CheckCPDBQuality.R")
  CheckCPDBQuality(fetusData.path = CheckCPDBQuality.fetusDataPath,
                   organoids = sample.ids,
                   organoidData.folder = CheckCPDBQuality.organoidFolder,
                   CPDB.resultsFolder = CheckCPDBQuality.CPDBResultsFolder,
                   writeTo = CheckCPDBQuality.writeTo)
}

rm(list = ls(pattern = "CheckCPDBQuality"))

## ---- IXN_ReformatCPDB

if (.ReformatCPDB == TRUE){
  source('code/IXN_ReformatCPDB.R')
  for (organoid in sample.ids){
    ReformatCPDB(organoid.id = organoid,
                 CPDB.resultsPath = ReformatCPDB.resultsFolder,
                 thresh = ReformatCPDB.threshold,
                 saveTo = paste0(ReformatCPDB.writeTo, organoid, ".csv"))
  }
  rm(organoid)
}

rm(list = ls(pattern = "ReformatCPDB"))

## ---- IXN_CountInteractions
#counts 

if (.CountInteractions == TRUE){
  source('code/IXN_CountInteractions.R')
  ixn.counts <- mapply(organoid.id = sample.ids, 
                       FUN = CountInteractions,
                       MoreArgs = list(threshold = CountInteractions.threshold,
                                       dataFolder = CountInteractions.dataFolder),
                       SIMPLIFY = FALSE)
  mapply(ixn.counts, sheetName=names(ixn.counts), FUN = write.xlsx, 
         MoreArgs = list(CountInteractions.writeTo, append=TRUE))
  rm(organoid); rm(ixn.counts)
}

rm(list = ls(pattern = "CountInteractions"))

## ---- IXN_InteractionsByCohort

if (.InteractionsByCohort == TRUE){
  source("code/IXN_InteractionsByCohort.R")
  InteractionsByCohort(sample.ids, threshold = InteractionsByCohort.threshold,
                       dataFolder = InteractionsByCohort.dataFolder,
                       writeTo = InteractionsByCohort.writeTo)
}

rm(list = ls(pattern = "InteractionsByCohort"))

## ---- IXN_CombineCPDBResults
#Combines the CPDB results from the different cohorts

if (.CombineCPDBResults == TRUE){
  source('code/IXN_CombineCPDBResults.R')
  CombineCPDBResults(organoid.ids = sample.ids, threshold = CombineCPDBResults.threshold,
                     CPDB.resultsFolder = CombineCPDBResults.CPDBResultsFolder,
                     fetus.weight = CombineCPDBResults.fw,
                     writeTo = CombineCPDBResults.writeTo)
}

rm(list = ls(pattern = "CombineCPDBResults"))

## ---- LabelNeuroDiseaseGenes
#annotates the LR interaction genes for relevant neurodegenerative diseases

if (.LabelNeuroDiseaseGenes == TRUE){
  source('code/IXN_LabelNeuroDiseaseGenes.R')
  LabelNeuroDiseaseGenes(interactionTable.path = LabelNeuroDiseaseGenes.interactionTablePath,
                         score_threshold = LabelNeuroDiseaseGenes.threshold,
                         writeTo = LabelNeuroDiseaseGenes.writeTo)
}

rm(list = ls(pattern = "LabelNeuroDiseaseGenes"))

## ---- NormalizeCounts
#divides LR interaction counts with total possible (necessary b/c not all cell types appear in all organoids)

if (.NormalizeCounts == TRUE){
  source('code/IXN_NormalizeCounts.R')
  NormalizeCounts(interactionTable.path = NormalizeCounts.interactionTablePath,
                  writeTo = NormalizeCounts.writeTo)
}

rm(list = ls(pattern = "NormalizeCounts"))

## ---- SubtypeSpecificityScores

if (.SubtypeSpecificityScores == TRUE){
  source('code/IXN_SubtypeSpecificityScores.R')
  SubtypeSpecificityScores(interactionTable.path = SubtypeSpecificityScores.interactionTablePath,
                         ps.count = SubtypeSpecificityScores.psCount,
                         writeTo = SubtypeSpecificityScores.writeTo)
}

rm(list = ls(pattern = "SubtypeSpecificityScores"))






