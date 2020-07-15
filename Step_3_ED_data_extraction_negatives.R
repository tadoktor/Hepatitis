#Author: Thomas Exner


compounds <- NULL
studyType <- NULL
organism <- NULL
organ <- NULL
dosing <- NULL
dose <- NULL
duration <- NULL
durationUnit <- NULL

compounds <- list("bromoethylamine","bucetin","chlormadinone","desmopressin acetate","fluoxetine hydrochloride","phenylanthranilic acid","
tannic acid","caffeine","metformin","theophylline")
studyType <- list("in_vivo")
organism <- list("Rat")
organ <- list("Liver")
#dosing <- list("Single")
dose <- list("high")
#duration <- list(9)
duration <- NULL
#durationUnit <- list("hr")

TGG_metadata <- "/OutputFiles/TGG_metadata_negatives.txt"
TGG_foldChanges <- "/OutputFiles/TGG_foldchanges_negatives.txt"
library(magrittr)
library(jsonlite)
library(httr)
set_config(config(ssl_verifypeer = 0L))
library(data.table)
library(tidyr)
library(RCurl)
library(XML)
library(stringi)
Root_location <- "C:\\Users\\tatyana\\Documents\\Projects\\Advance project\\Advance"

url <- "https://api.kit.edelweissconnect.com/datasets"
req <- fromJSON(url, flatten = TRUE)
studySet <- req$results[grepl("TG-GATES_STUDY", req$results[, 'name'], fixed=TRUE),]
studySet

url <- paste("https://api.kit.edelweissconnect.com/datasets/", studySet$id.id, "/versions/", studySet$id.version, "/data?limit=None", sep="")
req <- fromJSON(url, flatten = TRUE)
#head(req$results,2)

if (!is.null(compounds)) {filtered <- req$results[req$results$data.Compound %in% compounds,]}
if (!is.null(studyType)) {filtered <- filtered[filtered$"data.Study type" %in% studyType,]}
if (!is.null(organism)) {filtered <- filtered[filtered$data.Organism %in% organism,]}
if (!is.null(organ)) {filtered <- filtered[filtered$data.Organ %in% organ,]}
if (!is.null(dosing)) {filtered <- filtered[filtered$data.Dosing %in% dosing,]}
if (!is.null(dose)) {filtered <- filtered[filtered$data.Dose %in% dose,]}
if (!is.null(duration)) {filtered <- filtered[filtered$data.Duration %in% duration,]}
if (!is.null(durationUnit)) {filtered <- filtered[filtered$"data.Duration unit" %in% durationUnit,]}
filtered = filtered[!duplicated(filtered$"data.Fold changes dataset ID"),]
filtered


sampleInfo <- data.frame(matrix(NA, nrow = nrow(filtered), ncol = 0))
sampleInfo$TGG_compoundName <- filtered$data.Compound
sampleInfo$TGG_ID <- 0
sampleInfo$Casrn <- filtered$data.CAS 
sampleInfo$"_id_" <- filtered$id 
sampleInfo$sampleid <- filtered$id
sampleInfo$timepointHr <- filtered$data.Duration
sampleInfo$organism <- filtered$data.Organism
sampleInfo$tissue <- filtered$data.Organ
sampleInfo$doseLevel <- filtered$data.Dose
sampleInfo$repeatType <- filtered$data.Dosage
sampleInfo$cellType <- filtered$"data.Study type"

sampleInfo

append = FALSE
j = 0
for(i in rownames(filtered)) {
  j = j+1
  if (file.exists(paste(Root_location,TGG_foldChanges,sep=''))) {
    metadata <- read.csv(paste(Root_location,TGG_foldChanges,sep=''),sep="\t", header=TRUE)
    if (i %in% metadata$sampleid) {
      next
    }
  }
  url = paste("https://api.kit.edelweissconnect.com/datasets/", filtered[i,"data.Fold changes dataset ID"], "/versions/", filtered[i, "data.Fold changes dataset version"], "/data", sep="")
  req <- fromJSON(url, flatten = TRUE)
  foldChanges <- data.frame(matrix(NA, nrow = nrow(req$results), ncol = 0))
  foldChanges$sampleid <- filtered[i,"id"]
  foldChanges$assayid <- req$results$data.PROBEID
  foldChanges$geneSymbols <- req$results$data.SYMBOL
  foldChanges$value <- req$results$data.logFC
  foldChanges$valueType <- "log2fold"
  foldChanges$pvalue <- req$results$data.P.Value
  if (j==1) {
    write.table(foldChanges, paste(Root_location,TGG_foldChanges,sep=''),sep="\t",row.names=F,col.names=T, append=FALSE)
    write.table(sampleInfo[j,], paste(Root_location,TGG_metadata,sep=''),sep="\t",row.names=F,col.names=T, append=FALSE)
    append = TRUE
  } else {
    write.table(foldChanges, paste(Root_location,TGG_foldChanges,sep=''),sep="\t",row.names=F,col.names=F, append=TRUE)
    write.table(sampleInfo[j,], paste(Root_location,TGG_metadata,sep=''),sep="\t",row.names=F,col.names=F, append=TRUE)
  }
}

