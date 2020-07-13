OutDir <- paste(Root_location, "/IntermediateData/", sep='')

# Install libraries
#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library(data.table)
library(reshape2)
library("biomaRt")
library("dplyr")
library(collapsibleTree)

#Generation of GTX list of genes


#Generation of GTX list of genes for rat in vivo
#The files can be derived either from the previous Step_3B1, Step_3B2, Step3_B3 or directly downloaded from this link: https://mega.nz/#F!zDBmiKyZ

#Generation of PLP list of genes for rat in vivo
#The files can be derived either from the previous Step_3B1, Step_3B2, Step3_B3 or directly downloaded from this link: https://mega.nz/#F!zDBmiKyZ

TG_Gates_rat_log2fold_GTX<- read.table("C:/Users/tatyana/Documents/Projects/Advance project/Advance/OutputFiles/TGG_foldchanges_hepatitis.txt", header=TRUE)
GTX_sampleID<- read.table("C:/Users/tatyana/Documents/Projects/Advance project/Advance/OutputFiles/TGG_metadata_hepatitis.txt", header=TRUE)

GTX_sampleID_filter_2<-TG_Gates_rat_log2fold_GTX

GTX_sampleID$TGG_ID<- NULL
GTX_sampleID$Casrn<- NULL
GTX_sampleID$X_id_<- NULL


TG_Gates_rat_log2fold_GTX$assayId <- NULL
TG_Gates_rat_log2fold_GTX$pvalue <- NULL
TG_Gates_rat_log2fold_GTX$valueType <- NULL

GTX_gene_search<-merge(TG_Gates_rat_log2fold_GTX,GTX_sampleID, by.x="sampleid",by.y="sampleid")

TG_Gates_rat_log2fold_NGTX<- read.table("file:///C:/Users/tatyana/Documents/Projects/Advance project/Advance/OutputFiles/TGG_foldchanges_negatives.txt", header=TRUE)
NGTX_sampleID<- read.table("file:///C:/Users/tatyana/Documents/Projects/Advance project/Advance/OutputFiles/TGG_metadata_negatives.txt", header=TRUE)

NGTX_sampleID$TGG_ID<- NULL
NGTX_sampleID$Casrn<- NULL
NGTX_sampleID$X_id_<- NULL

TG_Gates_rat_log2fold_NGTX$assayId <- NULL
TG_Gates_rat_log2fold_NGTX$pvalue <- NULL
TG_Gates_rat_log2fold_NGTX$valueType <- NULL



NGTX_gene_search<-merge(TG_Gates_rat_log2fold_NGTX,NGTX_sampleID, by.x="sampleid",by.y="sampleid")

NGTX_gene_search$sampleId <- NULL
NGTX_gene_search$cellType <- NULL
NGTX_gene_search$repeatType <- NULL
NGTX_gene_search$tissue <- NULL
NGTX_gene_search$organism <- NULL

GTX_gene_search$sampleId <- NULL
GTX_gene_search$cellType <- NULL
GTX_gene_search$repeatType <- NULL
GTX_gene_search$tissue <- NULL
GTX_gene_search$organism <- NULL

#add new column
GTX_gene_search$group<-"PC"
NGTX_gene_search$group<-"NC"

write.csv(GTX_gene_search, file="TOTAL_PC.csv")
write.csv(NGTX_gene_search, file="TOTAL_NC.csv")


TOTAL<-rbind(GTX_gene_search, NGTX_gene_search)

colnames(TOTAL)[4]<-"FC"

write.csv(TOTAL, file="TOTAL_hepatitis.csv")



