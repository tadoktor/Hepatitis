####################################################################
####Author: Noffisat Oki, Tatyana Doktorova
####Version: 1.0 08-28-2018
####Description: GTX Subset script:
####             The script uses the information extracted from TG Gates and 
###              filtering of interesting genes according to fold change and p-value.
####             Subsets chemical-gene data to filter for only genes active
####             in at least 3 chemicals 
###              Reannotation to Ensemb Identifier and visualizations
####             
####Notes: The path to the 'ADVANCE_Project' project folder should be modified to user specified 
####       locations for proper execution of the script. The 'Root_location' object on line 20
####       should be used for this purpose.
####Potential issues:

##The root location folder should be specified here. 
#This is the main folder where any subfolders from which data will be read from or written to are located
Root_location <- "C:\\Users\\tatyana\\Documents\\Projects\\Advance project\\Advance"

##Output directory
outDir <- paste(Root_location, "/IntermediateData/", sep='')

# Install libraries
library(data.table)
library(reshape2)
library("biomaRt")
library("dplyr")
library(collapsibleTree)


#Generation of Disease list of genes for rat in vivo
#The files can be derived either from the previous Step_3B1, Step_3B2, Step3_B3 or directly downloaded from this link: https://mega.nz/#F!zDBmiKyZ

TG_Gates_rat_log2fold_GTX<- read.table("file:///C:/Users/tatyana/Documents/Projects/Advance project/Advance/OutputFiles/TGG_foldchanges_hepatitis.txt", header=TRUE)
GTX_sampleID<- read.table("file:///C:/Users/tatyana/Documents/Projects/Advance project/Advance/OutputFiles/TGG_metadata_hepatitis.txt", header=TRUE)
HCC_Toxcast_CTD<-read.csv("file:///C:/Users/tatyana/Documents/Projects/KIT_phospholipidosis/Hepatitis/HCC_Toxcast_CTD_hep.csv", header=TRUE)
Unique_genes = HCC_Toxcast_CTD %>% distinct(lhs)
GTX_sampleID_filter<-filter(GTX_sampleID, timepointHr== "29")

GTX_all<-merge(TG_Gates_rat_log2fold_GTX,GTX_sampleID_filter, by.x="sampleid",by.y="sampleid")

#Filtering according to p value and fold change


GTX_filtering<-GTX_all[!is.na(GTX_all$pvalue),]

GTX_pval0.05<- GTX_filtering[GTX_filtering$pvalue<0.05,]
GTX_pval0.05_FC_up_1<-GTX_pval0.05[GTX_pval0.05$ value > 1,]
GTX_pval0.05_FC_down_1<-GTX_pval0.05[GTX_pval0.05$ value < -1,]
GTX_specific<- rbind(GTX_pval0.05_FC_down_1,GTX_pval0.05_FC_up_1)
#write.csv(GTX_significant, file="GTX_specific_rat.csv")

#Generation of NGTX list of genes
TG_Gates_rat_log2fold_NGTX<- read.table("file:///C:/Users/tatyana/Documents/Projects/Advance project/Advance/OutputFiles/TGG_foldchanges_negatives.txt", header=TRUE)
NGTX_sampleID<- read.table("file:///C:/Users/tatyana/Documents/Projects/Advance project/Advance/OutputFiles/TGG_metadata_negatives.txt", header=TRUE)

NGTX_sampleID_filter<-filter(NGTX_sampleID, timepointHr== "29")

NGTX_all<-merge(TG_Gates_rat_log2fold_NGTX,NGTX_sampleID_filter, by.x="sampleid",by.y="sampleid")

NGTX_filtering<-NGTX_all[!is.na(NGTX_all$pvalue),]

NGTX_pval0.05<- NGTX_filtering[NGTX_filtering$pvalue<0.05,]
NGTX_pval0.05_FC_up_1<-NGTX_pval0.05[NGTX_pval0.05$ value > 1,]
NGTX_pval0.05_FC_down_1<-NGTX_pval0.05[NGTX_pval0.05$ value < -1,]
NGTX_specific<- rbind(NGTX_pval0.05_FC_down_1,NGTX_pval0.05_FC_up_1)
#write.csv(NGTX_significant, file="NGTX_specific.csv")





GTX_chemgene_data.long <- as.data.frame (GTX_specific)
##reformatting the data
GTX_chemgene_data.wide <- dcast(GTX_chemgene_data.long, TGG_compoundName ~ assayid , value.var="geneSymbols")
recoding_GTX_chemgene_data.wide <- (GTX_chemgene_data.wide)
recoding_GTX_chemgene_data.wide[,2:ncol(recoding_GTX_chemgene_data.wide)][recoding_GTX_chemgene_data.wide[,2:ncol(recoding_GTX_chemgene_data.wide)]>1] <- 1
rownames(recoding_GTX_chemgene_data.wide) <- recoding_GTX_chemgene_data.wide[,1]
recoding_GTX_chemgene_data.wide <- ((recoding_GTX_chemgene_data.wide[,-1]))
#write.table(recoding_GTX_chemgene_data.wide,file=file.path(outDir, "GTX_chem_gene_data_wide.csv"), sep=',', col.names=TRUE, row.names=TRUE,quote=FALSE)

Active_GTX_dataset<- recoding_GTX_chemgene_data.wide
Active_GTX_dataset[is.na(Active_GTX_dataset)] <- 0

Active_GTX_dataset<-data.matrix(Active_GTX_dataset)

#remove columns which are all '0' indicating no effect for an assay, gene or disease
Active_GTX_dataset<- Active_GTX_dataset[,!apply((Active_GTX_dataset),2,function(x) sum(abs(x), na.rm=TRUE) < 2)]
chemicals <- as.data.frame(rownames(Active_GTX_dataset))
colnames(chemicals) <- "chemical"
Probes <- as.data.frame(colnames(Active_GTX_dataset))
Active_GTX_dataset <- cbind(chemicals, Active_GTX_dataset)
Active_GTX_dataset_T<-t(Active_GTX_dataset)
Active_GTX_dataset<-cbind(Row.Names=rownames(Active_GTX_dataset_T),Active_GTX_dataset_T)

GTX_GL_final <-Probes
colnames(GTX_GL_final) <- "assayid"

#Selecting the more than 3 chemicals deregulated genes for the Negative group

NGTX_chemgene_data.long <- as.data.frame (NGTX_specific)
##reformatting the data
NGTX_chemgene_data.wide <- dcast(NGTX_chemgene_data.long, TGG_compoundName ~ assayid , value.var="geneSymbols")
recoding_NGTX_chemgene_data.wide <- (NGTX_chemgene_data.wide)
recoding_NGTX_chemgene_data.wide[,2:ncol(recoding_NGTX_chemgene_data.wide)][recoding_NGTX_chemgene_data.wide[,2:ncol(recoding_NGTX_chemgene_data.wide)]>1] <- 1
rownames(recoding_NGTX_chemgene_data.wide) <- recoding_NGTX_chemgene_data.wide[,1]
recoding_NGTX_chemgene_data.wide <- ((recoding_NGTX_chemgene_data.wide[,-1]))
#write.table(recoding_GTX_chemgene_data.wide,file=file.path(outDir, "GTX_chem_gene_data_wide.csv"), sep=',', col.names=TRUE, row.names=TRUE,quote=FALSE)

Active_NGTX_dataset<- recoding_NGTX_chemgene_data.wide

Active_NGTX_dataset[is.na(Active_NGTX_dataset)] <- 0

Active_NGTX_dataset<-data.matrix(Active_NGTX_dataset)

#remove columns which are all '0' indicating no effect for an assay, gene or disease
Active_NGTX_dataset<- Active_NGTX_dataset[,!apply((Active_NGTX_dataset),2,function(x) sum(abs(x), na.rm=TRUE) < 2)]
chemicals <- as.data.frame(rownames(Active_NGTX_dataset))
colnames(chemicals) <- "chemical"
NGTX_Probes <- as.data.frame(colnames(Active_NGTX_dataset))
Active_NGTX_dataset <- cbind(chemicals, Active_NGTX_dataset)
Active_NGTX_dataset_T<-t(Active_NGTX_dataset)
Active_NGTX_dataset<-cbind(Row.Names=rownames(Active_NGTX_dataset_T),Active_NGTX_dataset_T)

NGTX_GL_final <-NGTX_Probes
colnames(NGTX_GL_final) <- "assayid"

#Attach the gene symbols

GTX_GL_final<-merge(GTX_GL_final,GTX_specific)
GTX_GL_genes<- unique(subset(GTX_GL_final, select=c("geneSymbols")))

NGTX_GL_final<-merge(NGTX_GL_final,NGTX_specific)
NGTX_GL_genes<- unique(subset(NGTX_GL_final, select=c("geneSymbols")))




##Common genees between GTX, NGTX
Genes_common<-intersect(GTX_GL_genes, NGTX_GL_genes)

## Toxic group specific genes
TC_specific<-setdiff(GTX_GL_genes, NGTX_GL_genes)

##write.csv(NGTX_specific, file="NGTX_specific_rat_15Days_Middle.csv")

##write.csv(Unique_genes, file="Unique_genes_rat_15Days_Middle.csv")

##Overlap between the ToxCast_CTD HCC-specific genes and TG Gates
TC_specific<-as.data.frame(apply(TC_specific,2,toupper))
names(TC_specific)<-"lhs"

# Overlapping genes
Genes_TGGATEs_CTD_Toxcast<-as.data.frame(intersect(TC_specific$lhs, Unique_genes$lhs))

names(Genes_TGGATEs_CTD_Toxcast)<-"lhs"

                                              
##Overlap between the ToxCast_CTD HCC-specific genes and TG Gates


Disease<-merge(Genes_TGGATEs_CTD_Toxcast, HCC_Toxcast_CTD, by.x="lhs", by.y="lhs")

#Visualization tree
collapsibleTree(df=Disease, c("lhs", "Major"),  fill = "transparent")

#adding the fold changes to the genes
GTX_all<-as.data.frame(apply(GTX_all,2,toupper))
NGTX_all<-as.data.frame(apply(NGTX_all,2,toupper))

Heatmap_NGTX_Middle_GTX<-merge (Genes_TGGATEs_CTD_Toxcast,GTX_all, by.x="lhs", by.y="geneSymbols")

Heatmap_NGTX_Middle_NGTX<-merge (Genes_TGGATEs_CTD_Toxcast,NGTX_all, by.x="lhs", by.y="geneSymbols")


HEATMAP<- rbind(Heatmap_NGTX_Middle_GTX,Heatmap_NGTX_Middle_NGTX) 


#Create a mean value for all genes with several identifiers

HEATMAP$name<-paste (HEATMAP$lhs, HEATMAP$TGG_compoundName, sep=" ")
HCA_merge<-HEATMAP[,c(4)]
HCA_merge<-apply(as.matrix.noquote (HCA_merge),2,as.numeric)
HCA_1<-aggregate (HCA_merge, by=list(HEATMAP$name), FUN=mean,na.rm=TRUE)

HCA_1$Gene = as.character(lapply(strsplit(as.character(HCA_1$Group.1), split=" "), "[", 1))
HCA_1$Compound = as.character(lapply(strsplit(as.character(HCA_1$Group.1), split=" "), "[", 2))
HCA_1<-HCA_1[,c(3,4,2)]

library(reshape)
HCA<-cast(HCA_1, Gene ~ Compound)

write.csv(HCA, file="HCA_29DAYS_High.csv")


m = data.matrix( HCA[,-1])
##colnames(m) <- c("ACB", "ACE","ACF","AMI","ASP",
##               "CAR","CHP","CIS","COL","CYC",
##           "DIC","DOX","ETH","ETY","GRI",
##          "HEX","IBU","MEL","MET","NAP",
##               "NIF","NIT","PHE","ROS","TAN","THI","WY1")
rownames(m) = as.character( HCA[,1])
m[m == 0] <- NA

library(gplots)
colors = c(seq(-3,-2,length=100),seq(-1,0.5,length=100),seq(0.4,6,length=100))
my_palette <- colorRampPalette(c("green", "gray", "red")) (n=100000)

heatmap.2 (m, trace="none", na.color = "gray", scale="none", breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5),
           col = greenred(10), cexRow=1, cexCol = 1, ColSideColors= c("red","red","darkgreen","darkgreen","darkgreen","darkgreen","red","red","darkgreen",
                                                                      "red","red","darkgreen","red",
                                                                      "red","red","darkgreen","red","darkgreen","red",
                                                                      "darkgreen"))

