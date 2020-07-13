#if (!require(ggplot2)) {
# install.packages("ggplot2")}
library(ggplot2)


#if (!require(plyr)) {
#install.packages("plyr")}
library(plyr)


#if (!require(dplyr)) {
# install.packages("dplyr")
library(dplyr)


#if (!require(tidyr)) {install.packages("tidyr")}
library(tidyr)


#if (!require(tidyverse)) {install.packages("tidyverse")}
library(tidyverse)


#if (!require(reshape2)) { install.packages("reshape2")}
library(reshape2)

TOTAL_all<-read.csv("file:///C:/Users/tatyana/Documents/projects/KIT_phospholipidosis/TOTAL_hepatitis.csv",sep=",", header=TRUE)
TOTAL<- filter (TOTAL_all, geneSymbols== "Hbb")

TOTAL<- data.frame(TOTAL)

TOTAL<-ddply(TOTAL, .(TGG_compoundName,timepointHr,group,doseLevel), summarize, FC=mean(FC))

TOTAL$hour <- "h"

exp_time <- c("4", "8","15","29")


TOTAL$exposure_time_factor = factor(TOTAL$timepointHr, levels=c( "4","8","15", "29"))


# estimate limits
min <- min(TOTAL[TOTAL$timepointHr%in% exp_time , ]$FC)
print (min)
max <- max(TOTAL[TOTAL$timepointHr %in% exp_time , ]$FC)
print (max)
#Center the title

#Create a plot

ggplot(TOTAL[TOTAL$timepointHr %in% exp_time , ], aes(y=FC, x=TGG_compoundName)) +
  scale_y_continuous(name="log2 Hbb", limits=c(-2,2), breaks=c(-2,0,2), labels=c("-2", "0", "2")) +
  geom_bar(stat="identity", colour="white", aes(fill=group)) + 
  coord_flip() +
  facet_grid(rows=vars(group), cols=vars(exposure_time_factor), scale="free_y") +
  scale_fill_brewer(palette="Set1") +
  theme(axis.title.y=element_blank())

ggsave("Hbb.tiff")
