#!/usr/bin/env Rscript

#Written by Jake Gockley, Yale University, 2017
#This script identifies median fold activities of fragment tag distribution and quantile normalizes the polulation 
#It outputs the correlation of PermA and PermB quantile normalized median fragment scores
#This is utalized by rand_sampl.py and must be in the called run directory of that script

#Example Call: R CMD BATCH '--args <Number of Tags per-fragment>'

library(corrplot)
library(Hmisc)
library(readr)
library(plyr)
library(reshape2)
library(devtools)
library(Biobase)
library(preprocessCore)

args <- commandArgs(trailingOnly = TRUE)
fileN <- paste(args[1],"_Barcodes.txt",sep="")
print(fileN)
print(args)
wd<-getwd()
FOO_culled_data_A <- read.table(file=paste0(wd,"/Permute_A.txt"), sep='\t', header = T)
FOO_culled_data_B <- read.table(file=paste0(wd,"/Permute_B.txt"), sep='\t', header = T)

culled_data2 <- FOO_culled_data_B
culled_data <- FOO_culled_data_A

#Process data normally for subsample 1
temp <- normalize.quantiles(as.matrix(culled_data[, 44:50]))
temp <- temp - median(temp)
culled_data.melt <-
  melt(
    culled_data[c(
      "Tag",
      "Alignment",
      "Ortholog_Seq",
      "Chr",
      "Start",
      "Stop",
      "Species",
      "eVar_BPChange",
      "eVar_Pos",
      "Activity_Rep1_2",
      "Activity_Rep1_3",
      "Activity_Rep2_2",
      "Activity_Rep2_3",
      "Activity_Lin1",
      "Activity_Lin2",
      "Activity_Total"
    )],
    id = c(
      "Tag",
      "Alignment",
      "Ortholog_Seq",
      "Chr",
      "Start",
      "Stop",
      "Species",
      "eVar_BPChange",
      "eVar_Pos"
    )
  )

#Setup final data frame
culled_data.melt$value <- melt(temp[, 1:7])$value
RATIO.temp <- culled_data.melt[!duplicated(culled_data.melt$Alignment),!(names(culled_data.melt) %in% c("value", "variable"))]

#Take the median value of Rep 1.2 permuted activity by fragment
Activity_Rep1_2.temp <- culled_data.melt[grep("Activity_Rep1_2", culled_data.melt$variable),]
Activity_Rep1_2.value <- tapply(Activity_Rep1_2.temp$value, factor(Activity_Rep1_2.temp$Alignment), median)
#Aggreagte data back to data frame
culled_data.ratio <- merge(RATIO.temp, data.frame(Alignment = names(Activity_Rep1_2.value), Activity_Rep1_2.median = Activity_Rep1_2.value), by = "Alignment")

#Take the median value of Rep 1.3 permuted activity by fragment
Activity_Rep1_3.temp <- culled_data.melt[grep("Activity_Rep1_3", culled_data.melt$variable),]
Activity_Rep1_3.value <- tapply(Activity_Rep1_3.temp$value, factor(Activity_Rep1_3.temp$Alignment), median)
#Aggreagte data back to data frame
culled_data.ratio <- merge(culled_data.ratio, data.frame(Alignment = names(Activity_Rep1_3.value), Activity_Rep1_3.median = Activity_Rep1_3.value), by = "Alignment")

#Take the median value of Rep 2.2 permuted activity by fragment
Activity_Rep2_2.temp <- culled_data.melt[grep("Activity_Rep2_2", culled_data.melt$variable),]
Activity_Rep2_2.value <- tapply(Activity_Rep2_2.temp$value, factor(Activity_Rep2_2.temp$Alignment), median)
#Aggreagte data back to data frame
culled_data.ratio <- merge(culled_data.ratio, data.frame(Alignment = names(Activity_Rep2_2.value), Activity_Rep2_2.median = Activity_Rep2_2.value), by = "Alignment")

#Take the median value Rep 2.3 permuted activity by fragment
Activity_Rep2_3.temp <- culled_data.melt[grep("Activity_Rep2_3", culled_data.melt$variable),]
Activity_Rep2_3.value <- tapply(Activity_Rep2_3.temp$value, factor(Activity_Rep2_3.temp$Alignment), median)
#Aggreagte data back to data frame
culled_data.ratio <- merge(culled_data.ratio, data.frame(Alignment = names(Activity_Rep2_3.value), Activity_Rep2_3.median = Activity_Rep2_3.value), by = "Alignment")

#Take the median value of Lin1 permuted activity by fragment
Activity_Lin1.temp <- culled_data.melt[grep("Activity_Lin1", culled_data.melt$variable),]
Activity_Lin1.value <- tapply(Activity_Lin1.temp$value, factor(Activity_Lin1.temp$Alignment), median)
#Aggreagte data back to data frame
culled_data.ratio <- merge(culled_data.ratio, data.frame(Alignment = names(Activity_Lin1.value), Activity_Lin1.median = Activity_Lin1.value), by = "Alignment")

#Take the median value of Lin2 permuted activity by fragment
Activity_Lin2.temp <- culled_data.melt[grep("Activity_Lin2", culled_data.melt$variable),]
Activity_Lin2.value <- tapply(Activity_Lin2.temp$value, factor(Activity_Lin2.temp$Alignment), median)
#Aggreagte data back to data frame
culled_data.ratio <- merge(culled_data.ratio, data.frame(Alignment = names(Activity_Lin2.value), Activity_Lin2.median = Activity_Lin2.value), by = "Alignment")

#Take the median value of Lin2 permuted activity by fragment
Activity_Total.temp <- culled_data.melt[grep("Activity_Total", culled_data.melt$variable),]
Activity_Total.value <- tapply(Activity_Total.temp$value, factor(Activity_Total.temp$Alignment), median)
#Aggreagte data back to data frame
culled_data.ratio <- merge(culled_data.ratio, data.frame(Alignment = names(Activity_Total.value), Activity_Total.median = Activity_Total.value), by = "Alignment")

#Process data normally for subsample 2
temp <- normalize.quantiles(as.matrix(culled_data2[, 44:50]))
temp <- temp - median(temp)
culled_data2.melt <-
  melt(
    culled_data2[c(
      "Tag",
      "Alignment",
      "Ortholog_Seq",
      "Chr",
      "Start",
      "Stop",
      "Species",
      "eVar_BPChange",
      "eVar_Pos",
      "Activity_Rep1_2",
      "Activity_Rep1_3",
      "Activity_Rep2_2",
      "Activity_Rep2_3",
      "Activity_Lin1",
      "Activity_Lin2",
      "Activity_Total"
    )],
    id = c(
      "Tag",
      "Alignment",
      "Ortholog_Seq",
      "Chr",
      "Start",
      "Stop",
      "Species",
      "eVar_BPChange",
      "eVar_Pos"
    )
  )


#Setup final data frame
culled_data2.melt$value <- melt(temp[, 1:7])$value
RATIO.temp <- culled_data2.melt[!duplicated(culled_data2.melt$Alignment),!(names(culled_data2.melt) %in% c("value", "variable"))]

#Take the median value of Rep 1.2 permuted activity by fragment
Activity_Rep1_2.temp <- culled_data2.melt[grep("Activity_Rep1_2", culled_data2.melt$variable),]
Activity_Rep1_2.value <- tapply(Activity_Rep1_2.temp$value, factor(Activity_Rep1_2.temp$Alignment), median)
#Aggreagte data back to data frame
culled_data2.ratio <- merge(RATIO.temp, data.frame(Alignment = names(Activity_Rep1_2.value), Activity_Rep1_2.median = Activity_Rep1_2.value), by = "Alignment")

#Take the median value of Rep 1.3 permuted activity by fragment
Activity_Rep1_3.temp <- culled_data2.melt[grep("Activity_Rep1_3", culled_data2.melt$variable),]
Activity_Rep1_3.value <- tapply(Activity_Rep1_3.temp$value, factor(Activity_Rep1_3.temp$Alignment), median)
#Aggreagte data back to data frame
culled_data2.ratio <- merge(culled_data2.ratio, data.frame(Alignment = names(Activity_Rep1_3.value), Activity_Rep1_3.median = Activity_Rep1_3.value), by = "Alignment")

#Take the median value of Rep 2.2 permuted activity by fragment
Activity_Rep2_2.temp <- culled_data2.melt[grep("Activity_Rep2_2", culled_data2.melt$variable),]
Activity_Rep2_2.value <- tapply(Activity_Rep2_2.temp$value, factor(Activity_Rep2_2.temp$Alignment), median)
#Aggreagte data back to data frame
culled_data2.ratio <- merge(culled_data2.ratio, data.frame(Alignment = names(Activity_Rep2_2.value), Activity_Rep2_2.median = Activity_Rep2_2.value), by = "Alignment")

#Take the median value Rep 2.3 permuted activity by fragment
Activity_Rep2_3.temp <- culled_data2.melt[grep("Activity_Rep2_3", culled_data2.melt$variable),]
Activity_Rep2_3.value <- tapply(Activity_Rep2_3.temp$value, factor(Activity_Rep2_3.temp$Alignment), median)
#Aggreagte data back to data frame
culled_data2.ratio <- merge(culled_data2.ratio, data.frame(Alignment = names(Activity_Rep2_3.value), Activity_Rep2_3.median = Activity_Rep2_3.value), by = "Alignment")

#Take the median value of Lin1 permuted activity by fragment
Activity_Lin1.temp <- culled_data2.melt[grep("Activity_Lin1", culled_data2.melt$variable),]
Activity_Lin1.value <- tapply(Activity_Lin1.temp$value, factor(Activity_Lin1.temp$Alignment), median)
#Aggreagte data back to data frame
culled_data2.ratio <- merge(culled_data2.ratio, data.frame(Alignment = names(Activity_Lin1.value), Activity_Lin1.median = Activity_Lin1.value), by = "Alignment")

#Take the median value of Lin2 permuted activity by fragment
Activity_Lin2.temp <- culled_data2.melt[grep("Activity_Lin2", culled_data2.melt$variable),]
Activity_Lin2.value <- tapply(Activity_Lin2.temp$value, factor(Activity_Lin2.temp$Alignment), median)
#Aggreagte data back to data frame
culled_data2.ratio <- merge(culled_data2.ratio, data.frame(Alignment = names(Activity_Lin2.value), Activity_Lin2.median = Activity_Lin2.value), by = "Alignment")

#Take the median value of Lin2 permuted activity by fragment
Activity_Total.temp <- culled_data2.melt[grep("Activity_Total", culled_data2.melt$variable),]
Activity_Total.value <- tapply(Activity_Total.temp$value, factor(Activity_Total.temp$Alignment), median)
#Aggreagte data back to data frame
culled_data2.ratio <- merge(culled_data2.ratio, data.frame(Alignment = names(Activity_Total.value), Activity_Total.median = Activity_Total.value), by = "Alignment")

COR <- sapply(10:ncol(culled_data.ratio), function(i) cor(culled_data.ratio[,i], culled_data2.ratio[,i]))
write.table(data.frame(t(unlist(COR))), file = fileN, col.names = FALSE, row.names = FALSE, quote = FALSE)

q(save = "no")