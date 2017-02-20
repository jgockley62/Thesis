#!/usr/bin/env Rscript

#Written by Jake Gockley, Yale University, 2017
#This script loops through MPRA Fragments and calls one-sided Wilcoxon rank scores

##Example Call: Example Call: R CMD BATCH '--args Lin1_minP_filt.txt Lin1 Lin1_Index.txt 20' MPRA_Activity_Caller.R
#CMD BATCH '--args <Tag_Activity_File.txt> <Output File Desination> <Fragment_Index_File.txt> <Number of Cores to use>' MPRA_Activity_Caller.R

library(readr)
library(doSNOW)
library(foreach)

#Load Ooutput File Designation
args <- commandArgs(trailingOnly = TRUE)
Input <- paste0(args[1])
fileN <- paste0(args[2])
fileI <- paste0(args[3])
cores <- as.numeric(paste0(args[4]))
print(Input)
print(fileN)
print(fileI)
print(cores)

###Run This
Lin1_minP.filt <- data.frame(read_delim(
  file = paste0(Input),
  delim = "\t",
  col_names = T
))
Lin1_Index <- data.frame(read_delim(
  file = paste0(fileI),
  delim = "\t",
  col_names = T
))

Indexr <- function(n){
  Dummy<-data.frame(matrix(0,0,4))
  for(i in n){
    
    Dummy<-data.frame(matrix(0,0,4))
    SigVal <- wilcox.test(
      Lin1_minP.filt[Lin1_Index[i,2]:Lin1_Index[i,3], ]$value,
      Lin1_minP.filt[-( Lin1_Index[i,2]:Lin1_Index[i,3] ), ]$value,
      alternative = c("greater")
    )$p.value
    
    Dummy[1,]<-c(Lin1_Index[i,1],Lin1_Index[i,2],Lin1_Index[i,3],SigVal)
  }
  return(Dummy)
}

#Build the Cluster of N Cores
cluster = makeCluster(cores, type = "SOCK")
registerDoSNOW(cluster)
clusterExport(cluster, c("Lin1_minP.filt", "Indexr","Lin1_Index"))

#Run
Test = clusterApply( cluster, 1:dim(Lin1_Index)[1], Indexr)

#Merge and Save Index File
Index <-do.call("rbind", Test)
colnames(Index) <- c("Fragment", "RowStart", "RowEnd", "PVal")
write.table(Index, paste0(fileN,"_Wilcoxon.txt"), quote=FALSE, sep="\t", row.names=F, col.names=T )
stopCluster(cluster)
q(save = "no")
