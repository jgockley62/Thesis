#!/usr/bin/env Rscript

#Written by Jake Gockley, Yale University, 2017
#This script loops through MPRA Fragments and preps the Index Files

##Example Call: R CMD BATCH '--args Lin1_minP_filt.txt Lin1 20' MPRA_Activity_Indexer.R
#CMD BATCH '--args <Tag_Activity_File.txt> <Output File Desination> <Number of Cores to use>' MPRA_Activity_Indexer.R

library(readr)
library(doSNOW)
library(foreach)

#Load Ooutput File Designation
args <- commandArgs(trailingOnly = TRUE)
Input <- paste0(args[1])
fileN <- paste0(args[2])
cores <- as.numeric(paste0(args[3]))
print(Input)
print(fileN)
print(cores)

###Run This
Lin1_minP.filt <- data.frame(read_delim(
  file = paste0(Input),
  delim = "\t",
  col_names = T
))

foo <- levels(factor( Lin1_minP.filt$Alignment))
Temp <- data.frame(Lin1_minP.filt[,1:2])

Findr <- function(n){
  Dummy<-data.frame(matrix(0,0,4))
  for(i in n) {
    select <- as.numeric(rownames(subset(Temp, Alignment == foo[i] ) ) ) 
    Dummy[1,]<-c(foo[i],select[1],select[length(select)],0)
  }
  return(Dummy)
}

#Build the Cluster of N Cores
cluster = makeCluster(cores, type = "SOCK")
registerDoSNOW(cluster)
clusterExport(cluster, c("Temp", "Findr","foo"))

#Run 
Test = clusterApply( cluster, 1:length(foo), Findr)

#Merge and Save Index File
Index <-do.call("rbind", Test)
colnames(Index) <- c("Fragment", "RowStart", "RowEnd", "PVal")
write.table(Index, paste0(fileN,"_Index.txt"), quote=FALSE, sep="\t", row.names=F, col.names=T )
stopCluster(cluster)
q(save = "no")
