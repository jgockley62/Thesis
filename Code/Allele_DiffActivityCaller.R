#!/usr/bin/env Rscript

#Written by Jake Gockley, Yale University, 2017
#This script loops through MPRA Fragment ortholog pairs and calls wilcoxon two tailed differential activity

##Example Call: Example Call: R CMD BATCH '--args Lin1_minP_filt.txt Lin2_minP_filt.txt Lin1_Index.txt Lin2_Index.txt 20' Allele_DiffActivityCaller.R
#CMD BATCH '--args <Tag_Activity_FileLin1.txt> <Tag_Activity_FileLin2.txt> <Fragment_Index_FileLin1.txt> <Fragment_Index_FileLin2.txt> <Number of Cores to use>' MPRA_Activity_Caller.R
library(exactRankTests)
library(readr)
library(doSNOW)
library(foreach)

#Load Lineage Data File Designation
args <- commandArgs(trailingOnly = TRUE)
Inputa <- paste0(args[1])
Inputb <- paste0(args[2])
#Load Index File Designation
fileNa <- paste0(args[3])
fileNb <- paste0(args[4])
cores <- as.numeric(paste0(args[5]))

###Load in files
#Inputa <- "Lin1_minP_filt.txt"
Lin1_minP.filt <- data.frame(read_delim(
  file = paste0(Inputa),
  delim = "\t",
  col_names = T
))
#Inputb <- "Lin2_minP_filt.txt"
Lin2_minP.filt <- data.frame(read_delim(
  file = paste0(Inputb),
  delim = "\t",
  col_names = T
))
#fileNa<-"Lin1_Index.txt"
Lin1_Index <- data.frame(read_delim(
  file = paste0(fileNa),
  delim = "\t",
  col_names = T
))
#fileNb<-"Lin2_Index.txt"
Lin2_Index <- data.frame(read_delim(
  file = paste0(fileNb),
  delim = "\t",
  col_names = T
))
rownames(Lin1_Index) <- as.character(Lin1_Index[,1])
rownames(Lin2_Index) <- as.character(Lin2_Index[,1])

#Make Testable Allele Pairs from 
system(paste0("grep \'hg19\' ",Inputa," > Human_temp.txt"))
system("sort -k 3,3 Human_temp.txt > Human_srt_temp.txt") 
system("awk \'{ print $2\"\t\"$3 }\' Human_srt_temp.txt > Human_frag_srt_temp.txt")
system("uniq Human_frag_srt_temp.txt > Unique_Human_Ents_temp.txt")

system(paste0("grep \'PanTro2\' ",Inputa," > Chimp_temp.txt"))
system("sort -k 3,3 Chimp_temp.txt > Chimp_srt_temp.txt") 
system("awk \'{ print $3\"\t\"$2 }\' Chimp_srt_temp.txt > Chimp_frag_srt_temp.txt")
system("uniq Chimp_frag_srt_temp.txt > Unique_Chimp_Ents_temp.txt")

system("cat Unique_Human_Ents_temp.txt Unique_Chimp_Ents_temp.txt > Total_temp.txt")
system("sort -k 1,1 Total_temp.txt > Total_srt_temp.txt") 
system("uniq Total_srt_temp.txt > Testable_Alleles.txt")
system("rm *_temp.txt")

Alleles <- data.frame(read_delim(
  file = paste0("Testable_Alleles.txt"),
  delim = "\t",
  col_names = F
))
colnames(Alleles)<-c("Human", "Chimp")

AlleleCaller <- function(n){
  Dummy<-data.frame(matrix(0,0,8))
  for(i in n){
    #Test Lin1
    if(Alleles[i,2] %in% rownames(Lin1_Index) & Alleles[i,1] %in% rownames(Lin1_Index) ){
      Lin1Pval <- wilcox.exact(
        Lin1_minP.filt[Lin1_Index[Alleles[i,1],2]:Lin1_Index[Alleles[i,1],3], ]$value,
        Lin1_minP.filt[Lin1_Index[Alleles[i,2],2]:Lin1_Index[Alleles[i,2],3], ]$value
      )$p.value
    
      Lin1Pval <- Lin1Pval
      Lin1Status <- "Paired"
      
    }else{
      Lin1Pval <- NA
      Lin1Status <- "UnPaired"
    }
    if(Alleles[i,2] %in% rownames(Lin2_Index) & Alleles[i,1] %in% rownames(Lin2_Index) ){
      Lin2Pval <- wilcox.exact(
        Lin2_minP.filt[Lin2_Index[Alleles[i,1],2]:Lin2_Index[Alleles[i,1],3], ]$value,
        Lin2_minP.filt[Lin2_Index[Alleles[i,2],2]:Lin2_Index[Alleles[i,2],3], ]$value
      )$p.value
      
      Lin2Pval <- Lin2Pval
      Lin2Status <- "Paired"
      
    }else{
      Lin2Pval <- NA
      Lin2Status <- "UnPaired"
    }
    Dummy[1,]<-c(Alleles[i,1],Alleles[i,2],Lin1Pval,Lin1Status,"BLANK",Lin2Pval,Lin2Status,"BLANK")
    Dummy[2,]<-c(Alleles[i,2],Alleles[i,1],Lin1Pval,Lin1Status,"BLANK",Lin2Pval,Lin2Status,"BLANK")
  }
  return(Dummy)
}

#Build the Cluster of N Cores
cluster = makeCluster(cores, type = "SOCK")
registerDoSNOW(cluster)
clusterExport(cluster, c("Lin1_minP.filt","Lin2_minP.filt","Lin1_Index","Lin2_Index","Alleles","AlleleCaller"))
clusterEvalQ(cluster, library(exactRankTests))

#Run
Test = clusterApply( cluster, 1:dim(Alleles)[1], AlleleCaller)

#Aggregate and print to file
Index <-do.call("rbind", Test)
colnames(Index) <- c("Fragment", "Ortholog", "Lin1_Allele_pVal", "Lin1_Allele_Status","Lin1_Dif_Act","Lin2_Allele_pVal", "Lin2_Allele_Status","Lin2_Dif_Act")
write.table(Index, paste0("MPRA_AlleleActivity_Wilcoxon.txt"), quote=FALSE, sep="\t", row.names=F, col.names=T )
stopCluster(cluster)
q(save = "no")




