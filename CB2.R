#!/usr/local/bin/R

#Load required libraries
library("CB2")

#Read input arguments (counts table and sample map)
args = commandArgs(trailingOnly=TRUE)

#Process sample map according to PBNPA requirements
data<-read.table(args[1],header=TRUE)
design_mat<-read.table(args[2],header=TRUE)

rownames(data)<-data[,1]
data<-data[,-c(1:2)]

lvl<-levels(design_mat$group)

results<-run_estimation(data,design_mat,lvl[1],lvl[2])

write.table(results,file=paste("temp/cb2/results_",lvl[1],"_vs_",lvl[2],".txt",sep=""),quote = FALSE)
