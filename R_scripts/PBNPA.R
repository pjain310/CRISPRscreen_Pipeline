#!/usr/local/bin/R

#Load required libraries
library(PBNPA)

#Read input arguments (counts table and sample map)
args = commandArgs(trailingOnly=TRUE)

#Process sample map according to PBNPA requirements
data<-read.table(args[1],header=TRUE)
design_mat<-read.table(args[2],header=TRUE)

#get condition name for file naming
condition <- strsplit(args[2],"__|\\.")

dat<-list()
control<-list()
treatment<-list()


for (row in 1:nrow(design_mat)){
  Trt<-as.character(design_mat[row,"Treatment"])
  Ctl<-as.character(design_mat[row,"Control"])
  varname <- paste("R",row,sep="")

  assign(varname,data.frame(data[,1:2],data[,Ctl],data[,Trt]))
  temp<-get(varname)
  dat[[row]] = temp
  treatment[[row]] = as.character(Trt)
  control[[row]] = as.character(Ctl)
  names(dat[[row]]) <- c("sgRNA","Gene",control[[row]],treatment[[row]])
}

results<-PBNPA(dat)

write.table(results$pos.gene,file=paste("temp/pbnpa/pos_gene_pbnpa_",condition[[1]][2],".txt",sep=""),quote = FALSE)
write.table(results$neg.gene,file=paste("temp/pbnpa/neg_gene_pbnpa_",condition[[1]][2],".txt",sep=""),quote = FALSE)
write.table(results$final.result,file=paste("temp/pbnpa/pbnpa_",condition[[1]][2],".txt",sep=""),quote = FALSE)
