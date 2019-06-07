#!/usr/bin/env RScript

#Load required libraries
library(PBNPA)

#Read input arguments (counts table and sample map)
args = commandArgs(trailingOnly=TRUE)

#Process sample map according to PBNPA requirements
data<-read.table(args[1],header=TRUE)
design_mat<-read.table(args[2],header=TRUE)


dat<-list()
control<-list()
treatment<-list()


for (row in 1:nrow(design_mat)){
  Trt<-design_mat[row,"Treatment"]
  Ctl<-design_mat[row,"Control"]
  varname <- paste("R",row,sep="")

  assign(varname,data.frame(data[,1:2],data[,Ctl],data[,Trt]))
  temp<-get(varname)
  dat[[row]] = temp
  treatment[[row]] = as.character(Trt)
  control[[row]] = as.character(Ctl)
}

count=0

for (i in dat){
  count= count+1
  print(dim(i))
  print(names(i))
  names(i) <- c("sgRNA","Gene",control[count],treatment[count])
  print(names(i))
}

print(dim(dat[[1]]))
print(names(dat[[1]]))