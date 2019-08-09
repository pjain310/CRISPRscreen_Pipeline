#!/usr/local/bin/R
"
Script to generate confusion matrix for preedictions by different tools 

Parameters
----------
folder = path to results directory for simulated data 
condition_vs_cl = name of condition vs control 
fdr_threshold = fdr threshold to qualify sample as a hit or not 

Returns
-------
Conf: confusion matrix file 
"

#Load required libraries
library(reshape)

#Read input arguments (output directory path and condition vs control name)
args = commandArgs(trailingOnly=TRUE)

folder=args[1]
condition_vs_cl=args[2]
fdr_threshold = args[3]

####################################################################
#MageckMLE
####################################################################
mageckmle_file <- paste(folder,"mageckmle/mle_",condition_vs_cl,".gene_summary.txt",sep="")
print(mageckmle_file)
result_mle<-read.table(mageckmle_file,header=TRUE)

####################################################################
#MageckRRA
####################################################################
#Read in results file for mageck
mageckrra_file <- paste(folder,"mageckrra/mageck_",condition_vs_cl,".gene_summary.txt",sep="")
result_rra<-read.table(mageckrra_file,header=TRUE)

####################################################################
#PBNPA
####################################################################

#Read in results for pbnpa
pbnpa_file <- paste(folder,"pbnpa/pbnpa_",condition_vs_cl,".txt",sep="")
result_pbnpa<-read.table(pbnpa_file,header=TRUE)

####################################################################
#CB2
####################################################################
#Read in results for cb2
condition_vs_cl <- "Control_vs_Condition"
cb2_file <- paste(folder,"cb2/cb2_",condition_vs_cl,".txt",sep="")
result_cb2<-read.table(cb2_file,header=TRUE)

#Correcting gene names for cb2 simulated data
gnum <- colsplit(result_pbnpa$Gene, split="_",names=c('gnum','type'))
result_cb2<-merge(result_cb2, gnum, by.y = "gnum", by.x = "gene")
result_cb2$Gene <- paste(result_cb2$gene,"_",result_cb2$type,sep="")
print(result_cb2)


####################################################################
#Confusion matrix
####################################################################

cmat <- data.frame(Gene = result_cb2$Gene)
cmat$Ground_Truth <- as.numeric(grepl("inc|dec",cmat$Gene))
cmat$MLE <- as.numeric(cmat$Gene %in% result_mle[result_mle$Condition.wald.fdr < fdr_threshold,]$Gene)
cmat$RRA <- as.numeric(cmat$Gene %in% result_rra[result_rra$neg.fdr< fdr_threshold,]$id) + as.numeric(cmat$Gene %in% result_rra[result_rra$pos.fdr < fdr_threshold,]$id)
cmat$PBNPA <- as.numeric(cmat$Gene %in% result_pbnpa[result_pbnpa$neg.fdr < fdr_threshold,]$Gene) + as.numeric(cmat$Gene %in% result_pbnpa[result_pbnpa$pos.fdr < fdr_threshold,]$Gene)
cmat$CB2 <- as.numeric(cmat$Gene %in% result_cb2[(result_cb2$fdr_ts < fdr_threshold & abs(result_cb2$logFC) > 1.5),]$Gene)

write.table(file=paste("~/github/outputs/confusion_matrix","fname",".txt",sep=""),cmat,sep="\t",quote=F,row.names = F)

