#!/usr/local/bin/R

folder="../results_new/results_Evers/"
condition_vs_cl="t1RT112_vs_t0RT112"

####################################################################
#MageckMLE
####################################################################
mageckmle_file <- paste(folder,"mageckmle/results_",condition_vs_cl,".gene_summary.txt",sep="")

result_mle<-read.table(mageckmle_file,header=TRUE)

#Keep only required columns and sort according to fdr and beta
cols <- c(1,3,8)
req_mle<-data.frame(result_mle[,cols])
names(req_mle)<-names(result_mle)[cols]
req_mle <- req_mle[order(abs(req_mle[,2]),decreasing = T),]
req_mle <- req_mle[order(req_mle[,3],decreasing = F),]


#Add fdr cutoff
fdr_cutoff <- result_mle[,8] < 0.05
result_mle <- result_mle[fdr_cutoff,]

#Print top 20 ranked genes
print(head(req_mle, n=20))

top_mle<-head(req_mle, n=100)

####################################################################
#MageckRRA
####################################################################

#Read in results file for mageck
mageckrra_file <- paste(folder,"mageckrra/mageck_",condition_vs_cl,".gene_summary.txt",sep="")
file <- "results_XL27_sep_ctls/mageckrra/tumor_vs_pi_ot1.gene_summary.txt"

result_rra<-read.table(mageckrra_file,header=TRUE)
result_rra<-result_rra[order(result_rra$pos.fdr,decreasing = FALSE),]

#Add fdr cutoff
fdr_cutoff <- result_rra$neg.fdr < 0.05 | result_rra$pos.fdr < 0.05
result_rra <- result_rra[fdr_cutoff,]

#Print genes and fdrs
req <- data.frame(result_rra$id,result_rra$neg.fdr,result_rra$pos.fdr,result_rra$neg.lfc)
print(head(req, n=20))

top_rra<-head(req, n=100)

####################################################################
#PBNPA
####################################################################

#Read in results for pbnpa
pbnpa_file <- paste(folder,"pbnpa/pbnpa_",condition_vs_cl,".txt",sep="")
result_pbnpa<-read.table(pbnpa_file,header=TRUE)

#Add fdr cutoff
fdr_cutoff <- result_pbnpa$neg.fdr < 0.05 | result_pbnpa$pos.fdr < 0.05
result_pbnpa_filt <- result_pbnpa[fdr_cutoff,]

#Sort according to pos.fdr and print
sorted_data_pbnpa <- result_pbnpa_filt[order(result_pbnpa_filt$pos.fdr, decreasing = FALSE),]
#print(head(sorted_data_pbnpa,n=20))

top_pbnpa<-head(result_pbnpa, n=100)

####################################################################
#CB2
####################################################################

#Read in results for cb2
cb2_file <- paste(folder,"cb2/cb2_",condition_vs_cl,".txt",sep="")
result_cb2<-read.table(cb2_file,header=TRUE)

#Add fdr cutoff
fdr_cutoff <- result_cb2$fdr_ts < 0.05
result_cb2_filt <- result_cb2[fdr_cutoff,]

#Sort data according to fdr and logFC
sorted_data_cb2 <- result_cb2_filt[order(abs(result_cb2_filt$logFC), decreasing = TRUE),]
fc_sorted_data_cb2 <- sorted_data_cb2[order(sorted_data_cb2$fdr_ts, decreasing = FALSE),]

fc_sorted_data_cb2 <- data.frame(fc_sorted_data_cb2$gene,fc_sorted_data_cb2$logFC,fc_sorted_data_cb2$fdr_ts)
print(head(fc_sorted_data_cb2, n = 20))

####################################################################
#BAGEL
####################################################################

#Read in results for bagel
bagel_file <- paste(folder,"bagel/bagel_",condition_vs_cl,sep="")
result_bagel<-read.table(bagel_file,header=TRUE)

#Add fdr cutoff
fdr_cutoff <- abs(result_bagel$BF) > 15
result_bagel_filt <- result_bagel[fdr_cutoff,]

#Sort according to bayes factor and print
sorted_data_bagel <- result_bagel_filt[order(abs(result_bagel_filt$BF),decreasing = TRUE),]
print(head(sorted_data_bagel,n=20))

