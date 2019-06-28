#!/usr/local/bin/R

folder="../results_new/results_KBM/"
condition_vs_cl="KBM7_final_vs_KBM7_initial"

#Define historgram function
plot_hist <- function(d,header,bw = 0.025){

  x<-d[,1]
  ggplot2::ggplot(d, aes(x, fill = cut(x, 100))) +
    geom_histogram(show.legend = FALSE,binwidth = bw, color="darkslategray", fill = "skyblue") +
    scale_fill_discrete(h = c(240, 10), c = 120, l = 70) +
    labs(y = "Frequency") +
    ggtitle(header)
}

####################################################################
#MageckRRA
####################################################################

#Read in results file for mageck
mageckrra_file <- paste(folder,"mageckrra/mageck_",condition_vs_cl,".gene_summary.txt",sep="")
file <- "results_XL27_sep_ctls/mageckrra/tumor_vs_pi_ot1.gene_summary.txt"

result_rra<-read.table(mageckrra_file,header=TRUE)
result_rra<-result_rra[order(result_rra$pos.fdr,decreasing = FALSE),]

#Plot FDR Distribtuion
plot_file <- paste("../plots/",condition_vs_cl,"_mageckrra.jpeg",sep="")
jpeg(plot_file)
plot_hist(data.frame(result_rra$pos.fdr),"FDR Distribution for MageckRRA")
dev.off()

#Add fdr cutoff
fdr_cutoff <- result_rra$neg.fdr < 0.05 | result_rra$pos.fdr < 0.05
result_rra <- result_rra[fdr_cutoff,]

#Print genes and fdrs
req <- data.frame(result_rra$id,result_rra$neg.fdr,result_rra$pos.fdr,result_rra$neg.lfc)
print(head(req, n=20))

top_rra<-result_rra[,1]

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
#sorted_data_pbnpa <- result_pbnpa_filt[order(result_pbnpa_filt$pos.fdr, decreasing = FALSE),]
#print(head(sorted_data_pbnpa,n=20))

plot_file <- paste("../plots/",condition_vs_cl,"_mageckrra.png",sep="")
png(plot_file)
plot_hist(data.frame(result_pbnpa$pos.fdr),"FDR Distribution for PBNPA")
dev.off()

top_pbnpa<-result_pbnpa_filt[,1]

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
sorted_data_bagel <- result_bagel_filt[order(result_bagel_filt$BF,decreasing = TRUE),]
print(head(sorted_data_bagel,n=20))

plot_file <- paste("../plots/",condition_vs_cl,"_mageckrra.png",sep="")
png(plot_file)
plot_hist(data.frame(result_bagel$BF),"BF Distribution for BAGEL", 2)
dev.off()
top_bagel<-result_bagel_filt[,1]

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

plot_file <- paste("../plots/",condition_vs_cl,"_mageckrra.png",sep="")
png(plot_file)
plot_hist(data.frame(result_cb2$fdr_ts),"FDR Distribution for CB2")
dev.off()
top_cb2<-result_cb2_filt[,1]


