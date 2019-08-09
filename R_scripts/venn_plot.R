#Create Venn Diagrams from predicitons from different datasets
#Inputs; 1. Path to results directory 
#        2. Condition to check
#        3. Top number of genes to be checked 

#Load required libraries 
library(VennDiagram)

#Define input argument
folder="~/github/results_new/results_XL27/"
condition_vs_cl="tumor_vs_pi_cas9"
num_genes = 100

####################################################################
#MageckRRA
####################################################################

#Read in results file for mageck
mageckrra_file <- paste(folder,"mageckrra/mageck_",condition_vs_cl,".gene_summary.txt",sep="")
result_rra<-read.table(mageckrra_file,header=TRUE)

#Order according to mageck ranking
result_rra<-result_rra[order(pmin(result_rra$pos.rank, result_rra$neg.rank),decreasing = TRUE),]
top_rra<-head(result_rra[,1],n=num)

####################################################################
#PBNPA
####################################################################

#Read in results for pbnpa
pbnpa_file <- paste(folder,"pbnpa/pbnpa_",condition_vs_cl,".txt",sep="")
result_pbnpa<-read.table(pbnpa_file,header=TRUE)

#Sort according to pos.fdr 
result_pbnpa <- result_pbnpa[order(pmin(result_pbnpa$pos.fdr,result_pbnpa$neg.fdr), decreasing = FALSE),]
top_pbnpa<-head(result_pbnpa[,1], n=num)

####################################################################
#CB2
####################################################################

#Read in results for cb2
cb2_file <- paste(folder,"cb2/cb2_",condition_vs_cl,".txt",sep="")
result_cb2<-read.table(cb2_file,header=TRUE)

#Sort data according to fdr and logFC
sorted_data_cb2 <- result_cb2[order(abs(result_cb2$logFC), decreasing = TRUE),]
fc_sorted_data_cb2 <- sorted_data_cb2[order(sorted_data_cb2$fdr_ts, decreasing = FALSE),]

top_cb2 <- head(fc_sorted_data_cb2$gene,n=num)

####################################################################
#MageckMLE
####################################################################
mageckmle_file <- paste(folder,"mageckmle/mle_",condition_vs_cl,".gene_summary.txt",sep="")
mageckmle_file <- "~/github/results_new/results_XL27/mageckmle/results_tumor_vs_pi_ot1.gene_summary.txt"
result_mle<-read.table(mageckmle_file,header=TRUE)

#Keep only required columns and sort according to fdr and beta
cols <- c(1,3,8)
req_mle<-data.frame(result_mle[,cols])
names(req_mle)<-names(result_mle)[cols]

#Order according to lfc and fdr
req_mle <- req_mle[order(abs(req_mle[,2]),decreasing = T),]
req_mle <- req_mle[order(req_mle[,3],decreasing = F),]

top_mle<-head(req_mle[,1], n=num)

####################################################################
#VENN DIAGRAM
####################################################################

top_genes<-list(top_rra,top_pbnpa,top_mle,top_cb2)
names(top_genes)<-c("MageckRRA","PBNPA","MageckMLE","CB2")
top_genes <- lapply(as.list(top_genes), function(x) x[x != ""])

VENN.LIST <- top_genes
venn.plot <- venn.diagram(VENN.LIST , condition_vs_cl, fill=c( 'cornflowerblue','turquoise','magenta', 'darkorange'), alpha=c(0.5,0.5,0.5,0.5), cex = 1.5, fontface = 7, cat.fontface=4, category.names=names(top_genes), main=condition_vs_cl)
