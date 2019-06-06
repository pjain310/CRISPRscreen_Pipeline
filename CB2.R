library(CB2)
library(tibble)

df_design <- tribble(~group,~sample_name,"Base","DLD_T0","EtOH","DLD_ETOH_R1","EtOH","DLD_ETOH_R2","EtOH","DLD_ETOH_R3")
df<-read.table("readcount-DLD1-lib1",header = TRUE)
rownames(df)<-df[,1]
df<-df[,-c(1:2)]

results<-run_estimation(df,df_design,"Base","EtOH")

write.table(results,"results_CB2.txt")