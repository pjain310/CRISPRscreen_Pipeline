data<-read.table("readcount-DLD1-lib1",header=TRUE)

R1<-data.frame(data["GENE_CLONE"],data["GENE"],data["DLD_T0"],data["DLD_ETOH_R1"])
R2<-data.frame(data["GENE_CLONE"],data["GENE"],data["DLD_T0"],data["DLD_ETOH_R2"])
R3<-data.frame(data["GENE_CLONE"],data["GENE"],data["DLD_T0"],data["DLD_ETOH_R3"])

names(R1)<-c("sgRNA","Gene","Control","DLD_ETOH_R1")
names(R2)<-c("sgRNA","Gene","Control","DLD_ETOH_R2")
names(R3)<-c("sgRNA","Gene","Control","DLD_ETOH_R3")

dat<-list(R1,R2,R3)

results<-PBNPA(dat)