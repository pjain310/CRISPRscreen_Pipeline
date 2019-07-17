library(dplyr)

data <- read.table("counts_250g_p.25_4gd_std.05_mod.txt",sep="\t",header = T)
data$X <- NULL

df <- data
n <- round(100000)
m <- 2.5
opfile <- "counts_250g_mod_prolif_dirty_3r.txt"
#Set weights 
df$pro1 <- df$counts_bin1/sum(df$counts_bin1)*m + 0.000005 
df$pro2 <- df$counts_bin2/sum(df$counts_bin2)*m + 0.000005

#Noise generating function
makesomenoise <- function(counts,noise_percentage){
  corrupt <- rbinom(length(counts),1,noise_percentage)    # choose an average of 10% to corrupt at random
  corrupt <- as.logical(corrupt)
  noise <- rnorm(length(corrupt),40,20) # generate the noise to add
  counts[corrupt] <- counts[corrupt] + round(noise[corrupt])
  
  return(counts)
}

#Boostrap samples for each guide 
control_1 <- as.data.frame(table(sample(df$barcodeid,n,replace = T, prob = df$pro1)))
names(control_1) <- c("Var1","Control_1")
control_2 <- as.data.frame(table(sample(df$barcodeid,n,replace = T, prob = df$pro1)))
names(control_2) <- c("Var1","Control_2")
control_3 <- as.data.frame(table(sample(df$barcodeid,n,replace = T, prob = df$pro1)))
names(control_3) <- c("Var1","Control_3")
#control_4 <- as.data.frame(table(sample(df$barcodeid,n,replace = T, prob = df$pro1)))
#names(control_4) <- c("Var1","Control_4")
#control_5 <- as.data.frame(table(sample(df$barcodeid,n,replace = T, prob = df$pro1)))
#names(control_5) <- c("Var1","Control_5")

case_1 <- as.data.frame(table(sample(df$barcodeid,n,replace = T, prob = df$pro2)))
names(case_1) <- c("Var1","Condn_1")
case_2 <- as.data.frame(table(sample(df$barcodeid,n,replace = T, prob = df$pro2)))
names(case_2) <- c("Var1","Condn_2")
case_3 <- as.data.frame(table(sample(df$barcodeid,n,replace = T, prob = df$pro2)))
names(case_3) <- c("Var1","Condn_3")
#case_4 <- as.data.frame(table(sample(df$barcodeid,n,replace = T, prob = df$pro2)))
#names(case_4) <- c("Var1","Condn_4")
#case_5 <- as.data.frame(table(sample(df$barcodeid,n,replace = T, prob = df$pro2)))
#names(case_5) <- c("Var1","Condn_5")


genes <- data.frame(Var1=df$barcodeid,Gene=df$gene)
#counts<-Reduce(function(...) merge(..., all=TRUE, by = "Var1"), list(genes,control_1, control_2, control_3,control_4,control_5,case_1,case_2,case_3,case_4,case_5))
counts<-Reduce(function(...) merge(..., all=TRUE, by = "Var1"), list(genes,control_1, control_2, control_3, case_1,case_2,case_3))
#counts<-Reduce(function(...) merge(..., all=TRUE, by = "Var1"), list(genes,control_1, control_2, case_1,case_2))
#counts<-Reduce(function(...) merge(..., all=TRUE, by = "Var1"), list(genes,control_1, case_1))
counts <- cbind(sgRNA=0,counts)
counts$sgRNA <- paste(counts$Gene,"_",counts$Var1,sep="")
counts$Var1 <- NULL
row.names(counts) <- NULL

#Replace NA with median 
for(i in 3:ncol(counts)){
  counts[is.na(counts[,i]), i] <- median(counts[,i], na.rm = TRUE)
}

#Introduce noise
for(i in 3:ncol(counts)){
  counts[,i] <- makesomenoise(counts[,i],0.)
}
head(counts)

hist(log2(data$counts_bin2/data$counts_bin1),breaks = 20)
hist(log2(counts$Condn_1/counts$Control_1),breaks = 20)

# write.table(file=opfile,counts,sep="\t",quote=F,row.names = F)


