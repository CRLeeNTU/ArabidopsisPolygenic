# SS
rm(list=ls())
my.data <- read.table("01-SS_output.txt.gz",header=T)
my.data$freq.diff <- my.data$UPallele.freq.relict - my.data$UPallele.freq.nonrelict
my.data$relict.higher.freq <- NA
my.data$relict.higher.freq[which(my.data$freq.diff > 0)] <- T
my.data$relict.higher.freq[which(my.data$freq.diff < 0)] <- F
# Calculate minor allele freq
my.data$UPallele.freq.mean <- (my.data$UPallele.freq.relict + my.data$UPallele.freq.nonrelict)/2
my.data$minor.allele.freq <- 1 - (0.5+abs(0.5-my.data$UPallele.freq.mean))

my.data <- my.data[complete.cases(my.data),]
# 256245 to 241125

# Group SNPs into 10 bins of MAF
my.data$MAF.bin <- NA
bin.start <- seq(0,0.45,0.05)
for (i in 1:length(bin.start)) {
  my.data$MAF.bin[(my.data$minor.allele.freq > bin.start[i]) & (my.data$minor.allele.freq <= bin.start[i]+0.05)] <- i
}
my.data <- my.data[complete.cases(my.data),]  # 241112
my.data$top.0.0001 <- "BOT"
my.data$top.0.0001[which(my.data$p_score <= quantile(my.data$p_score,0.0001))] <- "TOP"
my.data$top.0.001 <- "BOT"
my.data$top.0.001[which(my.data$p_score <= quantile(my.data$p_score,0.001))] <- "TOP"
my.data$top.0.01 <- "BOT"
my.data$top.0.01[which(my.data$p_score <= quantile(my.data$p_score,0.01))] <- "TOP"
my.data$top.0.1 <- "BOT"
my.data$top.0.1[which(my.data$p_score <= quantile(my.data$p_score,0.1))] <- "TOP"
my.data$top.0.2 <- "BOT"
my.data$top.0.2[which(my.data$p_score <= quantile(my.data$p_score,0.2))] <- "TOP"
my.data$top.0.3 <- "BOT"
my.data$top.0.3[which(my.data$p_score <= quantile(my.data$p_score,0.3))] <- "TOP"

write.table(my.data,file="02-SS_data.table.txt",quote=F,col.names=T,row.names=F,sep="\t")




#################################################################################################
##### Modify here for directions. #

##### If relict has higher trait value, this is the right-direction SNP
##### Reverse it for traits where relict has lower trait value
#data.use <- my.data[(my.data$relict.higher.freq),]

##### If relict has higher larger trait value, this is the wrong-direction SNP
##### Reverse it for traits where relict has lower trait value
#data.use <- my.data[!(my.data$relict.higher.freq),]

##### Modify here to pick the top X proportion SNPs
##### This is not done automatically, need to change the proportion by hand and do it one my one
#MAF.bin.top.SNP <- data.use$MAF.bin[data.use$top.0.0001 == "TOP"]
#Fst.median.vec <- median(data.use$Fst.WEIR_AND_COCKERHAM_FST[data.use$top.0.0001 == "TOP"])
#################################################################################################
  
  
# Generate a list. data.use.list[[1]] contains all SNPs in MAF.bin == 1, etc.
data.use.list <- list()
for (i in 1:10) {
  data.use.list[[i]] <- data.use[data.use$MAF.bin == i,]
}
# Re-sample background SNPs
perm.num <- 1000
for (i in 1:perm.num) {
  temp.Fst.median <- c()
  for (j in 1:10) {
    data.use.this.bin <- data.use.list[[j]]
    TOP.SNP.num.this.bin <- sum(MAF.bin.top.SNP == j)
    if (TOP.SNP.num.this.bin > 0) {
      temp.Fst.median <- c(temp.Fst.median,sample(data.use.this.bin$Fst.WEIR_AND_COCKERHAM_FST,TOP.SNP.num.this.bin,replace=F))
    }
  }
  Fst.median.vec <- c(Fst.median.vec,median(temp.Fst.median))
}
sum(Fst.median.vec[2:length(Fst.median.vec)] > Fst.median.vec[1]) / (length(Fst.median.vec)-1)





# So for the section above, run it separately by hand one my one 

# I do it by hand, one by one
output.table <- data.frame("relictMORE.0.0001"=Fst.median.vec)
row.names(output.table) <- c("True",paste0("background",1:1000))
output.table$relictMORE.0.0001 <- Fst.median.vec
output.table$relictMORE.0.001 <- Fst.median.vec
output.table$relictMORE.0.01 <- Fst.median.vec

output.table$nonrelictMORE.0.0001 <- Fst.median.vec
output.table$nonrelictMORE.0.001 <- Fst.median.vec
output.table$nonrelictMORE.0.01 <- Fst.median.vec

write.table(output.table,file="02-SS_resample.txt",row.names=T,col.names=T,sep="\t",quote=F)
