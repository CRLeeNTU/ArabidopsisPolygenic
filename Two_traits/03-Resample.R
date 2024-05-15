rm(list=ls())
my.data <- read.table("02-Data.table.txt.gz",header=T,as.is=T)

# Use the data from step 2, but here we focus on the top SNPs in the "wrong direction"

# Now the top has to satisfy: 1. !(correct direction in both), 2. top X % SNP
top.num <- ceiling(nrow(my.data)*0.0001)
cutoff <- sort(my.data$composite.GWAS.score[!(my.data$correct.direction)],decreasing=T)[top.num]
my.data$top.0.0001 <- F
my.data$top.0.0001[which(!(my.data$correct.direction) & my.data$composite.GWAS.score >= cutoff)] <- T

top.num <- ceiling(nrow(my.data)*0.001)
cutoff <- sort(my.data$composite.GWAS.score[!(my.data$correct.direction)],decreasing=T)[top.num]
my.data$top.0.001 <- F
my.data$top.0.001[which(!(my.data$correct.direction) & my.data$composite.GWAS.score >= cutoff)] <- T

top.num <- ceiling(nrow(my.data)*0.01)
cutoff <- sort(my.data$composite.GWAS.score[!(my.data$correct.direction)],decreasing=T)[top.num]
my.data$top.0.01 <- F
my.data$top.0.01[which(!(my.data$correct.direction) & my.data$composite.GWAS.score >= cutoff)] <- T

# Check
#plot(my.data$SS.GWAS.score.signed,my.data$ONF.GWAS.score.signed)
#points(my.data$SS.GWAS.score.signed[my.data$top.0.0001],my.data$ONF.GWAS.score.signed[my.data$top.0.0001],pch=16,col="red")
#points(my.data$SS.GWAS.score.signed[my.data$top.0.001],my.data$ONF.GWAS.score.signed[my.data$top.0.001],pch=16,col="red")
#points(my.data$SS.GWAS.score.signed[my.data$top.0.01],my.data$ONF.GWAS.score.signed[my.data$top.0.01],pch=16,col="red")






# Let's resample background SNPs
# Which allele freq bin to resample from, for each top SNP
# Need to change here for other groups
####################################################################################
    MAF.bin.top.SNP <- my.data$MAF.bin[my.data$top.0.01]
    Fst.median.vec <- median(my.data$Fst[my.data$top.0.01])
####################################################################################

# No specific meaning for creating data.use, just a legacy of old codes
data.use <- my.data
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
      temp.Fst.median <- c(temp.Fst.median,sample(data.use.this.bin$Fst,TOP.SNP.num.this.bin,replace=F))
    }
  }
  Fst.median.vec <- c(Fst.median.vec,median(temp.Fst.median))
}
sum(Fst.median.vec[2:length(Fst.median.vec)] > Fst.median.vec[1]) / (length(Fst.median.vec)-1)

# I do it by hand, one by one
#output.table <- data.frame("AllSNP.0.0001"=Fst.median.vec)
#row.names(output.table) <- c("True",paste0("background",1:1000))

#output.table$AllSNP.0.001 <- Fst.median.vec
output.table$AllSNP.0.01 <- Fst.median.vec

write.table(output.table,file="03-Resample.txt",row.names=T,col.names=T,sep="\t",quote=F)
