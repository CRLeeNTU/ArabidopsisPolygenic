rm(list=ls())
ONF <- read.table("01-ONF_output.txt.gz",header=T)
SS <- read.table("01-SS_output.txt.gz",header=T)
my.data <- ONF[,1:4]
my.data$ONF.beta.relict.scale <- ONF$beta.relict.scale
my.data$SS.beta.relict.scale <- SS$beta.relict.scale
my.data$ONF.P <- ONF$p_score
my.data$SS.P <- SS$p_score
my.data$Fst <- ONF$Fst.WEIR_AND_COCKERHAM_FST

# Calculate minor allele freq, assuming the two populations have the same sample size
my.data$allele.freq.mean <- (ONF$allele1.freq.relict + ONF$allele1.freq.nonrelict)/2
my.data$minor.allele.freq <- 1 - (0.5+abs(0.5-my.data$allele.freq.mean))

my.data <- my.data[complete.cases(my.data),]
rm(ONF,SS)

# -log10P & their composite (geometric mean)
# The log of geometric mean of P simply equals the mean of their log
my.data$ONF.GWAS.score <- 0 - log10(my.data$ONF.P)
my.data$SS.GWAS.score <- 0 - log10(my.data$SS.P)
my.data$composite.GWAS.score <- (my.data$ONF.GWAS.score + my.data$SS.GWAS.score) / 2

# Get signed GWAS score
# This is the same as GWAS scoore, but the sign (+ or -) follows relict allele effect size
my.data$ONF.GWAS.score.signed <- my.data$ONF.GWAS.score
my.data$ONF.GWAS.score.signed[which(my.data$ONF.beta.relict.scale < 0)] <- 0 - my.data$ONF.GWAS.score[which(my.data$ONF.beta.relict.scale < 0)]
my.data$SS.GWAS.score.signed <- my.data$SS.GWAS.score
my.data$SS.GWAS.score.signed[which(my.data$SS.beta.relict.scale < 0)] <- 0 - my.data$SS.GWAS.score[which(my.data$SS.beta.relict.scale < 0)]
# Check
#plot(my.data$ONF.GWAS.score.signed,my.data$ONF.beta.relict.scale)
#plot(my.data$SS.GWAS.score.signed,my.data$SS.beta.relict.scale)


# Group SNPs into 10 bins of MAF
my.data$MAF.bin <- NA
bin.start <- seq(0,0.45,0.05)
for (i in 1:length(bin.start)) {
  my.data$MAF.bin[(my.data$minor.allele.freq > bin.start[i]) & (my.data$minor.allele.freq <= bin.start[i]+0.05)] <- i
}
my.data <- my.data[complete.cases(my.data),]


# Correct directrion
# Note that in the example here, relicts have higher SS but lower ONF
my.data$correct.direction <- F
my.data$correct.direction[which(my.data$SS.beta.relict.scale >= 0 & my.data$ONF.beta.relict.scale <= 0)] <- T
#plot(my.data$SS.GWAS.score.signed,my.data$ONF.GWAS.score.signed)
#points(my.data$SS.GWAS.score.signed[my.data$correct.direction],my.data$ONF.GWAS.score.signed[my.data$correct.direction],pch=16,col="red")



# Now the top SNPs have to satisfy: 1. correct direction in both (the correct quadrant), 2. top X % SNP
top.num <- ceiling(nrow(my.data)*0.0001)
cutoff <- sort(my.data$composite.GWAS.score[my.data$correct.direction],decreasing=T)[top.num]
my.data$top.0.0001 <- F
my.data$top.0.0001[which(my.data$correct.direction & my.data$composite.GWAS.score >= cutoff)] <- T

top.num <- ceiling(nrow(my.data)*0.001)
cutoff <- sort(my.data$composite.GWAS.score[my.data$correct.direction],decreasing=T)[top.num]
my.data$top.0.001 <- F
my.data$top.0.001[which(my.data$correct.direction & my.data$composite.GWAS.score >= cutoff)] <- T

top.num <- ceiling(nrow(my.data)*0.01)
cutoff <- sort(my.data$composite.GWAS.score[my.data$correct.direction],decreasing=T)[top.num]
my.data$top.0.01 <- F
my.data$top.0.01[which(my.data$correct.direction & my.data$composite.GWAS.score >= cutoff)] <- T

# Check
#plot(my.data$SS.GWAS.score.signed,my.data$ONF.GWAS.score.signed)
#points(my.data$SS.GWAS.score.signed[my.data$top.0.0001],my.data$ONF.GWAS.score.signed[my.data$top.0.0001],pch=16,col="red")
#points(my.data$SS.GWAS.score.signed[my.data$top.0.001],my.data$ONF.GWAS.score.signed[my.data$top.0.001],pch=16,col="red")
#points(my.data$SS.GWAS.score.signed[my.data$top.0.01],my.data$ONF.GWAS.score.signed[my.data$top.0.01],pch=16,col="red")

write.table(my.data,file="02-Data.table.txt",quote=F,col.names=T,row.names=F,sep="\t")








# Let's resample background SNPs
# Which allele freq bin to resample from, for each top SNP
# Need to change here for other groups, do it by hand for top 0.01, 0.001, 0.0001
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


#######################################################################################
# After doing the section avobe one by one by hand, for each step save the results here
#output.table <- data.frame("AllSNP.0.0001"=Fst.median.vec)
#row.names(output.table) <- c("True",paste0("background",1:1000))
#output.table$AllSNP.0.001 <- Fst.median.vec
#output.table$AllSNP.0.01 <- Fst.median.vec

write.table(output.table,file="02-Resample.txt",row.names=T,col.names=T,sep="\t",quote=F)
