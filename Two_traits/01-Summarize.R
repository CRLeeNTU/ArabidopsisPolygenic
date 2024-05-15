# For the two target traits, here we simply get their original results from step 01 of single-trant analysis
# and then get the SNP intersection for the next step.

ONF <- read.table("../Single_trait/01-ONF_output.txt.gz",header=T,as.is=T,sep="\t")
SS <- read.table("../Single_trait/01-SS_output.txt.gz",header=T,as.is=T,sep="\t")

ONF$SNP <- paste(ONF$chr,ONF$ps,sep=":")
SS$SNP <- paste(SS$chr,SS$ps,sep=":")

ONF <- ONF[ONF$SNP %in% SS$SNP,]
SS <- SS[SS$SNP %in% ONF$SNP,]

# 245,029 SNPs in total
sum(ONF$SNP == SS$SNP)

# The Relict.higher.freq field denotes whether relicts have higher allele 1 freq than non-relict
ONF$Relict.higher.freq <- F
ONF$Relict.higher.freq[ONF$allele1.freq.relict >= ONF$allele1.freq.nonrelict] <- T
SS$Relict.higher.freq <- F
SS$Relict.higher.freq[SS$allele1.freq.relict >= SS$allele1.freq.nonrelict] <- T

# The Relict.higher.freq field denotes whether relicts have higher allele 1 freq than non-relict
# If yes, the allele giving the effect size is correct
# If no, the effect size need to change sign
# So that the final effect size denotes the effect of the "relict allele"

ONF$beta.relict <- ONF$beta
ONF$beta.relict[!(ONF$Relict.higher.freq)] <- 0 - ONF$beta[!(ONF$Relict.higher.freq)]
SS$beta.relict <- SS$beta
SS$beta.relict[!(SS$Relict.higher.freq)] <- 0 - SS$beta[!(SS$Relict.higher.freq)]

ONF$beta.relict.scale <- ONF$beta.relict/sd(ONF$beta.relict)
SS$beta.relict.scale <- SS$beta.relict/sd(SS$beta.relict)

# OK, so the beta.relict.scale column denotes the effect of the "relict allele"
# "relict allele" = the allele with higher allele freq in relict

summary(lm(SS$beta.relict.scale ~ ONF$beta.relict.scale))
cor(SS$beta.relict.scale,ONF$beta.relict.scale) # -0.239596

write.table(ONF,file="01-ONF_output.txt",row.names=F,col.names=T,sep="\t",quote=F)
write.table(SS,file="01-SS_output.txt",row.names=F,col.names=T,sep="\t",quote=F)


