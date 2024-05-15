# The wrong direction SNPs have negative GWAS scores

# But here I change it to simply reflect the effect direction of relict alleles
# Simply swapping the results of ONF & Germ in panel A



# The relict.higher.freq column in the output table of step 02 is already polarized for the "UP" allele
# Note that in step 02, I separated alleles by whether they increase trait value,
# and then use the UP allele freq in relicts vs non-relicts to designated their "correct" vs "wrong".
# But when writing paper, I first separate alleles by relict vs non-relict and look at their directions

ONF.step2.output <- read.table("02-ONF_data.table.txt.gz",header=T,as.is=T,sep="\t")
# For the "UP" allele that has higher freq in relict, that's the wrong direction for ONF
ONF.step2.output$direction <- !(ONF.step2.output$relict.higher.freq)
if (F) {
  # Just to double check
  # To prove the "correct" vs "wrong" designaton in the "direction" column is really correct
  # First designate a relict allele, assume to be allele 0, unless allele 1 freq higher in relicts
  ONF.step2.output$Relict_allele <- 0
  ONF.step2.output$Relict_allele[ONF.step2.output$allele1.freq.relict >= ONF.step2.output$allele1.freq.nonrelict] <- 1
  # Then get the effect direction of relict allele
  # GWAS records effect of allele 1, so if relict allele is allele 1, unchange
  # But if relict allele is allele 0, reverse the effect size to 0 - beta
  ONF.step2.output$Relict_allele.effect <- ONF.step2.output$beta
  ONF.step2.output$Relict_allele.effect[ONF.step2.output$Relict_allele == 0] <- 0 - ONF.step2.output$beta[ONF.step2.output$Relict_allele == 0]
  # Then let's look at their association
  # For ONF, relict allele increasing trait value should be WRONG
  # So in the plot, the "direction" column's "FALSE" value should have positive reclit effect
  plot(ONF.step2.output$direction,ONF.step2.output$Relict_allele.effect)
  # Good, so this proves that the "direction" column is correct
  # This is the final designation of whether a SNP is the correct or wrong allele
}

SS.step2.output <- read.table("02-SS_data.table.txt.gz",header=T,as.is=T,sep="\t")
# For the "UP" allele that has higher freq in relict, that's the correct direction for SS
SS.step2.output$direction <- (SS.step2.output$relict.higher.freq)

SS <- read.table("02-SS_resample.txt",header=T,row.names=1,sep="\t",as.is=T)
ONF <- read.table("02-ONF_resample.txt",header=T,row.names=1,sep="\t",as.is=T)

SS[SS < 0] <- 0
ONF[ONF < 0] <- 0

# Get P value table
SS.vec <- c()
ONF.vec <- c()
for (i in 1:ncol(SS)) {
  SS.vec <- c(SS.vec,sum(SS[2:1001,i] >= SS[1,i]) / 1000)
  ONF.vec <- c(ONF.vec,sum(ONF[2:1001,i] >= ONF[1,i]) / 1000)
}

P.table <- data.frame("SS"=SS.vec,"ONF"=ONF.vec)
P.table$Type <- names(SS)

color.correct <- rgb(red = 148/255, green = 103/255, blue = 189/255, alpha = 0.5)
color.wrong <- rgb(red = 189/255, green = 189/255, blue = 34/255, alpha = 0.5)


# Same figures

# This is picking which proportion of SNPs for the top Fst SNPs
my.cutoff <- 0.99




tiff(file="03-Graph_abc.tif",width=1000,height=900,res=150)
par(oma=c(0.5,1.5,1.5,0.5),mar=c(5,5,1,1),cex.lab=1.8,cex.axis=1.4,mfrow=c(2,2))

# SS box
my.data <- SS.step2.output
my.GWAS.score <- -log10(my.data$p_score)
# If a SNP has wrong direction, then GWAS score becomes negative
my.GWAS.score[!(my.data$direction)] <- 0 - my.GWAS.score[!(my.data$direction)]
my.group <- rep(NA,length(my.GWAS.score))
my.group[(my.data$Fst.WEIR_AND_COCKERHAM_FST >= quantile(my.data$Fst.WEIR_AND_COCKERHAM_FST,my.cutoff))] <- "Top 1%"
my.group[(my.data$Fst.WEIR_AND_COCKERHAM_FST < quantile(my.data$Fst.WEIR_AND_COCKERHAM_FST,my.cutoff))] <- "Other SNPs"
final.table <- data.frame("GWAS"=my.GWAS.score,"group"=my.group)
final.table <- final.table[order(final.table$group,decreasing=T),]
boxplot(final.table$GWAS ~ final.table$group,horizontal=T,col=c(rgb(0.8,0.8,0.8),"orange"),xlab="SS signed GWAS score",ylab="SNP type")
my.test <- summary(lm(final.table$GWAS ~ final.table$group))
# Top Fst SNP higher in GWAS, P < 2.2e-16
text(x=5,y=0.65,labels=expression(italic(p) < 2.2e-16),cex=1.5)

mtext("a",side=3,outer=F,cex=2.5,at=-10)

# ONF box
my.data <- ONF.step2.output
my.GWAS.score <- -log10(my.data$p_score)
# If a SNP has wrong direction, then GWAS score becomes negative
my.GWAS.score[!(my.data$direction)] <- 0 - my.GWAS.score[!(my.data$direction)]
# For step 05-8:
# But here I want a graph with relict allele effect size to reflect the pleiotropy graph
# So I invert it again so that the relict allele decreasing trait (correct direction) has negative values
my.GWAS.score <- 0 - my.GWAS.score
my.group <- rep(NA,length(my.GWAS.score))
my.group[(my.data$Fst.WEIR_AND_COCKERHAM_FST >= quantile(my.data$Fst.WEIR_AND_COCKERHAM_FST,my.cutoff))] <- "Top 1%"
my.group[(my.data$Fst.WEIR_AND_COCKERHAM_FST < quantile(my.data$Fst.WEIR_AND_COCKERHAM_FST,my.cutoff))] <- "Other SNPs"
final.table <- data.frame("GWAS"=my.GWAS.score,"group"=my.group)
final.table <- final.table[order(final.table$group,decreasing=T),]
boxplot(final.table$GWAS ~ final.table$group,horizontal=T,col=c(rgb(0.8,0.8,0.8),"orange"),xlab="ONF signed GWAS score",ylab="SNP type")
my.test <- summary(lm(final.table$GWAS ~ final.table$group))
# Top Fst SNP higher in GWAS, P < 2.2e-16
text(x=2.5,y=0.65,labels=expression(italic(p) < 2.2e-16),cex=1.5)

# SS
plot(c(0.5,3.5),c(0,0.4),type="n",xaxt="n",xlab="GWAS Top SNPs",ylab=expression(Median ~ italic(F[ST])))
axis(side=1,at=1:3,labels=c("0.01%","0.1 %","1 %"))

target.col <- SS$relictMORE.0.0001
points(rnorm(n=1000,mean=0.8,sd=0.03),target.col[2:1001],pch=16,col=color.correct)
lines(c(0.7,0.9),c(target.col[1],target.col[1]),col="black",lwd=3)
target.col <- SS$relictMORE.0.001
points(rnorm(n=1000,mean=1.8,sd=0.03),target.col[2:1001],pch=16,col=color.correct)
lines(c(1.7,1.9),c(target.col[1],target.col[1]),col="black",lwd=3)
target.col <- SS$relictMORE.0.01
points(rnorm(n=1000,mean=2.8,sd=0.03),target.col[2:1001],pch=16,col=color.correct)
lines(c(2.7,2.9),c(target.col[1],target.col[1]),col="black",lwd=3)

target.col <- SS$nonrelictMORE.0.0001
points(rnorm(n=1000,mean=1.2,sd=0.03),target.col[2:1001],pch=16,col=color.wrong)
lines(c(1.1,1.3),c(target.col[1],target.col[1]),col="black",lwd=3)
target.col <- SS$nonrelictMORE.0.001
points(rnorm(n=1000,mean=2.2,sd=0.03),target.col[2:1001],pch=16,col=color.wrong)
lines(c(2.1,2.3),c(target.col[1],target.col[1]),col="black",lwd=3)
target.col <- SS$nonrelictMORE.0.01
points(rnorm(n=1000,mean=3.2,sd=0.03),target.col[2:1001],pch=16,col=color.wrong)
lines(c(3.1,3.3),c(target.col[1],target.col[1]),col="black",lwd=3)

points(c(0.8,1.8,2.8),c(0.32,0.185,0.145),pch="*",cex=3,col="black")
points(c(3.2),c(0.1),pch="*",cex=3,col="black")

text(x=0.7,y=0.385,labels="SS",cex=2)
#legend(x=1.9,y=0.4,legend=c("SNPs with correct direction","SNPs with wrong direction"),pch=16,col=c("#9467bd","#bcbd22"),cex=1.5)
legend(x=1.5,y=0.4,legend=c("SNPs with correct direction","SNPs with wrong direction"),pch=16,col=c("#9467bd","#bcbd22"),cex=1)

mtext("b",side=3,outer=F,cex=2.5,at=-0.45)



#ONF
plot(c(0.5,3.5),c(0,0.4),type="n",xaxt="n",xlab="GWAS Top SNPs",ylab=expression(Median ~ italic(F[ST])))
axis(side=1,at=1:3,labels=c("0.01%","0.1 %","1 %"))

target.col <- ONF$nonrelictMORE.0.0001
points(rnorm(n=1000,mean=0.8,sd=0.03),target.col[2:1001],pch=16,col=color.correct)
lines(c(0.7,0.9),c(target.col[1],target.col[1]),col="black",lwd=3)
target.col <- ONF$nonrelictMORE.0.001
points(rnorm(n=1000,mean=1.8,sd=0.03),target.col[2:1001],pch=16,col=color.correct)
lines(c(1.7,1.9),c(target.col[1],target.col[1]),col="black",lwd=3)
target.col <- ONF$nonrelictMORE.0.01
points(rnorm(n=1000,mean=2.8,sd=0.03),target.col[2:1001],pch=16,col=color.correct)
lines(c(2.7,2.9),c(target.col[1],target.col[1]),col="black",lwd=3)

target.col <- ONF$relictMORE.0.0001
points(rnorm(n=1000,mean=1.2,sd=0.03),target.col[2:1001],pch=16,col=color.wrong)
lines(c(1.1,1.3),c(target.col[1],target.col[1]),col="black",lwd=3)
target.col <- ONF$relictMORE.0.001
points(rnorm(n=1000,mean=2.2,sd=0.03),target.col[2:1001],pch=16,col=color.wrong)
lines(c(2.1,2.3),c(target.col[1],target.col[1]),col="black",lwd=3)
target.col <- ONF$relictMORE.0.01
points(rnorm(n=1000,mean=3.2,sd=0.03),target.col[2:1001],pch=16,col=color.wrong)
lines(c(3.1,3.3),c(target.col[1],target.col[1]),col="black",lwd=3)

points(c(2.2,3.2),c(0.02,0.02),pch="*",cex=3,col="black")
points(c(1.8,2.8),c(0.18,0.13),pch="*",cex=3,col="black")

text(x=0.8,y=0.385,labels="ONF",cex=2)
#legend(x=1.9,y=0.4,legend=c("SNPs with correct direction","SNPs with wrong direction"),pch=16,col=c("#9467bd","#bcbd22"),cex=1.5)
legend(x=1.5,y=0.4,legend=c("SNPs with correct direction","SNPs with wrong direction"),pch=16,col=c("#9467bd","#bcbd22"),cex=1)



dev.off()

