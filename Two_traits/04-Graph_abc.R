# Here plot the top 1% Fst

rm(list=ls())

color.correct <- rgb(red = 148/255, green = 103/255, blue = 189/255, alpha = 0.5)
color.wrong <- rgb(red = 189/255, green = 189/255, blue = 34/255, alpha = 0.5)

tiff(file="04-Graph_abc.tif",width=500,height=1500,res=150)
par(oma=c(0.5,1.5,1.5,0.5),mar=c(5,5,1,1),cex.lab=1.8,cex.axis=1.4,mfrow=c(3,1))

# SS-ONF scatter
my.data <- read.table("02-Data.table.txt.gz",header=T,as.is=T)
my.data <- my.data[order(my.data$Fst,decreasing=T),]
my.model <- summary(lm(my.data$ONF.GWAS.score.signed ~ my.data$SS.GWAS.score.signed))
plot(my.data$SS.GWAS.score.signed,my.data$ONF.GWAS.score.signed,type="n",main="",xlab="SS signed GWAS score",ylab="ONF signed GWAS score")
points(my.data$SS.GWAS.score.signed,my.data$ONF.GWAS.score.signed,pch=16,col=rgb(0.8,0.8,0.8))
# Top 1 % Fst SNP
points(my.data$SS.GWAS.score.signed[which(my.data$Fst >= quantile(my.data$Fst,0.99))],my.data$ONF.GWAS.score.signed[which(my.data$Fst >= quantile(my.data$Fst,0.99))],pch=16,col="orange")
abline(a=my.model$coefficients[1,1],b=my.model$coefficients[2,1],col="blue")
arrows(x0=0,y0=0,x1=4,y1=-4,col="red",lwd=2)
#legend(x=1.8,y=4.4,legend=c("Top 5% Fst SNPs","Other SNPs"),pch=16,col=c("orange",rgb(0.8,0.8,0.8)),cex=1)
#legend(x=1.5,y=4.4,legend=c("Top 1% Fst SNPs","Other SNPs"),pch=16,col=c("orange",rgb(0.8,0.8,0.8)),cex=1)
legend(x=1.5,y=4.4,legend=c(expression(paste('Top 1%') ~ italic(F[ST]) ~ SNPs),"Other SNPs"),pch=16,col=c("orange",rgb(0.8,0.8,0.8)),cex=1)




# SS-ONF box
my.data <- read.table("02-Data.table.txt.gz",header=T,as.is=T)
my.axis <- (my.data$SS.GWAS.score.signed - my.data$ONF.GWAS.score.signed)/2
my.axis.direction <- my.axis > 0
# First use the coordinates to calculate the distance to origin
my.axis <- ((my.axis^2) + (my.axis^2))^0.5
# Then use the direction to assign positive or negative
my.axis[!(my.axis.direction)] <- 0 - my.axis[!(my.axis.direction)]
my.group <- rep("Other SNPs",length(my.axis))
my.group[my.data$Fst >= quantile(my.data$Fst,0.99)] <- "Top 1%"
#my.group[my.data$Fst >= quantile(my.data$Fst,0.95)] <- "Top 5% Fst"
boxplot(my.axis ~ my.group,horizontal=T,col=c(rgb(0.8,0.8,0.8),"orange"),xlab="Axis of SS-ONF divergence",ylab="SNP type")
my.model <- summary(lm(my.axis ~ my.group))
# Top Fst SNPs have highest GWAS score, P < 2.2e-16
text(x=3.5,y=0.65,labels=expression(italic(p) < 2.2e-16),cex=1.5)




# SS-ONF Fst
data.right <- read.table("02-Resample.txt",header=T,as.is=T)
data.wrong <- read.table("03-Resample.txt",header=T,as.is=T)
data.right[data.right < 0] <- 0
data.wrong[data.wrong < 0] <- 0
plot(c(0.5,3.5),c(0,0.4),type="n",xaxt="n",xlab="GWAS Top SNPs",ylab=expression(Median ~ italic(F[ST])))
axis(side=1,at=1:3,labels=c("0.01%","0.1 %","1 %"))
target.col <- data.right$AllSNP.0.0001
points(rnorm(n=1000,mean=0.8,sd=0.03),target.col[2:1001],pch=16,col=color.correct)
lines(c(0.7,0.9),c(target.col[1],target.col[1]),col="black",lwd=3)
target.col <- data.right$AllSNP.0.001
points(rnorm(n=1000,mean=1.8,sd=0.03),target.col[2:1001],pch=16,col=color.correct)
lines(c(1.7,1.9),c(target.col[1],target.col[1]),col="black",lwd=3)
target.col <- data.right$AllSNP.0.01
points(rnorm(n=1000,mean=2.8,sd=0.03),target.col[2:1001],pch=16,col=color.correct)
lines(c(2.7,2.9),c(target.col[1],target.col[1]),col="black",lwd=3)
target.col <- data.wrong$AllSNP.0.0001
points(rnorm(n=1000,mean=1.2,sd=0.03),target.col[2:1001],pch=16,col=color.wrong)
lines(c(1.1,1.3),c(target.col[1],target.col[1]),col="black",lwd=3)
target.col <- data.wrong$AllSNP.0.001
points(rnorm(n=1000,mean=2.2,sd=0.03),target.col[2:1001],pch=16,col=color.wrong)
lines(c(2.1,2.3),c(target.col[1],target.col[1]),col="black",lwd=3)
target.col <- data.wrong$AllSNP.0.01
points(rnorm(n=1000,mean=3.2,sd=0.03),target.col[2:1001],pch=16,col=color.wrong)
lines(c(3.1,3.3),c(target.col[1],target.col[1]),col="black",lwd=3)
points(c(1.8,2.8),c(0.19,0.14),pch="*",cex=3,col="black")
points(c(2.2,3.2),c(0.025,0.04),pch="*",cex=3,col="black")
#text(x=0.7,y=0.385,labels="SS",cex=3)
legend(x=1.4,y=0.4,legend=c("SNPs with correct pleiotropy","SNPs with wrong pleiotropy"),pch=16,col=c("#9467bd","#bcbd22"),cex=1)
text(x=2.4,y=0.3,labels="SS-ONF divergence",cex=1.5)

dev.off()
