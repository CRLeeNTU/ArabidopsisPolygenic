# In theory the SNPs should be exactly the same as in GWAS result file

rm(list=ls())

relict.file <- "100_relict_SS-GWAS-SNP_freq.frq.gz"
nonrelict.file <- "99_non-relict_SS-GWAS-SNP_freq.frq.gz"
GWAS.file <- "SS_nonLDprune-kinship_LDprune.assoc_part.txt.gz"
output.file <- "01-SS_output.txt"
fst.file <- "199_SS-GWAS-SNP_fst.weir.fst.gz"

# There's no "," in the allele freq file, so I used this as field separator. This gets everything into one line.
# I do this because the field number is different for SNPs, so can't really get them into tables
relict <- read.table(file=relict.file,header=F,as.is=T,sep=",",skip=1)
# Get rid of the "Chr" in the beginning of chr names
relict[,1] <- substr(relict[,1],start=4,stop=1000)

nonrelict <- read.table(file=nonrelict.file,header=F,as.is=T,sep=",",skip=1)
nonrelict[,1] <- substr(nonrelict[,1],start=4,stop=1000)

GWAS <- read.table(file=GWAS.file,header=T,sep="\t",as.is=T)
rownames(GWAS) <- GWAS$rs

Fst <- read.table(file=fst.file,header=T,as.is=T)

freq <- cbind(relict,nonrelict)
names(freq) <- c("relict","nonrelict")
rownames(freq) <- paste(sapply(strsplit(freq[,1],"\t"),"[",1),sapply(strsplit(freq[,1],"\t"),"[",2),sep=":")

all.data <- merge(GWAS,freq,by="row.names")
all.data <- all.data[,-1]
all.data <- all.data[order(all.data$chr,all.data$ps),]

rm(GWAS,relict,nonrelict,GWAS.file,relict.file,nonrelict.file,freq)

# Check whether the row numbers are the same across the 3 files~!!!

# Frequency of allele 1
allele1.freq.relict.vec <- c()
allele1.freq.nonrelict.vec <- c()
# Frequency of trait-increasing allele
UPallele.freq.relict.vec <- c()
UPallele.freq.nonrelict.vec <- c()

for (i in 1:nrow(all.data)) {
  this.line <- all.data[i,]
  relict.line <- strsplit(this.line$relict,"\t")[[1]]
  nonrelict.line <- strsplit(this.line$nonrelict,"\t")[[1]]
  
  # Then get the element with allele1 in vector relict.line & nonrelict.line
  # There will be some lines that allele1 does not exist in allele freq file
  # At those sites, allele1.relict will be chracter(0)
  # In that case, create a new allele1.relict and give allele freq 0
  allele1.relict <- relict.line[grep(pattern=paste0("^",this.line$allele1,":"),x=relict.line)]
  if(length(allele1.relict) == 0) {
    allele1.relict <- paste0(this.line$allele1,":0")
  }
  allele1.nonrelict <- nonrelict.line[grep(pattern=paste0("^",this.line$allele1,":"),x=nonrelict.line)]
  if(length(allele1.nonrelict) == 0) {
    allele1.nonrelict <- paste0(this.line$allele1,":0")
  }
  
  # Remove the "A:" from "A:0.85"
  allele1.freq.relict <- as.numeric(substr(allele1.relict,3,9))
  allele1.freq.nonrelict <- as.numeric(substr(allele1.nonrelict,3,9))
  
  # Append to the vector
  allele1.freq.relict.vec <- c(allele1.freq.relict.vec,allele1.freq.relict)
  allele1.freq.nonrelict.vec <- c(allele1.freq.nonrelict.vec,allele1.freq.nonrelict)
  
  # Then calculate the trait-increasing allele
  if (this.line$beta >= 0) {
    # If beta >= 0, this means allele1 increases trait. Directly append.
    UPallele.freq.relict.vec <- c(UPallele.freq.relict.vec,allele1.freq.relict)
    UPallele.freq.nonrelict.vec <- c(UPallele.freq.nonrelict.vec,allele1.freq.nonrelict)
  } else {
    # Else if beta < 0, this means allele 0 increases the trait, Append 1 minus allele1.freq
    UPallele.freq.relict.vec <- c(UPallele.freq.relict.vec,(1-allele1.freq.relict))
    UPallele.freq.nonrelict.vec <- c(UPallele.freq.nonrelict.vec,(1-allele1.freq.nonrelict))
  }
}


output.table <- all.data[,c(1,3,4,5,6,7)]
output.table$allele1.freq.relict <- allele1.freq.relict.vec
output.table$allele1.freq.nonrelict <- allele1.freq.nonrelict.vec
output.table$UPallele.freq.relict <- UPallele.freq.relict.vec
output.table$UPallele.freq.nonrelict <- UPallele.freq.nonrelict.vec
output.table <- cbind(output.table,Fst$WEIR_AND_COCKERHAM_FST)

write.table(output.table,file=output.file,quote=F,col.names=T,row.names=F,sep="\t")
