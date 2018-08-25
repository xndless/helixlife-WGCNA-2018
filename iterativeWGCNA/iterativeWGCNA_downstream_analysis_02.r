library(WGCNA)

args <- commandArgs()

MEs <- read.table(args[6],header=T,row.names=1)
MMs <- read.table(args[7],header=T)
#MEs <- read.table("final-eigengenes.txt",header=T,row.names=1)
#MMs <- read.table("final-membership.txt",header=T)
expr <- read.table(args[8],header=T,check.names=F,row.names=1)

dir.create(args[9])
setwd(args[9])
#dir.create("02.relationship_analysis")
#setwd("02.relationship_analysis")

for(i in 1:(nrow(MEs)-1)) {
which.module <- rownames(MEs)[i]
pdf(file=paste("expression_ME_",which.module,".pdf",sep=""),25,10)
gene <- MMs[MMs$Module==which.module,]$Gene
expr1 <- expr[match(gene,rownames(expr)),]
ME <- as.matrix(MEs[which.module,])
colnames(ME) <- colnames(expr1)
layout(matrix(c(1,2)),heights=c(1.5,3))
par(mar=c(0.3,8.5,3,5.5))
plotMat(t(scale(t(expr1))),
	nrgcols=30,rlabels=F,rcols=rainbow(nrow(MEs))[i],
	main=paste(which.module),cex.main=1)
par(mar=c(5,4,0,1))
barplot(ME,col=rainbow(nrow(MEs))[i],main="",
	cex.names=1,cex.axis=1,ylab="module eigengene",las=3)
dev.off()
}
