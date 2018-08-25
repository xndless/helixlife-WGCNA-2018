library(WGCNA)
library(pheatmap)

args <- commandArgs()

MEs0 <- read.table(args[6],header=T,row.names=1)
MEs <- t(MEs0)

datExpr0 <- read.table(args[7],header=T,row.names=1,check.names=F,sep="\t")
datTraits0 <- read.table(args[8],header=T,row.names=1,check.names=F,sep="\t")

dir.create(args[9])
setwd(args[9])

MEDiss <- 1-cor(MEs)
METree <- hclust(as.dist(MEDiss),method="average")
pdf(file="modules_cluster_tree.pdf",width=7,height=5)
plot(METree,main="Clustering of module eigengenes",xlab="",sub="")
dev.off()

moduleCor <- corAndPvalue(MEs,use="p")
rowLabels <- colnames(MEs)
textMatrix <- paste(signif(moduleCor$cor,2),
				"\n(",signif(moduleCor$p,1),")",sep="")
dim(textMatrix) <- dim(moduleCor$cor)
pdf(file="modules_relationships.pdf",15,8)
par(mar=c(5,5,1,2))
labeledHeatmap(Matrix=moduleCor$cor,
			textMatrix=textMatrix,
			xLabels=rowLabels,yLabels=rowLabels,
			colorLabels=TRUE,colors=blueWhiteRed(50),
			setStdMargins=FALSE,xLabelsAngle=90,zlim=c(-1,1))
dev.off()

text <- paste("cor=",round(moduleCor$cor,4),
			";p-value=",round(moduleCor$p,4),sep="")
dim(text) <- dim(moduleCor$cor)
rownames(text) <- rowLabels
colnames(text) <- rowLabels
text <- cbind(rownames(text),text)
colnames(text)[1] <- "modules"
write.table(text,file="modules_relationships.xls",
		quote=F,sep="\t",row.names=F)

###################################################################################################

datExpr <- t(datExpr0)
sample_cor <- cor(t(datExpr),t(datExpr),
				use='pairwise.complete.obs')
moduleSampleCor <- cor(MEs,sample_cor,use="p")
nSamples <- nrow(datExpr)
moduleSamplePvalue <- corPvalueStudent(moduleSampleCor,nSamples)
textMatrix <- paste(signif(moduleSampleCor,2),
				"\n(",signif(moduleSamplePvalue,1),")",sep="")
dim(textMatrix) <- dim(moduleSampleCor)
rowLabels <- colnames(MEs)
pdf(file="modules_samples_relationships.pdf",20,6)
par(mar=c(5,6,1,1))
labeledHeatmap(Matrix=moduleSampleCor,
			xLabels=colnames(sample_cor),
			yLabels=rowLabels,ySymbols=names(MEs),
			colorLabels=TRUE,colors=blueWhiteRed(50),
			setStdMargins=FALSE,xLabelsAngle=90,zlim=c(-1,1))
dev.off()

text <- paste("cor=",round(moduleSampleCor,4),
			";p-value=",round(moduleSamplePvalue,4),sep="")
dim(text) <- dim(moduleSampleCor)
rownames(text) <- rownames(moduleSampleCor)
colnames(text) <- colnames(moduleSampleCor)
text <- cbind(rownames(text),text)
colnames(text)[1] <- "modules"
write.table(text,file="modules_samples_relationships.xls",
		quote=F,sep="\t",row.names=F)

###################################################################################################

datTraits <- datTraits0[rownames(MEs),]
moduleTraitCor <- cor(MEs,datTraits,use="p")
nSamples <- nrow(datExpr)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
textMatrix <- paste(signif(moduleTraitCor,2),
				"\n(",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix) <- dim(moduleTraitCor)
rowLabels <- colnames(MEs)
pdf(file="modules_traits_relationships.pdf",15,8)
par(mar=c(7,12,1,2))
labeledHeatmap(Matrix=moduleTraitCor,
			textMatrix=textMatrix,
			xLabels=colnames(datTraits),
			yLabels=rowLabels,ySymbols=names(MEs),
			colorLabels=TRUE,colors=blueWhiteRed(50),
			setStdMargins=FALSE,xLabelsAngle=90,zlim=c(-1,1))
dev.off()

text <- paste("cor=",round(moduleTraitCor,4),
			";p-value=",round(moduleTraitPvalue,4),sep="")
dim(text) <- dim(moduleTraitCor)
rownames(text) <- rownames(moduleTraitCor)
colnames(text) <- colnames(moduleTraitCor)
text <- cbind(rownames(text),text)
colnames(text)[1] <- "modules"
write.table(text,file="modules_traits_relationships.xls",
		quote=F,sep="\t",row.names=F)
