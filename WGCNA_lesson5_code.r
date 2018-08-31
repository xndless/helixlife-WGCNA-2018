dir.create("helix-WGCNA-lesson5")
setwd("helix-WGCNA-lesson5")
fileUrl <- "https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/JMiller/metaAnalysisFiles.zip"
download.file(fileUrl,destfile="metaAnalysisFiles.zip")
unzip("metaAnalysisFiles.zip")
list.files()
## [1] "__MACOSX" "metaAnalysisFiles" "metaAnalysisFiles.zip"

setwd("metaAnalysisFiles")
load("metaAnalysisData.RData")
ls()
## [1] "datExprA1" "datExprA2" "datExprB1" "datExprB2" "genesA"    "genesI"   
## [7] "probesA"   "probesI"

###################################################################################################

library(WGCNA)
library(reshape2)

commonProbesA <- intersect(rownames(datExprA1),rownames(datExprA2))
datExprA1p <- datExprA1[commonProbesA,]
datExprA2p <- datExprA2[commonProbesA,]

softPower <- 10
rankExprA1 <- rank(rowMeans(datExprA1p))
rankExprA2 <- rank(rowMeans(datExprA2p))
random5000 <- sample(commonProbesA,5000)
rankConnA1 <- rank(softConnectivity(t(datExprA1p[random5000,]),
				type="signed",power=softPower))
rankConnA2 <- rank(softConnectivity(t(datExprA2p[random5000,]),
				type="signed",power=softPower))

pdf("01.generalNetworkProperties.pdf",height=5,width=10)
par(mfrow=c(1,2))
verboseScatterplot(rankExprA1,rankExprA2,xlab="Ranked Expression (A1)",
				ylab="Ranked Expression (A2)")
verboseScatterplot(rankConnA1,rankConnA2,xlab="Ranked Connectivity (A1)",
				ylab="Ranked Connectivity (A2)")
dev.off()

###################################################################################################

networkType <- "signed"
netA1 <- blockwiseModules(t(datExprA1p),power=softPower,
				networkType=networkType,numericLabels=TRUE,
				mergeCutHeight=0.25,minModuleSize=30,
				maxBlockSize=30000,saveTOMs=TRUE,
				saveTOMFileBase="WGCNA_TOM_A1",verbose=5)
moduleLabelsA1 <- netA1$colors
moduleColorsA1 <- labels2colors(moduleLabelsA1)

netA2 <- blockwiseModules(t(datExprA2p),power=softPower,
                networkType=networkType,numericLabels=TRUE,
                mergeCutHeight=0.25,minModuleSize=30,
                maxBlockSize=30000,saveTOMs=TRUE,
                saveTOMFileBase="WGCNA_TOM_A2",verbose=5)
moduleLabelsA2 <- netA2$colors
moduleColorsA2 <- labels2colors(moduleLabelsA2)

pdf(file="02.comparison_modules_color.pdf",width=12,height=5)
plotDendroAndColors(netA1$dendrograms[[1]],
				cbind(moduleColorsA1[netA1$blockGenes[[1]]],
					moduleColorsA2[netA2$blockGenes[[1]]]),
				c("Modules A1","Modules A2"),
				dendroLabels=FALSE,addGuide=TRUE)
dev.off()

###################################################################################################

multiExpr <- list(A1=list(data=t(datExprA1p)),A2=list(data=t(datExprA2p)))
multiColor <- list(A1=moduleColorsA1[netA1$blockGenes[[1]]])
mp <- modulePreservation(multiExpr,multiColor,
				referenceNetworks=1,verbose=5,
				networkType="signed",
				nPermutations=30)
stats <- mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2
stats[order(-stats[,2]),c(1:2)]
##              moduleSize Zsummary.pres
## brown              1000      54.22666
## blue               1000      41.47798
## black               638      36.43791
## turquoise          1000      36.06894
## yellow             1000      33.92899

###################################################################################################

MEs0_A1 <- moduleEigengenes(t(datExprA1p)[,netA1$goodGenes],
						moduleColorsA1[netA1$blockGenes[[1]]])$eigengenes
MEs_A1 <- orderMEs(MEs0_A1)
geneModuleMembership1 <- signedKME(t(datExprA1p)[,netA1$goodGenes],MEs_A1)
name <- colnames(MEs_A1)
colnames(geneModuleMembership1) <- paste("PC",name,".cor",sep="")
geneModuleMembership1 <-
geneModuleMembership1[,order(colnames(geneModuleMembership1))]
MMPvalue1 <- corPvalueStudent(as.matrix(geneModuleMembership1),
			dim(datExprA1p)[[2]])
colnames(MMPvalue1) <- paste("PC",name[order(name)],".pval",sep="")
Gene1 <- rownames(datExprA1p)[netA1$goodGenes]
Gene1 <- as.data.frame(Gene1)
kMEtable1 <- cbind(Gene1,Gene1,moduleColorsA1[netA1$blockGenes[[1]]])
for (i in 1:nrow(t(MEs_A1))){
kMEtable1 <- cbind(kMEtable1,geneModuleMembership1[,i],MMPvalue1[,i])
}
colnames(kMEtable1) <- c("PSID","Gene","Module",
		sort(c(colnames(geneModuleMembership1),colnames(MMPvalue1))))
kMEtable1 <- kMEtable1[kMEtable1$Module!="grey",]
write.csv(na.omit(kMEtable1),"03.kMEtable_A1.csv",row.names=FALSE)

MEs0_A2 <- moduleEigengenes(t(datExprA2p)[,netA2$goodGenes],
						moduleColorsA2[netA2$blockGenes[[1]]])$eigengenes
MEs_A2 <- orderMEs(MEs0_A2)
geneModuleMembership2 <- signedKME(t(datExprA2p)[,netA2$goodGenes],MEs_A2)
name <- colnames(MEs_A2)
colnames(geneModuleMembership2) <- paste("PC",name,".cor",sep="")
geneModuleMembership2 <-
geneModuleMembership2[,order(colnames(geneModuleMembership2))]
MMPvalue2 <- corPvalueStudent(as.matrix(geneModuleMembership2),
			dim(datExprA2p)[[2]])
colnames(MMPvalue2) <- paste("PC",name[order(name)],".pval",sep="")
Gene2 <- rownames(datExprA2p)[netA2$goodGenes]
Gene2 <- as.data.frame(Gene2)
kMEtable2 <- cbind(Gene2,Gene2,moduleColorsA2[netA2$blockGenes[[1]]])
for (i in 1:nrow(t(MEs_A2))){
kMEtable2 <- cbind(kMEtable2,geneModuleMembership2[,i],MMPvalue2[,i])
}
colnames(kMEtable2) <- c("PSID","Gene","Module",
		sort(c(colnames(geneModuleMembership2),colnames(MMPvalue2))))
kMEtable2 <- kMEtable2[kMEtable2$Module!="grey",]
write.csv(na.omit(kMEtable2),"03.kMEtable_A2.csv",row.names=FALSE)

enrichments_A1_A2 <- userListEnrichment(rownames(datExprA1p)[netA1$goodGenes],
				moduleColorsA1[netA1$blockGenes[[1]]],
				"03.kMEtable_A2.csv",
				"","03.significant_enrichment_A1_A2.csv")
enrichments_A1_A2$pValues[,2] <- sub("__","",enrichments_A1_A2$pValues[,2])
head(enrichments_A1_A2$pValues)
##   InputCategories UserDefinedCategories Type NumOverlap      Pvalues
## 1           black                 black User        148 1.771052e-53
## 2           black                  blue User        260 6.720260e-61
## 3           black                 brown User          2 1.000000e+00
## 4           black                  cyan User          2 9.997414e-01
## 5           black                 green User          3 1.000000e+00
## 6           black           greenyellow User          1 1.000000e+00
##   CorrectedPvalues
## 1     6.021578e-51
## 2     2.284888e-58
## 3     1.000000e+00
## 4     1.000000e+00
## 5     1.000000e+00
## 6     1.000000e+00

N_A1_A2 <- dcast(enrichments_A1_A2$pValues,
			InputCategories~UserDefinedCategories,
			value.var="NumOverlap")
rownames(N_A1_A2) <- N_A1_A2[,1]
N_A1_A2 <- N_A1_A2[,-1]
N_A1_A2[1:5,1:5]
##       black blue brown cyan green
## black   148  260     2    2     3
## blue     21    9  1460    3   114
## brown     9    7   633   15   953
## cyan      6   48     0   74     5
## green     8    9   297   18   140

P_A1_A2 <- dcast(enrichments_A1_A2$pValues,
			InputCategories~UserDefinedCategories,
			value.var="Pvalues")
rownames(P_A1_A2) <- P_A1_A2[,1]
P_A1_A2 <- P_A1_A2[,-1]

CP_A1_A2 <- dcast(enrichments_A1_A2$pValues,
			InputCategories~UserDefinedCategories,
			value.var="CorrectedPvalues")
rownames(CP_A1_A2) <- CP_A1_A2[,1]
CP_A1_A2 <- CP_A1_A2[,-1]
LCP_A1_A2 <- -log10(CP_A1_A2)
LCP_A1_A2[LCP_A1_A2>=50] <- 50

textMatrix <- paste(as.matrix(N_A1_A2),
				"\n(",as.matrix(signif(CP_A1_A2,1)),")",sep="")

pdf(file="03.internetwork_modules_relationships.pdf",15,10)
par(mar=c(6,6,1,2))
labeledHeatmap(Matrix=t(LCP_A1_A2),
			textMatrix=textMatrix,
			xLabels=colnames(t(LCP_A1_A2)),yLabels=rownames(t(LCP_A1_A2)),
			xSymbols=paste("ME",colnames(t(LCP_A1_A2)),sep=""),
			ySymbols=paste("ME",rownames(t(LCP_A1_A2)),sep=""),
			colorLabels=TRUE,colors=blueWhiteRed(100)[51:100],
			setStdMargins=FALSE,xLabelsAngle=45,zlim=c(0,50))
dev.off()

###################################################################################################

dir.create("04.kMEtable1_vs_kMEtable2")
topGenesKME <- NULL
topNames <- NULL
for(i in 1:ncol(geneModuleMembership1)){
for(j in 1:ncol(geneModuleMembership2)){
M1 <- sub("PCME","",names(geneModuleMembership1)[i])
M1 <- sub(".cor","",M1)
M2 <- sub("PCME","",names(geneModuleMembership2)[j])
M2 <- sub(".cor","",M2)
inMod1 <- moduleColorsA1[netA1$blockGenes[[1]]]==M1
inMod2 <- moduleColorsA2[netA2$blockGenes[[1]]]==M2
inModp <- intersect(rownames(geneModuleMembership1)[inMod1],
rownames(geneModuleMembership2)[inMod2])

if(length(inModp)>10){
fn <- paste("04.kMEtable1_",M1,"_vs_kMEtable2_",M2,".pdf",sep="")
pdf(paste("04.kMEtable1_vs_kMEtable2/",fn,sep=""),height=8,width=8)
verboseScatterplot(geneModuleMembership1[inModp,i],
geneModuleMembership2[inModp,j],
main=paste("A1_",M1,"_vs_A2_",M2,sep=""),
xlab=paste("kME in A1 ",M1,sep=""),
ylab=paste("kME in A2 ",M2,sep=""))
dev.off()

kMErank1 <- rank(-geneModuleMembership1[inModp,i])
kMErank2 <- rank(-geneModuleMembership2[inModp,j])
maxKMErank <- rank(apply(cbind(kMErank1,kMErank2+.00001),1,max))
topGenesKME <- cbind(topGenesKME,inModp[maxKMErank<=10])
comp <- paste("A1_",M1,"_vs_A2_",M2,sep="")
topNames <- cbind(topNames,comp)
}
}
}
colnames(topGenesKME) <- topNames
write.table(topGenesKME,
	"04.shared_top10_hub_internetwork_modules.xls",
	quote=F,row.names=F,sep="\t")
