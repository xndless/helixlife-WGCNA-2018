library(WGCNA)

args <- commandArgs()

load("wgcna-blocks.RData")
load(paste(args[6],"-TOM-block.1.RData",sep=""))

eval(parse(text=paste(args[6]," <- as.matrix(TOM)",sep="")))

dir.create("module_result")
setwd("module_result")

modGenes <- rownames(expression[blocks$goodGenes,])
nTop <- length(modGenes)

eval(parse(text=paste("modTOM <- ",args[6],sep="")))
IMConn <- softConnectivity(t(expression[blocks$goodGenes,]))
top <- (rank(-IMConn) <= nTop)
cyt1 <- exportNetworkToCytoscape(modTOM[top,top],
		edgeFile=paste("CytoscapeInput-edges-",args[6],".txt",sep=""),
		nodeFile=paste("CytoscapeInput-nodes-",args[6],".txt",sep=""),
		weighted=TRUE,
		threshold=0.00,
		nodeNames=modGenes[top],
		altNodeNames=modGenes[top])
out <- cbind(modGenes[top],IMConn[top])
colnames(out) <- c("gene","connectivity")
out <- out[order(as.numeric(out[,2]),decreasing=T),]
write.table(out,paste(args[6],"-module-gene.txt",sep=""),sep="\t",quote=F,row.names=F)
