# TODO: Add comment
# 
# Author: cshughes
###############################################################################
#get the PCA analysis done for proteins
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
opca = readRDS('ch_OvC_wStd_PCA_Analysis.rds')
#order the data by most contribution to variance
opcaV = as.data.frame(opca$rotation[,1:2])
opcaV$Gene = row.names(opca$rotation)
row.names(opcaV) = NULL
opcaV = opcaV[order(-abs(opcaV$PC1)),]


#build the RNA set to be clustered
setwd(dir="/Users/cshughes/Documents/projects/OvC/RNAexpression") #change this to whatever directory you have stored the data in
orna = readRDS('ch_OvC_RNA_processedProbes.rds')
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput")
#get out the expression data
eSet = as.data.frame(exprs(orna$eset))
#calculate RLE
eSet$med = apply(eSet,1,function(x) median(x,na.rm=TRUE))
eRLE = as.data.frame(apply(eSet[,1:55],2,function(x) x - eSet$med))
#annotate the data
gene.symbols <- getSYMBOL(row.names(eRLE), "hgu133plus2")
eRLE$Gene = gene.symbols
#remove probes with NA annotation
eSub = subset(eRLE, !is.na(eRLE$Gene))
#aggregate into single gene measurements
eAgg = aggregate(eSub[,1:55],by=list(eSub$Gene), median, na.rm=TRUE)
colnames(eAgg)[1] = 'Gene'



###############################################################################
###clustering heat maps by subtype
###############################################################################
#build the RNA set to be clustered based on the top 1000 genes from the PCA
gn = opcaV[1:610,3]
hSet = eAgg[eAgg$Gene %in% gn,]


#built the heatmap
pdf('ch_OvC_TMT10_RNA_MarkerSet_Top500_HeatMap.pdf')
#make the plot labels and boundaries
xLabels<- c(rep('ccc',25),rep('emc',14),rep('hgs',16))
mybreaks = seq(-2,2,by=0.05)
ColSideColors = c(rep(brewer.pal(6,'Accent')[2],25),rep(brewer.pal(6,'Accent')[3],14),rep(brewer.pal(6,'Accent')[1],16))
#make the correlation heatmap
heatmap.2(
		as.matrix(hSet[2:56]),
		col= colorRampPalette(brewer.pal(6,"RdBu"))(length(mybreaks)-1),
		symkey=TRUE,
		Rowv=TRUE,
		Colv=TRUE,
		dendrogram="both",
		breaks=mybreaks,
		#cellnote = round(ovCor,2),
		labRow = '',
		labCol = xLabels,
		las=2,
		ColSideColors=ColSideColors,
		## labels
		main='OvC Subtype Clustering',
		## color key
		key = TRUE,
		keysize = 1,
		density.info = "none",
		scale = "none",
		trace = "none",
		mar=c(8,8),
		cexRow=1.5,
		cexCol=1.5
)
dev.off()







