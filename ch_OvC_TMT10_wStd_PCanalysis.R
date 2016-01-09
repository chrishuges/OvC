# TODO: Add comment
# 
# Author: cshughes
###############################################################################
#read in the data
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
a1 = readRDS('ch_OvC_wStd_processedPeptides_a1.rds')
a2 = readRDS('ch_OvC_wStd_processedPeptides_a2.rds')
a3 = readRDS('ch_OvC_wStd_processedPeptides_a3.rds')
b1 = readRDS('ch_OvC_wStd_processedPeptides_b1.rds')
b2 = readRDS('ch_OvC_wStd_processedPeptides_b2.rds')
b3 = readRDS('ch_OvC_wStd_processedPeptides_b3.rds')
#use this function for normalization
pepNorm = function(x,...){
	#mark the spike-in control peptides
	spikeins = grepl('EColi', x$Organism)
	#extract the expression values
	xnorm = as.matrix(x[,7:16])
	#make a noise model based on the spike ins
	spfit = vsn2(xnorm[spikeins,],lts.quantile=1)
	meanSdPlot(spfit)
	#apply it to the human peptides
	nkid = predict(spfit, newdata=xnorm)
	#extract only human expression values
	nkid.h = as.data.frame(nkid[!spikeins,])
	#normalize the values to the internal standard sample
	vnorm = apply(nkid.h[,1:9],2, function(x) x - nkid.h$aStd)
	snorm = scale(vnorm, center=TRUE,scale=TRUE)
	#recombine the replicates
	vnormFull = as.data.frame(cbind(x[!spikeins,c(1:6,17)],snorm))
	#output the data
	return(vnormFull)
}
#apply the function
a1n = pepNorm(a1)
a2n = pepNorm(a2)
a3n = pepNorm(a3)
b1n = pepNorm(b1)
b2n = pepNorm(b2)
b3n = pepNorm(b3)
#do a quick check of the data quality
intPlot = list(a1n,a2n,a3n,b1n,b2n,b3n)
pdf('ch_OvC_TMT10_wStd_Human_Peptides_Densities.pdf')
for (i in 1:6){
	z = as.data.frame(intPlot[[i]])
	#snorm = scale(z[,8:16],center=TRUE,scale=FALSE)
	plotDensities(z[,8:16],
		legend = FALSE,
		main = 'Expression Signal Density per Sample')
abline(v=0,col='red',lty=2,lwd=2)
}
dev.off()
#looks good after centering, implement into the pepNorm function
###############################################################################
#aggregate and build a protein set
###############################################################################
#merge replicates
z1 = rbind(a1n,a2n,a3n)
z2 = rbind(b1n,b2n,b3n)
#aggregate redundant peptides
agg1 = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8,a9)~Gene,data=z1,median,na.action=na.pass,na.rm=TRUE)
agg2 = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8,a9)~Gene,data=z2,median,na.action=na.pass,na.rm=TRUE)
#merge into a gene set
pro = merge(agg1,agg2,by='Gene')
#calculate mean expression for the subtypes
pro$hgs1 = rowMeans(pro[,c(2:4,11:13)],na.rm=TRUE)
pro$ccc1 = rowMeans(pro[,c(5:7,14:16)],na.rm=TRUE)
pro$emc1 = rowMeans(pro[,c(8:10,17:19)],na.rm=TRUE)
#calculate the variance and order by it
pro$var = apply(pro[,20:22], 1, function(x) var(x, na.rm=TRUE))
proV = pro[order(-pro$var),]
#remove rows that have any NA values
proS = subset(proV, rowSums(is.na(proV[,2:19]))==0)
#subset the top N genes by variance
proX = proS[1:500,]

###############################################################################
###PCA stuffs
###############################################################################
x = proX[,2:19]
row.names(x) = proX$Gene
hCols = brewer.pal(6,'Accent')
#pcCols = c(rep(hCols[1],3),rep(hCols[2],3),rep(hCols[3],3),rep(hCols[4],3),rep(hCols[5],3),rep(hCols[6],3))
pcCols = c(rep(hCols[1],6),rep(hCols[2],6),rep(hCols[3],6))
pcPCH = c(rep(15,3),rep(19,3),rep(8,3),rep(12,3),rep(2,3),rep(11,3))
pca <- prcomp(t(x))
#pca <- prcomp(t(xSub),scale=TRUE,center=TRUE)
plot(pca$x[,1:2],
		col=pcCols,
		pch=pcPCH)

###############################################################################
###clustering heat maps by subtype
###############################################################################
#specify the number of genes to keep for the clustering
proX = proS[1:6000,]
#built the heatmap
pdf('ch_OvC_TMT10_wStd_Human_Proteins_Subtypes_Top6000_HeatMap.pdf')
#make the plot labels and boundaries
xLabels<- c('hgs1','hgs2','hgs3','ccc1','ccc2','ccc3','emc1','emc2','emc3','hgs4','hgs5','hgs6','ccc4','ccc5','ccc6','emc4','emc5','emc6')
mybreaks = seq(-2,2,by=0.05)
ColSideColors = c(rep(brewer.pal(6,'Accent')[1],3),rep(brewer.pal(6,'Accent')[2],3),rep(brewer.pal(6,'Accent')[3],3),rep(brewer.pal(6,'Accent')[1],3),rep(brewer.pal(6,'Accent')[2],3),rep(brewer.pal(6,'Accent')[3],3))
#make the correlation heatmap
heatmap.2(
		as.matrix(proX[2:19]),
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











proX = proS[1:500,]
pOut = proX[proX$hgs1< -1 & proX$ccc1>0.5,c(1,20:23)]
#write out the data into a format usable in excel
write.table(pOut,'test.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)


























