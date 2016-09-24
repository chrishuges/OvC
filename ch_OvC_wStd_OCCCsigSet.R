# TODO: Add comment
# 
# Author: cshughes
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/signatures/")
ccc<-read.table("./OCCC_signature.txt", header=TRUE, sep='\t', na.strings='')



###############################################################################
#look in the HGS vs CCC comparison first
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
phvc = readRDS('ch_OvC_wStd_Proteins_HGSvCCC.rds')

#bring in the protein data
proh = phvc
#bring in marker data
setwd(dir="/Users/cshughes/Documents/projects/OvC/signatures/")
up = read.table("./OCCC_signatureUP.txt", header=TRUE, sep='\t')
dn = read.table("./OCCC_signatureDN.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput/")
#combine into a single set
ccUP = as.character(unique(up[,3]))
ccDN = as.character(unique(dn[,3]))
ccTOT = c(ccUP,ccDN)
#only keep markers that are present in the data
markUPDN = ccTOT[which(ccTOT %in% proh$Gene)]
markUP = ccDN[which(ccDN %in% proh$Gene)]
markDN = ccUP[which(ccUP %in% proh$Gene)]

#make a plot of the data
proSD<-sd(proh$PROexp, na.rm=TRUE)
#sort out colors
lnCols<-brewer.pal(6,"Blues")
cols<-rev(brewer.pal(6,"RdBu"))

###make an initial plot with all points
pdf('ch_OvC_TMT10_wStd_Human_Proteins_OCCCsig_Volcano.pdf')
xCol = col2rgb(ifelse(proh$Gene %in% markDN, cols[1], ifelse(proh$Gene %in% markUP, cols[6],'gray80')))
xCex = ifelse(proh$Gene %in% markUPDN, 2, 1)
plot(proh$PROexp,
		-log10(proh$score),
		col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),
		pch=20,
		cex = xCex,
		ylab = '-log10(Adjusted p-value)',
		xlab = 'log2(serous/clear cell)',
		main = 'OvC signature comparison',
		xlim = c(-4,4),
		ylim = c(0,5)
)
box(lwd=3)
abline(v = -proSD, col=lnCols[6],lwd=2,lty=2)
abline(v = proSD, col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.05), col=lnCols[6],lwd=2,lty=2)
text(3,0.5,paste('n=',nrow(proh),sep=""),cex=1.25)
text(3,2.3,paste('p<0.05'),cex=1.25)

###make a second plot with just the marker points
xCol = ifelse(proh$Gene %in% markDN, cols[1], ifelse(proh$Gene %in% markUP, cols[6],'gray60'))
gnRM = !grepl('gray', xCol)
proh.s = proh[gnRM,]
xCol = col2rgb(ifelse(proh.s$Gene %in% markDN, cols[1], ifelse(proh.s$Gene %in% markUP, cols[6],'gray80')))
plot(proh.s$PROexp,
		-log10(proh.s$score),
		col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),
		pch=20,
		cex = 2,
		ylab = '-log10(Adjusted p-value)',
		xlab = 'log2(serous/clear cell)',
		main = 'OvC signature comparison',
		xlim = c(-4,4),
		ylim = c(0,5)
)
box(lwd=3)
abline(v = -proSD, col=lnCols[6],lwd=2,lty=2)
abline(v = proSD, col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.05), col=lnCols[6],lwd=2,lty=2)
text(3,0.5,paste('n=',sum(gnRM),sep=""),cex=1.25)
text(3,2.3,paste('p<0.05'),cex=1.25)
#text(proh.s$PROexp, -log10(proh.s$score), proh.s$Gene, cex=0.5)
dev.off()











###############################################################################
###clustering heat maps by the marker set
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
#merge replicates
z1 = rbind(a1n,a2n,a3n)
z2 = rbind(b1n,b2n,b3n)
#aggregate redundant peptides
agg1 = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8,a9)~Gene,data=z1,median,na.action=na.pass,na.rm=TRUE)
agg2 = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8,a9)~Gene,data=z2,median,na.action=na.pass,na.rm=TRUE)
#merge into a gene set
pro = merge(agg1,agg2,by='Gene')
#for output only
#proOut = pro[,c(1:4,11:13,5:7,14:16,8:10,17:19)]
#colnames(proOut) = c('Gene','t1811D','t1768A','t376D','t1842B','t1630B','t1863B','t1215B','t1210A','t1246A','t1303A','t1493A','t867A','t782A','t1683A','t1911A','t1195A','t1328A','t614A')
#saveRDS(proOut,'ch_OvC_wStd_RLEset.rds')
#calculate mean expression for the subtypes
pro$hgs1 = rowMeans(pro[,c(2:4,11:13)],na.rm=TRUE)
pro$ccc1 = rowMeans(pro[,c(5:7,14:16)],na.rm=TRUE)
pro$emc1 = rowMeans(pro[,c(8:10,17:19)],na.rm=TRUE)
#calculate the variance and order by it
pro$var = apply(pro[,20:22], 1, function(x) var(x, na.rm=TRUE))
proV = pro[order(-pro$var),]
#remove rows that have any NA values
proS = subset(proV, rowSums(is.na(proV[,2:19]))==0)





#specify the number of genes to keep for the clustering
proX = proS[proS$Gene %in% ccTOT,]
proXs = proX[1:75,]
#built the heatmap
#pdf('ch_OvC_TMT10_wStd_Human_Proteins_OCCCsig_HeatMap.pdf')
pdf('ch_test.pdf')
#make the plot labels and boundaries
xLabels<- c('hgs1','hgs2','hgs3','ccc1','ccc2','ccc3','emc1','emc2','emc3','hgs4','hgs5','hgs6','ccc4','ccc5','ccc6','emc4','emc5','emc6')
mybreaks = seq(-2,2,by=0.05)
ColSideColors = c(rep(brewer.pal(6,'Accent')[1],3),rep(brewer.pal(6,'Accent')[2],3),rep(brewer.pal(6,'Accent')[3],3),rep(brewer.pal(6,'Accent')[1],3),rep(brewer.pal(6,'Accent')[2],3),rep(brewer.pal(6,'Accent')[3],3))
#make the correlation heatmap
heatmap.2(
		as.matrix(proXs[2:19]),
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


ovCor = cor(proXs[2:19], use='pairwise.complete.obs', method='pearson')
#make the plot
pdf('ch_test2.pdf')
#pdf('ch_test.pdf')
#make the plot labels and boundaries
xLabels<- c('hgs1','hgs2','hgs3','ccc1','ccc2','ccc3','emc1','emc2','emc3','hgs4','hgs5','hgs6','ccc4','ccc5','ccc6','emc4','emc5','emc6')
mybreaks = seq(0,0.9,by=0.05)
ColSideColors = c(rep(brewer.pal(6,'Accent')[1],3),rep(brewer.pal(6,'Accent')[2],3),rep(brewer.pal(6,'Accent')[3],3),rep(brewer.pal(6,'Accent')[1],3),rep(brewer.pal(6,'Accent')[2],3),rep(brewer.pal(6,'Accent')[3],3))
#make the correlation heatmap
heatmap.2(
		ovCor,
		col= colorRampPalette(brewer.pal(9,"YlOrRd"))(length(mybreaks)-1),
		symkey=FALSE,
		Rowv=TRUE,
		Colv=TRUE,
		dendrogram="both",
		breaks=mybreaks,
		cellnote = round(ovCor,2),
		labRow = xLabels,
		labCol = xLabels,
		notecol = 'black',
		notecex = 0.75,
		colsep = 1:12,
		rowsep = 1:12,
		sepwidth = c(0.03, 0.03),
		sepcolor = 'white',
		ColSideColors=ColSideColors,
		## labels
		main='OvC Subtype Correlation',
		## color key
		key = TRUE,
		keysize = 1,
		density.info = "none",
		scale = "none",
		trace = "none",
		mar=c(8,8),
		cexRow=1.2,
		cexCol=1.2
)
dev.off()



###############################################################################
###how many in the top 500 genes that describe variance between HGSC and CCC
###############################################################################
#read in the data
pro = readRDS('ch_OvC_wStd_RLEset.rds')
#calculate mean expression for the subtypes
pro$hgs1 = rowMeans(pro[,c(2:4,11:13)],na.rm=TRUE)
pro$ccc1 = rowMeans(pro[,c(5:7,14:16)],na.rm=TRUE)
#calculate the variance and order by it
pro$var = apply(pro[,20:21], 1, function(x) var(x, na.rm=TRUE))
proV = pro[order(-pro$var),]
#remove rows that have any NA values
proS = subset(proV, rowSums(is.na(proV[,2:13]))==0)
#do a quick PCA to look
hCols = brewer.pal(6,'Accent')
pcCols = c(rep(hCols[1],6),rep(hCols[2],6))
pcPCH = c(rep(15,6),rep(19,6))
pca <- prcomp(t(proS[,2:13]))
pdf('ch_OvC_TMT10_wStd_Human_Proteins_CCCSig_PCA.pdf')
plot(pca$x[,1:2],
		col=pcCols,
		pch=pcPCH,
		cex = 2
)
box(lwd=3)
dev.off()

#keep top 500 genes based on variance
proX = proS[1:500,]
#specify the number of genes to keep for the clustering
gin = proX[proX$Gene %in% ccTOT,]


