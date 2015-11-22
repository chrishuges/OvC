# TODO: Add comment
# 
# Author: cshughes
###############################################################################

##########################################################
##make a plot for peptide metrics
##########################################################
#subset out peptides with greater than 50 peptides assigned
vE = pro
vE$pepSum = ifelse(vE$pepNum.x > vE$pepNum.y, vE$pepNum.x, vE$pepNum.y)
nPeps<-data.frame(table(vE[,10]))
up<-subset(nPeps, as.numeric(as.character(nPeps$Var1))>25)
dn<-subset(nPeps, as.numeric(as.character(nPeps$Var1))<=25)
#sum the number of proteins id with more than 50 peptides
x<-sum(up$Freq)
dn[26,1]<-'26'
dn[26,2]<-x
#
mainCol<-brewer.pal(6,"RdBu")
colA<-col2rgb(mainCol[6])
xLabels<-seq(1,26,1)
xValues<-seq(1,26,1)
#
mp<-barplot(dn[,2])
pdf('ch_OvC-TMT10_PeptideIDMetrics.pdf')
barplot(dn[,2],
		col=rgb(colA[1,],colA[2,],colA[3,],100,maxColorValue=255),
		main="Distribution of Peptide Identification Numbers", 
		xlab="Number of Unique Peptides",
		ylim=c(0,max(dn[,2])+50),
		xaxt="n"
)
box(lwd=3)
axis(1,mp,xLabels,las=2,lwd=2)
dev.off()


##################################################
#check the normalization scheme in terms of abundance
##################################################
#use the ov.psm data frame from the processing sheet
pdf('ch_AbundanceNormalizationMethods.pdf')
boxplot(ov.psm[,c(7:16)],las=2,main='raw abundances')
boxplot(log10(ov.psm[,c(7:16)]),las=2,main='log-scale raw abundances')
boxplot(log10(ov.psm[grepl('ecoli',ov.psm$Source),c(7:16)]),las=2,main='e-coli abundances log scale')
boxplot(log10(ov.psm[grepl('human',ov.psm$Source),c(7:16)]),las=2,main='human abundances log scale')
boxplot(log10(normalize.quantiles(as.matrix(ov.psm[,c(7:16)]))),las=2,main='quantile normalized abundances log scale')
boxplot(normalizeVSN(as.matrix(ov.psm[,c(7:16)])),las=2,main='limma VSN abundances')
boxplot(justvsn(as.matrix(ov.psm[,c(7:16)])),las=2,main='VSN abundances')
dev.off()

##################################################
#do some investigation of the peptides
##################################################
###look at the peptide distributions
#assign replicate rows
#make the plots
xnorm = as.matrix(ov.psm[ovA,c(5:24)])
pdf('ch_preNorm_peptideDistributions.pdf')
ov.n = names(ov.psm[,5:24])
for (i in 1:10) {
	hist(log10(xnorm[,i]),
			main = paste('Replicate ',ov.psm$Rep[1],' ',ov.n[i], sep=''),
			breaks = 200,
			col = brewer.pal(6,'RdBu')[6],
			xlim = c(0,3.5)
	)	
}
boxplot(log10(xnorm),
		col = brewer.pal(6,'RdBu')[2],
		las = 2)
boxplot(justvsn(xnorm),
		col = brewer.pal(6,'RdBu')[2],
		las = 2)
dev.off()


##################################################
#plot the peptide replicates
##################################################
pA = cbind(ov.pepA,ov.qA)
pB = cbind(ov.pepB,ov.qB)
#subset anything with a small FC
pAs = subset(pA, logFC > 1 | logFC < -1)
ovCorA = cor(pAs[,c(7:16)], use='pairwise.complete.obs', method='pearson')
pBs = subset(pB, logFC > 1 | logFC < -1)
ovCorB = cor(pBs[,c(7:16)], use='pairwise.complete.obs', method='pearson')

pdf('ch_OvC_subsetReplicates_heatNorm.pdf')
#make the plot labels and boundaries
xLabels<- names(ovCor)
mybreaks = seq(0.2,1,by=0.05) 
#make the correlation heatmap
ovCorA[lower.tri(ovCorA)] <- ovCorB[lower.tri(ovCorB)]
heatmap.2(
		ovCorA,
		col= colorRampPalette(brewer.pal(9,"YlOrRd"))(length(mybreaks)-1),
		symkey=FALSE,
		Rowv=FALSE,
		Colv=FALSE,
		#hclust=hclustfunc,
		#distfun=distfunc,
		#na.color='black',
		#na.rm=TRUE,
		#symm=FALSE,
		dendrogram="none",
		breaks=mybreaks,
		cellnote = round(ovCorA,2),
		labRow = xLabels,
		labCol = xLabels,
		notecol = 'black',
		notecex = 0.7,
		colsep = 1:10,
		rowsep = 1:10,
		sepwidth = c(0.03, 0.03),
		sepcolor = 'white',
		#ColSideColors=ColSideColors,
		## labels
		main='Peptide reproducibility correlation map',
		#xlab='Time Point',
		#ylab='Protein',
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


##########################################################
##make a plot for interesting candidates
##########################################################

setwd(dir="/Users/cshughes/Documents/projects/OvC/")
hunt<-read.table("./Huntsman_GeneExpression.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/Routput/")
#need to reshape the hunt frame
hunt2 = melt(hunt, id='Gene')



pdf('ch_HuntsmanExpression.pdf')
barchart(hunt2[,3]~hunt2[,1],
		beside=FALSE,
		groups=hunt2[,2],
		las=2)
dev.off()

##########################################################
##make volcano plots of expression variance
##########################################################
#make some colors
lnCols<-brewer.pal(6,"Blues")
cols<-rev(brewer.pal(6,"RdBu"))
vE = pro
#get the list of candidate genes
setwd(dir="/Users/cshughes/Documents/projects/OvC/")
genesCC = read.table("./Schwartz_ClearCellTable_geneList.txt", header=FALSE, sep='\t')
genesS = read.table("./Schwartz_SerousTable_geneList.txt", header=FALSE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/Routput/")

geneCC = as.character(genesCC[,1])
geneS = as.character(genesS[,1])
geneTot = c(geneCC,geneS)

#make the plotting function 
pdf('ch_OvC_TMT10_expression_volcano_CC.pdf')
proSD<-sd(vE[,4], na.rm=TRUE)
xCol = col2rgb(ifelse(vE[,2] %in% geneCC, cols[1], ifelse(vE[,2] %in% geneS, cols[6],'gray30')))
xCex = ifelse(vE[,2] %in% geneTot, 2, 1)
#xCol = col2rgb(ifelse(-log10(vE[,5]) > -log10(0.01),ifelse(vE[,4] > proSD, cols[1],ifelse(vE[,4] < -proSD, cols[1], 'gray30')),'gray30'))
#plot with all points
plot(vE[,4],
			-log10(vE[,5]),
			col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),
			pch=20,
			cex = xCex,
			ylab = '-log10(Adjusted p-value)',
			xlab = 'log2(serous/clear cell)',
			main = 'OvC type comparison',
			xlim = c(-5,5)
	)
	box(lwd=3)
	abline(v = -proSD, col=lnCols[6],lwd=2,lty=2)
	abline(v = proSD, col=lnCols[6],lwd=2,lty=2)
	abline(h = -log10(0.05), col=lnCols[6],lwd=2,lty=2)
	abline(h = -log10(0.01), col=lnCols[6],lwd=2,lty=2)
	text(4,0.5,paste('n=',nrow(pro),sep=""),cex=1.25)
	text(4,1.6,paste('p<0.05'),cex=1.25)
	text(4,2.3,paste('p<0.01'),cex=1.25)
#plot with just subset points
xCol = ifelse(vE[,2] %in% geneCC, cols[1], ifelse(vE[,2] %in% geneS, cols[6],'gray30'))
gnRM = !grepl('gray', xCol)
vB = vE[gnRM,]
xCol = col2rgb(ifelse(vB[,2] %in% geneCC, cols[1], ifelse(vB[,2] %in% geneS, cols[6],'gray30')))
plot(vB[,4],
			-log10(vB[,5]),
			col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),
			pch=20,
			cex = 2,
			ylab = '-log10(Adjusted p-value)',
			xlab = 'log2(serous/clear cell)',
			main = 'OvC type comparison',
			xlim = c(-5,5)
	)
	box(lwd=3)
	abline(v = -proSD, col=lnCols[6],lwd=2,lty=2)
	abline(v = proSD, col=lnCols[6],lwd=2,lty=2)
	abline(h = -log10(0.05), col=lnCols[6],lwd=2,lty=2)
	abline(h = -log10(0.01), col=lnCols[6],lwd=2,lty=2)
	text(4,0.5,paste('n=',sum(gnRM),sep=""),cex=1.25)
	text(4,1.6,paste('p<0.05'),cex=1.25)
	text(4,2.3,paste('p<0.01'),cex=1.25)
dev.off()

##########################################################
##make a heat map just for markers
##########################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/")
genesCC = read.table("./Schwartz_ClearCellTable_geneList.txt", header=FALSE, sep='\t')
genesS = read.table("./Schwartz_SerousTable_geneList.txt", header=FALSE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/Routput/")
#make them into a searchable list
geneCC = as.character(genesCC[,1])
geneS = as.character(genesS[,1])
geneTot = c(geneCC,geneS)
#get the peptide data
pA = cbind(ov.pepA,ov.qA)
pB = cbind(ov.pepB,ov.qB)
#subset anything with a small FC
pAs = subset(pA, Gene %in% geneTot & pVal < 0.05)
ovCorA = cor(pAs[,c(7:16)], use='pairwise.complete.obs', method='pearson')
pBs = subset(pB, Gene %in% geneTot & pVal < 0.05)
ovCorB = cor(pBs[,c(7:16)], use='pairwise.complete.obs', method='pearson')

pdf('ch_OvC_Markers_heat.pdf')
#make the plot labels and boundaries
xLabels<- names(ovCor)
mybreaks = seq(0.2,1,by=0.05) 
#make the correlation heatmap
ovCorA[lower.tri(ovCorA)] <- ovCorB[lower.tri(ovCorB)]
heatmap.2(
		ovCorA,
		col= colorRampPalette(brewer.pal(9,"YlOrRd"))(length(mybreaks)-1),
		symkey=FALSE,
		Rowv=FALSE,
		Colv=FALSE,
		#hclust=hclustfunc,
		#distfun=distfunc,
		#na.color='black',
		#na.rm=TRUE,
		#symm=FALSE,
		dendrogram="none",
		breaks=mybreaks,
		cellnote = round(ovCorA,2),
		labRow = xLabels,
		labCol = xLabels,
		notecol = 'black',
		notecex = 0.7,
		colsep = 1:10,
		rowsep = 1:10,
		sepwidth = c(0.03, 0.03),
		sepcolor = 'white',
		#ColSideColors=ColSideColors,
		## labels
		main='Peptide reproducibility correlation map',
		#xlab='Time Point',
		#ylab='Protein',
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







