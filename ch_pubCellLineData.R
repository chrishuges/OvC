# TODO: using published data from The Cancer Cell Line Encyclopedia enables predictive modeling of anticancer drug sensitivity 
# 
# Author: cshughes
###############################################################################
#get the data
setwd(dir="/Users/cshughes/Documents/projects/OvC/pubData/")
ccle = read.table("./mRNA_expression_cshughes_expression_subset.txt", header=TRUE, sep='\t')
ccle.all = read.table("./CCLE_Expression_Entrez_2012-09-29.txt", header=TRUE, sep='\t')
ccle.ov = read.table("./CCLE_Expression_OVonly_2012-09-29.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/Routput/")


##########################################################
##compare with the tissue data
##########################################################
#use the pro data frame from the processing sheet
tis = pro[,c(1:4,7)]
#make new colnames for ccle
colnames(ccle)[2] = 'Gene'
#make a mean expression column
tis$meanFC = rowMeans(tis[,4:5], na.rm=TRUE)
ccle$meanSER = rowMeans(ccle[,3:4], na.rm=TRUE)
ccle$meanCC = rowMeans(ccle[,5:7], na.rm=TRUE)
#make a fold change column for the gene expression data
ccle$meanFC = log2(ccle$meanSER/ccle$meanCC)
#merge the data by gene name
tis.ccle = merge(tis,ccle,by='Gene')
#extract the expression data
ovex = tis.ccle[,c(1,6,15)]
#correlation with the protein data is 0.186!!! really shitty...but is expected to be shit

##########################################################
##compare with the cell line data
##########################################################
#use the pro data frame from the processing sheet
tis = proCell[,c(1:4)]
#make new colnames for ccle
colnames(ccle)[2] = 'Gene'
#make a mean expression column
ccle$meanSER = rowMeans(ccle[,3:4], na.rm=TRUE)
ccle$meanCC = rowMeans(ccle[,5:7], na.rm=TRUE)
#make a fold change column for the gene expression data
ccle$meanFC = log2(ccle$meanSER/ccle$meanCC)
#merge the data by gene name
tis.ccle = merge(tis,ccle,by='Gene')
#extract the expression data
ovex = tis.ccle[,c(1,4,13)]
#correlation with the protein data is 0.186!!! really shitty...but is expected to be shit
cor(ovex[,2],ovex[,3], use='pairwise.complete.obs')

##################################################
#plot the correlation of the two sets with the gene expression data
##################################################
#setup the tissue data
tis = pro[,c(1:4,7)]
tis$meanFC = rowMeans(tis[,4:5], na.rm=TRUE)
#setup the gene expression data
colnames(ccle)[2] = 'Gene'
ccle$meanSER = rowMeans(ccle[,3:4], na.rm=TRUE)
ccle$meanCC = rowMeans(ccle[,5:7], na.rm=TRUE)
ccle$meanFC = log2(ccle$meanSER/ccle$meanCC)
#setup the cell line data
cell = proCell[,c(1:4)]
#merge them all together
tis.cell = merge(tis,cell,by='Gene')
tis.cell.ccle = merge(tis.cell,ccle,by='Gene')
#rework the data frame
ovCor = cor(as.matrix(tis.cell.ccle[,c(6,9,18)]), use='pairwise.complete.obs', method='pearson')
#subset anything with a small FC
#pBs = subset(pB, logFC > 1 | logFC < -1)
#ovCorB = cor(pBs[,c(7:16)], use='pairwise.complete.obs', method='pearson')

pdf('ch_OvCall_ProteinandGeneExpression_Heat.pdf')
#make the plot labels and boundaries
xLabels<- c('Tissue','Cell Line','Gene')
mybreaks = seq(0.2,1,by=0.05) 
#make the correlation heatmap
heatmap.2(
		ovCor,
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
		cellnote = round(ovCor,2),
		labRow = xLabels,
		labCol = xLabels,
		notecol = 'black',
		notecex = 1.5,
		colsep = 1:3,
		rowsep = 1:3,
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
##make a volcano of the gene data and plot the protein markers on it
##########################################################
#need to go back to ccle and do the stats
pepQuan <- function(x,...){
	eset = as.matrix(x[,3:7])
	design <- cbind(SR=c(1,1,0,0,0),CC=c(0,0,1,1,1))
	fit <- lmFit(eset,design)
	cont.matrix <- makeContrasts(SRvsCC=SR-CC, levels=design)
	fit2 <- contrasts.fit(fit, cont.matrix)
	fit2 <- eBayes(fit2)
	x$logFC = fit2$coef
	x$pVal = fit2$p.value
	#output the data
	return(x[,c(1:2,11:12)])
}

ccle.q = pepQuan(ccle)

#output the data table
write.table(ccle.q,
		'ch_OvC_pubCellLine_CCLE_geneSet.txt',
		quote=FALSE,
		sep='\t',
		col.names=TRUE,
		row.names=FALSE)

###do the plotting
#make some colors
lnCols<-brewer.pal(6,"Blues")
cols<-rev(brewer.pal(6,"RdBu"))
#get the list of candidate genes
setwd(dir="/Users/cshughes/Documents/projects/OvC/")
genesCC = read.table("./Schwartz_ClearCellTable_geneList.txt", header=FALSE, sep='\t')
genesS = read.table("./Schwartz_SerousTable_geneList.txt", header=FALSE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/Routput/")
#make them into a searchable list
geneCC = as.character(genesCC[,1])
geneS = as.character(genesS[,1])
geneTot = c(geneCC,geneS)
#make the plotting function 
pdf('ch_OvC_TMT10_GENEexpression_volcano.pdf')
proSD<-sd(ccle.q[,3], na.rm=TRUE)
xCol = col2rgb(ifelse(ccle.q[,2] %in% geneCC, cols[1], ifelse(ccle.q[,2] %in% geneS, cols[6],'gray30')))
xCex = ifelse(ccle.q[,2] %in% geneTot, 2, 1)
#plot with all points
plot(ccle.q[,3],
		-log10(ccle.q[,4]),
		col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),
		pch=20,
		cex = xCex,
		ylab = '-log10(Adjusted p-value)',
		xlab = 'log2(serous/clear cell)',
		main = 'OvC type Cell Line comparison',
		xlim = c(-10,10)
)
box(lwd=3)
abline(v = -proSD, col=lnCols[6],lwd=2,lty=2)
abline(v = proSD, col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.05), col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.01), col=lnCols[6],lwd=2,lty=2)
text(8,0.5,paste('n=',nrow(ccle.q),sep=""),cex=1.25)
text(8,1.6,paste('p<0.05'),cex=1.25)
text(8,2.3,paste('p<0.01'),cex=1.25)
#plot with just subset points
xCol = ifelse(ccle.q[,2] %in% geneCC, cols[1], ifelse(ccle.q[,2] %in% geneS, cols[6],'gray30'))
gnRM = !grepl('gray', xCol)
vB = ccle.q[gnRM,]
xCol = col2rgb(ifelse(vB[,2] %in% geneCC, cols[1], ifelse(vB[,2] %in% geneS, cols[6],'gray30')))
plot(vB[,3],
		-log10(vB[,4]),
		col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),
		pch=20,
		cex = 2,
		ylab = '-log10(Adjusted p-value)',
		xlab = 'log2(serous/clear cell)',
		main = 'OvC type Cell Line comparison',
		xlim = c(-10,10)
)
box(lwd=3)
abline(v = -proSD, col=lnCols[6],lwd=2,lty=2)
abline(v = proSD, col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.05), col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.01), col=lnCols[6],lwd=2,lty=2)
text(8,0.5,paste('n=',sum(gnRM),sep=""),cex=1.25)
text(8,1.6,paste('p<0.05'),cex=1.25)
text(8,2.3,paste('p<0.01'),cex=1.25)
dev.off()


##########################################################
##make a volcano of the gene data and plot the tissue protein data on it
##########################################################
#use the pro data frame from the processing sheet
tis = pro[,c(1:5,7:8)]
#get the sd of the gene data
proSD<-sd(ccle.q[,3], na.rm=TRUE)
#make a mean expression column
tis$meanFC = rowMeans(tis[,c(4,6)], na.rm=TRUE)
tis$include = ifelse(tis[,8] > 1.25, 'yes', ifelse(tis[,8] < -1.25,'yes','no'))
#make vectors for color
upGNloc = tis$meanFC > 1.25
dnGNloc = tis$meanFC < -1.25
upGN = as.character(tis[upGNloc,2])
dnGN = as.character(tis[dnGNloc,2])
updnGN = c(upGN,dnGN)
#make the plotting function 
pdf('ch_OvC_TMT10_GENEexpressionPROTEINoverlay_volcano.pdf')
proSD<-sd(ccle.q[,3], na.rm=TRUE)
xCol = col2rgb(ifelse(ccle.q[,2] %in% dnGN, cols[1], ifelse(ccle.q[,2] %in% upGN, cols[6],'gray30')))
xCex = ifelse(ccle.q[,2] %in% updnGN, 2, 1)
#plot with all points
plot(ccle.q[,3],
		-log10(ccle.q[,4]),
		col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),
		pch=20,
		cex = xCex,
		ylab = '-log10(Adjusted p-value)',
		xlab = 'log2(serous/clear cell)',
		main = 'OvC type Cell Line comparison',
		xlim = c(-10,10)
)
box(lwd=3)
abline(v = -proSD, col=lnCols[6],lwd=2,lty=2)
abline(v = proSD, col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.05), col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.01), col=lnCols[6],lwd=2,lty=2)
text(8,0.5,paste('n=',nrow(ccle.q),sep=""),cex=1.25)
text(8,1.6,paste('p<0.05'),cex=1.25)
text(8,2.3,paste('p<0.01'),cex=1.25)
#plot with just subset points
xCol = ifelse(ccle.q[,2] %in% dnGN, cols[1], ifelse(ccle.q[,2] %in% upGN, cols[6],'gray30'))
gnRM = !grepl('gray', xCol)
vB = ccle.q[gnRM,]
xCol = col2rgb(ifelse(vB[,2] %in% dnGN, cols[1], ifelse(vB[,2] %in% upGN, cols[6],'gray30')))
plot(vB[,3],
		-log10(vB[,4]),
		col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),
		pch=20,
		cex = 2,
		ylab = '-log10(Adjusted p-value)',
		xlab = 'log2(serous/clear cell)',
		main = 'OvC type Cell Line comparison',
		xlim = c(-10,10)
)
box(lwd=3)
abline(v = -proSD, col=lnCols[6],lwd=2,lty=2)
abline(v = proSD, col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.05), col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.01), col=lnCols[6],lwd=2,lty=2)
text(8,0.5,paste('n=',sum(gnRM),sep=""),cex=1.25)
text(8,1.6,paste('p<0.05'),cex=1.25)
text(8,2.3,paste('p<0.01'),cex=1.25)
dev.off()

##########################################################
##make a volcano of the protein data and plot the cell line gene data on it
##########################################################
#use the pro data frame from the processing sheet
tis = pro[,c(1:5,7:8)]
#get the sd of the gene data
proSD<-sd(tis[,4], na.rm=TRUE)
#set up the gene data
cc = ccle.q
cc$include = ifelse(cc[,3] > 2, 'yes', ifelse(cc[,3] < -2,'yes','no'))
#make vectors for color
upGNloc = cc[,3] > 2
dnGNloc = cc[,3] < -2
upGN = as.character(cc[upGNloc,2])
dnGN = as.character(cc[dnGNloc,2])
updnGN = c(upGN,dnGN)
#make the plotting function 
pdf('ch_OvC_TMT10_PROTEINexpressionGENEoverlay_volcano.pdf')
proSD<-sd(tis[,4], na.rm=TRUE)
xCol = col2rgb(ifelse(tis[,2] %in% dnGN, cols[1], ifelse(tis[,2] %in% upGN, cols[6],'gray30')))
xCex = ifelse(tis[,2] %in% updnGN, 2, 1)
#plot with all points
plot(tis[,4],
		-log10(tis[,5]),
		col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),
		pch=20,
		cex = xCex,
		ylab = '-log10(Adjusted p-value)',
		xlab = 'log2(serous/clear cell)',
		main = 'OvC type Cell Line comparison',
		xlim = c(-5,5)
)
box(lwd=3)
abline(v = -proSD, col=lnCols[6],lwd=2,lty=2)
abline(v = proSD, col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.05), col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.01), col=lnCols[6],lwd=2,lty=2)
text(4,0.5,paste('n=',nrow(tis),sep=""),cex=1.25)
text(4,1.6,paste('p<0.05'),cex=1.25)
text(4,2.3,paste('p<0.01'),cex=1.25)
#plot with just subset points
xCol = ifelse(tis[,2] %in% dnGN, cols[1], ifelse(tis[,2] %in% upGN, cols[6],'gray30'))
gnRM = !grepl('gray', xCol)
vB = tis[gnRM,]
xCol = col2rgb(ifelse(vB[,2] %in% dnGN, cols[1], ifelse(vB[,2] %in% upGN, cols[6],'gray30')))
plot(vB[,4],
		-log10(vB[,5]),
		col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),
		pch=20,
		cex = 2,
		ylab = '-log10(Adjusted p-value)',
		xlab = 'log2(serous/clear cell)',
		main = 'OvC type Cell Line comparison',
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
##look at gene co-expression
##########################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/pubData/TCGA/")
wt1 = read.table("./coexpression_WT1_mRNA_expression_ov_tcga_pub.txt", header=TRUE, sep='\t')
k14 = read.table("./coexpression_KLHL14_mRNA_expression_ov_tcga_pub.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/Routput/")
vE = pro
#filter the data by 0.4 correlation
wt1.s = subset(wt1, Pearson.Score < -0.35 | Pearson.Score > 0.35)
k14.s = subset(k14, Pearson.Score < -0.35 | Pearson.Score > 0.35)
#make vectors for color
mkPlot <- function(x,...){
	upGNloc = x$Pearson.Score > 0.35
	dnGNloc = x$Pearson.Score < -0.35
	upGN = as.character(x[upGNloc,1])
	dnGN = as.character(x[dnGNloc,1])
	updnGN = c(upGN,dnGN)
	pdf('ch_OvC_TMT10_WT1_coxpress.pdf')
	proSD<-sd(vE[,4], na.rm=TRUE)
	xCex = ifelse(vE[,2] %in% updnGN, 2, 1)
	#plot with just subset points
	xCol = ifelse(vE[,2] %in% dnGN, cols[1], ifelse(vE[,2] %in% upGN, cols[6],'gray30'))
	gnRM = !grepl('gray', xCol)
	vB = vE[gnRM,]
	xCol = col2rgb(ifelse(vB[,2] %in% dnGN, cols[1], ifelse(vB[,2] %in% upGN, cols[6],'gray30')))
	plot(vB[,4],
		-log10(vB[,5]),
		col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),
		pch=20,
		cex = 2,
		ylab = '-log10(Adjusted p-value)',
		xlab = 'log2(serous/clear cell)',
		main = 'OvC CoExpression TCGA',
		xlim = c(-5,5),
		ylim = c(0,4)
)
	box(lwd=3)
	abline(v = -proSD, col=lnCols[6],lwd=2,lty=2)
	abline(v = proSD, col=lnCols[6],lwd=2,lty=2)
	abline(h = -log10(0.05), col=lnCols[6],lwd=2,lty=2)
	abline(h = -log10(0.01), col=lnCols[6],lwd=2,lty=2)
	text(4,0.5,paste('n=',sum(gnRM),sep=""),cex=1.25)
	text(4,1.6,paste('p<0.05'),cex=1.25)
	text(4,2.3,paste('p<0.01'),cex=1.25)
	text(vB[,4],-log10(vB[,5]),vB[,2],pos=3)
dev.off()
}

mkPlot(wt1.s)




##########################################################
##use the whole ccle data set and cluster based on consolidated markers
##########################################################
#import the rankings
setwd(dir="/Users/cshughes/Documents/projects/OvC/pubData/")
sr = read.table("./OvC_CellLineSuitabilityRank.txt", header=TRUE, sep='\t')
sc = read.table("./OvC_CellLineSuitabilityRank_wTCGA.txt", header=TRUE, sep='\t')
ccle.ov = read.table("./CCLE_Expression_OVonly_2012-09-29.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/Routput/")

##get the marker list out
ser = subset(tcg.q, -log10(pAdj) > -log10(0.05) & meanFC > 0)
ser$class = 'serous'
cc = subset(tcg.q, -log10(pAdj) > -log10(0.05) & meanFC < 0)
cc$class = 'clear'
ser.cc = rbind(ser,cc)

#modify the ccle data
colnames(ccle.ov)[2] = 'Gene'
exprs = merge(ser.cc[,c(1,7)],ccle.ov,by='Gene')
exprs = exprs[order(exprs$class),]
#remove the origin from the column names
colnames(exprs)[4:ncol(exprs)] = sapply(strsplit(colnames(exprs[,4:ncol(exprs)]), '_'),'[', 1)
#rework the data to include the cell line rankings
tex = data.frame(t(exprs[,4:ncol(exprs)]))
colnames(tex) = exprs[,1]
tex = cbind(Row.Names = rownames(tex), tex)
rownames(tex) <- NULL
colnames(tex)[1] = 'CellLine'
#merge with the ranks
tex.sr = merge(sr,tex,by='CellLine')
tex.sr = tex.sr[order(-tex.sr$Rank),]
cols = brewer.pal(6,'Set2')
Rowcols = ifelse(tex.sr$Type=='NS',cols[1],ifelse(tex.sr$Type=='Serous',cols[2],ifelse(tex.sr$Type=='Endometrioid',cols[3],ifelse(tex.sr$Type=='Mixed',cols[4],ifelse(tex.sr$Type=='Clear',cols[5],ifelse(tex.sr$Type=='Mucinous',cols[6],'black'))))))

#make a heatmap of the expression patterns
pdf('ch_OvCall_consolidatedMarkers_ccle.pdf')
#make the plot labels and boundaries
mybreaks = seq(3,10,by=0.1) 
#make the correlation heatmap
heatmap.2(
		as.matrix(tex.sr[,4:ncol(tex.sr)]),
		col= colorRampPalette(brewer.pal(9,"YlOrRd"))(length(mybreaks)-1),
		symkey=FALSE,
		Rowv=FALSE,
		Colv=FALSE,
		#hclust=hclustfunc,
		#distfun=distfunc,
		#na.color='black',
		na.rm=TRUE,
		#symm=FALSE,
		dendrogram="none",
		breaks=mybreaks,
		#cellnote = round(ovCor,2),
		labRow = tex.sr[,1],
		labCol = colnames(tex.sr[,4:ncol(tex.sr)]),
		notecol = 'black',
		notecex = 1.5,
		#colsep = 1:70,
		#rowsep = 1:70,
		sepwidth = c(0.03, 0.03),
		sepcolor = 'white',
		RowSideColors= Rowcols,
		## labels
		main='CCLE',
		#xlab='Time Point',
		#ylab='Protein',
		## color key
		key = TRUE,
		keysize = 1,
		density.info = "none",
		scale = "none",
		trace = "none",
		mar=c(8,8),
		cexRow=0.75,
		cexCol=0.75
)
dev.off()












