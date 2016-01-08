# TODO: Add comment
# 
# Author: cshughes
###############################################################################


##########################################################
##make a plot for peptide metrics
##########################################################
#subset out peptides with greater than 50 peptides assigned
vE = proCell
nPeps<-data.frame(table(vE[,6]))
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
pdf('ch_OvCcell-TMT10_PeptideIDMetrics.pdf')
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

##########################################################
##make volcano plots of marker gene expression
##########################################################
#make some colors
lnCols<-brewer.pal(6,"Blues")
cols<-rev(brewer.pal(6,"RdBu"))
vE = proCell
#get the list of candidate genes
setwd(dir="/Users/cshughes/Documents/projects/OvC/")
genesCC = read.table("./Schwartz_ClearCellTable_geneList.txt", header=FALSE, sep='\t')
genesS = read.table("./Schwartz_SerousTable_geneList.txt", header=FALSE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/cell/Routput/")

geneCC = as.character(genesCC[,1])
geneS = as.character(genesS[,1])
geneTot = c(geneCC,geneS)

#make the plotting function 
pdf('ch_OvCcell_TMT10_expression_volcano_CC.pdf')
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
		xlim = c(-7,7)
)
box(lwd=3)
abline(v = -proSD, col=lnCols[6],lwd=2,lty=2)
abline(v = proSD, col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.05), col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.01), col=lnCols[6],lwd=2,lty=2)
text(4,0.5,paste('n=',nrow(vE),sep=""),cex=1.25)
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
		xlim = c(-7,7)
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
##make a volcano of the cell line data and plot the tissue protein data on it
##########################################################
#use the pro data frame from the processing sheet
tis = pro[,c(1:5,7:8)]
#get the sd of the gene data
proSD<-sd(proCell[,4], na.rm=TRUE)
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
pdf('ch_OvC_TMT10_CELLexpressionTISSUEoverlay_volcano.pdf')
proSD<-sd(proCell[,4], na.rm=TRUE)
xCol = col2rgb(ifelse(proCell[,2] %in% dnGN, cols[1], ifelse(proCell[,2] %in% upGN, cols[6],'gray30')))
xCex = ifelse(proCell[,2] %in% updnGN, 2, 1)
#plot with all points
plot(proCell[,4],
		-log10(proCell[,5]),
		col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),
		pch=20,
		cex = xCex,
		ylab = '-log10(Adjusted p-value)',
		xlab = 'log2(serous/clear cell)',
		main = 'OvC type Cell Line comparison',
		xlim = c(-6,6)
)
box(lwd=3)
abline(v = -proSD, col=lnCols[6],lwd=2,lty=2)
abline(v = proSD, col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.05), col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.01), col=lnCols[6],lwd=2,lty=2)
text(4,0.5,paste('n=',nrow(proCell),sep=""),cex=1.25)
text(4,1.6,paste('p<0.05'),cex=1.25)
text(4,2.3,paste('p<0.01'),cex=1.25)
#plot with just subset points
xCol = ifelse(proCell[,2] %in% dnGN, cols[1], ifelse(proCell[,2] %in% upGN, cols[6],'gray30'))
gnRM = !grepl('gray', xCol)
vB = proCell[gnRM,]
xCol = col2rgb(ifelse(vB[,2] %in% dnGN, cols[1], ifelse(vB[,2] %in% upGN, cols[6],'gray30')))
plot(vB[,4],
		-log10(vB[,5]),
		col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),
		pch=20,
		cex = 2,
		ylab = '-log10(Adjusted p-value)',
		xlab = 'log2(serous/clear cell)',
		main = 'OvC type Cell Line comparison',
		xlim = c(-6,6)
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
##plot an MA of the tissue and cell line data
##########################################################
vC = proCell
vT = pro[,c(1:5,7:8)]
vT$meanFC = rowMeans(vT[,c(4,6)], na.rm=TRUE)
#make a single data set
vA = merge(vC,vT,by='Gene')
vA = vA[,c(4,6,13)]
#make the plot
cols = colorRampPalette(c("white", brewer.pal(9,'YlOrRd')))
pdf('ch_OvC_TissueANDCell_MAplot.pdf')
M= vA[,1]
A= vA[,3]
smoothScatter(A,M,colramp = cols, cex = 2, xlab = 'Average Intensity',ylab = 'Intensity Ratio',main='preNorm expression')
abline(h=0,col=brewer.pal(6,"RdBu")[6],lwd=3)
box(lwd=3)
dev.off()



##########################################################
##capture the markers present in all 3 sets
##########################################################




