# TODO: Add comment
# 
# Author: cshughes
###############################################################################
#use proCell and proTis for these analyses
##################################################
#read in the human protein atlas data
##################################################
setwd(dir="/Users/cshughes/Documents/projects/PaC/RNAseq/")
ens<-read.table("./HumanProteomeMap_Science_TableS1.txt", header=TRUE, sep='\t')
hpa<-read.table("./HumanProteinAtlas_FPKM_allsamples.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/Routput/")
#change the gene colname
colnames(ens)[2] = 'Gene'
colnames(hpa)[1] = 'ensg_id'
#merge the two human atlas sets
hpa.m = merge(ens,hpa, by='ensg_id')
#remove isoforms
isoforms = grepl('\\.', hpa.m$Gene)
isoforms2 = grepl('\\-', hpa.m$Gene)
hpa.mi = hpa.m[!isoforms,]
hpa.mi = hpa.mi[!isoforms2,]
#look at the FPKM distribution
pdf('ch_HumanProteinAtlas_FPKM_Ovary.pdf')
hist(log2(hpa.mi[,130]),
		col = 'blue',
		breaks=200,
		xlim = c(-10,15)
)
#hist(log2(hpa.mi[,130]),add=TRUE,breaks=200,col='red')
dev.off()
#combine with the protein data
ov.whpa = merge(proTis,hpa.mi,by='Gene',all=FALSE,sort=FALSE)
#plot the FPKM overlaid with the protein data
cols = brewer.pal(6,'Set2')
pdf('ch_HumanProteinAtlas_FPKMwProtein_Ovary_ctpTissue.pdf')
hist(log2(hpa.mi[,130]),
		col = cols[3],
		breaks=200,
		xlim = c(-10,15),
		main = 'Ovarian Tissue FPKM Distribution',
		xlab = 'log2(FPKM Pancreatic Tissue)',
		ylab = 'Frequency'
)
hist(log2(ov.whpa[,138]),add=TRUE,breaks=100,col=cols[6])
dev.off()


##########################################################
##make validation plots of significant trends to eliminate bias
##########################################################
pdf('ch_OvC_SignificanceVpepNumber_repAB_ctpTissue.pdf')
vE = proTis
vE$fcStd = rowMeans(vE[,c(4,7)],na.rm=TRUE)
vE$pvStd = rowMeans(vE[,c(5,8)],na.rm=TRUE)
vE$pepStd = rowMeans(vE[,c(6,9)])
lnCols<-brewer.pal(6,"Set3")
M= vE[,10]
A= log2(vE[,12])
heatscatter(A,M,
		cex = 1,
		xlab = 'log10(peptide number)',
		ylab = 'log2(Serous v Clear Cell)',
		main = 'fold change independence on peptide number',
		ylim = c(-4,4)
)
abline(h= sd(M,na.rm=TRUE), col=lnCols[1],lwd=3,lty=2)
abline(h= -sd(M,na.rm=TRUE), col=lnCols[1],lwd=3,lty=2)
#
M= -log10(vE[,11])
A= log10(vE[,12])
heatscatter(A,M,
		cex=1.15,
		xlab = 'log10(peptide number)',
		ylab = '-log10(p-value)',
		main = 'p-value independence on peptide number')
abline(h= -log10(0.05), col=lnCols[1],lwd=3,lty=2)
abline(h= -log10(0.01), col=lnCols[1],lwd=3,lty=2)
dev.off()

##########################################################
##make volcano plots of expression variance
##########################################################
#make some colors
lnCols<-brewer.pal(6,"Blues")
pogCols<-brewer.pal(6, "Greens")
cols<-rev(brewer.pal(6,"RdBu"))
vE = proTis
#get the list of candidate genes
setwd(dir="/Users/cshughes/Documents/projects/OvC/")
genesCC = read.table("./Schwartz_ClearCellTable_geneList.txt", header=FALSE, sep='\t')
genesS = read.table("./Schwartz_SerousTable_geneList.txt", header=FALSE, sep='\t')
pog = read.table("./POGMarkers.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/Routput/")
#
geneCC = as.character(genesCC[,1])
geneS = as.character(genesS[,1])
geneTot = c(geneCC,geneS)
colnames(pog) = c('ensg_id','Gene')
pogM = as.character(pog[,2])
#make the plotting function 
pdf('ch_OvC_Expression-POGmarkers_ctpTissue.pdf')
vE$tisStd = rowMeans(vE[,c(4,7)],na.rm=TRUE)
proSD<-sd(vE[,4], na.rm=TRUE)
xCol = col2rgb(ifelse(vE[,2] %in% pogM, pogCols[5], ifelse(vE[,2] %in% geneCC, cols[1],ifelse(vE[,2] %in% geneS, cols[6],'gray30'))))
xCex = ifelse(vE[,2] %in% pogM, 2, ifelse(vE[,2] %in% geneTot, 2,1))
#xCol = col2rgb(ifelse(-log10(vE[,5]) > -log10(0.01),ifelse(vE[,4] > proSD, cols[1],ifelse(vE[,4] < -proSD, cols[1], 'gray30')),'gray30'))
#plot with all points
plot(vE[,10],
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
xCol = ifelse(vE[,2] %in% pogM, pogCols[5], ifelse(vE[,2] %in% geneCC, cols[1],ifelse(vE[,2] %in% geneS, cols[6],'gray30')))
gnRM = !grepl('gray', xCol)
vB = vE[gnRM,]
xCol = col2rgb(ifelse(vB[,2] %in% pogM, pogCols[5], ifelse(vB[,2] %in% geneCC, cols[1],ifelse(vB[,2] %in% geneS, cols[6],'gray30'))))
plot(vB[,10],
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


##################################################
#using the human protein atlas data for proteins
##################################################
setwd(dir="/Users/cshughes/Documents/projects/PaC/RNAseq/")
ens<-read.table("./HumanProteomeMap_Science_TableS1.txt", header=TRUE, sep='\t')
hpa<-read.table("./HumanProteinAtlas_ProteinExp.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/Routput/")
####use the ov.a2 set from the ch_OvC_technicalReplicate.R script
ov.a2$pepNum = 1
ov.ap = aggregate(cbind(specNum,pepNum)~Accession+Gene+Descriptions, data=ov.a2, sum, na.action=na.pass, na.rm=TRUE)
ov.ap$ab = ov.ap$specNum/ov.ap$pepNum
#subset out only pancreatic tissue samples
hpa.p = hpa[grepl('ovar',hpa$Tumor),]
#x = aggregate(Count.patients~ensg_id,data=hpa.p,sum,na.action=na.pass,na.rm=TRUE)
#y = merge(x,pc.whpa,by='ensg_id')
#get it to a single measurement per gene based on a 60% threshold of expression
hpa.pe = hpa.p[hpa.p$Count.patients/hpa.p$Total.patients >= 0.6,]
#merge with the protein and gene expression data
ov.whpa = merge(ov.ap,hpa.mi,by='Gene',all=FALSE,sort=FALSE)
colnames(hpa.pe)[1] = 'ensg_id'
ov.all = merge(ov.whpa,hpa.pe, by='ensg_id')
ov.heat = ov.all[,c(2,145,7)]
#for getting numbers
x = aggregate(Count.patients~Gene, hpa.p, sum, na.action=na.pass, na.rm=TRUE)
colnames(x)[1] = 'ensg_id'
y = merge(ov.whpa,x,by='ensg_id')
#make the plot
pdf('ch_OvC_HumanProteinAtlas_Heat.pdf')
ov.heat = ov.heat[order(-ov.heat$ab),]
cols = brewer.pal(9,'YlOrRd')
ov.heat$sideCols = ifelse(grepl('High',ov.heat$Level), cols[9], ifelse(grepl('Medium',ov.heat$Level), cols[5], ifelse(grepl('Low',ov.heat$Level), brewer.pal(6,'RdBu')[6], 'gray30')))
#nd = grepl('Not',ov.heat$Level)
#ov.heatSub = ov.heat[!nd,]
x = as.matrix(log2(ov.heat[,3]))
y = cbind(x,x)
#make the plot labels and boundaries
xLabels<- names(ov.heat)[1]
mybreaks = seq(0,3,by=0.05) 
#make the correlation heatmap
heatmap.2(
		y,
		col= rev(colorRampPalette(brewer.pal(9,"RdBu"))(length(mybreaks)-1)),
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
		#cellnote = round(x,2),
		labRow = '',
		labCol = '',
		#notecol = 'black',
		#notecex = 0.75,
		colsep = 1,
		rowsep = 1,
		sepwidth = c(0.03, 0.03),
		sepcolor = 'white',
		RowSideColors=ov.heat$sideCols,
		## labels
		main='correlation map',
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

#try something else here....
####use the ov.psmh set from the ch_OvC_technicalReplicate.R script
ov.a2 = aggregate(cbind(specNum,Average.Reporter.SN)~Accession+Gene+Descriptions+Sequence, data=ov.psmh, sum, na.action=na.pass, na.rm=TRUE)
ov.a2$pepNum = 1
ov.ap = aggregate(cbind(specNum,pepNum,Average.Reporter.SN)~Accession+Gene+Descriptions, data=ov.a2, sum, na.action=na.pass, na.rm=TRUE)
ov.ap$ab = (ov.ap$specNum * ov.ap$Average.Reporter.SN) / ov.ap$pepNum
ov.whpa = merge(ov.ap,hpa.mi,by='Gene',all=FALSE,sort=FALSE)
colnames(hpa.pe)[1] = 'ensg_id'
ov.all = merge(ov.whpa,hpa.pe, by='ensg_id',all=TRUE)
ov.heat = ov.all[,c(1:2,146,8)]

#try to make a simple scatter of the tiered data
library(beeswarm)
pdf('ch_OvC_HumanProteinAtlas_VioPlot.pdf')
x = ov.heat[!is.na(ov.heat$ab),3:4]
x = x[x$ab>0,]
x$ab = log2(x$ab)
x$cat = 'Not in HPA'
high = grepl('High', x[,1])
low = grepl('Low', x[,1])
med = grepl('Medium', x[,1])
nd = grepl('Not', x[,1])
x[high,3]='High'
x[low,3]='Low'
x[med,3]='Med'
x[nd,3]='ND'
x1 <- x$ab[high]
x2 <- x$ab[low]
x3 <- x$ab[med]
x4 <- x$ab[nd]
x5 <- x$ab[nd]
boxplot(x[,2]~x[,3],
		horizontal=TRUE)
vioplot(x1, x2, x3, x4, x5,  
		col=brewer.pal(5,'Set3'),
		horizontal=TRUE,
		add=TRUE)
boxplot(x[,2]~x[,3],
		horizontal=TRUE,
		add=TRUE)
dev.off()





##################################################
#using the human protein atlas data for FPKM correlation
##################################################
setwd(dir="/Users/cshughes/Documents/projects/PaC/RNAseq/")
ens<-read.table("./HumanProteomeMap_Science_TableS1.txt", header=TRUE, sep='\t')
hpa<-read.table("./HumanProteinAtlas_FPKM_allsamples.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/Routput/")
#change the gene colname
colnames(ens)[2] = 'Gene'
colnames(hpa)[1] = 'ensg_id'
#merge the two human atlas sets
hpa.m = merge(ens,hpa, by='ensg_id')
#remove isoforms
isoforms = grepl('\\.', hpa.m$Gene)
isoforms2 = grepl('\\-', hpa.m$Gene)
hpa.mi = hpa.m[!isoforms,]
hpa.mi = hpa.mi[!isoforms2,]
####use the ov.psmh set from the ch_OvC_technicalReplicate.R script
ov.a2 = aggregate(cbind(specNum,Average.Reporter.SN)~Accession+Gene+Descriptions+Sequence, data=ov.psmh, sum, na.action=na.pass, na.rm=TRUE)
ov.a2$pepNum = 1
ov.ap = aggregate(cbind(specNum,pepNum,Average.Reporter.SN)~Accession+Gene+Descriptions, data=ov.a2, sum, na.action=na.pass, na.rm=TRUE)
ov.ap$ab = (ov.ap$specNum * ov.ap$Average.Reporter.SN) / ov.ap$pepNum
ov.whpa = merge(ov.ap,hpa.mi,by='Gene',all=FALSE,sort=FALSE)

pdf('ch_OvC_FPKM-ProteinExpressionCor_ctpTissue.pdf')
lnCols<-brewer.pal(6,"Blues")
x = subset(ov.whpa, ov.whpa[,135]>0)
x = subset(x, x$ab>0)
reg = lm(log2(x[,135])~log2(x$ab))
heatscatter(log2(x$ab),
		log2(x[,135]),
		ylim = c(-10,15),
		xlab = 'log2(HPA FPKM Ovary)',
		ylab = 'MS Expression Estimate')
abline(reg,col=lnCols[6],lwd=2,lty=2)
text(20,-10,paste('n=',nrow(x),sep=""))
text(20,-7,paste('cor=',round(cor(log2(x$ab),log2(x[,135]),method='spearman'),3),sep=""))
dev.off()



